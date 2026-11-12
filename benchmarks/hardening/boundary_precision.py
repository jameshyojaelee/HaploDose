"""Extended-precision replay of fixed posterior-boundary RNA profiles."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose import profile
from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('boundary_precision_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def logsum(x, axis=None):
    high = np.max(x, axis=axis, keepdims=True)
    value = high+np.log(np.sum(np.exp(x-high), axis=axis, keepdims=True,
                              dtype=np.longdouble))
    return value.squeeze() if axis is None else np.squeeze(value, axis=axis)


def runs(emission, prior, change):
    # Keep the original binary64 transition entries and normalized prior as
    # model inputs. Extend arithmetic, and accumulate each interval directly
    # instead of subtracting two cumulative sums.
    e = np.asarray(emission, np.longdouble)
    prior = np.asarray(prior, float)
    prior = prior/prior.sum()
    transition = load('shared_membership_hmm').transition(prior, change)
    T = np.log(transition.astype(np.longdouble))
    n, states = e.shape
    f, b = np.empty_like(e), np.empty_like(e)
    f[0] = np.log(prior.astype(np.longdouble))+e[0]
    for j in range(1,n):
        f[j] = e[j]+logsum(f[j-1,:,None]+T, axis=0)
    b[-1] = 0
    for j in range(n-2,-1,-1):
        b[j] = logsum(T+e[j+1,None,:]+b[j+1,None,:], axis=1)
    evidence = logsum(f[-1])
    off = T.copy()
    np.fill_diagonal(off, -np.inf)
    left, right = np.empty_like(e), np.empty_like(e)
    left[0], right[-1] = np.log(prior.astype(np.longdouble)), 0
    for j in range(1,n):
        left[j] = logsum(f[j-1,:,None]+off, axis=0)
    for j in range(n-1):
        right[j] = logsum(off+e[j+1,None,:]+b[j+1,None,:], axis=1)
    result = np.zeros((n,n), np.longdouble)
    for a in range(n):
        accumulated = np.zeros(states, np.longdouble)
        for z in range(a,n):
            accumulated += e[z]
            result[a,z] = np.exp(logsum(left[a]+accumulated+(z-a)*np.diag(T)+right[z]-evidence))
    mass_error = max(abs(result[:j+1,j:].sum(dtype=np.longdouble)-1) for j in range(n))
    return result, float(mass_error)


def checks():
    prior,change = np.array([.25,.5,.25]),.125
    likelihood = np.array([[.5,.25,.125],[.125,.5,.25],[.25,.125,.5]])
    T = load('shared_membership_hmm').transition(prior,change)
    paths = list(itertools.product(range(3),repeat=3))
    mass = []
    for path in paths:
        value = Fraction(float(prior[path[0]]))
        for j,state in enumerate(path):
            value *= Fraction(float(likelihood[j,state]))
            if j:value *= Fraction(float(T[path[j-1],state]))
        mass.append(value)
    total = sum(mass)
    expected = np.zeros((3,3),np.longdouble)
    for path,value in zip(paths,mass):
        p = value/total
        p = np.longdouble(str(p.numerator))/np.longdouble(str(p.denominator))
        for a in range(3):
            if a and path[a]==path[a-1]:continue
            b = a
            while b+1<3 and path[b+1]==path[a]:b+=1
            expected[a,b] += p
    actual,error = runs(np.log(likelihood.astype(np.longdouble)),prior,change)
    operations = 128*3*3*3+256
    eps = np.finfo(np.longdouble).eps
    budget = operations*eps/(1-operations*eps)*32
    discrepancy = float(np.max(abs(expected-actual)))
    assert max(discrepancy,error)<=budget
    return dict(exact_fraction_paths=27,maximum_probability_error=discrepancy,
                operation_rounding_budget=float(budget))


def replay(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    original = args.task/(args.mode+'_v1')
    inputs = {k:Path(v['path']) for k,v in
              json.loads((original/'run_manifest.json').read_text())['input']['files'].items()}
    inputs.update(script=Path(__file__), precision_design=args.task/'PRECISION_DESIGN.md',
                  original_predictions=original/'RNA_predictions.npz')
    rec = start_run(args.out, inputs, dict(mode=args.mode),
        dict(scope='Fixed model numerical replay; no model or threshold selection'),
        mode='boundary_extended_precision')
    try:
        assert np.finfo(np.longdouble).nmant > np.finfo(float).nmant
        numerical = checks()
        genes,bins,samples = [pd.read_csv(inputs[k],sep='\t') for k in ('genes','bins','samples')]
        effects = np.array([np.load(inputs[f'RNA{i}'])['corrected_gene'] for i in range(68)])
        with np.load(inputs['weights']) as weights:
            variance = weights['tumour_variance']
        decisions = {r['sample']:r for r in json.loads(inputs['decisions'].read_text())}
        reset = chromosome_starts(bins)|arm_starts(bins,genome_build='hg38')
        start = np.flatnonzero(reset)
        ends = np.r_[start[1:],len(bins)]
        params = json.loads(inputs['native_parameters'].read_text())['params']
        normal = samples.role.eq('normal').to_numpy()
        result = np.full((68,len(bins)), np.nan, np.longdouble)
        records = []
        average = load('posterior_region_average')
        for k,i in enumerate(np.flatnonzero(~normal)):
            effect,_,_ = load('reference_weighted_regions').normal_centered_effect(
                effects[i],effects[normal],genes.bin.to_numpy(int))
            fitted = decisions[str(samples.iloc[i]['sample'])]
            with np.load(inputs[f'fit{i}']) as z:
                y = z['bin_values'].astype(np.longdouble)
                total = z['total_variance'].astype(np.longdouble)
            means = np.longdouble(fitted['beta'])*np.asarray(profile.LOG_DOSE,np.longdouble)
            # State-independent Gaussian constants cancel from every path.
            emission = -(y[:,None]-means)**2/(2*total[:,None])
            errors = []
            for a,b in zip(start,ends):
                probability,error = runs(emission[a:b],profile.STATE_PRIOR,params['p_change'])
                keep = (genes.bin.to_numpy()>=a)&(genes.bin.to_numpy()<b)
                summaries = average.interval_summaries(effect[keep],variance[k,keep],
                    genes.bin.to_numpy(int)[keep]-a,b-a)[1].astype(np.longdouble)
                for j in range(b-a):
                    result[i,a+j] = np.sum(probability[:j+1,j:]*summaries[:j+1,j:],dtype=np.longdouble)
                errors.append(error)
            records.append(dict(sample=str(samples.iloc[i]['sample']),maximum_mass_error=max(errors)))
        name = ('averaged_gaussian_reference_median' if args.mode=='control'
                else 'decomposed_averaged_reference_median')
        with np.load(inputs['original_predictions']) as z:
            difference = float(np.nanmax(abs(result-z[name].astype(np.longdouble))))
        np.savez_compressed(args.out/'RNA_predictions.npz',primary=result)
        write_json(dict(mode=args.mode,mantissa_bits=np.finfo(np.longdouble).nmant,
            maximum_profile_change=difference,cases=records,model_parameters_changed=False),
            args.out/'report.json')
        write_json(numerical,args.out/'numerics.json')
        finish_run(args.out,rec,backend='CPU extended precision direct interval accumulation')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


def compare(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    original = args.task/'score_v1'
    inputs = {k:Path(v['path']) for k,v in
              json.loads((original/'run_manifest.json').read_text())['input']['files'].items()}
    inputs.update(script=Path(__file__),precision_design=args.task/'PRECISION_DESIGN.md',
        old_windows=original/'window_scores.npz',
        control=args.task/'precision_control_v1/RNA_predictions.npz',
        candidate=args.task/'precision_candidate_v1/RNA_predictions.npz')
    rec = start_run(args.out,inputs,dict(primary='decomposed_averaged_reference_median'),
        dict(scope='Numerical sign check only; no uncertainty/biological claim'),
        mode='boundary_precision_compare')
    try:
        genes,samples,mapping = [pd.read_csv(inputs[k],sep='\t')for k in ('genes','samples','mapping')]
        with np.load(inputs['old_windows']) as z:
            truth = {k:z[k]for k in z.files}
        row = pd.Index(samples['sample']).get_indexer(truth['samples'])
        assert np.all(row>=0)
        positions = genes.set_index('gene').bin.reindex(mapping.gene)
        has = positions.notna().to_numpy()
        assert has[mapping.shared.to_numpy(bool)].all()
        correlations, changes = {},{}
        for label,name in (('control','averaged_gaussian_reference_median'),
                           ('candidate','decomposed_averaged_reference_median')):
            with np.load(inputs[label]) as z:
                G = np.full((len(mapping),len(row)),np.nan,np.longdouble)
                G[has] = z['primary'][row][:,positions[has].to_numpy(int)].T
            G[~mapping.shared.to_numpy(bool)] = np.nan
            windows = np.full(truth['dna'].shape,np.nan,np.longdouble)
            for w in range(len(windows)):
                windows[w] = np.nanmedian(G[mapping.window.to_numpy(int)==w],axis=0)
            # Preserve the frozen scorer's binary32 window export separately.
            unchanged = np.array_equal(windows.astype(np.float32),truth[name],equal_nan=True)
            changes[label] = dict(original_float32_windows_exact=unchanged,
                maximum_window_change=float(np.nanmax(abs(windows-truth[name].astype(np.longdouble)))))
            values = []
            for j in np.flatnonzero(truth['split']=='evaluation'):
                mask = truth['common'][:,j]
                x,y = windows[mask,j],truth['dna'][mask,j].astype(np.longdouble)
                x,y = x-x.mean(),y-y.mean()
                values.append(np.sum(x*y)/np.sqrt(np.sum(x*x)*np.sum(y*y))
                    if len(x)>=3 and np.any(x)and np.any(y) else np.longdouble(np.nan))
            correlations[label] = np.array(values,np.longdouble)
        difference = correlations['candidate']-correlations['control']
        scores = pd.read_csv(original/'patient_scores.tsv',sep='\t',float_precision='round_trip')
        matrix = scores[scores.split.eq('evaluation')].pivot(index='sample',columns='method',values='pearson')
        before = float((matrix.decomposed_averaged_reference_median-matrix.averaged_gaussian_reference_median).mean())
        after = float(np.nanmean(difference))
        assert np.isfinite(difference).sum()==23
        assert math.copysign(1.,before)==math.copysign(1.,after)
        write_json(dict(defined_patients=23,original_mean_difference=before,
            extended_mean_difference=after,sign_unchanged=True,windows=changes,
            precision_establishes_scientific_superiority=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU extended Pearson and frozen binary32 window comparison')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('control','candidate','compare'),required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.mode=='compare':compare(args)
    else:replay(args)


if __name__=='__main__':main()
