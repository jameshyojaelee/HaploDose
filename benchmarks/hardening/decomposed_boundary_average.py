"""Research normal-noise decomposition composed with exact boundary averaging."""
import argparse
from functools import lru_cache
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr

from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.profile import ProfileParams
from haplodose.provenance import start_run, finish_run, file_record


NAMES = ('decomposed_averaged_uniform_median','decomposed_averaged_reference_median',
         'decomposed_averaged_reference_mean')
OLD = ('averaged_gaussian_uniform_median','averaged_gaussian_reference_median',
       'averaged_gaussian_reference_mean')


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('decomposed_average_'+name,Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checks():
    average,noise = load('posterior_region_average'),load('normal_noise_decomposition')
    y,levels = np.array([-.25,.2,.35]),np.array([-.4,0.,.4])
    observed,sampling,query = np.array([.25,.5,.25]),np.array([.125,.25,.125]),np.array([.125,.25,.125])
    prior,change = np.array([.2,.6,.2]),.1
    effects,weights,bins = np.array([-.5,0.,0.,.4,.2,.5]),np.array([1.,4.,1.,4.,1.,4.]),np.repeat(np.arange(3),2)
    summaries = average.interval_summaries(effects,weights,bins,3)
    results,errors,budgets = [],[],[]
    for variance in (observed+query,noise.component_variance(observed,sampling)+query):
        emission = np.array([[norm.logpdf(value,mean,math.sqrt(v)) for mean in levels]
                             for value,v in zip(y,variance)])
        probability,record = average.run_probabilities(emission,prior,change)
        paths = list(itertools.product(range(3),repeat=3))
        mass,profiles = [],[]
        for path in paths:
            p = prior[path[0]]
            for j,state in enumerate(path):
                p *= norm.pdf(y[j],levels[state],math.sqrt(variance[j]))
                if j:p *= 1-change if path[j-1]==state else change*prior[state]/(1-prior[path[j-1]])
            mass.append(p)
            profile = np.zeros((3,3))
            for a in range(3):
                if a and path[a]==path[a-1]:continue
                b=a
                while b+1<3 and path[b+1]==path[a]:b+=1
                keep=(bins>=a)&(bins<=b)
                med=load('reference_weighted_regions').weighted_median(effects[keep],1/np.sqrt(weights[keep]))
                mean=math.fsum(map(float,effects[keep]/weights[keep]))/math.fsum(map(float,1/weights[keep]))
                profile[:,a:b+1]=np.array([np.median(effects[keep]),med,mean])[:,None]
            profiles.append(profile)
        mass=np.array(mass)/math.fsum(mass)
        independent=sum(p*v for p,v in zip(mass,profiles))
        actual=np.array([[math.fsum(map(float,(probability[:j+1,j:]*summaries[k,:j+1,j:]).ravel()))
                         for j in range(3)]for k in range(3)])
        error=float(np.max(abs(actual-independent)))
        budget=record['mass_rounding_budget']*max(1.,float(np.max(abs(effects))))
        assert error<=budget
        results.append(actual);errors.append(error);budgets.append(budget)
    effect=float(np.max(abs(results[0]-results[1])))
    assert effect>sum(budgets)
    np.testing.assert_array_equal(noise.component_variance(observed,np.zeros(3))+query,observed+query)
    return dict(exhaustive_paths_each=27,full_density_errors=errors,rounding_budgets=budgets,
        variance_changes_averaged_profile_by=effect,no_normal_noise_limit_exact=True,
        existing_average_checks=average.checks(),existing_decomposition_checks=noise.checks(),
        calibration_established=False)


def predict(args):
    root=Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior=root/'dev/programme_gene_profile_20261109/refit_000_v1'
    decomp=root/'dev/normal_noise_decomposition_20261109/bulk_v1'
    old=root/'dev/posterior_region_average_20261109/bulk_v1'
    source=json.loads((prior/'run_manifest.json').read_text())['input']['files']
    native=Path(source['genes.tsv.gz']['path']).parent
    contexts=(root/'dev/reference_gaussian_regions_20261109/bulk_v1'if args.mode=='control'else decomp)
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',samples=Path(source['samples']['path']),
        genes=native/'genes.tsv.gz',bins=native/'bins.tsv',native_parameters=Path(source['run.json']['path']),
        decomposition=decomp/'normal_noise_decomposition.npz',weights=decomp/'reference_weights.npz',
        decisions=contexts/'Gaussian_decisions.json',old_RNA=old/'RNA_predictions.npz',
        averaging=Path(__file__).with_name('posterior_region_average.py'),
        noise=Path(__file__).with_name('normal_noise_decomposition.py'),
        gaussian=Path(__file__).with_name('reference_gaussian_regions.py'),
        weighting=Path(__file__).with_name('reference_weighted_regions.py'))
    for i in range(68):
        inputs[f'RNA{i}']=prior/f'context_{i:03d}.npz'
        if(contexts/f'Gaussian_context_{i:03d}.npz').exists():inputs[f'fit{i}']=contexts/f'Gaussian_context_{i:03d}.npz'
    rec=start_run(args.out,inputs,dict(mode=args.mode,primary=NAMES[1]),
        dict(scope='RNA-only fixed-fit boundary intervention; former evaluation regression-only'),mode='decomposed_boundary_predict')
    try:
        genes,bins,samples=[pd.read_csv(inputs[k],sep='\t')for k in ('genes','bins','samples')]
        effects,expected=[],[]
        for i in range(68):
            with np.load(inputs[f'RNA{i}'])as z:effects.append(z['corrected_gene']);expected.append(z['expected_gene'])
        effects,expected=np.asarray(effects),np.asarray(expected)
        normal=samples.role.eq('normal').to_numpy();assert len(samples)==68 and normal.sum()==20
        weighted,gaussian=load('reference_weighted_regions'),load('reference_gaussian_regions')
        training=gaussian.train(effects[normal],expected[normal],genes.bin.to_numpy(int))
        variance=load('normal_noise_decomposition').decompose(effects[normal],expected[normal],genes.bin.to_numpy(int),training)
        with np.load(inputs['decomposition'])as z:
            for key in z.files:np.testing.assert_array_equal(variance[key],z[key])
        with np.load(inputs['weights'])as z:
            spread=z['spread'];frozen_weights=z['tumour_variance']
            np.testing.assert_array_equal(z['sample_ids'],samples['sample'].to_numpy(str))
        reset=chromosome_starts(bins)|arm_starts(bins,genome_build='hg38')
        params=ProfileParams(**json.loads(inputs['native_parameters'].read_text())['params'])
        decisions={r['sample']:r for r in json.loads(inputs['decisions'].read_text())}
        names=OLD if args.mode=='control'else NAMES
        predictions={name:np.full((68,len(bins)),np.nan)for name in names};records=[]
        for k,i in enumerate(np.flatnonzero(~normal)):
            effect,_,_=weighted.normal_centered_effect(effects[i],effects[normal],genes.bin.to_numpy(int))
            weight=weighted.gene_variance(spread,expected[i],'reference')
            np.testing.assert_array_equal(weight,frozen_weights[k])
            fitted=dict(decisions[str(samples.iloc[i]['sample'])])
            with np.load(inputs[f'fit{i}'])as z:fitted.update({key:z[key]for key in z.files})
            reference=training['reference_variance']if args.mode=='control'else variance['intrinsic_variance']
            np.testing.assert_array_equal(fitted['total_variance'],reference+fitted['sampling_variance']+fitted['tau']**2)
            averaged,mass=load('posterior_region_average').average(effect,weight,genes.bin.to_numpy(int),reset,fitted,params)
            for j,name in enumerate(names):predictions[name][i]=averaged[j]
            records.append(dict(sample=str(samples.iloc[i]['sample']),beta=fitted['beta'],tau=fitted['tau'],mass=mass))
        np.savez_compressed(args.out/'RNA_predictions.npz',**predictions)
        write_json(dict(predictions=file_record(args.out/'RNA_predictions.npz'),cases=records,
            parameters_frozen=True,DNA_read=False,normal_variance_reconstruction_exact=True),args.out/'RNA_LOCKED.json')
        if args.mode=='control':
            with np.load(inputs['old_RNA'])as z:
                assert z.files==list(predictions)
                for name in z.files:np.testing.assert_array_equal(z[name],predictions[name])
        write_json(dict(mode=args.mode,queries=48,normal_variance_bins=len(reference),old_RNA_exact=args.mode=='control',
            DNA_read=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU exact posterior complete-run averages on frozen Gaussian fits')
    except BaseException as error:finish_run(args.out,rec,error=error);raise


def score(args):
    root=Path(os.environ['HAPLODOSE_DATA_ROOT']);base=load('profile')
    prior=root/'dev/programme_gene_profile_20261109/refit_000_v1'
    source=json.loads((prior/'run_manifest.json').read_text())['input']['files'];native=Path(source['genes.tsv.gz']['path']).parent
    previous=root/'dev/gene_laplace_profile_20261110/bulk_v1'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',samples=Path(source['samples']['path']),
        genes=native/'genes.tsv.gz',windows=previous/'window_scores.npz',old_scores=previous/'patient_scores.tsv',
        predictions=args.task/'candidate_v1/RNA_predictions.npz',RNA_lock=args.task/'candidate_v1/RNA_LOCKED.json',
        control=args.task/'control_v1/report.json',scorer=Path(base.__file__),transitive_scorer=Path(base.subarm.__file__),
        mapping=root/'dev/bulk_amplitude_reference_20261109/refit_000_v1/gene_map.tsv',
        split=root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv')
    assert json.loads(inputs['control'].read_text())['old_RNA_exact']
    rec=start_run(args.out,inputs,dict(primary=NAMES[1],comparators=42),
        dict(scope='Former evaluation regression-only; fixed coverage/threshold/scorer'),mode='decomposed_boundary_score')
    try:
        genes,samples=[pd.read_csv(inputs[k],sep='\t')for k in('genes','samples')]
        with np.load(inputs['windows'])as z:truth={k:z[k]for k in z.files}
        old=[k for k in truth if k not in('samples','common','loss','neutral','gain','dna','split')];assert len(old)==42
        windows={k:truth[k]for k in old};mapping=pd.read_csv(inputs['mapping'],sep='\t')
        with np.load(inputs['predictions'])as z:
            for name in NAMES:
                array=base.subarm.align_rows(z[name],samples['sample'],truth['samples'])
                mapped=base.subarm.gene_values(array,genes.set_index('gene').bin,mapping.gene.to_numpy())
                windows[name]=base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None],mapped,np.nan),mapping.window.to_numpy(int),truth['dna'].shape[0])
                assert np.isfinite(windows[name][truth['common']]).all()
        np.savez_compressed(args.out/'window_scores.npz',**truth,**{k:windows[k]for k in NAMES})
        split=pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any()and np.array_equal(split.split,truth['split'])
        calibration=np.flatnonzero(split.split.eq('calibration'));assert len(calibration)==24
        rows=[]
        for name,array in windows.items():
            threshold=base.loss_threshold(array,truth['neutral']&truth['common'],calibration)
            losses=base.score_loss(array,truth['loss'],truth['neutral'],truth['common'],threshold,negative=truth['gain']|truth['neutral'])
            for j,sample in enumerate(truth['samples']):
                x,y=array[truth['common'][:,j],j],truth['dna'][truth['common'][:,j],j]
                corr=float(pearsonr(x,y).statistic)if len(x)>=3 and np.isfinite(x).all()and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name in NAMES and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(corr):raise ArithmeticError('Candidate loses informative patient')
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=corr,**losses[j]))
        pd.DataFrame(rows).to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores=pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip')
        before=pd.read_csv(inputs['old_scores'],sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(before,scores[scores.method.isin(old)].reset_index(drop=True),check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        matrix=scores[scores.split.eq('evaluation')].pivot(index='sample',columns='method',values='pearson')
        contrasts=[]
        for name in old:
            difference=(matrix[NAMES[1]]-matrix[name]).dropna()
            contrasts.append(dict(baseline=name,n=len(difference),mean=float(difference.mean()),wins=int((difference>0).sum()),losses=int((difference<0).sum()),ties=int((difference==0).sum())))
        write_json(dict(methods=45,previous_methods_exact=42,primary=NAMES[1],means=matrix.mean().to_dict(),
            contrasts=contrasts,point_gate=all(c['mean']>0 for c in contrasts),independent_evaluation=False,
            production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU45 fixed bulk profiles, unchanged patient/threshold/DNA scoring')
    except BaseException as error:finish_run(args.out,rec,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('numerics','control','candidate','score'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.mode=='numerics':write_json(checks(),args.out/'report.json')
    elif args.mode=='score':score(args)
    else:predict(args)


if __name__=='__main__':main()
