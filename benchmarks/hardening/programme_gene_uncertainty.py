"""Exact-threshold nested uncertainty for the fixed bulk candidate family."""
import argparse
from collections import Counter
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom,pearsonr
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('programme_uncertainty_'+name,Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Threshold:
    """Cached empirical CDF with exact verification of the final quantile."""
    def __init__(self,scores,neutral,calibration,q):
        self.values = [np.sort(scores[neutral[:,j] & np.isfinite(scores[:,j]),j]) for j in calibration]
        self.n = np.array([len(v) for v in self.values],int)
        self.support = np.unique(np.concatenate([v for v in self.values if len(v)]))
        self.cdf = np.stack([np.searchsorted(v,self.support,side='right') for v in self.values])
        self.denominator = math.lcm(*(int(n) for n in self.n if n))
        self.weight = [self.denominator//int(n) if n else 0 for n in self.n]
        self.q = Fraction(str(q))

    def __call__(self,multiplicity):
        m = np.asarray(multiplicity,int)
        if m.shape != self.n.shape or np.any(m < 0):
            raise ValueError('Nonnegative calibration multiplicities required')
        total = int(m[self.n>0].sum())
        if not total:
            return np.nan
        mass = np.divide(m,self.n,out=np.zeros(len(m),float),where=self.n>0)@self.cdf
        k = min(int(np.searchsorted(mass,float(self.q)*total)),len(self.support)-1)
        target = self.q.numerator*self.denominator*total
        def enough(index):
            value = sum(int(m[j])*self.weight[j]*int(self.cdf[j,index]) for j in np.flatnonzero(m))
            return value*self.q.denominator >= target
        while k>0 and enough(k-1):
            k -= 1
        while not enough(k):
            k += 1
        if not enough(k) or (k>0 and enough(k-1)):
            raise ArithmeticError('Exact neutral threshold verification failed')
        return float(self.support[k])


def scores_at(values,positive,negative,threshold):
    result = []
    for mask in (positive,negative):
        result.append(np.array([float(np.searchsorted(np.sort(values[mask[:,j],j]),threshold,side='left')/mask[:,j].sum())
            if mask[:,j].any() else np.nan for j in range(values.shape[1])]))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    base,stats = load('profile'),load('statistics')
    methods = ['programme_gene_predictive','programme_gene_native_regions',
        'gene_predictive','naive_own_partition','hd_path_precision_mean']
    candidates,baselines = methods[:2],methods[2:]
    metrics = ('pearson','recall_loss','fpr_loss')
    family = [(a,b,e) for a in candidates for b in baselines for e in metrics]
    assert len(family) == 18
    inputs = dict(script=Path(__file__),design=args.task/'UNCERTAINTY_DESIGN.md',
        scorer=Path(base.__file__),statistics=Path(stats.__file__),
        transitive_scorer=Path(base.subarm.__file__),
        primary_predictions=args.task/'bulk_v1/RNA_predictions.npz',
        primary_scores=args.task/'bulk_v1/patient_scores.tsv')
    for i in range(9):
        inputs[f'windows{i}'] = args.task/f'refit_{i:03d}_v1/window_scores.npz'
        inputs[f'scores{i}'] = args.task/f'refit_{i:03d}_v1/patient_scores.tsv'
        inputs[f'predictions{i}'] = args.task/f'refit_{i:03d}_v1/RNA_predictions.npz'
    split_path = args.task.parents[1]/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv'
    inputs['split'] = split_path
    run = start_run(args.out,inputs,dict(draws=20000,seed=2026110914,family_size=18),
        dict(scope='Post-result finite-reference development characterization, not independent validation or coverage'),
        mode='programme_gene_uncertainty')
    try:
        with np.load(inputs['primary_predictions']) as old,np.load(inputs['predictions0']) as new:
            assert old.files == new.files
            for k in old.files:
                np.testing.assert_array_equal(old[k],new[k])
        pd.testing.assert_frame_equal(pd.read_csv(inputs['primary_scores'],sep='\t',float_precision='round_trip'),
            pd.read_csv(inputs['scores0'],sep='\t',float_precision='round_trip'),check_exact=True)
        split = pd.read_csv(split_path,sep='\t',dtype=str).set_index('sample')
        caches,threshold_checks,score_checks = [],0,0
        source_samples = None
        check_rng = np.random.default_rng(2026110915)
        for i in range(9):
            with np.load(inputs[f'windows{i}']) as z:
                d = {k:z[k] for k in z.files}
            if source_samples is None:
                source_samples = d['samples'].copy()
            np.testing.assert_array_equal(d['samples'],source_samples)
            meta = split.loc[d['samples']]
            cal,ev = np.flatnonzero(meta.split.eq('calibration')),np.flatnonzero(meta.split.eq('evaluation'))
            assert len(cal)==len(ev)==24 and not meta.patient.duplicated().any()
            methods_cache = {}
            for method in methods:
                values = d[method]
                threshold = Threshold(values,d['neutral']&d['common'],cal,base.FPR)
                allocations = [np.ones(len(cal),int)]+[np.bincount(check_rng.integers(len(cal),size=len(cal)),minlength=len(cal)) for _ in range(40)]
                for m in allocations:
                    columns = np.repeat(cal,m)
                    actual = threshold(m)
                    expected = base.loss_threshold(values,d['neutral']&d['common'],columns)
                    assert actual == expected or (np.isnan(actual) and np.isnan(expected))
                    threshold_checks += 1
                t = threshold(np.ones(len(cal),int))
                recall,fpr = scores_at(values,d['loss']&d['common'],d['neutral']&d['common'],t)
                ordinary = base.score_loss(values,d['loss'],d['neutral'],d['common'],t,negative=d['gain']|d['neutral'])
                np.testing.assert_array_equal(recall,np.array([r['recall_loss'] for r in ordinary]))
                np.testing.assert_array_equal(fpr,np.array([r['fpr_loss'] for r in ordinary]))
                score_checks += len(ordinary)*2
                cor = np.array([float(pearsonr(values[d['common'][:,j],j],d['dna'][d['common'][:,j],j]).statistic)
                    if d['common'][:,j].sum()>=3 and np.ptp(values[d['common'][:,j],j])>0 and np.ptp(d['dna'][d['common'][:,j],j])>0
                    else np.nan for j in range(len(d['samples']))])
                positive = [np.sort(values[d['loss'][:,j]&d['common'][:,j],j]) for j in range(len(d['samples']))]
                negative = [np.sort(values[d['neutral'][:,j]&d['common'][:,j],j]) for j in range(len(d['samples']))]
                methods_cache[method] = dict(threshold=threshold,pearson=cor,positive=positive,negative=negative,
                    point=dict(pearson=cor,recall_loss=recall,fpr_loss=fpr))
            caches.append(methods_cache)
        def mean_difference(a,b,e,score):
            sign = -1 if e=='fpr_loss' else 1
            delta = sign*(score[a][e][ev]-score[b][e][ev])
            valid = delta[np.isfinite(delta)]
            return float(valid.mean()) if len(valid) else np.nan
        point = {m:caches[0][m]['point'] for m in methods}
        estimates = {key:mean_difference(*key,point) for key in family}
        draws = np.empty((20000,len(family)))
        rng = np.random.default_rng(2026110914)
        for b in range(20000):
            refit = int(rng.integers(1,9))
            m = np.bincount(rng.integers(len(cal),size=len(cal)),minlength=len(cal))
            ev = rng.choice(np.flatnonzero(meta.split.eq('evaluation')),len(ev),replace=True)
            scored = {}
            for method in methods:
                cache = caches[refit][method]
                t = cache['threshold'](m)
                if np.isfinite(t):
                    recall,fpr = [np.array([float(np.searchsorted(v,t,side='left')/len(v)) if len(v) else np.nan for v in cache[k]])
                        for k in ('positive','negative')]
                else:
                    recall,fpr = [np.full(len(cache['pearson']),np.nan) for _ in range(2)]
                scored[method] = dict(pearson=cache['pearson'],recall_loss=recall,fpr_loss=fpr)
            for j,key in enumerate(family):
                draws[b,j] = mean_difference(*key,scored)
            if b<40:
                for j,(a,c,e) in enumerate(family):
                    sign = -1 if e=='fpr_loss' else 1
                    delta = sign*(scored[a][e][ev]-scored[c][e][ev])
                    finite = [float(v) for v in delta if np.isfinite(v)]
                    if finite:
                        expected = math.fsum(finite)/len(finite)
                        operations = 2*len(finite)+4
                        allowance = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,max(map(abs,finite)))
                        assert abs(draws[b,j]-expected) <= allowance
            if (b+1)%1000==0:
                print(json.dumps(dict(completed=b+1,total=20000)),flush=True)
        np.save(args.out/'bootstrap_improvement.npy',draws)
        records = []
        for j,(a,b,e) in enumerate(family):
            boot = draws[:,j]
            boot = boot[np.isfinite(boot)]
            delta = estimates[(a,b,e)]
            testable = bool(np.isfinite(delta) and len(boot)>=18000 and np.ptp(boot)>0)
            record = dict(a=a,b=b,endpoint=e,mean_improvement=delta,
                raw_mean_difference=-delta if e=='fpr_loss' else delta,
                family_size=18,testable=testable,p=1.,family_ci_low=np.nan,family_ci_high=np.nan,
                valid_bootstrap_draws=len(boot),independent_evaluation=False,
                uncertainty='finite normal-reference, calibration patient and comparison patient resampling')
            if testable:
                error = boot-delta
                record['p'] = float((1+np.sum(abs(error)>=abs(delta)))/(len(error)+1))
                record['family_ci_low'],record['family_ci_high'] = map(float,delta-np.quantile(error,[1-.05/36,.05/36]))
                endpoint_brackets = []
                for probability in (.05/36,1-.05/36):
                    ordered = np.sort(boot)
                    lo = max(0,int(binom.ppf(.05/72,len(boot),probability))-1)
                    hi = min(len(boot)-1,int(binom.ppf(1-.05/72,len(boot),probability)))
                    endpoint_brackets.append([float(ordered[lo]),float(ordered[hi])])
                record['Monte_Carlo_bootstrap_quantile_brackets'] = endpoint_brackets
            records.append(record)
        records = stats.adjust_family(records)
        write_json(dict(contrasts=records,reference_draws=8,biological_comparison_patients=24,
            threshold_exact_checks=threshold_checks,score_lookup_checks=score_checks,
            primary_prediction_and_score_replay_exact=True,family=18,draws=20000,
            biological_calibration=False,independent_evaluation=False,default_promoted=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU exact-threshold nested biological-unit bootstrap')
    except BaseException as error:
        finish_run(args.out,run,error=error)
        raise


if __name__ == '__main__':
    main()
