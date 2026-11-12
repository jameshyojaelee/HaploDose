"""Fixed reference-bootstrap aggregation with identical treatment of baselines."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


PRIMARY = 'reference_bag_averaged_gaussian_reference_median'
METADATA = ('samples','common','loss','neutral','gain','dna','split')


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('reference_bag_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mean(values):
    values = np.asarray(values,float)
    if values.shape[0]!=8:
        raise ValueError('All eight frozen resamples required')
    finite = np.isfinite(values)
    if not np.all(finite==finite[:1]):
        raise ValueError('Reference replicates must have identical finite support')
    result = np.full(values.shape[1:],np.nan)
    for position in zip(*np.where(finite[0])):
        result[position] = math.fsum(float(values[(j,*position)]) for j in range(8))/8
    return result


def aggregate(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inputs = dict(script=Path(__file__),design=args.task/'DESIGN.md')
    for i in range(1,9):
        inputs[f'refit{i}'] = root/f'dev/posterior_region_average_20261109/refit_{i:03d}_v1/window_scores.npz'
    rec = start_run(args.out,inputs,dict(primary=PRIMARY,refits=list(range(1,9)),operator='equal arithmetic mean'),
        dict(scope='RNA method arrays only; archives also contain frozen DNA arrays, not accessed'),
        mode='reference_bootstrap_profile_aggregation')
    try:
        toy = np.arange(16,dtype=float).reshape(8,2,1)/8
        expected = np.array([float(sum(Fraction(float(toy[i,j,0]))for i in range(8))/8)for j in range(2)])
        np.testing.assert_array_equal(mean(toy)[:,0],expected)
        np.testing.assert_array_equal(mean(toy[::-1]),mean(toy))
        bags,identifiers = {},None
        for i in range(1,9):
            with np.load(inputs[f'refit{i}']) as z:
                names = [k for k in z.files if k not in METADATA]
                if identifiers is None:
                    assert len(names)==39
                    identifiers = z['samples']
                    bags = {name:[] for name in names}
                assert names==list(bags)
                np.testing.assert_array_equal(identifiers,z['samples'])
                for name in names:bags[name].append(z[name])
        output = {'reference_bag_'+name:mean(value)for name,value in bags.items()}
        assert PRIMARY in output
        np.savez_compressed(args.out/'RNA_windows.npz',samples=identifiers,**output)
        write_json(dict(artifact=file_record(args.out/'RNA_windows.npz'),methods=39,
            primary=PRIMARY,DNA_arrays_accessed=False,resamples=list(range(1,9)),
            exact_fraction_and_permutation_checks=True,threshold_fitted=False),args.out/'RNA_LOCKED.json')
        finish_run(args.out,rec,backend='CPU compensated fixed eight-reference RNA mean')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


def score(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    previous = root/'dev/decomposed_boundary_average_20261110/score_v1'
    inputs = dict(script=Path(__file__),design=args.task/'DESIGN.md',
        windows=previous/'window_scores.npz',old_scores=previous/'patient_scores.tsv',
        RNA=args.task/'aggregate_v1/RNA_windows.npz',RNA_lock=args.task/'aggregate_v1/RNA_LOCKED.json',
        split=root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv',
        scorer=Path(__file__).with_name('profile.py'))
    base = load('profile')
    inputs['transitive_scorer'] = Path(base.subarm.__file__)
    rec = start_run(args.out,inputs,dict(primary=PRIMARY,old_methods=45,bagged_methods=39),
        dict(scope='Former evaluation regression-only; all methods retained'),mode='reference_bag_score')
    try:
        assert json.loads(inputs['RNA_lock'].read_text())['DNA_arrays_accessed']is False
        with np.load(inputs['windows'])as z:truth = {k:z[k]for k in z.files}
        old = [k for k in truth if k not in METADATA]
        assert len(old)==45
        methods = {name:truth[name]for name in old}
        with np.load(inputs['RNA'])as z:
            np.testing.assert_array_equal(z['samples'],truth['samples'])
            bags = [name for name in z.files if name!='samples']
            assert len(bags)==39
            for name in bags:
                methods[name] = z[name]
                assert np.isfinite(methods[name][truth['common']]).all()
        split = pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any()and np.array_equal(split.split,truth['split'])
        calibration = np.flatnonzero(split.split.eq('calibration'))
        assert len(calibration)==24
        rows = []
        for name,array in methods.items():
            threshold = base.loss_threshold(array,truth['neutral']&truth['common'],calibration)
            losses = base.score_loss(array,truth['loss'],truth['neutral'],truth['common'],threshold,
                negative=truth['gain']|truth['neutral'])
            for j,sample in enumerate(truth['samples']):
                x,y = array[truth['common'][:,j],j],truth['dna'][truth['common'][:,j],j]
                corr = float(pearsonr(x,y).statistic)if len(x)>=3 and np.isfinite(x).all()and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name in bags and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(corr):
                    raise ArithmeticError('Aggregated method loses informative patient')
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=corr,**losses[j]))
        pd.DataFrame(rows).to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores = pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip')
        before = pd.read_csv(inputs['old_scores'],sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(before,scores[scores.method.isin(old)].reset_index(drop=True),check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        matrix = scores[scores.split.eq('evaluation')].pivot(index='sample',columns='method',values='pearson')
        contrasts = []
        for name in methods:
            if name==PRIMARY:continue
            difference = (matrix[PRIMARY]-matrix[name]).dropna()
            contrasts.append(dict(baseline=name,n=len(difference),mean=float(difference.mean()),
                wins=int((difference>0).sum()),losses=int((difference<0).sum()),ties=int((difference==0).sum())))
        np.savez_compressed(args.out/'window_scores.npz',**truth,**{name:methods[name]for name in bags})
        write_json(dict(methods=84,old_methods_exact=45,primary=PRIMARY,means=matrix.mean().to_dict(),
            contrasts=contrasts,point_gate=all(c['mean']>0 for c in contrasts),
            independent_evaluation=False,biological_replication_added=0,production_changed=False,
            scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU fixed84-method patient-level paired bulk scoring')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('aggregate','score'),required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.mode=='aggregate':aggregate(args)
    else:score(args)


if __name__=='__main__':main()
