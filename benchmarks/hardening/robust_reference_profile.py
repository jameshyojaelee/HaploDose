"""Fixed robust reference aggregation to test attenuation of bulk loss signal."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


PRIMARY='reference_median_averaged_gaussian_reference_median'
METADATA=('samples','common','loss','neutral','gain','dna','split')


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('robust_reference_'+name,
        Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def median(values):
    values=np.asarray(values,float)
    if values.ndim!=3 or values.shape[0]!=8:
        raise ValueError('All eight fixed window-by-specimen reference fits required')
    support=np.isfinite(values)
    if not np.all(support==support[:1]):
        raise ValueError('Reference fits must have identical finite support')
    output=np.full(values.shape[1:],np.nan)
    # Same even-sample median for every method. Half first prevents overflow;
    # no trimming fraction, fit weight or DNA-selected reference is introduced.
    ordered=np.sort(values[:,support[0]],axis=0)
    output[support[0]]=ordered[3]/2+ordered[4]/2
    return output


def aggregate(args):
    root=Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md')
    for i in range(1,9):
        inputs[f'refit{i}']=root/f'dev/posterior_region_average_20261109/refit_{i:03d}_v1/window_scores.npz'
    rec=start_run(args.out,inputs,dict(primary=PRIMARY,refits=list(range(1,9)),operator='fixed even-sample median'),
        dict(scope='RNA arrays only; archived DNA arrays not accessed'),mode='robust_reference_aggregation')
    try:
        toy=np.array([-10.,-2.,-.5,.25,.75,3.,5.,100.]).reshape(8,1,1)
        exact=float((Fraction(float(toy[3,0,0]))+Fraction(float(toy[4,0,0])))/2)
        assert median(toy)[0,0]==exact
        np.testing.assert_array_equal(median(toy[::-1]),median(toy))
        # A minority reference outlier cannot attenuate a repeatedly supported
        # effect; this property is not a guarantee of biological validity.
        repeated=np.array([-1.]*7+[0.]).reshape(8,1,1)
        assert median(repeated)[0,0]==-1.
        bags,identifiers={},None
        for i in range(1,9):
            with np.load(inputs[f'refit{i}'])as z:
                names=[k for k in z.files if k not in METADATA]
                if identifiers is None:
                    assert len(names)==39
                    identifiers=z['samples'];bags={name:[]for name in names}
                assert names==list(bags)
                np.testing.assert_array_equal(identifiers,z['samples'])
                for name in names:bags[name].append(z[name])
        output={'reference_median_'+name:median(value)for name,value in bags.items()}
        assert PRIMARY in output
        np.savez_compressed(args.out/'RNA_windows.npz',samples=identifiers,**output)
        write_json(dict(artifact=file_record(args.out/'RNA_windows.npz'),primary=PRIMARY,
            methods=39,DNA_arrays_accessed=False,resamples=list(range(1,9)),
            exact_rational_fixture=True,permutation_exact=True,minority_outlier_fixture=True,
            threshold_fitted=False),args.out/'RNA_LOCKED.json')
        finish_run(args.out,rec,backend='CPU fixed eight-reference RNA median')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


def score(args):
    root=Path(os.environ['HAPLODOSE_DATA_ROOT'])
    previous=root/'dev/reference_bagged_profile_20261110/score_v1'
    base=load('profile')
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',
        windows=previous/'window_scores.npz',old_scores=previous/'patient_scores.tsv',
        RNA=args.task/'aggregate_v1/RNA_windows.npz',RNA_lock=args.task/'aggregate_v1/RNA_LOCKED.json',
        split=root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv',
        scorer=Path(base.__file__),transitive_scorer=Path(base.subarm.__file__))
    rec=start_run(args.out,inputs,dict(primary=PRIMARY,old_methods=84,added_methods=39),
        dict(scope='Exposed patient-level regression; no independent evaluation'),mode='robust_reference_score')
    try:
        assert json.loads(inputs['RNA_lock'].read_text())['DNA_arrays_accessed']is False
        with np.load(inputs['windows'])as z:truth={k:z[k]for k in z.files}
        old=[k for k in truth if k not in METADATA];assert len(old)==84
        methods={name:truth[name]for name in old}
        with np.load(inputs['RNA'])as z:
            np.testing.assert_array_equal(z['samples'],truth['samples'])
            added=[name for name in z.files if name!='samples'];assert len(added)==39
            for name in added:
                assert name not in methods
                methods[name]=z[name]
                assert np.isfinite(methods[name][truth['common']]).all()
        split=pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any()and np.array_equal(split.split,truth['split'])
        calibration=np.flatnonzero(split.split.eq('calibration'));assert len(calibration)==24
        rows=[]
        for name,array in methods.items():
            threshold=base.loss_threshold(array,truth['neutral']&truth['common'],calibration)
            losses=base.score_loss(array,truth['loss'],truth['neutral'],truth['common'],threshold,
                negative=truth['gain']|truth['neutral'])
            for j,sample in enumerate(truth['samples']):
                x,y=array[truth['common'][:,j],j],truth['dna'][truth['common'][:,j],j]
                corr=float(pearsonr(x,y).statistic)if len(x)>=3 and np.isfinite(x).all()and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name in added and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(corr):
                    raise ArithmeticError('Candidate loses informative patient')
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=corr,**losses[j]))
        scores=pd.DataFrame(rows);scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        replay=pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip')
        before=pd.read_csv(inputs['old_scores'],sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(before,replay[replay.method.isin(old)].reset_index(drop=True),check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        evaluation=scores[scores.split.eq('evaluation')]
        matrix=evaluation.pivot(index='sample',columns='method',values='pearson')
        contrasts=[]
        for name in methods:
            if name==PRIMARY:continue
            d=(matrix[PRIMARY]-matrix[name]).dropna()
            contrasts.append(dict(baseline=name,n=len(d),mean=float(d.mean()),wins=int((d>0).sum()),losses=int((d<0).sum())))
        comparator='averaged_gaussian_reference_median'
        recall=evaluation.pivot(index='sample',columns='method',values='recall_loss')
        fpr=evaluation.pivot(index='sample',columns='method',values='fpr_loss')
        recall_effect=float((recall[PRIMARY]-recall[comparator]).dropna().mean())
        fpr_effect=float((fpr[PRIMARY]-fpr[comparator]).dropna().mean())
        point=all(c['mean']>0 for c in contrasts)
        joint=point and recall_effect>=0 and fpr_effect<=0
        np.savez_compressed(args.out/'window_scores.npz',**truth,**{k:methods[k]for k in added})
        write_json(dict(primary=PRIMARY,methods=123,old_methods_exact=84,means=matrix.mean().to_dict(),
            contrasts=contrasts,point_gate=point,loss_comparator=comparator,
            paired_recall_effect=recall_effect,paired_fpr_effect=fpr_effect,
            joint_point_gate=joint,independent_evaluation=False,biological_replication_added=0,
            production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU fixed123-method paired bulk scoring')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--mode',choices=('aggregate','score'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    (aggregate if args.mode=='aggregate'else score)(args)


if __name__=='__main__':main()
