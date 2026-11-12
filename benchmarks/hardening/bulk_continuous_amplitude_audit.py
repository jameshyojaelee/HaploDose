"""Scalar count, Decimal amplitude, window and patient-score audit."""
import argparse
from decimal import Decimal,localcontext
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    spec=importlib.util.spec_from_file_location('continuous_audit_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--panel-root',default='panel_v1')
    parser.add_argument('--reference-refit',type=int,choices=range(9))
    a=parser.parse_args();a.out.mkdir(exist_ok=False);panel=a.task/a.panel_root;root=a.task.parents[1]
    hd=root/'dev/bulk_posterior_export_20261109/stage_v1/haplodose'
    previous=root/'dev/bulk_geometry_20261109/panel_v1'
    if a.reference_refit is not None:
        hd=root/f'dev/hardening_20261106_2119/profile_v6/raw_refit_{a.reference_refit:03d}/haplodose'
        previous=root/'dev/bulk_continuous_amplitude_20261109/panel_v2'
    quantile=load('bulk_posterior_stage_audit');correlation=load('pearson_reference');candidate=load('bulk_continuous_amplitude')
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',manifest=panel/'run_manifest.json',
                windows=panel/'window_scores.npz',report=panel/'report.json',scores=panel/'patient_scores.tsv',
                segments=panel/'segments.tsv',mapping=panel/'gene_map.tsv')
    manifest=json.loads(inputs['manifest'].read_text());assert manifest['status']=='complete'
    for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    rec=start_run(a.out,inputs,{},dict(scope='Independent numerical and scoring audit,not superiority inference'),mode='bulk_continuous_amplitude_audit')
    try:
        correlation.verify();u=np.finfo(float).eps
        fixture=np.array([0.,1.,3.]);variance=np.array([1.,2.,4.]);path=np.zeros(3,int);reset=np.array([1,0,0],bool)
        mean,precision,_=candidate.amplitudes(fixture,variance,path,reset)
        budget=32*u*3*(1+max(abs(fixture)))
        assert max(abs(mean-4/3))<=budget and max(abs(precision-5/7))<=budget
        shifted=candidate.amplitudes(fixture+.7,variance*17,path,reset)
        assert max(abs(shifted[0]-mean-.7))<=budget and max(abs(shifted[1]-precision-.7))<=budget
        report=json.loads(inputs['report'].read_text())
        if a.reference_refit is not None:assert report['reference_refit']==a.reference_refit
        with np.load(inputs['windows']) as z:w={k:z[k] for k in z.files}
        with np.load(previous/'window_scores.npz') as z:old={k:z[k] for k in z.files}
        for key in ('samples','common','loss','neutral','gain','dna','split'):assert np.array_equal(w[key],old[key],equal_nan=True) if w[key].dtype.kind not in 'US' else np.array_equal(w[key],old[key])
        old_methods=json.loads((previous/'report.json').read_text())['methods']
        if a.reference_refit is not None:assert set(report['methods'])==set(old_methods)
        if a.reference_refit is None or a.reference_refit==0:
            for method in old_methods:assert np.array_equal(w[method],old[method],equal_nan=True)
        genes=pd.read_csv(hd/'genes.tsv.gz',sep='\t');bins=pd.read_csv(hd/'bins.tsv',sep='\t')
        refit=a.reference_refit if a.reference_refit is not None else 0
        raw=pd.read_csv(root/f'dev/hardening_20261106_2119/profile_v6/raw_refit_{refit:03d}/counts.tsv',sep='\t',index_col=0)
        donor_checks=0
        if a.reference_refit is not None:
            original=root/'dev/hardening_20261106_2119/profile_v6/raw_refit_000'
            original_counts=pd.read_csv(original/'counts.tsv',sep='\t',index_col=0)
            original_meta=pd.read_csv(original/'samples.tsv',sep='\t').set_index('sample')
            input_meta=pd.read_csv(hd.parent/'samples.tsv',sep='\t')
            assert raw.index.is_unique and np.array_equal(raw.index,original_counts.index)
            assert input_meta['sample'].is_unique
            for row in input_meta.itertuples():
                source_sample=getattr(row,'source_sample',np.nan)
                source=source_sample if row.role=='normal' and pd.notna(source_sample) else row.sample
                assert original_meta.at[source,'role']==row.role
                assert original_meta.at[source,'patient']==row.patient
                assert np.array_equal(raw[row.sample].to_numpy(),original_counts[source].to_numpy())
                donor_checks+=int(row.role=='normal')
            assert donor_checks==20
        mapping=pd.read_csv(inputs['mapping'],sep='\t');segments=pd.read_csv(inputs['segments'],sep='\t',float_precision='round_trip')
        expected_reset=np.r_[True,(bins.chromosome.to_numpy()[1:]!=bins.chromosome.to_numpy()[:-1])|
                             (bins.genomic_block.to_numpy()[1:]!=bins.genomic_block.to_numpy()[:-1])]
        reference={name:[] for name in ('hd_path_mean','hd_path_precision_mean')};numerical=[];count_checks=0;decimal_checks=0
        for sample in w['samples']:
            z=np.load(panel/(str(sample)+'.npz'));assert np.array_equal(z['reset'],expected_reset)
            expected=np.array([sum(int(raw.at[g,str(sample)]) for g in genes.loc[genes.bin.eq(j),'gene']) for j in range(len(bins))])
            assert np.array_equal(z['Y'],expected);count_checks+=len(expected)
            Y,E=z['Y'].astype(np.longdouble),z['E'].astype(np.longdouble)
            track=np.log((Y+.5)/(E+.5))-np.longdouble(z['centre']);var=z['sd'].astype(np.longdouble)**2+1/(Y+.5)
            track_guard=64*u*(1+float(np.max(abs(track))))
            assert np.max(abs(track-z['track']))<=track_guard
            assert np.max(abs(var-z['variance']))<=64*u*(1+float(np.max(var)))
            path=z['path'];cuts=[0]+[j for j in range(1,len(path)) if expected_reset[j] or path[j]!=path[j-1]]+[len(path)]
            predicted={key:np.empty(len(path)) for key in reference}
            selected=segments[segments['sample'].eq(sample)]
            assert list(zip(selected.first_bin,selected.last_bin))==[(lo,hi-1) for lo,hi in zip(cuts[:-1],cuts[1:])]
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                r=np.asarray(track[lo:hi],float);v=np.asarray(var[lo:hi],float)
                first=math.fsum(r)/(hi-lo)
                left,right=min(r),max(r)
                second=left if left==right else brentq(lambda b:math.fsum((r-b)/v),left,right,
                    xtol=np.spacing(max(1.,abs(left),abs(right))),rtol=4*u)
                guard=64*u*(hi-lo+1)*(1+max(abs(r)))
                assert np.max(abs(z['mean'][lo:hi]-first))<=guard
                assert np.max(abs(z['precision_mean'][lo:hi]-second))<=guard
                if lo==0 or hi==len(path):
                    with localcontext() as ctx:
                        ctx.prec=80
                        values=[(Decimal.from_float(float(z['Y'][j]))+Decimal('.5'))/(Decimal.from_float(float(z['E'][j]))+Decimal('.5')) for j in range(lo,hi)]
                        dr=[value.ln()-Decimal.from_float(float(z['centre'])) for value in values]
                        dv=[Decimal.from_float(float(z['sd'][j]))**2+1/(Decimal.from_float(float(z['Y'][j]))+Decimal('.5')) for j in range(lo,hi)]
                        exact_mean=sum(dr)/len(dr);exact_weighted=sum(x/s for x,s in zip(dr,dv))/sum(1/s for s in dv)
                        assert abs(float(exact_mean)-first)<=guard and abs(float(exact_weighted)-second)<=guard
                        decimal_checks+=1
                predicted['hd_path_mean'][lo:hi]=first;predicted['hd_path_precision_mean'][lo:hi]=second
                numerical.append(dict(sample=str(sample),first_bin=lo,last_bin=hi-1,guard=guard,
                    mean_difference=float(max(abs(z['mean'][lo:hi]-first))),precision_difference=float(max(abs(z['precision_mean'][lo:hi]-second)))))
            for method,values in predicted.items():
                output=[]
                for window in range(w[method].shape[0]):
                    loc=mapping[(mapping.window==window)&mapping.shared&mapping.bin.notna()]
                    ordered=sorted(float(values[int(j)]) for j in loc.bin)
                    n=len(ordered);output.append(ordered[n//2] if n%2 else (ordered[n//2-1]+ordered[n//2])/2 if n else np.nan)
                reference[method].append(output)
        scores=pd.read_csv(inputs['scores'],sep='\t',float_precision='round_trip').set_index(['sample','method']);assert scores.index.is_unique
        ci=np.flatnonzero(w['split']=='calibration');assert len(ci)==24
        decisions=[];scoring_checks=0
        for method in report['methods']:
            values=w[method];threshold=quantile.threshold_reference(values,w['neutral']&w['common'],ci);assert threshold==report['thresholds'][method]
            for j,sample in enumerate(w['samples']):
                row=scores.loc[(sample,method)]
                for field,mask in (('recall_loss',w['loss']),('fpr_loss',w['neutral'])):
                    support=mask[:,j]&w['common'][:,j];expected=sum(float(x)<threshold for x in values[support,j])/int(sum(support)) if support.any() else np.nan
                    assert (pd.isna(row[field]) and np.isnan(expected)) or row[field]==expected;scoring_checks+=1
                selected=values[:,j]<threshold;eligible=w['common'][:,j];positives=w['loss'][:,j];negatives=w['gain'][:,j]|w['neutral'][:,j]
                calls=int(sum(selected&eligible&(positives|negatives)));true=int(sum(selected&eligible&positives))
                precision=true/calls if calls else np.nan
                assert (pd.isna(row.precision_loss) and np.isnan(precision)) or row.precision_loss==precision;scoring_checks+=1
                x,y=values[eligible,j],w['dna'][eligible,j]
                expected=correlation.reference(x,y)['value'] if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                assert (pd.isna(row.pearson) and np.isnan(expected)) or row.pearson==expected;scoring_checks+=1
            if method in reference:
                refined=np.array(reference[method]).T
                difference=np.abs(refined-values);valid=np.isfinite(difference)
                bound=max(r['guard'] for r in numerical)
                # The unchanged window scorer stores each median in float32.
                # Keep amplitude checks in binary64 and reproduce that declared
                # storage operation exactly rather than widening their guard.
                assert np.array_equal(refined.astype(np.float32),values,equal_nan=True)
                half_ulp=np.spacing(np.abs(values).astype(np.float32)).astype(float)/2
                assert np.all(difference[valid]<=bound+half_ulp[valid])
                new_threshold=quantile.threshold_reference(refined,w['neutral']&w['common'],ci)
                changed=((values<threshold)!=(refined<new_threshold))&w['common']
                decisions.append(dict(method=method,maximum_window_discrepancy=float(np.max(difference[valid])),
                    float32_window_storage_reproduces_exactly=True,amplitude_guard=bound,
                    maximum_half_ulp=float(np.max(half_ulp[valid])),
                    threshold=threshold,refined_threshold=new_threshold,callable_decision_disagreements=int(changed.sum())))
        paired=[]
        table=scores.reset_index();comparison=table[table.split.eq('evaluation')]
        for method in reference:
            controls=[m for m in report['methods'] if m not in reference] if a.reference_refit is not None else ('hd_posterior','naive_own_partition')
            for control in controls:
                one=comparison[comparison.method.eq(method)].set_index('sample');two=comparison[comparison.method.eq(control)].set_index('sample')
                for endpoint in ('pearson','recall_loss','fpr_loss','precision_loss'):
                    valid=one[endpoint].notna()&two[endpoint].notna();ids=one.index[valid]
                    effect=math.fsum(float(one.at[s,endpoint]-two.at[s,endpoint]) for s in ids)/len(ids) if len(ids) else None
                    paired.append(dict(method=method,control=control,endpoint=endpoint,paired_patients=len(ids),total_comparison_patients=24,
                                       mean_paired_effect=effect,patients=ids.tolist(),uncertainty='Not estimated;descriptive exposed-panel mechanism'))
        write_json(numerical,a.out/'amplitude_numerics.json');write_json(paired,a.out/'paired_effects.json');write_json(decisions,a.out/'scientific_boundaries.json')
        write_json(dict(all_checks_passed=True,integer_bin_sums=count_checks,segments=len(numerical),decimal80_checks=decimal_checks,
                        independent_patient_endpoints=scoring_checks,reference_decisions=decisions,
                        normal_source_correspondence_checks=donor_checks,method_promoted=False,new_hypothesis_tests=0),a.out/'report.json')
        finish_run(a.out,rec,backend='CPU scalar counts/roots,Decimal80,Fraction calibration and patient scoring')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
