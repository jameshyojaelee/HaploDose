"""Independent geometry, scientific-boundary and paired-population audit."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.provenance import start_run, finish_run, file_record


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False);panel=args.task/'panel_v1'
    hd=args.task.parents[1]/'dev/bulk_posterior_export_20261109/stage_v1/haplodose'
    path=Path(__file__).with_name('bulk_posterior_stage_audit.py')
    spec=importlib.util.spec_from_file_location('geometry_exact_quantile_reference',path)
    quantile=importlib.util.module_from_spec(spec);spec.loader.exec_module(quantile)
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',quantile_reference=path,
        manifest=panel/'run_manifest.json',geometry=panel/'geometry.npz',blocks=panel/'blocks.json',
        windows=panel/'window_scores.npz',scores=panel/'patient_scores.tsv',report=panel/'report.json',
        bins=hd/'bins.tsv',genes=hd/'genes.tsv.gz')
    manifest=json.loads(inputs['manifest'].read_text());assert manifest['status']=='complete'
    for record in manifest['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(args.out,inputs,{},dict(scope='Independent exposed-panel geometry audit; no promotion or new hypothesis tests'),mode='bulk_geometry_audit')
    try:
        with np.load(inputs['geometry']) as z:g={k:z[k] for k in z.files}
        with np.load(inputs['windows']) as z:w={k:z[k] for k in z.files}
        bins=pd.read_csv(inputs['bins'],sep='\t');genes=pd.read_csv(inputs['genes'],sep='\t')
        anchors=[]
        for j in range(len(bins)):
            values=sorted(int(x) for x in genes.loc[genes.bin.eq(j),'start']);n=len(values)
            assert n;anchors.append(values[n//2] if n%2 else (values[n//2-1]+values[n//2])/2)
        assert np.array_equal(anchors,g['anchors'])
        reset=np.r_[True,(bins.chromosome.to_numpy()[1:]!=bins.chromosome.to_numpy()[:-1])|
                    (bins.genomic_block.to_numpy()[1:]!=bins.genomic_block.to_numpy()[:-1])]
        assert np.array_equal(reset,g['reset'])
        q=np.array([.2,.5,.2,.1],np.float32).astype(float);q/=q.sum()
        stationary=q*(1-q);stationary/=stationary.sum();root=np.sqrt(stationary)
        value,vectors=np.linalg.eigh(root[:,None]*g['transition']/root[None]);assert (value>0).all()
        def spectral(time):return np.eye(4) if time==0 else ((vectors*value**time)@vectors.T)*root[None]/root[:,None]
        endpoints=[];u=np.finfo(float).eps
        for index in np.split(np.arange(len(bins)),np.flatnonzero(reset)[1:]):
            edges=index[1:];steps=len(edges)
            assert np.array_equal(g['distance_bp'][edges],np.diff(np.array(anchors)[index]))
            assert set(g['permutation'][edges])==set(edges)
            assert np.array_equal(g['shuffled_ratios'][edges],g['ratios'][g['permutation'][edges]])
            assert abs(math.fsum(g['ratios'][edges])-steps) <= (2*steps+3)*u*max(1,steps)
            product=np.eye(4);permuted=np.eye(4)
            for j in edges:
                product=product@spectral(g['ratios'][j]);permuted=permuted@spectral(g['shuffled_ratios'][j])
            target=np.linalg.matrix_power(g['transition'],steps)
            operations=4**3*max(1,steps)*4;budget=operations*u/(1-operations*u)
            error=float(max(np.max(np.abs(product-target)),np.max(np.abs(permuted-target))))
            assert error<=budget
            endpoints.append(dict(first_bin=int(index[0]),edges=steps,endpoint_discrepancy=error,binary64_budget=budget))
        scores=pd.read_csv(inputs['scores'],sep='\t',float_precision='round_trip')
        assert not scores[['sample','method']].duplicated().any()
        ci=np.flatnonzero(w['split']=='calibration');ei=np.flatnonzero(w['split']=='evaluation')
        assert len(ci)==len(ei)==24
        report=json.loads(inputs['report'].read_text());boundary=[];counts=[]
        for method in report['methods']:
            values=w[method];threshold=quantile.threshold_reference(values,w['neutral']&w['common'],ci)
            assert threshold==report['thresholds'][method]
            selected=scores[scores.method.eq(method)].set_index('sample').loc[w['samples']]
            assert len(selected)==48
            for j in range(48):
                for field,mask in (('recall_loss',w['loss']),('fpr_loss',w['neutral'])):
                    eligible=mask[:,j]&w['common'][:,j]
                    expected=sum(float(x)<threshold for x in values[eligible,j])/int(eligible.sum()) if eligible.any() else np.nan
                    actual=selected.iloc[j][field]
                    assert (pd.isna(actual) and np.isnan(expected)) or actual==expected
                hit=values[:,j]<threshold;eligible=w['common'][:,j];loss=w['loss'][:,j]
                true=int((hit&eligible&loss).sum());negative=w['gain'][:,j]|w['neutral'][:,j]
                calls=int((hit&eligible&(loss|negative)).sum())
                expected=true/calls if calls else np.nan;actual=selected.iloc[j].precision_loss
                assert (pd.isna(actual) and np.isnan(expected)) or actual==expected
            evaluation=selected.iloc[ei]
            counts.append(dict(method=method,patients=24,loss_patients=int(evaluation.recall_loss.count()),
                precision_patients=int(evaluation.precision_loss.count()),correlation_patients=int(evaluation.pearson.count())))
            if 'reference_'+method in w:
                refined=w['reference_'+method]
                refined_threshold=quantile.threshold_reference(refined,w['neutral']&w['common'],ci)
                mismatch=(values<threshold)!=(refined<refined_threshold)
                boundary.append(dict(method=method,threshold=threshold,reference_threshold=refined_threshold,
                    all_callable_decision_disagreements=int((mismatch&w['common']).sum()),
                    comparison_decision_disagreements=int((mismatch[:,ei]&w['common'][:,ei]).sum())))
        comparisons=[]
        # Fixed comparator family; paired precision uses the same patient population.
        for candidate in ('ideal_constant','distance','distance_shuffled'):
            a=scores[scores.method.eq(candidate)&scores.split.eq('evaluation')].set_index('sample')
            for comparator in ('hd_posterior','naive_own_partition','distance_shuffled'):
                if candidate==comparator:continue
                b=scores[scores.method.eq(comparator)&scores.split.eq('evaluation')].set_index('sample').loc[a.index]
                for metric in ('recall_loss','fpr_loss','precision_loss','pearson'):
                    valid=a[metric].notna()&b[metric].notna()
                    differences=(a.loc[valid,metric]-b.loc[valid,metric]).tolist()
                    comparisons.append(dict(candidate=candidate,comparator=comparator,metric=metric,
                        paired_patients=len(differences),candidate_defined=int(a[metric].notna().sum()),
                        comparator_defined=int(b[metric].notna().sum()),mean_paired_effect=math.fsum(differences)/len(differences) if differences else None,
                        uncertainty='descriptive; no biological CI, independent evaluation or superiority test'))
        pd.DataFrame(comparisons).to_csv(args.out/'paired_effects.tsv',sep='\t',index=False)
        pd.DataFrame(counts).to_csv(args.out/'defined_populations.tsv',sep='\t',index=False)
        (args.out/'report.json').write_text(json.dumps(dict(anchor_medians_exact=len(anchors),blocks=len(endpoints),
            endpoint_checks=endpoints,method_patient_metric_checks=len(report['methods'])*48*3,
            threshold_boundaries=boundary,defined_populations=counts,paired_effects=comparisons,
            method_promotion=False,scientific_thresholds_changed=False,new_hypothesis_tests=0),indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
