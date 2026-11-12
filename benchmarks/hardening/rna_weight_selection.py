"""One conditional RNA-risk selector between two frozen regional medians."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import statistics

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec=importlib.util.spec_from_file_location('rna_weight_selection_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def identities(values, variance, gene_ids, gene_bins, path, reset):
    values,variance,gene_ids,gene_bins,path,reset=map(np.asarray,(values,variance,gene_ids,gene_bins,path,reset))
    if (values.ndim!=1 or any(a.shape!=values.shape for a in (variance,gene_ids,gene_bins))
            or len(np.unique(gene_ids))!=len(gene_ids) or path.ndim!=1 or not len(path)
            or reset.shape!=path.shape or not reset[0] or not np.isfinite(values).all()
            or not np.isfinite(variance).all() or np.any(variance<=0)
            or not np.array_equal(np.unique(gene_bins),np.arange(len(path)))):
        raise ValueError('Unique aligned genes, complete bins, finite effects and positive variance required')
    groups=np.cumsum(np.r_[True,reset[1:]|(path[1:]!=path[:-1])])-1
    gene_groups=groups[gene_bins]
    folds=np.empty(len(values),int)
    for g in np.unique(gene_groups):
        ix=np.flatnonzero(gene_groups==g)
        order=ix[np.argsort(gene_ids[ix],kind='stable')]
        folds[order]=np.arange(len(order))%10
    return groups,gene_groups,folds


def training_predictions(values, variance, gene_groups, folds, excluded):
    weighted=load('reference_weighted_regions').weighted_median
    result=np.full((2,len(values)),np.nan)
    for g in np.unique(gene_groups):
        region=gene_groups==g
        train=region&(folds!=excluded)
        hold=region&(folds==excluded)
        if train.any() and hold.any():
            result[0,hold]=np.median(values[train])
            result[1,hold]=weighted(values[train],1/np.sqrt(variance[train]))
    return result


def select(values, variance, gene_ids, gene_bins, path, reset):
    values,variance,gene_ids,gene_bins,path,reset=map(np.asarray,(values,variance,gene_ids,gene_bins,path,reset))
    groups,gene_groups,folds=identities(values,variance,gene_ids,gene_bins,path,reset)
    predictions=np.full((2,len(values)),np.nan)
    for fold in range(10):
        pred=training_predictions(values,variance,gene_groups,folds,fold)
        predictions[:,folds==fold]=pred[:,folds==fold]
    eligible=np.isfinite(predictions).all(axis=0)
    if not eligible.any():
        raise ValueError('No independently held gene has training support in its fixed region')
    order=np.flatnonzero(eligible)[np.argsort(gene_ids[eligible],kind='stable')]
    risks=np.array([math.fsum(map(float,abs(p[order]-values[order]))) for p in predictions])
    operations=4*len(order)+32
    eps=np.finfo(float).eps
    magnitude=max(1.,float(risks.max()),math.fsum(map(float,abs(values[order]))))
    budget=operations*eps/(1-operations*eps)*magnitude
    chosen=int(risks[1]<risks[0]-2*budget)
    weighted=load('reference_weighted_regions').weighted_median
    fitted=np.empty(len(path))
    for group in np.unique(groups):
        keep=gene_groups==group
        fitted[groups==group]=(weighted(values[keep],1/np.sqrt(variance[keep])) if chosen else np.median(values[keep]))
    return fitted,dict(selected=('uniform','reference_weighted')[chosen],risk_sums=risks.tolist(),
        risk_rounding_budget=float(budget),testable_genes=int(eligible.sum()),excluded_genes=int((~eligible).sum())),folds


def checks():
    weighted=load('reference_weighted_regions')
    old=weighted.checks()
    y=np.r_[np.zeros(5),np.full(5,10.)]
    v=np.r_[np.full(5,.01),np.ones(5)]
    ids=np.array(list('abcdefghij'))
    bins=np.zeros(10,int);path=np.array([1]);reset=np.array([True])
    result,record,folds=select(y,v,ids,bins,path,reset)
    assert record['risk_sums']==[100.,50.] and record['selected']=='reference_weighted' and result[0]==0
    for f in range(10):
        pred=training_predictions(y,v,np.zeros(10,int),folds,f)
        hold=folds==f
        expected=statistics.median(map(float,y[~hold]))
        assert pred[0,hold][0]==expected
        changed=y.copy();changed[hold]+=1000
        np.testing.assert_array_equal(pred,training_predictions(changed,v,np.zeros(10,int),folds,f))
    order=np.array([9,0,5,2,6,4,1,8,3,7])
    alternate,other,other_folds=select(y[order],v[order],ids[order],bins[order],path,reset)
    np.testing.assert_array_equal(result,alternate);assert other==record
    np.testing.assert_array_equal(other_folds,folds[order])
    _,constant,_=select(np.ones(10),v,ids,bins,path,reset)
    assert constant['selected']=='uniform' and constant['risk_sums']==[0.,0.]
    fit,_,_=select(np.r_[np.zeros(5),np.full(5,100.)],np.ones(10),ids,np.repeat([0,1],5),
        np.array([1,1]),np.array([True,True]))
    np.testing.assert_array_equal(fit,[0.,100.])
    try:
        select(y,v,np.repeat('duplicate',10),bins,path,reset)
    except ValueError:
        pass
    else:
        raise AssertionError('Duplicate identities must fail')
    return dict(prior_weight_numerics=old,independent_scalar_risks=True,held_out_perturbation_exact=True,
        gene_permutation_exact=True,gap_resets=True,constant_selects_uniform=True,duplicate_identity_rejected=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true')
    parser.add_argument('--reference-refit',type=int,choices=range(9),default=0)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:
        write_json(checks(),args.out/'report.json');return
    root=args.task.parents[1]
    number=args.reference_refit
    prior=root/f'dev/programme_gene_profile_20261109/refit_{number:03d}_v1'
    suffix='bulk_v1' if number==0 else f'refit_{number:03d}_v1'
    ref=root/'dev/reference_weighted_regions_20261109'/suffix
    old=root/'dev/count_weighted_regions_20261109'/suffix
    manifest=json.loads((prior/'run_manifest.json').read_text())
    native=Path(manifest['input']['files']['genes.tsv.gz']['path']).parent
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',prior_manifest=prior/'run_manifest.json',
        weights=ref/'reference_weights.npz',uniform=prior/'RNA_predictions.npz',weighted=ref/'RNA_predictions.npz',
        previous_windows=old/'window_scores.npz',previous_scores=old/'patient_scores.tsv',
        genes=native/'genes.tsv.gz',bins=native/'bins.tsv',
        split=root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv',
        mapping=root/f'dev/bulk_amplitude_reference_20261109/refit_{number:03d}_v1/gene_map.tsv',
        samples=Path(manifest['input']['files']['samples']['path']))
    for i in range(68):inputs[f'context{i}']=prior/f'context_{i:03d}.npz'
    base=load('profile');inputs.update(scorer=Path(base.__file__),transitive_scorer=Path(base.subarm.__file__),
        weighted_reference=Path(__file__).with_name('reference_weighted_regions.py'))
    run=start_run(args.out,inputs,dict(folds=10,reference_refit=number),dict(scope='Conditional RNA selection, regression-only DNA characterization'),mode='rna_weight_selection')
    try:
        numerical=checks()
        genes=pd.read_csv(inputs['genes'],sep='\t');bins=pd.read_csv(inputs['bins'],sep='\t')
        samples=pd.read_csv(inputs['samples'],sep='\t',dtype=str)
        with np.load(inputs['weights']) as z:
            normal=z['normal'];variance=z['tumour_variance'];np.testing.assert_array_equal(z['sample_ids'],samples['sample'])
            np.testing.assert_array_equal(z['genes'],genes.gene)
        with np.load(inputs['uniform']) as z:uniform=z['programme_gene_native_regions']
        with np.load(inputs['weighted']) as z:weighted=z['reference_weighted_median']
        from haplodose.bins import chromosome_starts,arm_starts
        reset=chromosome_starts(bins)|arm_starts(bins,genome_build='hg38')
        candidate=np.full(uniform.shape,np.nan);records=[]
        for j,i in enumerate(np.flatnonzero(~normal)):
            with np.load(inputs[f'context{i}']) as z:
                fitted,record,folds=select(z['corrected_gene'],variance[j],genes.gene.to_numpy(str),
                    genes.bin.to_numpy(int),z['path'],reset)
            expected=uniform[i] if record['selected']=='uniform' else weighted[i]
            np.testing.assert_array_equal(fitted,expected)
            candidate[i]=fitted;records.append(dict(sample=str(samples.iloc[i]['sample']),**record))
        method='rna_selected_weighted_median'
        np.savez_compressed(args.out/'RNA_predictions.npz',**{method:candidate})
        write_json(dict(records=records,numerical=numerical,predictions=file_record(args.out/'RNA_predictions.npz')),
            args.out/'RNA_DECISIONS_LOCKED.json')
        with np.load(inputs['previous_windows']) as z:windows={k:z[k] for k in z.files}
        old_methods=[k for k in windows if k not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(old_methods)==23
        mapping=pd.read_csv(inputs['mapping'],sep='\t')
        aligned=base.subarm.align_rows(candidate,samples['sample'],windows['samples'])
        mapped=base.subarm.gene_values(aligned,genes.set_index('gene').bin,mapping.gene.to_numpy())
        windows[method]=base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None],mapped,np.nan),
            mapping.window.to_numpy(int),windows['dna'].shape[0])
        np.savez_compressed(args.out/'window_scores.npz',**windows)
        split=pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[windows['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split,windows['split'])
        cal=np.flatnonzero(split.split.eq('calibration'));rows=[]
        assert len(cal)==24 and split.split.eq('evaluation').sum()==24
        for name in old_methods+[method]:
            values=windows[name]
            threshold=base.loss_threshold(values,windows['neutral']&windows['common'],cal)
            losses=base.score_loss(values,windows['loss'],windows['neutral'],windows['common'],threshold,
                negative=windows['gain']|windows['neutral'])
            for j,sample in enumerate(windows['samples']):
                x,y=values[windows['common'][:,j],j],windows['dna'][windows['common'][:,j],j]
                r=float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name==method and len(x)>=3 and np.ptp(y)>0:assert np.isfinite(r)
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=r,**losses[j]))
        scores=pd.DataFrame(rows);scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores[scores.method.isin(old_methods)].to_csv(args.out/'old_method_replay.tsv',sep='\t',index=False)
        pd.testing.assert_frame_equal(pd.read_csv(inputs['previous_scores'],sep='\t',float_precision='round_trip'),
            pd.read_csv(args.out/'old_method_replay.tsv',sep='\t',float_precision='round_trip'),check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        selected=scores[scores.split.eq('evaluation')];contrasts=[]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix=selected.pivot(index='sample',columns='method',values=metric)
            for baseline in old_methods:
                delta=(matrix[method]-matrix[baseline]).dropna()
                signed=-delta if metric=='fpr_loss' else delta
                contrasts.append(dict(baseline=baseline,metric=metric,patients=len(delta),mean=float(delta.mean()),
                    wins=int((signed>0).sum()),losses=int((signed<0).sum()),ties=int((delta==0).sum())))
        choices={k:sum(r['selected']==k for r in records) for k in ('uniform','reference_weighted')}
        gate=all(choices.values()) and all(next(r['mean'] for r in contrasts if r['baseline']==baseline and r['metric']=='pearson')>0
            for baseline in ('programme_gene_native_regions','reference_weighted_median','naive_own_partition'))
        write_json(dict(choices=choices,contrasts=contrasts,continuation_gate_passed=bool(gate),old_methods_exact=23,
            methods=old_methods+[method],checks=numerical,independent_evaluation=False,production_changed=False,
            scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU conditional held-gene selection and unchanged DNA scorer')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
