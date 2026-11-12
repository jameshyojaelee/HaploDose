"""RNA-only segment selection against the complete frozen profile comparators."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path
import statistics

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('predictive_panel_'+name,path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checks(model):
    y = np.array([0.,0.,6.,6.])
    fitted,parts = model._model.partition(y,np.ones(4),np.array([True,False,False,False]))
    costs=[]
    for bits in itertools.product((False,True),repeat=3):
        edges = [0]+[i+1 for i,on in enumerate(bits) if on]+[4]
        value=0.
        for lo,hi in zip(edges,edges[1:]):
            centre=statistics.mean(float(x) for x in y[lo:hi])
            value+=sum((float(x)-centre)**2 for x in y[lo:hi])
        costs.append(value+3*np.log(4)*(len(edges)-2))
    reconstructed=sum(sum((float(x)-p['effect'])**2 for x in y[p['first_bin']:p['last_bin']+1]) for p in parts)
    reconstructed+=3*np.log(4)*(len(parts)-1)
    assert reconstructed==min(costs)
    np.testing.assert_array_equal(fitted,y)
    gbin=np.repeat(np.arange(7),4)
    values=np.repeat([0.,2.,0.,2.,0.,2.,0.],4)
    blocks=np.repeat('p',7)
    result,record,base=model.select(values,gbin,blocks)
    np.testing.assert_array_equal(result,values[::4])
    assert record['selected_model']=='running_median'
    constant,record,base=model.select(np.ones(len(values)),gbin,blocks)
    np.testing.assert_array_equal(constant,np.ones(7))
    assert record['status']=='segmentation_ineligible'
    folds=np.tile(np.arange(4),7)
    before=model.masked_segment(values,gbin,blocks,folds,0)[0]
    changed=values.copy();changed[folds==0]+=100.
    after=model.masked_segment(changed,gbin,blocks,folds,0)[0]
    np.testing.assert_array_equal(before,after)
    return dict(exhaustive_partition_optimum=True,known_noiseless_optimum=True,
                constant_scale_abstention=True,masked_value_perturbation_invariant=True)


def kinker(args,model,inputs):
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    base=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    prior=root/'dev/rna_predictive_smoothing_20261109/regression_v1'
    identities=json.loads((base/'predictions_locked.json').read_text())['cases']
    selected=identities if args.panel=='kinker' else [r for r in identities if r['case'] in
        ('2313287_STOMACH','42MGBA_CENTRAL_NERVOUS_SYSTEM')]
    assert len(selected)==(181 if args.panel=='kinker' else 2)
    inputs['identities']=base/'predictions_locked.json'
    inputs['truth']=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz'
    for r in selected:
        i=r['prediction_index'];inputs[f'gene{i}']=prior/f'predictions_case{i}.tsv.gz'
        inputs[f'bins{i}']=base/f'bins_case{i}.tsv'
    run=start_run(args.out,inputs,dict(panel=args.panel),dict(scope='Exposed Kinker regression; transformed RNA, no count calibration'),mode='predictive_partition_panel')
    try:
        numerical=checks(model)
        predictions,decisions=[],[]
        for r in selected:
            i=r['prediction_index'];gene=pd.read_csv(inputs[f'gene{i}'],sep='\t',float_precision='round_trip')
            bins=pd.read_csv(inputs[f'bins{i}'],sep='\t')
            result,decision,baseline=model.select(gene.rna_log_ratio.to_numpy(),gene.bin.to_numpy(),bins.smoothing_block.to_numpy())
            np.testing.assert_array_equal(baseline[gene.bin],gene.predictive_smoothed)
            gene['predictive_partition']=result[gene.bin.to_numpy()]
            gene.to_csv(args.out/f'predictions_case{i}.tsv.gz',sep='\t',index=False)
            decisions.append(dict(case=r['case'],prediction_index=i,**decision));predictions.append(gene)
        write_json(dict(cases=decisions,files=[file_record(args.out/f'predictions_case{r["prediction_index"]}.tsv.gz') for r in selected]),args.out/'RNA_DECISIONS_LOCKED.json')
        # DNA values enter only after all RNA decisions and their hashes exist.
        cn=pd.read_csv(inputs['truth'],sep='\t',index_col=0)
        methods=['pseudobulk','pseudobulk_sm5','pseudobulk_mdl','predictive_smoothed','predictive_partition']
        scores=[]
        for gene,r in zip(predictions,selected):
            assert len(r['query_DNA_line'])==1
            weights=pd.Series(r['reference_DNA_weights']);g=gene[gene.gene.isin(cn.index)].copy()
            g['truth']=cn.loc[g.gene,r['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene,weights.index].to_numpy()@weights.to_numpy()
            assert np.isfinite(g[['truth']+methods].to_numpy()).all()
            g['mb']=g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1_000_000).astype(str)
            mb=g.groupby('mb')[['truth']+methods].mean()
            mb.to_csv(args.out/f'windows_case{r["prediction_index"]}.tsv',sep='\t')
            for name in methods:
                scores.append(dict(case=r['case'],method=name,pearson=float(pearsonr(mb[name],mb.truth).statistic),
                                   windows=len(mb),common_genes=len(g)))
        frame=pd.DataFrame(scores);frame.to_csv(args.out/'scores.tsv',sep='\t',index=False)
        matrix=frame.pivot(index='case',columns='method',values='pearson');paired=[]
        for name in methods[:-1]:
            delta=matrix.predictive_partition-matrix[name]
            assert np.isfinite(delta).all()
            paired.append(dict(baseline=name,cases=len(delta),mean=float(delta.mean()),median=float(delta.median()),
                wins=int((delta>0).sum()),losses=int((delta<0).sum()),ties=int((delta==0).sum())))
        write_json(dict(panel=args.panel,checks=numerical,cases=len(selected),paired=paired,
            method_means=frame.groupby('method').pearson.mean().to_dict(),
            selected_partition=sum(r['selected_model']=='piecewise_mean_mdl' for r in decisions),
            independent_evaluation=False,default_promoted=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU conditional RNA model selection and fixed Kinker scorer')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


def bulk(args,model,inputs):
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    native=root/'dev/native_predictive_profile_20261109/raw_bulk_v1/haplodose'
    prior=root/'dev/bulk_amplitude_reference_20261109/refit_000_v1'
    panel=root/'dev/hardening_20261106_2119/profile_v6'
    base=load('profile')
    inputs.update(frozen=prior/'window_scores.npz',mapping=prior/'gene_map.tsv',split=panel/'split_manifest.tsv',
                  scorer=Path(base.__file__),transitive_scorer=Path(base.subarm.__file__),
                  statistics=Path(__file__).with_name('statistics.py'))
    for name in ('gene_rna_log_ratio.npy','rna_log_ratio_predictive.npy','genes.tsv.gz','bins.tsv','samples.tsv'):
        inputs[name]=native/name
    run=start_run(args.out,inputs,{},dict(scope='Exposed raw-bulk postprocessing; no new DNA or biological replication'),mode='predictive_partition_bulk')
    try:
        bins=pd.read_csv(inputs['bins.tsv'],sep='\t');genes=pd.read_csv(inputs['genes.tsv.gz'],sep='\t')
        samples=pd.read_csv(inputs['samples.tsv'],sep='\t',dtype={'sample':str})
        values=np.load(inputs['gene_rna_log_ratio.npy']);old=np.load(inputs['rna_log_ratio_predictive.npy'])
        blocks=np.cumsum(chromosome_starts(bins)|arm_starts(bins,genome_build='hg38'))-1
        fitted,records,median1,median5=[],[],[],[]
        from haplodose.smoothing import _running_median
        for row,sample in enumerate(samples['sample']):
            result,decision,baseline=model.select(values[row],genes.bin.to_numpy(),blocks)
            np.testing.assert_array_equal(baseline,old[row])
            full=model.medians(values[row],genes.bin.to_numpy(),np.ones(len(genes),bool),len(bins))
            median1.append(full);median5.append(_running_median(full,blocks,5));fitted.append(result)
            records.append(dict(sample=sample,**decision))
        fitted=np.array(fitted);np.save(args.out/'RNA_predictions.npy',fitted)
        write_json(dict(cases=records,predictions=file_record(args.out/'RNA_predictions.npy')),args.out/'RNA_DECISIONS_LOCKED.json')
        # All RNA decisions are frozen before loading the DNA/scoring archive.
        with np.load(inputs['frozen']) as z: source={key:z[key] for key in z.files}
        mapping=pd.read_csv(inputs['mapping'],sep='\t')
        methods=[key for key in source if key not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(methods)==12
        windows={name:source[name] for name in methods}
        gene_bin=genes.set_index('gene').bin
        for name,array in [('gene_predictive',old),('gene_median_w1',np.array(median1)),
                           ('gene_median_w5',np.array(median5)),('predictive_partition',fitted)]:
            aligned=base.subarm.align_rows(array,samples['sample'],source['samples'])
            mapped=base.subarm.gene_values(aligned,gene_bin,mapping.gene.to_numpy())
            windows[name]=base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None],mapped,np.nan),
                mapping.window.to_numpy(int),source['dna'].shape[0])
            assert np.isfinite(windows[name][source['common']]).all(),name
        split=pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[source['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split,source['split'])
        calibration=np.flatnonzero(split.split.eq('calibration'));assert len(calibration)==24
        scores=[]
        for name,array in windows.items():
            threshold=base.loss_threshold(array,source['neutral']&source['common'],calibration)
            loss=base.score_loss(array,source['loss'],source['neutral'],source['common'],threshold,
                                 negative=source['gain']|source['neutral'])
            for j,sample in enumerate(source['samples']):
                x,y=array[source['common'][:,j],j],source['dna'][source['common'][:,j],j]
                correlation=float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name=='predictive_partition' and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(correlation):
                    raise ArithmeticError('Candidate loses a truth-informative profile; cannot discard this patient')
                scores.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=correlation,**loss[j]))
        scores=pd.DataFrame(scores);scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        paired=[];evaluation=scores[scores.split.eq('evaluation')]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix=evaluation.pivot(index='sample',columns='method',values=metric)
            for name in [*methods,'gene_predictive','gene_median_w1','gene_median_w5']:
                delta=(matrix.predictive_partition-matrix[name]).dropna()
                improvement=-delta if metric=='fpr_loss' else delta
                paired.append(dict(metric=metric,baseline=name,patients=len(delta),mean=float(delta.mean()),
                    wins=int((improvement>0).sum()),losses=int((improvement<0).sum()),ties=int((delta==0).sum()),
                    higher_is_better=metric!='fpr_loss'))
        write_json(dict(panel='bulk',paired=paired,original_comparators=methods,
            selected_partition=sum(r['selected_model']=='piecewise_mean_mdl' for r in records),
            all_existing_callable_support_retained=True,independent_evaluation=False,default_promoted=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU conditional RNA segmentation and fixed raw-bulk window scorer')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--panel',choices=('checks','locked','kinker','bulk'),required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    model=load('rna_predictive_partition')
    inputs=dict(script=Path(__file__),model=Path(model.__file__),partition_source=model._path,
                design=args.task/'DESIGN.md')
    if args.panel=='checks':
        run=start_run(args.out,inputs,{},dict(scope='Numerical fixtures only'),mode='predictive_partition_checks')
        try:
            write_json(checks(model),args.out/'report.json');finish_run(args.out,run,backend='CPU exhaustive/scalar fixtures')
        except BaseException as error:
            finish_run(args.out,run,error=error);raise
    elif args.panel=='bulk':bulk(args,model,inputs)
    else:kinker(args,model,inputs)


if __name__=='__main__':main()
