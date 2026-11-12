"""Locked profile and bulk performance screens for the gene-Laplace expert."""
import argparse
from functools import lru_cache
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.stacking import stacked_rna_profile

OLD = ('pseudobulk','pseudobulk_sm5','pseudobulk_mdl','predictive_smoothed','predictive_partition','stacked_smoothed')
LOCKED = ('2313287_STOMACH','42MGBA_CENTRAL_NERVOUS_SYSTEM')


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('laplace_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sc(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    prior = root/'dev/rna_stacked_smoothing_20261110'
    identities = json.loads((base/'predictions_locked.json').read_text())['cases']
    selected = identities[args.start:args.stop] if args.all_lines else [r for r in identities if r['case'] in LOCKED]
    assert len(selected) == (args.stop-args.start if args.all_lines else 2)
    inputs = dict(script=Path(__file__), model=Path(__file__).with_name('gene_laplace_profile.py'),
        design=args.task/'DESIGN.md', identities=base/'predictions_locked.json',
        truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz',
        previous=prior/'summary_v3/scores.tsv')
    if args.all_lines:
        inputs['regression_design'] = args.task/'REGRESSION_DESIGN.md'
        inputs['locked_gate'] = args.task/'sc_v1/report.json'
        assert json.loads(inputs['locked_gate'].read_text())['gate']
    for r in selected:
        i = r['prediction_index']
        chunk = max(start for start in (0,46,92,138) if start <= i)
        inputs[f'RNA{i}'] = prior/f'regression_{chunk:03d}_v3/predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
    rec = start_run(args.out, inputs, dict(all_lines=args.all_lines, start=args.start, stop=args.stop),
        dict(scope='Former evaluation regression-only; prepared reconstructed RNA; no raw-count calibration'), mode='gene_laplace_sc')
    try:
        decisions, predictions = [], []
        for r in selected:
            i = r['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'],sep='\t',float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'],sep='\t')
            assert gene.gene.is_unique and gene.case.eq(r['case']).all()
            old, _ = stacked_rna_profile(gene.rna_log_ratio.to_numpy(),gene.bin.to_numpy(),bins.smoothing_block.to_numpy())
            np.testing.assert_array_equal(old[gene.bin.to_numpy()],gene.stacked_smoothed)
            new, decision, standalone = load('gene_laplace_profile').select(
                gene.rna_log_ratio.to_numpy(),gene.bin.to_numpy(),bins.smoothing_block.to_numpy())
            gene['gene_laplace_stack'] = new[gene.bin.to_numpy()]
            gene['gene_laplace_alone'] = standalone[gene.bin.to_numpy()] if standalone is not None else np.nan
            gene.to_csv(args.out/f'predictions_case{i}.tsv.gz',sep='\t',index=False)
            decisions.append(dict(case=r['case'],prediction_index=i,**decision))
            predictions.append(gene)
        write_json(dict(cases=decisions,outputs=[file_record(args.out/f"predictions_case{r['prediction_index']}.tsv.gz") for r in selected]),args.out/'RNA_LOCKED.json')
        cn = pd.read_csv(inputs['truth'],sep='\t',index_col=0)
        assert cn.index.is_unique and cn.columns.is_unique
        methods = OLD+('gene_laplace_alone','gene_laplace_stack')
        rows = []
        for r,gene in zip(selected,predictions):
            weights = pd.Series(r['reference_DNA_weights'])
            g = gene[gene.gene.isin(cn.index)].copy()
            g['truth'] = cn.loc[g.gene,r['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene,weights.index].to_numpy()@weights.to_numpy()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1_000_000).astype(str)
            mb = g.groupby('mb')[['truth']+list(methods)].mean()
            for name in methods:
                finite = np.isfinite(mb[[name,'truth']]).all().all() and mb[name].nunique()>1
                score = float(pearsonr(mb[name],mb.truth).statistic) if finite else np.nan
                rows.append(dict(case=r['case'],method=name,pearson=score,windows=len(mb),common_genes=len(g)))
        scores = pd.DataFrame(rows)
        scores.to_csv(args.out/'scores.tsv',sep='\t',index=False)
        actual = pd.read_csv(args.out/'scores.tsv',sep='\t',float_precision='round_trip')
        previous = pd.read_csv(inputs['previous'],sep='\t',float_precision='round_trip')
        before = previous[previous.case.isin([r['case'] for r in selected])].reset_index(drop=True)
        pd.testing.assert_frame_equal(before,actual[actual.method.isin(OLD)].reset_index(drop=True),check_exact=True)
        matrix = scores.pivot(index='case',columns='method',values='pearson')
        gate = len(selected)==2 and all((matrix.gene_laplace_stack>matrix[name]).all() for name in ('stacked_smoothed','pseudobulk_sm5'))
        write_json(dict(cases=len(selected),means=matrix.mean().to_dict(),gate=bool(gate),
            paired=[dict(case=c,current=float(matrix.loc[c,'stacked_smoothed']),candidate=float(matrix.loc[c,'gene_laplace_stack']),
                         sm5=float(matrix.loc[c,'pseudobulk_sm5']))for c in matrix.index],
            prior_scores_exact=len(before),candidate_selected=sum(r['selected']for r in decisions),
            standalone_is_ablation_only=True,independent_evaluation=False,biological_calibration=False,
            production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU exact integer gene-Laplace partitions and conditional sparse L1 stack')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


def bulk(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/programme_gene_profile_20261109/refit_000_v1'
    comparison = root/'dev/rna_stacked_smoothing_20261110/bulk_v3'
    source = json.loads((prior/'run_manifest.json').read_text())['input']['files']
    native = Path(source['genes.tsv.gz']['path']).parent
    inputs = dict(script=Path(__file__),model=Path(__file__).with_name('gene_laplace_profile.py'),
        design=args.task/'DESIGN.md',samples=Path(source['samples']['path']),
        genes=native/'genes.tsv.gz',bins=native/'bins.tsv',
        scorer=Path(__file__).with_name('profile.py'),transitive_scorer=Path(__file__).parents[1]/'tcga_bulk/subarm.py',
        normal_centering=Path(__file__).with_name('reference_weighted_regions.py'),
        previous_manifest=prior/'run_manifest.json',windows=comparison/'window_scores.npz',
        previous_scores=comparison/'patient_scores.tsv',
        mapping=root/'dev/bulk_amplitude_reference_20261109/refit_000_v1/gene_map.tsv',
        split=root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv')
    for i in range(68):
        inputs[f'context{i}'] = prior/f'context_{i:03d}.npz'
    rec = start_run(args.out,inputs,dict(experts='Original bank plus gene Laplace; no SC outcome dependence'),
        dict(scope='Raw-bulk prepared RNA effects; existing48 regression patients; no independent calibration'),mode='gene_laplace_bulk')
    try:
        genes,bins,samples = [pd.read_csv(inputs[k],sep='\t')for k in ('genes','bins','samples')]
        effects = []
        for i in range(68):
            with np.load(inputs[f'context{i}']) as z:
                effects.append(z['corrected_gene'])
        effects = np.array(effects)
        normal = samples.role.eq('normal').to_numpy()
        assert normal.sum()==20 and len(samples)==68
        blocks = np.cumsum(chromosome_starts(bins)|arm_starts(bins,genome_build='hg38'))-1
        predictions = {name:np.full((68,len(bins)),np.nan)for name in ('gene_laplace_stack','gene_laplace_alone')}
        decisions = []
        for i in np.flatnonzero(~normal):
            values,_,_ = load('reference_weighted_regions').normal_centered_effect(
                effects[i],effects[normal],genes.bin.to_numpy())
            fit,decision,alone = load('gene_laplace_profile').select(values,genes.bin.to_numpy(),blocks)
            predictions['gene_laplace_stack'][i] = fit
            if alone is not None:
                predictions['gene_laplace_alone'][i] = alone
            decisions.append(dict(sample=str(samples.iloc[i]['sample']),**decision))
        np.savez_compressed(args.out/'RNA_predictions.npz',**predictions)
        write_json(dict(cases=decisions,predictions=file_record(args.out/'RNA_predictions.npz')),args.out/'RNA_LOCKED.json')
        with np.load(inputs['windows']) as z:
            truth = {k:z[k]for k in z.files}
        old_names = [k for k in truth if k not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(old_names)==40
        windows = {k:truth[k]for k in old_names}
        base = load('profile')
        mapping = pd.read_csv(inputs['mapping'],sep='\t')
        for name,prediction in predictions.items():
            aligned = base.subarm.align_rows(prediction,samples['sample'],truth['samples'])
            projected = base.subarm.gene_values(aligned,genes.set_index('gene').bin,mapping.gene.to_numpy())
            windows[name] = base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None],projected,np.nan),
                mapping.window.to_numpy(int),truth['dna'].shape[0])
        assert np.isfinite(windows['gene_laplace_stack'][truth['common']]).all()
        np.savez_compressed(args.out/'window_scores.npz',**truth,**{k:windows[k]for k in predictions})
        split = pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split,truth['split'])
        calibration = np.flatnonzero(split.split.eq('calibration'))
        assert len(calibration)==24
        rows = []
        for name,array in windows.items():
            threshold = base.loss_threshold(array,truth['neutral']&truth['common'],calibration)
            losses = base.score_loss(array,truth['loss'],truth['neutral'],truth['common'],threshold,
                negative=truth['gain']|truth['neutral'])
            for j,sample in enumerate(truth['samples']):
                x,y = array[truth['common'][:,j],j],truth['dna'][truth['common'][:,j],j]
                score = float(pearsonr(x,y).statistic) if len(x)>=3 and np.isfinite(x).all() and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name=='gene_laplace_stack' and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(score):
                    raise ArithmeticError('Candidate loses an informative patient')
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=score,**losses[j]))
        scores = pd.DataFrame(rows)
        scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        after = pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip')
        before = pd.read_csv(inputs['previous_scores'],sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(before,after[after.method.isin(old_names)].reset_index(drop=True),check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        matrix = scores[scores.split.eq('evaluation')].pivot(index='sample',columns='method',values='pearson')
        paired = []
        for name in old_names:
            delta = (matrix.gene_laplace_stack-matrix[name]).dropna()
            paired.append(dict(baseline=name,n=len(delta),mean=float(delta.mean()),wins=int((delta>0).sum()),losses=int((delta<0).sum())))
        required = ('programme_gene_native_regions','naive_own_partition','averaged_gaussian_reference_median')
        gate = all(next(r['mean']for r in paired if r['baseline']==name)>0 for name in required)
        write_json(dict(methods=42,previous_methods_exact=40,primary='gene_laplace_stack',
            paired=paired,gate=gate,standalone_is_ablation_only=True,
            candidate_selected=sum(r['selected']for r in decisions),patients=48,
            independent_evaluation=False,biological_calibration=False,production_changed=False,
            scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU gene-Laplace expert with42 retained bulk profiles')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    p.add_argument('--mode',choices=('numerics','sc','bulk'),required=True)
    p.add_argument('--all-lines',action='store_true')
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--stop',type=int,default=181)
    args = p.parse_args()
    args.out.mkdir(exist_ok=False)
    if args.mode=='numerics':
        write_json(load('gene_laplace_profile').checks(),args.out/'report.json')
    elif args.mode=='sc':
        sc(args)
    else:
        bulk(args)


if __name__=='__main__':
    main()
