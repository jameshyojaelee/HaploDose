"""Same frozen RNA stack on raw-bulk normal-centered gene observations."""
import argparse
from functools import lru_cache
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from haplodose.bins import chromosome_starts, arm_starts
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('bulk_stack_'+name,Path(__file__).with_name(name+'.py'))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parents[1];prior=root/'dev/programme_gene_profile_20261109/refit_000_v1'
    compare=root/'dev/posterior_region_average_20261109/bulk_v1'
    source=json.loads((prior/'run_manifest.json').read_text())['input']['files']
    native=Path(source['genes.tsv.gz']['path']).parent
    inputs=dict(script=Path(__file__),design=args.task/'REGRESSION_DESIGN.md',
        stack=Path(__file__).with_name('rna_stacked_smoothing.py'),
        weighted=Path(__file__).with_name('reference_weighted_regions.py'),
        partition=Path(__file__).with_name('rna_predictive_partition.py'),
        partition_engine=Path(__file__).with_name('performance_candidates.py'),
        scorer=Path(__file__).with_name('profile.py'),transitive_scorer=Path(__file__).parents[1]/'tcga_bulk/subarm.py',
        previous_manifest=prior/'run_manifest.json',samples=Path(source['samples']['path']),
        genes=native/'genes.tsv.gz',bins=native/'bins.tsv',windows=compare/'window_scores.npz',
        previous_scores=compare/'patient_scores.tsv',
        mapping=root/'dev/bulk_amplitude_reference_20261109/refit_000_v1/gene_map.tsv',
        split=root/'dev/hardening_20261106_2119/profile_v6/split_manifest.tsv',
        locked_gate=args.task/'locked_v2/report.json')
    for i in range(68):inputs[f'context{i}']=prior/f'context_{i:03d}.npz'
    rec=start_run(args.out,inputs,dict(expert_bank='same as Kinker, original median widths plus MDL',
        reference_refit=0),dict(scope='Raw-bulk regression; no independent DNA or calibration claim'),mode='bulk_stacked_smoothing')
    try:
        assert json.loads(inputs['locked_gate'].read_text())['locked_two_case_gate']
        genes=pd.read_csv(inputs['genes'],sep='\t');bins=pd.read_csv(inputs['bins'],sep='\t')
        samples=pd.read_csv(inputs['samples'],sep='\t',dtype=str)
        effects=[]
        for i in range(len(samples)):
            with np.load(inputs[f'context{i}']) as z:effects.append(z['corrected_gene'])
        effects=np.asarray(effects);normal=samples.role.eq('normal').to_numpy()
        assert len(samples)==68 and normal.sum()==20
        reset=chromosome_starts(bins)|arm_starts(bins,genome_build='hg38');blocks=np.cumsum(reset)-1
        prediction=np.full((68,len(bins)),np.nan);decisions=[]
        for i in np.flatnonzero(~normal):
            effect,_,_=load('reference_weighted_regions').normal_centered_effect(effects[i],effects[normal],genes.bin.to_numpy())
            value,decision=load('rna_stacked_smoothing').select(effect,genes.bin.to_numpy(),blocks)
            prediction[i]=value;decisions.append(dict(sample=str(samples.iloc[i]['sample']),**decision))
        np.savez_compressed(args.out/'RNA_predictions.npz',stacked_normal_centered=prediction)
        write_json(dict(cases=decisions,predictions=file_record(args.out/'RNA_predictions.npz')),
            args.out/'RNA_DECISIONS_LOCKED.json')
        # All RNA model fitting is complete before corresponding DNA is opened.
        with np.load(inputs['windows']) as z:truth={k:z[k] for k in z.files}
        base=load('profile');mapping=pd.read_csv(inputs['mapping'],sep='\t')
        originals=[k for k in truth if k not in ('samples','common','loss','neutral','gain','dna','split')]
        assert len(originals)==39
        aligned=base.subarm.align_rows(prediction,samples['sample'],truth['samples'])
        mapped=base.subarm.gene_values(aligned,genes.set_index('gene').bin,mapping.gene.to_numpy())
        new=base.subarm.window_median(np.where(mapping.shared.to_numpy(bool)[:,None],mapped,np.nan),mapping.window.to_numpy(int),truth['dna'].shape[0])
        assert np.isfinite(new[truth['common']]).all()
        windows={k:truth[k] for k in originals};windows['stacked_normal_centered']=new
        np.savez_compressed(args.out/'window_scores.npz',**truth,stacked_normal_centered=new)
        split=pd.read_csv(inputs['split'],sep='\t',dtype=str).set_index('sample').loc[truth['samples']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split,truth['split'])
        calibration=np.flatnonzero(split.split.eq('calibration'));assert len(calibration)==24
        rows=[]
        for name,array in windows.items():
            threshold=base.loss_threshold(array,truth['neutral']&truth['common'],calibration)
            losses=base.score_loss(array,truth['loss'],truth['neutral'],truth['common'],threshold,negative=truth['gain']|truth['neutral'])
            for j,sample in enumerate(truth['samples']):
                x,y=array[truth['common'][:,j],j],truth['dna'][truth['common'][:,j],j]
                correlation=float(pearsonr(x,y).statistic) if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                if name=='stacked_normal_centered' and len(x)>=3 and np.ptp(y)>0 and not np.isfinite(correlation):
                    raise ArithmeticError('Stack loses an informative patient')
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=correlation,**losses[j]))
        scores=pd.DataFrame(rows);scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        old=pd.read_csv(inputs['previous_scores'],sep='\t',float_precision='round_trip')
        after=pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(old,after[after.method.isin(originals)].reset_index(drop=True),check_exact=True)
        scores.groupby(['split','method'])[['pearson','recall_loss','precision_loss','fpr_loss']].agg(['mean','count']).to_csv(args.out/'summary.tsv',sep='\t')
        paired=[]
        for metric in ('pearson','recall_loss','precision_loss','fpr_loss'):
            matrix=scores[scores.split.eq('evaluation')].pivot(index='sample',columns='method',values=metric)
            for name in originals:
                delta=(matrix.stacked_normal_centered-matrix[name]).dropna()
                sign=-delta if metric=='fpr_loss' else delta
                paired.append(dict(baseline=name,metric=metric,n=len(delta),mean=float(delta.mean()),
                    wins=int((sign>0).sum()),losses=int((sign<0).sum()),ties=int((sign==0).sum())))
        required=('programme_gene_native_regions','naive_own_partition','normal_centered_reference_median',
            'gaussian_region_reference_median','centered_gaussian_reference_median','decomposed_gaussian_reference_median',
            'averaged_gaussian_reference_median')
        gate=all(next(r['mean'] for r in paired if r['baseline']==name and r['metric']=='pearson')>0 for name in required)
        write_json(dict(methods=40,old_scores_exact=39,paired=paired,positive_gain_gate=gate,
            stack_selected=sum(r['selected'] for r in decisions),samples=48,
            independent_evaluation=False,biological_calibration=False,production_changed=False,
            scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU frozen conditional RNA stack,40 matched bulk comparators')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
