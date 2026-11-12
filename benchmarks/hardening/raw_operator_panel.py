"""Shared RNA-only segmentation operator on the already exposed raw bulk panel.

This reuses nine frozen fits. It neither refits HaploDose nor creates a new
independent evaluation cohort. Finite reference-refit bootstrap intervals remain
approximate, with no empirical coverage guarantee.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    path=Path(__file__).with_name(name+'.py')
    spec=importlib.util.spec_from_file_location('operator_'+name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def operator(naive,hd,arms,partition):
    """No truth, calibration labels, or other subjects enter operator decisions."""
    naive,hd=np.asarray(naive,float),np.asarray(hd,float)
    if naive.shape!=hd.shape or len(arms)!=len(naive):raise ValueError('Unaligned RNA windows')
    result=[np.full_like(naive,np.nan),np.full_like(hd,np.nan)];details=[]
    for j in range(naive.shape[1]):
        idx=np.flatnonzero(np.isfinite(naive[:,j])&np.isfinite(hd[:,j]))
        if not len(idx):
            details.append(dict(sample_column=j,failure='no_common_RNA_support'));continue
        breaks=np.r_[True,(np.diff(idx)!=1)|(arms[idx[1:]]!=arms[idx[:-1]])]
        blocks=np.split(idx,np.flatnonzero(breaks)[1:])
        differences=np.concatenate([np.diff(naive[b,j]) for b in blocks if len(b)>1])
        sigma=np.median(np.abs(differences-np.median(differences)))/(np.sqrt(2)*norm.ppf(.75))
        if not np.isfinite(sigma) or sigma<=0:
            raise ValueError('Unidentified RNA noise scale; no substitute epsilon')
        penalty=3*sigma*sigma*np.log(len(idx));boundaries=[]
        for method,values in enumerate((naive,hd)):
            regions=[]
            for b in blocks:
                fitted,parts,_=partition(values[b,j],penalty,3)
                result[method][b,j]=fitted
                regions.extend([(int(b[a]),int(b[z-1])) for a,z in parts])
            boundaries.append(regions)
        details.append(dict(sample_column=j,sigma=float(sigma),penalty=float(penalty),
                            minimum_windows=3,boundaries=boundaries))
    return result,details


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frozen-panel',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(exist_ok=False);base=load('profile');op=load('profile_operator_control');op.verify()
    precision=load('pearson_reference');precision.verify()
    root=Path(__file__).parent
    inputs={name:root/name for name in ('raw_operator_panel.py','profile.py','profile_operator_control.py','statistics.py','pearson_reference.py')}
    inputs.update(subarm=root.parent/'tcga_bulk/subarm.py',
                  split=a.frozen_panel/'split_manifest.tsv',genes=a.frozen_panel/'gene_universe.txt',
                  support=a.frozen_panel/'raw_support_and_truth.npz',scores=a.frozen_panel/'raw_patient_scores.tsv')
    for refit in range(9):
        directory=a.frozen_panel/f'raw_refit_{refit:03d}'
        for name in ('counts.tsv','samples.tsv','positions.tsv','haplodose/samples.tsv','haplodose/segments.tsv',
                     'haplodose/genes.tsv.gz','haplodose/log_ratio.npy','haplodose/run_manifest.json'):
            inputs[f'{refit}:{name}']=directory/name
    rec=start_run(a.out,inputs,dict(bootstrap=10000,seed=base.SEED),
        dict(role='already_exposed_development',truth_resolution='gene-anchored approximately 10Mb DNA windows',
             validation='not independent; conditional on finite normal-reference refits'),mode='raw_shared_operator')
    try:
        genes=(a.frozen_panel/'gene_universe.txt').read_text().splitlines()
        meta=pd.read_csv(a.frozen_panel/'split_manifest.tsv',sep='\t').set_index('sample')
        with np.load(a.frozen_panel/'raw_support_and_truth.npz',allow_pickle=True) as z:
            frozen={k:z[k] for k in z.files}
        g37=base.subarm.hg19_genes();arms=base.subarm.hg19_arms();tab,seg=base.subarm.absolute()
        truth=(g37,tab,seg,*arms);refits=[];before_scores=None
        ci=ei=None
        for i in range(9):
            directory=a.frozen_panel/f'raw_refit_{i:03d}'
            d=base.profiles(directory,directory/'haplodose',truth,genes)
            np.testing.assert_array_equal(d['tum'],frozen['samples'])
            for key,oldkey in [('D','dna'),('FL','loss'),('NEU','neutral')]:
                np.testing.assert_array_equal(d[key],frozen[oldkey])
            common=np.isfinite(d['D'])&np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS])
            np.testing.assert_array_equal(common,frozen['common'])
            if i==0:
                split=meta.reindex(d['tum']);ci=np.flatnonzero(split.split.eq('calibration'));ei=np.flatnonzero(split.split.eq('evaluation'))
                if len(ci)!=24 or len(ei)!=24 or split.patient.duplicated().any():raise ValueError('Biological split identity changed')
                before_scores=pd.DataFrame(base.score_profiles(d,ci)[2])
                old=pd.read_csv(a.frozen_panel/'raw_patient_scores.tsv',sep='\t').drop(columns='split')
                # Independent re-reading of the frozen fits must reproduce their
                # scorer before adding anything. Allow only text roundtrip error.
                pd.testing.assert_frame_equal(before_scores.drop(columns='pearson_window'),old.drop(columns='pearson_window'),
                                              check_dtype=False,rtol=0,atol=64*np.finfo(float).eps)
                numerical=[]
                for row,previous in zip(before_scores.to_dict('records'),old.to_dict('records')):
                    j=int(np.flatnonzero(d['tum']==row['sample'])[0]);ok=common[:,j]
                    if not np.isfinite(row['pearson_window']):
                        if np.isfinite(previous['pearson_window']):raise ValueError('Pearson missingness changed')
                        continue
                    x,y=d['W'][row['method']][ok,j],d['D'][ok,j]
                    if np.ptp(x)==0:
                        if row['pearson_window']!=0 or previous['pearson_window']!=0:raise ValueError('Constant-vector scoring changed')
                        continue
                    ref=precision.reference(x,y)
                    for value in (row['pearson_window'],previous['pearson_window']):
                        if not ref['float32_lower']<=value<=ref['float32_upper']:
                            raise ArithmeticError('Historical/reconstructed Pearson outside independent arithmetic enclosure')
                    numerical.append(dict(sample=row['sample'],method=row['method'],historical=previous['pearson_window'],
                        reconstructed=row['pearson_window'],**ref))
                (a.out/'pearson_precision.json').write_text(json.dumps(numerical,indent=2))
            for method in ('naive','naive_smooth'):d['W'][method]=d['W'][method].astype(float)*np.log(2.)
            (naive,hd),details=operator(d['W']['naive'],d['W']['hd_lr'],d['warm'],op.partition)
            d['W'].update(naive_own_partition=naive,hd_own_partition=hd)
            refits.append(d)
            np.savez_compressed(a.out/f'refit_{i:03d}.npz',samples=d['tum'],arms=d['warm'],**d['W'])
            (a.out/f'operator_{i:03d}.json').write_text(json.dumps(details,indent=2))
            print(json.dumps(dict(reconstructed_refit=i,original_support_identical=True)),flush=True)
        base.METHODS=base.METHODS+('naive_own_partition','hd_own_partition')
        original_correlations=base.correlations
        def refined_correlations(d,common):
            result=original_correlations(d,common)
            for method in base.METHODS:
                for j in range(len(d['tum'])):
                    ok=common[:,j];x,y=d['W'][method][ok,j],d['D'][ok,j]
                    if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:
                        result[method]['pearson_window'][j]=precision.reference(x,y)['value']
            return result
        base.correlations=refined_correlations
        # Retain the entire original family; adding competitors cannot hide a failure.
        base.CONTRASTS=base.CONTRASTS+(('hd_own_partition','naive_own_partition'),('hd_own_partition','naive_smooth'))
        _,thresholds,records,common=base.score_profiles(refits[0],ci)
        np.testing.assert_array_equal(common,frozen['common'])
        patient=pd.DataFrame(records);patient['split']=patient['sample'].map(meta.split)
        patient.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
        comparisons=base.nested_comparisons(refits[0],refits[1:],ci,ei,10000,base.SEED)
        pd.DataFrame(comparisons).to_csv(a.out/'comparisons.tsv',sep='\t',index=False)
        summary=dict(methods=base.METHODS,contrasts=base.CONTRASTS,thresholds=thresholds,
            biological_patients=48,calibration_patients=24,previously_inspected_reserved_patients=24,
            reference_donors=20,reference_bootstrap_refits=8,family_size=24,
            uncertainty='nested donor-reference, calibration-patient and evaluation-patient resampling; approximate',
            multiplicity='Holm 24 contrasts; Bonferroni basic-bootstrap intervals',
            units='all RNA outputs natural log; DNA relative log2 remains unchanged',
            minimum_segment_resolution='3 gene-anchored approximately 10Mb windows, not 3 original RNA bins',
            interpretation='development smoothing diagnostic; no independent superiority or coverage claim',
            comparisons=comparisons)
        (a.out/'summary.json').write_text(json.dumps(summary,indent=2))
        # Input immutability is checked after the full analysis too.
        for record in rec['input']['files'].values():
            if file_record(record['path'])['sha256']!=record['sha256']:raise ValueError('Input changed during analysis')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
