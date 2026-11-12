"""Fixed-boundary continuous RNA amplitudes, preserving all exposed comparators."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from haplodose.bins import arm_starts
from haplodose.bulk import BULK_PROFILE,fit_sample,reference_model
from haplodose.device import select
from haplodose.io import write_json
from haplodose.profile import log_ratio
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    spec=importlib.util.spec_from_file_location('continuous_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def amplitudes(track,variance,path,reset):
    boundaries=np.r_[0,np.flatnonzero(reset[1:]|(path[1:]!=path[:-1]))+1,len(path)]
    ordinary=np.empty(len(path));weighted=np.empty(len(path));rows=[]
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        r=track[lo:hi];w=1/variance[lo:hi]
        if not np.isfinite(r).all() or not np.isfinite(w).all() or (w<=0).any():
            raise ValueError('Finite RNA effects and positive precision required')
        a=float(np.mean(r));b=float(np.sum(w*r)/np.sum(w))
        ordinary[lo:hi]=a;weighted[lo:hi]=b
        rows.append(dict(first_bin=int(lo),last_bin=int(hi-1),state=int(path[lo]),mean=a,precision_mean=b))
    return ordinary,weighted,rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    a=parser.parse_args();a.out.mkdir(exist_ok=False)
    root=a.task.parents[1];stage=root/'dev/bulk_posterior_export_20261109/stage_v1';hd=stage/'haplodose'
    old=root/'dev/bulk_geometry_20261109/panel_v1';panel=root/'dev/hardening_20261106_2119/profile_v6';inp=panel/'raw_refit_000'
    base=load('profile');quantile=load('bulk_posterior_stage_audit');pearson=load('pearson_reference')
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',scorer=Path(base.__file__),
                transitive_scorer=Path(base.subarm.__file__),rank_reference=Path(quantile.__file__),pearson=Path(pearson.__file__),
                previous=old/'window_scores.npz',previous_scores=old/'patient_scores.tsv',previous_report=old/'report.json',
                split=panel/'split_manifest.tsv',gene_universe=panel/'gene_universe.txt',
                gtf=base.subarm.GTF37,cytoband=base.subarm.CYTO_HG19,centromere=base.subarm.CEN38,
                dna_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                dna_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for directory,prefix,names in ((inp,'input',('counts.tsv','samples.tsv','positions.tsv')),
        (hd,'hd',('genes.tsv.gz','bins.tsv','samples.tsv','rna_log_effect.npy','rna_state_posterior.npy','rna_state_index.npy'))):
        inputs.update({prefix+':'+name:directory/name for name in names})
    rec=start_run(a.out,inputs,dict(boundaries='unchanged_HMM_Viterbi',amplitude_models=['mean','precision_mean'],parameter_tuning=False),
                  dict(scope='Exposed-panel conditional RNA-amplitude mechanism;no superiority or biological calibration'),mode='bulk_continuous_amplitude')
    try:
        backend=select('cpu');assert backend=='cpu'
        bins=pd.read_csv(hd/'bins.tsv',sep='\t');genes=pd.read_csv(hd/'genes.tsv.gz',sep='\t')
        meta=pd.read_csv(hd/'samples.tsv',sep='\t',float_precision='round_trip')
        counts=pd.read_csv(inp/'counts.tsv',sep='\t',index_col=0).loc[genes.gene,meta['sample']]
        aggregation=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),genes.bin)),shape=(len(genes),len(bins)))
        Y=np.asarray((aggregation.T@counts.to_numpy(float)).T)
        assert np.array_equal(Y.sum(axis=1),meta.input_binned_counts)
        reset=np.r_[True,(bins.chromosome.to_numpy()[1:]!=bins.chromosome.to_numpy()[:-1])|
                    (bins.genomic_block.to_numpy()[1:]!=bins.genomic_block.to_numpy()[:-1])]
        ref=reference_model(Y[meta.role.eq('normal')],10,arm_starts(bins,genome_build='hg38'))
        original_post=np.load(hd/'rna_state_posterior.npy');original_effect=np.load(hd/'rna_log_effect.npy');original_path=np.load(hd/'rna_state_index.npy')
        effects={name:np.full_like(Y,np.nan) for name in ('hd_path_mean','hd_path_precision_mean')};segments=[]
        for i in np.flatnonzero(meta.role.eq('tumour')):
            sample=str(meta.iloc[i]['sample']);fit=fit_sample(Y[i],*ref,reset,BULK_PROFILE)
            assert np.array_equal(fit['post'],original_post[i]) and np.array_equal(fit['dose'],original_effect[i]) and np.array_equal(fit['path'],original_path[i])
            track=log_ratio(Y[i],fit['E'])-fit['centre'];variance=fit['sd']**2+1/(Y[i]+.5)
            ordinary,weighted,parts=amplitudes(track,variance,fit['path'],reset)
            effects['hd_path_mean'][i]=ordinary;effects['hd_path_precision_mean'][i]=weighted
            for part in parts:segments.append(dict(sample=sample,**part))
            np.savez_compressed(a.out/(sample+'.npz'),Y=Y[i],E=fit['E'],centre=fit['centre'],sd=fit['sd'],
                track=track,variance=variance,path=fit['path'],reset=reset,mean=ordinary,precision_mean=weighted)
        pd.DataFrame(segments).to_csv(a.out/'segments.tsv',sep='\t',index=False)
        # All RNA prediction decisions and arrays are complete before reading DNA scoring data.
        truth=base.subarm.hg19_genes();calls,dna_segments=base.subarm.absolute()
        d=base.profiles(inp,hd,(truth,calls,dna_segments,*base.subarm.hg19_arms()),(panel/'gene_universe.txt').read_text().splitlines())
        with np.load(old/'window_scores.npz') as z:frozen={k:z[k] for k in z.files}
        assert np.array_equal(d['tum'],frozen['samples'])
        for key,new in (('FL','loss'),('NEU','neutral'),('FG','gain'),('D','dna')):assert np.array_equal(d[key],frozen[new],equal_nan=True)
        common=np.isfinite(d['D'])&np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS])
        assert np.array_equal(common,frozen['common'])
        methods={name:frozen[name] for name in json.loads((old/'report.json').read_text())['methods']}
        mapping=genes.set_index('gene').bin
        gene_map=pd.DataFrame(dict(gene=d['genes'],bin=mapping.reindex(d['genes']).to_numpy(),
                                   window=d['window_id'],shared=d['shared_gene_mapping']))
        gene_map.to_csv(a.out/'gene_map.tsv',sep='\t',index=False)
        for name,values in effects.items():
            aligned=base.subarm.align_rows(values,meta['sample'],d['tum'])
            mapped=base.subarm.gene_values(aligned,mapping,d['genes'])
            mapped=np.where(d['shared_gene_mapping'][:,None],mapped,np.nan)
            methods[name]=base.subarm.window_median(mapped,d['window_id'],len(d['warm']))
            assert np.isfinite(methods[name][common]).all()
        split=pd.read_csv(panel/'split_manifest.tsv',sep='\t').set_index('sample').loc[d['tum']]
        assert np.array_equal(split.split.to_numpy(str),frozen['split'])
        ci=np.flatnonzero(split.split.eq('calibration'));assert len(ci)==24
        thresholds={};rows=[];pearson.verify()
        for name,values in methods.items():
            threshold=base.loss_threshold(values,d['NEU']&common,ci)
            assert threshold==quantile.threshold_reference(values,d['NEU']&common,ci);thresholds[name]=threshold
            loss=base.score_loss(values,d['FL'],d['NEU'],common,threshold,negative=d['FG']|d['NEU'])
            for j,sample in enumerate(d['tum']):
                x,y=values[common[:,j],j],d['D'][common[:,j],j]
                correlation=pearson.reference(x,y)['value'] if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=correlation,**loss[j]))
        scores=pd.DataFrame(rows);scores.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
        scores.groupby(['split','method']).agg(patients=('sample','size'),loss_patients=('recall_loss','count'),
            precision_patients=('precision_loss','count'),correlation_patients=('pearson','count'),mean_recall=('recall_loss','mean'),
            mean_fpr=('fpr_loss','mean'),mean_precision=('precision_loss','mean'),mean_pearson=('pearson','mean')).to_csv(a.out/'summary.tsv',sep='\t')
        old_scores=pd.read_csv(old/'patient_scores.tsv',sep='\t',float_precision='round_trip').set_index(['sample','method'])
        for row in scores[scores.method.isin(old_scores.index.get_level_values('method'))].itertuples():
            original=old_scores.loc[(row.sample,row.method)]
            for key in ('recall_loss','fpr_loss','precision_loss','pearson'):
                x,y=getattr(row,key),original[key];assert (pd.isna(x) and pd.isna(y)) or x==y
        np.savez_compressed(a.out/'window_scores.npz',samples=d['tum'].astype(str),common=common,loss=d['FL'],neutral=d['NEU'],
                            gain=d['FG'],dna=d['D'],split=split.split.to_numpy(str),**methods)
        write_json(dict(exact_tumour_refits=48,reference_donors=20,calibration_patients=24,comparison_patients=24,
                        methods=list(methods),thresholds=thresholds,segments=len(segments),old_methods_exact=True,
                        score_units='New amplitudes:natural-log conditional RNA effect,not DNA copy number',new_hypothesis_tests=0,
                        production_changed=False,method_promoted=False),a.out/'report.json')
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(a.out,rec,backend=backend)
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
