"""Matched fixed-fit ablation of HMM context on the exposed raw bulk panel.

Nuisance parameters remain the HMM estimates. This is not a complete workflow
alternative and cannot establish an independently fitted local model's quality.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import softmax

from haplodose.bins import arm_starts, chromosome_starts
from haplodose.bulk import BULK_PROFILE, fit_sample, reference_model
from haplodose.device import select
from haplodose.kernels.likelihood import pln_logpmf
from haplodose.profile import LOG_DOSE, STATE_PRIOR
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec=importlib.util.spec_from_file_location('local_ablation_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parent.parent; panel=root/'dev/hardening_20261106_2119/profile_v6'
    stage=root/'dev/bulk_posterior_export_20261109/stage_v1';hd=stage/'haplodose';inp=panel/'raw_refit_000'
    base=load('profile');audit=load('bulk_posterior_stage_audit')
    inputs=dict(script=Path(__file__),design=args.task/'LOCAL_ABLATION_DESIGN.md',
        profile=Path(__file__).with_name('profile.py'),quantile_reference=Path(audit.__file__),
        subarm=Path(base.subarm.__file__),statistics=Path(__file__).with_name('statistics.py'),
        split=panel/'split_manifest.tsv',genes=panel/'gene_universe.txt',
        dna_gtf=base.subarm.GTF37,cytoband=base.subarm.CYTO_HG19,centromeres=base.subarm.CEN38,
        absolute_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        absolute_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for directory,prefix,names in ((inp,'input',('counts.tsv','samples.tsv','positions.tsv')),
        (stage,'stage',('window_stages.npz','patient_stage_scores.tsv')),
        (hd,'hd',('bins.tsv','genes.tsv.gz','samples.tsv','rna_log_effect.npy','rna_state_posterior.npy'))):
        inputs.update({prefix+':'+name:directory/name for name in names})
    run=start_run(args.out,inputs,dict(parameter_selection=False,complete_workflow_ablation=False),
        dict(scope='Fixed-fit context ablation; exposed development panel; no independent evaluation'),mode='bulk_local_ablation')
    try:
        backend=select('cpu');assert backend=='cpu'
        bins=pd.read_csv(hd/'bins.tsv',sep='\t');gene=pd.read_csv(hd/'genes.tsv.gz',sep='\t')
        meta=pd.read_csv(hd/'samples.tsv',sep='\t',float_precision='round_trip')
        counts=pd.read_csv(inp/'counts.tsv',sep='\t',index_col=0).loc[gene.gene,meta['sample']]
        a=sparse.csr_matrix((np.ones(len(gene)),(np.arange(len(gene)),gene.bin)),shape=(len(gene),len(bins)))
        Y=np.asarray((a.T@counts.to_numpy(float)).T);assert np.array_equal(Y.sum(axis=1),meta.input_binned_counts)
        mu,v,W=reference_model(Y[meta.role.eq('normal')],10,arm_starts(bins,genome_build='hg38'))
        reset=chromosome_starts(bins);local_effect=np.full_like(Y,np.nan,dtype=np.float32)
        saved_post=np.load(hd/'rna_state_posterior.npy');saved_effect=np.load(hd/'rna_log_effect.npy')
        fitted=0
        # RNA computation is completed before any DNA truth value is loaded.
        for i in np.flatnonzero(meta.role.eq('tumour')):
            f=fit_sample(Y[i],mu,v,W,reset,BULK_PROFILE)
            assert np.array_equal(f['post'],saved_post[i]) and np.array_equal(f['dose'],saved_effect[i])
            means=np.log(f['E'].astype(np.float32))[:,None]+np.float32(f['centre'])+f['beta']*LOG_DOSE
            emit=np.asarray(pln_logpmf(Y[i].astype(np.float32)[:,None],means,f['sd'].astype(np.float32)[:,None]))
            local=softmax(emit.astype(float)+np.log(STATE_PRIOR.astype(float)),axis=1)
            local_effect[i]=local@(f['beta']*LOG_DOSE);fitted+=1
        np.save(args.out/'fixed_fit_local_effect.npy',local_effect)
        g37=base.subarm.hg19_genes();tab,seg=base.subarm.absolute()
        d=base.profiles(inp,hd,(g37,tab,seg,*base.subarm.hg19_arms()),(panel/'gene_universe.txt').read_text().splitlines())
        with np.load(stage/'window_stages.npz') as z:old={k:z[k] for k in z.files}
        assert np.array_equal(d['tum'],old['samples'])
        for key,k in (('FL','loss'),('NEU','neutral')):assert np.array_equal(d[key],old[k])
        mapped=base.subarm.gene_values(base.subarm.align_rows(local_effect,meta['sample'],d['tum']),gene.set_index('gene').bin,d['genes'])
        mapped=np.where(d['shared_gene_mapping'][:,None],mapped,np.nan)
        local=base.subarm.window_median(mapped,d['window_id'],len(d['warm']))
        common=old['common'];assert np.isfinite(local[common]).all()
        split=pd.read_csv(panel/'split_manifest.tsv',sep='\t').set_index('sample').loc[d['tum']]
        ci=np.flatnonzero(split.split.eq('calibration'));assert len(ci)==24
        methods=dict(hd_lr=old['ratio'],hd_posterior=old['posterior'],fixed_fit_local=local)
        rows=[];thresholds={}
        for method,value in methods.items():
            threshold=base.loss_threshold(value,d['NEU']&common,ci)
            assert threshold==audit.threshold_reference(value,d['NEU']&common,ci)
            thresholds[method]=threshold
            metrics=base.score_loss(value,d['FL'],d['NEU'],common,threshold,negative=d['FG']|d['NEU'])
            for j,sample in enumerate(d['tum']):
                eligible=common[:,j]&d['FL'][:,j];neutral=common[:,j]&d['NEU'][:,j]
                for label,mask in (('recall_loss',eligible),('fpr_loss',neutral)):
                    expected=int((value[mask,j]<threshold).sum())/int(mask.sum()) if mask.any() else np.nan
                    got=metrics[j][label];assert (np.isnan(expected) and np.isnan(got)) or got==expected
                rows.append(dict(sample=sample,split=split.iloc[j].split,method=method,**metrics[j]))
        frame=pd.DataFrame(rows);frame.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        frame.groupby(['split','method']).agg(patients=('sample','size'),loss_patients=('recall_loss','count'),
            mean_recall=('recall_loss','mean'),mean_fpr=('fpr_loss','mean'),mean_precision=('precision_loss','mean')).to_csv(args.out/'summary.tsv',sep='\t')
        np.savez_compressed(args.out/'window_scores.npz',samples=d['tum'].astype(str),common=common,loss=d['FL'],neutral=d['NEU'],**methods)
        reference=pd.read_csv(stage/'patient_stage_scores.tsv',sep='\t',float_precision='round_trip').set_index(['sample','method'])
        for row in frame[frame.method.ne('fixed_fit_local')].itertuples():
            previous=reference.loc[(row.sample,row.method)]
            for key in ('recall_loss','fpr_loss','precision_loss'):
                x=getattr(row,key);y=previous[key];assert (pd.isna(x) and pd.isna(y)) or x==y
        report=dict(exact_tumour_refits=fitted,reference_donors=20,calibration_patients=24,comparison_patients=24,
            thresholds=thresholds,independent_quantiles_exact=True,old_patient_scores_exact=True,
            inference_scope='Conditional fixed-nuisance ablation only; HMM used in fitting E, beta, centre and variance',
            new_inferential_tests=0,method_promotion=False,biological_calibration=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend=backend)
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
