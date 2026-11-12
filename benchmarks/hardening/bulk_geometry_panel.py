"""Matched transition-geometry ablation on the exposed raw-STAR bulk panel."""
import argparse
import importlib.util
import json
from pathlib import Path

import jax.numpy as jnp
import numpy as np
import pandas as pd
from scipy import sparse

from haplodose.bins import arm_starts
from haplodose.bulk import BULK_PROFILE, fit_sample, reference_model
from haplodose.device import select
from haplodose.kernels.hmm import forward_backward, sticky_transition
from haplodose.kernels.likelihood import pln_logpmf
from haplodose.profile import LOG_DOSE, STATE_PRIOR
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec=importlib.util.spec_from_file_location('geometry_panel_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parents[1];stage=root/'dev/bulk_posterior_export_20261109/stage_v1'
    hd=stage/'haplodose';panel=root/'dev/hardening_20261106_2119/profile_v6';inp=panel/'raw_refit_000'
    operator=root/'dev/joint_confirmation_20261108/raw_operator_v2/refit_000.npz'
    base=load('profile');geometry=load('genomic_geometry');ranks=load('bulk_posterior_stage_audit');pearson=load('pearson_reference')
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',geometry=Path(geometry.__file__),profile=Path(base.__file__),
        ranks=Path(ranks.__file__),pearson=Path(pearson.__file__),scorer=Path(base.subarm.__file__),
        statistics=Path(__file__).with_name('statistics.py'),operator=operator,
        preflight=args.task/'preflight_v1/checks.json',split=panel/'split_manifest.tsv',genes=panel/'gene_universe.txt',
        dna_gtf=base.subarm.GTF37,cytoband=base.subarm.CYTO_HG19,centromeres=base.subarm.CEN38,
        absolute_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        absolute_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for directory,prefix,names in ((inp,'input',('counts.tsv','samples.tsv','positions.tsv')),
        (stage,'stage',('window_stages.npz','patient_stage_scores.tsv')),
        (hd,'hd',('bins.tsv','genes.tsv.gz','samples.tsv','rna_log_effect.npy','rna_state_posterior.npy'))):
        inputs.update({prefix+':'+name:directory/name for name in names})
    run=start_run(args.out,inputs,dict(permutation_seed=20261109,parameter_tuning=False,nuisance_refit_for_geometry=False),
        dict(scope='Conditional exposed-panel geometry ablation; no independent evaluation or biological calibration'),mode='bulk_geometry_panel')
    try:
        backend=select('cpu');assert backend=='cpu'
        assert json.loads(inputs['preflight'].read_text())['enumeration_paths']==256
        bins=pd.read_csv(hd/'bins.tsv',sep='\t');genes=pd.read_csv(hd/'genes.tsv.gz',sep='\t')
        meta=pd.read_csv(hd/'samples.tsv',sep='\t',float_precision='round_trip')
        counts=pd.read_csv(inp/'counts.tsv',sep='\t',index_col=0).loc[genes.gene,meta['sample']]
        aggregation=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),genes.bin)),shape=(len(genes),len(bins)))
        Y=np.asarray((aggregation.T@counts.to_numpy(float)).T);assert np.array_equal(Y.sum(axis=1),meta.input_binned_counts)
        g=geometry.geometry(bins,genes);reset=g['reset'];J=len(bins)
        transition,generator,prior=geometry.transition_model(BULK_PROFILE.p_change)
        legacy=np.asarray(sticky_transition(4,BULK_PROFILE.p_change,STATE_PRIOR))
        logs=dict(hd_posterior=legacy,
            ideal_constant=np.broadcast_to(np.log(transition).astype(np.float32),(J,4,4)),
            distance=geometry.log_transitions(g['ratios'],generator)[0],
            distance_shuffled=geometry.log_transitions(g['shuffled_ratios'],generator)[0])
        np.savez_compressed(args.out/'geometry.npz',**{k:g[k] for k in
            ('reset','anchors','distance_bp','ratios','shuffled_ratios','permutation')},transition=transition,generator=generator,**logs)
        (args.out/'blocks.json').write_text(json.dumps(g['blocks'],indent=2)+'\n')
        pd.DataFrame(dict(bin=np.arange(J),chromosome=bins.chromosome,block=bins.genomic_block,
            anchor_start=g['anchors'],distance_bp=g['distance_bp'],process_time=g['ratios'],
            shuffled_time=g['shuffled_ratios'],shuffled_source_bin=g['permutation'])).to_csv(args.out/'geometry.tsv',sep='\t',index=False)
        ref=reference_model(Y[meta.role.eq('normal')],10,arm_starts(bins,genome_build='hg38'))
        original_post=np.load(hd/'rna_state_posterior.npy');original_effect=np.load(hd/'rna_log_effect.npy')
        effects={name:np.full_like(Y,np.nan,dtype=np.float32) for name in logs}
        reference_effect={name:np.full_like(Y,np.nan,dtype=float) for name in logs};numerics=[]
        fitted=0
        for i in np.flatnonzero(meta.role.eq('tumour')):
            sample=str(meta.iloc[i]['sample']);fit=fit_sample(Y[i],*ref,reset,BULK_PROFILE)
            assert np.array_equal(fit['post'],original_post[i]) and np.array_equal(fit['dose'],original_effect[i])
            means=np.log(fit['E'].astype(np.float32))[:,None]+np.float32(fit['centre'])+fit['beta']*LOG_DOSE
            emit=np.asarray(pln_logpmf(Y[i].astype(np.float32)[:,None],means,fit['sd'].astype(np.float32)[:,None]))
            for name,lt in logs.items():
                post=np.asarray(forward_backward(jnp.asarray(emit),jnp.asarray(lt),jnp.log(jnp.asarray(STATE_PRIOR)),jnp.asarray(reset))[0])
                reference_post=geometry.reference_posterior(emit,lt,reset)
                error=float(np.max(np.abs(post-reference_post)));assert error<1e-4,(sample,name,error)
                effects[name][i]=post@(fit['beta']*LOG_DOSE)
                reference_effect[name][i]=reference_post@(fit['beta']*LOG_DOSE)
                if name=='hd_posterior':assert np.array_equal(post,original_post[i]) and np.array_equal(effects[name][i],original_effect[i])
                numerics.append(dict(sample=sample,method=name,maximum_posterior_discrepancy=error,
                    maximum_effect_discrepancy=float(np.max(np.abs(effects[name][i]-reference_effect[name][i])))))
            np.savez_compressed(args.out/(sample+'.npz'),emission=emit,beta=fit['beta'],centre=fit['centre'],tau=fit['tau'])
            fitted+=1
        pd.DataFrame(numerics).to_csv(args.out/'numerical_references.tsv',sep='\t',index=False)
        # No DNA value enters transition construction, fitting or posterior projection.
        g37=base.subarm.hg19_genes();tab,seg=base.subarm.absolute()
        d=base.profiles(inp,hd,(g37,tab,seg,*base.subarm.hg19_arms()),(panel/'gene_universe.txt').read_text().splitlines())
        with np.load(stage/'window_stages.npz') as z:old={k:z[k] for k in z.files}
        assert np.array_equal(d['tum'],old['samples'])
        for key,k in (('FL','loss'),('NEU','neutral')):assert np.array_equal(d[key],old[k])
        common=np.isfinite(d['D'])&np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS])
        assert np.array_equal(common,old['common'])
        methods={name:d['W'][name] for name in base.METHODS}
        reference_windows={}
        mapping=genes.set_index('gene').bin
        for name in logs:
            for collection,target in ((effects,methods),(reference_effect,reference_windows)):
                matrix=base.subarm.align_rows(collection[name],meta['sample'],d['tum'])
                mapped=base.subarm.gene_values(matrix,mapping,d['genes'])
                mapped=np.where(d['shared_gene_mapping'][:,None],mapped,np.nan)
                target[name]=base.subarm.window_median(mapped,d['window_id'],len(d['warm']))
                assert np.isfinite(target[name][common]).all()
        assert np.array_equal(methods['hd_posterior'],old['posterior'],equal_nan=True)
        # Trusted local frozen NumPy archive contains the historical object sample axis.
        with np.load(operator,allow_pickle=True) as z:
            assert np.array_equal(z['samples'],d['tum'])
            methods['naive_own_partition']=z['naive_own_partition']
        assert np.isfinite(methods['naive_own_partition'][common]).all()
        split=pd.read_csv(panel/'split_manifest.tsv',sep='\t').set_index('sample').loc[d['tum']]
        ci=np.flatnonzero(split.split.eq('calibration'));assert len(ci)==24
        thresholds={};rows=[];pearson.verify()
        for name,value in methods.items():
            threshold=base.loss_threshold(value,d['NEU']&common,ci)
            assert threshold==ranks.threshold_reference(value,d['NEU']&common,ci);thresholds[name]=threshold
            loss=base.score_loss(value,d['FL'],d['NEU'],common,threshold,negative=d['FG']|d['NEU'])
            for j,sample in enumerate(d['tum']):
                mask=common[:,j];x,y=value[mask,j],d['D'][mask,j]
                correlation=pearson.reference(x,y)['value'] if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                rows.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=correlation,**loss[j]))
        scores=pd.DataFrame(rows);scores.to_csv(args.out/'patient_scores.tsv',sep='\t',index=False)
        scores.groupby(['split','method']).agg(patients=('sample','size'),loss_patients=('recall_loss','count'),
            mean_recall=('recall_loss','mean'),mean_fpr=('fpr_loss','mean'),mean_precision=('precision_loss','mean'),
            mean_pearson=('pearson','mean')).to_csv(args.out/'summary.tsv',sep='\t')
        np.savez_compressed(args.out/'window_scores.npz',samples=d['tum'].astype(str),common=common,loss=d['FL'],
            neutral=d['NEU'],gain=d['FG'],dna=d['D'],split=split.split.to_numpy(str),**methods,
            **{'reference_'+name:value for name,value in reference_windows.items()})
        old_scores=pd.read_csv(stage/'patient_stage_scores.tsv',sep='\t',float_precision='round_trip').set_index(['sample','method'])
        for row in scores[scores.method.isin([*base.METHODS,'hd_posterior'])].itertuples():
            old_row=old_scores.loc[(row.sample,row.method)]
            for key in ('recall_loss','fpr_loss','precision_loss'):
                a,b=getattr(row,key),old_row[key];assert (pd.isna(a) and pd.isna(b)) or a==b
        (args.out/'report.json').write_text(json.dumps(dict(exact_tumour_refits=fitted,reference_donors=20,
            thresholds=thresholds,methods=list(methods),conditional_nuisance_fits_unchanged=True,
            old_patient_scores_exact=True,independent_rank_thresholds_exact=True,new_tests=0,
            score_units={name:('log2_relative_RNA_rate' if name in ('naive','naive_smooth') else
                'natural_log_relative_RNA_rate') for name in methods},
            scope='Block-matched process time; not a universal per-bp breakpoint model',method_promoted=False),indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend=backend)
    except Exception as error:finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
