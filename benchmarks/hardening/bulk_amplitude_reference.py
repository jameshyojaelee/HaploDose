"""Frozen continuous-amplitude estimators under existing donor-reference draws."""
import argparse
import json
from pathlib import Path

import jax.numpy as jnp
import numpy as np
import pandas as pd
from scipy import sparse

from haplodose.bins import arm_starts
from haplodose.bulk import BULK_PROFILE,fit_sample,reference_model
from haplodose.device import select
from haplodose.kernels.hmm import forward_backward,sticky_transition
from haplodose.kernels.likelihood import pln_logpmf
from haplodose.profile import LOG_DOSE,STATE_PRIOR,log_ratio
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    import importlib.util
    spec=importlib.util.spec_from_file_location('amplitude_reference_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--refit',type=int,choices=range(9),required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    root=a.task.parents[1];panel=root/'dev/hardening_20261106_2119/profile_v6'
    inp=panel/f'raw_refit_{a.refit:03d}';hd=inp/'haplodose'
    previous=root/'dev/bulk_continuous_amplitude_20261109/panel_v2'
    operator=root/f'dev/joint_confirmation_20261108/raw_operator_v2/refit_{a.refit:03d}.npz'
    base=load('profile');geo=load('genomic_geometry');amp=load('bulk_continuous_amplitude')
    ranks=load('bulk_posterior_stage_audit');pearson=load('pearson_reference')
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',amplitudes=Path(amp.__file__),geometry=Path(geo.__file__),
        profile=Path(base.__file__),scorer=Path(base.subarm.__file__),ranks=Path(ranks.__file__),pearson=Path(pearson.__file__),
        operator=operator,split=panel/'split_manifest.tsv',genes=panel/'gene_universe.txt',
        original_windows=previous/'window_scores.npz',original_scores=previous/'patient_scores.tsv',
        dna_gtf=base.subarm.GTF37,cytoband=base.subarm.CYTO_HG19,centromeres=base.subarm.CEN38,
        absolute_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        absolute_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for directory,prefix,names in ((inp,'input',('counts.tsv','samples.tsv','positions.tsv')),
        (hd,'frozen',('bins.tsv','genes.tsv.gz','samples.tsv','segments.tsv','log_ratio.npy'))):
        inputs.update({prefix+':'+name:directory/name for name in names})
    rec=start_run(a.out,inputs,dict(reference_refit=a.refit,parameter_tuning=False,permutation_seed=20261109),
        dict(scope='Exposed-patient reference sensitivity; no independent validation or coverage claim'),mode='bulk_amplitude_reference')
    try:
        backend=select('cpu');assert backend=='cpu'
        bins=pd.read_csv(hd/'bins.tsv',sep='\t');genes=pd.read_csv(hd/'genes.tsv.gz',sep='\t')
        meta=pd.read_csv(hd/'samples.tsv',sep='\t',float_precision='round_trip')
        counts=pd.read_csv(inp/'counts.tsv',sep='\t',index_col=0).loc[genes.gene,meta['sample']]
        aggregation=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),genes.bin)),shape=(len(genes),len(bins)))
        Y=np.asarray((aggregation.T@counts.to_numpy(float)).T)
        assert np.array_equal(Y.sum(axis=1),meta.binned_counts)
        assert np.isfinite(Y).all() and np.all(Y>=0) and np.array_equal(Y,np.rint(Y))
        g=geo.geometry(bins,genes);reset=g['reset'];J=len(bins)
        transition,generator,_=geo.transition_model(BULK_PROFILE.p_change)
        logs=dict(ideal_constant=np.broadcast_to(np.log(transition).astype(np.float32),(J,4,4)),
            distance=geo.log_transitions(g['ratios'],generator)[0],
            distance_shuffled=geo.log_transitions(g['shuffled_ratios'],generator)[0])
        ref=reference_model(Y[meta.role.eq('normal')],10,arm_starts(bins,genome_build='hg38'))
        effects={name:np.full_like(Y,np.nan,dtype=np.float32 if name in logs or name=='hd_posterior' else float)
                 for name in ('hd_posterior',*logs,'hd_path_mean','hd_path_precision_mean')}
        frozen_ratio=np.load(hd/'log_ratio.npy');old_segments=pd.read_csv(hd/'segments.tsv',sep='\t',float_precision='round_trip')
        rows=[];numerics=[]
        for i in np.flatnonzero(meta.role.eq('tumour')):
            sample=str(meta.iloc[i]['sample']);fit=fit_sample(Y[i],*ref,reset,BULK_PROFILE)
            ratio=log_ratio(Y[i],fit['E']);ratio-=np.median(ratio)
            assert np.array_equal(ratio.astype(np.float32),frozen_ratio[i]),sample
            for key in ('centre','beta','tau'):assert fit[key]==meta.iloc[i][key],(sample,key)
            track=log_ratio(Y[i],fit['E'])-fit['centre'];variance=fit['sd']**2+1/(Y[i]+.5)
            ordinary,weighted,parts=amp.amplitudes(track,variance,fit['path'],reset)
            old=old_segments[old_segments['sample'].eq(sample)].reset_index(drop=True)
            assert len(parts)==len(old),(sample,'segment_count')
            for part,segment in zip(parts,old.itertuples()):
                lo,hi=part['first_bin'],part['last_bin']
                assert (lo,hi,part['state']+1)==(segment.first_bin,segment.last_bin,segment.cn_state)
                assert float(fit['dose'][lo:hi+1].mean())==segment.rna_log_effect
                rows.append(dict(sample=sample,**part))
            effects['hd_posterior'][i]=fit['dose'];effects['hd_path_mean'][i]=ordinary;effects['hd_path_precision_mean'][i]=weighted
            means=np.log(fit['E'].astype(np.float32))[:,None]+np.float32(fit['centre'])+fit['beta']*LOG_DOSE
            emit=np.asarray(pln_logpmf(Y[i].astype(np.float32)[:,None],means,fit['sd'].astype(np.float32)[:,None]))
            for name,lt in logs.items():
                post=np.asarray(forward_backward(jnp.asarray(emit),jnp.asarray(lt),jnp.log(jnp.asarray(STATE_PRIOR)),jnp.asarray(reset))[0])
                independent=geo.reference_posterior(emit,lt,reset)
                error=float(np.max(np.abs(post-independent)))
                assert error<1e-4,(sample,name,error)  # Existing independently studied HMM criterion, unchanged.
                effects[name][i]=post@(fit['beta']*LOG_DOSE)
                numerics.append(dict(sample=sample,method=name,maximum_posterior_discrepancy=error))
            np.savez_compressed(a.out/(sample+'.npz'),Y=Y[i],E=fit['E'],centre=fit['centre'],sd=fit['sd'],
                track=track,variance=variance,path=fit['path'],reset=reset,mean=ordinary,precision_mean=weighted)
        pd.DataFrame(rows).to_csv(a.out/'segments.tsv',sep='\t',index=False)
        pd.DataFrame(numerics).to_csv(a.out/'numerical_references.tsv',sep='\t',index=False)
        # All fitted RNA predictions precede loading DNA outcomes.
        truth=base.subarm.hg19_genes();calls,segments=base.subarm.absolute()
        d=base.profiles(inp,hd,(truth,calls,segments,*base.subarm.hg19_arms()),(panel/'gene_universe.txt').read_text().splitlines())
        with np.load(previous/'window_scores.npz') as z:frozen={k:z[k] for k in z.files}
        assert np.array_equal(d['tum'],frozen['samples'])
        common=np.isfinite(d['D'])&np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS])
        for new,old in ((common,'common'),(d['D'],'dna'),(d['FL'],'loss'),(d['NEU'],'neutral'),(d['FG'],'gain')):
            assert np.array_equal(new,frozen[old],equal_nan=True)
        methods={name:d['W'][name] for name in base.METHODS};mapping=genes.set_index('gene').bin
        pd.DataFrame(dict(gene=d['genes'],bin=mapping.reindex(d['genes']).to_numpy(),window=d['window_id'],
                          shared=d['shared_gene_mapping'])).to_csv(a.out/'gene_map.tsv',sep='\t',index=False)
        for name,values in effects.items():
            mapped=base.subarm.gene_values(base.subarm.align_rows(values,meta['sample'],d['tum']),mapping,d['genes'])
            methods[name]=base.subarm.window_median(np.where(d['shared_gene_mapping'][:,None],mapped,np.nan),d['window_id'],len(d['warm']))
            assert np.isfinite(methods[name][common]).all()
        with np.load(operator,allow_pickle=True) as z:
            assert np.array_equal(z['samples'],d['tum']);methods['naive_own_partition']=z['naive_own_partition']
        assert len(methods)==12
        if a.refit==0:
            for name,values in methods.items():assert np.array_equal(values,frozen[name],equal_nan=True),name
        split=pd.read_csv(panel/'split_manifest.tsv',sep='\t').set_index('sample').loc[d['tum']]
        assert not split.patient.duplicated().any() and np.array_equal(split.split.to_numpy(str),frozen['split'])
        ci=np.flatnonzero(split.split.eq('calibration'));assert len(ci)==24
        thresholds={};scores=[];pearson.verify()
        for name,values in methods.items():
            threshold=base.loss_threshold(values,d['NEU']&common,ci)
            assert threshold==ranks.threshold_reference(values,d['NEU']&common,ci);thresholds[name]=threshold
            loss=base.score_loss(values,d['FL'],d['NEU'],common,threshold,negative=d['FG']|d['NEU'])
            for j,sample in enumerate(d['tum']):
                x,y=values[common[:,j],j],d['D'][common[:,j],j]
                correlation=pearson.reference(x,y)['value'] if len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0 else np.nan
                scores.append(dict(sample=str(sample),split=split.iloc[j].split,method=name,pearson=correlation,**loss[j]))
        scores=pd.DataFrame(scores);scores.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
        if a.refit==0:
            old=pd.read_csv(previous/'patient_scores.tsv',sep='\t',float_precision='round_trip').set_index(['sample','method'])
            for row in scores.itertuples():
                for key in ('recall_loss','fpr_loss','precision_loss','pearson'):
                    x,y=getattr(row,key),old.loc[(row.sample,row.method),key]
                    assert (pd.isna(x) and pd.isna(y)) or x==y,(row.sample,row.method,key)
        np.savez_compressed(a.out/'window_scores.npz',samples=d['tum'].astype(str),common=common,loss=d['FL'],neutral=d['NEU'],
                            gain=d['FG'],dna=d['D'],split=split.split.to_numpy(str),**methods)
        (a.out/'report.json').write_text(json.dumps(dict(reference_refit=a.refit,methods=list(methods),thresholds=thresholds,
            exact_frozen_tumour_refits=48,reference_draws_are_not_independent_cohorts=True,
            original_reference_exact=a.refit==0,new_hypothesis_tests=0,method_promoted=False),indent=2)+'\n')
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(a.out,rec,backend=backend)
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
