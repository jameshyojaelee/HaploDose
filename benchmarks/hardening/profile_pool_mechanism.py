"""Reconstruct primary generating RNA means and isolate aggregation effects."""
import argparse
import importlib.util
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, spearmanr
from haplodose import simulate
from haplodose.profile import ProfileParams
from haplodose.provenance import start_run, finish_run, file_record


def centre(x):
    return x-np.median(x)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    prior=a.task.parent/'profile_pool_export_20261109/parity_v1';fit=prior/'fitted_new';data=prior/'data'
    op_path=Path(__file__).with_name('profile_operator_control.py')
    spec=importlib.util.spec_from_file_location('frozen_pool_operator',op_path)
    op=importlib.util.module_from_spec(spec);spec.loader.exec_module(op)
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',operator=op_path,
        simulator=Path(simulate.__file__),prior_manifest=prior/'run_manifest.json',
        counts=data/'counts.h5',input_cells=data/'cells.tsv',positions=data/'positions.tsv',
        cells=fit/'cells.tsv.gz',profiles=fit/'profiles.tsv.gz',clones=fit/'clones.tsv',
        genes=fit/'genes.tsv.gz',bins=fit/'bins.tsv',binned=fit/'bin_counts.npy',capture=prior/'fitted_new_capture.npz')
    previous=json.loads(inputs['prior_manifest'].read_text());assert previous['status']=='complete'
    for record in previous['output']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(a.out,inputs,dict(seed=20261109,specimens=1,new_fits=0),
        dict(scope='Conditional generating RNA target; oracle reference identity, no genomic/coverage claim'),
        mode='profile_pool_mechanism')
    try:
        observed=[]
        def trace(frame,event,arg):
            if event=='return' and frame.f_code.co_name=='draw' and frame.f_code.co_filename==simulate.__file__:
                v=frame.f_locals
                observed.append(dict(ids=list(v['ids']),mean=np.array(v['L'][:,None]*v['share'],copy=True),
                                     depth=np.array(v['L'],copy=True)))
        assert sys.getprofile() is None
        sys.setprofile(trace)
        try:
            X,genes,meta,positions,events,_=simulate.simulate(seed=20261109,n_chrom=4,genes_per_chrom=60,
                extra_genes=(),n_types=2,ref_cells_per_type=60,depth=1200.,
                specimens=[dict(name='A',n_cells=160,normal_frac=.1,
                    clones=[dict(frac=.5,events=[('chr3',0.,1.,3)]),dict(frac=.5,events=[('chr3',0.,1.,4)])])])
        finally:sys.setprofile(None)
        with h5py.File(inputs['counts']) as h:
            old=sparse.csr_matrix((h['X/data'][:],h['X/indices'][:],h['X/indptr'][:]),shape=X.shape)
            assert np.array_equal(h['genes'].asstr()[:],genes)
            assert np.array_equal(h['cell_ids'].asstr()[:],meta.cell_id)
        assert (old!=X).nnz==0
        M=np.vstack([v['mean'] for v in observed])
        assert np.array_equal(np.concatenate([v['ids'] for v in observed]),meta.cell_id)
        assert M.shape==X.shape and (M>0).all()
        cells=pd.read_csv(inputs['cells'],sep='\t');used=pd.read_csv(inputs['genes'],sep='\t')
        bins=pd.read_csv(inputs['bins'],sep='\t');profiles=pd.read_csv(inputs['profiles'],sep='\t',float_precision='round_trip')
        clones=pd.read_csv(inputs['clones'],sep='\t',float_precision='round_trip').sort_values('column')
        assert cells.cell_id.tolist()==meta.cell_id.tolist()
        gi=pd.Index(genes).get_indexer(used.gene);assert (gi>=0).all()
        A=sparse.csr_matrix((np.ones(len(used)),(np.arange(len(used)),used.bin)),shape=(len(used),len(bins)))
        Y=np.asarray((X[:,gi]@A).toarray());mu=np.asarray(M[:,gi]@A)
        assert np.array_equal(Y,np.load(inputs['binned']))
        q=cells.specimen.eq('A')&cells.status.eq('assigned')
        ref=meta.role.eq('reference')&meta.ref_split.eq('fit')&meta.ref_group.eq('type0')
        assert int(q.sum())==160 and int(ref.sum())==30
        y=Y[q].sum(axis=0);mean=mu[q].sum(axis=0);reference=Y[ref].sum(axis=0)
        reference_mean=mu[ref].sum(axis=0)
        truth=centre(np.log((mean/mean.sum())/(reference_mean/reference_mean.sum())))
        profile_y=profiles.pivot(index='bin',columns='clone',values='pooled_count')[clones.clone].to_numpy()
        profile_E=profiles.pivot(index='bin',columns='clone',values='expected_count')[clones.clone].to_numpy()
        effects=profiles.pivot(index='bin',columns='clone',values='log_ratio')[clones.clone].to_numpy()
        with np.load(inputs['capture']) as z:R=z['A__R']
        u=np.finfo(float).eps;operations=2*len(R)+2;gamma=operations*u/(1-operations*u)
        assert (np.abs(profile_y.sum(axis=1)-y)<=gamma*np.maximum(1,y)).all()
        n=cells.loc[q,'clone'].value_counts().reindex(clones.clone,fill_value=0).to_numpy()
        methods=dict(map_mean_log=centre(effects@n/n.sum()),
            responsibility_mean_log=centre(effects@R.sum(axis=0)/R.sum()),
            count_mass_mean_log=centre(effects@profile_y.sum(axis=0)/profile_y.sum()),
            pooled_model=centre(np.log((profile_y.sum(axis=1)+.5)/(profile_E.sum(axis=1)+.5))),
            naive_known_reference=centre(np.log((y+.5)/(reference/reference.sum()*y.sum()+.5))),
            oracle_reference_rate=centre(np.log((y+.5)/(reference_mean/reference_mean.sum()*y.sum()+.5))))
        reset=np.r_[True,bins.chromosome.to_numpy()[1:]!=bins.chromosome.to_numpy()[:-1]]
        blocks=np.split(np.arange(len(bins)),np.flatnonzero(reset)[1:])
        differences=np.concatenate([np.diff(methods['naive_known_reference'][b]) for b in blocks])
        sigma=np.median(np.abs(differences-np.median(differences)))/(np.sqrt(2)*norm.ppf(.75))
        assert sigma>0
        penalty=3*sigma*sigma*np.log(len(bins));minimum=ProfileParams().min_bins
        regions=[]
        for b in blocks:
            _,pieces,_=op.partition(methods['naive_known_reference'][b],penalty,minimum)
            regions.extend([[int(b[start]),int(b[end-1])+1] for start,end in pieces])
        metrics=[];stored=dict(truth=truth,observed=y,query_mean=mean,reference_count=reference,
            reference_mean=reference_mean,model_expected=profile_E.sum(axis=1),R=R,gene_mean=M,
            used_gene_mean=mu,query_mask=q.to_numpy(),reference_mask=ref.to_numpy(),gene_columns=gi)
        for method,value in methods.items():
            fitted=np.empty_like(value)
            for start,end in regions:fitted[start:end]=value[start:end].mean()
            for operator,v in [('none',value),('shared_naive_boundaries',fitted)]:
                metrics.append(dict(method=method,operator=operator,MSE=float(np.mean((v-truth)**2)),
                    Pearson=float(np.corrcoef(v,truth)[0,1]),Spearman=float(spearmanr(v,truth).statistic)))
                stored[method+'__'+operator]=v
        # Fragmentation cannot change a sum-based aggregate; per-component offsets can.
        fragmented_y=np.repeat(profile_y/2,2,axis=1);fragmented_E=np.repeat(profile_E/2,2,axis=1)
        fragmented=np.log((fragmented_y.sum(axis=1)+.5)/(fragmented_E.sum(axis=1)+.5))
        original=np.log((profile_y.sum(axis=1)+.5)/(profile_E.sum(axis=1)+.5))
        budget=32*u*(1+np.abs(original));assert (np.abs(fragmented-original)<=budget).all()
        # Independently calculable unequal-exposure example, no simulation/fitting.
        scalar=math.log((1*1+9*2)/(1+9));cell_mean=.5*math.log(2)
        assert scalar>cell_mean
        np.savez_compressed(a.out/'primary_means_and_profiles.npz',**stored)
        pd.DataFrame(metrics).to_csv(a.out/'metrics.tsv',sep='\t',index=False)
        (a.out/'regions.json').write_text(json.dumps(regions)+'\n')
        (a.out/'report.json').write_text(json.dumps(dict(count_matrix_exact=True,cells=int(X.shape[0]),
            query_cells=int(q.sum()),reference_cells=int(ref.sum()),retained_bins=len(bins),
            biological_specimens=0,draws_traced=len(observed),shared_sigma=float(sigma),shared_penalty=float(penalty),
            reference_identity='Known generating type0; not independently established from RNA',
            truth='Generating pooled relative RNA means, not DNA copy numbers',
            unequal_exposure_example=dict(pooled_log_effect=scalar,cell_mean_log_effect=cell_mean),
            max_fragmentation_discrepancy=float(np.max(np.abs(fragmented-original))),metrics=metrics,
            no_default_or_scientific_criteria_change=True),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU NumPy/SciPy, no new fit')


if __name__=='__main__':main()
