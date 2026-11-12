"""Bypass intermediate bins while preserving the existing pooled-reference effect."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    p=Path(__file__).with_name(name+'.py');s=importlib.util.spec_from_file_location('binning_'+name,p)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def median_bins(values,bins):
    values,bins=np.asarray(values),np.asarray(bins)
    if values.ndim!=2 or bins.shape!=(len(values),) or not np.array_equal(np.unique(bins),np.arange(bins.max()+1)):
        raise ValueError('Invalid aligned bin membership')
    return np.array([np.median(values[bins==b],axis=0) for b in range(bins.max()+1)])


def checks():
    values=np.array([0.,1.,100.,2.,3.,4.])[:,None];bins=np.array([0,0,0,1,1,1])
    med=median_bins(values,bins)
    if np.median(values)!=2.5 or np.median(med[bins])!=2.:raise AssertionError('Noncommutation counterexample failed')
    perm=np.array([2,4,0,3,5,1]);np.testing.assert_array_equal(median_bins(values[perm],bins[perm]),med)
    np.testing.assert_array_equal(median_bins(values+7,bins),med+7)
    counts=np.array([[10.,5.,20.],[30.,15.,25.],[5.,20.,10.],[15.,20.,5.]])
    normal=np.array([True,True,False]);bins=np.array([0,0,1,1])
    gene,pooled=gene_effects(counts,bins,normal)
    depths=[math.fsum(counts[:,s]) for s in range(3)]
    rates=[math.fsum(counts[g,:2])/math.fsum(depths[:2]) for g in range(4)]
    scalar=np.array([[math.log((counts[g,s]+.5)/(rates[g]*depths[s]+.5)) for s in range(3)] for g in range(4)])
    center=np.array([np.median([np.median(scalar[bins==b,s]) for b in range(2)]) for s in range(3)])
    allowance=32*np.finfo(float).eps*max(1.,float(np.max(np.abs(scalar))))
    np.testing.assert_allclose(gene,scalar-center,rtol=0,atol=allowance)
    np.testing.assert_allclose(pooled,median_bins(scalar,bins)-center,rtol=0,atol=allowance)
    return dict(independent_noncommutation_counterexample=True,within_bin_permutation=True,
                translation_invariance=True,independent_scalar_gene_effect=True,roundoff_allowance=allowance)


def gene_effects(counts,bins,normal):
    counts,bins,normal=np.asarray(counts,float),np.asarray(bins,int),np.asarray(normal,bool)
    depth=counts.sum(axis=0)
    if np.any(depth<=0) or not np.isfinite(counts).all() or np.any(counts<0):raise ValueError('Invalid counts')
    mu=counts[:,normal].sum(axis=1)/depth[normal].sum()
    if normal.sum()<2 or np.any(mu<=0):raise ValueError('Unsupported reference')
    gene=np.log((counts+.5)/(mu[:,None]*depth[None]+.5))
    pooled=median_bins(gene,bins);center=np.median(pooled,axis=0)
    return gene-center,pooled-center


def direct_gene(directory,d,subarm,programme_transfer=False):
    X=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0)
    X.index=subarm.canonical_symbols(X.index,'count genes')
    meta=pd.read_csv(directory/'samples.tsv',sep='\t').set_index('sample').loc[X.columns]
    mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t')
    mapping.gene=subarm.canonical_symbols(mapping.gene,'mapped genes')
    subarm.unique_ids(X.columns,'count samples');subarm.unique_ids(mapping.gene,'mapped genes')
    normal=meta.role.eq('normal').to_numpy()
    G=X.loc[mapping.gene].to_numpy(float);bin_ids=mapping.bin.to_numpy(int)
    gene,pooled=gene_effects(G,bin_ids,normal)
    # The previous gene-half-count median path used the same mapped-gene depth,
    # pooled reference and bin-median center. Verify all of them before bypassing.
    factor=load('aggregation_reference');previous,_=factor.effects(G,bin_ids,normal)
    np.testing.assert_array_equal(pooled.T,previous['pooled_median'])
    target=X.columns.get_indexer(d['tum'])
    if np.any(target<0):raise ValueError('Unknown tumour sample')
    original=subarm.gene_values(pooled[:,target].T,mapping.set_index('gene').bin,d['genes'])
    projected=pd.DataFrame(gene[:,target],index=mapping.gene,columns=d['tum']).reindex(d['genes']).to_numpy(float)
    values_by_gene=[projected,original]
    if programme_transfer:
        settings=json.loads((directory/'haplodose/run.json').read_text())
        if settings['bin_stat']!='sum' or settings['input_contract']['counts_kind']!='raw':
            raise ValueError('Programme transfer requires the frozen raw sum-bin reference')
        # This transfers the observed bin-level RNA-effect correction. It is not
        # a reconstruction of gene expected counts or a gene likelihood model:
        # the bin half-count and gene half-count are different measurements.
        A=sparse.csr_matrix((np.ones(len(G)),(np.arange(len(G)),bin_ids)),shape=(len(G),len(pooled)))
        Y=np.asarray((A.T@G).T);L=Y.sum(axis=1)
        np.testing.assert_array_equal(L,G.sum(axis=0))
        mu=Y[normal].sum(axis=0)/Y[normal].sum()
        raw=np.log((Y+.5)/(L[:,None]*mu[None]+.5));raw-=np.median(raw,axis=1,keepdims=True)
        saved_meta=pd.read_csv(directory/'haplodose/samples.tsv',sep='\t').set_index('sample')
        saved=np.load(directory/'haplodose/log_ratio.npy')
        order=saved_meta.index.get_indexer(X.columns)
        if np.any(order<0) or len(np.unique(order))!=len(order):raise ValueError('Unaligned corrected bin profiles')
        corrected=saved[order].astype(float)
        if corrected.shape!=raw.shape or not np.isfinite(corrected).all():raise ValueError('Invalid corrected profiles')
        delta=(corrected-raw).T
        for effects in (gene+delta[bin_ids],(pooled+delta)[bin_ids]):
            values_by_gene.append(pd.DataFrame(effects[:,target],index=mapping.gene,columns=X.columns[target])
                                  .reindex(d['genes']).to_numpy(float))
        # A common correction preserves within-bin contrasts and the original
        # bin path. Check those identities with an accumulation-scale allowance.
        allowance=16*np.finfo(float).eps*max(1.,float(np.max(np.abs(raw))),float(np.max(np.abs(delta))))
        np.testing.assert_allclose(raw+delta.T,corrected,rtol=0,atol=allowance)
        for b in range(len(pooled)):
            g=np.flatnonzero(bin_ids==b)
            np.testing.assert_allclose((gene[g]+delta[b])-(gene[g[0]]+delta[b]),
                                       gene[g]-gene[g[0]],rtol=0,atol=allowance)
    windows=[]
    for values in values_by_gene:
        values=np.where(d['shared_gene_mapping'][:,None],values,np.nan)
        windows.append(subarm.window_median(values,d['window_id'],len(d['warm'])))
    return windows,dict(genes_in_bins=len(mapping),bins=len(pooled),scoring_genes=len(d['genes']),
                       normal_samples=meta.index[normal].tolist(),sample_order=d['tum'].tolist(),
                       normalization='unchanged mapped-gene depth, pooled normal rate, gene half-count, bin-median center',
                       intervention='direct gene effects to evaluation-window medians; omit only intermediate bin medians',
                       programme_transfer=programme_transfer,
                       transfer_units='additive natural-log RNA effect, not gene expected counts or calibrated uncertainty')


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frozen-panel',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    p.add_argument('--programme-transfer',action='store_true')
    a=p.parse_args();a.out.mkdir(exist_ok=False);base=load('profile');op=load('raw_operator_panel')
    partition=load('profile_operator_control');precision=load('pearson_reference')
    (a.out/'checks.json').write_text(json.dumps(checks(),indent=2))
    root=Path(__file__).parent
    inputs={n:root/n for n in ('gene_aggregation_panel.py','aggregation_reference.py','profile.py','raw_operator_panel.py','profile_operator_control.py','pearson_reference.py','statistics.py')}
    inputs.update(subarm=root.parent/'tcga_bulk/subarm.py',split=a.frozen_panel/'split_manifest.tsv',
                  genes=a.frozen_panel/'gene_universe.txt',support=a.frozen_panel/'raw_support_and_truth.npz',
                  previous_summary=a.previous/'summary.json',previous_scores=a.previous/'patient_scores.tsv',
                  annotation=base.subarm.GTF37,arms19=base.subarm.CYTO_HG19,centromeres38=base.subarm.CEN38,
                  dna_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
                  dna_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for i in range(9):
        directory=a.frozen_panel/f'raw_refit_{i:03d}'
        for name in ('counts.tsv','samples.tsv','positions.tsv','haplodose/samples.tsv','haplodose/genes.tsv.gz',
                     'haplodose/log_ratio.npy','haplodose/segments.tsv'):
            inputs[f'{i}:{name}']=directory/name
        inputs[f'{i}:previous']=a.previous/f'refit_{i:03d}.npz'
    if a.programme_transfer:
        for i in range(9):
            inputs[f'{i}:settings']=a.frozen_panel/f'raw_refit_{i:03d}/haplodose/run.json'
    params=dict(bootstrap=20000,seed=base.SEED,family=93 if a.programme_transfer else 75,
                programme_transfer=a.programme_transfer)
    rec=start_run(a.out,inputs,params,dict(scope='Development/regression-only aggregation bypass; no new evaluation'),mode='gene_aggregation')
    try:
        with np.load(inputs['support'],allow_pickle=True) as z:frozen={k:z[k] for k in z.files}
        genes=inputs['genes'].read_text().splitlines();meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample')
        g37=base.subarm.hg19_genes();tab,seg=base.subarm.absolute();truth=(g37,tab,seg,*base.subarm.hg19_arms());refits=[]
        for i in range(9):
            directory=a.frozen_panel/f'raw_refit_{i:03d}';d=base.profiles(directory,directory/'haplodose',truth,genes)
            np.testing.assert_array_equal(d['tum'],frozen['samples'])
            for name,old in [('D','dna'),('FL','loss'),('NEU','neutral')]:np.testing.assert_array_equal(d[name],frozen[old])
            for name in ('naive','naive_smooth'):d['W'][name]=d['W'][name].astype(float)*np.log(2.)
            with np.load(a.previous/f'refit_{i:03d}.npz',allow_pickle=True) as old:
                np.testing.assert_array_equal(old['samples'],d['tum']);np.testing.assert_array_equal(old['arms'],d['warm'])
                for name in base.METHODS:np.testing.assert_array_equal(d['W'][name],old[name])
                for name in old.files:
                    if name not in ('samples','arms'):d['W'][name]=old[name]
            windows,identity=direct_gene(directory,d,base.subarm,a.programme_transfer)
            raw,control=windows[:2]
            (_,old_partition),_=op.operator(d['W']['naive'],control,d['warm'],partition.partition)
            np.testing.assert_array_equal(old_partition,d['W']['pooled_median_partition'])
            (naive,smoothed),details=op.operator(d['W']['naive'],raw,d['warm'],partition.partition)
            np.testing.assert_array_equal(naive,d['W']['naive_own_partition'])
            if a.programme_transfer:
                np.testing.assert_array_equal(smoothed,d['W']['pooled_gene_partition'])
                for name,values in zip(('programme_gene_partition','programme_median_partition'),windows[2:]):
                    (naive,modified),_=op.operator(d['W']['naive'],values,d['warm'],partition.partition)
                    np.testing.assert_array_equal(naive,d['W']['naive_own_partition'])
                    d['W'][name]=modified
            else:d['W']['pooled_gene_partition']=smoothed
            refits.append(d)
            np.savez_compressed(a.out/f'refit_{i:03d}.npz',samples=d['tum'],arms=d['warm'],**d['W'])
            (a.out/f'operator_{i:03d}.json').write_text(json.dumps(dict(identity=identity,operator=details),indent=2))
            print(json.dumps(dict(refit=i,original_pooled_median_exact=True)),flush=True)
        previous=json.loads(inputs['previous_summary'].read_text())
        if a.programme_transfer:
            added=('programme_gene_partition','programme_median_partition')
            contrasts=(('programme_gene_partition','pooled_gene_partition'),
                ('programme_median_partition','pooled_median_partition'),
                ('programme_gene_partition','programme_median_partition'),
                ('programme_gene_partition','hd_own_partition'),
                ('programme_gene_partition','naive_own_partition'),
                ('programme_median_partition','hd_own_partition'))
        else:
            added=('pooled_gene_partition',)
            contrasts=(('pooled_gene_partition','pooled_median_partition'),('pooled_gene_partition','naive_binned_partition'),
                ('pooled_gene_partition','naive_own_partition'),('pooled_gene_partition','hd_own_partition'))
        base.METHODS=tuple(previous['methods'])+added
        base.CONTRASTS=tuple(tuple(x) for x in previous['contrasts'])+contrasts
        if len(base.METHODS)!=(17 if a.programme_transfer else 15) or len(base.CONTRASTS)*3!=params['family']:
            raise AssertionError('Fixed comparator family changed')
        ordinary=base.correlations
        def refined(d,common):
            result=ordinary(d,common)
            for name in base.METHODS:
                for j in range(len(d['tum'])):
                    x,y=d['W'][name][common[:,j],j],d['D'][common[:,j],j]
                    if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:result[name]['pearson_window'][j]=precision.reference(x,y)['value']
            return result
        base.correlations=refined
        split=meta.reindex(refits[0]['tum']);ci=np.flatnonzero(split.split.eq('calibration'));ei=np.flatnonzero(split.split.eq('evaluation'))
        if len(ci)!=24 or len(ei)!=24 or split.patient.duplicated().any():raise AssertionError('Patient roles changed')
        _,thresholds,records,common=base.score_profiles(refits[0],ci);np.testing.assert_array_equal(common,frozen['common'])
        patient=pd.DataFrame(records);patient['split']=patient['sample'].map(meta.split)
        patient.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
        comparisons=base.nested_comparisons(refits[0],refits[1:],ci,ei,20000,base.SEED)
        pd.DataFrame(comparisons).to_csv(a.out/'comparisons.tsv',sep='\t',index=False)
        (a.out/'summary.json').write_text(json.dumps(dict(params=params,methods=base.METHODS,contrasts=base.CONTRASTS,
            thresholds=thresholds,comparisons=comparisons,
            limits='Same exposed development cohort; finite8reference ensemble and bootstrap tails do not certify interval coverage'),indent=2))
        for r in rec['input']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Frozen input changed')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
