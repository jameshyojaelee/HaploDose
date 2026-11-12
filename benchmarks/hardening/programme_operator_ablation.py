"""Remove only the bulk programme mean correction on frozen raw-count profiles.

This is a development mechanism test, not new evaluation or automatic promotion.
It reconstructs the no-programme RNA ratio algebraically; it does not refit HMM
states or claim that its corrected-profile SE has been calibrated.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.ndimage import median_filter
from haplodose.profile import log_ratio
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    path=Path(__file__).with_name(name+'.py');s=importlib.util.spec_from_file_location('ablation_'+name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def uncorrected(directory,d,subarm):
    """Same binned counts and pooled reference as bulk.fit_sample with empty W."""
    settings=json.loads((directory/'haplodose/run.json').read_text())
    if settings['bin_stat']!='sum' or settings['input_contract']['counts_kind']!='raw':
        raise ValueError('Only frozen raw sum-bin mode is supported')
    counts=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0)
    meta=pd.read_csv(directory/'samples.tsv',sep='\t')
    mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t')
    subarm.unique_ids(counts.index,'count genes');subarm.unique_ids(counts.columns,'count samples')
    subarm.unique_ids(meta['sample'],'sample identities');subarm.unique_ids(mapping.gene,'mapped genes')
    np.testing.assert_array_equal(counts.index.to_numpy()[mapping.column.to_numpy(int)],mapping.gene)
    values=counts.loc[:,meta['sample']].to_numpy(float)[mapping.column.to_numpy(int)]
    if not np.isfinite(values).all() or np.any(values<0) or not np.equal(values,np.floor(values)).all():
        raise ValueError('Raw count contract violated')
    bins=pd.read_csv(directory/'haplodose/bins.tsv',sep='\t')
    if np.any(mapping.bin<0) or mapping.bin.max()+1!=len(bins):raise ValueError('Invalid frozen bin map')
    onehot=sparse.csr_matrix((np.ones(len(mapping)),(np.arange(len(mapping)),mapping.bin)),shape=(len(mapping),len(bins)))
    Y=np.asarray((onehot.T@values).T)
    np.testing.assert_array_equal(Y.sum(axis=1),values.sum(axis=0))
    saved=pd.read_csv(directory/'haplodose/samples.tsv',sep='\t').set_index('sample').reindex(meta['sample'])
    np.testing.assert_array_equal(Y.sum(axis=1),saved.binned_counts)
    normal=meta.role.eq('normal').to_numpy();mu=Y[normal].sum(axis=0)/Y[normal].sum()
    target=meta.set_index('sample').index.get_indexer(d['tum'])
    if np.any(target<0) or normal[target].any():raise ValueError('Only tumour target profiles are reconstructed')
    yt=Y[target];expected=yt.sum(axis=1)[:,None]*mu
    if np.any(expected<=0):raise ValueError('Undefined raw reference ratio')
    r=log_ratio(yt,expected);r-=np.median(r,axis=1,keepdims=True)
    # Independently written scalar algebra for bounded representative entries.
    for sample in (0,len(yt)-1):
        ref=np.array([np.log((float(y)+.5)/(float(e)+.5)) for y,e in zip(yt[sample],expected[sample])])
        ref-=np.median(ref)
        np.testing.assert_allclose(r[sample],ref,rtol=0,atol=4*np.finfo(float).eps)
    by_gene=subarm.gene_values(r,mapping.set_index('gene').bin,d['genes'])
    by_gene=np.where(d['shared_gene_mapping'][:,None],by_gene,np.nan)
    smooth=np.full_like(by_gene,np.nan);chrom=d['gene_positions'].chromosome.to_numpy()
    for c in pd.unique(chrom):
        ix=chrom==c;smooth[ix]=median_filter(by_gene[ix],size=(subarm.SMOOTH,1),mode='nearest')
    n=len(d['warm']);raw=subarm.window_median(by_gene,d['window_id'],n)
    smoothed=subarm.window_median(smooth,d['window_id'],n)
    return raw,smoothed,dict(samples=d['tum'].tolist(),reference_samples=meta.loc[normal,'sample'].tolist(),
                            count_conservation=True,n_bins=len(bins),gene_mapping_identity=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frozen-panel',type=Path,required=True);p.add_argument('--operator-panel',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);base=load('profile');op=load('raw_operator_panel')
    partition=load('profile_operator_control');precision=load('pearson_reference');partition.verify();precision.verify()
    root=Path(__file__).parent;inputs={name:root/name for name in
        ('programme_operator_ablation.py','profile.py','raw_operator_panel.py','profile_operator_control.py','pearson_reference.py','statistics.py')}
    inputs.update(subarm=root.parent/'tcga_bulk/subarm.py',split=a.frozen_panel/'split_manifest.tsv',
        genes=a.frozen_panel/'gene_universe.txt',support=a.frozen_panel/'raw_support_and_truth.npz',
        old_operator_manifest=a.operator_panel/'run_manifest.json',annotation=base.subarm.GTF37,
        arms19=base.subarm.CYTO_HG19,centromeres38=base.subarm.CEN38,
        dna_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        dna_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for i in range(9):
        directory=a.frozen_panel/f'raw_refit_{i:03d}'
        for f in ('counts.tsv','samples.tsv','positions.tsv','haplodose/samples.tsv','haplodose/bins.tsv',
                  'haplodose/genes.tsv.gz','haplodose/log_ratio.npy','haplodose/segments.tsv','haplodose/run.json'):
            inputs[f'{i}:{f}']=directory/f
        inputs[f'{i}:old_operator']=a.operator_panel/f'refit_{i:03d}.npz'
    params=dict(bootstrap=10000,seed=base.SEED,primary_family=36)
    rec=start_run(a.out,inputs,params,dict(role='already_inspected_development',units='natural_log_RNA',
        intervention='remove programme mean correction only; fixed genes, bins, reference and scoring support',
        limitations='No HMM refit, no independently validated interval coverage or biological superiority'),mode='programme_operator_ablation')
    try:
        meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample');genes=inputs['genes'].read_text().splitlines()
        with np.load(inputs['support'],allow_pickle=True) as z:frozen={k:z[k] for k in z.files}
        g37=base.subarm.hg19_genes();tab,seg=base.subarm.absolute();truth=(g37,tab,seg,*base.subarm.hg19_arms())
        refits=[];identities=[]
        for i in range(9):
            directory=a.frozen_panel/f'raw_refit_{i:03d}';d=base.profiles(directory,directory/'haplodose',truth,genes)
            np.testing.assert_array_equal(d['tum'],frozen['samples'])
            for key,old in [('D','dna'),('FL','loss'),('NEU','neutral')]:np.testing.assert_array_equal(d[key],frozen[old])
            for name in ('naive','naive_smooth'):d['W'][name]=d['W'][name].astype(float)*np.log(2.)
            with np.load(a.operator_panel/f'refit_{i:03d}.npz',allow_pickle=True) as old:
                np.testing.assert_array_equal(d['tum'],old['samples']);np.testing.assert_array_equal(d['warm'],old['arms'])
                for name in base.METHODS:np.testing.assert_array_equal(d['W'][name],old[name])
                for name in ('naive_own_partition','hd_own_partition'):d['W'][name]=old[name]
            raw,smoothed,identity=uncorrected(directory,d,base.subarm)
            (naive,noprog),details=op.operator(d['W']['naive'],raw,d['warm'],partition.partition)
            np.testing.assert_array_equal(naive,d['W']['naive_own_partition'])
            d['W'].update(noprog_lr=raw,noprog_lr_smooth=smoothed,noprog_own_partition=noprog)
            identities.append(identity);refits.append(d)
            np.savez_compressed(a.out/f'refit_{i:03d}.npz',samples=d['tum'],arms=d['warm'],**d['W'])
            (a.out/f'operator_{i:03d}.json').write_text(json.dumps(details,indent=2))
            print(json.dumps(dict(reconstructed=i,counts_conserved=True,previous_operators_identical=True)),flush=True)
        base.METHODS=base.METHODS+('naive_own_partition','hd_own_partition','noprog_lr','noprog_lr_smooth','noprog_own_partition')
        base.CONTRASTS=base.CONTRASTS+(('hd_own_partition','naive_own_partition'),('hd_own_partition','naive_smooth'),
            ('hd_lr','noprog_lr'),('hd_lr_smooth','noprog_lr_smooth'),('hd_own_partition','noprog_own_partition'),
            ('noprog_own_partition','naive_own_partition'))
        original_correlations=base.correlations
        def correlations(d,common):
            result=original_correlations(d,common)
            for name in base.METHODS:
                for j in range(len(d['tum'])):
                    x,y=d['W'][name][common[:,j],j],d['D'][common[:,j],j]
                    if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:result[name]['pearson_window'][j]=precision.reference(x,y)['value']
            return result
        base.correlations=correlations
        split=meta.reindex(refits[0]['tum']);ci=np.flatnonzero(split.split.eq('calibration'));ei=np.flatnonzero(split.split.eq('evaluation'))
        if len(ci)!=24 or len(ei)!=24 or split.patient.duplicated().any():raise ValueError('Patient identity changed')
        _,thresholds,records,common=base.score_profiles(refits[0],ci);np.testing.assert_array_equal(common,frozen['common'])
        patient=pd.DataFrame(records);patient['split']=patient['sample'].map(meta.split);patient.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
        comparisons=base.nested_comparisons(refits[0],refits[1:],ci,ei,10000,base.SEED)
        pd.DataFrame(comparisons).to_csv(a.out/'comparisons.tsv',sep='\t',index=False)
        (a.out/'summary.json').write_text(json.dumps(dict(params=params,thresholds=thresholds,methods=base.METHODS,
            contrasts=base.CONTRASTS,comparisons=comparisons,identities=identities,
            interpretation='Component attribution on exposed development data; no default promotion or independently validated coverage'),indent=2))
        for record in rec['input']['files'].values():
            if file_record(record['path'])['sha256']!=record['sha256']:raise ValueError('Frozen input changed')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
