"""Match the naive baseline's intermediate gene bins before common segmentation."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
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
    return dict(independent_noncommutation_counterexample=True,within_bin_permutation=True,translation_invariance=True)


def binned_naive(directory,d,subarm):
    X=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0)
    X.index=subarm.canonical_symbols(X.index,'count genes')
    meta=pd.read_csv(directory/'samples.tsv',sep='\t').set_index('sample').loc[X.columns]
    mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t')
    mapping.gene=subarm.canonical_symbols(mapping.gene,'mapped genes')
    subarm.unique_ids(X.columns,'count samples');subarm.unique_ids(mapping.gene,'mapped genes')
    normal=meta.role.eq('normal').to_numpy();cpm=X/X.sum(axis=0)*1e6
    log=np.log2(cpm+1);R=log.sub(log.loc[:,normal].mean(axis=1),axis=0)
    # Preserve the frozen naive centering universe. Its membership is fixed by
    # RNA/DNA positional correspondence, never DNA dosage or method scores.
    R=R-R.loc[d['genes']].median(axis=0)
    original=R.loc[d['genes'],d['tum']].to_numpy(np.float32)
    original=np.where(d['shared_gene_mapping'][:,None],original,np.nan)
    np.testing.assert_array_equal(original,d['gene_values']['naive'])
    full=R.loc[mapping.gene,d['tum']].to_numpy(np.float32)
    pooled=median_bins(full,mapping.bin.to_numpy(int))
    projected=subarm.gene_values(pooled.T,mapping.set_index('gene').bin,d['genes'])
    projected=np.where(d['shared_gene_mapping'][:,None],projected,np.nan)
    result=subarm.window_median(projected,d['window_id'],len(d['warm'])).astype(float)*np.log(2.)
    return result,dict(genes_in_bins=len(mapping),bins=len(pooled),scoring_genes=len(d['genes']),
                       normal_samples=meta.index[normal].tolist(),sample_order=d['tum'].tolist(),
                       normalization='original full-library CPM+1, mean normal logCPM, original centering',
                       intervention='median within frozen HD bins, expand back to genes, original window median')


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frozen-panel',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);base=load('profile');op=load('raw_operator_panel')
    partition=load('profile_operator_control');precision=load('pearson_reference')
    (a.out/'checks.json').write_text(json.dumps(checks(),indent=2))
    root=Path(__file__).parent
    inputs={n:root/n for n in ('naive_binning_panel.py','profile.py','raw_operator_panel.py','profile_operator_control.py','pearson_reference.py','statistics.py')}
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
    params=dict(bootstrap=20000,seed=base.SEED,family=63)
    rec=start_run(a.out,inputs,params,dict(scope='Development-only intermediate-binning attribution; no new evaluation'),mode='naive_binning')
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
            raw,identity=binned_naive(directory,d,base.subarm)
            (naive,smoothed),details=op.operator(d['W']['naive'],raw,d['warm'],partition.partition)
            np.testing.assert_array_equal(naive,d['W']['naive_own_partition'])
            d['W']['naive_binned_partition']=smoothed;refits.append(d)
            np.savez_compressed(a.out/f'refit_{i:03d}.npz',samples=d['tum'],arms=d['warm'],**d['W'])
            (a.out/f'operator_{i:03d}.json').write_text(json.dumps(dict(identity=identity,operator=details),indent=2))
            print(json.dumps(dict(refit=i,original_gene_naive_exact=True)),flush=True)
        previous=json.loads(inputs['previous_summary'].read_text());base.METHODS=tuple(previous['methods'])+('naive_binned_partition',)
        base.CONTRASTS=tuple(tuple(x) for x in previous['contrasts'])+(
            ('naive_binned_partition','naive_own_partition'),('pooled_median_partition','naive_binned_partition'),
            ('hd_own_partition','naive_binned_partition'))
        if len(base.METHODS)!=14 or len(base.CONTRASTS)!=21:raise AssertionError('Fixed comparator family changed')
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
