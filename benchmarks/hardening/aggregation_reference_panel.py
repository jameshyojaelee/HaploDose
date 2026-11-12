"""Fixed aggregation/reference factorial on previously exposed raw bulk inputs."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    path=Path(__file__).with_name(name+'.py');s=importlib.util.spec_from_file_location('factor_'+name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--frozen-panel',type=Path,required=True);p.add_argument('--previous',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);base=load('profile');factor=load('aggregation_reference')
    op=load('raw_operator_panel');partition=load('profile_operator_control');precision=load('pearson_reference')
    checks=factor.verify();partition.verify();precision.verify()
    (a.out/'checks.json').write_text(json.dumps(checks,indent=2))
    root=Path(__file__).parent;inputs={name:root/name for name in ('aggregation_reference_panel.py','aggregation_reference.py',
        'profile.py','raw_operator_panel.py','profile_operator_control.py','pearson_reference.py','statistics.py')}
    inputs.update(subarm=root.parent/'tcga_bulk/subarm.py',split=a.frozen_panel/'split_manifest.tsv',
        genes=a.frozen_panel/'gene_universe.txt',support=a.frozen_panel/'raw_support_and_truth.npz',
        previous_manifest=a.previous/'run_manifest.json',previous_summary=a.previous/'summary.json',
        previous_scores=a.previous/'patient_scores.tsv',annotation=base.subarm.GTF37,arms19=base.subarm.CYTO_HG19,
        centromeres38=base.subarm.CEN38,dna_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        dna_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
    for i in range(9):
        directory=a.frozen_panel/f'raw_refit_{i:03d}'
        for f in ('counts.tsv','samples.tsv','positions.tsv','haplodose/samples.tsv','haplodose/bins.tsv',
                  'haplodose/genes.tsv.gz','haplodose/log_ratio.npy','haplodose/segments.tsv','haplodose/run.json'):
            inputs[f'{i}:{f}']=directory/f
        inputs[f'{i}:previous']=a.previous/f'refit_{i:03d}.npz'
    params=dict(bootstrap=20000,seed=base.SEED,family=54)
    rec=start_run(a.out,inputs,params,dict(role='development_only',estimands='count-weighted or median gene RNA effects',
        intervention='pooled-depth versus equal-donor rates; sum counts versus median gene log ratios',
        count_calibration='not claimed for median RNA effects',references='same frozen donors and roles'),mode='aggregation_reference')
    try:
        meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample');genes=inputs['genes'].read_text().splitlines()
        with np.load(inputs['support'],allow_pickle=True) as z:frozen={k:z[k] for k in z.files}
        g37=base.subarm.hg19_genes();tab,seg=base.subarm.absolute();truth=(g37,tab,seg,*base.subarm.hg19_arms());refits=[]
        for i in range(9):
            directory=a.frozen_panel/f'raw_refit_{i:03d}';d=base.profiles(directory,directory/'haplodose',truth,genes)
            np.testing.assert_array_equal(d['tum'],frozen['samples'])
            for key,old in [('D','dna'),('FL','loss'),('NEU','neutral')]:np.testing.assert_array_equal(d[key],frozen[old])
            for key in ('naive','naive_smooth'):d['W'][key]=d['W'][key].astype(float)*np.log(2.)
            with np.load(a.previous/f'refit_{i:03d}.npz',allow_pickle=True) as old:
                np.testing.assert_array_equal(d['tum'],old['samples']);np.testing.assert_array_equal(d['warm'],old['arms'])
                for key in base.METHODS:np.testing.assert_array_equal(d['W'][key],old[key])
                for key in old.files:
                    if key not in ('samples','arms'):d['W'][key]=old[key]
            X=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0);samples=pd.read_csv(directory/'samples.tsv',sep='\t')
            mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t')
            base.subarm.unique_ids(X.index,'counts genes');base.subarm.unique_ids(samples['sample'],'samples')
            base.subarm.unique_ids(mapping.gene,'mapped genes')
            np.testing.assert_array_equal(X.index.to_numpy()[mapping.column.to_numpy(int)],mapping.gene)
            G=X.loc[:,samples['sample']].to_numpy(float)[mapping.column.to_numpy(int)]
            effects,identity=factor.effects(G,mapping.bin.to_numpy(int),samples.role.eq('normal'))
            saved=pd.read_csv(directory/'haplodose/samples.tsv',sep='\t').set_index('sample').reindex(samples['sample'])
            np.testing.assert_array_equal(identity['depth'],saved.binned_counts)
            target=pd.Index(samples['sample']).get_indexer(d['tum'])
            if np.any(target<0):raise ValueError('Unknown target')
            details={};raw_store={}
            for name,values in effects.items():
                gene_values=base.subarm.gene_values(values[target],mapping.set_index('gene').bin,d['genes'])
                gene_values=np.where(d['shared_gene_mapping'][:,None],gene_values,np.nan)
                raw=base.subarm.window_median(gene_values,d['window_id'],len(d['warm']));raw_store[name]=raw
                (naive,fitted),details[name]=op.operator(d['W']['naive'],raw,d['warm'],partition.partition)
                np.testing.assert_array_equal(naive,d['W']['naive_own_partition'])
                if name=='pooled_sum':
                    np.testing.assert_array_equal(raw,d['W']['noprog_lr'])
                    np.testing.assert_array_equal(fitted,d['W']['noprog_own_partition'])
                else:d['W'][name+'_partition']=fitted
            refits.append(d)
            np.savez_compressed(a.out/f'refit_{i:03d}.npz',samples=d['tum'],arms=d['warm'],**d['W'])
            np.savez_compressed(a.out/f'raw_factorial_{i:03d}.npz',**raw_store)
            (a.out/f'operators_{i:03d}.json').write_text(json.dumps(details,indent=2))
            print(json.dumps(dict(refit=i,count_conservation=True,pooled_sum_control_exact=True)),flush=True)
        previous=json.loads((a.previous/'summary.json').read_text())
        base.METHODS=tuple(previous['methods'])+('equal_sum_partition','pooled_median_partition','equal_median_partition')
        base.CONTRASTS=tuple(tuple(x) for x in previous['contrasts'])+(
            ('equal_sum_partition','noprog_own_partition'),('equal_median_partition','pooled_median_partition'),
            ('pooled_median_partition','noprog_own_partition'),('equal_median_partition','equal_sum_partition'),
            ('pooled_median_partition','naive_own_partition'),('equal_median_partition','naive_own_partition'))
        if len(base.METHODS)!=13 or len(base.CONTRASTS)!=18:raise ValueError('Frozen family configuration changed')
        ordinary=base.correlations
        def refined(d,common):
            result=ordinary(d,common)
            for key in base.METHODS:
                for j in range(len(d['tum'])):
                    x,y=d['W'][key][common[:,j],j],d['D'][common[:,j],j]
                    if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:result[key]['pearson_window'][j]=precision.reference(x,y)['value']
            return result
        base.correlations=refined
        split=meta.reindex(refits[0]['tum']);ci=np.flatnonzero(split.split.eq('calibration'));ei=np.flatnonzero(split.split.eq('evaluation'))
        if len(ci)!=24 or len(ei)!=24 or split.patient.duplicated().any():raise ValueError('Patient roles changed')
        _,thresholds,records,common=base.score_profiles(refits[0],ci);np.testing.assert_array_equal(common,frozen['common'])
        patient=pd.DataFrame(records);patient['split']=patient['sample'].map(meta.split);patient.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
        comparisons=base.nested_comparisons(refits[0],refits[1:],ci,ei,20000,base.SEED)
        pd.DataFrame(comparisons).to_csv(a.out/'comparisons.tsv',sep='\t',index=False)
        (a.out/'summary.json').write_text(json.dumps(dict(params=params,methods=base.METHODS,contrasts=base.CONTRASTS,
            thresholds=thresholds,comparisons=comparisons,
            limits='Development mechanism test; finite reference-refit approximation, bootstrap tail Monte Carlo uncertainty and biological applicability remain unvalidated'),indent=2))
        for record in rec['input']['files'].values():
            if file_record(record['path'])['sha256']!=record['sha256']:raise ValueError('Frozen source/input changed')
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
