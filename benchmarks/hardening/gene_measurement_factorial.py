"""Frozen eight-way measurement diagnostic; no inferential promotion."""
import argparse,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record
import importlib.util

def load(name):
 p=Path(__file__).with_name(name+'.py');s=importlib.util.spec_from_file_location('factor_'+name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=Path('/gpfs/commons/home/jameslee/Haplodose')
 panel=root/'dev/hardening_20261106_2119/profile_v6';previous=root/'dev/gene_programme_20261108/gene_v1'
 base=load('profile');gene=load('gene_aggregation_panel');op=load('raw_operator_panel');partition=load('profile_operator_control')
 inputs=dict(design=a.task/'DESIGN.md',script=Path(__file__),genes=panel/'gene_universe.txt',split=panel/'split_manifest.tsv',support=panel/'raw_support_and_truth.npz',previous=previous/'summary.json')
 for f in Path(__file__).parent.glob('*.py'):inputs['source:'+f.name]=f
 inputs.update(subarm=Path(__file__).parent.parent/'tcga_bulk/subarm.py',annotation=base.subarm.GTF37,centromeres=base.subarm.CEN38,arms=base.subarm.CYTO_HG19,dna_calls=base.subarm.CN/'TCGA_mastercalls.abs_tables_JSedit.fixed.txt',dna_segments=base.subarm.CN/'TCGA_mastercalls.abs_segtabs.fixed.txt')
 for i in range(9):
  for name in ['counts.tsv','samples.tsv','positions.tsv','haplodose/genes.tsv.gz','haplodose/samples.tsv','haplodose/log_ratio.npy','haplodose/segments.tsv']:inputs[f'{i}:{name}']=panel/f'raw_refit_{i:03d}'/name
  inputs[f'{i}:old']=previous/f'refit_{i:03d}.npz'
 rec=start_run(a.out,inputs,dict(factors=['depth','reference','offset']),dict(scope='Exposed development diagnostic only; no new inferential claim'),mode='gene_factorial')
 try:
  genes=inputs['genes'].read_text().splitlines();meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample')
  with np.load(inputs['support'],allow_pickle=True) as z:frozen={k:z[k] for k in z.files}
  tab,seg=base.subarm.absolute();truth=(base.subarm.hg19_genes(),tab,seg,*base.subarm.hg19_arms())
  old_methods=json.loads(inputs['previous'].read_text())['methods'];records=[]
  for i in range(9):
   directory=panel/f'raw_refit_{i:03d}';d=base.profiles(directory,directory/'haplodose',truth,genes)
   np.testing.assert_array_equal(d['tum'],frozen['samples'])
   for key,saved in [('D','dna'),('FL','loss'),('NEU','neutral')]:np.testing.assert_array_equal(d[key],frozen[saved])
   with np.load(previous/f'refit_{i:03d}.npz',allow_pickle=True) as z:
    for name in old_methods:d['W'][name]=z[name].copy()
   X=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0);X.index=base.subarm.canonical_symbols(X.index,'genes')
   sm=pd.read_csv(directory/'samples.tsv',sep='\t').set_index('sample').loc[X.columns];normal=sm.role.eq('normal').to_numpy()
   mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t');G=X.loc[mapping.gene].to_numpy(float);bins=mapping.bin.to_numpy(int)
   target=X.columns.get_indexer(d['tum']);assert np.all(target>=0)
   names=[]
   for depth,ref,offset in itertools.product(['mapped','full'],['pooled','logmean'],['half','cpm']):
    L=G.sum(axis=0) if depth=='mapped' else X.to_numpy(float).sum(axis=0)
    pseudo=np.full(len(L),.5) if offset=='half' else L/1e6
    if ref=='pooled':
     mu=G[:,normal].sum(axis=1)/L[normal].sum();effect=np.log((G+pseudo)/(mu[:,None]*L+pseudo))
    else:
     log=np.log((G+pseudo)/L);effect=log-log[:,normal].mean(axis=1,keepdims=True)
    effect-=np.median(gene.median_bins(effect,bins),axis=0)
    if (depth,ref,offset)==('mapped','pooled','half'):np.testing.assert_array_equal(effect,gene.gene_effects(G,bins,normal)[0])
    projected=pd.DataFrame(effect[:,target],index=mapping.gene,columns=d['tum']).reindex(d['genes']).to_numpy(float)
    projected=np.where(d['shared_gene_mapping'][:,None],projected,np.nan)
    window=base.subarm.window_median(projected,d['window_id'],len(d['warm']))
    (naive,result),_=op.operator(d['W']['naive'],window,d['warm'],partition.partition)
    np.testing.assert_array_equal(naive,d['W']['naive_own_partition'])
    if (depth,ref,offset)==('mapped','pooled','half'):np.testing.assert_array_equal(result,d['W']['pooled_gene_partition'])
    name=f'factor_{depth}_{ref}_{offset}';d['W'][name]=result;names.append(name)
   base.METHODS=tuple(old_methods)+tuple(names)
   split=meta.reindex(d['tum']);assert split.patient.nunique()==48
   ci=np.flatnonzero(split.split.eq('calibration'));assert len(ci)==24
   _,thresholds,rows,common=base.score_profiles(d,ci);np.testing.assert_array_equal(common,frozen['common'])
   for row in rows:row.update(refit=i,split=meta.loc[row['sample'],'split'])
   records.extend(rows);np.savez_compressed(a.out/f'refit_{i:03d}.npz',samples=d['tum'],arms=d['warm'],**d['W'])
   print(json.dumps(dict(refit=i,control_exact=True,common_support_exact=True)),flush=True)
  pd.DataFrame(records).to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
  for item in rec['input']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
