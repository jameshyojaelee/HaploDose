"""Extended-precision independent reconstruction of factorial measurements."""
import argparse,importlib.util,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

def load(n):
 p=Path(__file__).with_name(n+'.py');s=importlib.util.spec_from_file_location('formula_'+n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 base=load('profile');op=load('raw_operator_panel');partition=load('profile_operator_control')
 root=a.task.parents[1];panel=root/'dev/hardening_20261106_2119/profile_v6';output=a.task/'panel_v1'
 manifest=json.loads((output/'run_manifest.json').read_text());assert manifest['status']=='complete'
 inputs={k:Path(v['path']) for k,v in manifest['input']['files'].items()};inputs['independent_script']=Path(__file__)
 for i in range(9):inputs[f'factorial{i}']=output/f'refit_{i:03d}.npz'
 rec=start_run(a.out,inputs,dict(arithmetic='longdouble; weighted-rate and log-difference formulas'),dict(scope='Computational measurement identity only'),mode='factorial_formula_audit')
 try:
  assert np.finfo(np.longdouble).eps<np.finfo(float).eps
  genes=(panel/'gene_universe.txt').read_text().splitlines();tab,seg=base.subarm.absolute();truth=(base.subarm.hg19_genes(),tab,seg,*base.subarm.hg19_arms());checks=[]
  for i in range(9):
   directory=panel/f'raw_refit_{i:03d}';d=base.profiles(directory,directory/'haplodose',truth,genes)
   X=pd.read_csv(directory/'counts.tsv',sep='\t',index_col=0);X.index=base.subarm.canonical_symbols(X.index,'genes')
   mapping=pd.read_csv(directory/'haplodose/genes.tsv.gz',sep='\t');sm=pd.read_csv(directory/'samples.tsv',sep='\t').set_index('sample').loc[X.columns]
   normal=sm.role.eq('normal').to_numpy();G=X.loc[mapping.gene].to_numpy(np.longdouble);bins=mapping.bin.to_numpy(int);tum=X.columns.get_indexer(d['tum']);assert (tum>=0).all()
   with np.load(inputs[f'factorial{i}'],allow_pickle=True) as saved:
    for depth,ref,offset in itertools.product(['mapped','full'],['pooled','logmean'],['half','cpm']):
     L=(G if depth=='mapped' else X.to_numpy(np.longdouble)).sum(axis=0)
     rate=G/L;shift=np.longdouble(.5)/L if offset=='half' else np.full(len(L),np.longdouble(1)/1000000)
     observed=np.log(rate+shift)
     if ref=='logmean':expected=np.mean(observed[:,normal],axis=1)[:,None]
     else:
      weights=L[normal]/L[normal].sum();reference=np.sum(rate[:,normal]*weights,axis=1)
      expected=np.log(reference[:,None]+shift)
     effect=observed-expected
     centers=np.median(np.stack([np.median(effect[bins==b],axis=0) for b in np.unique(bins)]),axis=0)
     effect-=centers
     projected=pd.DataFrame(effect[:,tum],index=mapping.gene,columns=d['tum']).reindex(d['genes']).to_numpy(np.longdouble)
     projected=np.where(d['shared_gene_mapping'][:,None],projected,np.nan)
     window=base.subarm.window_median(projected,d['window_id'],len(d['warm']))
     (_,result),_=op.operator(saved['naive'],window,d['warm'],partition.partition)
     name=f'factor_{depth}_{ref}_{offset}';actual=saved[name]
     finite=np.isfinite(actual)&np.isfinite(result);difference=float(np.max(np.abs(actual[finite]-result[finite])))
     checks.append(dict(refit=i,method=name,max_difference=difference,exact=bool(np.array_equal(actual,result,equal_nan=True))))
     (a.out/'checks.json').write_text(json.dumps(checks,indent=2))
     np.testing.assert_array_equal(result,actual)
   print(json.dumps(dict(refit=i,eight_transforms_exact=True)),flush=True)
  for item in rec['input']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
  (a.out/'report.json').write_text(json.dumps(dict(transforms=len(checks),all_exported_profiles_exact=True,biological_validation=False)))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
