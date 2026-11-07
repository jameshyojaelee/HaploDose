"""Post-fit oracle measurement diagnostic; no truth supplied to HaploDose fitting."""
import argparse,json,sys
from pathlib import Path
import h5py,numpy as np,pandas as pd
from scipy import special,stats
from haplodose.simulate import simulate
p=argparse.ArgumentParser();p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--trace-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False);summaries=[]
for seed in (20261108,20261109,20261110):
 captured=[]
 def observer(frame,event,result):
  if event=='return' and frame.f_code.co_name=='draw' and frame.f_code.co_filename==simulate.__code__.co_filename:
   x=frame.f_locals
   if x.get('cols',{}).get('specimen')=='S1':
    captured.append(dict(label=x['cols']['truth_clone'],rate=x['rate'].copy(),L=x['L'].copy(),gamma=x['gamma'].copy(),phi=x['phi'],depth=x['depth']))
 sys.setprofile(observer)
 try: X,genes,meta,*_=simulate(seed=seed)
 finally:sys.setprofile(None)
 source=a.panel_root/f'seed_{seed}'/'data/counts.h5'
 with h5py.File(source) as h:
  for name,values in [('data',X.data),('indices',X.indices),('indptr',X.indptr)]:
   if not np.array_equal(values,h['X/'+name][:]):raise RuntimeError('Generator does not exactly reproduce frozen counts')
  if not np.array_equal(meta.cell_id.to_numpy(),h['cell_ids'].asstr()[:]):raise RuntimeError('Cell identity mismatch')
 rows=np.flatnonzero(meta.specimen.eq('S1'));Y=X[rows].toarray();L=np.concatenate([x['L'] for x in captured]);phi=captured[0]['phi'];gamma=captured[0]['gamma'];depth=captured[0]['depth']
 if [x['label'] for x in captured]!=['normal','clone1','clone2']:raise RuntimeError('Unexpected population order')
 kept=pd.read_csv(a.trace_root/f'trace_{seed}'/'fit/genes.tsv.gz',sep='\t');idx=kept.loc[kept.bin>=0,'column'].to_numpy(int)
 ll_all=[];ll_kept=[]
 for group in captured:
  rate=group['rate'][None]*(L[:,None]/depth)**gamma[None];rate/=rate.sum(axis=1,keepdims=True);m=L[:,None]*rate
  ll=stats.nbinom.logpmf(Y,1/phi,1/(1+phi*m))
  ll_all.append(ll.sum(axis=1));ll_kept.append(ll[:,idx].sum(axis=1))
 prior=np.array([len(x['L']) for x in captured])/len(L);truth=meta.iloc[rows].truth_clone.to_numpy();names=np.array([x['label'] for x in captured]);table=pd.DataFrame(dict(cell_id=meta.iloc[rows].cell_id.to_numpy(),truth=truth,true_depth=L))
 record=dict(seed=seed,exact_frozen_count_identity=True,oracle_prior=prior.tolist(),gene_count=len(genes),retained_gene_count=len(idx))
 for name,values in [('all_genes',ll_all),('retained_genes',ll_kept)]:
  ll=np.array(values).T;post=special.softmax(ll+np.log(prior),axis=1);pred=names[post.argmax(axis=1)];normal=truth=='normal';tum=~normal
  table[name+'_p_normal']=post[:,0];table[name+'_assignment']=pred
  record[name]=dict(normal_false_calls=int((post[normal,0]<.5).sum()),normal_n=int(normal.sum()),tumour_detected=int((post[tum,0]<.5).sum()),tumour_n=int(tum.sum()),tumour_assignment_accuracy=float((pred[tum]==truth[tum]).mean()))
 table.to_csv(a.out/f'cells_{seed}.tsv',sep='\t',index=False);summaries.append(record)
(a.out/'summary.json').write_text(json.dumps(dict(cases=summaries,scope='Oracle true gene rates, dispersion, per-cell depth and population proportions; post-fit diagnostic only. Uses same frozen measurements; not deployable performance or an independent evaluation.'),indent=2))
