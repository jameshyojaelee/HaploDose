"""Indexed exact-threshold loss bootstrap draws with validated recall lookup."""
import argparse,importlib.util,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

def load(n):
 p=Path(__file__).with_name(n+'.py');s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);p.add_argument('--start',type=int,required=True);p.add_argument('--stop',type=int,required=True);a=p.parse_args()
 assert 0<=a.start<a.stop<=1000000;a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];factor=root/'dev/gene_factorial_20261108';panel=root/'dev/hardening_20261106_2119/profile_v6';prior=root/'dev/gene_programme_20261108/gene_v1/summary.json'
 base=load('profile');cache=load('loss_recall_cache')
 inputs=dict(script=Path(__file__),scorer=Path(__file__).with_name('profile.py'),cache=Path(__file__).with_name('loss_recall_cache.py'),cache_audit=root/'dev/loss_recall_cache_20261108/audit_v1/run_manifest.json',design=a.task/'DESIGN.md',prior=prior,support=panel/'raw_support_and_truth.npz',split=panel/'split_manifest.tsv')
 for i in range(9):inputs[f'fit{i}']=factor/f'panel_v1/refit_{i:03d}.npz'
 check=json.loads(inputs['cache_audit'].read_text());assert check['status']=='complete';assert file_record(inputs['cache'])['sha256']==check['input']['files']['cache']['sha256']
 rec=start_run(a.out,inputs,dict(start=a.start,stop=a.stop,total=1000000,seed=20261106,family=153),dict(scope='Numerical precision only; no biological validation'),mode='factor_loss_draws')
 try:
  previous=json.loads(prior.read_text());levels=[['mapped','full'],['pooled','logmean'],['half','cpm']];combinations=list(itertools.product(*levels));name=lambda v:'factor_'+'_'.join(v)
  methods=list(previous['methods'])+list(map(name,combinations));contrasts=list(map(tuple,previous['contrasts']))
  for v in combinations:
   for axis in range(3):
    if v[axis]==levels[axis][0]:
     other=list(v);other[axis]=levels[axis][1];contrasts.append((name(other),name(v)))
  contrasts.extend((name(v),'naive_own_partition') for v in combinations);assert len(methods)==25 and len(contrasts)==51
  with np.load(inputs['support'],allow_pickle=True) as z:L=z['loss'];N=z['neutral'];C=z['common'];samples=z['samples']
  meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample').loc[samples];ci=np.flatnonzero(meta.split.eq('calibration'));ei=np.flatnonzero(meta.split.eq('evaluation'));assert len(ci)==len(ei)==24 and meta.patient.nunique()==48
  data=[];curves=[]
  for i in range(9):
   with np.load(inputs[f'fit{i}'],allow_pickle=True) as z:
    np.testing.assert_array_equal(z['samples'],samples);d={m:z[m].copy() for m in methods}
   data.append(d);curves.append({m:cache.RecallCurve(d[m],L,C) for m in methods})
  def effects(scores,chosen):
   out=[]
   for x,y in contrasts:
    v=(scores[x]-scores[y])[chosen];v=v[np.isfinite(v)];out.append(float(v.mean()) if len(v) else np.nan)
   return out
  point={m:curves[0][m].at(base.loss_threshold(data[0][m],N&C,ci)) for m in methods}
  (a.out/'definition.json').write_text(json.dumps(dict(methods=methods,contrasts=contrasts,point_effects=effects(point,ei)),indent=2))
  draws=np.lib.format.open_memmap(a.out/'draws.npy',mode='w+',dtype=float,shape=(a.stop-a.start,51));indices=np.lib.format.open_memmap(a.out/'indices.npy',mode='w+',dtype=np.int16,shape=(a.stop-a.start,49))
  rng=np.random.default_rng(20261106);audited=0
  for iteration in range(a.stop):
   ri=int(rng.integers(8))+1;cal=rng.choice(ci,len(ci),replace=True);chosen=rng.choice(ei,len(ei),replace=True)
   if iteration<a.start:continue
   scored={}
   for m in methods:
    threshold=base.loss_threshold(data[ri][m],N&C,cal);scored[m]=curves[ri][m].at(threshold)
    if iteration in (0,1,9999,19999):
     reference=np.array([r['recall_loss'] for r in base.score_loss(data[ri][m],L,N,C,threshold)]);np.testing.assert_array_equal(scored[m],reference);audited+=1
   draws[iteration-a.start]=effects(scored,chosen);indices[iteration-a.start]=np.r_[ri,cal,chosen]
   if (iteration+1)%1000==0:print(json.dumps(dict(completed=iteration+1,start=a.start,stop=a.stop)),flush=True)
  draws.flush();indices.flush()
  for item in rec['input']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
  (a.out/'report.json').write_text(json.dumps(dict(start=a.start,stop=a.stop,full_scorer_vectors_checked=audited,scientific_target_changed=False)))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
