"""Independent direct-comparison audit of cached loss recall."""
import argparse,importlib.util,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

def load(n):
 p=Path(__file__).with_name(n+'.py');s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];factor=root/'dev/gene_factorial_20261108';base=load('profile');cache=load('loss_recall_cache')
 inputs=dict(script=Path(__file__),cache=Path(__file__).with_name('loss_recall_cache.py'),scorer=Path(__file__).with_name('profile.py'),support=root/'dev/hardening_20261106_2119/profile_v6/raw_support_and_truth.npz',thresholds=root/'dev/loss_quantile_identity_20261108/replay_v1/thresholds.tsv')
 for i in range(9):inputs[f'fit{i}']=factor/f'panel_v1/refit_{i:03d}.npz'
 rec=start_run(a.out,inputs,dict(seed=20261108,fixtures=100),dict(scope='Computational recall identity; no new statistical or biological claim'),mode='recall_cache_audit')
 try:
  rng=np.random.default_rng(20261108);synthetic=0;actual=0;maxbytes=0;slow=fast=0.
  def check(Y,L,C,thresholds):
   nonlocal synthetic,actual,maxbytes,slow,fast
   curve=cache.RecallCurve(Y,L,C);maxbytes=max(maxbytes,curve.counts.nbytes+curve.support.nbytes+curve.denominator.nbytes)
   for threshold in thresholds:
    begin=time.monotonic();reference=np.array([row['recall_loss'] for row in base.score_loss(Y,L,~L,C,threshold)]);slow+=time.monotonic()-begin
    begin=time.monotonic();observed=curve.at(threshold);fast+=time.monotonic()-begin
    np.testing.assert_array_equal(observed,reference)
   return len(thresholds)*Y.shape[1]
  for _ in range(100):
   Y=rng.integers(-5,6,(17,5)).astype(float);Y[0]=[np.nan,np.inf,-np.inf,np.nan,0]
   L=rng.random(Y.shape)<.5;C=rng.random(Y.shape)<.8;L[:,0]=False
   support=np.unique(Y[np.isfinite(Y)]);thresholds=[np.nan,np.inf,-np.inf]
   for v in support:thresholds.extend([np.nextafter(v,-np.inf),v,np.nextafter(v,np.inf)])
   synthetic+=check(Y,L,C,thresholds)
  with np.load(inputs['support'],allow_pickle=True) as z:L=z['loss'];C=z['common'];samples=z['samples']
  thresholds=pd.read_csv(inputs['thresholds'],sep='\t',float_precision='round_trip')
  for i in range(9):
   with np.load(inputs[f'fit{i}'],allow_pickle=True) as z:
    np.testing.assert_array_equal(z['samples'],samples)
    for row in thresholds[thresholds.refit.eq(i)].itertuples():
     Y=z[row.method];finite=np.unique(Y[np.isfinite(Y)])
     selected=[finite[0],finite[len(finite)//2],finite[-1],row.before,row.after]
     targets=[t for v in selected for t in [np.nextafter(v,-np.inf),v,np.nextafter(v,np.inf)]]
     actual+=check(Y,L,C,targets)
  for item in rec['input']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
  (a.out/'report.json').write_text(json.dumps(dict(synthetic_patient_threshold_checks=synthetic,real_patient_threshold_checks=actual,exact_equality=True,largest_curve_bytes=maxbytes,baseline_seconds=slow,lookup_seconds=fast,timing_scope='Evaluation only; initialization excluded; no universal speed guarantee'),indent=2))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
