"""Independent gate before expanding exact-threshold loss bootstrap draws."""
import argparse,importlib.util,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

def load(n):
 p=Path(__file__).with_name(n+'.py');s=importlib.util.spec_from_file_location(n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];prefix=a.task/'prefix_v1';full=root/'dev/gene_factorial_20261108/exact_uncertainty_v1'
 inputs=dict(script=Path(__file__),reference=Path(__file__).with_name('profile.py'),interval_reference=Path(__file__).with_name('profile_precision_audit.py'),prefix_manifest=prefix/'run_manifest.json',full_manifest=full/'run_manifest.json',summary=full/'summary.json',definition=prefix/'definition.json',draws=prefix/'draws.npy',indices=prefix/'indices.npy')
 pm=json.loads(inputs['prefix_manifest'].read_text());fm=json.loads(inputs['full_manifest'].read_text());assert pm['status']==fm['status']=='complete'
 for manifest,names in [(pm,['definition.json','draws.npy','indices.npy']), (fm,['summary.json'])]:
  for name in names:
   item=manifest['output']['files'][name];assert file_record(item['path'])['sha256']==item['sha256']
 for manifest in [pm,fm]:
  for v in manifest['input']['files'].values():assert file_record(v['path'])['sha256']==v['sha256']
 for name in ['support','split']:
  inputs[name]=Path(pm['input']['files'][name]['path'])
 for i in range(9):inputs[f'fit{i}']=Path(pm['input']['files'][f'fit{i}']['path'])
 rec=start_run(a.out,inputs,dict(draws=20000,scientific_family=153),dict(scope='Computational identity gate; no biological claim'),mode='loss_prefix_audit')
 try:
  base=load('profile');reference=load('profile_precision_audit');definition=json.loads(inputs['definition'].read_text());summary=json.loads(inputs['summary'].read_text());rows=[r for r in summary['comparisons'] if r['endpoint']=='recall_loss']
  assert len(rows)==51 and definition['contrasts']==[[r['a'],r['b']] for r in rows]
  draws=np.load(inputs['draws'],mmap_mode='r');indices=np.load(inputs['indices'],mmap_mode='r');assert draws.shape==(20000,51) and indices.shape==(20000,49)
  with np.load(inputs['support'],allow_pickle=True) as z:L=z['loss'];N=z['neutral'];C=z['common'];samples=z['samples']
  meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample').loc[samples];ci=np.flatnonzero(meta.split.eq('calibration'));ei=np.flatnonzero(meta.split.eq('evaluation'));assert meta.patient.nunique()==48 and len(ci)==len(ei)==24
  rng=np.random.default_rng(20261106);scalar_checks=0
  for i in range(20000):
   ri=int(rng.integers(8))+1;cal=rng.choice(ci,len(ci),replace=True);chosen=rng.choice(ei,len(ei),replace=True)
   np.testing.assert_array_equal(indices[i],np.r_[ri,cal,chosen])
   if i not in (0,1,9999,19999):continue
   with np.load(inputs[f'fit{ri}'],allow_pickle=True) as z:
    scored={}
    for method in definition['methods']:
     threshold=base.loss_threshold(z[method],N&C,cal)
     scored[method]=np.array([r['recall_loss'] for r in base.score_loss(z[method],L,N,C,threshold)])
   for k,(x,y) in enumerate(definition['contrasts']):
    values=(scored[x]-scored[y])[chosen];values=values[np.isfinite(values)]
    if not len(values):assert np.isnan(draws[i,k]);continue
    scalar=math.fsum(values)/len(values);eps=np.finfo(float).eps/2;gamma=(len(values)+1)*eps/(1-(len(values)+1)*eps);budget=gamma*math.fsum(abs(float(v)) for v in values)/len(values)
    assert abs(scalar-draws[i,k])<=budget;scalar_checks+=1
  endpoint_checks=0
  expected_n=int(np.sum((L&C).sum(axis=0)[ei]>0))
  for j,row in enumerate(rows):
   delta=definition['point_effects'][j];assert delta==row['mean_diff'] and row['n']==expected_n and row['n_missing']==24-expected_n
   boot=np.asarray(draws[:,j]);boot=boot[np.isfinite(boot)];assert len(boot)==row['valid_bootstrap_replicates']
   testable=bool(np.isfinite(delta) and len(boot)>=18000 and np.ptp(boot)>0);assert testable==row['testable']
   if not testable:assert row['p']==1.;continue
   error=boot-delta;k=int(np.count_nonzero(np.abs(error)>=abs(delta)));assert float(Fraction(k+1,len(boot)+1))==row['p'];ordered=np.sort(error)
   for field,q in [('ci_low',.975),('ci_high',.025),('family_ci_low',1-.05/306),('family_ci_high',.05/306)]:
    direct=float(delta-np.quantile(error,q));assert direct==row[field]
    scalar,budget=reference.basic_endpoint(ordered,delta,q);assert abs(scalar-direct)<=budget;endpoint_checks+=1
  (a.out/'report.json').write_text(json.dumps(dict(resampling_rows_exact=20000,loss_contrasts_exact=51,independent_scalar_draw_effects=scalar_checks,interval_endpoints_checked=endpoint_checks,continuation_gate_passed=True,biological_validation=False),indent=2))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
