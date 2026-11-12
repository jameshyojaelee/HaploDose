"""Prespecified 153-endpoint paired analysis; initial Monte Carlo precision only."""
import argparse,importlib.util,itertools,json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record

def load(n):
 p=Path(__file__).with_name(n+'.py');s=importlib.util.spec_from_file_location('uncertainty_'+n,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];panel=root/'dev/hardening_20261106_2119/profile_v6';old=root/'dev/gene_programme_20261108/gene_v1'
 formula=a.task/'formula_v1/report.json';assert json.loads(formula.read_text())['all_exported_profiles_exact']
 inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',formula=formula,audit=a.task/'audit_v2/report.json',old=old/'summary.json',split=panel/'split_manifest.tsv',support=panel/'raw_support_and_truth.npz')
 for f in Path(__file__).parent.glob('*.py'):inputs['source:'+f.name]=f
 for i in range(9):inputs[f'refit{i}']=a.task/f'panel_v1/refit_{i:03d}.npz'
 # Capture the original panel's transitive measurement and truth dependencies.
 manifest=json.loads((a.task/'panel_v1/run_manifest.json').read_text())
 for k,v in manifest['input']['files'].items():inputs['panel:'+k]=Path(v['path'])
 run=start_run(a.out,inputs,dict(bootstrap=20000,seed=20261106,family=153),dict(scope='Exposed development; finite reference ensemble and numerical tails unvalidated'),mode='factorial_uncertainty')
 try:
  base=load('profile');precision=load('pearson_reference');prior=json.loads(inputs['old'].read_text())
  levels=[['mapped','full'],['pooled','logmean'],['half','cpm']];comb=list(itertools.product(*levels));name=lambda v:'factor_'+'_'.join(v)
  edges=[]
  for v in comb:
   for axis in range(3):
    if v[axis]==levels[axis][0]:
     other=list(v);other[axis]=levels[axis][1];edges.append((name(other),name(v)))
  assert len(edges)==12
  base.METHODS=tuple(prior['methods'])+tuple(map(name,comb))
  base.CONTRASTS=tuple(map(tuple,prior['contrasts']))+tuple(edges)+tuple((name(v),'naive_own_partition') for v in comb)
  assert len(base.METHODS)==25 and len(base.CONTRASTS)==51
  ordinary=base.correlations
  def refined(d,common):
   scores=ordinary(d,common)
   for m in base.METHODS:
    for j in range(len(d['tum'])):
     x,y=d['W'][m][common[:,j],j],d['D'][common[:,j],j]
     if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:scores[m]['pearson_window'][j]=precision.reference(x,y)['value']
   return scores
  base.correlations=refined
  genes=(panel/'gene_universe.txt').read_text().splitlines();tab,seg=base.subarm.absolute();truth=(base.subarm.hg19_genes(),tab,seg,*base.subarm.hg19_arms())
  template=base.profiles(panel/'raw_refit_000',panel/'raw_refit_000/haplodose',truth,genes)
  with np.load(inputs['support'],allow_pickle=True) as z:common=z['common'];samples=z['samples']
  refs=[]
  for i in range(9):
   with np.load(inputs[f'refit{i}'],allow_pickle=True) as z:
    np.testing.assert_array_equal(z['samples'],samples);d=dict(template);d['W']={m:z[m].copy() for m in base.METHODS}
   actual=np.isfinite(d['D'])&np.logical_and.reduce([np.isfinite(d['W'][m]) for m in base.METHODS]);np.testing.assert_array_equal(actual,common);refs.append(d)
  meta=pd.read_csv(inputs['split'],sep='\t').set_index('sample').loc[samples];assert meta.patient.nunique()==48
  ci=np.flatnonzero(meta.split.eq('calibration'));ei=np.flatnonzero(meta.split.eq('evaluation'));assert len(ci)==len(ei)==24
  records=base.nested_comparisons(refs[0],refs[1:],ci,ei,20000,20261106);assert len(records)==153
  indexed={(r['endpoint'],r['a'],r['b']):r for r in records}
  for oldrow in prior['comparisons']:
   row=indexed[oldrow['endpoint'],oldrow['a'],oldrow['b']]
   for key in ['mean_diff','n','n_missing','p','testable','valid_bootstrap_replicates']:
    assert row[key]==oldrow[key],(key,row,oldrow)
  pd.DataFrame(records).to_csv(a.out/'comparisons.tsv',sep='\t',index=False)
  (a.out/'summary.json').write_text(json.dumps(dict(methods=base.METHODS,contrasts=base.CONTRASTS,comparisons=records,params=dict(seed=20261106,bootstrap=20000,family=153),limit='Do not promote before Monte Carlo precision audit; biological coverage unvalidated'),indent=2))
  for item in run['input']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
  finish_run(a.out,run,backend='cpu')
 except Exception as e:finish_run(a.out,run,error=e);raise
if __name__=='__main__':main()
