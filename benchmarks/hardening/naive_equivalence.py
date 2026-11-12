"""Describe translation/partition differences between algebraically equivalent RNA transforms."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run,file_record
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
root=a.task.parents[1];source=root/'dev/gene_factorial_20261108/panel_v1'
inputs=dict(script=Path(__file__))
for i in range(9):inputs[f'fit{i}']=source/f'refit_{i:03d}.npz'
rec=start_run(a.out,inputs,dict(a='factor_full_logmean_cpm',b='naive_own_partition'),dict(scope='Numerical attribution only; no accuracy selection or biological claim'),mode='naive_equivalence')
try:
 rows=[]
 for i in range(9):
  with np.load(inputs[f'fit{i}'],allow_pickle=True) as z:
   A=z['factor_full_logmean_cpm'];B=z['naive_own_partition'];arms=z['arms'];samples=z['samples']
   assert A.shape==B.shape
   for j,sample in enumerate(samples):
    idx=np.flatnonzero(np.isfinite(A[:,j])&np.isfinite(B[:,j]))
    if not len(idx):rows.append(dict(refit=i,sample=str(sample),windows=0,status='no_common_RNA_support'));continue
    delta=A[idx,j]-B[idx,j];translation=float(np.median(delta));adjacent=(np.diff(idx)==1)&(arms[idx[1:]]==arms[idx[:-1]])
    ac=np.diff(A[idx,j])!=0;bc=np.diff(B[idx,j])!=0
    rows.append(dict(refit=i,sample=str(sample),windows=len(idx),status='descriptive_only',translation=translation,maximum_residual=float(np.max(np.abs(delta-translation))),residual_range=float(np.ptp(delta)),boundary_disagreements=int(np.sum((ac!=bc)&adjacent)),boundaries_a=int(np.sum(ac&adjacent)),boundaries_b=int(np.sum(bc&adjacent))))
 table=pd.DataFrame(rows);table.to_csv(a.out/'profiles.tsv',sep='\t',index=False)
 (a.out/'report.json').write_text(json.dumps(dict(profiles=len(rows),missing_profiles=int(table.windows.eq(0).sum()),profiles_with_boundary_disagreements=int(table.boundary_disagreements.fillna(0).gt(0).sum()),maximum_residual=float(table.maximum_residual.max()),median_residual=float(table.maximum_residual.median()),scope='No tolerance-based equivalence or superiority claim; algebraic transform identity is distinct from finite-precision export identity'),indent=2))
 for v in rec['input']['files'].values():assert file_record(v['path'])['sha256']==v['sha256']
 finish_run(a.out,rec,backend='cpu')
except Exception as e:finish_run(a.out,rec,error=e);raise
