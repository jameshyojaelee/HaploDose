"""Audit retained arrays/support and independently refine Pearson scores."""
import argparse,hashlib,importlib.util,json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import start_run,finish_run

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 root=a.task.parents[1];panel=a.task/'panel_v1';previous=root/'dev/gene_programme_20261108/gene_v1'
 manifest=json.loads((panel/'run_manifest.json').read_text());assert manifest['status']=='complete'
 inputs=dict(script=Path(__file__),manifest=panel/'run_manifest.json',scores=panel/'patient_scores.tsv',old_scores=previous/'patient_scores.tsv',precision=Path(__file__).with_name('pearson_reference.py'),support=root/'dev/hardening_20261106_2119/profile_v6/raw_support_and_truth.npz')
 for i in range(9):inputs[f'new{i}']=panel/f'refit_{i:03d}.npz';inputs[f'old{i}']=previous/f'refit_{i:03d}.npz'
 run=start_run(a.out,inputs,dict(scope='Array identity, common support and Pearson reference'),dict(no_biological_validation=True),mode='gene_factorial_audit')
 try:
  for record in manifest['input']['files'].values():assert hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()==record['sha256']
  spec=importlib.util.spec_from_file_location('precision',inputs['precision']);precision=importlib.util.module_from_spec(spec);spec.loader.exec_module(precision);precision.verify()
  with np.load(inputs['support'],allow_pickle=True) as z:truth=z['dna'];common=z['common'];samples=z['samples']
  scores=pd.read_csv(inputs['scores'],sep='\t',float_precision='round_trip');assert len(scores)==9*25*48 and not scores.duplicated(['refit','sample','method']).any()
  scores=scores.set_index(['refit','sample','method']);maximum=0.;arrays=0
  for i in range(9):
   with np.load(inputs[f'new{i}'],allow_pickle=True) as new,np.load(inputs[f'old{i}'],allow_pickle=True) as old:
    assert set(old.files)<=set(new.files);np.testing.assert_array_equal(new['samples'],samples)
    for name in old.files:np.testing.assert_array_equal(new[name],old[name]);arrays+=1
    methods=[k for k in new.files if k not in ('samples','arms')];assert len(methods)==25
    actual=np.isfinite(truth)&np.logical_and.reduce([np.isfinite(new[k]) for k in methods]);np.testing.assert_array_equal(actual,common)
    np.testing.assert_array_equal(new['factor_mapped_pooled_half'],new['pooled_gene_partition'])
    for name in methods:
     for j,sample in enumerate(samples):
      mask=common[:,j];x=new[name][mask,j];y=truth[mask,j]
      if len(x)>5 and np.ptp(x)>0 and np.ptp(y)>0:
       value=precision.reference(x,y)['value'];key=(i,sample,name)
       before=scores.loc[key,'pearson_window'];maximum=max(maximum,abs(value-before));scores.loc[key,'pearson_window']=value
  corrected=scores.reset_index();corrected.to_csv(a.out/'patient_scores.tsv',sep='\t',index=False)
  # Historical non-Pearson measurements must remain exactly unchanged in value.
  old=pd.read_csv(inputs['old_scores'],sep='\t',float_precision='round_trip').set_index(['sample','method']).sort_index()
  point=corrected[corrected.refit.eq(0)&corrected.method.isin(old.index.get_level_values('method'))].drop(columns='refit').set_index(['sample','method']).sort_index()
  pd.testing.assert_frame_equal(point.drop(columns='pearson_window'),old.drop(columns='pearson_window'),check_exact=True)
  pd.testing.assert_series_equal(point.pearson_window,old.pearson_window,check_exact=True)
  summary=[]
  for (refit,split),g in corrected.groupby(['refit','split']):
   for endpoint in ['pearson_window','rho_window','recall_loss_failure_aware']:
    wide=g.pivot(index='sample',columns='method',values=endpoint)
    for name in sorted(wide.columns):
     delta=wide[name]-wide.naive_own_partition
     summary.append(dict(refit=int(refit),split=split,method=name,endpoint=endpoint,n=int(delta.notna().sum()),missing=int(delta.isna().sum()),paired_mean=delta.mean()))
  pd.DataFrame(summary).to_csv(a.out/'descriptive_effects.tsv',sep='\t',index=False)
  (a.out/'report.json').write_text(json.dumps(dict(retained_arrays_checked=arrays,reference_fits=9,patient_method_rows=len(scores),max_Pearson_refinement=maximum,old_non_Pearson_point_scores_exact=True,old_Pearson_exact_with_round_trip_parser=True,scientific_criteria_changed=False,scope='Descriptive only; factorial formula audit and153-endpoint uncertainty remain required'),indent=2))
  finish_run(a.out,run,backend='cpu')
 except Exception as e:finish_run(a.out,run,error=e);raise
if __name__=='__main__':main()
