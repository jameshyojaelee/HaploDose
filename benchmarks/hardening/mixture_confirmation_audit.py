"""Stage outcome, unchanged control and unit-correction audit."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
a=p.parse_args();a.out.mkdir(exist_ok=False);report={}
root=a.task.parents[1]
for name,path in [('baseline',a.task/'baseline'),('prior_final',root/'dev/attribution_20261108/final_source')]:
 manifest=json.loads((path/'source_manifest.json').read_text())
 for f,h in manifest.items():
  if hashlib.sha256((path/f).read_bytes()).hexdigest()!=h:raise ValueError('Frozen source changed')
 report[name+'_files_preserved']=len(manifest)
for f in ('tests/test_split_mode.py','python/haplodose/simulate.py','python/haplodose/clones.py','python/haplodose/profile.py'):
 if (root/f).read_bytes()!=(a.task/'baseline'/f).read_bytes():raise ValueError('Original criterion/generator/package changed')
report['original_criteria_generators_and_package_unchanged']=True
cases=[]
for seed in (20261031,20261107,20261108,20261109,20261110):
 previous=root/f'dev/confirmation_20261107/panel_merge_v1/seed_{seed}/fit'
 control=a.task/f'panel_merge_control_v1/seed_{seed}/fit';candidate=a.task/f'panel_mixture_v1/seed_{seed}/fit'
 if not np.array_equal(np.load(previous/'posteriors.npy'),np.load(control/'posteriors.npy'),equal_nan=True):raise ValueError('Control changed')
 for name in ('cells.tsv.gz','differences.tsv','segments.tsv'):
  pd.testing.assert_frame_equal(pd.read_csv(previous/name,sep='\t'),pd.read_csv(control/name,sep='\t'))
 same=np.array_equal(np.load(candidate/'posteriors.npy'),np.load(control/'posteriors.npy'),equal_nan=True)
 summary=json.loads((candidate.parent/'summary.json').read_text())
 cases.append(dict(seed=seed,control_identical=True,candidate_posteriors_identical=same,original_assertion_passed=summary['original_assertion_passed']))
report['cases']=cases
prior=root/'dev/attribution_20261108/profile_v1';current=a.task/'profile_units_v1'
profiles=json.loads((current/'summary.json').read_text());report['profile_units']={}
for r in profiles:
 name=r['case'];old=pd.read_csv(prior/f'{name}_profiles.tsv.gz',sep='\t');new=pd.read_csv(current/f'{name}_profiles.tsv.gz',sep='\t')
 if name!='TCGA-05-4389-01':pd.testing.assert_frame_equal(old,new)
 for col in (c for c in old if not c.endswith('_own_partition')):
  pd.testing.assert_series_equal(old[col],new[col])
 report['profile_units'][name]=dict(scores=r['scores'],input_log_units=r['input_log_units'],output_units=r['partition_output_units'])
report['decision']='Mixture confirmation fails recovery-improvement gate; no default promotion or large null sweep. TCGA shared-unit result supersedes the earlier mixed-unit penalty comparison.'
(a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
