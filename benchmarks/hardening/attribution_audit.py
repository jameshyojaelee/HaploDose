"""Independent numeric discrepancies and frozen source/input preservation."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
a=p.parse_args();a.out.mkdir(exist_ok=False);report={}
for name,root in [('initial',a.task/'baseline'),('previous',a.task.parent/'confirmation_20261107/final_source')]:
 manifest=json.loads((root/'source_manifest.json').read_text())
 for f,h in manifest.items():
  if hashlib.sha256((root/f).read_bytes()).hexdigest()!=h:raise ValueError('Frozen source changed')
 report[name+'_source_files_preserved']=len(manifest)
rows=json.loads((a.task/'results_v1/summary.json').read_text())
if len(rows)!=18 or len({(r['seed'],r['exposure'],r['membership']) for r in rows})!=18:raise ValueError('Incomplete factorial')
report['optimizer_warnings']=sum(not r['optimizer_success'] for r in rows)
report['max_objective_discrepancy']=max(abs(r['independent_objective_discrepancy']) for r in rows)
report['max_scalar_root_parameter_difference']=max(float(np.max(np.abs(np.asarray(r['parameters'])-r['independent_scalar_roots']))) for r in rows if r['independent_scalar_roots'] is not None)
report['independent_optimizer_all_success']=all(r['independent_optimizer_success'] for r in rows)
report['numerical_interpretation']='Independent objective/parameter agreement does not erase convergence warnings or certify a global optimum.'
root=a.task.parent/'hardening_20261106_2119/artifact_diagnostics_v8'
profiles=json.loads((a.task/'profile_v1/summary.json').read_text())
if len(profiles)!=5:raise ValueError('Missing locked profile case')
for r in profiles:
 if hashlib.sha256((root/(r['case']+'_profiles.tsv.gz')).read_bytes()).hexdigest()!=r['input_sha256']:raise ValueError('Frozen profile changed')
report['frozen_profile_inputs_preserved']=len(profiles)
(a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)
