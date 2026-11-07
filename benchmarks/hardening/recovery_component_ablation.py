"""One-component diagnostic ablations on all five frozen development inputs."""
import argparse,importlib.util,json,shutil,time
from dataclasses import asdict
from pathlib import Path
import pandas as pd
from haplodose.clones import CloneParams
from haplodose.pipeline import run
from haplodose.provenance import start_run,finish_run
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True);p.add_argument('--component',choices=['reference_depth','no_programmes'],required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
sp=importlib.util.spec_from_file_location('original_assertions',a.assertion_source);mod=importlib.util.module_from_spec(sp);sp.loader.exec_module(mod)
cp=importlib.util.spec_from_file_location('unchanged_scoring',Path(__file__).with_name('recovery_baseline.py'));score=importlib.util.module_from_spec(cp);cp.loader.exec_module(score)
params=CloneParams(split_mode='thin',**({'depth_source':'reference'} if a.component=='reference_depth' else {'n_programmes':0}))
rec=start_run(a.out,{'script':__file__,'assertions':a.assertion_source,'scorer':Path(__file__).with_name('recovery_baseline.py')},asdict(params),{'scope':'paired development diagnostic; one component, no default promotion'},mode='component_ablation');results=[]
try:
 for seed in (20261031,20261107,20261108,20261109,20261110):
  src=a.panel_root/f'seed_{seed}';out=a.out/f'seed_{seed}';out.mkdir();start=time.monotonic()
  run(src/'data/counts.h5',src/'model_cells.tsv',src/'data/positions.tsv',out/'fit',chunk=97,params=params,**mod.INPUT_CONTRACT)
  before=pd.read_csv(src/'data/cells.tsv',sep='\t');fit=pd.read_csv(out/'fit/cells.tsv.gz',sep='\t');truth=[c for c in before if c.startswith('truth_')];joined=fit.merge(before[['cell_id']+truth],on='cell_id',validate='one_to_one');view=out/'assertion_view';view.mkdir()
  joined.to_csv(view/'cells.tsv.gz',sep='\t',index=False)
  for name in ('qc.json','differences.tsv'):shutil.copyfile(out/'fit'/name,view/name)
  ok=True
  try:mod.test_thin_records_mode_and_recovers_clones({'thin':view})
  except AssertionError:ok=False
  criteria=score.criteria(joined,pd.read_csv(out/'fit/differences.tsv',sep='\t'))
  row=dict(seed=seed,original_assertion_passed=ok,criteria=criteria,seconds=time.monotonic()-start,
           s3_groups=int(joined.loc[joined.specimen.eq('S3') & joined.malignant_call,'clone'].nunique()))
  results.append(row);(out/'summary.json').write_text(json.dumps(row,indent=2));print(json.dumps(row),flush=True)
 (a.out/'summary.json').write_text(json.dumps(dict(component=a.component,params=asdict(params),cases=results,scope='Same inspected simulation inputs, regression/development only. No parameter grid, new independent evaluation or default change.'),indent=2));finish_run(a.out,rec,backend='cpu')
except Exception as err:
 finish_run(a.out,rec,error=err);raise
