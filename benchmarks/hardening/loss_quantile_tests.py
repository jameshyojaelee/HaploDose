"""Run focused scorer regressions against the frozen source on a compute node."""
import argparse,json,subprocess,sys
from pathlib import Path
from haplodose.provenance import start_run,finish_run
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
root=Path(__file__).parents[2]
r=start_run(a.out,dict(script=Path(__file__),test=root/'tests/test_hardening_scoring.py',scorer=Path(__file__).with_name('profile.py')),dict(scope='Focused scoring regression'),dict(biological_validation=False),mode='loss_quantile_tests')
try:
 command=[sys.executable,'-I','-m','pytest','--import-mode=importlib','-p','no:cacheprovider','--basetemp',str(a.out/'pytest_tmp'),'-q',str(root/'tests/test_hardening_scoring.py')]
 with (a.out/'pytest.txt').open('w') as h:subprocess.run(command,check=True,stdout=h,stderr=subprocess.STDOUT)
 finish_run(a.out,r,backend='cpu')
except Exception as e:finish_run(a.out,r,error=e);raise
