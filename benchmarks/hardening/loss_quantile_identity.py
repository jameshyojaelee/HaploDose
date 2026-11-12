"""Equal-patient quantiles must not depend on redundant identical windows."""
import argparse,importlib.util,json
from pathlib import Path
from fractions import Fraction
import numpy as np
from haplodose.provenance import start_run,finish_run

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
 source=Path(__file__).with_name('profile.py');s=importlib.util.spec_from_file_location('profile',source);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
 rec=start_run(a.out,dict(script=Path(__file__),scorer=source),dict(patients=20,q='1/20',duplicated_windows=list(range(1,101))),dict(scope='Exact numerical counterexample; no biological specimens or performance selection'),mode='loss_quantile_identity')
 try:
  results=[]
  for n in range(1,101):
   values=np.tile(np.arange(20,dtype=float),(n,1));mask=np.ones_like(values,dtype=bool)
   observed=m.loss_threshold(values,mask,np.arange(20))
   # Each patient has one constant score and total mass1/20, independently
   # of how many times that score is represented as identical windows.
   mass=Fraction(n,20*n);assert mass==Fraction(1,20)
   expected=0.
   results.append(dict(windows_per_patient=n,expected=expected,observed=observed,equal=observed==expected))
  (a.out/'report.json').write_text(json.dumps(dict(cases=results,failures=sum(not r['equal'] for r in results),scope='Computational invariance only; historical results unchanged'),indent=2))
  finish_run(a.out,rec,backend='cpu')
 except Exception as e:finish_run(a.out,rec,error=e);raise
if __name__=='__main__':main()
