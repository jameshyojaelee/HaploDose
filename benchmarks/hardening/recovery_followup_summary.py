"""Before/after same-input traces and independent saved-proposal optimization."""
import argparse,json,subprocess,sys
from pathlib import Path
import numpy as np
from scipy import optimize,special,stats
from haplodose.clones import _fit_event_mixture
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
for name,root in [('before',a.root),('after',a.root/'after_retry')]:
 subprocess.run([sys.executable,'-I',str(Path(__file__).with_name('recovery_cause_report.py')),'--root',str(root),'--out',str(a.out/name)],check=True)
r=a.root/'trace_20261108';events=json.loads((r/'trace.json').read_text());e=next(x for x in events if x['stage']=='event_mixtures');beta=next(x['beta'] for x in events if x['stage']=='difference_segments')
with np.load(r/e['arrays']) as z:
 keep=(z['cols']>=209)&(z['cols']<=217)
 O=z['Yp'][:,keep];E=z['Ep'][:,keep];R=np.broadcast_to(z['Rnb'],z['Yp'].shape)[:,keep]
fit=_fit_event_mixture(O,E,R,beta*np.log(1.5))
def component(eta):return stats.nbinom.logpmf(O,R,R/(R+E*np.exp(eta))).sum(axis=1)
null=optimize.minimize_scalar(lambda x:-component(x).sum(),bounds=(-10,10),method='bounded',options={'xatol':np.sqrt(np.finfo(float).eps)})
null_ll=-null.fun
def objective(x):return -special.logsumexp(np.stack([component(x[0])-np.logaddexp(0,x[2]),component(x[1])-np.logaddexp(0,-x[2])]),axis=0).sum()
ref=optimize.minimize(objective,fit['parameters'],method='Nelder-Mead',bounds=[(-10,10),(-10,10),(-np.log(999),np.log(999))],options={'maxiter':3000,'xatol':np.sqrt(np.finfo(float).eps),'fatol':np.finfo(float).eps*abs(null_ll)})
allowance=fit['numerical_gain_budget']+16*O.size*np.finfo(float).eps*max(1,abs(null_ll),abs(ref.fun))
if not ref.success or abs(fit['gain']-(-ref.fun-null_ll))>allowance:raise RuntimeError('Independent proposal reference disagrees')
(a.out/'optimizer.json').write_text(json.dumps(dict(solver=fit['solver'],gain=fit['gain'],reference_gain=-ref.fun-null_ll,difference=fit['gain']+ref.fun+null_ll,roundoff_allowance=allowance,local_gain_error=fit['local_quadratic_gain_error'],numerical_budget=fit['numerical_gain_budget'],scope='same saved counts/exposures; independent full NB pmf/Nelder-Mead local optimum, no global-optimum claim'),indent=2))
