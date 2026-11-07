"""Conditional two-component NB fit with frozen proposal probabilities.

This is a diagnostic, not independent confirmation: fitted nuisance and proposal
probability uncertainty are omitted. Truth is used only to summarize assignments
AFTER fitting. No pipeline output or scientific gate is changed.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import nbinom


def objective(x,O,E,r,q):
    mean=E[None]*np.exp(x)[:,None,None]
    ll=nbinom.logpmf(O[None],r[None],r[None]/(r[None]+mean)).sum(axis=2)
    with np.errstate(divide='ignore'):
        joint=ll+np.stack([np.log1p(-q),np.log(q)])
    total=logsumexp(joint,axis=0)
    w=np.exp(joint[1]-total)
    score=(r[None]*(O[None]-mean)/(r[None]+mean)).sum(axis=2)
    curv=((O[None]+r[None])*r[None]*mean/(r[None]+mean)**2).sum(axis=2)
    grad=np.array([np.dot(1-w,score[0]),np.dot(w,score[1])])
    t=w*(1-w)
    H=np.array([[np.dot(1-w,curv[0])-np.dot(t,score[0]**2),np.dot(t,score[0]*score[1])],
                [np.dot(t,score[0]*score[1]),np.dot(w,curv[1])-np.dot(t,score[1]**2)]])
    return -total.sum(),-grad,H,w


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--trace-root',type=Path,required=True);p.add_argument('--panel-root',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);rows=[]
    for seed in [20261108,20261109,20261110]:
        root=a.trace_root/f'trace_{seed}'
        for e in json.loads((root/'trace.json').read_text()):
            if e['stage']!='event_fit':continue
            z=np.load(root/e['arrays']);O,E,r,q=(z[k] for k in ['Ytest','E','r','q'])
            fun=lambda x:objective(x,O,E,r,q)[:2]
            result=minimize(fun,z['parameters'][:2],jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':1000})
            ref=minimize(lambda x:objective(x,O,E,r,q)[0],z['parameters'][:2],method='Nelder-Mead',
                         options={'xatol':1e-10,'fatol':1e-10,'maxiter':2000})
            value,gradient,H,w=objective(result.x,O,E,r,q)
            delta=np.array([-1.,1.]);eig=np.linalg.eigvalsh(H)
            se=float(np.sqrt(delta@np.linalg.solve(H,delta))) if np.all(eig>0) else None
            # Check analytic derivatives with step refinement, independent of optimizer stopping.
            derivative_checks=[]
            for j in range(2):
                h=np.cbrt(np.finfo(float).eps)*max(1.,abs(result.x[j]))
                ds=[]
                for step in [h,h/2]:
                    shift=np.zeros(2);shift[j]=step
                    ds.append((fun(result.x+shift)[0]-fun(result.x-shift)[0])/(2*step))
                allowance=4*abs(ds[1]-ds[0])/3+16*np.finfo(float).eps*max(1.,abs(value))/h
                derivative_checks.append(dict(error=float(abs(ds[1]-gradient[j])),allowance=float(allowance)))
                if abs(ds[1]-gradient[j])>allowance:raise RuntimeError('Independent derivative check failed')
            truth=pd.read_csv(a.panel_root/f'seed_{seed}'/'data/cells.tsv',sep='\t')
            labels=truth.iloc[z['rows'][z['cells']]].truth_clone.to_numpy();tumour=np.isin(labels,['clone1','clone2'])
            before=q>.5;after=w>.5;y=labels=='clone2'
            rows.append(dict(seed=seed,trace_index=e['index'],log_rates=result.x.tolist(),
                latent_log_contrast=float(delta@result.x),conditional_SE=se,hessian_eigenvalues=eig.tolist(),
                score=gradient.tolist(),optimizer_success=bool(result.success),optimizer_message=str(result.message),
                independent_optimizer_success=bool(ref.success),independent_objective_difference=float(value-ref.fun),
                derivative_checks=derivative_checks,
                proposal_accuracy=float(np.mean(before[tumour]==y[tumour])),
                combined_assignment_accuracy=float(np.mean(after[tumour]==y[tumour])),
                scope='Postfit diagnostic; truth-free estimation; conditional SE excludes fitted q/nuisance/selection; no calibration or pipeline recovery claim'))
    (a.out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows))


if __name__=='__main__':main()
