"""Balanced-clone versus continuous-state oracle; no HaploDose fitting or tuning."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import norm
from haplodose.provenance import start_run,finish_run


def log_ratio(x,delta,sigma,dtype=np.float64):
    x=np.asarray(x,dtype=dtype);n=x.shape[1];k=n//2
    if n%2 or sigma<=0:raise ValueError('Positive sigma and even cell count required')
    delta,sigma=dtype(delta),dtype(sigma);half=delta/2;variance=sigma*sigma+half*half
    logw=delta*x/(sigma*sigma)
    dp=np.full((len(x),k+1),-np.inf,dtype=dtype);dp[:,0]=0
    for i in range(n):
        stop=min(i+1,k)
        dp[:,1:stop+1]=np.logaddexp(dp[:,1:stop+1],dp[:,:stop]+logw[:,i,None])
    constant=dtype(n)/2*np.log1p(half*half/(sigma*sigma))
    first=np.sum((x+half)**2,axis=1)/(2*sigma*sigma)
    second=np.sum(x*x,axis=1)/(2*variance)
    combinatorial=np.log(dtype(str(math.comb(n,k))))
    values=constant-first+second+dp[:,k]-combinatorial
    # logaddexp is 1-Lipschitz in its input sup norm. Count normal-density,
    # ratio and recurrence operations along a path, rather than using a fixed
    # epsilon independent of the number/scale of cells.
    unit=np.finfo(dtype).eps/2;operations=20*n+32
    gamma=operations*unit/(1-operations*unit)
    budget=gamma*(1+abs(constant)+abs(first)+abs(second)+abs(dp[:,k])+abs(combinatorial)+np.sum(abs(logw),axis=1))
    return values,budget


def checks(delta,sigma):
    rng=np.random.default_rng([20261108,700]);reports=[]
    for n in (4,6,8):
        x=rng.normal(0,sigma,size=(3,n));computed,budget=log_ratio(x,delta,sigma)
        for i,row in enumerate(x):
            densities=[]
            for chosen in itertools.combinations(range(n),n//2):
                mean=np.full(n,-delta/2);mean[list(chosen)]=delta/2
                densities.append(float(norm.logpdf(row,mean,sigma).sum()))
            independent=float(logsumexp(densities)-math.log(len(densities))-
                norm.logpdf(row,0,math.sqrt(sigma*sigma+delta*delta/4)).sum())
            if abs(computed[i]-independent)>budget[i]:raise AssertionError('Exhaustive density disagreement')
            reports.append(dict(cells=n,exhaustive_difference=float(computed[i]-independent),budget=float(budget[i])))
    x=rng.normal(0,sigma,size=(4,600))
    zero,budget=log_ratio(x,0,sigma)
    if np.any(abs(zero)>budget):raise AssertionError('Zero-effect likelihood ratio is not zero')
    value,budget=log_ratio(x,delta,sigma)
    other,_=log_ratio(x[:,rng.permutation(600)],delta,sigma)
    if np.any(abs(value-other)>2*budget):raise AssertionError('Balanced density not permutation invariant')
    if np.finfo(np.longdouble).eps>=np.finfo(float).eps:raise RuntimeError('Extended precision unavailable')
    extended,_=log_ratio(x,delta,sigma,np.longdouble)
    if np.any(abs(np.asarray(value,np.longdouble)-extended)>budget):raise ArithmeticError('Extended precision outside operation budget')
    return dict(exhaustive=reports,zero_effect_max_error=float(abs(zero).max()),
        full_size_permutation_max_difference=float(abs(value-other).max()),
        full_size_extended_max_difference=float(abs(np.asarray(value,np.longdouble)-extended).max()),
        full_size_max_budget=float(budget.max()),
        scope='Operation-count error estimates checked against independent densities and extended precision; not formal interval arithmetic')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--condition',choices=('null','alternative'),required=True);a=p.parse_args()
    n,total=600,20000;delta=.6*np.log(1.5);sigma=.25
    a.out.mkdir(exist_ok=False)
    rec=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md'),
        dict(cells=n,simulated_specimens=total,delta=delta,sigma=sigma,
            condition=a.condition,seed=[20261108,701 if a.condition=='null' else 702]),
        dict(scope='Noiseless latent amplitude diagnostic only; no biological data or method fitting'),mode='state_information')
    try:
        (a.out/'checks.json').write_text(json.dumps(checks(delta,sigma),indent=2))
        rng=np.random.default_rng([20261108,701 if a.condition=='null' else 702])
        amplitude=np.lib.format.open_memmap(a.out/'latent.npy',mode='w+',dtype=np.float64,shape=(total,n))
        labels=np.lib.format.open_memmap(a.out/'truth_labels.npy',mode='w+',dtype=np.int8,shape=(total,n))
        scores=np.empty(total);budgets=np.empty(total)
        for start in range(0,total,128):
            stop=min(start+128,total);x=rng.normal(size=(stop-start,n));g=np.zeros(x.shape,np.int8)
            if a.condition=='null':x*=np.sqrt(sigma*sigma+delta*delta/4)
            else:
                x*=sigma
                for i in range(len(x)):g[i]=(rng.permutation(n)<n//2).astype(np.int8)
                x+=(2*g-1)*delta/2
                if not np.all(g.sum(axis=1)==n//2):raise AssertionError('Changed balanced alternative')
            amplitude[start:stop]=x;labels[start:stop]=g
            scores[start:stop],budgets[start:stop]=log_ratio(x,delta,sigma)
            if not np.isfinite(scores[start:stop]).all():raise ArithmeticError('Nonfinite oracle score')
            if stop%2048==0 or stop==total:print(json.dumps(dict(condition=a.condition,completed=stop,total=total)),flush=True)
        amplitude.flush();labels.flush()
        np.savez(a.out/'scores.npz',log_likelihood_ratio=scores,rounding_allowance=budgets,
            lower=scores-budgets,upper=scores+budgets)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
