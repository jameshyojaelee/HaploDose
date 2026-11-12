"""Numerically checked oracle power envelope for the fixed-population simulator."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.stats import beta, binom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def statistic(values, k, sigma, dtype=np.float64):
    """Fixed-k exchangeable mixture density, not an iid Bernoulli mixture."""
    x=np.asarray(values,dtype=dtype); m,n=x.shape
    variance=dtype(sigma)**2; p=dtype(k)/n
    null_variance=variance+p*(1-p)
    log0=-dtype(.5)*(np.log(2*dtype(np.pi)*variance)+x*x/variance)
    lognull=-dtype(.5)*(np.log(2*dtype(np.pi)*null_variance)+(x-p)**2/null_variance)
    odds=(x-dtype(.5))/variance
    coefficient=np.full((m,k+1),-np.inf,dtype=dtype); coefficient[:,0]=0
    # One logaddexp update per observation. The recurrence is 1-Lipschitz
    # in its inputs in max norm; accumulate a magnitude/operation allowance.
    allowance=np.zeros(m,dtype=dtype);eps=np.finfo(dtype).eps
    for i in range(n):
        high=min(i+1,k)
        if high:
            a=coefficient[:,1:high+1]
            b=coefficient[:,:high]+odds[:,i,None]
            updated=np.logaddexp(a,b)
            allowance+=8*eps*(1+abs(odds[:,i])+np.max(abs(updated),axis=1))
            coefficient[:,1:high+1]=updated
    constant=dtype(math.log(math.comb(n,k)))
    value=np.sum(log0-lognull,axis=1)+coefficient[:,k]-constant
    allowance+=16*eps*(1+np.sum(abs(log0)+abs(lognull),axis=1)+abs(coefficient[:,k])+abs(constant))
    assert np.isfinite(value).all() and np.isfinite(allowance).all()
    return np.asarray(value,float),np.asarray(allowance,float)


def numerical_checks(sigma):
    rng=np.random.default_rng(2026111017); x=rng.normal(size=(7,6)); k=2
    actual,budget=statistic(x,k,sigma)
    reference=[]
    for row in x:
        states=[]
        for chosen in combinations(range(6),k):
            mu=np.zeros(6);mu[list(chosen)]=1
            states.append(float(np.sum(-.5*(np.log(2*np.pi*sigma**2)+(row-mu)**2/sigma**2))))
        null_variance=sigma**2+(k/6)*(1-k/6)
        null=float(np.sum(-.5*(np.log(2*np.pi*null_variance)+(row-k/6)**2/null_variance)))
        reference.append(float(logsumexp(states)-math.log(math.comb(6,k))-null))
    discrepancy=np.abs(actual-reference);assert np.all(discrepancy<=budget)
    permuted,other=statistic(x[:,::-1],k,sigma)
    assert np.all(abs(actual-permuted)<=budget+other)
    mirrored,other=statistic(1-x,6-k,sigma)
    assert np.all(abs(actual-mirrored)<=budget+other)
    for count in (0,6):
        value,error=statistic(x,count,sigma)
        assert np.all(abs(value)<=error)
    assert np.finfo(np.longdouble).eps<np.finfo(float).eps
    large=rng.normal(size=(8,600))
    double,allowance=statistic(large,240,sigma)
    extended,extended_allowance=statistic(large,240,np.longdouble(sigma),np.longdouble)
    gap=abs(double-extended)
    assert np.all(gap<=allowance+extended_allowance)
    return dict(exhaustive_mixture_max_error=float(discrepancy.max()),
                full600_extended_precision_max_error=float(gap.max()),
                full600_max_arithmetic_allowance=float(allowance.max()),
                permutations_and_label_complement_passed=True,
                arithmetic_allowance_scope='Magnitude/operation diagnostic checked independently, not directed-rounding certification')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    generation=args.task.parent/'phase8_membership_trace_20261110/trace_v1/replay/seed_202611098001/generation.npz'
    source=Path(__file__).with_name('gene_veto_fresh.py')
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md',
                  generation=generation,generator=source),
                  dict(calibration_draws=10000,power_draws=20000,chunk=128,threads=4,seed=2026111018),
                  dict(scope='Fixed-population latent-rate oracle, no HaploDose or biological calibration'),
                  mode='fixed_population_oracle')
    try:
        with np.load(generation) as z:
            truth=z['truth'];delta=float(z['delta']);sigma=float(z['sigma'])/delta
        n=len(truth);k=int(truth.sum());assert (n,k)==(600,240)
        checks=numerical_checks(sigma)
        def batch(item):
            arm,start,count=item
            rng=np.random.default_rng(np.random.SeedSequence([2026111018,arm,start]))
            x=rng.normal(size=(count,n))
            if arm==0:x=x*np.sqrt(sigma**2+(k/n)*(1-k/n))+k/n
            else:
                x*=sigma;x[:,:k]+=1
            return statistic(x,k,sigma)
        populations=[]
        for arm,total in ((0,10000),(1,20000)):
            jobs=[(arm,start,min(128,total-start))for start in range(0,total,128)]
            with ThreadPoolExecutor(max_workers=4)as executor:
                parts=list(executor.map(batch,jobs))
            values=np.concatenate([a for a,b in parts]);errors=np.concatenate([b for a,b in parts])
            np.savez_compressed(args.out/('null.npz'if arm==0 else'alternative.npz'),statistic=values,arithmetic_allowance=errors)
            populations.append((values,errors))
        (cal,cal_error),(power,power_error)=populations
        # 97.5% order-statistic bracket for the true95th null percentile.
        # Remaining alpha=.025 supplies a one-sided exact upper power bound.
        alpha_quantile=alpha_power=.025;q=.95
        lower_rank=int(binom.ppf(alpha_quantile/2,len(cal),q))
        upper_rank=int(binom.ppf(1-alpha_quantile/2,len(cal),q))+1
        ordered=np.sort(cal)
        lower=ordered[lower_rank-1]-cal_error.max()
        upper=ordered[upper_rank-1]+cal_error.max()
        possible=int(np.sum(power+power_error>=lower))
        bound=1. if possible==len(power)else float(beta.ppf(1-alpha_power,possible+1,len(power)-possible))
        independent=1. if possible==len(power)else brentq(lambda p:binom.cdf(possible,len(power),p)-alpha_power,
            possible/len(power),1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
        inverse_gap=abs(bound-independent)
        assert inverse_gap<=32*np.finfo(float).eps*(1+abs(bound))
        report=dict(n=n,fixed_clone_cells=k,delta=delta,state_sigma_scaled=sigma,
            null='one Gaussian RNA-state population, matching marginal mean/variance of the alternative',
            alternative='exactly240/600 shifted cells, with independent Gaussian state noise; exchangeable labels',
            observation='ideal per-cell regional log-rate before count sampling, not observed RNA or DNA',
            calibration_draws=len(cal),power_draws=len(power),nominal_FPR=.05,
            order_statistic_ranks=[lower_rank,upper_rank],critical_value_bracket=[float(lower),float(upper)],
            possible_power_events=possible,optimal_oracle_power_upper=bound,
            simultaneous_Monte_Carlo_confidence=.95,alpha_quantile=alpha_quantile,alpha_power=alpha_power,
            independent_binomial_inverse_gap=inverse_gap,numerical_checks=checks,
            calibration_max_arithmetic_allowance=float(cal_error.max()),
            alternative_max_arithmetic_allowance=float(power_error.max()),
            biological_specimens=0,package_thresholds_changed=False,
            original_recovery_failures_excused=False,
            interpretation='Conditional ceiling for this specified matched-moment Gaussian competitor and fixed simulator, subject to the recorded numerical checks; not a universal biological limit or achieved pipeline power')
        write_json(report,args.out/'report.json')
        finish_run(args.out,run,backend='CPU exact fixed-count mixture recurrence and independent Monte Carlo populations')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
