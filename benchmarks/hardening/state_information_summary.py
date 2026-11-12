"""Conservative upper bound on the specified oracle's level-.05 detection power."""
import argparse
from decimal import Decimal,localcontext
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import binom
from haplodose.provenance import start_run,finish_run,file_record


def power_upper(k,n,precision):
    if k==n:return 1.
    with localcontext() as context:
        context.prec=precision;one=Decimal(1);target=one/40;lo=Decimal(0);hi=one
        for _ in range(100):
            p=(lo+hi)/2;term=(one-p)**n;total=term;odds=p/(one-p)
            for j in range(k):
                term=term*(n-j)/(j+1)*odds;total+=term
            if total>target:lo=p
            else:hi=p
        return float((lo+hi)/2)


def polynomial_cdf(k,n,p):
    """Independent exact-integer coefficients and high-precision scalar powers."""
    with localcontext() as context:
        context.prec=60;p=Decimal(p);q=1-p
        return sum((Decimal(math.comb(n,j))*p**j*q**(n-j) for j in range(k+1)),Decimal(0))


def convolution_ratio(x,delta,sigma):
    """Independent balanced density via posterior-count polynomial convolution."""
    x=np.asarray(x,np.longdouble);delta=np.longdouble(delta);sigma=np.longdouble(sigma)
    n=len(x);half=delta/2;variance=sigma*sigma+half*half
    logw=delta*x/(sigma*sigma);q=1/(1+np.exp(-logw))
    polynomials=[np.array([1-v,v],np.longdouble) for v in q]
    while len(polynomials)>1:
        polynomials=[np.convolve(polynomials[i],polynomials[i+1]) if i+1<len(polynomials)
                     else polynomials[i] for i in range(0,len(polynomials),2)]
    # Normalising constants involving pi cancel between numerator/denominator.
    lp=-np.log(sigma)-(x-half)**2/(2*sigma*sigma)
    lm=-np.log(sigma)-(x+half)**2/(2*sigma*sigma)
    lq=-np.log(variance)/2-x*x/(2*variance)
    iid=np.sum(np.logaddexp(lp,lm)-np.log(np.longdouble(2))-lq)
    prior=np.log(np.longdouble(str(math.comb(n,n//2))))-n*np.log(np.longdouble(2))
    return iid+np.log(polynomials[0][n//2])-prior


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',
        original_generator=a.task.parents[1]/'dev/recovery_next_20261107/rank_source_v1/benchmarks/hardening/global_rank_generator.py')
    data={}
    for condition in ('null','alternative'):
        directory=a.task/(condition+'_v1');manifest=json.loads((directory/'run_manifest.json').read_text())
        if manifest['status']!='complete':raise AssertionError('Incomplete oracle sample')
        settings=manifest['execution']['parameters']
        if (settings['cells'],settings['simulated_specimens'],settings['condition'])!=(600,20000,condition):
            raise AssertionError('Oracle population changed')
        if settings['delta']!=.6*np.log(1.5) or settings['sigma']!=.25:raise AssertionError('Oracle parameters changed')
        for r in manifest['output']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Oracle output changed')
        for r in manifest['input']['files'].values():
            if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Oracle input changed')
        for name in ('run_manifest.json','scores.npz','checks.json','latent.npy','truth_labels.npy'):
            inputs[f'{condition}:{name}']=directory/name
        with np.load(directory/'scores.npz') as z:data[condition]={key:z[key] for key in z.files}
        labels=np.load(directory/'truth_labels.npy',mmap_mode='r')
        if labels.shape!=(20000,600) or not np.all(labels.sum(axis=1)==(0 if condition=='null' else 300)):
            raise AssertionError('Balanced-label identity changed')
    rec=start_run(a.out,inputs,dict(alpha=.05,calibration_failure=.025,power_bound_failure=.025,
        calibration_specimens=20000,alternative_specimens=20000),
        dict(scope='Model-conditional information diagnostic, not biological validation or a changed recovery criterion'),
        mode='state_information_summary')
    try:
        n=20000;rank=int(binom.ppf(.025,n,.95))
        # Independent binomial law for failure count, avoiding a large upper-tail
        # sum. These inequalities prove the largest admissible lower quantile rank.
        below=Decimal(1)-polynomial_cdf(n-rank,n,Decimal(1)/20)
        next_below=Decimal(1)-polynomial_cdf(n-rank-1,n,Decimal(1)/20)
        if not below<=Decimal(1)/40<next_below:raise AssertionError('Calibration quantile rank invalid')
        threshold=float(np.nextafter(np.sort(data['null']['lower'])[rank-1],-np.inf))
        events=np.nextafter(data['alternative']['upper'],np.inf)>threshold;k=int(events.sum())
        upper=power_upper(k,n,80)
        if upper!=power_upper(k,n,60):raise ArithmeticError('Power bound has not converged')
        enclosure=None
        if k<n:
            lo=float(np.nextafter(upper,-np.inf));hi=float(np.nextafter(upper,np.inf))
            cdf_lo=polynomial_cdf(k,n,Decimal.from_float(lo));cdf_hi=polynomial_cdf(k,n,Decimal.from_float(hi))
            if not cdf_hi<Decimal(1)/40<cdf_lo:raise AssertionError('Independent power-bound polynomial bracket failed')
            enclosure=dict(lower=lo,upper=hi,cdf_lower=str(cdf_lo),cdf_upper=str(cdf_hi))
        density_checks=[]
        for condition in ('null','alternative'):
            values=data[condition]['log_likelihood_ratio'];indices={int(np.argmin(values)),int(np.argmax(values)),
                int(np.argmin(abs(values-threshold)))}
            if condition=='null':indices.add(int(np.argsort(data[condition]['lower'])[rank-1]))
            latent=np.load(inputs[f'{condition}:latent.npy'],mmap_mode='r')
            for index in sorted(indices):
                independent=convolution_ratio(latent[index],.6*np.log(1.5),.25)
                difference=abs(independent-np.longdouble(values[index]))
                if difference>data[condition]['rounding_allowance'][index]:
                    raise ArithmeticError('Independent convolution density differs at a decisive sample')
                density_checks.append(dict(condition=condition,index=index,
                    independent_log_ratio=float(independent),recorded_log_ratio=float(values[index]),
                    difference=float(difference),allowance=float(data[condition]['rounding_allowance'][index])))
        ambiguous=int(np.count_nonzero((np.nextafter(data['alternative']['lower'],-np.inf)<=threshold)&events))
        report=dict(calibration_rank=rank,calibration_failure_probability=str(below),
            next_rank_failure_probability=str(next_below),threshold_lower=threshold,
            alternative_conservative_rejections=k,alternative_specimens=n,
            conservative_rejection_fraction=k/n,optimal_level_005_power_upper=1. if enclosure is None else enclosure['upper'],
            power_root_nearest_float=upper,independent_density_checks=density_checks,
            joint_confidence_at_least=.95,numerical_boundary_ambiguous=ambiguous,
            independent_power_enclosure=enclosure,
            maximum_score_rounding_allowance={c:float(d['rounding_allowance'].max()) for c,d in data.items()},
            assumptions='Specified balanced-label P and moment-matched Gaussian-state Q; unobserved exchangeable cell identities; '
                'conditional count law adds no class information given amplitudes; checked numerical enclosures',
            interpretation='Upper bound on the most powerful level-.05 specimen test in this specified oracle experiment, '
                'not a universal biological bound, not cell-membership accuracy, and not a change to earlier criteria',
            biological_specimens=0,earlier_empirical_failures_overridden=False)
        (a.out/'report.json').write_text(json.dumps(report,indent=2));finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
