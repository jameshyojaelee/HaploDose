"""Locked dependence-aware interval diagnostic on all frozen simulation outputs."""
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats,special
from haplodose.provenance import start_run,finish_run


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    path=Path(__file__).with_name('gaussian_interval_maximum.py')
    spec=importlib.util.spec_from_file_location('max_reference',path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    previous=json.loads((a.previous/'summary.json').read_text());design=previous['design']
    if design['coverage_target']!=.95 or design['replicates_per_condition']!=1000:raise AssertionError('Original criteria changed')
    inputs={n:a.previous/n for n in ('summary.json','intervals.tsv','simulations.npz','run_manifest.json')}
    inputs.update(script=Path(__file__),method=path)
    params=dict(calibration_n=20000,seed_namespace=[20261108,771],alpha_count=.025,
                alpha_latent_tail=.0125,alpha_calibration_total=.0125,calibration_covariances=3,
                primary_coverage_family=12,original_coverage_target=.95)
    rec=start_run(a.out,inputs,params,dict(scope='Oracle covariance/reference only; no biological validation'),mode='dependent_uncertainty')
    try:
        checks=m.checks();(a.out/'checks.json').write_text(json.dumps(checks,indent=2))
        data=np.load(a.previous/'simulations.npz');J=design['bins'];reset=np.zeros(J,bool);reset[::J//4]=True
        original=pd.read_csv(a.previous/'intervals.tsv',sep='\t',float_precision='round_trip')
        old=original[original.method.eq('simultaneous')].copy()
        if len(old)!=10800:raise AssertionError('Incomplete frozen interval panel')
        calibration={};geometry={};maximums={}
        for ci,condition in enumerate(design['conditions']):
            S=data[f'coverage_{condition}_effect0_covariance']
            report,values,g=m.calibrate(S,reset,n=20000,seed=[20261108,771,ci],tail=.0125,failure=.0125/3)
            # Smaller prefixes are convergence diagnostics only. The final bank
            # size is fixed in advance; no prefix chooses the operating cutoff.
            report['prefix_diagnostics']=[]
            for n in (5000,10000):
                k=m.tolerance_rank(n,.0125,.0125/3)
                report['prefix_diagnostics'].append(dict(n=n,rank=k,critical=float(np.partition(values[:n],k-1)[k-1])))
            calibration[condition]=report;geometry[condition]=g;maximums[condition]=values
        (a.out/'calibration.json').write_text(json.dumps(calibration,indent=2))
        np.savez_compressed(a.out/'calibration_maxima.npz',**maximums)
        records=[];max_old_error=0.
        for (phase,condition,effect),frame in old.groupby(['phase','condition','effect'],sort=False):
            ei=0 if phase=='coverage' else previous['effect_grid'].index(effect)
            key=f'{phase}_{condition}_effect{ei}';Y=data[key+'_counts'];S=data[key+'_covariance']
            if not np.array_equal(S,data[f'coverage_{condition}_effect0_covariance']):raise AssertionError('Covariance changed across phases')
            g=geometry[condition];lookup={(int(b),int(e)-1):i for i,(b,e) in enumerate(zip(g['first'],g['stop']))}
            tail=.025/(2*J);low=np.zeros_like(Y);positive=Y>0
            low[positive]=stats.gamma.ppf(tail,Y[positive]);high=stats.gamma.isf(tail,Y+1)
            with np.errstate(divide='ignore'):lowlog=np.log(low)-np.log(2000.)
            highlog=np.log(high)-np.log(2000.)
            original_critical=stats.norm.isf(.025/(2*len(lookup)))
            for r in frame.to_dict('records'):
                b,e,i=int(r['first_bin']),int(r['last_bin']),int(r['replicate']);n=e-b+1
                sd=g['sum_sd_upper'][lookup[(b,e)]]/n
                count_low=float(lowlog[i,b:e+1].mean());count_high=float(highlog[i,b:e+1].mean())
                error=max(abs(count_low-original_critical*sd-r['lower']),abs(count_high+original_critical*sd-r['upper']))
                allowance=8*J*J*np.finfo(float).eps*max(1.,abs(r['upper']),abs(r['lower']),np.log(2000.))
                if error>allowance:raise ArithmeticError('Original count/latent interval failed reconstruction')
                max_old_error=max(max_old_error,error)
                h=calibration[condition]['critical']*sd;lo,hi=count_low-h,count_high+h
                r.update(method='dependent_maximum',lower=lo,upper=hi,width=hi-lo,
                         covered=bool(lo<=r['true_mu_I']<=hi),excludes_zero=bool(lo>0 or hi<0))
                records.append(r)
        candidate=pd.DataFrame(records)
        joined=pd.concat([old,candidate],ignore_index=True)
        joined.to_csv(a.out/'intervals.tsv',sep='\t',index=False)
        coverage=[]
        for (method,condition,interval),frame in joined[joined.phase.eq('coverage')].groupby(['method','condition','interval'],sort=False):
            n=len(frame);k=int(frame.covered.sum())
            if n!=1000 or len(set(frame.replicate))!=n:raise AssertionError('Coverage endpoint population changed')
            # A complete probability sum can round slightly above one. Its
            # exact log probability cannot exceed zero; retain finite log tails.
            logp=min(0.,float(special.logsumexp(stats.binom.logpmf(np.arange(k+1),n,.95))))
            ci=stats.binomtest(k,n).proportion_ci(confidence_level=1-.05/12,method='exact')
            coverage.append(dict(method=method,condition=condition,interval=interval,n=n,covered=k,
                                 coverage=k/n,simultaneous_cp_low=float(ci.low),simultaneous_cp_high=float(ci.high),log_p=logp))
        if len(coverage)!=12:raise AssertionError('Coverage family incomplete')
        running=-np.inf
        for rank,i in enumerate(sorted(range(12),key=lambda i:coverage[i]['log_p'])):
            running=max(running,min(0.,coverage[i]['log_p']+np.log(12-rank)))
            coverage[i].update(log_p_holm=float(running),falsified_undercoverage=bool(running<np.log(.05)))
        curves=[]
        for keys,frame in joined[joined.phase.eq('power')].groupby(['method','condition','effect','interval'],sort=False):
            if len(frame)!=200:raise AssertionError('Power population changed')
            curves.append(dict(zip(['method','condition','effect','interval'],keys),n=len(frame),
                               covered=int(frame.covered.sum()),zero_exclusions=int(frame.excludes_zero.sum()),
                               median_width=float(frame.width.median()),infinite_width=int(np.isinf(frame.width).sum())))
        report=dict(params=params,coverage=coverage,power=curves,calibration=calibration,
                    original_reconstruction_max_error=max_old_error,
                    error_budget='Count failure .025 + Gaussian maximum tail .0125 + simultaneous calibration-bank failure .0125 <= .05',
                    limit='Known covariance and independent expectation only; conditional bank assurance differs from unconditional calibration-and-data coverage; no DNA confidence claim')
        (a.out/'summary.json').write_text(json.dumps(report,indent=2))
        print(json.dumps(dict(coverage=coverage,critical={c:r['critical'] for c,r in calibration.items()})),flush=True)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
