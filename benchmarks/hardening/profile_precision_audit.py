"""Independent reconstruction of fixed-budget bootstrap reports and MC decisions."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta,binom
from haplodose.provenance import start_run,finish_run,file_record


def holm(values):
    return [min(1,max((len(values)-sum(v<u for v in values))*u
        for u in values if u<=p)) for p in values]


def basic_endpoint(ordered,delta,q):
    # Exact rational interpolation of the recorded binary64 values, independent
    # of numpy.quantile. The budget covers rank multiplication, lerp and the
    # final subtraction (at most eight rounded operations).
    exact_h=Fraction.from_float(q)*(len(ordered)-1);j=exact_h.numerator//exact_h.denominator
    weight=exact_h-j;a=float(ordered[j]);b=float(ordered[min(j+1,len(ordered)-1)])
    value=Fraction.from_float(delta)-(1-weight)*Fraction.from_float(a)-weight*Fraction.from_float(b)
    eps=np.finfo(float).eps/2;gamma=8*eps/(1-8*eps)
    rounded_h=float(q*(len(ordered)-1));rank_error=abs(float(Fraction.from_float(rounded_h)-exact_h))
    budget=gamma*(abs(delta)+abs(a)+abs(b)+abs(b-a))+rank_error*abs(b-a)
    return float(value),budget


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--refined',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    mp=a.refined/'run_manifest.json';manifest=json.loads(mp.read_text())
    if manifest['status']!='complete':raise AssertionError('Refinement incomplete')
    for r in manifest['input']['files'].values():
        if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Refinement input changed')
    for name in ('summary.json','bootstrap_draws.npy','paired_patient_effects.npy'):
        r=manifest['output']['files'][name]
        if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Refinement output changed')
    inputs=dict(script=Path(__file__),manifest=mp,summary=a.refined/'summary.json',
        draws=a.refined/'bootstrap_draws.npy',effects=a.refined/'paired_patient_effects.npy',
        previous=Path(manifest['input']['files']['summary']['path']),mc_original=a.task/'audit_v2/report.json')
    rec=start_run(a.out,inputs,dict(scientific_family=93,final_mc_family=279),
        dict(scope='Numerical reproducibility and precision, not biological coverage'),mode='profile_precision_audit')
    try:
        summary=json.loads(inputs['summary'].read_text());old=json.loads(inputs['previous'].read_text())
        original=old['comparisons'];corr=[r for r in original if r['endpoint']!='recall_loss']
        draws=np.load(inputs['draws'],mmap_mode='r');effects=np.load(inputs['effects'])
        if draws.shape!=(1000000,62) or effects.shape!=(9,62,48):raise AssertionError('Changed numerical population')
        source_manifest=json.loads(Path(manifest['input']['files']['manifest']['path']).read_text())
        split=pd.read_csv(source_manifest['input']['files']['split']['path'],sep='\t').set_index('sample')
        with np.load(source_manifest['input']['files']['support']['path'],allow_pickle=True) as z:samples=z['samples']
        split=split.reindex(samples);ci=np.flatnonzero(split.split.eq('calibration'));ei=np.flatnonzero(split.split.eq('evaluation'))
        rng=np.random.default_rng(old['params']['seed']);selected={0,1,19999,20000,199999,200000,999999};checked=0
        eps=np.finfo(float).eps/2
        for i in range(1000000):
            ri=int(rng.integers(8));rng.choice(ci,len(ci),replace=True);idx=rng.choice(ei,len(ei),replace=True)
            if i not in selected:continue
            for k in range(62):
                values=effects[ri+1,k,idx];values=values[np.isfinite(values)]
                if not len(values):
                    if not np.isnan(draws[i,k]):raise AssertionError('Missing-draw identity differs')
                    continue
                expected=math.fsum(values)/len(values)
                gamma=(len(values)+1)*eps/(1-(len(values)+1)*eps)
                allowance=gamma*math.fsum(abs(float(v)) for v in values)/len(values)
                if abs(float(draws[i,k])-expected)>allowance:raise AssertionError('Independent scalar draw differs')
                checked+=1
        final_mc={};prefix_reports=[];tail=.05/(2*279)
        for prefix in summary['prefixes']:
            n=prefix['draws'];records=prefix['comparisons']
            if n not in (20000,200000,1000000) or len(records)!=93:raise AssertionError('Changed checkpoint/family')
            for r,previous in zip(records,original):
                for key in ('endpoint','a','b','n','n_missing','mean_diff'):
                    if r[key]!=previous[key]:raise AssertionError('Scientific estimand/population changed')
                if r['endpoint']=='recall_loss':
                    for key in ('p','ci_low','ci_high','family_ci_low','family_ci_high','valid_bootstrap_replicates'):
                        if r[key]!=previous[key]:raise AssertionError('Loss result changed')
            refined=[r for r in records if r['endpoint']!='recall_loss'];probabilities=[];endpoint_checks=0
            for j,r in enumerate(refined):
                values=np.asarray(draws[:n,j]);values=values[np.isfinite(values)];delta=r['mean_diff'];errors=values-delta
                k=sum(int(np.sum(np.abs(errors[start:start+10000])>=abs(delta))) for start in range(0,len(errors),10000))
                expected=Fraction(k+1,len(errors)+1)
                if not r['testable']:
                    expected=Fraction(1)
                if float(expected)!=r['p']:raise AssertionError('Independent probability count differs')
                probabilities.append(expected)
                if not r['testable']:continue
                ordered=np.sort(errors)
                for field,q in (('ci_low',.975),('ci_high',.025),
                    ('family_ci_low',1-.05/(2*93)),('family_ci_high',.05/(2*93))):
                    value,budget=basic_endpoint(ordered,delta,q)
                    if abs(value-r[field])>budget:raise AssertionError('Independent quantile interpolation differs')
                    endpoint_checks+=1
                if n==1000000:
                    regions=[]
                    for q in (.05/(2*93),1-.05/(2*93)):
                        lo=int(binom.ppf(tail,len(errors),q));hi=int(binom.isf(tail,len(errors),q))+1
                        if lo<1 or hi>len(errors):raise AssertionError('Final numerical quantile bound unbounded')
                        if binom.cdf(lo-1,len(errors),q)>tail or binom.sf(hi-1,len(errors),q)>tail:
                            raise AssertionError('Invalid numerical quantile rank')
                        regions.append((float(ordered[lo-1]),float(ordered[hi-1])))
                    final_mc[(r['endpoint'],r['a'],r['b'])]=dict(exceedances=k,n=len(errors),
                        p_lower=0. if k==0 else float(beta.ppf(tail,k,len(errors)-k+1)),
                        p_upper=1. if k==len(errors) else float(beta.isf(tail,k+1,len(errors)-k)),
                        family_ci_low_mc=[delta-regions[1][1],delta-regions[1][0]],
                        family_ci_high_mc=[delta-regions[0][1],delta-regions[0][0]])
            pmap={(r['endpoint'],r['a'],r['b']):v for r,v in zip(refined,probabilities)}
            exact=holm([pmap.get((r['endpoint'],r['a'],r['b']),Fraction.from_float(r['p'])) for r in records])
            for r,reference in zip(records,exact):
                # One rounding of p, multiplication by at most93, and output conversion.
                allowance=3*eps*93*max(float(reference),r['p'])
                if abs(r['p_holm']-float(reference))>allowance:raise AssertionError('Independent Holm differs')
                if (float(reference)<.05)!= (r['p_holm']<.05):raise AssertionError('Holm decision differs')
            prefix_reports.append(dict(draws=n,independent_probability_counts=62,independent_interval_endpoints=endpoint_checks))
        original_mc=json.loads(inputs['mc_original'].read_text())['results'];bounds=[]
        final=summary['prefixes'][-1]['comparisons']
        for r,previous_mc in zip(final,original_mc):
            key=(r['endpoint'],r['a'],r['b'])
            if key in final_mc:b=final_mc[key]
            else:b=dict(p_lower=previous_mc.get('bootstrap_exceedance_probability_mc_lower',1.),
                        p_upper=previous_mc.get('bootstrap_exceedance_probability_mc_upper',1.))
            bounds.append(dict(endpoint=r['endpoint'],a=r['a'],b=r['b'],**b))
        lows=holm([r['p_lower'] for r in bounds]);highs=holm([r['p_upper'] for r in bounds])
        for r,lo,hi,scientific in zip(bounds,lows,highs,final):
            r.update(holm_mc_lower=lo,holm_mc_upper=hi,
                positive_claim_numerically_stable=bool(scientific['superior'] and hi<.05 and
                    r.get('family_ci_low_mc',[-np.inf])[0]>0))
        (a.out/'report.json').write_text(json.dumps(dict(prefixes=prefix_reports,
            independent_scalar_draw_endpoints=checked,final_mc=bounds,scientific_family=93,
            numerical_family=279,biological_units_added=0,
            scope='Numerical precision only; biological bootstrap coverage and reference adequacy not established'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
