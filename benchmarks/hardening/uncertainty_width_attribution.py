"""Frozen interval-width decomposition and infinite-depth Gaussian reference.

No new confidence procedure, fit, seed, operating threshold or primary family.
The infinite-depth reference is conditional on known covariance and expectation.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import integrate, stats
from haplodose.provenance import start_run,finish_run


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(exist_ok=False)
    s=json.loads((a.previous/'summary.json').read_text());design=s['design']
    if design['bins']!=120 or design['coverage_target']!=.95:
        raise AssertionError('Frozen design changed')
    inputs={name:a.previous/name for name in ('summary.json','intervals.tsv','simulations.npz','run_manifest.json')}
    inputs['script']=Path(__file__)
    rec=start_run(a.out,inputs,dict(scope='Width attribution only',original_design=design),
                  dict(oracle='Known generating covariance and independent expectation; no biological guarantee'),
                  mode='uncertainty_width_attribution')
    try:
        frame=pd.read_csv(a.previous/'intervals.tsv',sep='\t',float_precision='round_trip')
        selected=frame[frame.method.eq('simultaneous')]
        if len(selected)!=10800: raise AssertionError('Incomplete original coverage/power panel')
        arrays=np.load(a.previous/'simulations.npz')
        J=design['bins'];length=J//4;family=4*length*(length+1)//2
        alpha=1-design['coverage_target']
        count_tail=(alpha/2)/(2*J);critical=float(stats.norm.isf((alpha/2)/(2*family)))
        eps=np.finfo(float).eps
        output=[];limits=[];maximum_error=0.
        effects=s['effect_grid']
        # Independent dense scalar sums are the reference for covariance and
        # count-log bounds; compare every width to the preserved output.
        for (phase,condition,effect),table in selected.groupby(['phase','condition','effect'],sort=False):
            ei=0 if phase=='coverage' else effects.index(effect)
            key=f'{phase}_{condition}_effect{ei}'
            Y=arrays[key+'_counts'];S=arrays[key+'_covariance'];mu=arrays[key+'_systematic_mu']
            required=design['replicates_per_condition'] if phase=='coverage' else design['power_replicates_per_condition_effect']
            if len(Y)!=required or len(table)!=2*required: raise AssertionError('Missing interval population')
            low=np.zeros_like(Y);positive=Y>0
            low[positive]=stats.gamma.ppf(count_tail,Y[positive])
            high=stats.gamma.isf(count_tail,Y+1)
            with np.errstate(divide='ignore'): logwidth=np.log(high)-np.log(low)
            for r in table.to_dict('records'):
                i,b,e=int(r['replicate']),int(r['first_bin']),int(r['last_bin']);n=e-b+1
                if b//length!=e//length: raise AssertionError('Interval crossed genomic reset')
                terms=S[b:e+1,b:e+1].ravel().tolist()
                variance=math.fsum(terms)/(n*n)
                gamma=len(terms)*eps/(1-len(terms)*eps)
                padded=(math.fsum(terms)+gamma*math.fsum(abs(t) for t in terms))/(n*n)
                count_width=math.fsum(logwidth[i,b:e+1].tolist())/n
                latent_width=2*critical*math.sqrt(padded)
                width=count_width+latent_width
                # Quantiles are the same library, but combination/summation is
                # independent. Account for log(E) subtraction in saved bounds.
                allowance=8*(J*J)*eps*max(1.,abs(width),math.log(2000.))
                discrepancy=abs(width-r['width'])
                if discrepancy>allowance: raise ArithmeticError('Saved interval width does not reconstruct')
                maximum_error=max(maximum_error,discrepancy)
                output.append(dict(phase=phase,condition=condition,effect=effect,interval=r['interval'],
                                   replicate=i,first_bin=b,last_bin=e,count_width=count_width,
                                   latent_width=latent_width,latent_mean_sd=math.sqrt(variance),
                                   width=width,latent_fraction=latent_width/width))
            if phase=='power':
                b=(length-design['fixed_bins'])//2;e=b+design['fixed_bins']
                variance=math.fsum(S[b:e,b:e].ravel().tolist())/(e-b)**2
                sd=math.sqrt(variance)
                # Infinite depth removes Poisson uncertainty, not latent effects.
                # Fixed-interval two-sided Gaussian exclusion probabilities are
                # exact under the original model. No claim for selected intervals.
                for label,z in [('current_all_interval_latent_bound',critical),
                                ('known_covariance_fixed_interval_reference',float(stats.norm.isf(alpha/2)))]:
                    half=z*sd
                    power=float(stats.norm.cdf((-half-effect)/sd)+stats.norm.sf((half-effect)/sd))
                    lo,elo=integrate.quad(lambda x:stats.norm.pdf(x,loc=effect,scale=sd),-np.inf,-half)
                    hi,ehi=integrate.quad(lambda x:stats.norm.pdf(x,loc=effect,scale=sd),half,np.inf)
                    if abs(power-lo-hi)>elo+ehi+16*eps: raise ArithmeticError('Independent Gaussian power discrepancy')
                    limits.append(dict(condition=condition,effect=effect,reference=label,latent_sd=sd,
                                       infinite_depth_width=2*half,infinite_depth_power=power,
                                       quadrature_error_estimate=elo+ehi))
        result=pd.DataFrame(output);result.to_csv(a.out/'decomposition.tsv',sep='\t',index=False)
        summaries=[]
        for key,g in result.groupby(['phase','condition','effect','interval'],sort=False):
            summaries.append(dict(zip(['phase','condition','effect','interval'],key),n=len(g),
                median_count_width=float(g.count_width.median()),median_latent_width=float(g.latent_width.median()),
                median_latent_fraction=float(g.latent_fraction.median()),median_total_width=float(g.width.median())))
        report=dict(reconstructed_intervals=len(output),maximum_width_discrepancy=maximum_error,
                    original_interval_family=family,latent_critical=critical,decomposition=summaries,
                    infinite_depth_fixed_interval_references=limits,
                    interpretation='Mechanism attribution, not a new calibrated interval. Fixed-interval oracle cannot replace selection-aware bounds; known covariance/expectation remain unsupported in routine biological use.')
        (a.out/'summary.json').write_text(json.dumps(report,indent=2))
        print(json.dumps(dict(reconstructed=len(output),maximum_width_discrepancy=maximum_error)),flush=True)
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error)
        raise


if __name__=='__main__':main()
