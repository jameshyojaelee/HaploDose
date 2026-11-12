"""Research continuous RNA-contrast HMM for robust regional boundaries."""
import argparse
from functools import lru_cache
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp, ndtri
from scipy.stats import norm
from haplodose import bulk, profile
from haplodose.io import write_json


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('reference_gaussian_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def bin_observations(effects,expected,variance,gene_bins):
    y,e,v,bins=map(np.asarray,(effects,expected,variance,gene_bins))
    if (y.ndim!=1 or not len(y) or any(x.shape!=y.shape for x in (e,v,bins))
            or bins.dtype.kind not in 'iu' or not np.array_equal(np.unique(bins),np.arange(int(bins.max())+1))
            or not all(np.isfinite(x).all() for x in (y,e,v)) or np.any(e<0) or np.any(v<=0)):
        raise ValueError('Finite aligned effects, expectations, positive scales and complete bin IDs required')
    median=load('reference_weighted_regions').weighted_median
    values=[];sampling=[]
    for j in range(int(bins.max())+1):
        keep=bins==j;a=1/np.sqrt(v[keep]);density=a*np.sqrt(e[keep]+.5)
        values.append(median(y[keep],a))
        sampling.append((np.pi/2)*math.fsum(map(float,a*a))/math.fsum(map(float,density))**2)
    return np.asarray(values),np.asarray(sampling)


def train(normal_effects,normal_expected,gene_bins):
    effects,expected=np.asarray(normal_effects,float),np.asarray(normal_expected,float)
    if effects.ndim!=2 or expected.shape!=effects.shape or len(effects)<3:
        raise ValueError('At least three aligned normal profiles required')
    weighted=load('reference_weighted_regions');rows=[]
    for i in range(len(effects)):
        other=np.delete(effects,i,axis=0)
        adjusted,centre,_=weighted.normal_centered_effect(effects[i],other,gene_bins)
        spread=(np.median(abs(other-centre),axis=0)/ndtri(.75))**2
        variance=weighted.gene_variance(spread,expected[i],'reference')
        rows.append(bin_observations(adjusted,expected[i],variance,gene_bins)[0])
    rows=np.asarray(rows);location=np.median(rows,axis=0)
    variance=(np.median(abs(rows-location),axis=0)/ndtri(.75))**2
    return dict(normal_bins=rows,reference_variance=variance,
        status='Conditional working spread; no calibrated sampling or DNA-state probabilities')


def chain_blocks(emission,reset,prior,change):
    emission,reset,prior=np.asarray(emission,float),np.asarray(reset,bool),np.asarray(prior,float)
    if emission.ndim!=2 or reset.shape!=(len(emission),) or not len(reset) or not reset[0]:
        raise ValueError('Aligned emissions and genomic resets required')
    engine=load('shared_membership_hmm');prior=prior/prior.sum()
    logT=np.log(engine.transition(prior,change));initial=np.log(prior)
    starts=np.flatnonzero(reset);stops=np.r_[starts[1:],len(reset)]
    post=np.empty_like(emission);path=np.empty(len(reset),int);scores=[]
    for a,b in zip(starts,stops):
        post[a:b],value,path[a:b]=engine.chain(emission[a:b],logT,initial,np.zeros(len(prior)))
        scores.append(value)
    return post,math.fsum(scores),path


def fit(effects,expected,gene_variance,gene_bins,reference_variance,reset,params,*,centre_mode='zero'):
    if centre_mode not in ('zero','bin-median'):raise ValueError('Unsupported boundary intercept')
    values,sampling=bin_observations(effects,expected,gene_variance,gene_bins)
    variance=np.asarray(reference_variance,float)+sampling;reset=np.asarray(reset,bool)
    if (variance.shape!=values.shape or reset.shape!=values.shape or not reset[0]
            or not np.isfinite(variance).all() or np.any(variance<=0)):
        raise ValueError('Finite positive bin scales and complete genomic reset correspondence required')
    inside=~reset[1:];difference=np.diff(values)[inside]
    if len(difference)<10:tau=params.tau_floor
    else:
        mad=np.median(abs(difference-np.median(difference)))/ndtri(.75)
        tau=math.sqrt(max(params.tau_floor**2,(mad**2-float(np.median((variance[:-1]+variance[1:])[inside])))/2))
    variance=variance+tau**2;fits=[]
    centre=0. if centre_mode=='zero' else float(np.median(values))
    for beta in bulk.BETA_GRID:
        if not params.beta_bounds[0]<=beta<=params.beta_bounds[1]:continue
        means=beta*np.asarray(profile.LOG_DOSE,float)
        emission=-.5*(np.log(2*np.pi*variance)[:,None]+(values[:,None]-centre-means)**2/variance[:,None])
        post,score,path=chain_blocks(emission,reset,profile.STATE_PRIOR,params.p_change)
        fits.append(dict(beta=beta,score=score,path=path,posterior=post))
    if not fits:raise ValueError('No existing bulk beta grid value satisfies the declared bounds')
    best=max(range(len(fits)),key=lambda j:fits[j]['score']);selected=fits[best]
    operations=128*len(values)*len(profile.CN)+256
    allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,max(abs(x['score']) for x in fits))
    near=[j for j,x in enumerate(fits) if abs(x['score']-selected['score'])<=2*allowance]
    return dict(path=selected['path'],posterior=selected['posterior'],bin_values=values,
        sampling_variance=sampling,total_variance=variance,beta=selected['beta'],tau=tau,centre=centre,centre_mode=centre_mode,
        grid_scores=[dict(beta=x['beta'],score=x['score']) for x in fits],
        near_tied_grid_indices=near if len(near)>1 else [],score_rounding_budget=allowance,
        calibration='Unvalidated conditional RNA boundary model')


def checks():
    prior=np.array([.2,.5,.2,.1]);change=.02;levels=np.array([-.4,0.,.2,.4])
    y=np.array([-.3,.05,.3]);variance=np.array([.2,.3,.4])
    emission=-.5*(np.log(2*np.pi*variance)[:,None]+(y[:,None]-levels)**2/variance[:,None])
    post,evidence,path=chain_blocks(emission,[True,False,False],prior,change)
    states=list(itertools.product(range(4),repeat=3));joint=[]
    for state in states:
        value=math.log(prior[state[0]])
        for j in range(3):
            value+=float(norm.logpdf(y[j],levels[state[j]],math.sqrt(variance[j])))
            if j:value+=math.log(1-change if state[j-1]==state[j] else change*prior[state[j]]/(1-prior[state[j-1]]))
        joint.append(value)
    independent=float(logsumexp(joint));weights=np.exp(np.asarray(joint)-independent);expected=np.zeros_like(post)
    for state,weight in zip(states,weights):
        for j,s in enumerate(state):expected[j,s]+=weight
    operations=128*len(states)*len(y)+128
    allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(independent))
    assert abs(evidence-independent)<=allowance and np.max(abs(post-expected))<=allowance
    np.testing.assert_array_equal(path,states[int(np.argmax(joint))])
    split=chain_blocks(emission,[True,False,True],prior,change)
    first=chain_blocks(emission[:2],[True,False],prior,change);last=chain_blocks(emission[2:],[True],prior,change)
    np.testing.assert_array_equal(split[0],np.vstack([first[0],last[0]]));np.testing.assert_array_equal(split[2],np.r_[first[2],last[2]])
    assert split[1]==math.fsum([first[1],last[1]])
    values,sampling=bin_observations([0.,2.,4.,8.],[3.5]*4,[1.]*4,np.array([0,0,1,1]))
    np.testing.assert_array_equal(values,[1.,6.])
    scalar=(math.pi/2)*(1.+1.)/(2.+2.)**2
    np.testing.assert_array_equal(sampling,[scalar,scalar])
    normals=np.array([[0.,1.,2.,3.],[1.,2.,4.,5.],[2.,3.,6.,7.]])
    training=train(normals,np.full_like(normals,3.5),np.array([0,0,1,1]))
    scalar_normal=[]
    weighted=load('reference_weighted_regions')
    for i in range(3):
        other=normals[[j for j in range(3) if j!=i]]
        centred,_,_=weighted.normal_centered_effect(normals[i],other,np.array([0,0,1,1]))
        # Within each two-gene bin the two scales are equal in this fixture.
        scalar_normal.append([float(np.median(centred[:2])),float(np.median(centred[2:]))])
    np.testing.assert_array_equal(training['normal_bins'],scalar_normal)
    permutation=np.array([3,0,2,1]);alternate=train(normals[:,permutation],np.full_like(normals,3.5)[:,permutation],np.array([0,0,1,1])[permutation])
    for key in ('normal_bins','reference_variance'):np.testing.assert_array_equal(training[key],alternate[key])
    return dict(independent_scipy_path_enumeration=True,maximum_posterior_error=float(np.max(abs(post-expected))),
        enumeration_allowance=allowance,independent_sampling_variance=True,normal_loo_identity=True,
        gene_permutation=True,gap_blocks_exact=True,working_variance_only=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--reference-refit',type=int,choices=range(9),default=0)
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:
        load('reference_weighted_regions').panel(argparse.Namespace(out=args.out,task=args.task,
            reference_refit=args.reference_refit,comparison_root=args.task.parent/'normal_centered_regions_20261109',
            weighting='reference',normal_centering=True,reference_gaussian_regions=True))
        if args.reference_refit==0 and args.out.name!='bulk_v1':
            original=args.task/'bulk_v1'
            with np.load(original/'RNA_predictions.npz') as old,np.load(args.out/'RNA_predictions.npz') as new:
                assert old.files==new.files
                for key in old.files:np.testing.assert_array_equal(old[key],new[key])
            pd.testing.assert_frame_equal(pd.read_csv(original/'patient_scores.tsv',sep='\t',float_precision='round_trip'),
                pd.read_csv(args.out/'patient_scores.tsv',sep='\t',float_precision='round_trip'),check_exact=True)
            write_json(dict(primary_prediction_and_all30_scores_exact=True),args.out/'refit_control.json')


if __name__=='__main__':main()
