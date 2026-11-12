"""Deterministic reference-depth counterexample and conditional intercept repair."""
import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.optimize import minimize,root
from scipy.special import logsumexp
from haplodose.clones import SpecimenModel,depth_moments
from haplodose.provenance import start_run,finish_run


def rates(intercept,gamma,depth,origin):
    z=intercept[None]+gamma[None]*(np.log(depth)-origin)[:,None]
    return np.exp(z-logsumexp(z,axis=1,keepdims=True))


def calibrate(mu,gamma,depth,origin,tolerance):
    mu=np.asarray(mu,float);gamma=np.asarray(gamma,float);depth=np.asarray(depth,float)
    if mu.ndim!=1 or gamma.shape!=mu.shape or not np.isfinite(mu).all() or not (mu>0).all():
        raise ValueError('Positive identified reference rates and aligned slopes required')
    if not np.isfinite(depth).all() or not (depth>0).all():raise ValueError('Positive reference exposures required')
    mu=mu/mu.sum();w=depth/depth.sum();x=np.log(depth)-origin
    def evaluate(z):
        b=np.r_[z,0.];logits=b[None]+x[:,None]*gamma[None]
        p=np.exp(logits-logsumexp(logits,axis=1,keepdims=True));mean=w@p
        gradient=mean-mu
        hessian=sum(float(weight)*(np.diag(q)-q[:,None]*q[None]) for weight,q in zip(w,p))
        value=float(w@logsumexp(logits,axis=1)-mu@b)
        return value,gradient[:-1],hessian[:-1,:-1]
    start=np.log(mu[:-1]/mu[-1])
    fit=minimize(lambda z:evaluate(z)[:2],start,jac=True,method='L-BFGS-B',
        options=dict(gtol=tolerance,ftol=0.,maxiter=2000,maxls=50))
    polished=root(lambda z:evaluate(z)[1],fit.x,jac=lambda z:evaluate(z)[2],
                  options=dict(xtol=tolerance,maxfev=2000))
    b=np.r_[polished.x,0.];prediction=w@rates(b,gamma,depth,origin)
    residual=float(np.max(np.abs(prediction-mu)))
    u=np.finfo(float).eps;operations=64*len(mu)*len(depth)
    budget=operations*u/(1-operations*u)
    if residual>max(tolerance,budget):raise ArithmeticError('Reference intercept calibration unresolved')
    return b,dict(tolerance=tolerance,optimizer_success=bool(fit.success),root_success=bool(polished.success),
        reference_residual=residual,operation_budget=budget,optimizer_message=str(fit.message),root_message=str(polished.message))


def current_expected(base,gamma,depth,origin):
    model=SimpleNamespace(base=base,gamma=gamma)
    model.depth_log_rate=lambda L,XL:SpecimenModel.depth_log_rate(model,L,XL)
    moments=depth_moments(depth,np.log(depth)-origin).sum(axis=0)
    return SpecimenModel.expected(model,depth.sum(),moments)[0]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md')
    run=start_run(a.out,inputs,{},dict(scope='Deterministic known-nuisance forward model, no observed counts or biological claim'),
                  mode='depth_reference_anchor')
    try:
        baseline=np.full(40,1/40);slope=np.r_[np.full(10,-.5),np.zeros(20),np.full(10,.5)]
        base_intercept=np.log(baseline);true_origin=math.log(9000.);delta=.6*math.log(1.5)
        conditions=dict(matched=(np.array([9000.]),np.full(8,9000.),slope),
            shifted=(np.array([9000.]),np.full(8,36000.),slope),
            reference_spread=(np.array([4500.,9000.,18000.]),np.full(8,36000.),slope),
            zero_slope=(np.array([9000.]),np.full(8,36000.),np.zeros(40)))
        rows=[];solvers=[];stored={};checks=[]
        for condition,(reference,query,gamma) in conditions.items():
            ref_rates=rates(base_intercept,gamma,reference,true_origin)
            mu=(reference[:,None]*ref_rates).sum(axis=0)/reference.sum()
            origin=float(np.log(reference).mean())
            null=(query[:,None]*rates(base_intercept,gamma,query,true_origin)).sum(axis=0)
            estimates={};intercepts={}
            for tolerance in (1e-10,1e-12):
                b,record=calibrate(mu,gamma,reference,origin,tolerance)
                solvers.append(dict(condition=condition,**record));intercepts[str(tolerance)]=b
                estimates[str(tolerance)]=(query[:,None]*rates(b,gamma,query,origin)).sum(axis=0)
            b=intercepts['1e-12']
            shifted_origin=origin+math.log(2)
            b2,record=calibrate(mu,gamma,reference,shifted_origin,1e-12)
            solvers.append(dict(condition=condition,origin_shift=True,**record))
            alternate=(query[:,None]*rates(b2,gamma,query,shifted_origin)).sum(axis=0)
            expected=dict(current=current_expected(mu,gamma,query,float(np.log(query).mean())),
                origin_shift_only=current_expected(mu,gamma,query,origin),anchored=estimates['1e-12'])
            # Independently compute exact forward expectations with scalar arithmetic.
            scalar=np.zeros(40)
            for depth in query:
                weights=[float(q)*math.exp(float(g)*math.log(float(depth)/9000.)) for q,g in zip(baseline,gamma)]
                total=math.fsum(weights)
                scalar+=np.array([float(depth)*weight/total for weight in weights])
            u=np.finfo(float).eps;operations=64*40*len(query);budget=operations*u/(1-operations*u)
            assert (np.abs(scalar-null)<=budget*np.maximum(1,null)).all()
            known_b=base_intercept+gamma*(origin-true_origin);known_b-=known_b[-1]
            # Error bounds here use the known intercept, not a scientific cutoff.
            parameter_error=float(np.max(np.abs(b-known_b)))
            # For softmax, max pointwise log-rate error <=2*max intercept error.
            log_budget=2*parameter_error+budget
            assert np.max(np.abs(np.log(expected['anchored']/null)))<=log_budget
            shifted_known=known_b+gamma*math.log(2);shifted_known-=shifted_known[-1]
            shifted_error=float(np.max(np.abs(b2-shifted_known)))
            assert np.max(np.abs(np.log(alternate/null)))<=2*shifted_error+budget
            checks.append(dict(condition=condition,intercept_error=parameter_error,shifted_intercept_error=shifted_error,
                max_origin_prediction_change=float(np.max(np.abs(np.log(alternate/expected['anchored'])))),
                max_refinement_prediction_change=float(np.max(np.abs(np.log(estimates['1e-10']/estimates['1e-12'])))),
                scalar_forward_discrepancy=float(np.max(np.abs(scalar-null)))))
            for event in ('null','gain'):
                effect=np.r_[np.full(10,delta),np.zeros(30)] if event=='gain' else np.zeros(40)
                observed=(query[:,None]*rates(base_intercept+effect,gamma,query,true_origin)).sum(axis=0)
                stored[condition+'__'+event+'__mean']=observed
                for method,E in expected.items():
                    residual=np.log(observed/E);residual-=np.median(residual)
                    error=residual-effect
                    if method=='anchored':assert np.max(np.abs(error))<=2*log_budget
                    rows.append(dict(condition=condition,event=event,method=method,
                        max_abs_centred_log_error=float(np.max(np.abs(error))),
                        first_region_mean=float(residual[:10].mean()),last_region_mean=float(residual[-10:].mean())))
                    stored[condition+'__'+method+'__expected']=E
            stored[condition+'__reference_rate']=mu;stored[condition+'__reference_depth']=reference
            stored[condition+'__query_depth']=query;stored[condition+'__gamma']=gamma
            stored[condition+'__b']=b;stored[condition+'__origin']=np.array(origin)
        pd.DataFrame(rows).to_csv(a.out/'contrast_errors.tsv',sep='\t',index=False)
        np.savez_compressed(a.out/'forward_models.npz',**stored)
        (a.out/'report.json').write_text(json.dumps(dict(contrast_checks=rows,solver_records=solvers,
            numerical_checks=checks,known_gamma_only=True,biological_specimens=0,
            no_inference_default_or_scientific_criteria_change=True),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU SciPy/current installed expected() and scalar forward reference')


if __name__=='__main__':main()
