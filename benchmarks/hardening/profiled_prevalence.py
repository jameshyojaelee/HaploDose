"""Profile the concave mixture-prevalence likelihood, including its boundaries."""
from fractions import Fraction
from functools import lru_cache
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('profiled_prevalence_'+name,
        Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def prevalence(log0,log1):
    difference=np.asarray(log1)-np.asarray(log0)
    a=np.exp(-np.maximum(difference,0.))
    b=np.exp(np.minimum(difference,0.))
    # expm1 retains the contrast when two component likelihoods nearly agree.
    contrast=np.where(difference>=0.,-np.expm1(-np.maximum(difference,0.)),
                      np.expm1(np.minimum(difference,0.)))
    def score(p):
        with np.errstate(divide='ignore'):
            values=contrast/((1-p)*a+p*b)
        return math.fsum(map(float,values))
    if np.all(contrast==0):
        return .5,dict(identified=False,boundary=False,score=0.,iterations=0)
    if score(0.)<=0:
        return 0.,dict(identified=True,boundary=True,score=score(0.),iterations=0)
    if score(1.)>=0:
        return 1.,dict(identified=True,boundary=True,score=score(1.),iterations=0)
    lo,hi,iterations=0.,1.,0
    while True:
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        value=score(mid);iterations+=1
        if value==0:
            lo=hi=mid;break
        if value>0:lo=mid
        else:hi=mid
    def likelihood(p):
        return math.fsum(map(float,np.log((1-p)*a+p*b)))
    chosen=max((lo,hi),key=likelihood)
    return chosen,dict(identified=True,boundary=False,score=score(chosen),
        bracket=[lo,hi],iterations=iterations)


def objective(x,y,e,r,strata,order,*,details=False,shared=False):
    k=len(x)-2
    delta,variance=x[-2:]
    sigma=math.sqrt(variance)
    gradient,total=np.zeros_like(x),0.
    records=[]
    predictor,guard=load('state_predictive_confirmation'),load('nb_state_guard')
    components=[]
    for s in range(k):
        use=strata==s
        mu0=x[s]
        log0,_=predictor.marginal(y[use],e[use],r[use],[mu0,sigma],False,order)
        log1,_=predictor.marginal(y[use],e[use],r[use],[mu0+delta,sigma],False,order)
        components.append((use,mu0,log0,log1))
    if shared:
        pi,shared_record=prevalence(np.concatenate([c[2]for c in components]),
                                    np.concatenate([c[3]for c in components]))
    for use,mu0,log0,log1 in components:
        s=len(records)
        if shared:record=dict(shared_record)
        else:pi,record=prevalence(log0,log1)
        eta=-np.inf if pi==0 else np.inf if pi==1 else math.log(pi)-math.log1p(-pi)
        value,g=guard.variance_objective([mu0,mu0+delta,eta,variance],y[use],e[use],r[use],order,True)
        weight=use.sum()/len(y)
        gradient[s]=weight*(g[0]+g[1])
        gradient[-2]+=weight*g[1]
        gradient[-1]+=weight*g[-1]
        total+=weight*value
        records.append(dict(prevalence=pi,**record))
    if details:return float(total),gradient,records
    return float(total),gradient


def fit(D,E,r,strata,q,old,*,shared=False):
    levels,index=np.unique(strata,return_inverse=True)
    assert levels.tolist()==old['levels']
    k,order=len(levels),old['order']
    pi=np.array(old['prevalence'])
    np.testing.assert_array_equal(pi,np.array([q[index==s].mean()for s in range(k)]))
    delta=old['alternative'][-2]
    centres=np.logaddexp(np.log1p(-pi),np.log(pi)+delta)
    anchor=np.r_[np.asarray(old['alternative'][:-2])-centres,delta,old['alternative'][-1]]
    def function(x,*args,**kwargs):
        return objective(x,*args,shared=shared,**kwargs)
    anchor_value,anchor_gradient=function(anchor,D,E,r,index,order)
    operations=512*D.size*order
    eps=np.finfo(float).eps
    budget=operations*eps/(1-operations*eps)*max(1.,abs(anchor_value),abs(old['alternative_objective']))
    if not shared:assert anchor_value<=old['alternative_objective']+budget
    z=np.log((D.sum(axis=1)+.5)/(E.sum(axis=1)+.5))
    spread=float(np.diff(np.quantile(z,[.2,.8]))[0])
    null=np.asarray(old['null'])
    starts=[anchor]+[np.r_[null[:-1]-np.logaddexp(np.log1p(-pi),np.log(pi)+d),d,v]
        for d in (0.,spread,-spread)for v in (0.,null[-1])]
    options=dict(ftol=eps,gtol=np.sqrt(eps),maxiter=1000,maxls=50)
    fits=[minimize(function,x,args=(D,E,r,index,order),jac=True,
        bounds=[(None,None)]*(k+1)+[(0.,None)],method='L-BFGS-B',options=options)for x in starts]
    good=[f for f in fits if f.success and np.isfinite(f.fun)and np.isfinite(f.x).all()]
    if not good:raise ArithmeticError('Profiled prevalence optimization unresolved')
    best=min(good,key=lambda f:f.fun)
    if best.fun>anchor_value+budget:
        raise ArithmeticError('Profiled optimizer loses feasible likelihood: '+str(dict(
            old=old['alternative_objective'],profiled_anchor=anchor_value,best=best.fun,budget=budget)))
    value,gradient,records=function(best.x,D,E,r,index,order,details=True)
    projected=gradient.copy()
    if best.x[-1]==0:projected[-1]=min(0.,projected[-1])
    return dict(order=order,levels=levels.tolist(),null=old['null'],
        base_mean=best.x[:k].tolist(),delta=float(best.x[-2]),variance=float(best.x[-1]),
        prevalence=[r['prevalence']for r in records],profiled_prevalence=records,
        alternative_objective=value,fixed_prevalence_objective=old['alternative_objective'],
        profiled_anchor_objective=anchor_value,
        discovery_objective_gain=float(old['alternative_objective']-value),
        prevalence_shared_across_depth=shared,
        maximum_projected_gradient=float(np.max(abs(projected))),successful_starts=len(good),
        starts=[dict(success=bool(f.success),objective=float(f.fun),iterations=int(f.nit),
                     message=str(f.message))for f in fits],
        optimization='Exact concave prevalence profile; best successful local location/contrast/variance fit; global joint optimum not certified',
        prevalence_source='Discovery marginal likelihood profiled to adjacent floats or exact boundary')


def predict(D,C,ED,EC,rD,rC,strata,model):
    per_cell,qd,qj=[np.zeros(len(D))for _ in range(3)]
    rounding=0.
    for s,level in enumerate(model['levels']):
        use=strata==level
        pi,mean=model['prevalence'][s],model['base_mean'][s]
        sigma=math.sqrt(model['variance'])
        null=[model['null'][s],math.sqrt(model['null'][-1])]
        if pi in (0.,1.):
            mu=mean+pi*model['delta']
            alternative=[mu,mu,0.,sigma]
        else:
            alternative=[mean,mean+model['delta'],math.log(pi)-math.log1p(-pi),sigma]
        result=load('state_predictive_confirmation').predictive(
            D[use],C[use],ED[use],EC[use],rD[use],rC[use],null,alternative,model['order'])
        per_cell[use]=result['per_cell']
        qd[use]=pi if pi in (0.,1.)else result['mixture_membership_discovery']
        qj[use]=pi if pi in (0.,1.)else result['mixture_membership_joint']
        rounding+=result['rounding_budget']
    def effective(x):return float(x.sum()**2/np.sum(x*x))if np.any(x)else 0.
    return dict(score=float(math.fsum(map(float,per_cell))),rounding_budget=rounding,
        per_cell=per_cell,discovery_membership=qd,joint_membership=qj,
        discovery_expected_members=[float((1-qd).sum()),float(qd.sum())],
        discovery_effective_members=[effective(1-qd),effective(qd)])


def checks(*,shared=False):
    a,b=np.array([.75,.25,.75]),np.array([.25,.75,.25])
    pi,record=prevalence(np.log(a),np.log(b))
    # Exact solution of -2/(3/2-p)+1/(1/2+p)=0 is p=1/6.
    error=abs(pi-float(Fraction(1,6)))
    eps=np.finfo(float).eps
    budget=256*len(a)*eps/(1-256*len(a)*eps)
    assert error<=budget and record['bracket'][1]-record['bracket'][0]<=np.spacing(pi)
    assert prevalence(np.log(a),np.log(a))[1]['identified']is False
    assert prevalence(np.log(a),np.log(a/2))[0]==0
    assert prevalence(np.log(a/2),np.log(a))[0]==1
    reverse,_=prevalence(np.log(b),np.log(a))
    assert abs(reverse-(1-pi))<=budget
    # Envelope gradient checked against the independently verified original
    # full-PMF parameterization, with prevalence fixed at its profiled value.
    D=np.array([[0,0],[4,6],[0,1],[5,7]])
    E=np.tile([1.,2.],(4,1));r=np.full(D.shape,2.)
    strata=np.array([0,0,1,1]);x=np.array([-.4,-.3,.7,.04])
    value,g,records=objective(x,D,E,r,strata,96,details=True,shared=shared)
    p=np.array([rec['prevalence']for rec in records])
    assert np.all((p>0)&(p<1))
    eta=np.log(p)-np.log1p(-p)
    centre=np.logaddexp(np.log1p(-p),np.log(p)+x[-2])
    old_x=np.r_[x[:2]+centre,x[-2],eta,x[-1]]
    old_value,old_g=load('likelihood_prevalence').objective(old_x,D,E,r,strata,96)
    transformed=np.r_[old_g[:2],old_g[2]+np.dot(expit(eta+x[-2]),old_g[:2]),old_g[-1]]
    # At the profiled solution the nuisance score is zero to rounding.
    reference=load('likelihood_prevalence').checks()
    allowance=reference['quad_and_rounding_budget']
    assert abs(value-old_value)<=allowance and np.max(abs(g-transformed))<=allowance
    nuisance=old_g[3:5]+(expit(eta+x[-2])-p)*old_g[:2]
    nuisance_error=abs(nuisance.sum())if shared else np.max(abs(nuisance))
    assert nuisance_error<=allowance
    # Independent five-point differences across a convergence sequence.
    differences=[]
    for step in (2.**-9,2.**-11,2.**-13,2.**-15):
        gradient=np.empty_like(x)
        for j in range(len(x)):
            move=np.zeros_like(x);move[j]=step
            f=[objective(x+k*move,D,E,r,strata,96,shared=shared)[0]for k in (-2,-1,1,2)]
            gradient[j]=(f[0]-8*f[1]+8*f[2]-f[3])/(12*step)
        differences.append(dict(step=step,error=float(np.max(abs(gradient-g)))))
    # A fixed independent full-PMF reference sets the gradient allowance; the
    # finite-difference sequence additionally must converge towards it.
    assert differences[-1]['error']<=allowance
    return dict(exact_fraction_prevalence='1/6',exact_prevalence_error=error,
        adjacent_float_or_exact_root=True,boundaries=True,flat_unidentified=True,
        component_swap=True,profiled_full_PMF_reference_error=float(np.max(abs(g-transformed))),
        fixed_component_prevalence_score_error=float(nuisance_error),
        prevalence_shared_across_depth=shared,
        full_PMF_reference_allowance=allowance,five_point_differences=differences,
        unprofiled_independent_checks=reference,scientific_thresholds_changed=False)
