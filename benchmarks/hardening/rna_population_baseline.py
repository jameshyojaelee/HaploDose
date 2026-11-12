"""Genome-wide RNA population baseline without generator-specific information.

PCA of observed library-normalized log ratios followed by a tied-variance
two-Gaussian mixture. Empirical complete-workflow calibration is separate.
This descriptive RNA working model does not identify genomic clone identity.
"""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import norm

from haplodose.provenance import start_run, finish_run, file_record


ROLES = {'cal_nb':('matched_nb',110,199), 'cal_state':('correlated_state',111,199),
         'eval_nb':('matched_nb',120,200), 'eval_state':('correlated_state',120,200),
         'eval_two':('two_clones',120,200)}


def mixture_objective(parameters,z):
    means=np.asarray(parameters[:2]);sd=np.exp(parameters[2]);logit=parameters[3]
    log_weights=-np.logaddexp(0.,np.array([logit,-logit]))
    residual=(z-means[None])/sd
    component=log_weights[None]-.5*np.log(2*np.pi)-parameters[2]-.5*residual**2
    total=logsumexp(component,axis=1);r=np.exp(component-total[:,None])
    gradient=np.r_[np.sum(r*(z-means[None])/sd**2,axis=0),
        np.sum(r*(residual**2-1)),np.sum(r[:,1])-len(z)*expit(logit)]
    return -float(total.mean()),-gradient/len(z)


def fit_population(counts, reference):
    """Only observed counts and externally fixed reference rates enter discovery."""
    y=np.asarray(counts,float);reference=np.asarray(reference,float)
    assert y.ndim==2 and reference.shape==(y.shape[1],)
    assert np.isfinite(y).all() and (y>=0).all() and (y.sum(axis=1)>0).all()
    assert np.isfinite(reference).all() and (reference>0).all()
    # Canonical observation ordering keeps randomized starts invariant to row order.
    order=np.lexsort(y[:,::-1].T);ys=y[order]
    x=np.log((ys+.5)/(ys.sum(axis=1)[:,None]*reference[None]+.5))
    x-=x.mean(axis=0)
    _,singular,vt=np.linalg.svd(x,full_matrices=False)
    if len(singular)<2 or singular[0]<=0:raise ValueError('Unidentified principal component')
    loading=vt[0].copy()
    if loading[np.argmax(np.abs(loading))]<0:loading*=-1
    score=x@loading;variance=float(np.mean((score-score.mean())**2))
    if variance<=0:raise ValueError('Unidentified component variance')
    # Standardization is invertible and common to both likelihoods.
    z=((score-score.mean())/np.sqrt(variance))[:,None]
    ll1=float(np.sum(-.5*(np.log(2*np.pi)+z[:,0]**2)))
    # Likelihood optima have means in the observed range, and shared variance
    # between the best hard two-means residual variance and total variance.
    # These data-derived bounds replace no model parameter by a numerical ridge.
    sorted_z=np.sort(z[:,0]);n=len(z);c=np.cumsum(sorted_z);s=np.cumsum(sorted_z**2)
    split=np.arange(1,n)
    within=s[:-1]-c[:-1]**2/split+(s[-1]-s[:-1])-(c[-1]-c[:-1])**2/(n-split)
    min_variance=float(within.min()/n)
    if min_variance<=0:raise ValueError('Unbounded or numerically unidentified tied-mixture variance')
    bounds=[(float(z.min()),float(z.max()))]*2+[(.5*np.log(min_variance),0.),(None,None)]
    def objective(parameters):return mixture_objective(parameters,z)
    starts=[]
    for a,b in ((.1,.9),(.2,.8),(.25,.75),(.3,.7),(.4,.6)):
        for weight in (.3,.7):
            means=np.quantile(z[:,0],[a,b]);between=weight*(1-weight)*np.diff(means)[0]**2
            sd=np.sqrt(max(min_variance,1-between))
            starts.append(np.r_[means,np.log(sd),np.log(weight/(1-weight))])
    results=[]
    for tolerance in (1e-8,1e-10):
        candidates=[minimize(objective,start,method='L-BFGS-B',jac=True,bounds=bounds,
            options=dict(maxiter=2000,ftol=tolerance,gtol=tolerance,maxls=50)) for start in starts]
        valid=[model for model in candidates if model.success and np.isfinite(model.fun)]
        if not valid:raise ArithmeticError('No mixture optimization converged')
        model=min(valid,key=lambda result:result.fun);idx=np.argsort(model.x[:2])
        means=model.x[:2][idx];weights=np.array([expit(-model.x[3]),expit(model.x[3])])[idx]
        sd=float(np.exp(model.x[2]))
        components=np.log(weights)[None]-.5*np.log(2*np.pi*sd**2)-.5*((z-means[None])/sd)**2
        scalar=math.fsum(float(logsumexp(np.log(weights)+norm.logpdf(value,means,sd))) for value in z[:,0])
        installed=-float(model.fun)*n
        # Operation-count scaled binary64 allowance for independent summation.
        budget=64*len(z)*np.finfo(float).eps*(1+float(np.max(np.abs(components))))
        if abs(scalar-installed)>budget:raise ArithmeticError('Independent mixture likelihood disagrees')
        probability=np.exp(components-logsumexp(components,axis=1,keepdims=True))
        restored=np.empty_like(probability);restored[order]=probability
        results.append(dict(tolerance=tolerance,iterations=model.nit,converged_starts=len(valid),score=max(0.,2*(scalar-ll1)),
            means=means.tolist(),weights=weights.tolist(),sd=sd,probability=restored,
            independent_loglik_difference=abs(scalar-installed),roundoff_budget=budget))
    fine=results[-1];coarse=results[0]
    literal_changes=int(np.count_nonzero(coarse['probability'].argmax(axis=1) !=
                                         fine['probability'].argmax(axis=1)))
    return dict(score=fine['score'],coarse_score=coarse['score'],probability=fine['probability'],
        loading=loading,projection_variance=variance,singular_values=singular[:2],
        refinement_score_change=abs(fine['score']-coarse['score']),
        refinement_MAP_changes=min(literal_changes,len(y)-literal_changes),
        refinement_literal_label_changes=literal_changes,
        refinement_pairwise_disagreements=literal_changes*(len(y)-literal_changes),
        fits=[{k:v for k,v in r.items() if k!='probability'} for r in results])


def checks():
    rng=np.random.default_rng(20261109)
    y=rng.poisson(20,size=(80,16)).astype(float);y[40:,:4]+=15
    reference=np.full(16,1/16)
    a=fit_population(y,reference);permutation=rng.permutation(len(y))
    b=fit_population(y[permutation],reference)
    assert a['score']==b['score']
    assert np.array_equal(a['probability'][permutation],b['probability'])
    assert np.array_equal(a['loading'],b['loading'])
    partition=np.array([0,0,1,1,1])
    assert min(np.count_nonzero(partition!=(1-partition)),np.count_nonzero(partition==1-partition))==0
    from scipy.optimize._numdiff import approx_derivative
    z=np.linspace(-2.,2.,31)[:,None]
    for parameters in (np.array([-.7,.9,-.3,.2]),np.array([-.1,.1,0.,-1.])):
        value,gradient=mixture_objective(parameters,z)
        fd=approx_derivative(lambda x:mixture_objective(x,z)[0],parameters,method='3-point').ravel()
        allowance=64*np.finfo(float).eps**(2/3)*(1+abs(value)+np.max(np.abs(gradient)))
        assert np.max(np.abs(fd-gradient))<allowance
    return dict(cell_permutation_score_exact=True,cell_permutation_membership_exact=True,
        independent_mixture_loglik_checked=True,analytic_gradient_finite_difference_checked=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--part',choices=('preflight','calibration','evaluation'),required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parents[1];prior=root/'dev/joint_fresh_null_20261108'
    roles=() if args.part=='preflight' else tuple(r for r in ROLES if r.startswith('cal_' if args.part=='calibration' else 'eval_'))
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md')
    records={}
    for role in roles:
        records[role]=[]
        for start in (0,50,100,150):
            directory=prior/'chunks_v1/joint_full'/f'{role}_{start}'
            m=json.loads((directory/'run_manifest.json').read_text());assert m['status']=='complete'
            record=m['output']['files']['records.jsonl'];assert file_record(record['path'])['sha256']==record['sha256']
            inputs[f'{role}:{start}:manifest']=directory/'run_manifest.json'
            inputs[f'{role}:{start}:records']=Path(record['path'])
            with Path(record['path']).open() as stream:
                for line in stream:
                    previous=json.loads(line)
                    records[role].append({key:previous[key] for key in ('seed_key','input_sha256')})
    generator_path=root/'dev/recovery_next_20261107/rank_source_v1/benchmarks/hardening/global_rank_generator.py'
    inputs['generator']=generator_path
    if roles:assert file_record(generator_path)['sha256']==m['input']['files']['global_rank_generator.py']['sha256']
    spec=importlib.util.spec_from_file_location('fixed_population_generator',generator_path)
    generator=importlib.util.module_from_spec(spec);spec.loader.exec_module(generator)
    run=start_run(args.out,inputs,dict(part=args.part,n_components=[1,2],reference='fixed known reference only',
        seed=20261109,mixture_starts=10,iteration_budget=2000,refinement_tolerances=[1e-8,1e-10]),
        dict(scope='Already exposed binned simulations; RNA components not validated CNA clones',
        excluded_inputs=['true_depth','event_coordinates','CNA_effect','nuisance_distribution','labels','balanced_class_sizes']),
        mode='rna_population_baseline')
    try:
        if args.part=='preflight':
            (args.out/'checks.json').write_text(json.dumps(checks(),indent=2)+'\n')
        for role in roles:
            condition,seed_role,n=ROLES[role];assert len(records[role])==n
            for i,previous in enumerate(records[role]):
                key=[20261107,seed_role,generator.CONDITIONS.index(condition),i]
                assert previous['seed_key']==key
                y,reference,phi,truth,latent=generator.generate(np.random.default_rng(np.random.SeedSequence(key)),condition,600,240,return_latent=True)
                reset=np.zeros(240,bool);reset[::20]=True
                for name,value in dict(counts=y,reference=reference,phi=phi,truth=truth,reset=reset,
                    state_effect=latent['state_effect'],mean=latent['mean']).items():
                    assert generator.array_hash(value)==previous['input_sha256'][name],(role,i,name)
                fitted=fit_population(y,reference)
                probabilities=fitted.pop('probability');labels=probabilities.argmax(axis=1)
                arrays=args.out/f'{role}_{i:03d}.npz'
                np.savez_compressed(arrays,probability=probabilities,simulated_CNA_labels=truth,
                    loading=fitted.pop('loading'),singular_values=fitted.pop('singular_values'))
                # Truth enters only this postfit descriptive assignment score.
                correctness=max(np.mean(labels==truth),np.mean((1-labels)==truth)) if condition=='two_clones' else None
                row=dict(role=role,replicate=i,seed_key=key,**fitted,all_cell_two_group_correctness=correctness,
                    minimum_MAP_group=int(np.bincount(labels,minlength=2).min()),
                    input_sha256=previous['input_sha256'],arrays=file_record(arrays))
                with (args.out/'records.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
                if i%25==0:print(json.dumps(dict(role=role,completed=i+1,total=n)),flush=True)
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
