"""Fit-only reference-rate candidates, then heldout scoring and truth reveal."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logsumexp
from scipy.stats import nbinom
from haplodose.provenance import start_run, finish_run, file_record


def objective(z,Y,L,r,initial):
    logits=np.r_[z,0.];lograte=logits-logsumexp(logits)
    logmean=np.log(L)[:,None]+lograte[None]
    normalizer=np.logaddexp(np.log(r)[None],logmean)
    base=np.log(L)[:,None]+np.log(initial)[None]
    base_normalizer=np.logaddexp(np.log(r)[None],base)
    scale=Y.sum()
    value=np.sum((Y+r)*(normalizer-base_normalizer)-Y*(logmean-base))/scale
    g=((Y+r)*np.exp(logmean-normalizer)-Y).sum(axis=0)
    gradient=(g-np.exp(lograte)*g.sum())/scale
    return float(value),gradient[:-1]


def scalar_complex_objective(z,Y,L,r,initial):
    weights=np.exp(np.r_[z,0.]);rate=weights/weights.sum();total=0j
    for i,library in enumerate(L):
        for j,q in enumerate(rate):
            mean=library*q;old=library*initial[j]
            total+=(Y[i,j]+r[j])*(np.log(r[j]+mean)-np.log(r[j]+old))-Y[i,j]*(np.log(mean)-np.log(old))
    return total/Y.sum()


def fit_rates(counts):
    Y=np.asarray(counts,float)
    if Y.ndim!=2 or not np.isfinite(Y).all() or (Y<0).any() or not np.equal(Y,np.round(Y)).all():
        raise ValueError('Finite nonnegative observed count matrix required')
    Y=Y[sorted(range(len(Y)),key=lambda i:tuple(Y[i]))]
    L=Y.sum(axis=1)
    if not (L>0).all():raise ValueError('Positive reference library totals required')
    pooled=np.maximum(Y.sum(axis=0),.1);pooled/=pooled.sum()
    equal=np.mean(Y/L[:,None],axis=0);equal/=equal.sum()
    if not (equal>0).all():raise ValueError('Zero equal-cell rate; no invented floor')
    mean=L[:,None]*pooled
    phi=np.clip(np.sum((Y-mean)**2-mean,axis=0)/np.sum(mean**2,axis=0),1e-3,10.)
    r=1/phi;fits=[];rates={};checks=[]
    for tolerance in (1e-8,1e-10):
        candidates=[]
        for label,start in [('pooled',pooled),('equal',equal)]:
            z=np.log(start[:-1]/start[-1])
            # Complex-step references avoid subtractive finite-difference cancellation.
            actual=objective(z,Y,L,r,pooled)[1]
            references=[]
            for step in (1e-15,1e-25):
                derivative=[]
                for j in range(len(z)):
                    point=z.astype(complex);point[j]+=step*1j
                    derivative.append(scalar_complex_objective(point,Y,L,r,pooled).imag/step)
                references.append(np.array(derivative))
            eps=np.finfo(float).eps;operations=32*Y.size
            budget=operations*eps/(1-operations*eps)*max(1.,float(np.max(np.abs(actual))))
            discrepancy=max(float(np.max(np.abs(actual-v))) for v in references)
            assert discrepancy<=budget
            checks.append(dict(start=label,tolerance=tolerance,gradient_discrepancy=discrepancy,operation_budget=budget))
            result=minimize(objective,z,args=(Y,L,r,pooled),jac=True,method='L-BFGS-B',
                options=dict(maxiter=2000,gtol=tolerance,ftol=0.,maxls=50))
            q=np.exp(np.r_[result.x,0.]-logsumexp(np.r_[result.x,0.]))
            assert np.isfinite(q).all() and (q>0).all()
            # Independently constructed full NB log-PMF difference.
            new=nbinom.logpmf(Y,r,r/(r+L[:,None]*q))
            old=nbinom.logpmf(Y,r,r/(r+L[:,None]*pooled))
            independent=float((old-new).sum()/Y.sum())
            allowance=operations*eps/(1-operations*eps)*max(1.,float((np.abs(old)+np.abs(new)).sum()/Y.sum()))
            assert abs(independent-result.fun)<=allowance
            record=dict(tolerance=tolerance,start=label,success=bool(result.success),message=str(result.message),
                iterations=int(result.nit),objective=float(result.fun),scaled_gradient_max=float(np.max(np.abs(result.jac))),
                full_PMF_discrepancy=abs(independent-result.fun),full_PMF_budget=allowance)
            fits.append(record);candidates.append((result,q))
        successful=[(result,q) for result,q in candidates if result.success]
        if not successful:raise ArithmeticError('No converged NB reference fit')
        result,q=min(successful,key=lambda pair:pair[0].fun)
        if result.fun>allowance:raise ArithmeticError('Selected fit worse than pooled starting rate')
        rates[str(tolerance)]=q
    return dict(pooled=pooled,equal=equal,nb=rates['1e-10']),phi,rates,fits,checks


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--means',type=Path,help='Explicit generating-mean source for information-leakage challenges')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    prior=a.task.parent/'profile_pool_export_20261109/parity_v1/fitted_new'
    truth_path=a.means or a.task.parent/'profile_pool_mechanism_20261109/panel_v1/primary_means_and_profiles.npz'
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',cells=prior/'cells.tsv.gz',counts=prior/'bin_counts.npy',
                means=truth_path,prior_audit=a.task.parent/'profile_pool_mechanism_20261109/audit_v1/run_manifest.json')
    assert json.loads(inputs['prior_audit'].read_text())['status']=='complete'
    run=start_run(a.out,inputs,dict(methods=['pooled','equal','nb']),
        dict(scope='Exposed fit/holdout RNA reference prediction; no DNA or calibration claim'),mode='reference_rate_prediction')
    try:
        equal_library=np.array([[2.,3.,5.],[5.,2.,3.],[3.,5.,2.]])
        analytic,_,_,_,_=fit_rates(equal_library)
        assert np.max(np.abs(analytic['nb']-np.full(3,1/3)))<=32*np.finfo(float).eps
        meta=pd.read_csv(inputs['cells'],sep='\t');Y=np.load(inputs['counts'])
        groups=sorted(meta.loc[meta.role.eq('reference'),'ref_group'].unique())
        predictions={};scores=[];optimizations=[];checks=[];selection=[]
        for group in groups:
            fit=meta.role.eq('reference')&meta.ref_group.eq(group)&meta.ref_split.eq('fit')
            hold=meta.role.eq('reference')&meta.ref_group.eq(group)&meta.ref_split.eq('holdout')
            assert fit.sum()==hold.sum()==30
            rates,phi,refinements,records,gradient_checks=fit_rates(Y[fit])
            reordered=fit_rates(Y[fit][::-1])[0]
            assert all(np.array_equal(value,reordered[method]) for method,value in rates.items())
            predictions[group+'__phi']=phi;predictions[group+'__fit']=fit.to_numpy();predictions[group+'__holdout']=hold.to_numpy()
            LH=Y[hold].sum(axis=1);r=1/phi
            for method,q in rates.items():
                predictions[group+'__'+method]=q
                values=nbinom.logpmf(Y[hold],r,r/(r+LH[:,None]*q))
                scores.append(dict(group=group,method=method,fit_cells=int(fit.sum()),holdout_cells=int(hold.sum()),
                    holdout_composite_log_score=float(values.sum()),fit_minimum_phi=int(np.sum(phi==1e-3)),
                    fit_maximum_phi=int(np.sum(phi==10))))
            for key,q in refinements.items():predictions[group+'__nb_'+key]=q
            optimizations.extend([dict(group=group,**record) for record in records])
            checks.extend([dict(group=group,**record) for record in gradient_checks])
            rankings=[]
            for tolerance in ('1e-08','1e-10'):
                alternative={**rates,'nb':refinements[str(float(tolerance))]}
                ss={method:float(nbinom.logpmf(Y[hold],r,r/(r+LH[:,None]*q)).sum()) for method,q in alternative.items()}
                rankings.append(sorted(ss,key=lambda method:-ss[method]))
            selection.append(dict(group=group,coarse_ranking=rankings[0],fine_ranking=rankings[1],
                refinement_choice_stable=rankings[0][0]==rankings[1][0]))
        # Freeze all predictions and selection before any generating truth is read.
        with (a.out/'predictions.npz').open('xb') as stream:np.savez_compressed(stream,**predictions)
        pd.DataFrame(scores).to_csv(a.out/'heldout_scores.tsv',sep='\t',index=False)
        (a.out/'selection_before_truth.json').write_text(json.dumps(selection,indent=2)+'\n')
        frozen={name:file_record(a.out/name) for name in ('predictions.npz','heldout_scores.tsv','selection_before_truth.json')}
        (a.out/'prediction_freeze.json').write_text(json.dumps(frozen,indent=2)+'\n')
        with np.load(truth_path) as z:M=z['used_gene_mean']
        errors=[]
        for group in groups:
            for split in ('fit','holdout'):
                mask=predictions[group+'__'+split];target=M[mask].sum(axis=0);target/=target.sum()
                for method in ('pooled','equal','nb'):
                    q=predictions[group+'__'+method]
                    errors.append(dict(group=group,target_split=split,method=method,
                        centred_log_rate_MSE=float(np.mean(((np.log(q)-np.median(np.log(q)))-
                            (np.log(target)-np.median(np.log(target))))**2))))
        pd.DataFrame(errors).to_csv(a.out/'post_freeze_truth_errors.tsv',sep='\t',index=False)
        for record in frozen.values():assert file_record(record['path'])['sha256']==record['sha256']
        (a.out/'optimizations.json').write_text(json.dumps(optimizations,indent=2)+'\n')
        (a.out/'numerical_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
        (a.out/'report.json').write_text(json.dumps(dict(biological_specimens=0,conditional_fixture=1,
            reference_types=groups,scores=scores,selection=selection,truth_errors=errors,
            predictions_frozen_before_truth=True,production_changed=False,scientific_criteria_changed=False),indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out,run,error=error);raise
    finish_run(a.out,run,backend='CPU NumPy/SciPy, no clone workflow refit')


if __name__=='__main__':main()
