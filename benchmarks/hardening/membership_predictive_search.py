"""Experimental whole-search statistic integrating uncertain proposal memberships.

No clone labels or pipeline decisions are changed. Every eligible event mixture
encountered by the executed fit contributes a confirmation predictive likelihood
relative to its profiled single-rate NB null. The statistic is their maximum.
Only complete-workflow simulation calibrates this statistic here; the fitted
nuisance quantities and adaptive reuse prohibit an independent e-value claim.
"""
from contextlib import contextmanager
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp
from haplodose import clones


def conditional_gain(O, E, r, probability, log_rates, *, return_hard=False):
    O, E, r = np.broadcast_arrays(np.asarray(O,float),np.asarray(E,float),np.asarray(r,float))
    q=np.asarray(probability,float)
    if O.ndim!=2 or q.shape!=(len(O),):raise ValueError('Expected cells x bins and per-cell probabilities')
    if not (np.isfinite(O).all() and np.isfinite(E).all() and np.isfinite(r).all()
            and np.isfinite(q).all() and np.isfinite(log_rates).all()):raise ValueError('Nonfinite inputs')
    if np.any(O<0) or np.any(E<=0) or np.any(r<=0) or np.any((q<0)|(q>1)):
        raise ValueError('Invalid count, exposure, size or probability')
    if not O.sum():return (0., 0.) if return_hard else 0.
    logE=np.log(E);logr=np.log(r)
    def score(x):return float((O-(O+r)*expit(logE+x-logr)).sum())
    lo,hi=-1.,1.
    while score(lo)<=0:lo*=2
    while score(hi)>=0:hi*=2
    eta=brentq(score,lo,hi,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
    base=logE+eta;values=[]
    for x in log_rates:
        d=x-eta
        term=(np.log1p(expit(base-logr)*np.expm1(d)) if abs(d)<.5
              else np.logaddexp(logr,base+d)-np.logaddexp(logr,base))
        values.append((O*d-(O+r)*term).sum(axis=1))
    with np.errstate(divide='ignore'):
        joint=np.stack(values)+np.stack([np.log1p(-q),np.log(q)])
    integrated=float(logsumexp(joint,axis=0).sum())
    if return_hard:
        hard=float(np.stack(values)[(q>.5).astype(int),np.arange(len(q))].sum())
        return integrated,hard
    return integrated


@contextmanager
def capture_search():
    if sys.getprofile() is not None:raise RuntimeError('Existing profiler')
    collected=[]
    def observe(frame,event,result):
        if event!='return' or frame.f_code is not clones._fit_event_mixture.__code__ or result is None:return
        outer=frame.f_back
        if outer.f_code is not clones.SpecimenModel.event_mixtures.__code__:return
        loc=outer.f_locals;m=loc['self']
        if m.graph is not None or m.alleles is not None or m.genes is not None:
            raise ValueError('Experimental statistic only supports nonspatial binned expression')
        q=result['responsibility'];groups=q>.5
        if result['gain']<=result.get('numerical_gain_budget',0.) or min(groups.sum(),(~groups).sum())<m.params.min_cells:return
        cols=loc['cols'][loc['part']]
        test=np.flatnonzero(m.test_mask & (np.arange(m.J)>=cols[0]) & (np.arange(m.J)<=cols[-1]))
        if len(test)<2:return
        A=loc['lr'][None]+loc['xc'][:,None]*m.gamma[None]+loc['uw']
        scale=np.log(loc['Lp'])+np.log(np.maximum(loc['s'],1e-12))-logsumexp(A[:,loc['cols']],axis=1)
        E=np.exp(scale[:,None]+A[:,test])
        O=m.rows.take(loc['cells'],m.test_fold)[:,test]
        share=1. if m.test_fold is None else m.fold_scale[m.test_fold]
        r=share/np.maximum(m.phi[test],1e-3)
        value,hard=conditional_gain(O,E,r,q,result['parameters'][:2],return_hard=True)
        if not np.isfinite([value,hard]).all():raise FloatingPointError('Nonfinite predictive statistic')
        collected.append(dict(score=value,hard_score=hard,first_bin=int(test[0]),last_bin=int(test[-1]),
            cells=len(q),component_sizes=[int((~groups).sum()),int(groups.sum())]))
    sys.setprofile(observe)
    try:yield collected
    finally:sys.setprofile(None)


def instrument_rank(module):
    original=module.run_specimen
    def run(*args,**kwargs):
        with capture_search() as records:
            row=original(*args,**kwargs)
        row['original_final_difference_score']=row['score']
        row['score']=max([0.]+[r['score'] for r in records])
        row['hard_score']=max([0.]+[r['hard_score'] for r in records])
        row['conditional_proposal_count']=len(records)
        row['conditional_predictive_diagnostics']=records
        row['statistic_scope']='Maximum over executed event proposals; RNA-model evidence only; no independent confirmation or membership validation'
        return row
    module.run_specimen=run
