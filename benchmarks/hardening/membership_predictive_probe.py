"""Conditional-likelihood diagnostic on saved proposal fits, never an acceptance gate.

A latent component is integrated once per cell over its entire region. The
proposal probabilities, rates and exposure are frozen. Null rate alone is fitted
on confirmation counts, giving a conservative denominator for a known-nuisance
conditional e-value. Here nuisance and thinning sizes are fitted and reused, so
these values DO NOT carry an e-value calibration claim.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit, logsumexp
from scipy.stats import nbinom


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--trace-root',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    rows=[]
    for seed in [20261108,20261109,20261110]:
        root=a.trace_root/f'trace_{seed}'
        for event in json.loads((root/'trace.json').read_text()):
            if event['stage']!='event_fit':continue
            z=np.load(root/event['arrays']);O,E,r,q=(z[k] for k in ['Ytest','E','r','q'])
            theta=np.exp(z['parameters'][:2]);g=q>.5
            def ll(t):return nbinom.logpmf(O,r,r/(r+E*t)).sum(axis=1)
            logE=np.log(E);logr=np.log(r)
            def score(t):return float((O-(O+r)*expit(logE+t-logr)).sum())
            lo,hi=-1.,1.
            while score(lo)<=0:lo*=2
            while score(hi)>=0:hi*=2
            eta=brentq(score,lo,hi,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
            null=ll(np.exp(eta));components=np.stack([ll(theta[0]),ll(theta[1])],axis=1)
            with np.errstate(divide='ignore'):
                prediction=logsumexp(components+np.stack([np.log1p(-q),np.log(q)],axis=1),axis=1)
            hard=components[np.arange(len(q)),g.astype(int)]
            weights=E.sum(axis=1)
            purity=np.array([np.average(q[g==j],weights=weights[g==j]) for j in [0,1]])
            expected_log_pool_contrast=np.log((1-purity[1])*theta[0]+purity[1]*theta[1])-np.log((1-purity[0])*theta[0]+purity[0]*theta[1])
            rows.append(dict(seed=seed,trace_index=event['index'],proposal_gain=event['gain'],
                proposal_latent_log_contrast=float(np.log(theta[1]/theta[0])),
                expected_assigned_pool_log_contrast=float(expected_log_pool_contrast),
                hard_pool_purities=purity.tolist(),
                conditional_predictive_log_gain=float((prediction-null).sum()),
                hard_assignment_predictive_log_gain=float((hard-null).sum()),
                scope='Fitted-nuisance diagnostic, not a calibrated test or new pipeline output'))
    (a.out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows))


if __name__=='__main__':main()
