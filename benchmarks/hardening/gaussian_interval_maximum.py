"""Research Gaussian all-interval tolerance bound; no automatic package use.

Known covariance, fixed eligible intervals and independent calibration randomness
are essential. Statistical error allocation includes calibration failure.
"""
import math
import numpy as np
from scipy import stats


def tolerance_rank(n, tail, failure):
    """Smallest k with P[Bin(n,1-tail) >= k] <= failure, one-based rank.

    F(M_(k)) ~ Beta(k,n+1-k) for continuous maxima. For calibration failure
    probability <= failure, P(new M > M_(k) | bank) <= tail. A larger computed
    maximum is conservative. No fitted quantile approximation is used.
    """
    if n<1 or not 0<tail<1 or not 0<failure<1:raise ValueError('Invalid tolerance design')
    left,right=1,n+1
    while left<right:
        mid=(left+right)//2
        if stats.binom.sf(mid-1,n,1-tail)<=failure:right=mid
        else:left=mid+1
    if left>n:raise ValueError('Calibration sample too small for finite tolerance bound')
    return left


def geometry(covariance,reset):
    S=np.asarray(covariance,float);reset=np.asarray(reset)
    J=len(reset)
    if reset.dtype.kind!='b' or not J or not reset[0] or S.shape!=(J,J):raise ValueError('Invalid interval geometry')
    if not np.isfinite(S).all() or not np.array_equal(S,S.T):raise ValueError('Finite symmetric covariance required')
    boundaries=np.r_[np.flatnonzero(reset),J];first=[];last=[];lower=[];upper=[]
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        for start in range(a,b):
            for stop in range(start+1,b+1):
                values=S[start:stop,start:stop].ravel().tolist();terms=len(values)
                total=math.fsum(values);gamma=terms*np.finfo(float).eps/(1-terms*np.finfo(float).eps)
                error=gamma*math.fsum(abs(x) for x in values)
                if total-error<=0:raise ValueError('Zero/unresolved interval variance: no automatic PSD repair')
                first.append(start);last.append(stop);lower.append(math.sqrt(total-error));upper.append(math.sqrt(total+error))
    return dict(first=np.array(first),stop=np.array(last),sum_sd_lower=np.array(lower),sum_sd_upper=np.array(upper))


def maxima(Z,g):
    """Upper-rounded standardized sum maxima over every eligible interval."""
    Z=np.asarray(Z,float);J=Z.shape[1];eps=np.finfo(float).eps
    prefix=np.column_stack([np.zeros(len(Z)),np.cumsum(Z,axis=1)])
    absolute=np.column_stack([np.zeros(len(Z)),np.cumsum(abs(Z),axis=1)])
    gamma=J*eps/(1-J*eps)
    first,stop=g['first'],g['stop']
    sums=prefix[:,stop]-prefix[:,first]
    # Prefix accumulation, subtraction, division and final comparison roundoff.
    error=gamma*(absolute[:,stop]+absolute[:,first])+eps*(abs(prefix[:,stop])+abs(prefix[:,first]))
    value=(abs(sums)+error)/g['sum_sd_lower']
    return np.nextafter(value.max(axis=1)*(1+2*eps),np.inf)


def calibrate(S,reset,*,n,seed,tail,failure,batch=512):
    g=geometry(S,reset);L=np.linalg.cholesky(S)
    eps=np.finfo(float).eps;J=len(reset);gamma=(2*J)*eps/(1-2*J*eps)
    residual=float(np.max(abs(L@L.T-S)))
    bound=float(np.max(gamma*(abs(L)@abs(L).T+abs(S))))
    if residual>bound:raise ArithmeticError('Covariance factor reconstruction unresolved')
    rng=np.random.default_rng(seed);values=[]
    for start in range(0,n,batch):
        Z=rng.standard_normal((min(batch,n-start),J))@L.T
        values.append(maxima(Z,g))
    values=np.concatenate(values);k=tolerance_rank(n,tail,failure)
    critical=float(np.partition(values,k-1)[k-1])
    report=dict(n=n,rank=k,tail=tail,calibration_failure_budget=failure,
                exact_binomial_failure=float(stats.binom.sf(k-1,n,1-tail)),
                previous_rank_failure=float(stats.binom.sf(k-2,n,1-tail)),
                beta_reference_failure=float(stats.beta.cdf(1-tail,k,n+1-k)),
                conditional_tail_upper=float(1-stats.beta.ppf(failure,k,n+1-k)),
                critical=critical,family_size=len(g['first']),seed=seed,
                factor_residual=residual,factor_error_allowance=bound,
                numerical_scope='Float64 covariance-factor and sum checks; theoretical guarantee assumes exact Gaussian calibration law')
    return report,values,g


def checks():
    # Exhaustive scalar sums differ from prefix evaluation; reset excludes [1,2].
    S=np.array([[1.,.2,0.],[.2,2.,0.],[0.,0.,.5]])
    reset=np.array([True,False,True]);g=geometry(S,reset)
    if list(zip(g['first'],g['stop']))!=[(0,1),(0,2),(1,2),(2,3)]:raise AssertionError('Gap/interval family mismatch')
    Z=np.array([[1.,-2.,.5],[-3.,1.,2.],[1e5,-1e5,0.]])
    ref=np.array([max(abs(math.fsum(row[a:b]))/math.sqrt(math.fsum(S[a:b,a:b].ravel()))
                     for a,b in zip(g['first'],g['stop'])) for row in Z])
    observed=maxima(Z,g)
    allowance=32*Z.shape[1]*np.finfo(float).eps*np.maximum(1,np.max(abs(Z),axis=1))/g['sum_sd_lower'].min()
    if np.any(observed<ref) or np.any(observed-ref>allowance):raise ArithmeticError('Independent maximum disagreement')
    # Scaling Z and its covariance must leave the standardized statistic fixed.
    scaled=maxima(7*Z,geometry(49*S,reset))
    if np.any(abs(scaled-observed)>allowance):raise ArithmeticError('Unit invariance failed')
    # Independent finite binomial enumeration and beta integration identity.
    for n,tail,failure in [(59,.05,.05),(199,.05,.01),(20000,.0125,.0125/3)]:
        k=tolerance_rank(n,tail,failure)
        beta=float(stats.beta.cdf(1-tail,k,n+1-k));binomial=float(stats.binom.sf(k-1,n,1-tail))
        if abs(beta-binomial)>8*n*np.finfo(float).eps:raise ArithmeticError('Order-statistic reference disagreement')
        if beta>failure or stats.binom.sf(k-2,n,1-tail)<=failure:raise AssertionError('Tolerance rank not minimal')
    # A single interval has M=|N(0,1)|; its exact quantile transformation provides
    # a reference without Monte Carlo fitting or another maximum implementation.
    u=np.linspace(.01,.99,99);z=stats.norm.ppf((1+u)/2)
    single=maxima(z[:,None],geometry(np.ones((1,1)),np.array([True])))
    if np.max(abs(single-z))>16*np.finfo(float).eps*np.max(z):raise ArithmeticError('One-dimensional reference failed')
    return dict(exhaustive_maximum=True,gap_family=True,unit_invariance=True,
                binomial_beta_identity=True,single_interval_normal_reference=True)
