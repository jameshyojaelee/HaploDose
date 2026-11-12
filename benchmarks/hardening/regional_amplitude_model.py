"""Research semi-Markov difference model with a shared signed regional effect.

The discrete-state path prior is unchanged. A nonzero run has one truncated
normal RNA effect, integrated jointly across its bins, rather than independent
per-bin draws or a finer effect grid. This is not a DNA-identification model.
"""
import numpy as np
from scipy.special import log_ndtr,logsumexp


def segment_evidence(track,variance,step,scale):
    """Integrated likelihood ratio to the zero-effect model for every run."""
    d,v=np.asarray(track,float),np.asarray(variance,float)
    if d.ndim!=1 or v.shape!=d.shape or not np.isfinite(d).all() or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('Finite aligned log-rate contrasts and positive variances required')
    if not np.isfinite(step) or step<=0 or not np.isfinite(scale) or scale<0:raise ValueError('Invalid prior scale')
    n=len(d);out=np.full((n,n,3),-np.inf)
    precision=np.r_[0.,np.cumsum(1/v)];weighted=np.r_[0.,np.cumsum(d/v)]
    for lo in range(n):
        A=precision[lo+1:]-precision[lo];B=weighted[lo+1:]-weighted[lo]
        out[lo,lo:,1]=0.
        for state,sign in ((0,-1),(2,1)):
            mu=sign*step
            if scale==0:
                value=mu*B-.5*mu**2*A
            else:
                H=A+1/scale**2;M=B+mu/scale**2
                value=(-.5*(mu/scale)**2+.5*M*M/H-.5*np.log(scale**2*H)
                       +log_ndtr(sign*M/np.sqrt(H))-log_ndtr(sign*mu/scale))
            out[lo,lo:,state]=value
    return out


def block(track,variance,transition,initial,step,scale):
    """Exact sum and max over state paths after integrating each run effect."""
    n=len(track);T=np.asarray(transition,float);I=np.asarray(initial,float)
    if n<1 or T.shape!=(3,3) or I.shape!=(3,):raise ValueError('Nonempty three-state block required')
    evidence=segment_evidence(track,variance,step,scale)
    scores=evidence.copy()
    for lo in range(n):scores[lo,lo:]+=np.arange(n-lo)[:,None]*np.diag(T)
    alpha=np.full((n,3),-np.inf);entry=np.full((n,3),-np.inf);entry[0]=I
    best=np.full((n,3),-np.inf);best_entry=entry.copy();previous=np.full((n,3),-1,int)
    first=np.full((n,3),-1,int);last_state=np.full((n,3),-1,int)
    for hi in range(n):
        for s in range(3):
            if hi:
                other=np.flatnonzero(np.arange(3)!=s)
                entry[hi,s]=logsumexp(alpha[hi-1,other]+T[other,s])
                values=best[hi-1,other]+T[other,s];which=int(np.argmax(values))
                best_entry[hi,s]=values[which];previous[hi,s]=other[which]
            alpha[hi,s]=logsumexp(entry[:hi+1,s]+scores[:hi+1,hi,s])
            values=best_entry[:hi+1,s]+scores[:hi+1,hi,s];lo=int(np.argmax(values))
            best[hi,s]=values[lo];first[hi,s]=lo;last_state[hi,s]=previous[lo,s]
    normalizer=float(logsumexp(alpha[-1]));back=np.full((n,3),-np.inf);future=np.zeros((n+1,3))
    for lo in range(n-1,-1,-1):
        for s in range(3):
            if lo<n-1:
                other=np.flatnonzero(np.arange(3)!=s)
                future[lo+1,s]=logsumexp(T[s,other]+back[lo+1,other])
            back[lo,s]=logsumexp(scores[lo,lo:,s]+future[lo+1:,s])
    reverse=float(logsumexp(I+back[0]))
    # Both recursions must agree within a magnitude/operation accumulation bound.
    u=np.finfo(float).eps/2;gamma=(64*n)*u/(1-64*n*u)
    guard=gamma*(1+abs(normalizer)+float(np.max(np.abs(scores[np.isfinite(scores)]))))
    assert abs(normalizer-reverse)<=guard
    mass=np.zeros_like(scores)
    for lo in range(n):mass[lo,lo:]=np.exp(entry[lo]+scores[lo,lo:]+future[lo+1:]-normalizer)
    cumulative=np.cumsum(mass,axis=0)
    post=np.stack([cumulative[j,j:].sum(axis=0) for j in range(n)])
    assert np.isfinite(post).all() and np.all(post>=0) and np.max(abs(post.sum(1)-1))<=guard
    path=np.empty(n,int);hi=n-1;s=int(np.argmax(best[-1]))
    while hi>=0:
        lo=int(first[hi,s]);path[lo:hi+1]=s;s=int(last_state[hi,s]);hi=lo-1
    return dict(post=post,path=path,log_evidence_ratio=normalizer,forward_backward_guard=guard,
                maximum_row_sum_error=float(np.max(abs(post.sum(1)-1))))


def infer(track,variance,reset,transition,initial,step,scale):
    track=np.asarray(track,float);reset=np.asarray(reset,bool)
    if reset.shape!=track.shape or not len(track) or not reset[0]:raise ValueError('Aligned genomic resets required')
    post=np.empty((len(track),3));path=np.empty(len(track),int);records=[]
    for idx in np.split(np.arange(len(track)),np.flatnonzero(reset)[1:]):
        result=block(track[idx],np.asarray(variance)[idx],transition,initial,step,scale)
        post[idx]=result['post'];path[idx]=result['path']
        records.append(dict(first=int(idx[0]),last=int(idx[-1]),**{k:v for k,v in result.items() if k not in ('post','path')}))
    return dict(post=post,path=path,blocks=records)
