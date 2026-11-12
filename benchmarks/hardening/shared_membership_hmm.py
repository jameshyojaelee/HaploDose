"""Research NB model with a shared cell membership and a dosage-state chain.

Exact state enumeration is a reference for short intervals. Mean-field inference
is an approximation, not a calibrated confidence procedure. Neither is enabled
in the package's inference path.
"""
import numpy as np
from scipy.special import expit, logit, logsumexp, rel_entr


def transition(prior, change):
    prior=np.asarray(prior,float)
    if np.any(prior<=0) or not np.isclose(prior.sum(),1.) or not 0<change<1:
        raise ValueError('Invalid state prior/transition probability')
    prior=prior/prior.sum()
    T=change*prior[None,:]/(1-prior[:,None])
    np.fill_diagonal(T,1-change)
    return T


def chain(emission, logT, initial, terminal):
    """Double-precision forward/backward and MAP path for one interval."""
    e=np.asarray(emission,float);n,k=e.shape
    forward=np.empty_like(e);forward[0]=initial+e[0]
    best=forward[0].copy();back=np.zeros((n,k),int)
    for j in range(1,n):
        forward[j]=e[j]+logsumexp(forward[j-1][:,None]+logT,axis=0)
        scores=best[:,None]+logT
        back[j]=scores.argmax(axis=0);best=e[j]+scores.max(axis=0)
    evidence=float(logsumexp(forward[-1]+terminal))
    backward=np.empty_like(e);backward[-1]=terminal
    for j in range(n-2,-1,-1):
        backward[j]=logsumexp(logT+e[j+1][None,:]+backward[j+1][None,:],axis=1)
    posterior=np.exp(forward+backward-evidence)
    path=np.empty(n,int);path[-1]=np.argmax(best+terminal)
    for j in range(n-2,-1,-1):path[j]=back[j+1,path[j+1]]
    return posterior,evidence,path


def build(counts, exposure, size, probability, effects, prior, change,
          *, left_neutral, right_neutral, mixture_proportion=None):
    O,E,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(exposure,float),np.asarray(size,float))
    q=np.asarray(probability,float);effects=np.asarray(effects,float)
    if O.ndim!=2 or q.shape!=(len(O),) or effects.shape!=(3,):raise ValueError('Invalid shapes')
    if not all(np.isfinite(x).all() for x in (O,E,r,q,effects)):raise ValueError('Nonfinite inputs')
    if np.any(O<0) or np.any(E<=0) or np.any(r<=0) or np.any((q<0)|(q>1)):
        raise ValueError('Invalid count/exposure/size/probability')
    # A cell subset retains the original component rates. Recentring on its
    # selected membership composition would silently change the fitted model.
    pi=(np.clip((E.T@q)/E.sum(axis=0),0.,1.) if mixture_proportion is None
        else np.asarray(mixture_proportion,float))
    if pi.shape!=(O.shape[1],) or not np.isfinite(pi).all() or np.any((pi<0)|(pi>1)):
        raise ValueError('Invalid component-normalization proportions')
    centre=np.log1p(pi[:,None]*np.expm1(effects)[None,:])
    shifts=np.stack([-centre,effects[None,:]-centre])
    ratio=expit(np.log(E)-np.log(r))
    ll=[]
    for shift in shifts:
        ll.append(O[:,:,None]*shift[None]-(O+r)[:,:,None]*np.log1p(ratio[:,:,None]*np.expm1(shift)[None]))
    logT=np.log(transition(prior,change));initial=logT[1].copy() if left_neutral else np.log(np.asarray(prior,float))
    terminal=logT[:,1].copy() if right_neutral else np.zeros(3)
    normalizer=chain(np.zeros((O.shape[1],3)),logT,initial,terminal)[1]
    return dict(ll=np.asarray(ll),q=q,logT=logT,initial=initial,terminal=terminal,
                prior_normalizer=normalizer,shifts=shifts,mixture_proportion=pi)


def enumerate_states(model, max_states=1000000, batch_size=128):
    """Integrate every Z independently conditional on each shared state path."""
    ll,q=model['ll'],model['q'];_,cells,length,states=ll.shape
    number=states**length
    if number>max_states:raise ValueError('Exact state enumeration exceeds declared computational cap')
    paths=((np.arange(number)[:,None]//states**np.arange(length-1,-1,-1))%states).astype(np.int8)
    with np.errstate(divide='ignore'):lq0,lq1=np.log1p(-q),np.log(q)
    def evaluate(paths):
        a=np.broadcast_to(lq0,(len(paths),cells)).copy();b=np.broadcast_to(lq1,a.shape).copy()
        for j in range(length):
            a+=ll[0,:,j,:][:,paths[:,j]].T;b+=ll[1,:,j,:][:,paths[:,j]].T
        joint=np.logaddexp(a,b)
        value=model['initial'][paths[:,0]]+model['terminal'][paths[:,-1]]
        for j in range(1,length):value+=model['logT'][paths[:,j-1],paths[:,j]]
        return value+joint.sum(axis=1)-model['prior_normalizer'],expit(b-a)
    log_joint=np.empty(number)
    for start in range(0,number,batch_size):
        stop=min(start+batch_size,number);log_joint[start:stop]=evaluate(paths[start:stop])[0]
    evidence=float(logsumexp(log_joint));weights=np.exp(log_joint-evidence)
    post=np.stack([np.bincount(paths[:,j],weights=weights,minlength=states) for j in range(length)])
    membership=np.zeros(cells)
    for start in range(0,number,batch_size):
        stop=min(start+batch_size,number)
        membership+=weights[start:stop]@evaluate(paths[start:stop])[1]
    return dict(log_evidence=evidence,posterior=post,membership=membership,
                path=paths[np.argmax(log_joint)].copy(),states_enumerated=number)


def variational(model, tolerance=None, max_iterations=2000):
    """Coordinate ascent in Q(Z)Q(S); retain all starts and convergence status."""
    if tolerance is None:tolerance=np.sqrt(np.finfo(float).eps)
    ll,q=model['ll'],model['q'];logodds=logit(q)
    starts=[q.copy()]+[expit(logodds+(ll[1,:,:,s]-ll[0,:,:,s]).sum(axis=1)) for s in (0,2)]
    allowance=64*np.finfo(float).eps*max(1.,float(np.abs(ll).sum()))
    fits=[]
    for start in starts:
        w=start;history=[];converged=False
        for iteration in range(max_iterations):
            emit=((1-w)[:,None,None]*ll[0]+w[:,None,None]*ll[1]).sum(axis=0)
            post,logZ,path=chain(emit,model['logT'],model['initial'],model['terminal'])
            elbo=float(logZ-model['prior_normalizer']-np.sum(rel_entr(w,q)+rel_entr(1-w,1-q)))
            if history and elbo<history[-1]-allowance:raise ArithmeticError('Coordinate ascent decreased ELBO')
            new=expit(logodds+np.einsum('cjs,js->c',ll[1]-ll[0],post))
            residual=float(np.max(np.abs(new-w)))
            history.append(elbo)
            if residual<=tolerance and (len(history)==1 or abs(history[-1]-history[-2])<=allowance):
                converged=True;break
            if iteration+1<max_iterations:w=new
        fits.append(dict(elbo=elbo,posterior=post,membership=w,path=path,
                         converged=converged,iterations=iteration+1,residual=residual,history=history))
    winner=max(range(len(fits)),key=lambda i:fits[i]['elbo'])
    return dict(**fits[winner],starts=fits,selected_start=winner,numerical_elbo_allowance=allowance)
