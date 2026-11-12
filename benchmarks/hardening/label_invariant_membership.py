"""Exact posterior on unordered two-component trajectories, for research only.

The higher-rate component at the first nonneutral bin is named component one.
One label swap applies to the ENTIRE trajectory and every cell, not each bin.
This is not applicable to confirmation with anchored cell-specific memberships.
"""
import numpy as np
from scipy.special import expit, logsumexp


def canonical(paths):
    paths=np.asarray(paths,int)
    nonneutral=paths!=1;identified=nonneutral.any(axis=1)
    first=np.argmax(nonneutral,axis=1)
    flip=identified & (paths[np.arange(len(paths)),first]==0)
    return np.where(flip[:,None],2-paths,paths),flip,identified,first


def enumerate_unordered(model,max_states=1000000,batch_size=128):
    ll,q=model['ll'],np.asarray(model['q']);_,cells,length,states=ll.shape
    if states!=3 or not np.array_equal(q,np.full_like(q,q[0])) or not 0<q[0]<1:
        raise ValueError('Unordered inference requires a common unanchored, interior membership prior')
    if (not np.array_equal(model['logT'],model['logT'][::-1,::-1]) or
        not np.array_equal(model['initial'],model['initial'][::-1]) or
        not np.array_equal(model['terminal'],model['terminal'][::-1])):
        raise ValueError('State prior must be symmetric under global component exchange')
    number=states**length
    if number>max_states:raise ValueError('Exact enumeration exceeds declared cap')
    place=states**np.arange(length-1,-1,-1)
    paths=((np.arange(number)[:,None]//place)%states).astype(np.int8)
    named,flip,identified,first=canonical(paths)
    code=named@place
    def evaluate(part):
        a=np.full((len(part),cells),np.log1p(-q[0]));b=np.full_like(a,np.log(q[0]))
        for j in range(length):
            a+=ll[0,:,j,:][:,part[:,j]].T;b+=ll[1,:,j,:][:,part[:,j]].T
        prior=model['initial'][part[:,0]]+model['terminal'][part[:,-1]]
        for j in range(1,length):prior+=model['logT'][part[:,j-1],part[:,j]]
        return prior+np.logaddexp(a,b).sum(axis=1)-model['prior_normalizer'],expit(b-a)
    log_joint=np.empty(number)
    for start in range(0,number,batch_size):
        stop=min(start+batch_size,number);log_joint[start:stop]=evaluate(paths[start:stop])[0]
    evidence=float(logsumexp(log_joint));mass=np.exp(log_joint-evidence)
    # Sum full trajectory orbits BEFORE selecting a MAP representative. Taking
    # abs(state) independently per bin would destroy relative sign information.
    collapsed=np.full(number,-np.inf);np.logaddexp.at(collapsed,code,log_joint)
    present=np.flatnonzero(np.isfinite(collapsed));probability=np.exp(collapsed[present]-evidence)
    representatives=paths[present]
    posterior=np.array([np.bincount(representatives[:,j],weights=probability,minlength=3) for j in range(length)])
    membership=np.zeros(cells)
    for start in range(0,number,batch_size):
        stop=min(start+batch_size,number);conditional=evaluate(paths[start:stop])[1]
        conditional=np.where(flip[start:stop,None],1-conditional,conditional)
        # A completely neutral trajectory has no CNA-defined component labels.
        conditional=np.where(identified[start:stop,None],conditional,.5)
        membership+=mass[start:stop]@conditional
    anchors=np.bincount(first[identified],weights=mass[identified],minlength=length)
    return dict(log_evidence=evidence,posterior=posterior,membership=membership,
        path=paths[np.argmax(collapsed)].copy(),neutral_trajectory_probability=float(mass[~identified].sum()),
        anchor_probability=anchors,states_enumerated=number,unordered_trajectories=len(present),
        label_convention='component1 higher at first nonneutral bin; undefined on all-neutral trajectory')
