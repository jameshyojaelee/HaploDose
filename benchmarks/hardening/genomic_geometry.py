"""Research transition geometry with matched block process time, not a default."""
import argparse
import itertools
import json
from pathlib import Path

import jax.numpy as jnp
import numpy as np
from scipy.linalg import expm, logm
from scipy.special import logsumexp

from haplodose.kernels.hmm import forward_backward, sticky_transition
from haplodose.profile import STATE_PRIOR
from haplodose.provenance import start_run, finish_run


def transition_model(p_change=.002):
    prior=STATE_PRIOR.astype(float);prior/=prior.sum()
    transition=p_change*prior[None]/(1-prior[:,None])
    np.fill_diagonal(transition,1-p_change)
    generator=logm(transition)
    if np.iscomplexobj(generator):raise ValueError('Transition logarithm is not real')
    off=generator.copy();np.fill_diagonal(off,0.)
    if np.any(off<0):raise ValueError('Transition is not embeddable in a Markov generator')
    return transition,generator,prior


def log_transitions(ratios,generator):
    matrices=np.stack([np.eye(4) if ratio==0 else expm(generator*ratio) for ratio in ratios])
    if np.any(matrices<0):raise ArithmeticError('Negative transition probability; do not clip it')
    with np.errstate(divide='ignore'):logs=np.log(matrices).astype(np.float32)
    return logs,matrices


def geometry(bins,genes,seed=20261109):
    from haplodose.bins import chromosome_starts
    reset=chromosome_starts(bins)
    anchors=genes.groupby('bin').start.median().reindex(range(len(bins))).to_numpy(float)
    if not np.isfinite(anchors).all():raise ValueError('Missing measured-locus anchor')
    ratios=np.zeros(len(bins));shuffled=np.zeros(len(bins));distance=np.zeros(len(bins))
    permutation=np.arange(len(bins));rng=np.random.default_rng(seed);blocks=[]
    for index in np.split(np.arange(len(bins)),np.flatnonzero(reset)[1:]):
        edges=index[1:];d=np.diff(anchors[index])
        if np.any(d<0):raise ValueError('Nonmonotone genomic anchors within block')
        if len(edges) and d.sum()<=0:raise ValueError('No physical extent identifies this block scale')
        if len(edges):
            ratios[edges]=d/d.mean();distance[edges]=d
            perm=rng.permutation(edges);shuffled[edges]=ratios[perm];permutation[edges]=perm
        blocks.append(dict(first_bin=int(index[0]),last_bin=int(index[-1]),edges=len(edges),
            anchor_extent_bp=float(d.sum()),scale_bp=float(d.mean()) if len(d) else None,
            process_time=float(ratios[edges].sum())))
    return dict(reset=reset,anchors=anchors,distance_bp=distance,ratios=ratios,
                shuffled_ratios=shuffled,permutation=permutation,blocks=blocks)


def reference_posterior(emit,log_transition,reset,prior=STATE_PRIOR):
    """Independent float64 recursion, resetting both messages for every block."""
    log_transition=np.asarray(log_transition,float);emit=np.asarray(emit,float)
    init=np.log(np.asarray(prior,np.float32)).astype(float);out=np.empty_like(emit)
    for index in np.split(np.arange(len(emit)),np.flatnonzero(reset)[1:]):
        e=emit[index];a=np.empty_like(e);b=np.zeros_like(e)
        a[0]=init+e[0];a[0]-=logsumexp(a[0])
        for k in range(1,len(index)):
            t=log_transition if log_transition.ndim==2 else log_transition[index[k]]
            a[k]=e[k]+logsumexp(a[k-1,:,None]+t,axis=0);a[k]-=logsumexp(a[k])
        for k in range(len(index)-2,-1,-1):
            t=log_transition if log_transition.ndim==2 else log_transition[index[k+1]]
            b[k]=logsumexp(t+e[k+1,None,:]+b[k+1,None,:],axis=1);b[k]-=max(b[k])
        total=a+b;out[index]=np.exp(total-logsumexp(total,axis=1,keepdims=True))
    return out


def checks():
    transition,generator,prior=transition_model()
    ratios=np.array([0.,.2,1.,2.,5.,8.3])
    logs,matrices=log_transitions(ratios,generator)
    # Independent spectral construction follows detailed balance.
    stationary=prior*(1-prior);stationary/=stationary.sum();root=np.sqrt(stationary)
    symmetric=root[:,None]*transition/root[None]
    eigenvalue,vectors=np.linalg.eigh(symmetric)
    if np.any(eigenvalue<=0):raise ValueError('Positive eigenvalues required')
    differences=[]
    for ratio,matrix in zip(ratios,matrices):
        spectral=np.eye(4) if ratio==0 else ((vectors*eigenvalue**ratio)@vectors.T)*root[None]/root[:,None]
        differences.append(float(np.max(np.abs(matrix-spectral))))
    # Error allowance scales with dense4x4 matrix products, not a scientific cutoff.
    operations=4**3*len(ratios);u=np.finfo(float).eps
    budget=operations*u/(1-operations*u)
    assert max(differences)<budget
    assert np.max(np.abs(matrices.sum(axis=2)-1))<budget
    for steps in (0,1,2,5):
        assert np.max(np.abs(expm(generator*steps)-np.linalg.matrix_power(transition,steps)))<budget
    assert np.max(np.abs(expm(generator*.2)@expm(generator*8.3)-expm(generator*8.5)))<budget
    # Exhaustive state paths include a zero-distance identity and a chain reset.
    rng=np.random.default_rng(20261109);emit=rng.normal(size=(4,4)).astype(np.float32)
    reset=np.array([True,False,True,False]);lt,_=log_transitions([0.,0.,0.,2.],generator)
    enumerated=np.zeros((4,4));weights=[];paths=[]
    for path in itertools.product(range(4),repeat=4):
        value=sum(float(emit[t,s]) for t,s in enumerate(path))
        for t,s in enumerate(path):value+=float(np.log(STATE_PRIOR[s])) if reset[t] else float(lt[t,path[t-1],s])
        weights.append(value);paths.append(path)
    weights=np.exp(np.array(weights)-logsumexp(weights))
    for path,weight in zip(paths,weights):
        for t,s in enumerate(path):enumerated[t,s]+=weight
    reference=reference_posterior(emit,lt,reset)
    got=np.asarray(forward_backward(jnp.asarray(emit),jnp.asarray(lt),jnp.log(jnp.asarray(STATE_PRIOR)),jnp.asarray(reset))[0])
    # Retain the existing1e-4 HMM posterior criterion, with a binary64 enumeration check too.
    assert np.max(np.abs(reference-enumerated))<budget
    assert np.max(np.abs(got-reference))<1e-4
    return dict(spectral_transition_max_difference=max(differences),binary64_matrix_budget=budget,
        enumeration_paths=len(paths),reference_enumeration_max_difference=float(np.max(np.abs(reference-enumerated))),
        JAX_reference_max_difference=float(np.max(np.abs(got-reference))),existing_HMM_criterion=1e-4)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md'),{},
        dict(scope='Independent numerical preflight; no biological claim'),mode='genomic_geometry_preflight')
    try:
        (args.out/'checks.json').write_text(json.dumps(checks(),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
