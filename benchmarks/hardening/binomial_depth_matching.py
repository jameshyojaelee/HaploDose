"""Research physical count-depth matching with exact conditional LR averaging."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import itertools
import math
from pathlib import Path

import numpy as np
from scipy.stats import binom, nbinom
from haplodose.io import write_json


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('binomial_depth_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def model(counts,expected,size,q,gene_bins,effects,prior,change,weights,
          *,left_neutral,right_neutral,mixture_proportion):
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    w=np.asarray(weights,float)
    if (y.ndim!=2 or w.shape!=(len(y),) or not np.isfinite(w).all()
            or np.any((w<0)|(w>1)) or np.any(y!=np.floor(y))):
        raise ValueError('Integer counts and aligned finite matching weights in[0,1] required')
    shared,engine=load('gene_shared_confirmation'),load('shared_membership_hmm')
    opts=dict(left_neutral=left_neutral,right_neutral=right_neutral,mixture_proportion=mixture_proportion)
    if np.all(w==1):return shared.gene_model(y,e,r,q,gene_bins,effects,prior,change,**opts)
    local=np.asarray(gene_bins,int);length=len(mixture_proportion)
    if local.shape!=(y.shape[1],) or set(local)!=set(range(length)):
        raise ValueError('Every event bin requires identified gene composition')
    result=engine.build(y,e,r,q,effects,prior,change,left_neutral=left_neutral,
        right_neutral=right_neutral,mixture_proportion=np.asarray(mixture_proportion)[local])
    old=result['ll'];ll=[];mean=e*w[:,None]
    for component,shift in enumerate(result['shifts']):
        v=np.log1p((mean/(mean+r))[:,:,None]*np.expm1(shift)[None])
        averaged=y[:,:,None]*np.log1p(w[:,None,None]*np.expm1(shift[None]-v))-r[:,:,None]*v
        averaged[w==1]=old[component,w==1];averaged[w==0]=0.
        ll.append(averaged)
    ll=np.asarray(ll)
    result['ll']=np.stack([ll[:,:,local==j,:].sum(axis=2) for j in range(length)],axis=2)
    result['prior_normalizer']=engine.chain(np.zeros((length,3)),result['logT'],result['initial'],result['terminal'])[1]
    return result


def checks():
    # Exact marginal identity for physical downsampling of NB(2,mean2).
    # Remaining NB tail encloses every omitted binomial probability.
    tail=Fraction(243,2**242);marginal=[]
    for k in range(3):
        partial=sum(Fraction((n+1)*math.comb(n,k),2**(2*n+2)) for n in range(k,241))
        exact=Fraction(4*(k+1),3**(k+2))
        assert 0<=exact-partial<=tail
        marginal.append(dict(count=k,probability=str(exact),omitted_upper=float(tail)))
    covariance=Fraction(1,2)*Fraction(1,2)*Fraction(2)**2/Fraction(2)
    assert covariance==Fraction(1,2)
    Y=np.array([[2,1],[1,0]]);E=np.array([[2.,3.],[4.,2.]])
    sizes=np.full(2,4.);weights=np.array([.25,.5]);q=np.array([.35,.7])
    bins=np.array([0,1]);effects=np.array([-math.log(1.5),0.,math.log(1.5)])
    prior=np.array([.15,.7,.15]);change=.02;pi=np.array([.4,.6])
    opts=dict(left_neutral=False,right_neutral=False,mixture_proportion=pi)
    engine=load('shared_membership_hmm');resolved=model(Y,E,sizes,q,bins,effects,prior,change,weights,**opts)
    fitted=engine.enumerate_states(resolved)
    paths=list(itertools.product(range(3),repeat=2));labels=list(itertools.product(range(2),repeat=2))
    retained=list(itertools.product(*(range(int(n)+1) for n in Y.ravel())))
    joint=np.zeros(len(paths));member=np.zeros(2)
    for j,path in enumerate(paths):
        state=prior[path[0]]*((1-change) if path[0]==path[1] else change*prior[path[1]]/(1-prior[path[0]]))
        for allocation in labels:
            pz=math.prod(q[i] if allocation[i] else 1-q[i] for i in range(2));value=0.
            for values in retained:
                K=np.array(values).reshape(Y.shape)
                kernel=math.prod(float(binom.pmf(K[i,g],Y[i,g],weights[i])) for i in range(2) for g in range(2))
                ratios=[]
                for i in range(2):
                    for g in range(2):
                        delta=effects[path[g]];multiplier=math.exp(allocation[i]*delta)/(1-pi[g]+pi[g]*math.exp(delta))
                        shape=sizes[g];mean=weights[i]*E[i,g]
                        ratios.append(math.exp(float(nbinom.logpmf(K[i,g],shape,shape/(shape+mean*multiplier))
                            -nbinom.logpmf(K[i,g],shape,shape/(shape+mean)))))
                value+=kernel*math.prod(ratios)
            mass=state*pz*value;joint[j]+=mass;member+=mass*np.array(allocation)
    evidence=math.fsum(map(float,joint));post=np.zeros((2,3))
    for path,mass in zip(paths,joint):
        for g,s in enumerate(path):post[g,s]+=mass/evidence
    operations=32*len(paths)*len(labels)*len(retained)*Y.size+256
    allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(math.log(evidence)))
    errors=[abs(fitted['log_evidence']-math.log(evidence)),float(np.max(abs(post-fitted['posterior']))),
        float(np.max(abs(member/evidence-fitted['membership'])))]
    assert max(errors)<=allowance
    # Independently evaluate closed forms at extended precision, including rare
    # retained reads and a high count. Actual extended precision is recorded.
    for count,mean,size,weight in [(0,3.,4.,.2),(4,3.,4.,.2),(1000,50.,2.,.999),(61,7.,1.2,.01)]:
        yy,ee,rr,ww=map(np.longdouble,(count,mean,size,weight))
        tiny=model([[count]],[[mean]],[size],[.35],[0],effects,prior,change,[weight],
            left_neutral=False,right_neutral=False,mixture_proportion=[.4])
        for component in range(2):
            for state,delta in enumerate(effects):
                ss=np.longdouble(component)*np.longdouble(delta)-np.log1p(np.longdouble(.4)*np.expm1(np.longdouble(delta)))
                vv=np.log1p(ww*ee/(rr+ww*ee)*np.expm1(ss))
                reference=float(yy*np.log1p(ww*np.expm1(ss-vv))-rr*vv)
                actual=float(tiny['ll'][component,0,0,state]);budget=256*np.finfo(float).eps*max(1.,abs(reference),count)
                assert abs(actual-reference)<=budget;errors.append(abs(actual-reference))
    unchanged=model(Y,E,sizes,q,bins,effects,prior,change,np.ones(2),**opts)
    old=load('gene_shared_confirmation').gene_model(Y,E,sizes,q,bins,effects,prior,change,**opts)
    for key in ('ll','q','logT','initial','terminal'):np.testing.assert_array_equal(unchanged[key],old[key])
    np.testing.assert_array_equal(model(Y,E,sizes,q,bins,effects,prior,change,np.zeros(2),**opts)['ll'],np.zeros_like(old['ll']))
    swapped=engine.enumerate_states(model(Y[::-1,::-1],E[::-1,::-1],sizes[::-1],q[::-1],bins[::-1],
        effects,prior,change,weights[::-1],**opts))
    assert np.max(abs(swapped['posterior']-fitted['posterior']))<=allowance
    assert np.max(abs(swapped['membership'][::-1]-fitted['membership']))<=allowance
    return dict(exact_marginal_enclosures=marginal,retained_discarded_covariance=str(covariance),
        complete_count_label_path_enumeration=True,enumeration_allowance=allowance,
        maximum_reference_discrepancy=max(errors),extended_precision_bits=np.finfo(np.longdouble).nmant,
        endpoints_exact=True,cell_gene_permutation=True,
        scope='Fixed NB observation law; no adaptive or biological false-split guarantee')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:load('averaged_count_matching').study(args.out,args.task,matching_model=model,numerical_checks=checks,
        extra_inputs=dict(binomial_kernel=Path(__file__)))


if __name__=='__main__':main()
