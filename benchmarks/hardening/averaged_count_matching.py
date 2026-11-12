"""Conditional averaging of the existing beta-binomial count-matching kernel."""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.special import logsumexp
from scipy.stats import betabinom, nbinom
from haplodose import profile
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('averaged_matching_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def log_moment(counts, size, weights, exponent):
    """log E[exp(t*K)|Y] over the COMPLETE beta-binomial support0..Y."""
    y,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(size,float))
    w=np.broadcast_to(np.asarray(weights,float),y.shape)
    t=np.asarray(exponent,float)
    if (y.ndim!=2 or t.shape[:2]!=y.shape or t.ndim!=3 or not np.isfinite(y).all()
            or np.any(y<0) or np.any(y!=np.floor(y)) or not np.isfinite(r).all()
            or np.any(r<=0) or not np.isfinite(w).all() or np.any((w<0)|(w>1))
            or not np.isfinite(t).all()):
        raise ValueError('Integer counts, positive sizes, weights in[0,1] and aligned finite exponents required')
    out=np.zeros(t.shape)
    one=w==1
    out[one]=y[one][:,None]*t[one]
    interior=(w>0)&(w<1)
    for n in np.unique(y[interior]).astype(np.int64):
        keep=interior&(y==n)
        alpha,beta=(r*w)[keep],(r*(1-w))[keep]
        tt=t[keep]
        # Relative probability at K=0 is1. Normalization cancels its unknown
        # constant; the recurrence never needs subtractive gamma functions.
        logp=np.zeros(len(tt));normalizer=np.zeros(len(tt));moment=np.zeros(tt.shape)
        for k in range(1,int(n)+1):
            logp+=math.log(int(n)-k+1)-math.log(k)+np.log(alpha+(k-1))-np.log(beta+(int(n)-k))
            normalizer=np.logaddexp(normalizer,logp)
            moment=np.logaddexp(moment,logp[:,None]+k*tt)
        out[keep]=moment-normalizer[:,None]
    return out


def model(counts, expected, size, q, gene_bins, effects, prior, change, weights,
          *, left_neutral, right_neutral, mixture_proportion):
    shared,engine=load('gene_shared_confirmation'),load('shared_membership_hmm')
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    w=np.asarray(weights,float)
    if w.shape!=(len(y),):raise ValueError('One executed matching weight per cell required')
    options=dict(left_neutral=left_neutral,right_neutral=right_neutral,mixture_proportion=mixture_proportion)
    if np.all(w==1):
        # Exact original calculation, including its floating-point operation order.
        return shared.gene_model(y,e,r,q,gene_bins,effects,prior,change,**options)
    local=np.asarray(gene_bins,int);length=len(mixture_proportion)
    if local.shape!=(y.shape[1],) or set(local)!=set(range(length)):
        raise ValueError('Complete gene-to-event-bin correspondence required')
    result=engine.build(y,e,r,q,effects,prior,change,left_neutral=left_neutral,
        right_neutral=right_neutral,mixture_proportion=np.asarray(mixture_proportion)[local])
    old=result['ll']
    ll=[]
    for k,shift in enumerate(result['shifts']):
        v=np.log1p((e/(e+r))[:,:,None]*np.expm1(shift)[None])
        averaged=log_moment(y,r,w[:,None],shift[None]-v)-(w[:,None]*r)[:,:,None]*v
        averaged[w==1]=old[k,w==1]
        averaged[w==0]=0.
        ll.append(averaged)
    ll=np.asarray(ll)
    result['ll']=np.stack([ll[:,:,local==j,:].sum(axis=2) for j in range(length)],axis=2)
    result['prior_normalizer']=engine.chain(np.zeros((length,3)),result['logT'],result['initial'],result['terminal'])[1]
    return result


def decimal_moment(n, alpha, beta, t, precision):
    with localcontext() as context:
        context.prec=precision
        a,b,tt=map(Decimal.from_float,(float(alpha),float(beta),float(t)))
        p=Decimal(1);mass=p;moment=p
        for k in range(1,n+1):
            p*=Decimal(n-k+1)/Decimal(k)*(a+k-1)/(b+n-k)
            mass+=p;moment+=p*(k*tt).exp()
        return float((moment/mass).ln())


def checks():
    # Y=4, r=3, w=1/3; K~BetaBinomial(4,1,2).
    # Baseline matched mean1 and alternative2 have LR=(2/3)*(4/3)^K.
    exact=sum(Fraction(2*(5-k),30)*Fraction(2,3)*Fraction(4,3)**k for k in range(5))
    y=np.array([[4.]]);r=np.array([[3.]]);w=np.array([[1/3]])
    t=np.array([[[math.log(4/3)]]])
    observed=float(log_moment(y,r,w,t)[0,0,0]+math.log(2/3))
    reference=float(math.log(float(exact)))
    operations=64*5+128
    fraction_allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(reference))
    assert abs(observed-reference)<=fraction_allowance
    errors=[abs(observed-reference)];references=[]
    for n,rr,ww,tt in [(0,2.5,.3,-.4),(1,3.,.2,.3),(6,2.5,.3,-.4),(23,7.5,.8,.2),(61,1.2,.01,.6)]:
        actual=float(log_moment([[n]],[[rr]],[[ww]],[[[tt]]])[0,0,0])
        ref60=decimal_moment(n,rr*ww,rr*(1-ww),tt,60)
        ref90=decimal_moment(n,rr*ww,rr*(1-ww),tt,90)
        assert ref60==ref90
        operations=64*(n+1)+128
        allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(ref90),n*abs(tt))
        assert abs(actual-ref90)<=allowance
        errors.append(abs(actual-ref90));references.append(dict(count=n,actual=actual,reference=ref90,allowance=allowance))
    # A rigorous rational enclosure for the fixed-null e-value expectation.
    # Null Y=K+U: K~NB(1,mean1), U~NB(2,mean2). Under LR weighting,
    # K instead has mean2. A union bound encloses the omitted N>240 tail.
    partial=Fraction(0)
    for n in range(241):
        kernel=sum(Fraction(2*(n-k+1),(n+1)*(n+2))*Fraction(2,3)*Fraction(4,3)**k for k in range(n+1))
        partial+=Fraction(math.comb(n+2,2),2**(n+3))*kernel
    tail=Fraction(2,3)**121+Fraction(123,2**122)
    assert 0<=1-partial<=tail
    # Independent full SciPy enumeration of retained counts, cell labels and paths.
    Y=np.array([[2,1],[1,0]]);E=np.array([[2.,3.],[4.,2.]])
    sizes=np.full(2,4.);weight=np.array([.25,.5]);q=np.array([.35,.7])
    bins=np.array([0,1]);effects=np.array([-math.log(1.5),0.,math.log(1.5)])
    prior=np.array([.15,.7,.15]);change=.02;pi=np.array([.4,.6])
    opts=dict(left_neutral=False,right_neutral=False,mixture_proportion=pi)
    fitted=load('shared_membership_hmm').enumerate_states(model(Y,E,sizes,q,bins,effects,prior,change,weight,**opts))
    paths=list(itertools.product(range(3),repeat=2));allocations=list(itertools.product(range(2),repeat=2))
    thin=list(itertools.product(*(range(int(n)+1) for n in Y.ravel())))
    joint=np.zeros(len(paths));member=np.zeros(2)
    for j,path in enumerate(paths):
        pstate=prior[path[0]]*((1-change) if path[0]==path[1] else change*prior[path[1]]/(1-prior[path[0]]))
        for labels in allocations:
            pz=math.prod(q[i] if labels[i] else 1-q[i] for i in range(2))
            value=0.
            for retained in thin:
                K=np.array(retained).reshape(Y.shape)
                kernel=math.prod(float(betabinom.pmf(K[i,g],Y[i,g],weight[i]*sizes[g],(1-weight[i])*sizes[g]))
                    for i in range(2) for g in range(2))
                ratios=[]
                for i in range(2):
                    for g in range(2):
                        delta=effects[path[g]]
                        multiplier=math.exp(labels[i]*delta)/(1-pi[g]+pi[g]*math.exp(delta))
                        shape=weight[i]*sizes[g];mean=weight[i]*E[i,g]
                        ratios.append(math.exp(float(nbinom.logpmf(K[i,g],shape,shape/(shape+mean*multiplier))
                            -nbinom.logpmf(K[i,g],shape,shape/(shape+mean)))))
                value+=kernel*math.prod(ratios)
            mass=pstate*pz*value;joint[j]+=mass;member+=mass*np.array(labels)
    evidence=math.fsum(map(float,joint));post=np.zeros((2,3))
    for path,mass in zip(paths,joint):
        for g,s in enumerate(path):post[g,s]+=mass/evidence
    operations=32*len(paths)*len(allocations)*len(thin)*Y.size+256
    allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(math.log(evidence)))
    for error in [abs(fitted['log_evidence']-math.log(evidence)),float(np.max(abs(post-fitted['posterior']))),
                  float(np.max(abs(member/evidence-fitted['membership'])))]:
        assert error<=allowance;errors.append(error)
    old=load('gene_shared_confirmation').gene_model(Y,E,sizes,q,bins,effects,prior,change,**opts)
    unchanged=model(Y,E,sizes,q,bins,effects,prior,change,np.ones(2),**opts)
    for key in ('ll','q','logT','initial','terminal'):np.testing.assert_array_equal(old[key],unchanged[key])
    zero=model(Y,E,sizes,q,bins,effects,prior,change,np.zeros(2),**opts)
    np.testing.assert_array_equal(zero['ll'],np.zeros_like(zero['ll']))
    swapped=model(Y[::-1,::-1],E[::-1,::-1],sizes[::-1],q[::-1],bins[::-1],effects,prior,change,weight[::-1],**opts)
    other=load('shared_membership_hmm').enumerate_states(swapped)
    assert np.max(abs(other['posterior']-fitted['posterior']))<=allowance
    assert np.max(abs(other['membership'][::-1]-fitted['membership']))<=allowance
    rejected=0
    for yy,ww in [([[.5]],[[.5]]),([[1]],[[1.1]]),([[1]],[[-.1]])]:
        try:log_moment(yy,[[3.]],ww,[[[.2]]])
        except ValueError:rejected+=1
    assert rejected==3
    return dict(exact_fraction_ratio=str(exact),fraction_log_error=abs(observed-reference),
        fraction_allowance=fraction_allowance,invalid_inputs_rejected=rejected,
        decimal_convergence=references,maximum_reference_discrepancy=max(errors),
        complete_count_label_path_enumeration=True,enumeration_allowance=allowance,
        one_endpoint_exact=True,zero_endpoint_no_evidence=True,cell_gene_permutation=True,
        null_expectation_fraction_enclosed=True,null_expectation_tail_upper=float(tail),
        scope='Fixed correctly specified null only; no adaptive or biological error guarantee')


def study(out,task,*,matching_model=None,numerical_checks=None,extra_inputs=None,weight_provider=None):
    matching_model=model if matching_model is None else matching_model
    numerical_checks=checks if numerical_checks is None else numerical_checks
    root=task.parents[1];prior=root/'dev/gene_shared_confirmation_20261109/count_gene_v1'
    trace=json.loads((prior/'stage_trace.json').read_text())
    splits=json.loads((prior/'candidate/split_records.json').read_text())
    inputs=dict(script=Path(__file__),design=task/'DESIGN.md',trace=prior/'stage_trace.json',
        splits=prior/'candidate/split_records.json',scorer=Path(profile.__file__),
        engine=Path(__file__).with_name('shared_membership_hmm.py'),
        gene=Path(__file__).with_name('gene_shared_confirmation.py'))
    for i,row in enumerate(splits):inputs[f'fold{i}']=prior/f"candidate/folds/model_{row['model']}.npz"
    for i,entry in enumerate(trace):
        for j,call in enumerate(entry.get('calls',[])):
            if call.get('arrays'):inputs[f'comparison{i}_{j}']=Path(call['arrays'])
    inputs.update(extra_inputs or {})
    run=start_run(out,inputs,dict(posterior=.95,min_bins=5),
        dict(scope='All saved exact comparisons, conditional evidence only; not a recovery run'),mode='averaged_count_matching')
    try:
        numerical=numerical_checks();folds={};records=[];new_core=False;replay_max=0.
        for row in splits:
            path=prior/f"candidate/folds/model_{row['model']}.npz"
            with np.load(path) as z:
                data={k:z[k] for k in z.files}
            def csr(name):return sparse.csr_matrix((data[name+'_data'],data[name+'_indices'],data[name+'_indptr']),shape=tuple(data[name+'_shape']))
            source_counts,proposal_counts=csr('source'),csr('prop')
            folds.setdefault(row['case_seed'],[]).append(dict(rows=data['global_rows'],counts=source_counts-proposal_counts,
                source=source_counts,proposal=proposal_counts,gene_bins=data['gbin']))
        engine=load('shared_membership_hmm');diff=profile.ProfileParams(posterior=.95,min_bins=5)
        for i,entry in enumerate(trace):
            if entry.get('stage')!='gene_shared_confirmation':
                records.append(dict(index=i,stage=entry.get('stage'),source=entry.get('source'),
                    skipped='Different recorded pipeline stage; original record retained in input trace'))
                continue
            initial=next((c for c in entry.get('calls',[]) if c.get('context')=='initial' and c.get('arrays')),None)
            if initial is None:
                records.append(dict(index=i,case_seed=entry['case_seed'],source=entry['source'],
                    skipped='No saved identified initial gene model; original abstention/fallback retained'))
                continue
            with np.load(initial['arrays']) as z:a={k:z[k] for k in z.files}
            candidates=[f for f in folds[entry['case_seed']] if set(a['global_rows']).issubset(set(f['rows']))]
            if len(candidates)!=1:raise RuntimeError('Ambiguous original specimen fold')
            source=candidates[0];lookup={int(row):j for j,row in enumerate(source['rows'])}
            selected=np.array([lookup[int(row)] for row in a['global_rows']])
            original=source['counts'][selected][:,a['genes']].toarray()
            np.testing.assert_array_equal(original,a['counts'])
            for j,call in enumerate(entry['calls']):
                if not call.get('arrays'):continue
                with np.load(call['arrays']) as z:b={k:z[k] for k in z.files}
                for key in ('q','cells','global_rows','genes','local_bin','mixture_proportion','effects','prior','cols'):
                    np.testing.assert_array_equal(a[key],b[key])
                w=b['matching_weights'];choose=b['selected']
                np.testing.assert_array_equal(b['expected'],a['expected']*w[:,None])
                if call['context']=='depth_matched':
                    np.testing.assert_array_equal(b['size'],np.broadcast_to(a['size'],original.shape)*w[:,None])
                assert np.all(b['counts']<=original)
                weight_record={}
                if weight_provider is not None and call['context']=='depth_matched':
                    assert len(lookup)==len(source['rows']) and len(np.unique(a['global_rows']))==len(selected)
                    w,weight_record=weight_provider(a,b,source,selected)
                    choose=np.flatnonzero(w>0)
                    if not weight_record['eligible']:
                        records.append(dict(index=i,case_seed=entry['case_seed'],source=entry['source'],context=call['context'],
                            skipped='Unchanged effective-cell/group-support criterion fails',weight_record=weight_record,
                            new_cores=[],decisive_case_qualifies=False))
                        continue
                cols=b['cols'];opts=dict(left_neutral=bool(np.any(b['initial']!=np.log(b['prior']))),
                    right_neutral=bool(np.any(b['terminal']!=0)),mixture_proportion=b['mixture_proportion'])
                resolved=matching_model(original[choose],a['expected'][choose],a['size'],b['q'][choose],b['local_bin'],
                    b['effects'],b['prior'],float(b['change']),w[choose],**opts)
                # Actual boundary messages are authoritative, including their
                # normalization; replace reconstructed boundary conditions.
                resolved['initial']=b['initial'];resolved['terminal']=b['terminal']
                resolved['prior_normalizer']=engine.chain(np.zeros((len(cols),3)),resolved['logT'],
                    b['initial'],b['terminal'])[1]
                fitted=engine.enumerate_states(resolved)
                reset=np.r_[True,np.zeros(len(cols)-1,bool)]
                cores=profile.credible_difference_runs(fitted['path'],fitted['posterior'],reset,diff)
                if call['context']=='initial':
                    error=float(np.max(abs(fitted['posterior']-b['posterior'])));replay_max=max(replay_max,error)
                    operations=64*fitted['states_enumerated']*len(choose)*len(cols)+128
                    budget=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(fitted['log_evidence']))
                    assert error<=budget
                decisive=entry['case_seed']==20261109 and entry['source'].startswith('event_chrom12_') and call['context']=='depth_matched' and bool(cores)
                new_core|=decisive
                record=dict(index=i,case_seed=entry['case_seed'],source=entry['source'],context=call['context'],
                    original_core=call.get('exact_segments',[]),new_cores=[list(c) for c in cores],
                    old_max_nonzero=float(b['posterior'][:,[0,2]].max()),new_max_nonzero=float(fitted['posterior'][:,[0,2]].max()),
                    cells=len(choose),max_count=int(original.max()),finite_support_terms=int(np.sum(original+1)),
                    log_evidence=fitted['log_evidence'],decisive_case_qualifies=bool(decisive))
                if weight_record:record['weight_record']=weight_record
                records.append(record)
                np.savez_compressed(out/f'comparison_{i:03d}_{j:02d}.npz',posterior=fitted['posterior'],path=fitted['path'],
                    membership=fitted['membership'],ll=resolved['ll'],initial=b['initial'],terminal=b['terminal'],
                    matching_weights=w,selected=choose)
        write_json(dict(checks=numerical,comparisons=records,initial_posterior_max_replay_error=replay_max,
            initial_posterior_replay_literal_exact=(replay_max==0.),
            offline_continuation_gate_passed=bool(new_core),full_recovery_run=False,
            old_historical_gate_replaced=False,production_changed=False,scientific_criteria_changed=False),out/'report.json')
        finish_run(out,run,backend='CPU exact finite-kernel averaging and complete state enumeration')
    except BaseException as error:
        finish_run(out,run,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:study(args.out,args.task)


if __name__=='__main__':main()
