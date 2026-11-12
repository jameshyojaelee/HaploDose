"""Research predictive NB confirmation with discovery-only depth strata."""
import argparse
from functools import lru_cache
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.stats import nbinom
from haplodose import profile
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('conditioned_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def prepare(counts,expected,size,q,depth,*,n_strata=8,min_cells=30):
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    q,depth=np.asarray(q,float),np.asarray(depth,float)
    if (y.ndim!=2 or q.shape!=(len(y),) or depth.shape!=q.shape or n_strata<1 or min_cells<1
            or not all(np.isfinite(v).all() for v in (y,e,r,q,depth))
            or np.any(y<0) or np.any(y!=np.floor(y)) or np.any(e<=0) or np.any(r<=0)
            or np.any((q<0)|(q>1))):
        raise ValueError('Aligned integer counts, positive exposure/size and finite discovery metadata required')
    edges=np.quantile(depth,np.linspace(0,1,n_strata+1)[1:-1])
    strata=np.searchsorted(edges,depth);labels=(q>.5).astype(int)
    offsets=np.zeros(len(y));pi=np.zeros_like(e);selected=np.zeros(len(y),bool);records=[]
    for s in range(n_strata):
        use=np.flatnonzero(strata==s);groups=np.bincount(labels[use],minlength=2)
        if not len(use) or groups.min()==0:
            records.append(dict(stratum=s,groups=groups.tolist(),included=False,reason='no_hard_group_overlap'));continue
        fit=load('predictive_gene_effect').weighted_eta(y[use],e[use],r[use],np.ones(len(use)))
        if not fit['identified']:
            records.append(dict(stratum=s,groups=groups.tolist(),included=False,reason=fit['reason']));continue
        offsets[use]=fit['eta'];pi[use]=(e[use].T@q[use])/e[use].sum(axis=0);selected[use]=True
        records.append(dict(stratum=s,groups=groups.tolist(),included=True,fit=fit))
    groups=np.bincount(labels[selected],minlength=2)
    return dict(identified=bool(groups.min()>=min_cells),selected=np.flatnonzero(selected).tolist(),
        offsets=offsets.tolist(),mixture_proportion=pi.tolist(),edges=edges.tolist(),strata=strata.tolist(),
        strata_records=records,contributing_groups=groups.tolist(),effective_cells=groups.tolist(),
        minimum_cells=min_cells,reason='discovery_only_stratum_offsets' if groups.min()>=min_cells else 'insufficient_common_support')


def model(counts,expected,size,q,bins,effects,prior,change,parameters,initial,terminal):
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    if np.any(y!=np.floor(y)) or np.any(y<0):raise ValueError('Raw integer confirmation counts required')
    use=np.asarray(parameters['selected'],int);local=np.asarray(bins,int);effects=np.asarray(effects,float)
    y,e,r=y[use],e[use]*np.exp(np.asarray(parameters['offsets'])[use,None]),r[use]
    pi=np.asarray(parameters['mixture_proportion'])[use];probability=np.asarray(q,float)[use]
    centre=np.log1p(pi[:,:,None]*np.expm1(effects)[None,None])
    shifts=np.stack([-centre,effects[None,None,:]-centre])
    ratio=e/(e+r)
    ll=y[None,:,:,None]*shifts-(y+r)[None,:,:,None]*np.log1p(ratio[None,:,:,None]*np.expm1(shifts))
    length=int(local.max())+1
    if set(local)!=set(range(length)):raise ValueError('Complete event gene/bin mapping required')
    grouped=np.stack([ll[:,:,local==j,:].sum(axis=2) for j in range(length)],axis=2)
    engine=load('shared_membership_hmm');logT=np.log(engine.transition(prior,change))
    return dict(ll=grouped,q=probability,logT=logT,initial=initial,terminal=terminal,
        prior_normalizer=engine.chain(np.zeros((length,3)),logT,initial,terminal)[1])


def checks():
    y=np.array([[2,3],[4,1]]);e=y.astype(float);r=np.array([3.,5.]);q=np.array([.1,.8])
    effects=np.array([-.4,0.,.4]);prior=np.array([.15,.7,.15]);change=.02
    initial=np.log(prior);terminal=np.zeros(3);bins=np.array([0,1])
    fit=prepare(y,e,r,q,[1.,2.],n_strata=1,min_cells=1)
    assert fit['identified']
    assert max(abs(np.asarray(fit['offsets'])))<=fit['strata_records'][0]['fit']['numerical_eta_error']
    # The old-model equivalence fixture specifies the analytically exact zero
    # offset; the estimated root is checked separately against its error bound.
    fit['offsets']=[0.,0.]
    new=model(y,e,r,q,bins,effects,prior,change,fit,initial,terminal)
    old=load('gene_shared_confirmation').gene_model(y,e,r,q,bins,effects,prior,change,
        left_neutral=False,right_neutral=False,mixture_proportion=(e.T@q)/e.sum(axis=0))
    allowance=512*y.size*np.finfo(float).eps*max(1.,float(np.abs(old['ll']).sum()))
    error=float(np.max(abs(old['ll']-new['ll'])));assert error<=allowance
    # Independent NB probability assembly, integrating all labels and paths.
    engine=load('shared_membership_hmm');actual=engine.enumerate_states(new)
    T=engine.transition(prior,change);joint=[];paths=list(itertools.product(range(3),repeat=2))
    pi=np.asarray(fit['mixture_proportion']);null=nbinom.logpmf(y,r,r/(r+e)).sum()
    for path in paths:
        mass=0.
        for z in itertools.product(range(2),repeat=2):
            multiplier=np.exp(np.asarray(z)[:,None]*effects[list(path)][None,:])/(1-pi+pi*np.exp(effects[list(path)])[None,:])
            alternative=nbinom.logpmf(y,r,r/(r+e*multiplier)).sum()
            pz=np.prod(np.where(z,q,1-q));mass+=pz*np.exp(alternative-null)
        joint.append(prior[path[0]]*T[path[0],path[1]]*mass)
    total=sum(joint);reference=np.zeros((2,3))
    for path,p in zip(paths,joint):
        for j,s in enumerate(path):reference[j,s]+=p/total
    operations=512*y.size*len(paths)*4
    budget=operations*np.finfo(float).eps*max(1.,abs(null),abs(np.log(total)))
    posterior_error=float(np.max(abs(reference-actual['posterior'])))
    evidence_error=abs(np.log(total)-actual['log_evidence'])
    assert max(posterior_error,evidence_error)<=budget
    pi=np.asarray(fit['mixture_proportion']);d=.4
    mixture=((1-q)[:,None]+q[:,None]*np.exp(d))/(1-pi+pi*np.exp(d))
    assert np.max(abs((e*mixture).sum(axis=0)-e.sum(axis=0)))<=budget
    permutation=np.array([1,0]);gene=np.array([1,0])
    perm=prepare(y[permutation][:,gene],e[permutation][:,gene],r[gene],q[permutation],[2.,1.],n_strata=1,min_cells=1)
    reordered=engine.enumerate_states(model(y[permutation][:,gene],e[permutation][:,gene],r[gene],q[permutation],bins[gene],
        effects,prior,change,perm,initial,terminal))
    assert np.max(abs(reordered['posterior']-actual['posterior']))<=budget
    test=y.copy();locked=prepare(y,e,r,q,[1.,2.],n_strata=1,min_cells=1);test+=100
    assert locked==prepare(y,e,r,q,[1.,2.],n_strata=1,min_cells=1)
    assert not prepare(y,e,r,q,[1.,2.],n_strata=2,min_cells=1)['identified']
    return dict(old_single_stratum_log_ratio_error=error,old_model_allowance=allowance,
        independent_scipy_posterior_error=posterior_error,independent_evidence_error=evidence_error,
        enumeration_allowance=budget,stratum_mean_conservation=True,cell_gene_permutation=True,
        no_overlap_abstains=True,confirmation_perturbation_parameters_exact=True,
        NB_score_reference=load('predictive_gene_effect').checks())


def study(args):
    root=args.task.parents[1];prior=root/'dev/gene_shared_confirmation_20261109/count_gene_v1'
    trace_path=prior/'stage_trace.json';split_path=prior/'candidate/split_records.json'
    trace=json.loads(trace_path.read_text());splits=json.loads(split_path.read_text())
    entries=[(i,t) for i,t in enumerate(trace) if t.get('stage')=='gene_shared_confirmation'
        and any(c.get('context')=='depth_matched' and c.get('arrays') for c in t.get('calls',[]))]
    assert len(entries)==5
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',trace=trace_path,splits=split_path,
        engine=Path(__file__).with_name('shared_membership_hmm.py'),
        gene=Path(__file__).with_name('gene_shared_confirmation.py'),
        rate=Path(__file__).with_name('predictive_gene_effect.py'),
        kernel=Path(__file__).with_name('averaged_count_matching.py'),
        matching=Path(__file__).with_name('outside_event_depth.py'),scorer=Path(profile.__file__))
    for row in splits:inputs[f"fold{row['model']}"]=prior/f"candidate/folds/model_{row['model']}.npz"
    for i,t in entries:
        for j,c in enumerate(t['calls']):
            if c.get('arrays'):inputs[f'comparison{i}_{j}']=Path(c['arrays'])
    run=start_run(args.out,inputs,dict(posterior=.95,min_bins=5,min_cells=30,strata=8),
        dict(scope='Fixed conditional depth-model experiment, not full recovery or biological validation'),mode='depth_conditioned_genes')
    try:
        numerical=checks();reader=load('predictive_gene_effect');folds=[]
        for row in splits:
            with np.load(inputs[f"fold{row['model']}"]) as z:
                folds.append(dict(path=inputs[f"fold{row['model']}"],case_seed=row['case_seed'],rows=z['global_rows'],
                    proposal=reader.read_csr(z,'prop'),gene_bins=z['gbin'],share=float(z['share']),size=z['gene_size']))
        records=[];prepared=[]
        metadata=('expected','size','q','genes','global_rows','cols','local_bin','mixture_proportion','effects','prior','change','initial','terminal')
        for i,t in entries:
            initial=next(c for c in t['calls'] if c['context']=='initial' and c.get('arrays'))
            with np.load(initial['arrays']) as z:a={k:z[k] for k in metadata}
            matches=[f for f in folds if f['case_seed']==t['case_seed'] and set(a['global_rows']).issubset(set(f['rows']))]
            if len(matches)!=1:raise RuntimeError('Ambiguous proposal fold')
            f=matches[0];lookup={int(v):j for j,v in enumerate(f['rows'])}
            assert len(lookup)==len(f['rows']) and len(np.unique(a['global_rows']))==len(a['global_rows'])
            selected=np.array([lookup[int(v)] for v in a['global_rows']]);share=f['share']
            sizes=np.broadcast_to(f['size'],(f['proposal'].shape[1],))[a['genes']]
            np.testing.assert_array_equal(a['size'],sizes*(1-share))
            y=f['proposal'][selected][:,a['genes']].toarray();e=a['expected']*share/(1-share)
            outside=~np.isin(f['gene_bins'],a['cols'])
            if not outside.any():raise RuntimeError('No outside-event discovery depth available')
            depth=np.log(np.maximum(np.asarray(f['proposal'][selected][:,outside].sum(axis=1)).ravel(),1.))
            fitted=prepare(y,e,sizes*share,a['q'],depth)
            records.append(dict(index=i,case_seed=t['case_seed'],source=t['source'],rows=a['global_rows'].tolist(),**fitted))
            prepared.append((i,t,a,f,selected,fitted))
        write_json(dict(proposals=records,confirmation_count_arrays_opened=False,full_source_CSR_opened=False),
            args.out/'PROPOSAL_PARAMETERS_LOCKED.json')
        engine=load('shared_membership_hmm');kernel=load('averaged_count_matching');matcher=load('outside_event_depth')
        comparisons=[];passed=[];diff=profile.ProfileParams(posterior=.95,min_bins=5)
        for i,t,a,f,selected,fitted in prepared:
            with np.load(f['path']) as z:full=reader.read_csr(z,'source')
            counts=(full-f['proposal'])[selected][:,a['genes']].toarray()
            controls=[]
            for j,c in enumerate(t['calls']):
                if not c.get('arrays'):continue
                with np.load(c['arrays']) as z:b={k:z[k] for k in z.files}
                for key in ('q','genes','global_rows','cols','effects','mixture_proportion'):np.testing.assert_array_equal(a[key],b[key])
                if c['context']=='initial':
                    np.testing.assert_array_equal(counts,b['counts']);use=np.arange(len(counts))
                else:
                    source=dict(source=full,proposal=f['proposal'],gene_bins=f['gene_bins'])
                    # This asserts every original matching weight exactly.
                    matcher.provider('uniform_full')(a,b,source,selected);use=b['selected']
                sizes=b['size'][use] if b['size'].ndim==2 else b['size']
                control=kernel.load('gene_shared_confirmation').gene_model(
                    b['counts'][use],b['expected'][use],sizes,b['q'][use],b['local_bin'],b['effects'],b['prior'],float(b['change']),
                    left_neutral=False,right_neutral=False,mixture_proportion=a['mixture_proportion'])
                control['initial']=a['initial'];control['terminal']=a['terminal']
                control['prior_normalizer']=engine.chain(np.zeros((len(a['cols']),3)),control['logT'],a['initial'],a['terminal'])[1]
                replay=engine.enumerate_states(control);error=float(np.max(abs(replay['posterior']-b['posterior'])))
                operations=128*replay['states_enumerated']*len(counts)*len(a['cols'])+256
                budget=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(replay['log_evidence']))
                assert error<=budget
                controls.append(dict(context=c['context'],max_error=error,literal_exact=(error==0.)))
            cores=[];maximum=None
            if fitted['identified']:
                candidate=model(counts,a['expected'],a['size'],a['q'],a['local_bin'],a['effects'],a['prior'],float(a['change']),fitted,a['initial'],a['terminal'])
                result=engine.enumerate_states(candidate)
                cores=profile.credible_difference_runs(result['path'],result['posterior'],np.r_[True,np.zeros(len(a['cols'])-1,bool)],diff)
                maximum=float(result['posterior'][:,[0,2]].max())
                np.savez_compressed(args.out/f'comparison_{i:03d}.npz',posterior=result['posterior'],path=result['path'],
                    membership=result['membership'],selected=fitted['selected'])
            original_initial=next(c for c in t['calls'] if c['context']=='initial' and c.get('arrays'))
            with np.load(original_initial['arrays']) as z:
                old_cores=profile.credible_difference_runs(z['path'],z['posterior'],np.r_[True,np.zeros(len(a['cols'])-1,bool)],diff)
            qualified=bool(old_cores and cores);passed.append(qualified)
            comparisons.append(dict(index=i,case_seed=t['case_seed'],source=t['source'],context='depth_conditioned',
                controls=controls,identified=fitted['identified'],contributing_groups=fitted['contributing_groups'],
                original_initial_core=[list(x) for x in old_cores],new_cores=[list(x) for x in cores],
                maximum_nonzero_support=maximum,conditional_gate_passed=qualified))
        write_json(dict(checks=numerical,comparisons=comparisons,selected_pairs=5,original_trace_entries_preserved=len(trace),
            proposal_lock=file_record(args.out/'PROPOSAL_PARAMETERS_LOCKED.json'),offline_continuation_gate_passed=all(passed),
            full_recovery_run=False,scientific_criteria_changed=False,production_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU discovery-conditioned NB likelihood and exact shared-state enumeration')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:study(args)


if __name__=='__main__':main()
