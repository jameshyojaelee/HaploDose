"""Discovery-only regional RNA-rate separation for conditional confirmation."""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import nbinom
from haplodose import profile
from haplodose.clones import _fit_event_mixture
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('predictive_gene_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def weighted_eta(counts,expected,size,weights):
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    w=np.asarray(weights,float)
    if (y.ndim!=2 or w.shape!=(len(y),) or not all(np.isfinite(x).all() for x in (y,e,r,w))
            or np.any(y<0) or np.any(y!=np.floor(y)) or np.any(e<=0) or np.any(r<=0) or np.any((w<0)|(w>1))):
        raise ValueError('Supported integer gene counts, positive expectations/sizes and aligned responsibilities required')
    mass=float(np.sum(w[:,None]*y))
    if mass==0:return dict(identified=False,reason='zero_weighted_observed_count')
    loge,logr=np.log(e),np.log(r)
    def score(eta):return float(np.sum(w[:,None]*(y-(y+r)*expit(loge+eta-logr))))
    lo,hi=-1.,1.
    while score(lo)<=0:lo*=2
    while score(hi)>=0:hi*=2
    eps=np.finfo(float).eps
    eta=float(brentq(score,lo,hi,xtol=eps,rtol=4*eps))
    p=expit(loge+eta-logr)
    information=float(np.sum(w[:,None]*(y+r)*p*(1-p)))
    if information<=0:return dict(identified=False,reason='zero_numerical_information')
    operations=32*y.size+128
    gamma=operations*eps/(1-operations*eps)
    magnitude=float(np.sum(w[:,None]*(y+(y+r)*p)))
    eta_error=4*eps*max(1.,abs(eta))+gamma*max(1.,magnitude)/information
    residual=abs(score(eta))
    assert residual<=information*eta_error
    return dict(identified=True,eta=eta,numerical_eta_error=eta_error,score_residual=residual,
        conditional_information=information,weighted_observed_count=mass)


def parameters(counts,expected,size,q):
    q=np.asarray(q,float)
    first=weighted_eta(counts,expected,size,1-q);second=weighted_eta(counts,expected,size,q)
    if not first['identified'] or not second['identified']:
        return dict(identified=False,components=[first,second],reason='unidentified_component_rate')
    delta=abs(second['eta']-first['eta']);error=first['numerical_eta_error']+second['numerical_eta_error']
    return dict(identified=bool(delta>error),components=[first,second],delta=delta,
        numerical_delta_error=error,reason='discovery_weighted_NB_rates' if delta>error else 'numerically_unresolved_separation')


def joint_parameters(counts,expected,size,step):
    try:fitted=_fit_event_mixture(counts,expected,size,step)
    except RuntimeError as error:
        if not str(error).startswith('Event mixture did not converge:'):raise
        return dict(identified=False,reason='gene_mixture_optimizer_not_certified',failure=str(error))
    if fitted['parameters'] is None or fitted['gain']<=fitted['numerical_gain_budget']:
        return dict(identified=False,reason='unresolved_or_collapsed_gene_mixture')
    q=fitted['responsibility'];eta=fitted['parameters'][:2]
    conditional=parameters(counts,expected,size,q)
    if not conditional['identified']:return conditional
    discrepancy=sum(abs(float(eta[k])-conditional['components'][k]['eta']) for k in (0,1))
    delta=abs(float(eta[1]-eta[0]));error=conditional['numerical_delta_error']+discrepancy
    return dict(identified=bool(delta>error),delta=delta,numerical_delta_error=error,
        reason='joint_discovery_gene_mixture' if delta>error else 'numerically_unresolved_joint_separation',
        responsibility=q.tolist(),mixture_parameters=fitted['parameters'].tolist(),
        log_likelihood_gain=float(fitted['gain']),gain_error_budget=float(fitted['numerical_gain_budget']),
        conditional_rate_discrepancy=discrepancy,solver=fitted['solver'],
        conditional_components=conditional['components'])


def decimal_eta(counts,expected,size,weights,precision):
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    with localcontext() as context:
        context.prec=precision
        rows=[tuple(Decimal.from_float(float(x)) for x in (y[i,g],e[i,g],r[i,g],weights[i]))
            for i in range(len(y)) for g in range(y.shape[1])]
        def score(rate):return sum(w*(yy-(yy+rr)*rate*ee/(rr+rate*ee)) for yy,ee,rr,w in rows)
        lo,hi=Decimal(0),Decimal(1)
        while score(hi)>0:hi*=2
        for _ in range(320):
            mid=(lo+hi)/2
            if score(mid)>0:lo=mid
            else:hi=mid
        return float(((lo+hi)/2).ln())


def checks():
    y=np.array([[0,2],[1,3]]);expected=np.ones((2,2));size=np.full(2,4.);q=np.array([.25,.75])
    fitted=parameters(y,expected,size,q)
    exact=[Fraction(5,4),Fraction(7,4)]
    errors=[]
    for i,w in enumerate((1-q,q)):
        result=fitted['components'][i];reference=math.log(float(exact[i]))
        assert abs(result['eta']-reference)<=result['numerical_eta_error'];errors.append(abs(result['eta']-reference))
    for yy,ee,rr,ww in [(y,expected,size,q),(y,np.array([[.5,3.],[2.,.75]]),np.array([2.5,7.]),q)]:
        actual=weighted_eta(yy,ee,rr,ww);ref80=decimal_eta(yy,ee,rr,ww,80);ref100=decimal_eta(yy,ee,rr,ww,100)
        assert ref80==ref100 and abs(actual['eta']-ref100)<=actual['numerical_eta_error']
        errors.append(abs(actual['eta']-ref100))
    swapped=parameters(y[::-1,::-1],expected[::-1,::-1],size[::-1],q[::-1])
    assert abs(swapped['delta']-fitted['delta'])<=swapped['numerical_delta_error']+fitted['numerical_delta_error']
    reversed_labels=parameters(y,expected,size,1-q)
    assert reversed_labels['delta']==fitted['delta']
    # The actual parameter interface takes no confirmation observation. Hold
    # its saved arguments fixed while perturbing the disjoint test observation.
    discovery=dict(counts=y,expected=expected,size=size,q=q)
    test=y.copy();first=parameters(**discovery);test+=1000
    assert first==parameters(**discovery)
    assert not weighted_eta(np.zeros_like(y),expected,size,q)['identified']
    fixture=np.vstack([np.tile([1,2],(40,1)),np.tile([12,16],(40,1))])
    exposure=np.full(fixture.shape,5.);shape=np.full(2,3.)
    joint=joint_parameters(fixture,exposure,shape,.3)
    assert joint['identified']
    # Independently assemble the full mixture probability with SciPy NB PMFs.
    eta0,eta1,logit=joint['mixture_parameters'];pi=expit(logit)
    null_eta=math.log(float(fixture.mean()/5))
    first=nbinom.logpmf(fixture,shape,shape/(shape+exposure*np.exp(eta0))).sum(axis=1)
    second=nbinom.logpmf(fixture,shape,shape/(shape+exposure*np.exp(eta1))).sum(axis=1)
    alternative=float(np.logaddexp(np.log1p(-pi)+first,np.log(pi)+second).sum())
    null=float(nbinom.logpmf(fixture,shape,shape/(shape+exposure*np.exp(null_eta))).sum())
    operations=256*fixture.size+128
    allowance=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(null)+abs(alternative))
    assert abs((alternative-null)-joint['log_likelihood_gain'])<=allowance
    return dict(independent_fraction_rates=list(map(str,exact)),decimal_80_100_convergence=True,
        maximum_reference_eta_error=max(errors),component_swap_and_permutation=True,
        confirmation_not_a_parameter_argument=True,confirmation_perturbation_parameters_exact=True,
        zero_mass_abstains=True,independent_scipy_joint_gain_error=abs((alternative-null)-joint['log_likelihood_gain']),
        joint_gain_allowance=allowance,kernel_reference=load('averaged_count_matching').checks())


def read_csr(z,prefix):
    return sparse.csr_matrix((z[prefix+'_data'],z[prefix+'_indices'],z[prefix+'_indptr']),shape=tuple(z[prefix+'_shape']))


def study(args):
    root=args.task.parents[1];prior=root/'dev/gene_shared_confirmation_20261109/count_gene_v1'
    trace_path=prior/'stage_trace.json';split_path=prior/'candidate/split_records.json'
    trace=json.loads(trace_path.read_text());splits=json.loads(split_path.read_text())
    entries=[(i,t) for i,t in enumerate(trace) if t.get('stage')=='gene_shared_confirmation'
        and any(c.get('context')=='depth_matched' and c.get('arrays') for c in t.get('calls',[]))]
    assert len(entries)==5
    uniform=root/'dev/outside_event_depth_20261109/uniform_full_v1'
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',trace=trace_path,splits=split_path,
        kernel=Path(__file__).with_name('averaged_count_matching.py'),matching=Path(__file__).with_name('outside_event_depth.py'),
        engine=Path(__file__).with_name('shared_membership_hmm.py'),scorer=Path(profile.__file__))
    if args.polish_membership:inputs['membership_design']=args.task/'MEMBERSHIP_DESIGN.md'
    for row in splits:inputs[f"fold{row['model']}"]=prior/f"candidate/folds/model_{row['model']}.npz"
    for i,t in entries:
        for j,c in enumerate(t['calls']):
            if c.get('arrays'):
                inputs[f'comparison{i}_{j}']=Path(c['arrays'])
                if c['context']=='depth_matched':inputs[f'uniform{i}_{j}']=uniform/f'comparison_{i:03d}_{j:02d}.npz'
    run=start_run(args.out,inputs,dict(posterior=.95,min_bins=5,min_effective_cells=30,polish_membership=args.polish_membership),
        dict(scope='Five fixed proposal/confirmation pairs; no full recovery or genomic validation'),mode='predictive_gene_effect')
    try:
        numerical=checks();folds=[]
        # Proposal phase never accesses full source counts or held-out count
        # arrays. NPZ entries are loaded individually, not unpacked wholesale.
        for row in splits:
            path=inputs[f"fold{row['model']}"]
            with np.load(path) as z:
                folds.append(dict(path=path,case_seed=row['case_seed'],rows=z['global_rows'],proposal=read_csr(z,'prop'),
                    gene_bins=z['gbin'],share=float(z['share']),size=z['gene_size']))
        proposals={};locked=[]
        metadata=('expected','size','q','genes','global_rows','cells','cols','local_bin','mixture_proportion','effects','prior','change','initial','terminal')
        for i,t in entries:
            initial=next(c for c in t['calls'] if c['context']=='initial' and c.get('arrays'))
            with np.load(initial['arrays']) as z:a={k:z[k] for k in metadata}
            candidates=[f for f in folds if f['case_seed']==t['case_seed'] and set(a['global_rows']).issubset(set(f['rows']))]
            if len(candidates)!=1:raise RuntimeError('Ambiguous original specimen fold')
            source=candidates[0];lookup={int(r):j for j,r in enumerate(source['rows'])}
            assert len(lookup)==len(source['rows']) and len(np.unique(a['global_rows']))==len(a['global_rows'])
            selected=np.array([lookup[int(r)] for r in a['global_rows']]);share=source['share']
            sizes=np.broadcast_to(source['size'],(source['proposal'].shape[1],))[a['genes']]
            np.testing.assert_array_equal(a['size'],sizes*(1-share))
            counts=source['proposal'][selected][:,a['genes']].toarray()
            expected=a['expected']*(share/(1-share))
            fitted=(joint_parameters(counts,expected,sizes*share,float(abs(a['effects'][2])))
                if args.polish_membership else parameters(counts,expected,sizes*share,a['q']))
            if fitted['identified'] and args.polish_membership:
                probabilities=np.asarray(fitted['responsibility'])
                bin_expected=np.stack([a['expected'][:,a['local_bin']==j].sum(axis=1) for j in range(len(a['cols']))],axis=1)
                fitted['mixture_proportion']=((bin_expected.T@probabilities)/bin_expected.sum(axis=0)).tolist()
            record=dict(index=i,case_seed=t['case_seed'],source=t['source'],old_step=float(abs(a['effects'][2])),
                rows=a['global_rows'].tolist(),genes=a['genes'].tolist(),**fitted)
            locked.append(record);proposals[i]=(a,source,selected,fitted)
        write_json(dict(proposals=locked,confirmation_count_arrays_opened=False,
            parameter_interface='proposal counts, proposal expectations, proposal NB sizes, discovery responsibilities'),
            args.out/'PROPOSAL_PARAMETERS_LOCKED.json')
        engine=load('shared_membership_hmm');kernel=load('averaged_count_matching');matcher=load('outside_event_depth')
        diff=profile.ProfileParams(posterior=.95,min_bins=5);records=[];gate=False
        for i,t in entries:
            a,source,selected,prediction=proposals[i]
            with np.load(source['path']) as z:source['source']=read_csr(z,'source')
            source['counts']=source['source']-source['proposal']
            original=source['counts'][selected][:,a['genes']].toarray()
            initial_core=False
            for j,c in enumerate(t['calls']):
                if not c.get('arrays'):continue
                with np.load(c['arrays']) as z:b={k:z[k] for k in z.files}
                if c['context']=='initial':np.testing.assert_array_equal(original,b['counts'])
                for key in ('q','genes','global_rows','cols','effects','mixture_proportion'):np.testing.assert_array_equal(a[key],b[key])
                if c['context']=='depth_matched':w,support=matcher.provider('uniform_full')(a,b,source,selected)
                else:
                    w=np.ones(len(original));groups=np.bincount((a['q']>.5).astype(int),minlength=2)
                    support=dict(eligible=bool(groups.min()>=30),effective_cells=groups.tolist(),minimum_cells=30)
                baseline_w=w.copy()
                probabilities=(np.asarray(prediction['responsibility']) if args.polish_membership and prediction['identified'] else a['q'])
                proportions=(np.asarray(prediction['mixture_proportion']) if args.polish_membership and prediction['identified'] else a['mixture_proportion'])
                if args.polish_membership and prediction['identified']:
                    labels=(probabilities>.5).astype(int);groups=np.bincount(labels,minlength=2)
                    if groups.min()<30:
                        support=dict(support,eligible=False,new_hard_group_cells=groups.tolist())
                    elif c['context']=='depth_matched':
                        counts=np.asarray(source['source'].sum(axis=1)).ravel();depth=np.log(np.maximum(counts,1.))
                        depth-=depth.mean();w,edges=matcher.matched_weights(depth[selected],labels,np.ones(len(labels)))
                        effective=[];positive=[]
                        for k in (0,1):
                            v=w[(labels==k)&(w>0)].astype(np.longdouble);positive.append(len(v))
                            effective.append(0. if not len(v) else float(len(v)) if np.all(v==v[0]) else float(v.sum()**2/(v*v).sum()))
                        support=dict(support,eligible=bool(min(positive)>=30 and min(effective)>=30),
                            baseline_uniform_mass=float(baseline_w.sum()),new_mass=float(w.sum()),new_edges=edges.tolist(),
                            effective_cells=effective,positive_cells=positive,new_hard_group_cells=groups.tolist())
                    else:support=dict(support,eligible=True,effective_cells=groups.tolist(),new_hard_group_cells=groups.tolist())
                if not prediction['identified'] or not support['eligible']:
                    records.append(dict(index=i,case_seed=t['case_seed'],source=t['source'],context=c['context'],
                        skipped='Unidentified discovery effect or unchanged cell support',support=support,new_cores=[]));continue
                choose=np.flatnonzero(w>0);cols=a['cols']
                def evaluate(effects,weights,probabilities,proportions):
                    use=np.flatnonzero(weights>0)
                    model=kernel.model(original[use],a['expected'][use],a['size'],probabilities[use],a['local_bin'],
                        effects,a['prior'],float(a['change']),weights[use],left_neutral=False,right_neutral=False,
                        mixture_proportion=proportions)
                    model['initial']=a['initial'];model['terminal']=a['terminal']
                    model['prior_normalizer']=engine.chain(np.zeros((len(cols),3)),model['logT'],a['initial'],a['terminal'])[1]
                    return engine.enumerate_states(model)
                baseline=evaluate(a['effects'],baseline_w,a['q'],a['mixture_proportion'])
                if c['context']=='depth_matched':
                    with np.load(inputs[f'uniform{i}_{j}']) as old:
                        np.testing.assert_array_equal(baseline_w,old['matching_weights']);baseline_post=old['posterior']
                else:baseline_post=b['posterior']
                error=float(np.max(abs(baseline['posterior']-baseline_post)))
                operations=64*baseline['states_enumerated']*int(np.sum(baseline_w>0))*len(cols)+128
                budget=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)*max(1.,abs(baseline['log_evidence']))
                assert error<=budget
                delta=prediction['delta'];fitted=evaluate(np.array([-delta,0.,delta]),w,probabilities,proportions)
                cores=profile.credible_difference_runs(fitted['path'],fitted['posterior'],np.r_[True,np.zeros(len(cols)-1,bool)],diff)
                if c['context']=='initial':initial_core=bool(cores)
                decisive=bool(t['case_seed']==20261109 and t['source'].startswith('event_chrom12_')
                    and c['context']=='depth_matched' and cores and initial_core)
                gate|=decisive
                records.append(dict(index=i,case_seed=t['case_seed'],source=t['source'],context=c['context'],
                    old_step=float(abs(a['effects'][2])),new_step=delta,support=support,
                    baseline_replay_max_error=error,baseline_literal_exact=(error==0.),
                    old_max_nonzero=float(baseline['posterior'][:,[0,2]].max()),
                    new_max_nonzero=float(fitted['posterior'][:,[0,2]].max()),new_cores=[list(x) for x in cores],
                    paired_initial_core=initial_core,decisive_case_qualifies=decisive))
                np.savez_compressed(args.out/f'comparison_{i:03d}_{j:02d}.npz',posterior=fitted['posterior'],
                    path=fitted['path'],membership=fitted['membership'],effects=[-delta,0.,delta],matching_weights=w,selected=choose)
        write_json(dict(checks=numerical,comparisons=records,proposal_lock=file_record(args.out/'PROPOSAL_PARAMETERS_LOCKED.json'),
            original_trace_entries_preserved=len(trace),selected_matched_pairs=5,offline_continuation_gate_passed=bool(gate),
            full_recovery_run=False,old_failed_exact_gates_replaced=False,production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU discovery-only NB rate scores and exact conditional state enumeration')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--polish-membership',action='store_true')
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:study(args)


if __name__=='__main__':main()
