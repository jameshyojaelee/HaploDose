"""Research NB gene-concordance scores for frozen depth strata and membership."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import expit
from scipy.stats import binom, nbinom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('conditioned_gene_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def efficient_scores(counts,expected,size,q,strata):
    y,e,r=np.broadcast_arrays(np.asarray(counts,float),np.asarray(expected,float),np.asarray(size,float))
    q,s=np.asarray(q,float),np.asarray(strata)
    if (y.ndim!=2 or q.shape!=(len(y),) or s.shape!=q.shape or s.dtype.kind not in 'iu'
            or not all(np.isfinite(v).all() for v in (y,e,r,q)) or np.any(y<0)
            or np.any(y!=np.floor(y)) or np.any(e<=0) or np.any(r<=0) or np.any((q<0)|(q>1))):
        raise ValueError('Aligned raw gene counts, positive NB means/sizes, fixed q and integer strata required')
    scores=np.zeros(y.shape[1]);information=np.zeros_like(scores);errors=np.zeros_like(scores);records=[]
    eta_solver=load('predictive_gene_effect').weighted_eta
    for level in np.unique(s):
        use=np.flatnonzero(s==level);yy,ee,rr=y[use],e[use],r[use]
        for g in range(y.shape[1]):
            fitted=eta_solver(yy[:,g,None],ee[:,g,None],rr[:,g,None],np.ones(len(use)))
            if not fitted['identified']:
                records.append(dict(stratum=int(level),gene=g,identified=False,reason=fitted['reason']));continue
            p=expit(np.log(ee[:,g])+fitted['eta']-np.log(rr[:,g]));weight=rr[:,g]*p
            if weight.sum()<=0:raise ArithmeticError('Positive count stratum has no representable NB information')
            qbar=float(np.dot(weight,q[use])/weight.sum());centred=q[use]-qbar
            residual=yy[:,g]-(yy[:,g]+rr[:,g])*p
            value=float(np.dot(centred,residual));info=float(np.dot(weight,centred*centred))
            # The derivative of this efficient score with respect to eta is
            # bounded by sum(y+r): NB score derivatives and the weighted-q
            # covariance each contribute at most a quarter of that scale.
            # Add a conservative elementary-operation accumulation bound.
            operations=128*len(use)+256;eps=np.finfo(float).eps
            gamma=operations*eps/(1-operations*eps);magnitude=float(np.sum(yy[:,g]+rr[:,g]))
            error=(fitted['numerical_eta_error']+8*gamma)*magnitude
            scores[g]+=value;information[g]+=info;errors[g]+=error
            records.append(dict(stratum=int(level),gene=g,identified=True,score=value,information=info,
                numerical_score_error=error,null_eta=fitted['eta'],qbar=qbar,
                null_score_residual=fitted['score_residual']))
    return dict(score=scores,information=information,numerical_score_error=errors,null_fits=records)


def concordance(counts,expected,size,q,strata,gene_bins,callable_genes,segment,*,min_genes=8,alpha=.01):
    counts=np.asarray(counts);bins=np.asarray(gene_bins);valid=np.asarray(callable_genes,bool)
    if counts.ndim!=2 or not np.isfinite(counts).all() or np.any(counts<0) or np.any(counts!=np.floor(counts)):
        raise ValueError('Finite nonnegative integer gene observations required')
    if bins.shape!=(counts.shape[1],) or valid.shape!=bins.shape or bins.dtype.kind not in 'iu':
        raise ValueError('Explicit gene/bin/callable correspondence required')
    if min_genes<1 or not 0<alpha<1 or segment['direction'] not in ('lower','higher'):
        raise ValueError('Explicit gene support, alpha and oriented contrast required')
    expected=np.asarray(expected,float);r=np.broadcast_to(np.asarray(size,float),counts.shape)
    if expected.shape!=counts.shape or not np.isfinite(expected).all() or np.any(expected<=0) or not np.isfinite(r).all() or np.any(r<=0):
        raise ValueError('Positive finite aligned NB means and sizes required')
    informative=(counts.sum(axis=0)>=10)&valid&(bins>=segment['first_bin'])&(bins<=segment['last_bin'])
    selected=np.flatnonzero(informative);n=len(selected)
    if n<min_genes:
        return dict(available=False,passed=False,reason='insufficient_informative_genes',n=n,k=0,p=None,
            min_genes=min_genes,alpha=alpha,selected_genes=selected.tolist())
    result=efficient_scores(counts[:,selected],np.asarray(expected)[:,selected],r[:,selected],q,strata)
    signs=np.where(abs(result['score'])>result['numerical_score_error'],np.sign(result['score']),0.)
    signs[result['information']<=0]=0.
    # Existing difference orientation is group0 minus1; positive score is1>0.
    wanted=1 if segment['direction']=='lower' else -1
    k=int(np.sum(signs==wanted));p=float(binom.sf(k-1,n,.5))
    return dict(available=True,passed=bool(p<alpha),reason='nominal_gene_score_sign_test',n=n,k=k,p=p,
        min_genes=min_genes,alpha=alpha,selected_genes=selected.tolist(),score=result['score'].tolist(),
        information=result['information'].tolist(),numerical_score_error=result['numerical_score_error'].tolist(),
        unresolved_signs=int(np.sum(signs==0)),null_fits=result['null_fits'])


def checks():
    y=np.array([[2],[4],[6],[8]]);e=np.full_like(y,2.,dtype=float);r=np.full(1,4.)
    q=np.array([.25,.75,.25,.75]);s=np.array([0,0,1,1]);fitted=efficient_scores(y,e,r,q,s)
    score,info=float(Fraction(36,77)),float(Fraction(41,77))
    assert abs(fitted['score'][0]-score)<=fitted['numerical_score_error'][0]
    operations=256*y.size+128;eps=np.finfo(float).eps;allowance=operations*eps/(1-operations*eps)*max(1.,info)
    assert abs(fitted['information'][0]-info)<=allowance
    theta=np.array([1.5,1.5,3.5,3.5])[:,None]
    def loglik(beta):
        mean=e*theta*(1-q[:,None]+q[:,None]*np.exp(beta))
        return float(nbinom.logpmf(y,r,r/(r+mean)).sum())
    references=[]
    for h in (2.**-6,2.**-8,2.**-10,2.**-12):
        a,b=loglik(h),loglik(-h);gradient=(a-b)/(2*h)
        truncation=h*h*float(np.sum(2*y+r))/24
        rounding=operations*eps/(1-operations*eps)*max(1.,abs(a)+abs(b))/h
        assert abs(gradient-score)<=truncation+rounding
        references.append(dict(step=h,gradient=gradient,error=abs(gradient-score),
            analytic_third_derivative_truncation=truncation,rounding_allowance=rounding))
    other=efficient_scores(y,e,r,1-q,s)
    assert abs(other['score'][0]+fitted['score'][0])<=sum((other['numerical_score_error']+fitted['numerical_score_error']))
    permutation=np.array([3,0,2,1]);permuted=efficient_scores(y[permutation],e[permutation],r,q[permutation],s[permutation])
    assert abs(permuted['score'][0]-fitted['score'][0])<=sum((permuted['numerical_score_error']+fitted['numerical_score_error']))
    constant=efficient_scores(y,e,r,np.array([.1,.1,.9,.9]),s)
    assert abs(constant['score'][0])<=constant['numerical_score_error'][0]
    assert constant['information'][0]<=allowance
    nuisance=efficient_scores([[2],[2],[8],[8]],e,r,q,s)
    assert abs(nuisance['score'][0])<=nuisance['numerical_score_error'][0]
    matrix=np.tile(y,(1,8));expected=np.tile(e,(1,8));sizes=np.full(8,4.)
    segment=dict(first_bin=0,last_bin=0,direction='lower');bins=np.zeros(8,int)
    passing=concordance(matrix,expected,sizes,q,s,bins,np.ones(8,bool),segment)
    assert passing['k']==passing['n']==8 and passing['p']==float(Fraction(1,256)) and passing['passed']
    matrix[:,0]=[4,2,8,6];failing=concordance(matrix,expected,sizes,q,s,bins,np.ones(8,bool),segment)
    assert failing['k']==7 and failing['p']==float(Fraction(9,256)) and not failing['passed']
    gp=np.array([7,0,5,2,4,1,6,3]);reordered=concordance(matrix[:,gp],expected[:,gp],sizes[gp],q,s,bins[gp],np.ones(8,bool),segment)
    assert (reordered['k'],reordered['p'],reordered['passed'])==(failing['k'],failing['p'],failing['passed'])
    assert not concordance(matrix[:,:7],expected[:,:7],sizes[:7],q,s,bins[:7],np.ones(7,bool),segment)['passed']
    zero=efficient_scores(np.zeros_like(y),e,r,q,s)
    np.testing.assert_array_equal(zero['score'],[0.]);np.testing.assert_array_equal(zero['information'],[0.])
    return dict(exact_fraction_score=str(Fraction(36,77)),exact_fraction_information=str(Fraction(41,77)),
        maximum_score_reference_error=float(abs(fitted['score'][0]-score)),numerical_score_budget=float(fitted['numerical_score_error'][0]),
        information_reference_error=float(abs(fitted['information'][0]-info)),information_budget=allowance,
        independent_scipy_logpmf_gradient=references,exact_binomial_tails=True,component_swap=True,
        cell_gene_permutation=True,no_within_stratum_information=True,pure_depth_nuisance_score_zero=True,
        zero_count_boundary_no_invented_rate=True,insufficient_genes_abstain=True,whole_specimen_calibration=False)


def panel(args):
    prior=args.task.parent/'depth_conditioned_genes_20261109/full_candidate_v2';trace_path=prior/'stage_trace.json'
    trace=json.loads(trace_path.read_text());selected=[]
    for i,t in enumerate(trace):
        if t.get('stage')!='gene_shared_confirmation':continue
        for j,c in enumerate(t.get('calls',[])):
            if c.get('context')=='depth_matched' and c.get('arrays') and c.get('exact_segments'):
                selected.append((i,j,t,c))
    if not selected:raise RuntimeError('No frozen matched cores')
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',trace=trace_path,
        rate=Path(__file__).with_name('predictive_gene_effect.py'),previous_veto=prior/'candidate/weight_records.json')
    for i,j,t,c in selected:
        inputs[f'core{i}_{j}']=Path(c['arrays']);inputs[f'lock{i}_{j}']=Path(t['proposal_lock'])
    run=start_run(args.out,inputs,dict(alpha=.01,min_genes=8,informative_count=10),
        dict(scope='Frozen conditional gene-score panel; not full recovery or biological validation'),mode='conditioned_gene_score')
    try:
        numerical=checks()
        write_json(dict(entries=[dict(index=i,call=j,case_seed=t['case_seed'],source=t['source'],
            original_accepted=t['accepted'],segments=c['exact_segments'],arrays=c['arrays'],proposal_lock=t['proposal_lock'])
            for i,j,t,c in selected],source_trace_entries=len(trace)),args.out/'PANEL_LOCKED.json')
        records=[];required=[];decisive=[]
        for i,j,t,c in selected:
            lock_path=Path(t['proposal_lock']);lock_record=file_record(lock_path);lock=json.loads(lock_path.read_text())
            with np.load(c['arrays']) as z:a={k:z[k] for k in z.files}
            params=lock['parameters'];use=np.asarray(params['selected'],int)
            np.testing.assert_array_equal(use,a['selected']);np.testing.assert_array_equal(lock['genes'],a['genes'])
            np.testing.assert_array_equal(lock['rows'],a['global_rows'])
            size=np.broadcast_to(a['size'],a['counts'].shape)[use]
            expected=a['expected'][use]*np.exp(np.asarray(params['offsets'])[use,None])
            bins=a['cols'][a['local_bin']];q=a['q'][use];strata=np.asarray(params['strata'],int)[use]
            for segment in c['exact_segments']:
                result=concordance(a['counts'][use],expected,size,q,strata,bins,np.ones(len(bins),bool),segment)
                key=dict(index=i,call=j,case_seed=t['case_seed'],source=t['source'],segment=segment,
                    original_accepted=t['accepted'],contributing_groups=params['contributing_groups'])
                records.append(dict(**key,result=result))
                if t['accepted']:required.append(result['passed'])
                if t['case_seed']==20261109 and t['source']=='event_chrom12_bins209-217':decisive.append(result['passed'])
            assert file_record(lock_path)==lock_record
        gate=bool(decisive and all(decisive) and all(required))
        write_json(dict(checks=numerical,comparisons=records,source_trace_entries=len(trace),matched_entries=len(selected),
            required_originally_accepted_cores=len(required),decisive109=decisive,offline_continuation_gate_passed=gate,
            event_callable_contract='Executed adapter rejects any event with a noncallable bin; only retained event genes used',
            full_recovery_run=False,production_changed=False,scientific_criteria_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU stratum NB efficient gene scores and fixed nominal sign threshold')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--checks-only',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    if args.checks_only:write_json(checks(),args.out/'report.json')
    else:panel(args)


if __name__=='__main__':main()
