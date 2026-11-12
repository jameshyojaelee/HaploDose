"""Conditional mechanism probe: read-fold versus whole-cell population prediction."""
import argparse
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.stats import nbinom,norm
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


@lru_cache(None)
def load(name):
    spec=importlib.util.spec_from_file_location('cell_population_'+name,
        Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def half(ids,strata,seed):
    ids,strata=np.asarray(ids),np.asarray(strata)
    assert len(np.unique(ids))==len(ids) and ids.shape==strata.shape
    train=np.zeros(len(ids),bool)
    for level in np.unique(strata):
        use=np.flatnonzero(strata==level);assert len(use)>=2
        keys=[hashlib.sha256(f'20261110|{seed}|{int(ids[i])}'.encode()).digest()for i in use]
        assert len(set(keys))==len(keys)
        choose=sorted(range(len(use)),key=lambda j:keys[j])[:len(use)//2]
        train[use[choose]]=True
    assert train.any()and(~train).any()
    return train


def population(y,e,r,strata,model):
    """Unconditional marginal RNA density for cells excluded from model fitting."""
    values=np.zeros(len(y));rounding=0.;predictor=load('state_predictive_confirmation')
    for s,level in enumerate(model['levels']):
        use=strata==level
        if not use.any():continue
        null=[model['null'][s],math.sqrt(model['null'][-1])]
        mean,pi=model['base_mean'][s],model['prevalence'][s]
        if pi in(0.,1.):
            alt=[mean+pi*model['delta'],math.sqrt(model['variance'])];mixture=False
        else:
            alt=[mean,mean+model['delta'],math.log(pi)-math.log1p(-pi),
                math.sqrt(model['variance'])];mixture=True
        a,_=predictor.marginal(y[use],e[use],r[use],alt,mixture,model['order'])
        n,_=predictor.marginal(y[use],e[use],r[use],null,False,model['order'])
        values[use]=a-n
        operations=128*y[use].size*model['order'];eps=np.finfo(float).eps
        assert operations*eps<1
        magnitude=max(1.,float(np.max(abs(np.r_[a,n]))))
        rounding+=operations*eps/(1-operations*eps)*magnitude*int(use.sum())
    return dict(score=float(math.fsum(map(float,values))),rounding_budget=rounding,
        per_cell=values)


def numerics():
    ids=np.array([5,3,8,1,9,2]);strata=np.array([0,0,0,1,1,1]);perm=np.array([4,2,0,5,1,3])
    a=half(ids,strata,17);b=half(ids[perm],strata[perm],17)
    assert set(ids[a])==set(ids[perm][b])and not set(ids[a])&set(ids[~a])
    assert set(ids[a])|set(ids[~a])==set(ids)
    y=np.array([[0,1],[3,2],[2,0],[5,4]]);e=np.tile([1.,2.],(4,1));r=np.full(y.shape,2.)
    strata=np.array([0,0,1,1]);differences=[];budgets=[]
    for pi in([.4,.4],[0.,1.]):
        model=dict(order=96,levels=[0,1],null=[-.04,.08,.03],base_mean=[-.17,-.12],
            delta=.34,variance=.01,prevalence=pi)
        result=population(y,e,r,strata,model);direct=[];allowance=0.
        for i in range(4):
            s=strata[i];densities=[]
            for means,weights,variance in(([model['null'][s]],[1.],model['null'][-1]),
                    ([model['base_mean'][s],model['base_mean'][s]+model['delta']],
                     [1-pi[s],pi[s]],model['variance'])):
                mass,error=0.,0.
                for mean,weight in zip(means,weights):
                    if not weight:continue
                    def integrand(z):
                        mu=e[i]*math.exp(mean+math.sqrt(variance)*z)
                        return float(np.exp(nbinom.logpmf(y[i],r[i],r[i]/(r[i]+mu)).sum())*norm.pdf(z))
                    value,err=quad(integrand,-12.,12.,epsabs=64*np.finfo(float).eps,
                        epsrel=64*np.finfo(float).eps,limit=200)
                    mass+=weight*value;error+=weight*(err+2*norm.sf(12.))
                assert mass>error
                densities.append(mass);allowance+=error/(mass-error)
            direct.append(math.log(densities[1]/densities[0]))
        error=abs(math.fsum(direct)-result['score']);budget=allowance+result['rounding_budget']
        assert error<=budget;differences.append(error);budgets.append(budget)
    return dict(independent_full_PMF_quad_errors=differences,quad_and_rounding_budgets=budgets,
        stable_identity_partition=True,disjoint_exhaustive_cells=True,zero_one_prevalence=True,
        direct_training_input_excludes_test_cells=True,upstream_selection_independent=False,
        existing_predictor_checks=load('state_predictive_confirmation').numerical_checks())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--case',type=int);parser.add_argument('--numerics',action='store_true')
    parser.add_argument('--summary',action='store_true');args=parser.parse_args();args.out.mkdir(exist_ok=False)
    cases=json.loads((args.task/'CASES.json').read_text())
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md',cases=args.task/'CASES.json',
        fit=Path(__file__).with_name('profiled_prevalence.py'),
        original_fit=Path(__file__).with_name('depth_state_predictive.py'),
        predictor=Path(__file__).with_name('state_predictive_confirmation.py'))
    prepared=None
    if args.case is not None:
        case=cases[args.case];inputs.update(models=Path(case['models']),trace=Path(case['trace']))
        prepared=load('state_predictive_confirmation').read_proposal(case)
        inputs.update({name:Path(prepared[5][name])for name in('fold','confirmation','proposal_lock')})
        inputs.update({f'baseline_{order}':Path(case['predictions'])/f'prediction_{order}.npz'for order in(24,48,96)})
    if not args.numerics:
        inputs['numerics']=args.task/'numerics_v1/report.json'
        inputs['numeric_manifest']=args.task/'numerics_v1/run_manifest.json'
    if args.summary:
        inputs.update({f'case{i}':args.task/f'case_{i:02d}_v1/report.json'for i in range(len(cases))})
    rec=start_run(args.out,inputs,dict(case=args.case,orders=[24,48,96],allocation_key='20261110|seed|global_row'),
        dict(scope='Conditional exposed-case mechanism probe; upstream selection fixed,biological n=0'),
        mode='cell_population_confirmation')
    try:
        if args.numerics:write_json(numerics(),args.out/'report.json')
        elif args.summary:
            rows=[json.loads(inputs[f'case{i}'].read_text())for i in range(len(cases))]
            passed=all(all(d.get('resolved',False)and(d['score']>0 if c['kind']=='original_recovery'else d['score']<0)
                for d in row['directions'])for c,row in zip(cases,rows))
            write_json(dict(cases=rows,mechanism_separation_passed=passed,
                full_workflow_run=False,independent_confirmation_established=False,
                false_split_control_established=False,production_changed=False,
                decision='Permit a separately specified full cell-discovery workflow experiment'if passed else
                    'Reject this correction or resolve declared numerics; no pipeline expansion'),args.out/'report.json')
        else:
            assert json.loads(inputs['numerics'].read_text())['stable_identity_partition']
            assert json.loads(inputs['numeric_manifest'].read_text())['status']=='complete'
            predictor=load('state_predictive_confirmation');method=load('profiled_prevalence')
            D,ED,rD,_,metadata,provenance,EC,rC=prepared
            with np.load(provenance['confirmation'])as z:C=z['counts'][metadata['selected']]
            with np.load(provenance['fold'])as z:
                full=sparse.csr_matrix((z['source_data'],z['source_indices'],z['source_indptr']),shape=z['source_shape'])
            Y=full[metadata['cells']][:,metadata['genes']].toarray()[metadata['selected']]
            np.testing.assert_array_equal(Y,D+C)
            lock=json.loads(Path(provenance['proposal_lock']).read_text())
            strata=np.asarray(lock['parameters']['strata'],int)[metadata['selected']]
            q=metadata['q'][metadata['selected']];ids=metadata['global_rows'][metadata['selected']]
            old_models=json.loads(inputs['models'].read_text())['models'];baseline=[]
            assert [m['order']for m in old_models]==[24,48,96]
            for old in old_models:
                value=method.predict(D,C,ED,EC,rD,rC,strata,old)
                with np.load(Path(case['predictions'])/f"prediction_{old['order']}.npz")as z:
                    for key in('per_cell','discovery_membership','joint_membership'):
                        np.testing.assert_array_equal(value[key],z[key])
                baseline.append(dict(order=old['order'],score=value['score']))
            E,R=ED+EC,rD+rC;mask=half(ids,strata,case['seed']);locked=[]
            for direction,train in enumerate((mask,~mask)):
                row=dict(direction=direction,train_rows=ids[train].tolist(),test_rows=ids[~train].tolist(),models=[])
                try:
                    for order in(24,48,96):
                        print(json.dumps(dict(case=args.case,direction=direction,order=order,status='training_only_fit')),flush=True)
                        old=load('depth_state_predictive').fit(Y[train],E[train],R[train],strata[train],q[train],order)
                        row['models'].append(method.fit(Y[train],E[train],R[train],strata[train],q[train],old,shared=True))
                except(ArithmeticError,ValueError)as error:row['unresolved_fit']=str(error)
                locked.append(row)
            write_json(dict(directions=locked,upstream_selection_fixed=True,
                full_pipeline_cell_independence=False),args.out/'TRAINING_MODELS_LOCKED.json')
            directions=[]
            for item,train in zip(locked,(mask,~mask)):
                row=dict(direction=item['direction'],training_cells=int(train.sum()),test_cells=int((~train).sum()))
                if 'unresolved_fit'in item:row.update(resolved=False,reason=item['unresolved_fit'])
                else:
                    values=[]
                    for model in item['models']:
                        value=population(Y[~train],E[~train],R[~train],strata[~train],model)
                        np.save(args.out/f"population_{item['direction']}_{model['order']}.npy",value.pop('per_cell'))
                        values.append(dict(order=model['order'],**value))
                    diagnostic=2*max(abs(values[i]['score']-values[i-1]['score'])for i in(1,2))+sum(v['rounding_budget']for v in values)
                    row.update(score=values[-1]['score'],numerical_diagnostic=diagnostic,
                        resolved=bool(abs(values[-1]['score'])>diagnostic),quadrature=values)
                directions.append(row)
            write_json(dict(case=case,baseline_read_fold=baseline,baseline_prediction_arrays_exact=9,
                directions=directions,provenance=provenance,upstream_selection_fixed=True,
                population_error_estimate=False,full_workflow_run=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU independent-cell marginal RNA probe; frozen upstream proposal')
    except BaseException as error:finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
