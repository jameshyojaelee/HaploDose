"""Discovery-likelihood prevalence for the depth-conditioned predictive screen."""
import argparse
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.integrate import quad
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import nbinom, norm
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


@lru_cache(None)
def load(name):
    spec = importlib.util.spec_from_file_location('likelihood_prevalence_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parameters(alpha, delta, eta, variance):
    centre = np.logaddexp(-np.logaddexp(0.,eta),-np.logaddexp(0.,-eta)+delta)
    return np.array([alpha-centre,alpha+delta-centre,eta,variance])


def objective(x, y, e, r, strata, order):
    k = (len(x)-2)//2
    assert len(x)==2*k+2
    gradient,total = np.zeros_like(x),0.
    for s in range(k):
        use = strata==s
        alpha,delta,eta,variance = x[s],x[k],x[k+1+s],x[-1]
        value,g = load('nb_state_guard').variance_objective(
            parameters(alpha,delta,eta,variance),y[use],e[use],r[use],order,True)
        weight = use.sum()/len(y)
        tilted,pi = expit(eta+delta),expit(eta)
        gradient[s] = weight*(g[0]+g[1])
        gradient[k] += weight*(-tilted*g[0]+(1-tilted)*g[1])
        gradient[k+1+s] = weight*(g[2]-(tilted-pi)*(g[0]+g[1]))
        gradient[-1] += weight*g[-1]
        total += weight*value
    return float(total),gradient


def fit(D,E,r,strata,q,old):
    levels,index = np.unique(strata,return_inverse=True)
    assert levels.tolist()==old['levels']
    k,order = len(levels),old['order']
    pi = np.array(old['prevalence'])
    np.testing.assert_array_equal(pi,np.array([q[index==s].mean()for s in range(k)]))
    eta = np.log(pi)-np.log1p(-pi)
    anchor = np.r_[old['alternative'][:-2],old['alternative'][-2],eta,old['alternative'][-1]]
    old_value,anchor_gradient = objective(anchor,D,E,r,index,order)
    operations = 256*D.size*order
    eps = np.finfo(float).eps
    gamma = operations*eps/(1-operations*eps)
    budget = gamma*max(1.,abs(old_value),abs(old['alternative_objective']))
    assert abs(old_value-old['alternative_objective'])<=budget
    z = np.log((D.sum(axis=1)+.5)/(E.sum(axis=1)+.5))
    spread = float(np.diff(np.quantile(z,[.2,.8]))[0])
    null = np.asarray(old['null'])
    starts = [anchor]+[np.r_[null[:-1],d,eta,v]
        for d in (0.,spread,-spread)for v in (0.,null[-1])]
    options = dict(ftol=eps,gtol=np.sqrt(eps),maxiter=1000,maxls=50)
    fits = [minimize(objective,x,args=(D,E,r,index,order),jac=True,
        bounds=[(None,None)]*(2*k+1)+[(0.,None)],method='L-BFGS-B',options=options)
        for x in starts]
    good = [f for f in fits if f.success and np.isfinite(f.fun)and np.isfinite(f.x).all()]
    if not good:
        raise ArithmeticError('Discovery-only free-prevalence optimization unresolved')
    best = min(good,key=lambda f:f.fun)
    if best.fun>old_value+budget:
        raise ArithmeticError('Expanded model fails its feasible fixed-prevalence anchor')
    gradient = objective(best.x,D,E,r,index,order)[1]
    projected = gradient.copy()
    if best.x[-1]==0:projected[-1]=min(0.,projected[-1])
    return dict(order=order,levels=levels.tolist(),null=old['null'],
        alpha=best.x[:k].tolist(),delta=float(best.x[k]),eta=best.x[k+1:-1].tolist(),
        prevalence=expit(best.x[k+1:-1]).tolist(),variance=float(best.x[-1]),
        alternative_objective=float(best.fun),fixed_prevalence_objective=old_value,
        discovery_objective_gain=float(old_value-best.fun),
        anchor_eta_gradient=anchor_gradient[k+1:-1].tolist(),
        maximum_projected_gradient=float(np.max(abs(projected))),
        successful_starts=len(good),starts=[dict(success=bool(f.success),
            objective=float(f.fun),iterations=int(f.nit),message=str(f.message))for f in fits],
        optimization='Best finite successful local L-BFGS-B fit; global optimum not certified',
        prevalence_source='Discovery marginal likelihood; old q used only to initialize')


def predict(D,C,ED,EC,rD,rC,strata,model):
    per_cell,qd,qj = [np.zeros(len(D))for _ in range(3)]
    rounding = 0.
    for s,level in enumerate(model['levels']):
        use = strata==level
        null = [model['null'][s],math.sqrt(model['null'][-1])]
        alternative = parameters(model['alpha'][s],model['delta'],model['eta'][s],model['variance'])
        alternative[-1] = math.sqrt(alternative[-1])
        result = load('state_predictive_confirmation').predictive(
            D[use],C[use],ED[use],EC[use],rD[use],rC[use],null,alternative,model['order'])
        per_cell[use] = result['per_cell']
        qd[use],qj[use] = result['mixture_membership_discovery'],result['mixture_membership_joint']
        rounding += result['rounding_budget']
    def effective(x):
        return float(x.sum()**2/np.sum(x*x))if np.any(x)else 0.
    return dict(score=float(math.fsum(map(float,per_cell))),rounding_budget=rounding,
        per_cell=per_cell,discovery_membership=qd,joint_membership=qj,
        discovery_expected_members=[float((1-qd).sum()),float(qd.sum())],
        discovery_effective_members=[effective(1-qd),effective(qd)])


def checks():
    D = np.array([[1,2],[3,0],[0,1],[2,3]])
    E = np.tile([1.,2.],(4,1))
    r = np.full(D.shape,2.)
    strata = np.array([0,0,1,1])
    x = np.array([.1,-.2,.3,-1.1,.4,.04])
    value,gradient = objective(x,D,E,r,strata,96)
    logdensity,direct = 0.,np.zeros_like(x)
    allowance = 0.
    for i in range(4):
        s = strata[i]
        alpha,delta,eta,v = x[s],x[2],x[3+s],x[-1]
        pi,tilted = expit(eta),expit(eta+delta)
        centre = math.log((1-pi)+pi*math.exp(delta))
        means = [alpha-centre,alpha+delta-centre]
        mass,derivative,errors = 0.,np.zeros(4),np.zeros(5)
        M = float(np.sum(D[i]+r[i]))
        for component,weight in enumerate((1-pi,pi)):
            def integrand(z,kind):
                mu = E[i]*math.exp(means[component]+math.sqrt(v)*z)
                probability = mu/(r[i]+mu)
                score = float(np.sum(D[i]-(D[i]+r[i])*probability))
                hessian = -float(np.sum((D[i]+r[i])*probability*(1-probability)))
                density = float(np.exp(nbinom.logpmf(D[i],r[i],r[i]/(r[i]+mu)).sum())*norm.pdf(z))
                factors = (1.,score,(component-tilted)*score,
                    component-pi-(tilted-pi)*score,.5*(score*score+hessian))
                return density*factors[kind]
            bounds = (1.,M,M,1+M,.5*(M*M+M/4))
            for kind in range(5):
                part,error = quad(lambda z:integrand(z,kind),-12.,12.,
                    epsabs=256*np.finfo(float).eps,epsrel=256*np.finfo(float).eps,limit=200)
                errors[kind] += weight*(error+2*norm.sf(12.)*bounds[kind])
                if kind==0:mass+=weight*part
                else:derivative[kind-1]+=weight*part
        assert mass>errors[0]
        logdensity += math.log(mass)
        direct[s] -= derivative[0]/mass/4
        direct[2] -= derivative[1]/mass/4
        direct[3+s] -= derivative[2]/mass/4
        direct[-1] -= derivative[3]/mass/4
        allowance += errors[0]/(mass-errors[0])/4
        allowance += float(np.sum(errors[1:]/(mass-errors[0])+
            abs(derivative)*errors[0]/(mass*(mass-errors[0]))))/4
    baseline = float(nbinom.logpmf(D,r,r/(r+E)).sum())
    expected = -(logdensity-baseline)/4
    operations = 512*96*D.size
    eps = np.finfo(float).eps
    rounding = operations*eps/(1-operations*eps)*max(1.,abs(baseline),abs(logdensity))
    budget = allowance+rounding
    assert max(abs(expected-value),np.max(abs(direct-gradient)))<=budget
    # At fixed component locations, prevalence score is pi minus posterior
    # mean membership. Alpha's chain term must cancel to recover that identity.
    from fractions import Fraction
    eta = math.log(1/3)
    delta = math.log(2.)
    alpha = math.log(1.25)
    fixture = np.array([alpha,delta,eta,0.])
    v,g = objective(fixture,D,E,r,np.zeros(4,int),96)
    responsibilities = []
    for row in D:
        probabilities = []
        for multiplier in (1,2):
            p = Fraction(1)
            for count,expected_mean in zip(row,(1,2)):
                q = Fraction(2,2+expected_mean*multiplier)
                p *= math.comb(int(count)+1,int(count))*q**2*(1-q)**int(count)
            probabilities.append(p)
        responsibilities.append(probabilities[1]/(3*probabilities[0]+probabilities[1]))
    exact = float(Fraction(1,4)-sum(responsibilities)/4)
    pi,tilted = expit(eta),expit(eta+delta)
    identity_error = abs(g[2]+(tilted-pi)*g[0]-exact)
    assert identity_error<=budget and abs(exact)>budget
    # Exact nesting, label and stratum/cell symmetries, including variance0.
    for variance in (0.,.04):
        trial = x.copy();trial[-1]=variance
        a = objective(trial,D,E,r,strata,96)
        swap = trial.copy();swap[2]*=-1;swap[3:5]*=-1
        b = objective(swap,D,E,r,strata,96)
        changed = b[1].copy();changed[2]*=-1;changed[3:5]*=-1
        assert abs(a[0]-b[0])<=budget and np.max(abs(a[1]-changed))<=budget
        perm = [3,0,2,1]
        b = objective(trial,D[perm],E[perm],r[perm],strata[perm],96)
        assert abs(a[0]-b[0])<=budget and np.max(abs(a[1]-b[1]))<=budget
        exchange = trial[[1,0,2,4,3,5]]
        b = objective(exchange,D,E,r,1-strata,96)
        assert abs(a[0]-b[0])<=budget and np.max(abs(a[1]-b[1][[1,0,2,4,3,5]]))<=budget
    nested = x.copy();nested[2]=0
    local = load('depth_state_predictive')
    old_pi = expit(x[3:5])
    old_x = np.r_[x[:2],0.,x[-1]]
    a = local.objective(old_x,D,E,r,strata,old_pi,96,True)
    b = objective(nested,D,E,r,strata,96)
    assert abs(a[0]-b[0])<=budget and np.max(abs(b[1][3:5]))<=budget
    return dict(independent_full_PMF_quad_objective_error=abs(expected-value),
        independent_full_PMF_quad_gradient_error=float(np.max(abs(direct-gradient))),
        quad_and_rounding_budget=budget,exact_fraction_prevalence_score=exact,
        exact_fraction_prevalence_score_error=identity_error,component_swap=True,
        cell_permutation=True,stratum_permutation=True,nested_null=True,zero_variance_boundary=True,
        old_numerical_evidence='Reused unchanged depth_state_predictive/state_predictive implementations; new prevalence derivative tested here',
        biological_calibration_established=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--case',type=int)
    parser.add_argument('--numerics',action='store_true')
    parser.add_argument('--summary',action='store_true')
    parser.add_argument('--profile-prevalence',action='store_true')
    parser.add_argument('--case-paths',type=Path)
    parser.add_argument('--numerical-report',type=Path)
    parser.add_argument('--shared-prevalence',action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    cases = json.loads((args.task/'CASES.json').read_text())
    prior = args.task.parent/'adjacent_depth_state_20261110'
    inputs = dict(script=Path(__file__),design=args.task/'DESIGN.md',cases=args.task/'CASES.json',
        optimizer=Path(__file__).with_name('nb_state_guard.py'),
        predictor=Path(__file__).with_name('state_predictive_confirmation.py'),
        fixed_prevalence=Path(__file__).with_name('depth_state_predictive.py'))
    method = load('profiled_prevalence')if args.profile_prevalence else None
    if method is not None:
        inputs['profiled_optimizer']=Path(method.__file__)
        inputs['optimizer_correction_design']=args.task/'PROFILE_DESIGN.md'
    if args.case is not None:
        inputs['old_models'] = prior/f'case_{args.case:02d}_v1/DISCOVERY_MODELS_LOCKED.json'
        inputs['old_report'] = prior/f'case_{args.case:02d}_v1/report.json'
        inputs['numerics'] = args.numerical_report or args.task/('profile_numerics_v1/report.json'if method else'numerics_v1/report.json')
        inputs['trace'] = Path(cases[args.case]['trace'])
        prepared = load('state_predictive_confirmation').read_proposal(cases[args.case])
        for name in ('fold','confirmation','proposal_lock'):inputs[name]=Path(prepared[5][name])
        for order in (24,48,96):inputs[f'old_scores_{order}']=prior/f'case_{args.case:02d}_v1/per_cell_{order}.npy'
    if args.summary:
        paths=json.loads(args.case_paths.read_text())if args.case_paths else None
        if paths is not None:
            assert len(paths)==len(cases)
            inputs['case_paths']=args.case_paths
        for i in range(len(cases)):
            inputs[f'case{i}']=Path(paths[i])/'report.json'if paths else args.task/f'case_{i:02d}_v1/report.json'
    rec = start_run(args.out,inputs,dict(case=args.case,orders=[24,48,96]),
        dict(scope='Fixed exposed mechanism screen; no whole-specimen error estimate'),
        mode='discovery_likelihood_prevalence')
    try:
        if args.numerics:write_json(method.checks(shared=args.shared_prevalence)if method else checks(),args.out/'report.json')
        elif args.summary:
            results = [json.loads(inputs[f'case{i}'].read_text())for i in range(len(cases))]
            passed = all(r['resolved']and(r['score']>0 if c['kind']=='original_recovery'else r['score']<0)
                for c,r in zip(cases,results))
            write_json(dict(cases=results,screen_passed=passed,full_workflow_run=False,
                false_split_control_established=False,decision='Permit separately frozen full-workflow and fresh error/power tests'if passed else
                'Reject expansion, threshold tuning and production promotion'),args.out/'report.json')
        else:
            D,ED,rD,g,metadata,provenance,EC,rC = prepared
            lock = json.loads(inputs['proposal_lock'].read_text())
            strata = np.asarray(lock['parameters']['strata'],int)[metadata['selected']]
            q = metadata['q'][metadata['selected']]
            previous = json.loads(inputs['old_models'].read_text())
            assert previous['provenance']==provenance
            old_models = previous['models']
            models = []
            for old in old_models:
                print(json.dumps(dict(case=args.case,order=old['order'],status='fitting')),flush=True)
                models.append(method.fit(D,ED,rD,strata,q,old,shared=args.shared_prevalence)if method else fit(D,ED,rD,strata,q,old))
                print(json.dumps(dict(case=args.case,order=old['order'],status='discovery_fit_locked')),flush=True)
            write_json(dict(models=models,provenance=provenance,strata=strata.tolist(),
                proposal_q=q.tolist(),confirmation_counts_opened=False),args.out/'DISCOVERY_MODELS_LOCKED.json')
            with np.load(provenance['confirmation'])as z:C=z['counts'][metadata['selected']]
            with np.load(provenance['fold'])as z:
                full=sparse.csr_matrix((z['source_data'],z['source_indices'],z['source_indptr']),shape=z['source_shape'])
            total=full[metadata['cells']][:,metadata['genes']].toarray()[metadata['selected']]
            np.testing.assert_array_equal(total,D+C)
            rows,old_rows = [],[]
            for old,model in zip(old_models,models):
                control = load('depth_state_predictive').predict(D,C,ED,EC,rD,rC,strata,old)
                np.testing.assert_array_equal(control['per_cell'],np.load(inputs[f"old_scores_{old['order']}"]))
                old_rows.append(dict(order=old['order'],score=control['score']))
                result = (method.predict if method else predict)(D,C,ED,EC,rD,rC,strata,model)
                arrays = {key:result.pop(key)for key in ('per_cell','discovery_membership','joint_membership')}
                np.savez_compressed(args.out/f"prediction_{model['order']}.npz",**arrays)
                rows.append(dict(order=model['order'],**result))
            uncertainty = 2*max(abs(rows[j]['score']-rows[j-1]['score'])for j in (1,2))+sum(r['rounding_budget']for r in rows)
            assert old_rows[-1]['score']==json.loads(inputs['old_report'].read_text())['score']
            write_json(dict(case=cases[args.case],cells=len(D),genes=D.shape[1],
                quadrature=rows,old_scores=old_rows,old_prediction_arrays_exact=3,
                score=rows[-1]['score'],numerical_uncertainty_diagnostic=uncertainty,
                resolved=bool(abs(rows[-1]['score'])>uncertainty),
                convergence_is_not_formal_bound=True,fit_reads_confirmation=False,
                fold_count_conservation_exact=True,full_workflow_run=False,
                false_split_control_established=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU discovery-only free-prevalence NB shared-state prediction')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


if __name__=='__main__':main()
