"""Research shared-state NB prediction conditional on locked depth strata."""
import argparse
from fractions import Fraction
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
    spec = importlib.util.spec_from_file_location('depth_state_'+name, Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def component_parameters(x, pi, mixture):
    if not mixture:
        return [x[0], math.sqrt(x[-1])]
    delta = x[1]
    centre = np.logaddexp(np.log1p(-pi), np.log(pi)+delta)
    return [x[0]-centre, x[0]+delta-centre, np.log(pi)-np.log1p(-pi), math.sqrt(x[-1])]


def objective(x, y, e, r, strata, prevalence, order, mixture):
    k = len(prevalence)
    gradient, total = np.zeros(len(x)), 0.
    for level, pi in enumerate(prevalence):
        use = strata == level
        local = [x[level], x[-2], x[-1]] if mixture else [x[level], x[-1]]
        parameters = np.array(component_parameters(local, pi, mixture))
        parameters[-1] = x[-1]  # Corrected optimizer accepts variance, not SD.
        value, derivative = load('nb_state_guard').variance_objective(
            parameters, y[use], e[use], r[use], order, mixture)
        weight = use.sum()/len(y)
        total += weight*value
        if mixture:
            tilted = expit(np.log(pi)-np.log1p(-pi)+x[-2])
            gradient[level] = weight*(derivative[0]+derivative[1])
            gradient[-2] += weight*(-tilted*derivative[0]+(1-tilted)*derivative[1])
        else:
            gradient[level] = weight*derivative[0]
        gradient[-1] += weight*derivative[-1]
    return float(total), gradient


def fit(y, e, r, strata, q, order):
    y, e, r = np.broadcast_arrays(np.asarray(y, float), np.asarray(e, float), np.asarray(r, float))
    strata, q = np.asarray(strata), np.asarray(q, float)
    if (y.ndim != 2 or strata.shape != (len(y),) or strata.dtype.kind not in 'iu'
            or q.shape != strata.shape or not all(np.isfinite(v).all() for v in (y,e,r,q))
            or np.any(y < 0) or np.any(y != np.floor(y)) or np.any(e <= 0) or np.any(r <= 0)
            or np.any((q < 0)|(q > 1))):
        raise ValueError('Integer counts, supported NB means/sizes and locked strata/q required')
    levels, index = np.unique(strata, return_inverse=True)
    prevalence = np.array([q[index == s].mean() for s in range(len(levels))])
    if np.any((prevalence <= 0)|(prevalence >= 1)):
        raise ValueError('Both discovery mixture components need support in every retained stratum')
    z = np.log((y.sum(axis=1)+.5)/(e.sum(axis=1)+.5))
    alpha = np.array([np.median(z[index == s]) for s in range(len(levels))])
    variance = float(np.mean((z-alpha[index])**2))
    options = dict(ftol=np.finfo(float).eps, gtol=np.sqrt(np.finfo(float).eps), maxiter=1000, maxls=50)
    def optimize(starts, mixture):
        bounds = [(None, None)]*(len(levels)+int(mixture))+[(0., None)]
        fits = [minimize(objective, x, args=(y,e,r,index,prevalence,order,mixture),
                         jac=True, bounds=bounds, method='L-BFGS-B', options=options) for x in starts]
        good = [f for f in fits if f.success and np.isfinite(f.fun)]
        if not good:
            raise ArithmeticError('Depth-conditioned state optimization unresolved')
        return min(good, key=lambda f:f.fun), len(good)
    null, nulls = optimize([np.r_[alpha,v] for v in (0.,variance)], False)
    spread = float(np.diff(np.quantile(z, [.2,.8]))[0])
    alternative, alternatives = optimize([np.r_[null.x[:-1],d,v]
        for d in (0.,spread,-spread) for v in (0.,null.x[-1])], True)
    return dict(order=order, levels=levels.tolist(), prevalence=prevalence.tolist(),
        null=null.x.tolist(), alternative=alternative.x.tolist(),
        null_objective=float(null.fun), alternative_objective=float(alternative.fun),
        null_solutions=nulls, alternative_solutions=alternatives,
        optimization='Finite successful multistart L-BFGS-B; global optimum not certified')


def predict(D,C,ED,EC,rD,rC,strata,model):
    score, rounding = np.zeros(len(D)), 0.
    for s, pi in enumerate(model['prevalence']):
        use = strata == model['levels'][s]
        null = component_parameters([model['null'][s],model['null'][-1]], pi, False)
        alt = component_parameters([model['alternative'][s],*model['alternative'][-2:]], pi, True)
        result = load('state_predictive_confirmation').predictive(
            D[use],C[use],ED[use],EC[use],rD[use],rC[use],null,alt,model['order'])
        score[use] = result['per_cell']
        rounding += result['rounding_budget']
    return dict(score=float(math.fsum(map(float,score))), rounding_budget=rounding, per_cell=score)


def checks():
    D = np.array([[1,2],[3,0],[0,1],[2,3]])
    C = D[::-1].copy()
    E = np.tile([1.,2.], (4,1))
    rD, rC = np.full(D.shape,2.),np.full(D.shape,3.)
    strata, pi = np.array([0,0,1,1]),np.array([.25,.75])
    model = dict(levels=[0,1],prevalence=pi,order=96,null=[0.,0.,0.],
                 alternative=[math.log(1.25),math.log(1.75),math.log(2.),0.])
    actual = predict(D,C,E,E,rD,rC,strata,model)
    def pmf(y,e,r,t):
        p = Fraction(r)/(Fraction(r)+Fraction(e)*t)
        return math.comb(y+r-1,y)*p**r*(1-p)**y
    independent = []
    for i in range(4):
        weight = Fraction(float(pi[strata[i]]))
        mixture_d = mixture_j = Fraction(0)
        null_d = null_j = Fraction(1)
        for t,w in ((Fraction(1),1-weight),(Fraction(2),weight)):
            d = j = Fraction(1)
            for g in range(2):
                d *= pmf(int(D[i,g]),int(E[i,g]),2,t)
                j *= pmf(int(D[i,g]),int(E[i,g]),2,t)*pmf(int(C[i,g]),int(E[i,g]),3,t)
            mixture_d += w*d
            mixture_j += w*j
        for g in range(2):
            null_d *= pmf(int(D[i,g]),int(E[i,g]),2,Fraction(1))
            null_j *= pmf(int(D[i,g]),int(E[i,g]),2,Fraction(1))*pmf(int(C[i,g]),int(E[i,g]),3,Fraction(1))
        independent.append(math.log(float(mixture_j/mixture_d*null_d/null_j)))
    error = abs(actual['score']-math.fsum(independent))
    assert error <= actual['rounding_budget']
    nested = dict(model,alternative=[.1,-.2,0.,.04],null=[.1,-.2,.04])
    zero = predict(D,C,E,E,rD,rC,strata,nested)
    assert abs(zero['score']) <= zero['rounding_budget']
    flipped = dict(model, prevalence=1-pi, alternative=[math.log(1.25),math.log(1.75),-math.log(2.),0.])
    swapped = predict(D,C,E,E,rD,rC,strata,flipped)
    assert abs(swapped['score']-actual['score']) <= swapped['rounding_budget']+actual['rounding_budget']
    x = np.array([.1,-.2,.3,.04])
    value, gradient = objective(x,D,E,rD,strata,pi,96,True)
    logdensity, direct_gradient, quadrature_error = 0.,np.zeros(4),0.
    # Independently integrate full SciPy NB PMFs and the Gaussian heat-equation
    # derivative, rather than differentiating Gauss-Hermite nodes in the model.
    for i in range(4):
        s = strata[i];p = pi[s]
        means = component_parameters([x[s],x[-2],x[-1]],p,True)[:2]
        tilted = expit(np.log(p)-np.log1p(-p)+x[-2])
        mass, derivatives, errors = 0.,np.zeros(3),np.zeros(4)
        M = float(np.sum(D[i]+rD[i]))
        for component,weight in enumerate((1-p,p)):
            def integrand(z,kind):
                eta = means[component]+math.sqrt(x[-1])*z
                mu = E[i]*np.exp(eta)
                probability = mu/(rD[i]+mu)
                score = float(np.sum(D[i]-(D[i]+rD[i])*probability))
                hessian = -float(np.sum((D[i]+rD[i])*probability*(1-probability)))
                base = float(np.exp(nbinom.logpmf(D[i],rD[i],rD[i]/(rD[i]+mu)).sum())*norm.pdf(z))
                factor = (1.,score,score*(component-tilted),.5*(score*score+hessian))[kind]
                return base*factor
            for kind in range(4):
                part,err = quad(lambda z:integrand(z,kind),-12.,12.,
                    epsabs=256*np.finfo(float).eps,epsrel=256*np.finfo(float).eps,limit=200)
                bound = (1.,M,M,.5*(M*M+M/4))[kind]
                errors[kind] += weight*(err+2*norm.sf(12.)*bound)
                if kind == 0:mass += weight*part
                else:derivatives[kind-1] += weight*part
        assert mass > errors[0]
        logdensity += math.log(mass)
        direct_gradient[s] -= derivatives[0]/mass/4
        direct_gradient[-2:] -= derivatives[1:]/mass/4
        quadrature_error += errors[0]/(mass-errors[0])/4
        quadrature_error += float(np.sum(errors[1:]/(mass-errors[0])+
            abs(derivatives)*errors[0]/(mass*(mass-errors[0]))))/4
    baseline = float(nbinom.logpmf(D,rD,rD/(rD+E)).sum())
    reference_value = -(logdensity-baseline)/4
    operations = 512*96*D.size
    eps = np.finfo(float).eps
    rounding = operations*eps/(1-operations*eps)*max(1.,abs(baseline),abs(logdensity))
    budget = quadrature_error+rounding
    objective_error = abs(reference_value-value)
    gradient_error = float(np.max(abs(direct_gradient-gradient)))
    assert max(objective_error,gradient_error) <= budget
    original = objective(x,D,E,rD,strata,pi,96,True)
    perm = np.array([3,0,2,1])
    other = objective(x,D[perm],E[perm],rD[perm],strata[perm],pi,96,True)
    assert abs(original[0]-other[0]) <= budget and np.max(abs(original[1]-other[1])) <= budget
    exchanged = objective(x[[1,0,2,3]],D,E,rD,1-strata,pi[::-1],96,True)
    assert abs(original[0]-exchanged[0]) <= budget
    assert np.max(abs(original[1]-exchanged[1][[1,0,2,3]])) <= budget
    return dict(exact_fraction_score_error=error,exact_fraction_budget=actual['rounding_budget'],
        independent_full_PMF_quad_objective_error=objective_error,
        independent_heat_equation_quad_gradient_error=gradient_error,quad_and_rounding_budget=budget,
        nested_null=True,component_swap=True,cell_permutation=True,stratum_permutation=True,
        existing_predictive_integrals=load('state_predictive_confirmation').numerical_checks())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--case',type=int)
    parser.add_argument('--numerics',action='store_true')
    parser.add_argument('--summary',action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__),design=args.task/'DESIGN.md',cases=args.task/'CASES.json',
        predictive=Path(__file__).with_name('state_predictive_confirmation.py'),
        optimizer=Path(__file__).with_name('nb_state_guard.py'))
    cases = json.loads(inputs['cases'].read_text())
    prepared = None
    if args.case is not None:
        prepared = load('state_predictive_confirmation').read_proposal(cases[args.case])
        for name in ('fold','confirmation','proposal_lock'):
            inputs[name] = Path(prepared[5][name])
        inputs['trace'] = Path(cases[args.case]['trace'])
        inputs['numerical_checks'] = args.task/'numerics_v1/report.json'
    if args.summary:
        for i in range(len(cases)):
            inputs[f'case{i}'] = args.task/f'case_{i:02d}_v1/report.json'
    rec = start_run(args.out,inputs,dict(orders=[24,48,96],case=args.case),
        dict(scope='Exposed depth-conditioned predictive screen; biological n=0; no5% error claim'),
        mode='depth_state_predictive')
    try:
        if args.numerics:
            write_json(checks(),args.out/'report.json')
        elif args.summary:
            results = [json.loads(inputs[f'case{i}'].read_text()) for i in range(len(cases))]
            passed = all(r['resolved'] and (r['score']>0 if c['kind']=='original_recovery' else r['score']<0)
                         for c,r in zip(cases,results))
            write_json(dict(cases=results,screen_passed=passed,independent_evaluation=False,
                false_split_control_established=False,full_workflow_run=False,
                decision='Permit full original recovery and fresh whole-workflow error/power' if passed else
                'Reject expansion and threshold tuning; retain all original scientific criteria'),args.out/'report.json')
        else:
            D,ED,rD,g,metadata,provenance,EC,rC = prepared
            lock = json.loads(inputs['proposal_lock'].read_text())
            strata = np.asarray(lock['parameters']['strata'],int)[metadata['selected']]
            q = metadata['q'][metadata['selected']]
            models = [fit(D,ED,rD,strata,q,order) for order in (24,48,96)]
            write_json(dict(models=models,provenance=provenance,strata=strata.tolist(),
                q=q.tolist(),confirmation_counts_opened=False),args.out/'DISCOVERY_MODELS_LOCKED.json')
            with np.load(provenance['confirmation']) as z:
                C = z['counts'][metadata['selected']]
            with np.load(provenance['fold']) as z:
                full = sparse.csr_matrix((z['source_data'],z['source_indices'],z['source_indptr']),shape=z['source_shape'])
            total = full[metadata['cells']][:,metadata['genes']].toarray()[metadata['selected']]
            np.testing.assert_array_equal(total,D+C)
            rows = []
            for model in models:
                result = predict(D,C,ED,EC,rD,rC,strata,model)
                np.save(args.out/f"per_cell_{model['order']}.npy",result.pop('per_cell'))
                rows.append(dict(order=model['order'],**result))
            discrepancy = max(abs(rows[j]['score']-rows[j-1]['score']) for j in (1,2))
            uncertainty = 2*discrepancy+sum(r['rounding_budget'] for r in rows)
            write_json(dict(case=cases[args.case],cells=len(D),genes=D.shape[1],strata=len(models[-1]['levels']),
                quadrature=rows,score=rows[-1]['score'],numerical_uncertainty_diagnostic=uncertainty,
                resolved=bool(abs(rows[-1]['score'])>uncertainty),convergence_is_not_formal_bound=True,
                fit_reads_confirmation=False,fold_count_conservation_exact=True,
                full_workflow_run=False,false_split_control_established=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU stratum-conditioned shared-state NB integration and prediction')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


if __name__ == '__main__':
    main()
