"""Independent proposal-only optimization diagnostic on three locked events.

This is a finite multistart search, not a certificate of a global optimum.
Confirmation and simulated labels cannot select a starting point or solution.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, minimize
from scipy.special import expit, logsumexp
from scipy.stats import nbinom

from haplodose import clones, profile
from haplodose.provenance import start_run, finish_run


def canonical(parameters):
    x = np.asarray(parameters).copy()
    if x[0] > x[1]:
        x = x[[1, 0, 2]]
        x[2] *= -1
    return x


def reference(counts, exposure, size):
    """Full SciPy NB probabilities and separately derived component scores."""
    def objective(x):
        means = exposure[None] * np.exp(np.asarray(x[:2]))[:, None, None]
        components = nbinom.logpmf(counts[None], size[None], size[None] / (size[None] + means)).sum(axis=2)
        weighted = components - np.logaddexp(0., np.array([x[2], -x[2]]))[:, None]
        total = logsumexp(weighted, axis=0)
        q = np.exp(weighted[1] - total)
        scores = (size[None] * (counts[None] - means) / (size[None] + means)).sum(axis=2)
        gradient = -np.array([np.dot(1-q, scores[0]), np.dot(q, scores[1]),
                              np.sum(q-expit(x[2]))]) / len(counts)
        return -float(total.mean()), gradient, q
    return objective


def scalar_loglik(x, counts, exposure, size):
    terms = []
    for observed, expected, shape in zip(counts, exposure, size):
        components = []
        for eta in x[:2]:
            values = []
            for o, e, r in zip(observed, expected, shape):
                mean = e * math.exp(eta)
                values.append(math.lgamma(o+r)-math.lgamma(r)-math.lgamma(o+1)
                              +r*math.log(r/(r+mean))+o*math.log(mean/(r+mean)))
            components.append(math.fsum(values))
        a = components[0]-float(np.logaddexp(0., x[2]))
        b = components[1]-float(np.logaddexp(0., -x[2]))
        terms.append(max(a,b)+math.log1p(math.exp(-abs(a-b))))
    return math.fsum(terms)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--trace-root', type=Path, required=True)
    parser.add_argument('--panel-root', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    model_path = Path(__file__).with_name('shared_membership_hmm.py')
    spec = importlib.util.spec_from_file_location('exact_model', model_path)
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    inputs = {'runner': __file__, 'exact_model': model_path}
    events = []
    for seed in (20261108,20261109,20261110):
        directory = args.trace_root / f'trace_{seed}'
        for name in ('trace.json','fit/bins.tsv','fit/specimens.tsv'):
            inputs[f'{seed}:{name}'] = directory / name
        inputs[f'{seed}:truth'] = args.panel_root / f'seed_{seed}/data/cells.tsv'
        for event in json.loads((directory/'trace.json').read_text()):
            if event['stage'] == 'event_fit':
                path = directory/event['arrays']
                inputs[f'{seed}:event{event["index"]}'] = path
                events.append((seed,event,path,directory))
    if len(events) != 3:
        raise ValueError('Locked event population changed')
    rec = start_run(args.out, inputs, dict(starts=10, solvers=['L-BFGS-B','Nelder-Mead']),
                    dict(selection='Proposal counts only; numerical ties retain saved solution',
                         scope='Postfit conditional diagnostic, no calibration or global optimum claim'),
                    mode='membership_optimizer_probe')
    try:
        rows = []
        eps = np.finfo(float).eps
        pp = clones._thin_diff_params(clones.CloneParams(split_mode='thin').profile)
        bound = math.log(1-1e-3)-math.log(1e-3)
        for seed,event,path,directory in events:
            with np.load(path) as z:
                arrays = {key:z[key].copy() for key in z.files}
            O,E = arrays['O'],arrays['E']
            r = np.broadcast_to(arrays['r'],O.shape)
            saved = canonical(arrays['parameters'])
            objective = reference(O,E,r)
            # Null score is monotone and this bracket is expanded, not clipped.
            def score(eta):
                mu = E*np.exp(eta)
                return float(np.sum(r*(O-mu)/(r+mu)))
            lo,hi = -1.,1.
            while score(lo)<=0: lo *= 2
            while score(hi)>=0: hi *= 2
            null = brentq(score,lo,hi,xtol=eps,rtol=4*eps)
            beta = float(pd.read_csv(directory/'fit/specimens.tsv',sep='\t').set_index('specimen').loc['S1','beta'])
            step = beta*np.log(1.5)
            starts = [saved] + [np.array([null-d*step,null+d*step,math.log(p/(1-p))])
                               for d in (.25,1.,4.) for p in (.1,.5,.9)]
            saved_value = objective(saved)[0]
            # Full PMF subtraction has a larger rounding scale than the package's
            # null-relative objective. Preserve ties within this explicit budget.
            operations = O.size+32*len(O)+64
            gamma = operations*eps/(1-operations*eps)
            budget = gamma*max(1.,abs(saved_value)*len(O),float(O.sum()+r.sum()))
            fits = []
            best_x,best_value = saved,saved_value
            for solver in ('L-BFGS-B','Nelder-Mead'):
                for index,start in enumerate(starts):
                    options = (dict(maxiter=3000,ftol=eps,gtol=eps**.75,maxls=50)
                               if solver=='L-BFGS-B' else dict(maxiter=6000,xatol=eps**.5,fatol=eps))
                    fit = minimize((lambda x: objective(x)[:2]) if solver=='L-BFGS-B' else
                                   (lambda x: objective(x)[0]),start,method=solver,
                                   jac=True if solver=='L-BFGS-B' else False,
                                   bounds=[(None,None),(None,None),(-bound,bound)],options=options)
                    x = canonical(fit.x)
                    value,gradient,q = objective(x)
                    fits.append(dict(solver=solver,start=index,success=bool(fit.success),message=str(fit.message),
                                     parameters=x.tolist(),proposal_gain_vs_saved=(saved_value-value)*len(O),
                                     gradient_per_cell=gradient.tolist(),iterations=int(fit.nit)))
                    if fit.success and np.isfinite(value) and (best_value-value)*len(O)>budget:
                        best_x,best_value = x,value
            references = []
            for label,x in (('saved',saved),('selected',best_x)):
                value,gradient,q = objective(x)
                scalar = scalar_loglik(x,O,E,r)
                discrepancy = abs(scalar+value*len(O))
                if discrepancy>budget:
                    raise ArithmeticError('Scalar NB reference disagrees beyond rounding budget')
                null_value = objective(np.array([null,null,0.]))[0]
                fast_value,fast_gradient,fast_q = clones._event_mixture_objective(x,O,E,r,null)
                if abs(fast_value-(value-null_value)*len(O))>budget:
                    raise ArithmeticError('Package relative objective disagrees with full PMF reference')
                if np.max(np.abs(fast_q-q))>gamma*max(1.,abs(value)):
                    raise ArithmeticError('Package membership disagrees with full PMF reference')
                # Report a step-size convergence study, not a single convenient epsilon.
                derivatives = []
                for power in (8,12,16,20):
                    h = 2.**(-power)
                    fd = []
                    for axis in range(3):
                        delta = np.zeros(3);delta[axis] = h
                        fd.append((objective(x+delta)[0]-objective(x-delta)[0])/(2*h))
                    derivatives.append(dict(step=h,maximum_score_difference=float(np.max(np.abs(fd-gradient)))))
                references.append(dict(label=label,scalar_loglik_discrepancy=discrepancy,
                    relative_objective_discrepancy=abs(fast_value-(value-null_value)*len(O)),
                    analytic_gradient_discrepancy=float(np.max(np.abs(fast_gradient/len(O)-gradient))),
                    finite_difference_convergence=derivatives))
            # Proposal-only selection is complete before confirmation or truth is read.
            cols = arrays['cols']
            chrom = pd.read_csv(directory/'fit/bins.tsv',sep='\t').chromosome.to_numpy()
            confirmation = []
            for label,x in (('saved',saved),('selected',best_x)):
                q = objective(x)[2]
                fitted = model.build(arrays['Ytest'],E,r,q,profile.DIFF_LEVELS*step,
                    profile.DIFF_PRIOR,pp.p_change,left_neutral=bool(cols[0]>0 and chrom[cols[0]-1]==chrom[cols[0]]),
                    right_neutral=bool(cols[-1]+1<len(chrom) and chrom[cols[-1]+1]==chrom[cols[-1]]))
                result = model.enumerate_states(fitted)
                cores = [(int(cols[a]),int(cols[b]),int(s)) for a,b,s in
                         profile.credible_difference_runs(result['path'],result['posterior'],
                                                         np.r_[True,np.zeros(len(cols)-1,bool)],pp)]
                confirmation.append(dict(label=label,cores=cores,posterior=result['posterior'].tolist(),
                                         log_evidence=result['log_evidence']))
            truth = pd.read_csv(inputs[f'{seed}:truth'],sep='\t').iloc[arrays['rows'][arrays['cells']]]
            if not truth.specimen.eq('S1').all(): raise ValueError('Specimen mismatch')
            y = truth.truth_clone.eq('clone2').to_numpy(float)
            q_before,q_after = objective(saved)[2],objective(best_x)[2]
            row = dict(seed=seed,trace_index=event['index'],fits=fits,references=references,
                numerical_gain_budget=budget,selected_parameters=best_x.tolist(),
                selected_gain_vs_saved=(saved_value-best_value)*len(O),
                saved_retained=bool(np.array_equal(best_x,saved)),
                max_membership_change=float(np.max(np.abs(q_before-q_after))),
                postfit_brier_before=float(np.mean((q_before-y)**2)),postfit_brier_after=float(np.mean((q_after-y)**2)),
                confirmation=confirmation)
            rows.append(row)
            print(json.dumps({k:v for k,v in row.items() if k not in ('fits','references','confirmation')}),flush=True)
        base = json.loads((args.task/'baseline/source_manifest.json').read_text())
        unchanged = 0
        for relative,digest in base.items():
            if hashlib.sha256((args.task/'baseline'/relative).read_bytes()).hexdigest()!=digest:
                raise ValueError('Baseline changed')
            if relative.startswith(('python/','tests/','crates/')):
                if hashlib.sha256((args.task.parents[1]/relative).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Production source changed')
                unchanged += 1
        (args.out/'summary.json').write_text(json.dumps(dict(cases=rows,baseline_files_preserved=len(base),
            core_test_native_unchanged=unchanged,scope='Finite proposal-only multistart search; not a global optimum certificate; no inference changes'),indent=2))
        finish_run(args.out,rec,backend='cpu')
    except Exception as error:
        finish_run(args.out,rec,error=error)
        raise


if __name__ == '__main__':
    main()
