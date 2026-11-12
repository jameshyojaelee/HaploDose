"""Population dispersion moments under soft and MAP RNA mixture membership."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import nbinom
from haplodose.provenance import start_run, finish_run


def tail_bounds(means, prior, N, r=5.):
    p = r/(r+means)
    tail = nbinom.sf(N, r, p)
    first = means*nbinom.sf(N-1, r+1, p)
    factorial_second = means**2*(1+1/r)*nbinom.sf(N-2, r+2, p)
    largest = float(max(means))
    # Bound the absolute residual for either selected mean; conditional count law
    # remains the original component law when forming its factorial moments.
    numerator = float(prior@(factorial_second+(2*largest+2)*first+(largest**2+largest)*tail))
    denominator = float(largest**2*(prior@tail))
    return numerator, denominator


def calculate(means, prior, N, subtraction='observed'):
    y = np.arange(N+1)
    logjoint = nbinom.logpmf(y[:, None], 5, 5/(5+means))[...]+np.log(prior)
    logmix = logsumexp(logjoint, axis=1)
    probability = np.exp(logmix)
    responsibility = np.exp(logjoint-logmix[:, None])
    assignment = logjoint.argmax(axis=1)
    moment = (y[:, None]-means)**2-(y[:, None] if subtraction == 'observed' else means)
    rows = []
    for method, weight in [('soft', responsibility), ('MAP', np.eye(2)[assignment])]:
        numerator = float(np.sum(probability[:, None]*weight*moment))
        denominator = float(np.sum(probability[:, None]*weight*means**2))
        estimate = numerator/denominator
        num_tail, den_tail = tail_bounds(means, prior, N)
        bound = (num_tail+abs(estimate)*den_tail)/denominator
        rows.append(dict(method=method, numerator=numerator, denominator=denominator,
            estimated_phi=estimate, ratio_tail_bound=bound))
    return rows, assignment


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--subtraction', choices=['observed', 'mean'], default='observed')
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    design = a.task/('DESIGN.md' if a.subtraction == 'observed' else 'MATCHED_DESIGN.md')
    run = start_run(a.out, dict(script=Path(__file__), design=design), dict(subtraction=a.subtraction),
        dict(scope='Known-parameter RNA count-mixture identity; biological n=0'), mode='dispersion_membership')
    try:
        rows, cases, arrays = [], [], {}
        for high in (6, 12, 30):
            for pi in (.1, .5, .9):
                means, prior = np.array([4., high]), np.array([pi, 1-pi])
                N = 32
                while True:
                    coarse, _ = calculate(means, prior, N, a.subtraction)
                    if max(row['ratio_tail_bound'] for row in coarse) < 1e-12:
                        break
                    N *= 2
                    if N > 16384:
                        raise ArithmeticError('Unresolved NB moment tail')
                fine, assignment = calculate(means, prior, 2*N, a.subtraction)
                key = f'm{high}_p{int(pi*10)}'
                arrays[key+'__assignment'] = assignment
                for c, f in zip(coarse, fine):
                    operations = 256*(2*N+1)
                    unit = np.finfo(float).eps
                    budget = operations*unit/(1-operations*unit)*max(1., abs(f['estimated_phi']))
                    difference = abs(c['estimated_phi']-f['estimated_phi'])
                    assert difference <= c['ratio_tail_bound']+f['ratio_tail_bound']+budget
                    if f['method'] == 'soft':
                        assert abs(f['estimated_phi']-.2) <= f['ratio_tail_bound']+budget
                    rows.append(dict(case=key, mean_low=4, mean_high=high, prior_low=pi,
                        true_phi=.2, N=2*N, coarse_N=N, numerical_budget=budget,
                        refinement_change=difference, **f))
                cases.append(dict(case=key, means=means.tolist(), prior=prior.tolist(), N=2*N))
        pd.DataFrame(rows).to_csv(a.out/'moments.tsv', sep='\t', index=False)
        np.savez_compressed(a.out/'assignments.npz', **arrays)
        (a.out/'report.json').write_text(json.dumps(dict(cases=cases, biological_specimens=0, subtraction=a.subtraction,
            production_changed=False, scientific_criteria_changed=False,
            estimand='Known-mean common NB dispersion; no normalization or clipping'), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, run, error=error)
        raise
    finish_run(a.out, run, backend='CPU SciPy NB mixture moments and factorial-moment tail bounds')


if __name__ == '__main__':
    main()
