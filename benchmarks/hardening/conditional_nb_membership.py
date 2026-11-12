"""Conditional NB normalizers on saved RNA states; no inference replacement."""
import argparse
from fractions import Fraction
from functools import lru_cache
import importlib.util
import itertools
import json
import math
from pathlib import Path

import jax
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import logsumexp, roots_legendre
from scipy.stats import nbinom

from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


@lru_cache(None)
def load(name):
    path = Path(__file__).with_name(name+'.py')
    spec = importlib.util.spec_from_file_location('conditional_nb_'+name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def support(means, sizes, count):
    mu, r = np.broadcast_arrays(np.asarray(means, float), np.asarray(sizes, float))
    if (mu.ndim != 1 or not len(mu) or not np.isfinite(mu).all() or np.any(mu <= 0)
            or not np.isfinite(r).all() or np.any(r <= 0) or not np.isfinite(count)
            or count < 0 or count != int(count)):
        raise ValueError('Positive finite NB means/sizes and an integer total required')
    return mu, r


def log_characteristic(t, mu, r, count):
    # expm1 avoids loss in 1-exp(it) near zero. The -it*count phase supplies
    # the Fourier coefficient, not the probability of an unconditioned mean.
    t = np.asarray(t)
    return -np.sum(r[:, None]*np.log1p(-(mu/r)[:, None]*np.expm1(1j*t)), axis=0)-1j*count*t


def integration_setup(mu, r, count):
    variance = float(np.sum(mu+mu*mu/r))
    sd = math.sqrt(variance)
    epsilon = np.finfo(float).eps
    operations = 128*len(mu)+256
    gamma = operations*epsilon/(1-operations*epsilon)
    # |characteristic(t)| decreases on[0,pi]. The omitted real Fourier
    # integral is bounded by (pi-cut)*|characteristic(cut)|/pi.
    target = gamma/max(1., sd)
    lo, hi = 0., math.pi
    for _ in range(np.finfo(float).nmant+2):
        cut = (lo+hi)/2
        logmod = -.5*np.sum(r*np.log1p(4*(mu/r)*(1+mu/r)*math.sin(cut/2)**2))
        tail = (math.pi-cut)*math.exp(float(logmod))/math.pi
        if tail > target:
            lo = cut
        else:
            hi = cut
    cut = hi
    logmod = -.5*np.sum(r*np.log1p(4*(mu/r)*(1+mu/r)*math.sin(cut/2)**2))
    tail = (math.pi-cut)*math.exp(float(logmod))/math.pi
    return cut, tail, gamma, sd


def total_logpmf(count, means, sizes, orders=(32, 64, 128, 256, 512)):
    mu, r = support(means, sizes, count)
    cut, tail, gamma, sd = integration_setup(mu, r, count)
    studies = []
    previous = None
    for order in orders:
        nodes, weights = roots_legendre(order)
        t = .5*cut*(nodes+1)
        values = log_characteristic(t, mu, r, count)
        probability = .5*cut/math.pi*float(weights@np.exp(values).real)
        # Phase and log-PGF cancellation are explicitly represented. This is
        # an arithmetic estimate plus an analytic tail bound, not a quadrature
        # error enclosure or a statistical uncertainty interval.
        magnitude = 1+abs(count)*cut+float(np.max(abs(values)))
        arithmetic = gamma*magnitude*cut/math.pi
        if probability <= tail+arithmetic:
            previous = probability
            studies.append(dict(order=order, probability=probability, unresolved=True))
            continue
        difference = abs(probability-previous) if previous is not None else math.inf
        diagnostic = tail+arithmetic+2*difference
        studies.append(dict(order=order, probability=probability, tail_bound=tail,
                            arithmetic_estimate=arithmetic, consecutive_difference=difference))
        if previous is not None and difference <= tail+arithmetic:
            return dict(log_probability=math.log(probability), probability=probability,
                log_error_diagnostic=diagnostic/(probability-diagnostic) if probability > diagnostic else math.inf,
                order=order, cut=cut, sd=sd, studies=studies,
                uncertainty_scope='analytic tail plus arithmetic and convergence diagnostic; not coverage')
        previous = probability
    raise ArithmeticError('Fourier coefficient did not resolve at prespecified orders')


def adaptive_reference(count, mu, r):
    cut, tail, gamma, sd = integration_setup(mu, r, count)
    def integrand(u):
        return float(np.exp(log_characteristic(np.array([u/sd]), mu, r, count))[0].real)/math.pi/sd
    value, error = quad(integrand, 0., cut*sd, epsabs=gamma/max(1., sd),
                        epsrel=gamma, limit=200)
    if value <= error+tail:
        raise ArithmeticError('Independent adaptive integral unresolved')
    return math.log(value), (error+tail)/(value-error-tail)


def numerical_checks():
    recurrence = load('gene_bin_reference')
    rows = []
    for count, mu, r in [(0, [1., 4.], [1., 2.]), (5, [1., 4.], [1., 2.]),
                          (20, [.2, 3., 20.], [2., 10., .5]),
                          (120, [6., 12., 42.], [1., 2., 7.])]:
        mu, r = np.asarray(mu), np.asarray(r)
        actual = total_logpmf(count, mu, r)
        reference, budget = recurrence.summed_logpmf(np.array([count]), mu[None], r)
        independent, independent_error = adaptive_reference(count, mu, r)
        allowance = actual['log_error_diagnostic']+float(budget[0])+independent_error
        assert abs(actual['log_probability']-reference[0]) <= allowance
        assert abs(actual['log_probability']-independent) <= allowance
        permuted = total_logpmf(count, mu[::-1], r[::-1])
        assert abs(permuted['log_probability']-actual['log_probability']) <= 2*allowance
        rows.append(dict(count=count, fourier_log_probability=actual['log_probability'],
            recurrence_discrepancy=abs(actual['log_probability']-float(reference[0])),
            adaptive_discrepancy=abs(actual['log_probability']-independent), numerical_allowance=allowance))
    # Exact enumeration of all allocations conditional on total4, with two
    # distinct NB probabilities. The denominator is independently rational.
    total = 4
    probabilities = [Fraction(1, 2)**(a+1)*Fraction(1, 3)*Fraction(2, 3)**(total-a)
                     for a in range(total+1)]
    exact = sum(probabilities)
    normalizer = total_logpmf(total, [1., 2.], [1., 1.])
    assert abs(normalizer['probability']-float(exact)) <= normalizer['probability']*normalizer['log_error_diagnostic']
    conditional = [p/exact for p in probabilities]
    assert sum(conditional) == 1
    for a, expected in enumerate(conditional):
        raw = nbinom.logpmf(a, 1., .5)+nbinom.logpmf(total-a, 1., 1/3)
        actual = math.exp(float(raw)-normalizer['log_probability'])
        assert abs(actual-float(expected)) <= normalizer['log_error_diagnostic']
    return dict(passed=True, independent_recurrence_and_adaptive_checks=rows,
                exact_conditional_sum='1', conditional_allocations=list(map(str, conditional)),
                permutation_checked=True, scope='Numerical agreement, not genomic or coverage validation')


def matching(labels, truth):
    normal = truth == 'normal'
    components = sorted(set(labels)-{0})
    best = int(np.sum((labels == 0) & normal))
    # One-to-one mapping of tumour components, with extra groups receiving no
    # duplicate credit. All cells, including reference-like calls, remain.
    scores = []
    for assignments in itertools.permutations(components, min(2, len(components))):
        for names in itertools.permutations(('clone1', 'clone2'), len(assignments)):
            scores.append(sum(int(np.sum((labels == k) & (truth == name)))
                              for k, name in zip(assignments, names)))
    return best+max(scores, default=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if jax.default_backend() != 'cpu':
        raise RuntimeError('CPU compute required')
    args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    source = root/'dev/production_profile_trace_20261109/trace_v1'
    source_report = source/'report.json'
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DESIGN.md',
        source=source_report, source_manifest=source/'run_manifest.json',
        reference=Path(__file__).with_name('gene_bin_reference.py')),
        dict(seeds=[20261031, 20261107, 20261108, 20261109, 20261110],
             selection='first S1 state at maximum fitted component count per original case'),
        dict(scope='Frozen exposed RNA states, no inference rerun or genomic validation'),
        mode='conditional_nb_membership')
    try:
        assert json.loads((source/'run_manifest.json').read_text())['status'] == 'complete'
        write_json(numerical_checks(), args.out/'numerics.json')
        rows = json.loads(source_report.read_text())['records']
        selected = []
        for seed in (20261031, 20261107, 20261108, 20261109, 20261110):
            candidates = []
            for row in rows:
                if row['seed'] != seed:
                    continue
                with np.load(row['model']) as m, np.load(row['arrays']) as z:
                    if len(m['L']) == 600:
                        candidates.append((len(z['before']), row))
            maximum = max(c[0] for c in candidates)
            selected.append(next(row for k, row in candidates if k == maximum))
        write_json(dict(records=selected, selection_independent_of_truth=True), args.out/'STATES_LOCKED.json')
        outputs = []
        for index, row in enumerate(selected):
            with np.load(row['model']) as m, np.load(row['arrays']) as z:
                Y, L, x = m['Y'], m['L'], m['x']
                assert len(Y) == 600 and np.array_equal(Y.sum(axis=1), L)
                R, doses = z['R'].astype(float), z['before']
                R /= R.sum(axis=1, keepdims=True)
                offset = np.log(z['base'])+z['offset']
                nuisance = x[:, None]*z['gamma'][None]+z['U']@z['W']
                shapes = 1/z['phi']
                logs, errors = np.empty_like(R), np.empty_like(R)
                checks = []
                for k, dose in enumerate(doses):
                    lr = offset[None]+dose[None]+nuisance
                    expected = L[:, None]*np.exp(lr-logsumexp(lr, axis=1, keepdims=True))
                    for cell in range(len(Y)):
                        fitted = total_logpmf(int(L[cell]), expected[cell], shapes)
                        logs[cell, k], errors[cell, k] = fitted['log_probability'], fitted['log_error_diagnostic']
                        if cell in (0, len(Y)-1):
                            independent, independent_error = adaptive_reference(int(L[cell]), expected[cell], shapes)
                            assert abs(independent-logs[cell, k]) <= independent_error+errors[cell, k]
                            checks.append(dict(cell=cell, component=k, independent_log_probability=independent,
                                discrepancy=abs(independent-logs[cell, k]), budget=independent_error+errors[cell, k]))
                logR = np.full_like(R, -np.inf)
                np.log(R, out=logR, where=R > 0)
                corrected_log = logR-logs
                corrected = np.exp(corrected_log-logsumexp(corrected_log, axis=1, keepdims=True))
                labels = corrected.argmax(axis=1)
                sorted_logits = np.sort(corrected_log, axis=1)
                margins = sorted_logits[:, -1]-sorted_logits[:, -2]
                unresolved = margins <= 2*np.max(errors, axis=1)
                path = args.out/f'case_{index}.npz'
                np.savez_compressed(path, source_R=R, conditional_R=corrected, log_normalizers=logs,
                    log_error_diagnostics=errors, rows=m['rows'], numerical_decision_unresolved=unresolved)
                outputs.append(dict(seed=row['seed'], source=row, prediction=file_record(path),
                    components=R.shape[1], map_changes=int(np.sum(labels != R.argmax(axis=1))),
                    maximum_responsibility_change=float(np.max(abs(corrected-R))),
                    maximum_pairwise_log_odds_correction=float(np.max(np.ptp(logs, axis=1))),
                    stored_zero_responsibilities=int(np.sum(R == 0)),
                    unresolved_decisions=int(unresolved.sum()), independent_checks=checks))
                print(json.dumps({k:v for k,v in outputs[-1].items() if k not in ('source', 'independent_checks')}), flush=True)
        write_json(dict(cases=outputs, truth_opened=False), args.out/'RNA_LOCKED.json')
        for output in outputs:
            seed = output['seed']
            truth_path = root/f'dev/recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv'
            truth = pd.read_csv(truth_path, sep='\t')
            with np.load(output['prediction']['path']) as z:
                truth = truth.iloc[z['rows']]
                assert truth.specimen.eq('S1').all()
                before = matching(z['source_R'].argmax(axis=1), truth.truth_clone.to_numpy())
                after = matching(z['conditional_R'].argmax(axis=1), truth.truth_clone.to_numpy())
            output.update(truth=file_record(truth_path), correct_before=before, correct_after=after)
        advance = (any(o['correct_after'] > o['correct_before'] for o in outputs)
                   and all(o['correct_after'] >= o['correct_before'] and o['unresolved_decisions'] == 0 for o in outputs))
        write_json(dict(cases=outputs, continuation_passed=advance, inference_rerun=False,
            original_recovery_criteria_changed=False, original_recovery_improved=False,
            biological_specimens=0, production_changed=False, promotion=False,
            interpretation='Conditional reweighting of fixed saved states; not complete-workflow accuracy'),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU Fourier coefficients and independent quadrature; frozen states')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
