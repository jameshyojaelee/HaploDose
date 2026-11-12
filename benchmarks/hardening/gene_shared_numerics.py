"""Independent joint-state/gene-count checks for research confirmation."""
import argparse
import importlib.util
import itertools
from pathlib import Path

import numpy as np
from scipy.special import logsumexp
from scipy.stats import nbinom, betabinom
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('gene_shared_confirmation.py')
    spec = importlib.util.spec_from_file_location('checked_gene_model', path)
    method = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(method)
    engine = method.load('shared_membership_hmm')
    run = start_run(args.out, dict(script=Path(__file__), method=path,
        engine=Path(engine.__file__), design=args.task/'DESIGN.md'), {},
        dict(scope='Conditional exact numerical model, no biological validation'), mode='gene_shared_numerics')
    try:
        O = np.array([[1., 3., 2., 4.], [4., 2., 1., 5.], [2., 7., 3., 1.]])
        E = np.array([[1.5, 2., 3., 2.], [2., 3., 2., 4.], [4., 1., 2., 3.]])
        r = np.array([1.2, 2.3, 3.4, 4.5])
        bins, q = np.array([0, 0, 1, 1]), np.array([.15, .65, .8])
        effects, prior, change, pi = np.array([-.3, 0., .3]), np.array([.1, .8, .1]), .02, np.array([.4, .6])
        results = []
        for left, right in itertools.product((False, True), repeat=2):
            options = dict(left_neutral=left, right_neutral=right, mixture_proportion=pi)
            model = method.gene_model(O, E, r, q, bins, effects, prior, change, **options)
            fit = engine.enumerate_states(model)
            paths = list(itertools.product(range(3), repeat=2))
            def transition(a, b):
                return 1-change if a == b else change*prior[b]/(1-prior[a])
            masses = []
            for s in paths:
                weight = transition(1, s[0]) if left else prior[s[0]]
                weight *= transition(s[0], s[1])
                if right:
                    weight *= transition(s[-1], 1)
                masses.append(weight)
            norm, base, logs, assignments = sum(masses), nbinom.logpmf(O, r, r/(r+E)).sum(), [], []
            for s, weight in zip(paths, masses):
                for z in itertools.product((0, 1), repeat=3):
                    mu = np.empty_like(E)
                    for c in range(3):
                        for g in range(4):
                            delta = effects[s[bins[g]]]
                            mu[c, g] = E[c, g]*np.exp(delta*z[c])/(1-pi[bins[g]]+pi[bins[g]]*np.exp(delta))
                    prob = np.prod([q[c] if z[c] else 1-q[c] for c in range(3)])
                    logs.append(float(nbinom.logpmf(O, r, r/(r+mu)).sum()-base+np.log(weight/norm)+np.log(prob)))
                    assignments.append((s, z))
            evidence = float(logsumexp(logs))
            post, membership = np.zeros((2, 3)), np.zeros(3)
            for w, (s, z) in zip(np.exp(np.asarray(logs)-evidence), assignments):
                for j in range(2):
                    post[j, s[j]] += w
                membership += w*np.asarray(z)
            operations = 64*(len(logs)+O.size)
            unit = np.finfo(float).eps
            budget = operations*unit/(1-operations*unit)*max(1., abs(base), max(map(abs, logs)))
            errors = [abs(fit['log_evidence']-evidence), float(np.max(abs(fit['posterior']-post))),
                      float(np.max(abs(fit['membership']-membership)))]
            assert max(errors) <= budget
            perm = np.array([2, 0, 3, 1])
            other = engine.enumerate_states(method.gene_model(O[:, perm], E[:, perm], r[perm], q,
                bins[perm], effects, prior, change, **options))
            np.testing.assert_allclose(other['posterior'], post, rtol=0., atol=budget)
            cell_perm = np.array([2, 0, 1])
            reordered = engine.enumerate_states(method.gene_model(O[cell_perm], E[cell_perm], r, q[cell_perm],
                bins, effects, prior, change, **options))
            np.testing.assert_allclose(reordered['posterior'], post, rtol=0., atol=budget)
            np.testing.assert_allclose(reordered['membership'], membership[cell_perm], rtol=0., atol=budget)
            choose = np.array([2, 0])
            subset = method.gene_model(O[choose], E[choose], r, q[choose], bins, effects, prior, change, **options)
            np.testing.assert_array_equal(subset['ll'], model['ll'][:, choose])
            results.append(dict(left=left, right=right, errors=errors, operation_rounding_budget=budget))
        # Retain independently checked algebraic counterexamples from the proposal
        # experiment: equal summed counts do not imply equal gene NB evidence.
        k = np.arange(7)
        mass = betabinom.pmf(k, 6, 2/3, 4/3)
        # Independent exact rising-factorial probability and moments.
        from fractions import Fraction
        from math import comb
        def rising(x, n):
            result = Fraction(1)
            for j in range(n):
                result *= x+j
            return result
        exact_mass = [comb(6, j)*rising(Fraction(2,3),j)*rising(Fraction(4,3),6-j)/rising(Fraction(2),6) for j in range(7)]
        assert sum(exact_mass) == 1
        mean = sum(j*p for j,p in enumerate(exact_mass))
        variance = sum((j-mean)**2*p for j,p in enumerate(exact_mass))
        assert mean == 2 and variance == Fraction(32,9)
        # Four accumulation/conversion operations per factor across the two
        # degree-six products and denominator, for all seven probabilities.
        probability_operations = 4*len(k)*(2*int(k.max())+1)
        probability_budget = probability_operations*np.finfo(float).eps/(1-probability_operations*np.finfo(float).eps)
        probability_error = float(np.max(abs(mass-np.asarray(list(map(float, exact_mass))))))
        assert probability_error <= probability_budget
        fixture = np.arange(12).reshape(3,4)
        weights = np.array([0., 1/3, 1.])
        thinned = method.thin_matching(fixture, 2., weights, 2026110913)
        np.testing.assert_array_equal(thinned, method.thin_matching(fixture, 2., weights, 2026110913))
        np.testing.assert_array_equal(thinned[0], np.zeros(4))
        np.testing.assert_array_equal(thinned[2], fixture[2])
        assert np.all(thinned <= fixture)
        write_json(dict(all_checks_passed=True, independent_joint_enumeration=results,
            exact_beta_binomial_mass=True, conditional_mean=str(mean), conditional_variance=str(variance),
            beta_binomial_probability_discrepancy=probability_error,
            beta_binomial_probability_rounding_budget=probability_budget,
            thinning_endpoints_count_bounds_determinism=True,
            gene_permutation=True, fixed_rate_cell_subsetting=True,
            same_sum_information_reference='dev/gene_event_proposals_20261109/numerics_v1/report.json'), args.out/'report.json')
        finish_run(args.out, run, backend='CPU explicit SciPy joint state/cell enumeration')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
