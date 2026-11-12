"""Positive NB-sum coefficients; numerical evidence, not a deployable classifier."""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats
from haplodose.provenance import start_run, finish_run


def summed_logpmf(counts, means, sizes):
    counts = np.asarray(counts)
    means, sizes = np.asarray(means, np.longdouble), np.asarray(sizes, np.longdouble)
    if counts.ndim != 1 or not len(counts) or means.ndim != 2 or means.shape[0] != len(counts):
        raise ValueError('Aligned count vector and cell-by-gene means required')
    if sizes.ndim > 1 or sizes.size not in (1, means.shape[1]):
        raise ValueError('A scalar or gene-aligned size vector required')
    if not np.isfinite(means).all() or not np.isfinite(sizes).all() or np.any(means <= 0) or np.any(sizes <= 0):
        raise ValueError('Finite positive means and sizes required')
    if np.any(counts < 0) or np.any(counts != np.floor(counts)):
        raise ValueError('Nonnegative integer counts required')
    maximum = int(counts.max())
    q = means/(sizes+means)
    if np.any(q <= 0) or np.any(q >= 1):
        raise ArithmeticError('Unresolved NB generating probabilities')
    power = np.ones_like(means)
    b = np.zeros((len(counts), maximum+1), np.longdouble)
    h = np.zeros_like(b)
    h[:, 0] = 1
    for n in range(1, maximum+1):
        power *= q
        b[:, n] = np.sum(sizes*power, axis=1)
    for n in range(1, maximum+1):
        h[:, n] = np.sum(b[:, 1:n+1]*h[:, n-1::-1], axis=1)/n
    chosen = h[np.arange(len(counts)), counts.astype(int)]
    if not np.isfinite(chosen).all() or np.any(chosen <= 0):
        raise ArithmeticError('Unresolved positive NB-sum coefficient')
    logzero = -np.sum(sizes*np.log1p(means/sizes), axis=1)
    logcoefficient = np.log(chosen)
    value = logzero+logcoefficient
    if not np.isfinite(value).all():
        raise ArithmeticError('Nonfinite NB-sum log probability')
    # Conservative accumulation allowance at input float64 precision, retaining
    # the magnitudes of the two terms that cancel in the returned log PMF.
    operations = 16*(counts+1.)**2+16*means.shape[1]*(counts+1.)+32
    eps = np.finfo(float).eps
    budget = operations*eps/(1-operations*eps)*(1+np.abs(logzero)+np.abs(logcoefficient))
    return np.asarray(value, float), np.asarray(budget, float)


def convolved_logpmf(count, means, sizes):
    """Independent direct PMF polynomial multiplication, no PGF recurrence."""
    polynomial = np.array([1.])
    support = np.arange(int(count)+1)
    for mean, size in zip(means, np.broadcast_to(sizes, np.shape(means))):
        component = stats.nbinom.pmf(support, size, size/(size+mean))
        polynomial = np.convolve(polynomial, component)[:int(count)+1]
    probability = polynomial[int(count)]
    if not np.isfinite(probability) or probability <= 0:
        raise ArithmeticError('Independent convolution unresolved')
    return float(np.log(probability))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    rec = start_run(a.out, dict(script=Path(__file__), design=a.task/'DESIGN.md'), {},
        dict(scope='Independent numerical NB-sum gate, biological n=0'), mode='gene_bin_reference')
    try:
        checks = []
        counts = np.arange(21)
        means = np.tile([1., 4.], (len(counts), 1))
        values, budgets = summed_logpmf(counts, means, np.array([1., 2.]))
        for n, value, budget in zip(counts, values, budgets):
            probability = sum(Fraction(1, 2)*Fraction(1, 2)**k *
                (n-k+1)*Fraction(1, 3)**2*Fraction(2, 3)**(n-k) for k in range(n+1))
            expected = math.log(float(probability))
            assert abs(value-expected) <= budget
            checks.append(dict(kind='Fraction', count=int(n), error=abs(value-expected), budget=budget))
        for q in (.2, .7):
            sizes = np.array([1., 2., 7.])
            means = np.tile(sizes*q/(1-q), (7, 1))
            counts = np.array([0, 1, 5, 10, 30, 60, 120])
            values, budgets = summed_logpmf(counts, means, sizes)
            reference = stats.nbinom.logpmf(counts, sizes.sum(), 1-q)
            for n, value, expected, budget in zip(counts, values, reference, budgets):
                assert abs(value-expected) <= budget
                checks.append(dict(kind='NB_closure', count=int(n), error=abs(value-expected), budget=budget))
        counts = np.array([0, 1, 5, 20, 60, 120])
        means = np.tile([.2, 3., 20.], (len(counts), 1))
        sizes = np.array([2., 10., .5])
        values, budgets = summed_logpmf(counts, means, sizes)
        permuted, _ = summed_logpmf(counts[::-1], means[::-1, ::-1], sizes[::-1])
        assert np.all(np.abs(values-permuted[::-1]) <= budgets)
        for n, value, row, budget in zip(counts, values, means, budgets):
            expected = convolved_logpmf(n, row, sizes)
            assert abs(value-expected) <= budget
            checks.append(dict(kind='positive_convolution', count=int(n), error=abs(value-expected), budget=budget))
        # Two geometric genes, two global rate scales, equal sums. Exact
        # likelihood ratios differ, so sums are not sufficient for this model.
        ratios = []
        for y in ((2, 0), (0, 2)):
            pmfs = []
            for means in ((1, 4), (2, 8)):
                pmfs.append(math.prod(Fraction(1, 1+m)*Fraction(m, 1+m)**n for m, n in zip(means, y)))
            ratios.append(pmfs[1]/pmfs[0])
        assert ratios[0] != ratios[1]
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True, checks=checks,
            unequal_gene_NB_sufficiency_counterexample=dict(count_vectors=[[2, 0], [0, 2]],
                class_means=[[1, 4], [2, 8]], sizes=[1, 1], exact_likelihood_ratios=[str(v) for v in ratios]),
            longdouble_epsilon=str(np.finfo(np.longdouble).eps), biological_specimens=0,
            production_changed=False, scientific_criteria_changed=False), indent=2)+'\n')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU Fraction, NB closure and independent positive PMF convolution')


if __name__ == '__main__':
    main()
