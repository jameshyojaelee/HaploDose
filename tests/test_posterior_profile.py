from fractions import Fraction
import itertools
import math

import numpy as np
import pytest

from haplodose import posterior


def test_covering_interval_sums_round_once_after_exact_cancellation():
    tiny = np.nextafter(0., 1.)
    intervals = np.array([[tiny, 2.**60, 1., -2.**60],
                          [0., 0., tiny, 0.],
                          [0., 0., -1., tiny],
                          [0., 0., 0., 2.**60]])
    exact = [float(sum(Fraction(float(intervals[a, b]))
                       for a in range(j+1) for b in range(j, len(intervals))))
             for j in range(len(intervals))]
    np.testing.assert_array_equal(posterior._sum_covering_intervals(intervals), exact)


def test_posterior_runs_against_independent_path_enumeration():
    prior = np.array([.2, .5, .3])
    change = .1
    emission = np.log(np.array([[.8, .2, .5], [.1, .6, .3], [.2, .5, .7]]))
    actual, diagnostic = posterior._run_probabilities(emission, prior, change)
    # Compute each full path without the forward/backward recurrence.
    paths = list(itertools.product(range(3), repeat=3))
    weights = []
    for path in paths:
        p = prior[path[0]]
        for j, state in enumerate(path):
            p *= math.exp(emission[j, state])
            if j:
                p *= (1-change if state == path[j-1]
                      else change*prior[state]/(1-prior[path[j-1]]))
        weights.append(p)
    normalizer = math.fsum(weights)
    expected = np.zeros((3, 3))
    for path, weight in zip(paths, weights):
        lo = 0
        for hi in range(3):
            if hi == 2 or path[hi+1] != path[hi]:
                expected[lo, hi] += weight/normalizer
                lo = hi+1
    operations = 128*len(paths)*len(emission)
    gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    assert np.max(abs(actual-expected)) <= gamma
    assert diagnostic['maximum_mass_error'] <= diagnostic['mass_rounding_budget']


def test_dyadic_weighted_interval_medians_against_fraction():
    values = np.array([0., 2., 1., 4., -1., 5.])
    variance = np.array([1., 4., 1., 9., 4., 1.])
    bins = np.repeat(np.arange(3), 2)
    actual = posterior._interval_summaries(values, variance, bins, 3)
    for lo in range(3):
        for hi in range(lo, 3):
            keep = (bins >= lo) & (bins <= hi)
            pairs = sorted((float(y), Fraction(float(1/math.sqrt(v))))
                           for y, v in zip(values[keep], variance[keep]))
            total = sum(w for _, w in pairs)
            cumulative = Fraction(0)
            for j, (value, weight) in enumerate(pairs):
                cumulative += weight
                if 2*cumulative >= total:
                    expected = (value+pairs[j+1][0])/2 if 2*cumulative == total else value
                    break
            assert actual[1, lo, hi] == expected
            assert actual[0, lo, hi] == np.median(values[keep])
    np.testing.assert_array_equal(actual, posterior._interval_summaries(
        values[::-1], variance[::-1], bins[::-1], 3))
    np.testing.assert_array_equal(actual[1:2], posterior._interval_summaries(
        values, variance, bins, 3, weighted_only=True))


def prepared_inputs():
    return dict(gene_log_ratios=np.array([-.2, 0., .4, .8, -.1, .2]),
        gene_expected=np.full(6, 4.), gene_bins=np.repeat(np.arange(3), 2),
        reference_log_ratios=np.array([[0., .1, -.2, .3, .2, -.1],
            [.1, -.1, .1, -.2, .1, -.1], [-.1, .2, 0., .1, -.2, .2]]),
        reference_expected=np.full((3, 6), 4.), blocks=np.array(['p', 'p', 'q']),
        reference_id='prespecified-reference', counts_kind='raw')


def test_public_profile_preserves_gene_and_reference_correspondence():
    inputs = prepared_inputs()
    profile, diagnostic = posterior.posterior_rna_profile(**inputs)
    order = np.array([5, 1, 3, 0, 4, 2])
    permuted = dict(inputs)
    for key in ('gene_log_ratios', 'gene_expected', 'gene_bins'):
        permuted[key] = inputs[key][order]
    for key in ('reference_log_ratios', 'reference_expected'):
        permuted[key] = inputs[key][::-1, :][:, order]
    other, _ = posterior.posterior_rna_profile(**permuted)
    np.testing.assert_array_equal(profile, other)
    assert profile.shape == (3,) and np.isfinite(profile).all()
    assert diagnostic['units'] == 'natural_log_relative_rna_rate'
    assert diagnostic['calibration_status'] == 'unvalidated'
    reconstructed = dict(inputs, counts_kind='reconstructed')
    np.testing.assert_array_equal(profile, posterior.posterior_rna_profile(**reconstructed)[0])


@pytest.mark.parametrize('key,value', [
    ('gene_bins', np.array([0., 0., 1., 1., 2., 2.])),
    ('gene_bins', np.array([0, 0, 1, 1, 3, 3])),
    ('blocks', np.array(['p', 'q', 'p'])),
    ('blocks', np.array(['p', None, 'q'], dtype=object)),
    ('reference_id', ''), ('counts_kind', 'normalized'),
    ('gene_log_ratios', np.full(6, np.nan)),
    ('gene_expected', np.full(6, -1.)),
    ('reference_log_ratios', np.zeros((2, 6))),
    ('reference_expected', np.zeros((3, 5))),
])
def test_public_profile_rejects_unsupported_inputs(key, value):
    inputs = prepared_inputs()
    inputs[key] = value
    with pytest.raises(ValueError):
        posterior.posterior_rna_profile(**inputs)


def test_genomic_resets_factorize_under_fixed_working_parameters():
    e = np.log(np.array([[.8, .2], [.4, .6], [.1, .9]]))
    prior = np.array([.4, .6])
    both = posterior._chain_blocks(e, [True, False, True], prior, .1)
    first = posterior._chain_blocks(e[:2], [True, False], prior, .1)
    last = posterior._chain_blocks(e[2:], [True], prior, .1)
    np.testing.assert_array_equal(both[0], np.vstack([first[0], last[0]]))
    np.testing.assert_array_equal(both[2], np.r_[first[2], last[2]])
    assert both[1] == math.fsum([first[1], last[1]])


def test_reference_average_aligns_regions_and_uses_exact_sum():
    arrays = [[1e16, 2., np.nan], [-1e16, np.nan, 4.], [6., 1., np.nan]]
    axes = [['a', 'b', 'c'], ['a', 'c', 'b'], ['b', 'a', 'c']]
    result, diagnostic = posterior.average_rna_profiles(arrays, axes, reference_ids=['r1', 'r2', 'r3'])
    expected = [float((Fraction(10**16)-10**16+1)/3), 4., np.nan]
    np.testing.assert_array_equal(result, expected)
    assert diagnostic['region_ids'] == ['a', 'b', 'c']
    # A permutation of fits returns that first fit's region order, with the
    # same identity-aligned estimator and exact compensated mean.
    result2, _ = posterior.average_rna_profiles(arrays[1:]+arrays[:1], axes[1:]+axes[:1],
        reference_ids=['r2', 'r3', 'r1'])
    np.testing.assert_array_equal(result2, np.asarray(expected)[[0, 2, 1]])


@pytest.mark.parametrize('arrays,axes,ids', [
    ([[1., 2.], [3., np.nan]], [['a', 'b'], ['a', 'b']], ['x', 'y']),
    ([[1., 2.], [3., 4.]], [['a', 'b'], ['a', 'c']], ['x', 'y']),
    ([[1., 2.], [3., 4.]], [['a', 'a'], ['a', 'b']], ['x', 'y']),
    ([[1., 2.], [3., 4.]], [['a', 'b'], ['a', 'b']], ['x', 'x']),
    ([[1., 2.], [3., 4.]], [['a', 'b'], ['a', 'b']], ['x', None]),
    ([[1., np.inf]], [['a', 'b']], ['x']),
])
def test_reference_average_never_drops_regions_or_failed_fits(arrays, axes, ids):
    with pytest.raises(ValueError):
        posterior.average_rna_profiles(arrays, axes, reference_ids=ids)
