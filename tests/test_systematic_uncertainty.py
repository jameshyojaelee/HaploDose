"""Independent checks of the optional model-conditional simultaneous bounds."""
import numpy as np
import pytest
from scipy import optimize, stats

from haplodose.uncertainty import systematic_rna_intervals


def calculate(y, E, S, reset, intervals, **kwargs):
    return systematic_rna_intervals(y, E, S, reset, intervals,
        counts_kind='raw', count_provenance='raw synthetic Poisson counts without transformation',
        expectation_provenance='independent deterministic test expectation',
        covariance_provenance='known test generator covariance',
        bin_universe_provenance='fixed test bin universe and reset boundaries', **kwargs)


def test_garwood_components_against_independent_poisson_cdf_inversion():
    # Include the experiment's count scale and its 120-bin tail allocation.
    y = np.resize(np.array([0, 1, 3, 9, 2000]), 120)
    result = calculate(y, np.ones(120), np.zeros((120, 120)), np.r_[True, np.zeros(119, bool)],
                       [[j, j] for j in range(5)])
    tail = result['count_two_sided_tail']
    for count, row in zip(y, result['intervals']):
        bracket = max(1000., 2. * (count + 1))
        lower = 0. if count == 0 else optimize.brentq(
            lambda rate: stats.poisson.sf(count - 1, rate) - tail, 0., bracket,
            xtol=np.finfo(float).eps, rtol=4 * np.finfo(float).eps)
        upper = optimize.brentq(lambda rate: stats.poisson.cdf(count, rate) - tail, 0., bracket,
                               xtol=np.finfo(float).eps, rtol=4 * np.finfo(float).eps)
        # Independent CDF inversion and gamma quantiles have separate rounding.
        np.testing.assert_allclose(np.exp(row['count_log_lower']), lower,
                                   rtol=64 * np.finfo(float).eps, atol=np.finfo(float).eps)
        np.testing.assert_allclose(np.exp(row['count_log_upper']), upper,
                                   rtol=64 * np.finfo(float).eps)
        if count == 0:
            assert row['lower'] == -np.inf and row['width'] == np.inf


def test_covariance_quadratic_form_family_and_error_allocation():
    factor = np.array([[1., 0., 0.], [.5, 1., 0.], [-.5, .25, 1.]]) / 10
    S = factor @ factor.T
    pairs = [[0, 0], [0, 1], [0, 2], [1, 1], [1, 2], [2, 2]]
    result = calculate([10, 20, 30], [15, 15, 15], S, np.array([True, False, False]), pairs)
    assert result['count_family_size'] == 3
    assert result['latent_family_size'] == 6
    assert result['alpha_count'] + result['alpha_latent'] == 1 - result['confidence']
    critical = stats.norm.isf(result['alpha_latent'] / 12)
    assert result['latent_critical_value'] == critical
    for pair, row in zip(pairs, result['intervals']):
        a = np.zeros(3)
        a[pair[0]:pair[1] + 1] = 1 / row['n_bins']
        direct = float(a @ S @ a)
        padded = S + np.eye(3) * result['covariance_diagonal_padding']
        sum_roundoff = (row['n_bins'] ** 2 * np.finfo(float).eps
                        / (1 - row['n_bins'] ** 2 * np.finfo(float).eps))
        upper = float(a @ padded @ a) + sum_roundoff * float(np.abs(padded).sum())
        assert direct <= row['latent_mean_variance_upper'] <= upper
        assert row['lower'] == row['count_log_lower'] - row['latent_halfwidth']
        assert row['upper'] == row['count_log_upper'] + row['latent_halfwidth']


def test_selection_and_gaps_do_not_shrink_the_family():
    y, E = np.array([1, 30, 3, 80]), np.full(4, 10.)
    reset = np.array([True, False, True, False])
    S = np.eye(4) * .1
    pairs = [[0, 0], [0, 1], [1, 1], [2, 2], [2, 3], [3, 3]]
    all_bounds = calculate(y, E, S, reset, pairs)
    selected = max(pairs, key=lambda pair: abs(np.log((y[pair[0]:pair[1] + 1] + .5) / 10.5).mean()))
    one = calculate(y, E, S, reset, [selected])
    assert one['latent_family_size'] == all_bounds['latent_family_size'] == 6
    assert one['intervals'][0] == all_bounds['intervals'][pairs.index(selected)]
    with pytest.raises(ValueError, match='gap reset'):
        calculate(y, E, S, reset, [[1, 2]])


def test_single_bin_count_component_exact_discrete_coverage():
    # Enumerate the Poisson law rather than checking simulated execution success.
    eps = np.finfo(float).eps
    for rate in (.1, 1., 5., 20.):
        maximum = int(stats.poisson.isf(eps, rate))
        covered_probability = 0.
        for count in range(maximum + 1):
            result = calculate([count], [1.], [[0.]], np.array([True]), [[0, 0]])
            row = result['intervals'][0]
            if row['lower'] <= np.log(rate) <= row['upper']:
                covered_probability += stats.poisson.pmf(count, rate)
        roundoff = (maximum + 1) * eps / (1 - (maximum + 1) * eps)
        assert covered_probability >= 1 - result['alpha_count'] - eps - roundoff


def test_covariance_validation_without_automatic_repair():
    reset = np.array([True, False])
    with pytest.raises(ValueError, match='semidefinite'):
        calculate([1, 2], [1, 1], [[1., 2.], [2., 1.]], reset, [[0, 1]])
    with pytest.raises(ValueError, match='symmetric'):
        calculate([1, 2], [1, 1], [[1., .1], [.2, 1.]], reset, [[0, 1]])
    with pytest.raises(ValueError, match='independent resolution'):
        calculate([1, 2], [1, 1], np.diag([1., -np.finfo(float).eps]), reset, [[0, 1]])
    # Rank-deficient PSD covariance is permitted; an exactly zero covariance
    # adds no latent uncertainty or arbitrary numerical floor.
    singular = calculate([1, 2], [1, 1], [[1., -1.], [-1., 1.]], reset, [[0, 1]])
    assert singular['intervals'][0]['latent_mean_variance_upper'] >= 0
    zero = calculate([1, 2], [1, 1], np.zeros((2, 2)), reset, [[0, 1]])
    assert zero['intervals'][0]['latent_halfwidth'] == 0
    assert zero['covariance_diagonal_padding'] == 0


@pytest.mark.parametrize('counts,expectation,covariance,reset,intervals', [
    ([1.5], [1.], [[0.]], [True], [[0, 0]]),
    ([np.nan], [1.], [[0.]], [True], [[0, 0]]),
    ([1], [0.], [[0.]], [True], [[0, 0]]),
    ([1], [np.inf], [[0.]], [True], [[0, 0]]),
    ([1], [1.], [[np.nan]], [True], [[0, 0]]),
    ([1], [1.], [[0.]], [False], [[0, 0]]),
    ([1], [1.], [[0.]], [True], [[0., 0.]]),
    ([1], [1.], [[0.]], [True], [[0, 1]]),
])
def test_invalid_inputs(counts, expectation, covariance, reset, intervals):
    with pytest.raises(ValueError):
        calculate(counts, expectation, covariance, np.asarray(reset), intervals)


def test_explicit_external_input_provenance_required():
    with pytest.raises(ValueError, match='expectation_provenance'):
        systematic_rna_intervals([1], [1.], [[0.]], np.array([True]), [[0, 0]],
                                 counts_kind='raw', count_provenance='raw test counts',
                                 expectation_provenance='', covariance_provenance='known generator',
                                 bin_universe_provenance='fixed bins')
    with pytest.raises(TypeError):
        systematic_rna_intervals([1], [1.], [[0.]], np.array([True]), [[0, 0]])


def test_transformed_integers_and_undeclared_measurements_cannot_receive_intervals():
    args = ([10], [10.], [[0.]], np.array([True]), [[0, 0]])
    external = dict(expectation_provenance='fixed E', covariance_provenance='known covariance',
                    bin_universe_provenance='fixed bin universe')
    with pytest.raises(TypeError):
        systematic_rna_intervals(*args, **external)
    for kind in ('normalized', 'reconstructed', 'unknown', None):
        with pytest.raises(ValueError, match='explicitly be raw'):
            systematic_rna_intervals(*args, counts_kind=kind, count_provenance='integer-valued transformed data',
                                     **external)
    with pytest.raises(ValueError, match='count_provenance'):
        systematic_rna_intervals(*args, counts_kind='raw', count_provenance='', **external)
    result = systematic_rna_intervals(*args, counts_kind='raw', count_provenance='raw test counts', **external)
    assert result['counts_kind'] == 'raw' and result['count_provenance'] == 'raw test counts'
