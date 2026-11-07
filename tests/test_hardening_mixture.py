"""Independent simplex-MLE, optimality and invariance checks for references."""
import numpy as np
import pytest
from scipy import optimize

from haplodose.clones import mixture_weights


def normalized(values):
    return values / values.sum()


def test_reference_absorption_boundary_reaches_independent_poisson_optimum():
    J = 120
    base = np.full(J, 1 / J)
    alternate = normalized(base * np.exp(np.r_[np.full(J // 4, .6 * np.log(1.5)), np.zeros(3 * J // 4)]))
    mu = np.stack([base, alternate], axis=1)
    counts = np.rint(500000 * alternate)
    counts[:J // 4] += 10
    counts[-J // 4:] -= 10
    objective = lambda weight: -float(counts @ np.log(mu @ np.array([1 - weight, weight])))
    independent = optimize.minimize_scalar(objective, bounds=(0, 1), method='bounded')
    endpoint = min([0., independent.x, 1.], key=objective)
    assert endpoint == 1.
    weights, diagnostic = mixture_weights(counts, mu, return_diagnostics=True)
    actual_gap = objective(weights[1]) - objective(endpoint)
    reduction = (J + 16) * np.finfo(float).eps
    objective_roundoff = reduction * abs(objective(endpoint))
    assert actual_gap <= counts.sum() * diagnostic['float64_roundoff_budget'] + objective_roundoff
    assert diagnostic['converged']
    assert diagnostic['normalized_objective_gap_upper'] <= diagnostic['float64_roundoff_budget']
    assert weights[0] == 0.


def test_known_interior_solution_and_common_count_scale():
    mu = np.eye(3)
    counts = np.array([20., 30., 50.])
    for scale in (1e-200, 1., 1e200):
        weights, diagnostic = mixture_weights(counts * scale, mu, return_diagnostics=True)
        # Here the likelihood gap equals KL(p||weights); Pinsker gives this
        # probability-error budget directly from the objective certificate.
        probability_budget = np.sqrt(2 * diagnostic['float64_roundoff_budget'])
        assert np.linalg.norm(weights - normalized(counts), ord=1) <= probability_budget
        assert diagnostic['normalized_objective_gap_upper'] <= diagnostic['float64_roundoff_budget']


def test_identical_references_are_symmetric_not_identified():
    column = normalized(np.array([1., 2., 3., 4.]))
    mu = np.tile(column[:, None], (1, 4))
    weights, diagnostic = mixture_weights(np.array([10., 5., 30., 20.]), mu, return_diagnostics=True)
    np.testing.assert_array_equal(weights, np.full(4, .25))
    assert diagnostic['iterations'] == 0
    np.testing.assert_array_equal(mu @ weights, column)


def test_reference_permutation_and_redundancy_preserve_fitted_composition():
    first = np.array([.6, .3, .1])
    second = np.array([.1, .3, .6])
    mu = np.stack([first, second, .5 * (first + second)], axis=1)
    truth = .3 * first + .7 * second
    counts = truth * 10000
    weights, diagnostic = mixture_weights(counts, mu, return_diagnostics=True)
    order = np.array([2, 0, 1])
    permuted, other = mixture_weights(counts, mu[:, order], return_diagnostics=True)
    budget = np.sqrt(2 * diagnostic['float64_roundoff_budget']) + np.sqrt(2 * other['float64_roundoff_budget'])
    assert np.linalg.norm(mu @ weights - mu[:, order] @ permuted, ord=1) <= budget


def test_zero_count_bins_zero_reference_rates_and_extreme_support():
    weights, diagnostic = mixture_weights(np.array([10., 0.]), np.eye(2), return_diagnostics=True)
    # At this boundary h(optimum)-h(w)=-log(w0). Certification allows a
    # roundoff-sized remainder; exact bitwise zero is not an optimizer contract.
    assert -np.log(weights[0]) <= diagnostic['float64_roundoff_budget']
    assert weights[1] <= -np.expm1(-diagnostic['float64_roundoff_budget'])
    mu = np.array([[1., 1.], [1e-300, 2e-300]])
    weights, diagnostic = mixture_weights(np.array([1., 1.]), mu, return_diagnostics=True)
    # The boundary directional gradient margin is 1/4, so concavity bounds
    # residual weight on the inferior component by gap/(1/4).
    assert weights[0] <= 4 * diagnostic['float64_roundoff_budget']
    actual_gap = .5 * np.log1p((mu[:, 1] - mu @ weights) / (mu @ weights)).sum()
    assert actual_gap <= diagnostic['float64_roundoff_budget']
    assert np.isfinite(diagnostic['normalized_log_likelihood'])
    # Normalization must not overflow when the finite individual counts sum
    # above float64's representable range.
    np.testing.assert_array_equal(mixture_weights(np.full(3, 1e308), np.eye(3)), np.full(3, 1 / 3))


@pytest.mark.parametrize('counts,mu,message', [
    (np.zeros(2), np.eye(2), 'positive total'),
    (np.array([1., -1.]), np.eye(2), 'finite, nonnegative'),
    (np.array([1., np.nan]), np.eye(2), 'finite, nonnegative'),
    (np.ones(2), np.array([[1., 1.], [0., 0.]]), 'no support'),
    (np.ones(2), np.eye(2) * 2, 'sum to one'),
    (np.ones(2), np.array([[1., 0.], [0., 0.]]), 'sum to one'),
    (np.ones(2), np.ones((3, 2)), 'aligned'),
])
def test_mixture_rejects_invalid_measurements(counts, mu, message):
    with pytest.raises(ValueError, match=message):
        mixture_weights(counts, mu)


def test_iteration_safeguard_is_not_reported_as_convergence():
    mu = np.array([[.7, .05, .1], [.2, .8, .2], [.1, .15, .7]])
    with pytest.raises(RuntimeError, match='converge'):
        mixture_weights(np.array([31., 27., 42.]), mu, iters=1)
