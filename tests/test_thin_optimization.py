"""Independent numerical references for the proposal-region optimizer."""
import numpy as np
import pytest
from scipy.optimize import minimize, minimize_scalar
from scipy.special import expit, logsumexp
from scipy.stats import nbinom

from haplodose.clones import _event_mixture_objective, _fit_event_mixture


def fixture():
    rng = np.random.default_rng(71821)
    E = rng.uniform(8., 35., (80, 6))
    r = np.broadcast_to(np.linspace(3., 19., 6), E.shape)
    theta = np.r_[np.full(40, .45), np.full(40, 1.8)]
    O = rng.negative_binomial(r, r / (r + E * theta[:, None])).astype(float)
    return O, E, r


def independent_component(eta, O, E, r):
    return nbinom.logpmf(O, r, r / (r + E * np.exp(eta))).sum(axis=1)


def test_event_gain_matches_full_pmf_independent_optimization():
    O, E, r = fixture()
    fitted = _fit_event_mixture(O, E, r, .25)
    null = minimize_scalar(lambda eta: -independent_component(eta, O, E, r).sum(),
                           bounds=(-10., 10.), method='bounded',
                           options={'xatol': np.sqrt(np.finfo(float).eps)})
    null_ll = independent_component(null.x, O, E, r).sum()

    def objective(x):
        ll = np.stack([independent_component(x[0], O, E, r) - np.logaddexp(0., x[2]),
                       independent_component(x[1], O, E, r) - np.logaddexp(0., -x[2])])
        return -logsumexp(ll, axis=0).sum()

    bound = np.log(999.)
    reference = minimize(objective, [np.log(.45), np.log(1.8), 0.], method='Nelder-Mead',
                         bounds=[(-10., 10.), (-10., 10.), (-bound, bound)],
                         options={'maxiter': 3000, 'xatol': np.sqrt(np.finfo(float).eps),
                                  'fatol': np.finfo(float).eps * abs(null_ll)})
    assert reference.success
    # The full-pmf reference subtracts large quantities. Its roundoff budget is
    # deliberately distinct from the production, cancellation-free gain budget.
    eps = np.finfo(float).eps
    reference_budget = 16 * O.size * eps * max(1., abs(null_ll), abs(reference.fun))
    assert abs(fitted['gain'] - (-reference.fun - null_ll)) <= (
        reference_budget + fitted['numerical_gain_budget'])
    assert np.array_equal(fitted['responsibility'] > .5,
                          independent_component(reference.x[1], O, E, r)
                          - np.logaddexp(0., -reference.x[2])
                          > independent_component(reference.x[0], O, E, r)
                          - np.logaddexp(0., reference.x[2]))


def test_event_gradient_and_constant_cancellation():
    O, E, r = fixture()
    x, null = np.array([-.6, .5, -.3]), .1
    value, analytic, responsibility = _event_mixture_objective(x, O, E, r, null)
    baseline = independent_component(null, O, E, r)
    a = independent_component(x[0], O, E, r) - baseline - np.logaddexp(0., x[2])
    b = independent_component(x[1], O, E, r) - baseline - np.logaddexp(0., -x[2])
    eps = np.finfo(float).eps
    full_budget = 16 * O.size * eps * max(1., np.abs(baseline).sum())
    assert abs(value + np.logaddexp(a, b).sum()) <= full_budget
    np.testing.assert_allclose(responsibility, expit(b - a), rtol=O.size * eps,
                               atol=O.size * eps)
    for axis in range(3):
        h = np.cbrt(eps) * max(1., abs(x[axis]))

        def derivative(step):
            delta = np.zeros(3)
            delta[axis] = step
            plus = _event_mixture_objective(x + delta, O, E, r, null)[0]
            minus = _event_mixture_objective(x - delta, O, E, r, null)[0]
            return (plus - minus) / (2 * step), max(abs(plus), abs(minus), 1.)

        coarse, magnitude = derivative(h)
        fine, _ = derivative(h / 2)
        # Richardson's estimate of central-difference truncation plus floating
        # subtraction error; no common tolerance shared with the optimizer.
        allowance = 4 * abs(fine - coarse) / 3 + 16 * eps * magnitude / h
        assert abs(fine - analytic[axis]) <= allowance


def test_event_invariant_to_exposure_scale_and_cell_order():
    O, E, r = fixture()
    original = _fit_event_mixture(O, E, r, .25)
    scaled = _fit_event_mixture(O[::-1], E[::-1] * 1e8, r[::-1], .25)
    allowance = original['numerical_gain_budget'] + scaled['numerical_gain_budget']
    assert abs(original['gain'] - scaled['gain']) <= allowance
    assert np.array_equal(original['responsibility'] > .5,
                          scaled['responsibility'][::-1] > .5)


def test_event_empty_counts_invalid_support_and_iteration_safeguard():
    O, E, r = fixture()
    assert _fit_event_mixture(np.zeros_like(O), E, r, .25)['gain'] == 0
    with pytest.raises(ValueError, match='supported'):
        _fit_event_mixture(O, np.zeros_like(E), r, .25)
    for invalid in (0, -1, 1.5, True):
        with pytest.raises(ValueError, match='positive integer'):
            _fit_event_mixture(O, E, r, .25, maxiter=invalid)


def test_line_search_failure_retries_without_bypassing_numerical_checks(monkeypatch):
    import scipy.optimize
    actual = scipy.optimize.minimize
    methods = []
    def reported_failure(*args, **kwargs):
        result = actual(*args, **kwargs)
        methods.append(kwargs.get('method'))
        if kwargs.get('method') == 'L-BFGS-B':
            result.success = False
            result.message = 'Injected line-search status failure at the computed solution'
        return result
    O, E, r = fixture()
    reference = _fit_event_mixture(O, E, r, .25)
    monkeypatch.setattr(scipy.optimize, 'minimize', reported_failure)
    fitted = _fit_event_mixture(O, E, r, .25)
    assert methods == ['L-BFGS-B', 'SLSQP']
    assert fitted['solver'] == 'SLSQP_retry'
    assert abs(reference['gain'] - fitted['gain']) <= reference['numerical_gain_budget'] + fitted['numerical_gain_budget']
    assert fitted['local_quadratic_gain_error'] <= fitted['numerical_gain_budget']
    assert np.all(np.abs(fitted['projected_score']) <= fitted['score_budget'])
