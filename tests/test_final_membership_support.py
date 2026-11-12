"""A last profile update must not export components below the existing support."""
from types import SimpleNamespace

import numpy as np
from scipy.special import softmax

from haplodose.clones import SpecimenModel


def controller(minimum=1, final_loglik=(0., np.log(2.), -np.log(16.))):
    # The initial three equal-rate components each have mass4/3. The sole
    # profile update changes the likelihood; its *returned* responsibilities
    # have an independently calculable, different support. No convergence or
    # biological calibration is assumed by this deterministic controller.
    m = SpecimenModel.__new__(SpecimenModel)
    m.params = SimpleNamespace(outer_rounds=1, normal_prior=.05, min_cells=minimum)
    m.J = 1
    m.alleles = m.B = m.R_prev = None
    m.fixed_gamma = True
    m.log = []
    m.lam = 0.
    m.loglik = lambda doses: np.tile(
        [final_loglik[int(d[0])] for d in doses], (4, 1))
    m.posterior = lambda LL, pi: softmax(LL+np.log(pi)[None], axis=1)
    observations = np.array([[3.], [5.], [7.], [11.]])
    depth_moments = np.column_stack([np.arange(4), np.arange(4)**2])
    m.pooled = lambda R: (R.T@observations, R.sum(axis=0), R.T@depth_moments)
    m.shared_profiles = lambda Y, *args, **kw: dict(
        dose=np.arange(len(Y))[:, None], E=np.ones_like(Y),
        post=np.ones((len(Y), 1, 1)), path=np.zeros_like(Y), tau=np.ones(len(Y)))
    m.has_events = lambda *args: True
    m.allele_events = lambda k: False
    m._update_offset = m._update_phi = m.update_scores = lambda *args: None
    return m


def fit(m):
    return m.em([np.zeros(1) for _ in range(3)], np.full(3, 1/3), stop_on_map=False)


def test_returned_responsibilities_meet_existing_component_minimum():
    m = controller()
    doses, pi, R, LL, Y, L, XL = fit(m)
    assert len(doses) == 2
    assert np.all(R.sum(axis=0)[1:] >= m.params.min_cells)
    # Deleting the unsupported component leaves a1:2 likelihood ratio and
    # equal retained priors. This is exact independently of a refinement loop.
    expected = np.tile([1/3, 2/3], (4, 1))
    allowance = 32*np.finfo(float).eps
    assert np.max(abs(R-expected)) <= allowance
    np.testing.assert_array_equal(LL, m.loglik(doses))
    np.testing.assert_array_equal(R, softmax(LL+np.log(pi)[None], axis=1))
    for actual, independent in zip((Y, L, XL), m.pooled(R)):
        np.testing.assert_array_equal(actual, independent)
    event = next(x for x in m.log if x['action'] == 'drop_small_final')
    assert event['clone'] == 2 and 0 < event['n_eff'] < 1
    assert event['minimum'] == 1


def test_terminal_support_can_retain_only_the_normal_component():
    m = controller(final_loglik=(0., -np.log(16.), -np.log(32.)))
    doses, pi, R, LL, Y, L, XL = fit(m)
    assert len(doses) == 1
    np.testing.assert_array_equal(doses, [[0.]])
    np.testing.assert_array_equal(pi, [1.])
    np.testing.assert_array_equal(R, np.ones((4, 1)))
    np.testing.assert_array_equal(LL, m.loglik(doses))
    for actual, expected in zip((Y, L, XL), m.pooled(R)):
        np.testing.assert_array_equal(actual, expected)


def test_zero_minimum_does_not_remove_positive_uncertain_components():
    m = controller(minimum=0)
    doses, pi, R, LL, *_ = fit(m)
    assert len(doses) == 3
    assert 0 < R[:, 2].sum() < 1
    assert not any(x['action'] == 'drop_small_final' for x in m.log)
    np.testing.assert_array_equal(R, softmax(LL+np.log(pi)[None], axis=1))


def test_normal_component_is_retained_even_below_tumour_minimum():
    m = controller(final_loglik=(-np.log(16.), np.log(2.), 0.))
    doses, pi, R, LL, *_ = fit(m)
    assert len(doses) == 3 and R[:, 0].sum() < 1
    assert np.all(R[:, 1:].sum(axis=0) >= 1)
    np.testing.assert_array_equal(doses[0], [0.])
    np.testing.assert_array_equal(R, softmax(LL+np.log(pi)[None], axis=1))


def test_exact_dyadic_support_at_minimum_remains_eligible():
    m = controller()
    m.posterior = lambda LL, pi: np.tile([.25, .25, .5], (4, 1))
    doses, _, R, *_ = fit(m)
    assert len(doses) == 3
    assert R[:, 1].sum() == m.params.min_cells
    assert not any(x['action'] == 'drop_small_final' for x in m.log)
