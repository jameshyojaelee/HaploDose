"""A finite refinement result must describe one stored allele parameter state."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.special import softmax
from scipy.stats import binom

from haplodose.clones import SpecimenModel


def allele_controller(rounds):
    # Isolate the controller with an exactly specified binomial allele model.
    # No convergence or biological calibration is presumed. Refitting allele
    # fractions after scoring changes the observations' likelihoods.
    a = np.array([[9.], [8.], [2.], [1.]])
    d = np.full_like(a, 10.)
    m = SpecimenModel.__new__(SpecimenModel)
    m.params = SimpleNamespace(outer_rounds=rounds, normal_prior=.05, min_cells=0)
    m.J = 1
    m.R_prev = None
    m.fixed_gamma = True
    m.log = []
    m.lam = 0.
    m.B = np.array([[.5], [.8]])
    m.Ak = m.Dk = m.allele_post = None
    m.alleles = SimpleNamespace(pool=lambda R: (R.T @ a, R.T @ d),
        fit=lambda Ak, Dk: (np.zeros((len(Ak), 1, 9)), Ak / Dk))

    def loglik(doses):
        b = m.B[:, 0]
        return a * np.log(b[None] / .5) + (d-a) * np.log((1-b[None]) / .5)

    m.loglik = loglik
    m.posterior = lambda LL, pi: softmax(LL + np.log(pi)[None], axis=1)
    m.pooled = lambda R: (R.T @ d, R.sum(axis=0), np.zeros((2, 2)))
    m.shared_profiles = lambda *args, **kw: dict(dose=np.zeros((2, 1)),
        E=np.ones((2, 1)), post=np.ones((2, 1, 1)), path=np.zeros((2, 1)), tau=np.ones(2))
    m.has_events = lambda *args: True
    m.allele_events = lambda k: False
    m._update_offset = m._update_phi = m.update_scores = lambda *args: None
    return m, a, d


@pytest.mark.parametrize('rounds,stop_on_map', [(1, False), (8, False), (8, True)])
def test_returned_allele_scores_responsibilities_and_pools_share_snapshot(rounds, stop_on_map):
    m, a, d = allele_controller(rounds)
    doses, pi, R, LL, Yk, Lk, XLk = m.em(
        [np.zeros(1), np.zeros(1)], np.array([.5, .5]), stop_on_map=stop_on_map)
    fresh = m.loglik(doses)
    np.testing.assert_array_equal(LL, fresh)
    np.testing.assert_array_equal(R, softmax(fresh + np.log(pi)[None], axis=1))
    for actual, expected in zip((Yk, Lk, XLk), m.pooled(R)):
        np.testing.assert_array_equal(actual, expected)
    independent = binom.logpmf(a, d, m.B[:, 0][None]) - binom.logpmf(a, d, .5)
    budget = 64*np.finfo(float).eps*float(d.sum())
    assert np.max(abs(LL-independent)) <= budget
    # Frozen allele parameters come from their recorded fit pools. They need
    # not be a fixed point of the returned responsibilities after finite work.
    np.testing.assert_array_equal(m.B[1:], (m.Ak/m.Dk)[1:])
    assert m.log[-1]['convergence_status'].startswith('not_certified_')
