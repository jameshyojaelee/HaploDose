"""Exact expression mixture redundancy must not depend on hard MAP sizes."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.special import logsumexp

from haplodose.clones import SpecimenModel, posteriors


def model(alleles=None, graph=None):
    m = SpecimenModel.__new__(SpecimenModel)
    m.alleles, m.graph, m.Ak = alleles, graph, None
    m.thin, m.prop_fold = True, 'prop'
    m.params, m.log = SimpleNamespace(min_cells=30), []
    # No pair has enough MAP cells in these small fixtures. No profile refit
    # should occur for a purely algebraic collapse.
    m.loglik = lambda doses, fold: np.zeros((4, len(doses)))
    m.pooled = lambda R: pytest.fail('Exact collapse must preserve the fitted rates')
    return m


def test_identical_profiles_collapse_even_without_two_eligible_map_groups():
    m = model()
    doses = [np.zeros(2)] + [np.array([.25, -.25]) for _ in range(4)]
    pi = np.array([.125, .25, .25, .25, .125])
    # An independent likelihood table with four identical tumour columns.
    LL = np.array([[-4., 0., 0., 0., 0.], [0., -2., -2., -2., -2.],
                   [-2., 1., 1., 1., 1.], [1., -1., -1., -1., -1.]])
    R = posteriors(LL, np.log(pi))
    before_R, before_pi = R.copy(), pi.copy()
    result, mass, merged = m.merge_unsupported(R, doses, pi)
    assert merged and len(result) == 2
    np.testing.assert_array_equal(result, doses[:2])
    np.testing.assert_array_equal(mass, [.125, .875])
    budget = 64 * LL.size * np.finfo(float).eps
    # Marginal likelihood and aggregate posterior are independently computed.
    np.testing.assert_allclose(logsumexp(LL + np.log(pi), axis=1),
                               logsumexp(LL[:, :2] + np.log(mass), axis=1), rtol=0, atol=budget)
    expected_R = np.column_stack((R[:, 0], R[:, 1:].sum(axis=1)))
    np.testing.assert_allclose(posteriors(LL[:, :2], np.log(mass)), expected_R,
                               rtol=0, atol=budget)
    counts = np.array([[2, 7], [11, 3], [5, 13], [17, 19]])
    np.testing.assert_allclose(expected_R.T @ counts,
                               np.vstack(((R.T @ counts)[0], (R.T @ counts)[1:].sum(axis=0))),
                               rtol=0, atol=budget * counts.sum())
    np.testing.assert_array_equal(R, before_R)
    np.testing.assert_array_equal(pi, before_pi)
    assert len(m.log) == 3
    assert all(r['action'] == 'merge_identical_profiles' for r in m.log)


@pytest.mark.parametrize('mode', ['allele', 'spatial', 'near', 'different', 'nonfinite'])
def test_no_collapse_when_expression_equivalence_is_insufficient_or_absent(mode):
    m = model(alleles=object() if mode == 'allele' else None,
              graph=object() if mode == 'spatial' else None)
    a = np.array([.25, -.25])
    b = a.copy()
    if mode == 'near':
        b[0] = np.nextafter(b[0], np.inf)
    elif mode == 'different':
        b *= 2
    elif mode == 'nonfinite':
        a[0] = b[0] = np.inf
    doses = [np.zeros(2), a, b]
    result, mass, merged = m.merge_unsupported(np.full((4, 3), 1/3), doses, np.array([.25, .5, .25]))
    assert not merged and len(result) == 3 and not m.log


def test_normal_component_is_never_collapsed_into_a_tumour():
    m = model()
    doses, pi, merged = m.merge_unsupported(np.full((4, 3), 1/3),
                                           [np.zeros(2) for _ in range(3)], np.array([.25, .5, .25]))
    assert merged and len(doses) == 2
    np.testing.assert_array_equal(pi, [.25, .75])


def test_exact_collapse_then_statistical_merge_still_refits_aggregated_pools():
    m = model()
    m.thin, m.J, m.allele_post = False, 1, None
    m.matched_differences = lambda labels, pairs: [dict(segments=[], distance=0., testable=True)]
    m.supported = lambda *args, **kwargs: False
    seen = []
    def pool(R):
        seen.append(R.copy())
        return np.zeros((2, 1)), np.ones(2), np.zeros((2, 2))
    m.pooled = pool
    m.profiles = lambda *args: dict(dose=np.array([[0.], [.75]]))
    R = np.tile([.1, .6, .2, .1], (60, 1))
    R[30:] = [.1, .1, .2, .6]
    doses, _, merged = m.merge_unsupported(R, [np.zeros(1), np.ones(1), np.ones(1), np.full(1, 2.)],
                                            np.array([.1, .3, .3, .3]))
    assert merged and len(doses) == 2 and len(seen) == 1
    np.testing.assert_array_equal(doses[1], [.75])
    budget = 64 * R.size * np.finfo(float).eps
    np.testing.assert_allclose(seen[0], np.column_stack((R[:, 0], R[:, 1:].sum(axis=1))),
                               rtol=0, atol=budget)
    assert [e['action'] for e in m.log] == ['merge_identical_profiles', 'merge_matched']
