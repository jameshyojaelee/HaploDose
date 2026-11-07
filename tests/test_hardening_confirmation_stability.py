"""Confirmation must abstain when all frozen candidate likelihoods underflow."""
import numpy as np

from haplodose.confirmation import FrozenDiscovery, score_confirmation


def artifact(rates):
    return FrozenDiscovery(log_rates=np.asarray(rates, dtype=float), phi=np.ones(4),
                           proposal_mask=np.array([True, True, False, False]),
                           test_mask=np.array([False, False, True, True]), contrasts=(),
                           discovery_cell_ids=('discovery',))


def test_positive_proposal_depth_with_no_finite_candidate_abstains():
    frozen = artifact([[0, -1000, 0, 0], [0, -2000, 0, 0]])
    scored = score_confirmation(np.array([[0, 1, 20, 0]]), ['confirmation'], frozen)
    assert scored['candidate'].tolist() == [-1]
    assert scored['unresolved_reason'].tolist() == ['no_finite_candidate_score']
    assert not np.isfinite(scored['proposal_log_scores']).any()


def test_one_underflowed_candidate_does_not_hide_finite_candidate():
    frozen = artifact([[0, -1000, 0, 0], [0, 0, 0, 0]])
    scored = score_confirmation(np.array([[0, 1, 20, 0]]), ['confirmation'], frozen)
    assert scored['candidate'].tolist() == [1]
    assert scored['unresolved_reason'].tolist() == ['']


def test_test_counts_cannot_rescue_unsupported_proposal_membership():
    frozen = artifact([[0, -1000, 0, 0], [0, -2000, 0, 0]])
    first = score_confirmation(np.array([[0, 1, 0, 0]]), ['a'], frozen)
    second = score_confirmation(np.array([[0, 1, 100000, 100000]]), ['b'], frozen)
    np.testing.assert_array_equal(first['candidate'], second['candidate'])
    np.testing.assert_array_equal(first['proposal_log_scores'], second['proposal_log_scores'])
