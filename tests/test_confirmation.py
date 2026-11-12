"""Role and feature leakage invariants for frozen candidate evaluation."""
import numpy as np
import pytest

from haplodose.confirmation import Contrast, FrozenDiscovery, assign_roles, score_confirmation


def artifact():
    return FrozenDiscovery(log_rates=np.log([[.4, .1, .4, .1], [.1, .4, .1, .4]]),
                           phi=np.full(4, .1), proposal_mask=[True, True, False, False],
                           test_mask=[False, False, True, True],
                           contrasts=(Contrast(0, 1, (2,), 1, 'region'),), discovery_cell_ids=('train',),
                           programme_loadings=np.array([[.1, -.1, .2, -.2]]), programme_precision=[10.])


def test_roles_are_id_stable_and_order_invariant():
    ids = np.array([f'cell{i}' for i in range(200)])
    kwargs = dict(discovery_fraction=.5, seed=42, min_discovery=20, min_confirmation=20)
    a, b = assign_roles(ids, **kwargs), assign_roles(ids[::-1], **kwargs)
    np.testing.assert_array_equal(a['discovery'], b['discovery'][::-1])
    extended = assign_roles(np.concatenate([ids, ['additional']]), **kwargs)
    np.testing.assert_array_equal(a['discovery'], extended['discovery'][:-1])
    assert not (a['discovery'] & a['confirmation']).any()
    with pytest.raises(ValueError, match='minimum cell counts'):
        assign_roles(ids[:2], **kwargs)


def test_confirmation_counts_cannot_change_membership_or_nuisance_scores():
    frozen = artifact()
    Y = np.array([[80, 20, 80, 20], [20, 80, 20, 80], [70, 30, 70, 30], [30, 70, 30, 70]])
    ids = ['a', 'b', 'c', 'd']
    first = score_confirmation(Y, ids, frozen)
    changed = Y.copy()
    changed[:, 2:] = np.array([[1, 999], [999, 1], [1, 999], [999, 1]])
    second = score_confirmation(changed, ids, frozen)
    np.testing.assert_array_equal(first['candidate'], second['candidate'])
    np.testing.assert_array_equal(first['proposal_log_scores'], second['proposal_log_scores'])
    np.testing.assert_array_equal(first['programme_scores'], second['programme_scores'])
    assert first['contrasts'][0]['log_odds_contrast'] > 0
    assert second['contrasts'][0]['log_odds_contrast'] < 0
    assert not first['evidence']['specimen_error_control']


def test_frozen_scoring_is_row_order_invariant_and_inputs_are_copied():
    frozen = artifact()
    Y = np.array([[80, 20, 80, 20], [20, 80, 20, 80]])
    a = score_confirmation(Y, ['a', 'b'], frozen)
    b = score_confirmation(Y[::-1], ['b', 'a'], frozen)
    np.testing.assert_array_equal(a['candidate'], b['candidate'][::-1])
    np.testing.assert_array_equal(a['proposal_log_scores'], b['proposal_log_scores'][::-1])
    assert a['contrasts'] == b['contrasts']
    with pytest.raises(ValueError):
        frozen.log_rates[0, 0] = 1.


def test_discovery_rows_and_overlapping_feature_roles_are_rejected():
    with pytest.raises(ValueError, match='discovery cells'):
        score_confirmation(np.ones((1, 4)), ['train'], artifact())
    with pytest.raises(ValueError, match='disjoint'):
        FrozenDiscovery(log_rates=np.zeros((2, 4)), phi=np.ones(4), proposal_mask=[1, 1, 0, 0],
                        test_mask=[0, 1, 1, 1], contrasts=(), discovery_cell_ids=('train',))


def test_zero_proposal_depth_abstains_even_with_test_signal():
    result = score_confirmation(np.array([[0, 0, 10000, 1]]), ['heldout'], artifact())
    assert result['candidate'].tolist() == [-1]
    assert not result['contrasts'][0]['testable']


@pytest.mark.parametrize('permutation', [[0, 1, 2], [2, 1, 0], [1, 0, 2]])
def test_identical_top_proposal_profiles_abstain_under_candidate_permutation(permutation):
    # The test profiles disagree, but cannot supply proposal membership evidence.
    rates = np.array([[.4, .1, .4, .1], [.4, .1, .1, .4], [.1, .4, .1, .4]])
    frozen = FrozenDiscovery(log_rates=np.log(rates[permutation]), phi=np.full(4, .1),
        proposal_mask=[True, True, False, False], test_mask=[False, False, True, True],
        contrasts=(), discovery_cell_ids=('train',))
    for test_counts in ([80, 20], [20, 80], [0, 0]):
        result = score_confirmation([[80, 20, *test_counts]], ['heldout'], frozen)
        assert result['candidate'].tolist() == [-1]
        assert result['unresolved_reason'].tolist() == ['tied_candidate_scores']


def test_symmetric_evidence_ties_do_not_enter_either_contrast_pool():
    frozen = FrozenDiscovery(log_rates=np.log([[.4, .1, .4, .1], [.1, .4, .1, .4]]),
        phi=np.full(4, .1), proposal_mask=[True, True, False, False],
        test_mask=[False, False, True, True],
        contrasts=(Contrast(0, 1, (2,), 1, 'region'),), discovery_cell_ids=('train',))
    Y = [[80, 20, 8, 2], [20, 80, 2, 8], [10, 10, 10000, 1]]
    result = score_confirmation(Y, ['a', 'b', 'tie'], frozen)
    assert result['candidate'].tolist() == [0, 1, -1]
    assert result['contrasts'][0]['n_cells'] == [1, 1]
    assert result['contrasts'][0]['test_counts'] == [[8., 2.], [2., 8.]]
    np.testing.assert_array_equal(result['proposal_log_scores'][2, 0], result['proposal_log_scores'][2, 1])


def test_single_finite_candidate_is_not_a_tie():
    frozen = FrozenDiscovery(log_rates=np.zeros((1, 4)), phi=np.ones(4),
        proposal_mask=[True, True, False, False], test_mask=[False, False, True, True],
        contrasts=(), discovery_cell_ids=('train',))
    result = score_confirmation([[1, 1, 0, 0]], ['a'], frozen)
    assert result['candidate'].tolist() == [0]
    assert result['unresolved_reason'].tolist() == ['']
