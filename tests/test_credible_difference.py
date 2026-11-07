"""Difference evidence survives uncertain flanks without accepting weak bins."""
import numpy as np

from haplodose.profile import ProfileParams, credible_difference_runs


def posterior(support):
    support = np.asarray(support)
    return np.column_stack([support, 1 - support, np.zeros(len(support))])


def test_uncertain_flanks_do_not_erase_a_supported_interior():
    params = ProfileParams(min_bins=5, posterior=.95)
    support = [.99] * 9 + [.826, .534]
    reset = np.r_[True, np.zeros(10, bool)]
    assert np.mean(support) < params.posterior
    core = credible_difference_runs(np.zeros(11, int), posterior(support), reset, params)
    assert core == [(0, 8, 0)]
    # Adding another uncertain flank changes no retained evidence or threshold.
    extended = credible_difference_runs(np.zeros(12, int), posterior([.60, *support]),
                                        np.r_[True, np.zeros(11, bool)], params)
    assert extended == [(1, 9, 0)]


def test_uncertain_interior_and_genomic_gaps_cannot_be_bridged():
    params = ProfileParams(min_bins=5, posterior=.95)
    support = [.999] * 4 + [.90] + [.999] * 4
    # High average confidence cannot license two individually too-short pieces.
    assert np.mean(support) > params.posterior
    assert credible_difference_runs(np.zeros(9, int), posterior(support),
                                    np.r_[True, np.zeros(8, bool)], params) == []
    reset = np.zeros(9, bool)
    reset[[0, 4]] = True
    assert credible_difference_runs(np.zeros(9, int), posterior([.999] * 9), reset, params) == [(4, 8, 0)]


def test_neutral_or_low_support_bins_do_not_become_differences():
    params = ProfileParams(min_bins=5, posterior=.95)
    reset = np.r_[True, np.zeros(9, bool)]
    assert credible_difference_runs(np.ones(10, int), posterior([.001] * 10), reset, params) == []
    assert credible_difference_runs(np.zeros(10, int), posterior([.949] * 10), reset, params) == []
