"""Descriptive RNA prediction has known optima and respects genomic boundaries."""
import statistics

import numpy as np
import pytest

from haplodose.smoothing import _running_median, predictive_rna_profile


def test_noiseless_contrasts_and_ties_have_known_optima():
    bins = np.repeat(np.arange(7), 4)
    truth = np.array([0., 2., 0., 2., 0., 2., 0.])
    fitted, decision = predictive_rna_profile(np.repeat(truth, 4), bins, ['p']*7)
    np.testing.assert_array_equal(fitted, truth)
    assert decision['selected_width'] == 1
    assert decision['mean_prediction_error'][0] == 0.
    fitted, decision = predictive_rna_profile(np.ones(len(bins)), bins, ['p']*7)
    assert decision['selected_width'] == 1
    assert decision['numerically_tied_widths'] == decision['widths']
    np.testing.assert_array_equal(fitted, np.ones(7))


def test_rolling_medians_match_scalar_reference_and_stop_at_gaps():
    values = [2., 9., 1., 4., 100., 110., 101., 109.]
    expected = []
    for block in [values[:4], values[4:]]:
        expected.extend(statistics.median(block[max(0, j-1):min(4, j+2)]) for j in range(4))
    np.testing.assert_array_equal(_running_median(values, np.repeat(['a', 'b'], 4), 3), expected)
    fitted, _ = predictive_rna_profile(np.repeat([0.]*4+[100.]*4, 3),
                                      np.repeat(np.arange(8), 3), np.repeat(['a', 'b'], 4))
    np.testing.assert_array_equal(fitted, [0.]*4+[100.]*4)


def test_singletons_abstain_and_remain_available_as_training_neighbours():
    fitted, record = predictive_rna_profile([2., 7., 4.], np.arange(3), ['p']*3)
    np.testing.assert_array_equal(fitted, [2., 7., 4.])
    assert record['selected_width'] is None
    assert record['status'] == 'insufficient_masked_gene_support'
    assert record['calibration_status'] == 'unvalidated'
    # In fold zero, the singleton middle bin remains observed. The width-three
    # predictions are the medians of [0,10] and [10,20], hence errors 5 and 5.
    _, record = predictive_rna_profile([0., 0., 10., 20., 20.],
                                      np.array([0, 0, 1, 2, 2]), ['p']*3)
    assert record['fold_error_sums'][0] == [0., 10.]
    assert sum(record['held_out_genes']) == 4


@pytest.mark.parametrize('values,bins,blocks', [
    ([0., np.nan], [0, 1], ['p', 'p']),
    ([0., 1.], [0, 2], ['p', 'p']),
    ([0., 1.], [0., 1.], ['p', 'p']),
    ([0., 1., 2.], [0, 1, 2], ['p', 'q', 'p']),
    ([0., 1.], [0, 1], ['p', None]),
])
def test_invalid_observations_and_coordinates_fail(values, bins, blocks):
    with pytest.raises(ValueError):
        predictive_rna_profile(values, np.asarray(bins), blocks)
