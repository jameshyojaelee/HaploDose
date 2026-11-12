"""Count-unit and genomic-boundary invariants for descriptive regional effects."""
from fractions import Fraction
import statistics

import numpy as np
import pytest

from haplodose.bulk import _gene_expectation, run_bulk
from haplodose.smoothing import rna_profile_on_regions


def test_gene_expectation_conservation_and_reference_permutation():
    reference = np.array([3., 7., 0., 2., 6.])
    bins = np.array([0, 0, 0, 1, 1])
    expected = _gene_expectation([20., 16.], reference, bins)
    exact = [Fraction(20)*Fraction(3, 10), Fraction(20)*Fraction(7, 10),
             Fraction(0), Fraction(16)*Fraction(2, 8), Fraction(16)*Fraction(6, 8)]
    np.testing.assert_array_equal(expected, list(map(float, exact)))
    np.testing.assert_array_equal(np.bincount(bins, weights=expected), [20., 16.])
    permutation = np.array([4, 0, 2, 1, 3])
    np.testing.assert_array_equal(_gene_expectation([20., 16.], reference[permutation], bins[permutation]),
                                  expected[permutation])
    with pytest.raises(ValueError, match='positive reference support'):
        _gene_expectation([20., 16.], [3., 7., 0., 0., 0.], bins)


def test_region_gene_median_is_not_median_of_bin_medians():
    # Three genes vs one makes the two possible aggregation rules disagree.
    values = np.array([0., 2., 4., 20., 10., 12.])
    bins = np.array([0, 0, 0, 1, 2, 3])
    expected = [statistics.median([0., 2., 4., 20.])]*2+[11., 11.]
    result, metadata = rna_profile_on_regions(values, bins, np.array([2, 2, 2, 2]), [0, 0, 1, 1])
    np.testing.assert_array_equal(result, expected)
    assert metadata['n_regions'] == 2
    permutation = np.array([3, 0, 5, 2, 1, 4])
    permuted, _ = rna_profile_on_regions(values[permutation], bins[permutation],
                                        np.array([9, 9, 9, 9]), [0, 0, 1, 1])
    np.testing.assert_array_equal(result, permuted)
    split, _ = rna_profile_on_regions(values, bins, np.array([0, 1, 1, 2]), [0, 0, 1, 1])
    np.testing.assert_array_equal(split, [2., 20., 10., 12.])


@pytest.mark.parametrize('values,bins,path,blocks', [
    ([0., np.nan], [0, 1], [0, 0], [0, 0]),
    ([0., 1.], [0, 2], [0, 0], [0, 0]),
    ([0., 1.], [0., 1.], [0, 0], [0, 0]),
    ([0., 1.], [0, 1], [-1, 0], [0, 0]),
    ([0., 1., 2.], [0, 1, 2], [0, 0, 0], [0, 1, 0]),
])
def test_unsupported_region_inputs_fail(values, bins, path, blocks):
    with pytest.raises(ValueError):
        rna_profile_on_regions(values, bins, path, blocks)


def test_median_pseudocount_expectations_rejected_before_output_creation(tmp_path):
    out = tmp_path/'output'
    with pytest.raises(ValueError, match="requires bin_stat='sum'"):
        run_bulk('unused', 'unused', 'unused', out, counts_kind='raw', genome_build='hg38',
                 coordinate_system='zero-based-half-open', gene_namespace='custom',
                 bin_stat='median', rna_profile_smoothing='native-regions')
    assert not out.exists()
