"""Independent exact partition scores and optional profile contracts."""
from fractions import Fraction
from itertools import product
import math

import numpy as np
import pandas as pd
import pytest

from haplodose.gene_segments import _block_partition, _partition_profile
from haplodose.stacking import stacked_count_profile, stacked_rna_profile


def test_exact_gene_partition_matches_every_partition():
    groups = [np.array([0., 0.]), np.array([1., 7.]),
              np.array([2., 3.]), np.array([9., 11.])]
    scale = .75/math.log(2.)
    fitted, record = _block_partition(groups, scale)
    candidates = []
    for cuts in product((False, True), repeat=3):
        edges = [0]+[i+1 for i, cut in enumerate(cuts) if cut]+[4]
        loss = Fraction(record['penalty'])*(len(edges)-2)
        profile = np.empty(4)
        for lo, hi in zip(edges[:-1], edges[1:]):
            observed = sorted(Fraction(float(v)) for g in groups[lo:hi] for v in g)
            n = len(observed)
            median = (observed[(n-1)//2]+observed[n//2])/2
            loss += 2*sum(abs(v-median) for v in observed)/Fraction(scale)
            profile[lo:hi] = float(median)
        candidates.append((loss, profile))
    loss, expected = min(candidates, key=lambda item: item[0])
    assert Fraction(int(record['exact_objective_numerator']),
                    int(record['exact_objective_denominator'])) == loss
    np.testing.assert_array_equal(fitted, expected)
    translated, other = _block_partition([g+2**20 for g in groups], scale)
    np.testing.assert_array_equal(translated, fitted+2**20)
    assert other['exact_objective_numerator'] == record['exact_objective_numerator']
    multiplied, other = _block_partition([g*8 for g in groups], scale*8)
    np.testing.assert_array_equal(multiplied, fitted*8)
    assert Fraction(int(other['exact_objective_numerator']),
                    int(other['exact_objective_denominator'])) == loss


def test_masking_resets_and_unidentified_scale():
    values = np.concatenate([np.arange(4.)+shift for shift in (0., 1., 10., 11.)])
    bins, blocks = np.repeat(np.arange(4), 4), np.array(['p', 'p', 'q', 'q'])
    hold = np.tile([True, False, False, False], 4)
    first, record = _partition_profile(values[~hold], bins[~hold], blocks)
    changed = values.copy()
    changed[hold] += 100.
    second, other = _partition_profile(changed[~hold], bins[~hold], blocks)
    np.testing.assert_array_equal(first, second)
    assert record == other
    assert first[:2].max() < first[2:].min()
    fitted, record = _partition_profile(np.ones(16), bins, blocks)
    assert fitted is None and record['status'] == 'positive_working_scale_not_identified'


def test_count_api_forwards_optional_mode_and_preserves_effects():
    counts = np.array([[2, 8, 3, 7, 20, 60, 30, 50], [5]*8])
    genes = [f'g{i}' for i in range(8)]
    bins, blocks = np.repeat(np.arange(4), 2), ['p', 'p', 'q', 'q']
    args = (counts, ['q', 'r'], genes, pd.DataFrame(dict(gene=genes, bin=bins)),
            blocks, ['q'], ['r'])
    original, effects, previous = stacked_count_profile(*args, counts_kind='raw')
    same, same_effects, same_record = stacked_count_profile(*args, counts_kind='raw', gene_segments=False)
    np.testing.assert_array_equal(original, same)
    np.testing.assert_array_equal(effects, same_effects)
    assert previous == same_record
    fitted, other_effects, record = stacked_count_profile(*args, counts_kind='raw', gene_segments=True)
    expected, selection = stacked_rna_profile(effects, bins, blocks, gene_segments=True)
    np.testing.assert_array_equal(other_effects, effects)
    np.testing.assert_array_equal(fitted, expected)
    for key, value in selection.items():
        assert record[key] == value
    assert record['units'] == 'natural_log_relative_rna_rate'
    assert record['calibration_status'] == 'unvalidated'


@pytest.mark.parametrize('flag', ['False', 0, 1, None])
def test_nonboolean_mode_rejected(flag):
    with pytest.raises(ValueError, match='boolean'):
        stacked_rna_profile(np.ones(4), [0, 0, 1, 1], ['p', 'p'], gene_segments=flag)


def test_optional_mode_abstains_and_rejects_invalid_identity():
    fitted, record = stacked_rna_profile([1., 2.], [0, 1], ['p', 'p'], gene_segments=True)
    np.testing.assert_array_equal(fitted, [1., 2.])
    assert not record['selected']
    with pytest.raises(ValueError):
        stacked_rna_profile([1., 2.], [0., .5], ['p'], gene_segments=True)
