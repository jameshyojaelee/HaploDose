"""Independent partition sums and public averaging/masking contracts."""
from fractions import Fraction
from itertools import product
import math
import sys

import numpy as np
import pandas as pd
import pytest

from haplodose.average_segments import _average, _fit
from haplodose.gene_segments import _block_partition
from haplodose.stacking import stacked_count_profile, stacked_rna_profile


def test_raw_gene_exhaustive_partition_reference():
    groups = [np.array([0., 0.]), np.array([1., 7.]),
              np.array([2., 3.]), np.array([9., 11.])]
    scale = .75/math.log(2.)
    captured = []
    _block_partition(groups, scale, capture=lambda *args: captured.append(_average(*args)))
    actual, record = captured[0]
    paths = []
    penalty = Fraction(3*math.log(sum(map(len, groups))))
    for cuts in product((False, True), repeat=3):
        edges = [0]+[i+1 for i, cut in enumerate(cuts) if cut]+[4]
        objective, profile = penalty*(len(edges)-1), [None]*4
        for lo, hi in zip(edges[:-1], edges[1:]):
            observed = sorted(Fraction(float(v)) for g in groups[lo:hi] for v in g)
            n = len(observed)
            median = (observed[(n-1)//2]+observed[n//2])/2
            objective += 2*sum(abs(v-median) for v in observed)/Fraction(scale)
            profile[lo:hi] = [float(median)]*(hi-lo)
        paths.append((math.exp(-.5*float(objective)), profile))
    total = math.fsum(weight for weight, _ in paths)
    expected = [math.fsum(weight*p[j] for weight, p in paths)/total for j in range(4)]
    assert np.max(abs(actual-expected)) <= record['effect_arithmetic_budget']


def test_uniform_rational_partition_average_and_singleton():
    costs = [[0]*4 for _ in range(4)]
    medians = [[2*(lo+hi) for hi in range(4)] for lo in range(4)]
    actual, record = _average(costs, medians, 1, 1, 0)
    expected = [Fraction(0)]*4
    for cuts in product((False, True), repeat=3):
        edges = [0]+[i+1 for i, cut in enumerate(cuts) if cut]+[4]
        for lo, hi in zip(edges[:-1], edges[1:]):
            for j in range(lo, hi):
                expected[j] += Fraction(medians[lo][hi-1], 16)
    assert np.max(abs(actual-np.array(list(map(float, expected))))) <= record['effect_arithmetic_budget']
    singleton, record = _average([[0]], [[6]], 1, 1, 0)
    np.testing.assert_array_equal(singleton, [3.])


def test_averaging_invariances():
    groups = [np.array([0., 0.]), np.array([1., 7.]),
              np.array([2., 3.]), np.array([9., 11.])]
    def evaluate(values, scale):
        captured = []
        _block_partition(values, scale, capture=lambda *args: captured.append(_average(*args)))
        return captured[0]
    first, a = evaluate(groups, .75)
    shifted, b = evaluate([g+8 for g in groups], .75)
    scaled, c = evaluate([g*4 for g in groups], 3.)
    assert np.max(abs(shifted-first-8)) <= a['effect_arithmetic_budget']+b['effect_arithmetic_budget']
    assert np.max(abs(scaled-first*4)) <= 4*a['effect_arithmetic_budget']+c['effect_arithmetic_budget']


def test_masked_targets_resets_and_no_global_hook():
    values = np.tile([-.25, .125, .375, -.125, .25, -.375, .5, -.5, .625, -.625], 4)
    values += np.repeat([0., .5, 10., 10.5], 10)
    bins, blocks = np.repeat(np.arange(4), 10), ['p', 'p', 'q', 'q']
    hook = sys.getprofile()
    _, standalone, _, first = _fit(values, bins, blocks)
    assert sys.getprofile() is hook
    changed = values.copy(); changed[0] += 64
    _, _, _, second = _fit(changed, bins, blocks)
    assert first['average_predictions'][0] == second['average_predictions'][0]
    assert standalone[:2].max() < standalone[2:].min()
    # Other chromosome values cannot change this chromosome's full expert.
    changed = values.copy(); changed[bins >= 2] += 64
    # Translation preserves the common within-bin residual scale in this fixture.
    _, other, _, _ = _fit(changed, bins, blocks)
    np.testing.assert_array_equal(standalone[:2], other[:2])


def test_public_count_forwarding_and_default_preservation():
    counts = np.array([[2, 8, 3, 7, 20, 60, 30, 50], [5]*8])
    genes = [f'g{i}' for i in range(8)]
    bins, blocks = np.repeat(np.arange(4), 2), ['p', 'p', 'q', 'q']
    args = (counts, ['q', 'r'], genes, pd.DataFrame(dict(gene=genes, bin=bins)),
            blocks, ['q'], ['r'])
    old, effects, old_record = stacked_count_profile(*args, counts_kind='raw', gene_segments=True)
    same, same_effects, same_record = stacked_count_profile(*args, counts_kind='raw',
        gene_segments=True, average_segments=False)
    np.testing.assert_array_equal(old, same)
    np.testing.assert_array_equal(effects, same_effects)
    assert old_record == same_record
    new, new_effects, record = stacked_count_profile(*args, counts_kind='raw',
        gene_segments=True, average_segments=True)
    expected, selection = stacked_rna_profile(effects, bins, blocks,
        gene_segments=True, average_segments=True)
    np.testing.assert_array_equal(new_effects, effects)
    np.testing.assert_array_equal(new, expected)
    assert all(record[k] == v for k, v in selection.items())
    assert record['calibration_status'] == 'unvalidated'
    assert record['units'] == 'natural_log_relative_rna_rate'


@pytest.mark.parametrize('flag', ['False', 0, 1, None])
def test_nonboolean_average_mode_rejected(flag):
    with pytest.raises(ValueError, match='average_segments must be boolean'):
        stacked_rna_profile(np.ones(4), [0, 0, 1, 1], ['p', 'p'], average_segments=flag)


def test_average_mode_requires_gene_segments_and_preserves_abstention():
    with pytest.raises(ValueError, match='requires gene_segments=True'):
        stacked_rna_profile(np.ones(4), [0, 0, 1, 1], ['p', 'p'], average_segments=True)
    profile, record = stacked_rna_profile([1., 2.], [0, 1], ['p', 'p'],
        gene_segments=True, average_segments=True)
    np.testing.assert_array_equal(profile, [1., 2.])
    assert not record['selected'] and record['reason'] == 'old_expert_support_unavailable'
