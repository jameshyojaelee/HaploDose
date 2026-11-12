"""Independent numerical optima and masking contracts for RNA stacking."""
from fractions import Fraction
from itertools import product
from types import SimpleNamespace

import numpy as np
import pytest

from haplodose import stacking
from haplodose.stacking import _convex_l1, _partition, _training_predictions, stacked_rna_profile


@pytest.mark.parametrize('solver', ['primal', 'dual'])
def test_simplex_l1_matches_exact_breakpoints_and_feasible_dual(solver):
    P = np.array([[0., 2.], [0., 2.], [2., 4.]])
    y = np.array([1., 2., 3.])
    candidates = {Fraction(0), Fraction(1)}
    for a, b, t in zip(P[:, 0], P[:, 1], y):
        candidates.add(min(Fraction(1), max(Fraction(0),
            (Fraction(float(t))-Fraction(float(a)))/(Fraction(float(b))-Fraction(float(a))))))
    losses = [sum(abs((1-z)*Fraction(float(a))+z*Fraction(float(b))-Fraction(float(t)))
                  for a, b, t in zip(P[:, 0], P[:, 1], y))/len(y) for z in candidates]
    assert min(losses) == Fraction(1, 3)
    weights, record = _convex_l1(P, y, solver=solver)
    np.testing.assert_array_equal(weights, [.5, .5])
    assert abs(record['mean_absolute_error']-1/3) <= record['rounding_budget']
    assert abs(record['feasible_dual_lower']-1/3) <= record['rounding_budget']
    permuted, _ = _convex_l1(P[::-1], y[::-1], solver=solver)
    np.testing.assert_array_equal(weights, permuted)


def test_failed_dual_falls_back_to_the_same_primal_fit(monkeypatch):
    P = np.array([[0., 2.], [0., 2.], [2., 4.]])
    y = np.array([1., 2., 3.])
    reference, reference_record = _convex_l1(P, y)
    failed = SimpleNamespace(success=False, status=4)
    monkeypatch.setattr(stacking, '_l1_dual', lambda *_: (failed, None, None))
    weights, record = _convex_l1(P, y, solver='dual')
    np.testing.assert_array_equal(weights, reference)
    assert record == dict(reference_record, solver='primal_fallback', dual_solver_status=4)


def test_dual_objective_can_be_identified_with_nonunique_expert_weights():
    P = np.array([[0., 0., 2.], [0., 0., 2.], [2., 2., 4.]])
    y = np.array([1., 2., 3.])
    for solver in ('primal', 'dual'):
        w, record = _convex_l1(P, y, solver=solver)
        np.testing.assert_array_equal(P @ w, [1., 1., 3.])
        assert abs(record['feasible_dual_lower']-1/3) <= record['rounding_budget']


def test_segment_score_matches_exhaustive_partitions_and_stops_at_gaps():
    # Every partition scored independently with direct residuals, no prefix sums.
    y = np.array([0., .2, 4., 3.8, -9., -8.9])
    var = np.array([.1, .2, .15, .1, .1, .2])
    reset = np.array([True, False, False, False, True, False])
    fitted, _ = _partition(y, var, reset)
    for begin, end in [(0, 4), (4, 6)]:
        expected = []
        for cuts in product([False, True], repeat=end-begin-1):
            edges = [begin]+[begin+i+1 for i, cut in enumerate(cuts) if cut]+[end]
            loss, candidate = 3*np.log(end-begin)*(len(edges)-2), np.empty(end-begin)
            for lo, hi in zip(edges[:-1], edges[1:]):
                weight = 1/var[lo:hi]
                mean = np.average(y[lo:hi], weights=weight)
                loss += np.sum(weight*(y[lo:hi]-mean)**2)
                candidate[lo-begin:hi-begin] = mean
            expected.append((loss, candidate))
        reference = min(expected, key=lambda item: item[0])[1]
        np.testing.assert_array_equal(fitted[begin:end], reference)


def test_base_predictions_exclude_targets_and_respect_genomic_blocks():
    bins = np.repeat(np.arange(7), 4)
    blocks = np.array(['a']*4+['b']*3)
    values = np.repeat([0., .2, 1., 1.3, 10., 10.2, 11.], 4)
    folds = np.tile(np.arange(4), 7)
    first, segments, hold, _ = _training_predictions(values, bins, blocks, [1, 3], folds, 0)
    changed = values.copy()
    changed[hold] += 100.
    second, other, _, _ = _training_predictions(changed, bins, blocks, [1, 3], folds, 0)
    np.testing.assert_array_equal(first, second)
    np.testing.assert_array_equal(segments, other)
    assert first[:, :4].max() < first[:, 4:].min()


def test_constant_fallback_singleton_abstention_and_units():
    fitted, record = stacked_rna_profile(np.ones(28), np.repeat(np.arange(7), 4), ['a']*7)
    np.testing.assert_array_equal(fitted, np.ones(7))
    assert not record['selected']
    fitted, record = stacked_rna_profile([1., 2., 3.], np.arange(3), ['a']*3)
    np.testing.assert_array_equal(fitted, [1., 2., 3.])
    assert not record['selected'] and record['status'] == 'insufficient_masked_gene_support'
    assert record['units'] == 'natural_log_relative_rna_rate'
    assert record['calibration_status'] == 'unvalidated'


@pytest.mark.parametrize('values,bins,blocks', [
    ([0., 1.], [0., 1.], ['a', 'a']),
    ([0., 1.], [False, True], ['a', 'a']),
    ([0., 1.], [0, 2], ['a', 'a']),
    ([0., np.nan], [0, 1], ['a', 'a']),
    ([0., 1., 2.], [0, 1, 2], ['a', 'b', 'a']),
    ([0., 1.], [0, 1], ['a', None]),
])
def test_invalid_bin_identity_is_rejected_before_conversion(values, bins, blocks):
    with pytest.raises(ValueError):
        stacked_rna_profile(values, bins, blocks)
