"""Arithmetic of the unpromoted conditional candidate, not the production default.

These six checks stay intact after the full candidate failed primary recovery.
The original full-workflow assertion remains unchanged and is not replaced here.
"""
import importlib.util
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from haplodose.clones import Rows

spec = importlib.util.spec_from_file_location('conditional_offset',
    Path(__file__).parents[1]/'benchmarks/hardening/offset_conditional.py')
conditional = importlib.util.module_from_spec(spec)
spec.loader.exec_module(conditional)


def model(chunk=2):
    # Base ratio1, dose ratio2, depth ratio2 and programme ratio4 give rational
    # cell probabilities. A fourfold exposure ratio gives centred log-depth
    # +/-log2, so the covariate also agrees with the constructor's actual depth.
    Y = np.array([[10880, 1360], [340, 2720], [11520, 720], [612, 2448]])
    L = np.array([12240., 3060.]*2)
    return SimpleNamespace(rows=Rows(Y, np.arange(4), chunk), J=2, base=np.full(2, .5),
        L=L, x=np.log(L)-np.log(L).mean(),
        gamma=np.array([.5, -.5]), U=np.array([[math.log(2)], [-math.log(2)]]*2),
        W=np.array([[1., -1.]]), offset=np.array([5., -5.]))


@pytest.mark.parametrize('keep', [[0, 1], [0]])
def test_represented_nuisance_cannot_create_a_shared_offset(keep):
    m = model()
    R = np.eye(2)[[0, 0, 1, 1]]
    doses = np.array([[0., 0.], [math.log(2), 0.]])[keep]
    conditional.update_offset(m, R, doses, keep)
    # Integer rational counts are exact; only exp/log and pooled accumulation
    # contribute floating error. Bound64 operations per cell/component/bin.
    budget = 64*4*2*2*np.finfo(float).eps
    assert np.max(np.abs(m.offset)) <= budget


def test_kept_population_residual_preserves_extra_counts_and_chunk_invariance():
    R = np.eye(2)[[0, 0, 1, 1]]
    results = []
    # Only the first two cells have kept mass. Their zero-offset pool is(11220,4080).
    expected = math.log((11220+34+.5)/(11220+.5))-math.log((4080-34+.5)/(4080+.5))
    for chunk in (1, 2, 7):
        m = model(chunk)
        m.rows.Y = m.rows.Y+np.array([17, -17])
        conditional.update_offset(m, R, [[0., 0.]], [0])
        results.append(m.offset)
        budget = 64*4*2*2*np.finfo(float).eps
        assert abs((m.offset[0]-m.offset[1])-expected) <= budget
    for result in results[1:]:
        np.testing.assert_allclose(result, results[0], rtol=0., atol=budget)


def test_offset_is_invariant_to_component_labels_and_duplicate_representation():
    R = np.array([[.5, .5]]*4)
    doses = np.array([[0., 0.], [math.log(2), 0.]])
    first = model()
    first.rows.Y = np.array([[11200, 1040], [476, 2584]]*2)
    conditional.update_offset(first, R, doses, [0, 1])
    permutation = model()
    permutation.rows.Y = first.rows.Y
    conditional.update_offset(permutation, R[:, ::-1], doses[::-1], [0, 1])
    duplicated = model()
    duplicated.rows.Y = first.rows.Y
    conditional.update_offset(duplicated, np.repeat(R/2, 2, axis=1), np.repeat(doses, 2, axis=0), list(range(4)))
    budget = 64*4*4*2*np.finfo(float).eps
    assert np.max(np.abs(first.offset)) <= budget
    assert np.max(np.abs(permutation.offset-first.offset)) <= budget
    assert np.max(np.abs(duplicated.offset-first.offset)) <= budget


def test_partial_soft_mass_uses_the_same_observed_count_population():
    m = model()
    R = np.tile([.25, .75], (4, 1))
    conditional.update_offset(m, R, [[0., 0.]], [0])
    # The retained mass weights the observed totals(23352,7248) by1/4;
    # the retained baseline pools are(5610,2040). The original half-count applies
    # once to each of these pools, not to an unweighted whole-specimen numerator.
    raw = np.log((np.array([5838., 1812.])+.5)/(np.array([5610., 2040.])+.5))
    expected = raw-np.median(raw)
    budget = 64*4*2*2*np.finfo(float).eps
    assert np.max(np.abs(m.offset-expected)) <= budget


def test_normalized_mean_is_invariant_to_logit_gauges():
    m = model()
    m.W += 7.
    doses = np.array([[3., 3.], [math.log(2)-2., -2.]])
    conditional.update_offset(m, np.eye(2)[[0, 0, 1, 1]], doses, [0, 1])
    budget = 64*4*2*2*np.finfo(float).eps
    assert np.max(np.abs(m.offset)) <= budget
