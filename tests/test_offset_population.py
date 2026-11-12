"""The default legacy residual must compare the same retained count population."""
from fractions import Fraction
import math
from types import SimpleNamespace

import numpy as np
import pytest
from haplodose.clones import Rows, SpecimenModel


@pytest.mark.parametrize('chunk', [1, 3])
def test_dropped_population_cannot_change_its_retained_null_offset(chunk):
    Y = np.array([[50, 50], [400, 100]])
    m = SimpleNamespace(rows=Rows(Y, np.arange(2), chunk), base=np.array([.5, .5]), L=Y.sum(axis=1),
        Ytot=Y.sum(axis=0), J=2)
    SpecimenModel._update_offset(m, np.eye(2), [[0., 0.]], [0])
    assert np.array_equal(m.offset, [0., 0.])


def test_partial_soft_mass_is_not_replaced_by_whole_specimen_counts():
    Y = np.array([[60, 40], [400, 100]])
    R = np.array([[.5, .5], [.25, .75]])
    m = SimpleNamespace(rows=Rows(Y, np.arange(2), 1), base=np.array([.5, .5]), L=Y.sum(axis=1),
        Ytot=Y.sum(axis=0), J=2)
    SpecimenModel._update_offset(m, R, [[0., 0.]], [0])
    # Independent retained counts(130,45), exposure175, null expectation(87.5,87.5).
    ratios = [(Fraction(130)+Fraction(1, 2))/(Fraction(175, 2)+Fraction(1, 2)),
              (Fraction(45)+Fraction(1, 2))/(Fraction(175, 2)+Fraction(1, 2))]
    raw = np.array([math.log(float(r)) for r in ratios])
    expected = raw-np.median(raw)
    budget = 64*Y.size*np.finfo(float).eps
    assert np.max(np.abs(m.offset-expected)) <= budget


def test_component_duplication_preserves_retained_exposure_and_half_count():
    Y = np.array([[60, 40], [400, 100]])
    R = np.array([[.5, .5], [.25, .75]])
    def score(weights):
        m = SimpleNamespace(rows=Rows(Y, np.arange(2), 1), base=np.array([.5, .5]), L=Y.sum(axis=1),
            Ytot=Y.sum(axis=0), J=2)
        SpecimenModel._update_offset(m, weights, np.zeros((weights.shape[1], 2)), list(range(weights.shape[1])))
        return m.offset
    np.testing.assert_array_equal(score(R), score(np.repeat(R/2, 2, axis=1)))
