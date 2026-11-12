"""Scientific invariants for exact arm aggregation, independent of input dtype."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


def scorer():
    path = Path(__file__).resolve().parents[1] / 'benchmarks/hardening/arm_metrics.py'
    spec = importlib.util.spec_from_file_location('arm_scoring_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.aggregate_metrics


@pytest.mark.parametrize('aggregation', ['mean', 'median'])
@pytest.mark.parametrize('dtype', [np.float32, np.float64])
def test_equal_arm_summaries_have_exact_tied_ranks(aggregation, dtype):
    sizes = [3, 7, 5, 9, 11]
    values = np.repeat(np.array([.1, .1, .2, .3, .4], dtype=dtype), sizes)
    truth = np.repeat(np.array([.1, .2, .2, .3, .5], dtype=dtype), sizes)
    groups = np.repeat(['a', 'b', 'c', 'd', 'e'], sizes)
    score = scorer()(values, truth, groups, aggregation)
    # Ranks centered at 3 have covariance 35/4 and both sums of squares 19/2.
    assert score['spearman'] == 35 / 38
    permutation = np.random.default_rng(2026111109).permutation(len(groups))
    assert scorer()(values[permutation], truth[permutation], groups[permutation], aggregation) == score


@pytest.mark.parametrize('aggregation', ['mean', 'median'])
def test_float32_promotion_preserves_the_same_observations(aggregation):
    values = np.array([1., np.nextafter(np.float32(1), np.float32(2)), 2., 2., 3., 3., 4., 4.], dtype=np.float32)
    truth = np.repeat([4., 1., 3., 2.], 2)
    groups = np.repeat(['a', 'b', 'c', 'd'], 2)
    assert scorer()(values, truth, groups, aggregation) == scorer()(values.astype(np.float64), truth, groups, aggregation)


def test_fixed_support_is_not_silently_reduced():
    with pytest.raises(ValueError, match='Finite'):
        scorer()([1., 2., np.nan, 4.], [1., 2., 3., 4.], ['a', 'b', 'c', 'd'], 'mean')
    with pytest.raises(ValueError, match='Aligned'):
        scorer()([1.], [1., 2.], ['a'], 'mean')


def test_ucan_wrapper_from_an_independent_directory():
    path = Path(__file__).resolve().parents[1] / 'benchmarks/hardening/ucan.py'
    spec = importlib.util.spec_from_file_location('arm_wrapper_test', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    values, truth, groups = [1., 2., 3., 4.], [4., 3., 2., 1.], ['a', 'b', 'c', 'd']
    assert module.arm_metrics(values, truth, groups, 'median') == dict(n=4, pearson=-1., spearman=-1.)
