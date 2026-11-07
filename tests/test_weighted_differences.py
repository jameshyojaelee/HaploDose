"""Independent conditional-Poisson and weight-unit invariants for pooled contrasts."""
import numpy as np
import pytest

from haplodose.clones import Rows, SpecimenModel
from haplodose.profile import ProfileParams, difference_segments, pooled_lograte


def test_unit_weights_preserve_original_formula():
    Y = np.array([[0., 3., 100.]])
    lr, v = pooled_lograte(Y, [1000.], Y, [1.])
    np.testing.assert_array_equal(lr, np.log((Y + .5) / 1000.))
    np.testing.assert_array_equal(v, 1 / (Y + .5))


def test_independent_poisson_variance_and_weight_units():
    # Independent Poisson means 100 and 200, weighted by 1/4 and 3/4:
    # mean=175; conditional sampling variance=100/16+9*200/16=118.75.
    S, Q, c = np.array([[175.]]), np.array([[118.75]]), np.array([.625])
    lr, v = pooled_lograte(S, [300.], Q, c)
    assert v[0, 0] == (118.75 + .5 * .625**2) / (175 + .5 * .625)**2
    for scale in [.125, 2., 16.]:
        actual = pooled_lograte(scale*S, [scale*300], scale**2*Q, scale*c)
        for expected, observed in zip((lr, v), actual):
            np.testing.assert_allclose(observed, expected, rtol=32*np.finfo(float).eps,
                                       atol=32*np.finfo(float).eps)


def test_matched_pooling_forwards_squared_weights_and_ignores_unassigned_cells():
    Y = np.arange(18, dtype=float).reshape(9, 2) + 1
    model = SpecimenModel.__new__(SpecimenModel)
    model.rows = Rows(Y, np.arange(9))
    model.L = np.full(9, 100.)
    model.x = np.array([0., 0., 1., 0., 0., 0., 1., 1., 100.])
    model.test_fold = None
    model.differences = lambda *args, **kwargs: (args, kwargs)
    args, kwargs = model.matched_differences(np.array([0,0,0,1,1,1,1,1,-1]), [(0,1)], n_strata=2)
    np.testing.assert_array_equal(kwargs['count_variance'][0], Y[:3].sum(axis=0, keepdims=True))
    expected = (4/9)*Y[3:6].sum(axis=0) + .25*Y[6:8].sum(axis=0)
    np.testing.assert_allclose(kwargs['count_variance'][1][0], expected,
                               rtol=32*np.finfo(float).eps, atol=0)
    np.testing.assert_allclose(kwargs['weight_scale'][1], [11/18], rtol=32*np.finfo(float).eps)
    np.testing.assert_allclose([args[1][0],args[4][0]], [300.,300.],
                               rtol=32*np.finfo(float).eps, atol=0)


def test_masked_reset_contrasts_are_invariant_to_pool_weight_scale():
    Ya = np.full((1, 20), 10000.)
    Yb = np.full((1, 20), 10000.)
    Ya[:, :10] *= 1.5
    reset = np.zeros(20, bool); reset[[0,10]] = True
    mask = np.arange(20) % 2 == 1
    params = ProfileParams(min_bins=3)
    base = difference_segments(Ya, [250000.], Yb, [200000.], reset, .6, params, mask=mask,
                               count_variance=(Ya*.6,Yb*.7), weight_scale=([.6],[.7]))
    scaled = difference_segments(Ya*2, [500000.], Yb*.25, [50000.], reset, .6, params, mask=mask,
        count_variance=(Ya*.6*4,Yb*.7/16), weight_scale=([1.2],[.175]))
    assert len(base[0]['segments']) > 0
    assert [(s['first_bin'],s['last_bin']) for s in base[0]['segments']] == [
        (s['first_bin'],s['last_bin']) for s in scaled[0]['segments']]
    for segment in scaled[0]['segments']:
        assert segment['first_bin'] % 2 == segment['last_bin'] % 2 == 1
        assert segment['first_bin']//10 == segment['last_bin']//10
    for a,b in zip(base[0]['segments'],scaled[0]['segments']):
        np.testing.assert_allclose(a['mean_diff'],b['mean_diff'],rtol=128*np.finfo(float).eps,
                                   atol=128*np.finfo(float).eps)


def test_zero_exposure_is_finite_but_untestable_and_invalid_variance_rejected():
    lr,v = pooled_lograte(np.zeros((1,6)), [0.], np.zeros((1,6)), [0.])
    assert np.isfinite(lr).all() and np.isfinite(v).all()
    result = difference_segments(np.zeros(6), 0., np.ones(6), 6., np.r_[True,np.zeros(5,bool)],
                                 .6, ProfileParams(), count_variance=(np.zeros(6),np.ones(6)),
                                 weight_scale=([0.],[1.]))[0]
    assert result['testable'] is False and result['segments'] == []
    for bad in [-1., np.nan, np.inf]:
        with pytest.raises(ValueError):
            pooled_lograte([[1.]], [1.], [[bad]], [1.])
        with pytest.raises(ValueError):
            pooled_lograte([[1.]], [1.], [[1.]], [bad])
