"""Independent mass, variance and effective-support invariants for weighted pools."""
import numpy as np
import pytest
from haplodose.clones import Rows, SpecimenModel


def fixture():
    m=SpecimenModel.__new__(SpecimenModel)
    m.rows=Rows(np.arange(16,dtype=float).reshape(8,2)+1,np.arange(8))
    m.L=np.full(8,100.);m.x=np.array([0.,0.,1.,1.,0.,0.,1.,1.]);m.test_fold=None
    m.differences=lambda *args,**kwargs:(args,kwargs)
    return m,np.repeat([0,1],4)


def test_fractional_mass_matching_and_squared_weights():
    m,labels=fixture();w=np.array([.1,.3,.2,.4,.2,.2,.1,.1])
    args,kw=m.matched_differences(labels,[(0,1)],n_strata=2,cell_weights=w)
    # Stratum masses .4/.6 versus .4/.2. Common masses .4 and .2;
    # first group's high-depth weights are divided by three.
    expected=np.array([.1,.3,.2/3,.4/3,.2,.2,.1,.1])
    Y=m.rows.Y
    for group,position in [(0,0),(1,3)]:
        take=labels==group
        np.testing.assert_allclose(args[position][0],expected[take]@Y[take],rtol=32*np.finfo(float).eps)
        np.testing.assert_allclose(kw['count_variance'][group][0],expected[take]**2@Y[take],rtol=32*np.finfo(float).eps)
        np.testing.assert_allclose(args[position+1],[60.],rtol=32*np.finfo(float).eps)


def test_unit_weights_reproduce_legacy_pooling_and_scale_conserves_information():
    m,labels=fixture()
    old=m.matched_differences(labels,[(0,1)],n_strata=2)
    explicit=m.matched_differences(labels,[(0,1)],n_strata=2,cell_weights=np.ones(8))
    for a,b in zip(old[0],explicit[0]):np.testing.assert_array_equal(a,b)
    for key in old[1]:
        for a,b in zip(old[1][key],explicit[1][key]):np.testing.assert_array_equal(a,b)
    w=np.array([.1,.3,.2,.4,.2,.2,.1,.1]);base=m.matched_differences(labels,[(0,1)],n_strata=2,cell_weights=w)
    for scale in [.125,8.]:
        scaled=m.matched_differences(labels,[(0,1)],n_strata=2,cell_weights=scale*w)
        for a,b in zip(base[0],scaled[0]):np.testing.assert_allclose(scale*a,b,rtol=64*np.finfo(float).eps)
        for a,b in zip(base[1]['count_variance'],scaled[1]['count_variance']):np.testing.assert_allclose(scale**2*a,b,rtol=64*np.finfo(float).eps)
        for a,b in zip(base[1]['weight_scale'],scaled[1]['weight_scale']):np.testing.assert_allclose(scale*a,b,rtol=64*np.finfo(float).eps)


def test_zero_confidence_abstains_and_support_is_scale_invariant():
    m,labels=fixture();m.differences=lambda *args,**kwargs:[dict(segments=[{'dummy':True}],distance=1.,testable=True)]
    for scale in [1.,.01,10.]:
        result=m.matched_differences(labels,[(0,1)],n_strata=2,cell_weights=np.full(8,scale),min_effective_cells=4)[0]
        assert result['testable'] and result['effective_cells']==[4.,4.]
    result=m.matched_differences(labels,[(0,1)],cell_weights=np.zeros(8),min_effective_cells=4)[0]
    assert not result['testable'] and result['segments']==[] and result['effective_cells']==[0.,0.]
    for bad in [np.ones(7),np.full(8,-1.),np.full(8,np.nan)]:
        with pytest.raises(ValueError):m.matched_differences(labels,[(0,1)],cell_weights=bad)
