"""Independent numerical checks for the experimental conditional-mixture statistic."""
import importlib.util
from pathlib import Path
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp
from scipy.stats import nbinom


def load():
    path=Path(__file__).resolve().parents[1]/'benchmarks/hardening/membership_predictive_search.py'
    spec=importlib.util.spec_from_file_location('membership_statistic',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.conditional_gain


def test_conditional_gain_matches_independent_full_pmf_and_label_permutation():
    gain=load()
    O=np.array([[2.,6.],[1.,0.],[8.,9.],[3.,4.]])
    E=np.array([[2.,3.],[.5,1.],[4.,6.],[2.,4.]])
    r=np.broadcast_to([2.,4.],O.shape)
    q=np.array([0.,.25,.7,1.]);rates=np.log([.7,1.4])
    ll=lambda eta: nbinom.logpmf(O,r,r/(r+E*np.exp(eta))).sum(axis=1)
    null=minimize_scalar(lambda x:-ll(x).sum(),bounds=(-5.,5.),method='bounded',options={'xatol':1e-12})
    with np.errstate(divide='ignore'):
        alt=logsumexp(np.stack([ll(rates[0])+np.log1p(-q),ll(rates[1])+np.log(q)]),axis=0).sum()
    expected=alt+null.fun
    # Full PMF subtraction has a reduction-sized float64 error; the optimizer
    # is independently bracketed and its objective is locally quadratic.
    allowance=256*O.size*np.finfo(float).eps*(abs(alt)+abs(null.fun))
    assert abs(gain(O,E,r,q,rates)-expected)<=allowance
    assert abs(gain(O,E,r,1-q,rates[::-1])-expected)<=allowance
    assert abs(gain(O[::-1],E[::-1],r[::-1],q[::-1],rates)-expected)<=allowance


def test_equal_components_ignore_membership_uncertainty():
    gain=load();O=np.full((8,3),4.);E=np.full_like(O,4.);r=np.full_like(O,3.)
    assert abs(gain(O,E,r,np.linspace(0,1,8),np.zeros(2)))<=128*O.size*np.finfo(float).eps


def test_matched_hard_control_shares_null_and_component_likelihoods():
    gain=load()
    O=np.array([[2.,6.],[1.,0.],[8.,9.],[3.,4.]])
    E=np.array([[2.,3.],[.5,1.],[4.,6.],[2.,4.]])
    r=np.broadcast_to([2.,4.],O.shape);q=np.array([.1,.25,.7,.9]);rates=np.log([.7,1.4])
    soft,hard=gain(O,E,r,q,rates,return_hard=True)
    assert soft==gain(O,E,r,q,rates)
    # Certain memberships in the same mixture reduce exactly to the hard rule;
    # it shares all nuisance parameters and the profiled null denominator.
    assert hard==gain(O,E,r,(q>.5).astype(float),rates)
