"""Independent conditional moments and exact pure-group emission regression."""
import argparse,importlib.util,itertools,json
from pathlib import Path
import numpy as np
from haplodose.provenance import start_run,finish_run
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
path=Path(__file__).with_name('mixture_confirmation.py');s=importlib.util.spec_from_file_location('mixture',path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
rec=start_run(a.out,{'runner':__file__,'model':path},{},{'scope':'Conditional on fixed membership probabilities'},mode='mixture_moments')
try:
 q=np.array([.1,.3,.7,.9]);L=np.array([2.,3.,5.,7.]);W=np.array([[1,0],[.5,0],[0,.4],[0,1.]])
 purity,h=m.allocation_moments(q,L,W);delta=np.array([-.6*np.log(1.5),0.,.6*np.log(1.5)])
 errors=[]
 for effect in delta:
    values=[];prob=[]
    for state in itertools.product((0,1),repeat=4):
        z=np.array(state);prob.append(np.prod(np.where(z,q,1-q)))
        values.append(W.T@(L*np.exp(effect*z)))
    values=np.asarray(values);prob=np.asarray(prob);expected=prob@values
    variance=prob@((values-expected)**2)
    exposure=W.T@L
    predicted=exposure*(1+purity*np.expm1(effect))
    predicted_var=exposure**2*h*np.expm1(effect)**2
    np.testing.assert_allclose(expected,predicted,rtol=64*np.finfo(float).eps,atol=64*np.finfo(float).eps)
    np.testing.assert_allclose(variance,predicted_var,rtol=128*np.finfo(float).eps,atol=128*np.finfo(float).eps)
    mean,var=m.observable_states(purity[:1],purity[1:],h[:1],h[1:],[effect])
    np.testing.assert_allclose(mean[0,0],np.log(expected[0]/exposure[0])-np.log(expected[1]/exposure[1]),atol=64*np.finfo(float).eps)
    np.testing.assert_allclose(var[0,0],np.sum(variance/expected**2),atol=64*np.finfo(float).eps)
    errors.append(float(np.max(np.abs(variance-predicted_var))))
 for factor in (.125,8.):
    for before,after in zip((purity,h),m.allocation_moments(q,L,factor*W)):
        np.testing.assert_array_equal(before,after)
 mean,var=m.observable_states([1],[0],[0],[0],delta)
 np.testing.assert_array_equal(mean[0],delta)
 np.testing.assert_array_equal(var,np.zeros((1,3)))
 (a.out/'summary.json').write_text(json.dumps(dict(exhaustive_states_per_effect=16,effects=3,moment_max_errors=errors,pure_membership_identity=True,weight_scale_invariance=True),indent=2))
 finish_run(a.out,rec,backend='cpu')
except Exception as error:
 finish_run(a.out,rec,error=error);raise
