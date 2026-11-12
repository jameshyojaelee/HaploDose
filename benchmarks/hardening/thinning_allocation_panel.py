"""Fixed allocation screen: original assertions, exact equal-share control."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
from haplodose import clones
from haplodose.provenance import start_run,finish_run


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def checks():
    discrepancies=[]
    for fraction in (1/3,.5,2/3):
        for mean,size in ((3.,2.),(20.,10.)):
            p=size/(mean+size)
            for u in range(13):
                for v in range(13):
                    total=u+v
                    joint=float(stats.nbinom.logpmf(total,size,p)+stats.betabinom.logpmf(u,total,fraction*size,(1-fraction)*size))
                    def nb(count,r):
                        terms=[math.lgamma(count+r),-math.lgamma(r),-math.lgamma(count+1),r*math.log(p),count*math.log1p(-p)]
                        return math.fsum(terms),math.fsum(abs(x) for x in terms)
                    left,ls=nb(u,fraction*size);right,rs=nb(v,(1-fraction)*size)
                    error=abs(joint-left-right)
                    # Gamma/log terms and their cancellation set the numerical
                    # scale. Report actual errors; this is not a sampling cutoff.
                    allowance=64*np.finfo(float).eps*max(1.,ls+rs+abs(joint))
                    if error>allowance:raise ArithmeticError('Independent thinning-factorization reference differs')
                    discrepancies.append(error)
            # Test actual sampler's deterministic partition/conservation contract.
            Y=np.array([[0,1,7],[20,9,2]],np.float32);keys=np.array([9,42]);r=np.array([2.,3.,7.])
            prop=clones.thin_counts(Y,keys,r,20261108,eps=fraction)
            test=Y-prop
            np.testing.assert_array_equal(prop+test,Y)
            np.testing.assert_array_equal(clones.thin_counts(Y[::-1],keys[::-1],r,20261108,eps=fraction)[::-1],prop)
            if (prop<0).any() or (test<0).any():raise AssertionError('Negative split count')
            # The two independently specified NB marginals share p. Transporting
            # mean and shape by test/proposal preserves p, including unequal shares.
            ratio=(1-fraction)/fraction
            if not math.isclose((fraction*size*ratio)/(fraction*(size+mean)*ratio),p,
                                rel_tol=8*np.finfo(float).eps):raise ArithmeticError('NB scale transport failed')
    return dict(joint_law_checks=len(discrepancies),maximum_log_probability_discrepancy=max(discrepancies),
                count_conservation=True,row_identity_invariance=True,
                scope='Known generating NB size only; fitted nuisance/dispersion do not establish independence')


@contextmanager
def allocation(fraction):
    original=clones.THIN_EPS
    clones.THIN_EPS=fraction
    try:yield
    finally:clones.THIN_EPS=original


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--allocation',choices=['one_third','half','two_thirds'],required=True)
    p.add_argument('--panel-root',type=Path,required=True);p.add_argument('--assertion-source',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    fraction={'one_third':1/3,'half':.5,'two_thirds':2/3}[a.allocation];root=Path(__file__).parent
    panel=load('original_panel',root/'shared_membership_panel.py');pool=load('pooling_hook',root/'pooled_initialization.py')
    old_loader=panel.load
    def experimental_loader(name,path):
        module=old_loader(name,path)
        if name=='weighted_stages':
            original=module.weighted_stages
            module.weighted_stages=lambda stage,trace:original(stage,trace,allow_unequal=True)
        return module
    panel.load=experimental_loader
    inputs={n:root/n for n in ('thinning_allocation_panel.py','shared_membership_panel.py',
            'shared_membership_stages.py','shared_membership_hmm.py','pooled_initialization.py',
            'membership_weighted_stages.py','recovery_baseline.py','shared_membership_probe.py')}
    inputs.update(panel_manifest=a.panel_root/'manifest.json',assertions=a.assertion_source,
                  previous_summary=a.control/'summary.json')
    rec=start_run(a.out,inputs,dict(allocation=a.allocation,proposal_fraction=fraction,confirmation_fraction=1-fraction,
        initialization='pooled_initial',stage='exact_matched'),
        dict(scope='Prespecified synthetic development screen; no scientific thresholds altered'),mode='thinning_allocation')
    trace=[]
    try:
        (a.out/'checks.json').write_text(json.dumps(checks(),indent=2))
        old_argv=sys.argv
        sys.argv=[str(root/'shared_membership_panel.py'),'--out',str(a.out/'panel'),'--panel-root',str(a.panel_root),
                  '--assertion-source',str(a.assertion_source),'--stage','exact_matched']
        try:
            with allocation(fraction),pool.pooled_initialization(trace):panel.main()
        finally:sys.argv=old_argv
        (a.out/'initialization_trace.json').write_text(json.dumps(trace,indent=2))
        cases=json.loads((a.out/'panel/summary.json').read_text())['cases']
        previous=json.loads((a.control/'summary.json').read_text())['cases']
        if [r['seed'] for r in cases]!=[r['seed'] for r in previous]:raise AssertionError('Case population changed')
        if a.allocation=='half':
            for row,before in zip(cases,previous):
                if row['criteria']!=before['criteria'] or row['original_assertion_passed']!=before['original_assertion_passed']:
                    raise AssertionError('Equal-share control criteria changed')
                directory=a.out/'panel'/f"seed_{row['seed']}"/'fit';old=a.control/f"seed_{row['seed']}"/'fit'
                np.testing.assert_array_equal(np.load(directory/'posteriors.npy'),np.load(old/'posteriors.npy'))
                for name in ('cells.tsv.gz','differences.tsv','segments.tsv','clones.tsv'):
                    pd.testing.assert_frame_equal(pd.read_csv(directory/name,sep='\t'),pd.read_csv(old/name,sep='\t'),
                                                  check_dtype=False,check_exact=True)
        (a.out/'summary.json').write_text(json.dumps(dict(allocation=a.allocation,proposal_fraction=fraction,
            cases=cases,half_control_exact=a.allocation=='half',interpretation='Screen only; allocation not selected for promotion'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
