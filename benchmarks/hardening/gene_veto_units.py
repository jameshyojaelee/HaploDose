"""Exact weight-unit counterexample and captured-pool numerical checks."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from haplodose.io import write_json
from haplodose.profile import pooled_lograte
from haplodose.provenance import start_run, finish_run


SCALES=(.125,.5,1.,2.,8.)


def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def check_moments(S,Q,L,c):
    base=pooled_lograte(S,L,Q,c);checks=[]
    for f in SCALES:
        actual=pooled_lograte(f*S,f*L,f*f*Q,f*c)
        direct=(np.log((f*S+.5*f*c[:,None])/(f*L[:,None])),
            (f*f*Q+.5*(f*c[:,None])**2)/(f*S+.5*f*c[:,None])**2)
        for a,b in zip(actual,direct):
            budget=64*np.finfo(float).eps*np.maximum(1.,abs(b))
            assert np.all(abs(a-b)<=budget)
        error=max(float(np.max(abs(a-b))) for a,b in zip(actual,base))
        assert all(np.all(abs(a-b)<=64*np.finfo(float).eps*np.maximum(1.,abs(b))) for a,b in zip(actual,base))
        checks.append(dict(scale=f,max_invariance_discrepancy=error))
    return checks


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    root=Path(__file__).parent
    veto=load('veto',root/'gene_veto_weights.py');reference=load('exact_signs',root/'gene_veto_weights_audit.py')
    prior=a.task.parent/'gene_veto_aligned_20261109/candidate_v1/report.json'
    calls=json.loads(prior.read_text())['calls']
    inputs=dict(script=Path(__file__),veto=root/'gene_veto_weights.py',reference=root/'gene_veto_weights_audit.py',
        design=a.task/'DESIGN.md',prior=prior)
    for row in calls:
        inputs['model:'+row['model']]=Path(row['model']);inputs['call:'+str(row['index'])]=Path(row['arrays'])
    rec=start_run(a.out,inputs,{},dict(scope='Research-veto weight-unit falsification; no pipeline or threshold change'),mode='gene_veto_units')
    try:
        ya=np.r_[np.full(4,20.),np.full(4,10.),np.full(16,20.)]
        yb=np.r_[np.full(4,10.),np.full(4,20.),np.full(16,20.)]
        gbin=np.r_[np.zeros(8,int),np.ones(16,int)];mask=np.ones(2,bool)
        seg=dict(first_bin=0,last_bin=0,direction='higher')
        fixture=[]
        for f in SCALES:
            got=veto.decision(f*ya,f*yb,gbin,mask,seg);exact=reference.exact_sign(f*ya,f*yb,gbin,mask,seg,8,.01)
            assert got['passed']==exact['passed']
            fixture.append(dict(scale=f,decision=got,exact=exact))
        assert fixture[0]['decision']['passed'] and not fixture[2]['decision']['passed']
        assert Fraction(fixture[2]['exact']['tail'])==Fraction(163,256)
        fixture_moments=check_moments(np.array([ya,yb]),np.array([ya,yb]),np.array([sum(ya),sum(yb)]),np.ones(2))
        actual,models=[],{}
        for row in calls:
            with np.load(row['arrays']) as z:
                if 'weights' not in z:continue
                if row['model'] not in models:
                    with np.load(row['model']) as m:
                        X=sparse.csr_matrix((m['data'],m['indices'],m['indptr']),shape=tuple(m['shape']))
                        models[row['model']]=(X,m['sub'].copy(),m['gbin'].copy(),m['mask'].copy(),int(m['seed']))
                X,sub,gbin,mask,seed=models[row['model']];W=z['weights'][sub]
                S=np.asarray(W.T@X);Q=np.asarray((W*W).T@X)
                L=S.sum(axis=1);c=(W*W).sum(axis=0)/W.sum(axis=0)
                checks=check_moments(S,Q,L,c)
                decisions=[]
                for f in SCALES:
                    got=veto.decision(f*S[0],f*S[1],gbin,mask,row['seg'],row['min_genes'],row['alpha'])
                    exact=reference.exact_sign(f*S[0],f*S[1],gbin,mask,row['seg'],row['min_genes'],row['alpha'])
                    assert got['passed']==exact['passed']
                    if 'n' in got:assert (got['n'],got['k'])==(exact['n'],exact['k'])
                    decisions.append(dict(scale=f,decision=got,exact=exact))
                actual.append(dict(seed=seed,index=row['index'],caller=row['caller'],
                    changed=len({d['decision']['passed'] for d in decisions})>1,
                    decisions=decisions,moment_checks=checks))
        write_json(dict(all_checks_passed=True,research_veto_unit_invariance_falsified=True,
            fixture=fixture,fixture_moments=fixture_moments,real_pools=actual,
            real_pools_with_decision_changes=sum(v['changed'] for v in actual),
            production_changed=False,scientific_criteria_changed=False,biological_specimens=0,
            scope='Unit falsification does not calibrate weighted sign tests or intervals'),a.out/'report.json')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finish_run(a.out,rec,backend='CPU exact rational sign tails and independent moment formulas')


if __name__=='__main__':main()
