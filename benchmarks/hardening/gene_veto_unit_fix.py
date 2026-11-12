"""Independent weight-unit gate for a bounded research gene-veto correction."""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path

import numpy as np
from haplodose.io import write_json
from haplodose.profile import pooled_lograte
from haplodose.provenance import start_run, finish_run


def load(name,path):
    sp=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m


def exact_units(W):
    return [sum(Fraction(float(w))**2 for w in W[:,j])/sum(Fraction(float(w)) for w in W[:,j]) for j in (0,1)]


def check(veto,reference,S,W,gbin,mask,seg,minimum=8,alpha=.01):
    a,b,c=veto.effective_counts(S[0],S[1],W)
    rational=exact_units(W)
    budget=64*len(W)*np.finfo(float).eps
    assert np.all(abs(c-np.array([float(x) for x in rational])) <= budget*np.maximum(1.,abs(c)))
    independent=np.array([[float(Fraction(float(x))/rational[j]) for x in S[j]] for j in (0,1)])
    assert np.all(abs(np.array([a,b])-independent) <= budget*np.maximum(1.,abs(independent)))
    L=S.sum(axis=1)
    actual=pooled_lograte(S,L,S,c)[0]
    normalized=np.log((np.array([a,b])+.5)/(L/c)[:,None])
    assert np.all(abs(actual-normalized)<=64*np.finfo(float).eps*np.maximum(1.,abs(actual)))
    expected=reference.exact_sign(independent[0],independent[1],gbin,mask,seg,minimum,alpha)
    baseline=veto.decision(a,b,gbin,mask,seg,minimum,alpha)
    assert baseline['passed']==expected['passed']
    if 'n' in baseline:assert (baseline['n'],baseline['k'])==(expected['n'],expected['k'])
    scales=(.125,.5,1.,2.,8.);checks=[]
    for fa in scales:
        for fb in scales:
            factors=np.array([fa,fb]);new=veto.effective_counts(fa*S[0],fb*S[1],W*factors)
            assert np.array_equal(new[0],a) and np.array_equal(new[1],b)
            assert veto.decision(new[0],new[1],gbin,mask,seg,minimum,alpha)==baseline
            checks.append([fa,fb])
    return dict(decision=baseline,exact=expected,weight_scale=c.tolist(),scale_pairs_checked=len(checks))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=Path(__file__).parent
    veto=load('unit_veto',root/'gene_veto_weights.py');reference=load('exact_reference',root/'gene_veto_weights_audit.py')
    prior=a.task.parent/'gene_veto_aligned_20261109/candidate_v1/report.json';calls=json.loads(prior.read_text())['calls']
    inputs=dict(script=Path(__file__),veto=root/'gene_veto_weights.py',reference=root/'gene_veto_weights_audit.py',design=a.task/'DESIGN.md',prior=prior)
    for row in calls:
        inputs['model:'+row['model']]=Path(row['model']);inputs['call:'+str(row['index'])]=Path(row['arrays'])
    rec=start_run(a.out,inputs,{},dict(scope='Weight-unit correction only; no gene-sign calibration'),mode='gene_veto_unit_fix')
    try:
        S=np.array([np.r_[np.full(4,20.),np.full(4,10.),np.full(16,20.)],
            np.r_[np.full(4,10.),np.full(4,20.),np.full(16,20.)]])
        gbin=np.r_[np.zeros(8,int),np.ones(16,int)];mask=np.ones(2,bool);seg=dict(first_bin=0,last_bin=0,direction='higher')
        fixture=check(veto,reference,S,np.eye(2),gbin,mask,seg)
        assert fixture['decision']['passed'] is False and Fraction(fixture['exact']['tail'])==Fraction(163,256)
        rng=np.random.default_rng(2026110914);binary=[]
        for i in range(32):
            W=np.zeros((20,2));W[:10,0]=rng.integers(0,2,10);W[10:,1]=rng.integers(0,2,10)
            W[0,0]=W[-1,1]=1
            X=rng.integers(0,40,(20,24));pooled=W.T@X
            ca,cb,c=veto.effective_counts(*pooled,W)
            assert np.array_equal(c,np.ones(2)) and np.array_equal(np.array([ca,cb]),pooled)
            assert veto.decision(ca,cb,gbin,mask,seg)==veto.decision(*pooled,gbin,mask,seg)
            binary.append(check(veto,reference,pooled,W,gbin,mask,seg))
        invalid=0
        for W in (np.zeros((2,2)),np.array([[1.,-1.],[0.,1.]]),np.array([[np.nan,1.],[1.,0.]]),np.ones((2,3))):
            try:veto.effective_counts(S[0],S[1],W)
            except ValueError:invalid+=1
            else:raise AssertionError('Invalid pool units accepted')
        real=[]
        for row in calls:
            with np.load(row['arrays']) as z:
                if 'weights' not in z:continue
                with np.load(row['model']) as m:
                    W=z['weights'][m['sub']];gbin=m['gbin'];mask=m['mask'];seed=int(m['seed'])
                result=check(veto,reference,np.array([z['weighted_a'],z['weighted_b']]),W,gbin,mask,
                    row['seg'],row['min_genes'],row['alpha'])
                real.append(dict(seed=seed,index=row['index'],**result))
        write_json(dict(all_checks_passed=True,fixture=fixture,binary_fixtures=len(binary),invalid_fixtures=invalid,
            real_pools=real,scale_pairs_per_case=25,biological_specimens=0,production_changed=False,
            scientific_criteria_changed=False,scope='Working weight units; inherited selected gene signs remain uncalibrated'),a.out/'report.json')
    except BaseException as error:finish_run(a.out,rec,error=error);raise
    finish_run(a.out,rec,backend='CPU exact-rational unit/sign references and invariants')


if __name__=='__main__':main()
