"""Independent threshold, endpoint, identity and indexed-bootstrap preflight."""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from haplodose.provenance import start_run,finish_run,file_record


def load(name):
    import importlib.util
    spec=importlib.util.spec_from_file_location('uncertainty_audit_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def scalar_scores(values,loss,neutral,gain,common,threshold):
    result=np.full((values.shape[1],3),np.nan)
    if not np.isfinite(threshold):return result
    for j in range(values.shape[1]):
        for k,mask in enumerate((loss,neutral)):
            selected=[float(values[i,j])<threshold for i in range(len(values)) if mask[i,j] and common[i,j]]
            if selected:result[j,k]=sum(selected)/len(selected)
        classified=[i for i in range(len(values)) if common[i,j] and np.isfinite(values[i,j])
                    and (loss[i,j] or neutral[i,j] or gain[i,j]) and values[i,j]<threshold]
        if classified:result[j,2]=sum(bool(loss[i,j]) for i in classified)/len(classified)
    return result


def fixtures(candidate,ranks):
    values=np.array([[0.,np.nan,2.],[1.,np.inf,2.],[2.,-np.inf,3.],[3.,4.,4.]])
    neutral=np.ones_like(values,bool);ci=np.arange(3);cache=candidate.Calibration(values,neutral,ci);checks=0
    for columns in itertools.product(ci,repeat=4):
        expected=ranks.threshold_reference(values,neutral,columns);actual=cache.threshold(columns)
        assert actual==expected or (np.isnan(actual) and np.isnan(expected));checks+=1
    empty=candidate.Calibration(values,np.zeros_like(neutral),ci)
    assert np.isnan(empty.threshold(ci));checks+=1
    rng=np.random.default_rng(20261109)
    for _ in range(24):
        scores=rng.choice([-np.inf,-1.,0.,1.,np.inf,np.nan],size=(11,4))
        labels=rng.integers(0,4,size=scores.shape);loss=labels==0;neutral=labels==1;gain=labels==2
        common=rng.random(scores.shape)>.2;curve=candidate.LossCurve(scores,loss,neutral,gain,common)
        for threshold in (-np.inf,-2.,-1.,0.,1.,2.,np.inf,np.nan):
            assert np.array_equal(curve.at(threshold),scalar_scores(scores,loss,neutral,gain,common,threshold),equal_nan=True)
            checks+=1
    return checks


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--draw-root',default='pilot_v1')
    a=p.parse_args();a.out.mkdir(exist_ok=False);directory=a.task/a.draw_root
    candidate=load('bulk_amplitude_uncertainty');ranks=load('bulk_posterior_stage_audit')
    manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    bank,corr,methods,inputs=candidate.read_bank(a.task.parents[1])
    inputs.update(script=Path(__file__),candidate=Path(candidate.__file__),ranks=Path(ranks.__file__),design=a.task/'DESIGN.md',
        manifest=directory/'run_manifest.json',definition=directory/'definition.json',draws=directory/'draws.npy',
        indices=directory/'indices.npy',paired=directory/'paired.npy',thresholds=directory/'thresholds.npy')
    assert file_record(inputs['candidate'])['sha256']==manifest['input']['files']['script']['sha256']
    rec=start_run(a.out,inputs,{},dict(scope='Computational audit only,not empirical confidence coverage'),mode='bulk_amplitude_uncertainty_audit')
    try:
        fixture_checks=fixtures(candidate,ranks);ci=np.flatnonzero(bank[0]['split']=='calibration');ei=np.flatnonzero(bank[0]['split']=='evaluation')
        definition=json.loads(inputs['definition'].read_text());contrasts=candidate.definitions(methods);assert definition['contrasts']==contrasts
        draws=np.load(inputs['draws'],mmap_mode='r');indices=np.load(inputs['indices'],mmap_mode='r')
        paired=np.load(inputs['paired'],mmap_mode='r');thresholds=np.load(inputs['thresholds'],mmap_mode='r')
        report=json.loads((directory/'report.json').read_text());start,stop=report['start'],report['stop']
        assert draws.shape==(stop-start,80) and paired.shape==draws.shape and indices.shape==(stop-start,49)
        assert thresholds.shape==(stop-start,12)
        cached=[{m:candidate.Calibration(w[m],w['neutral']&w['common'],ci) for m in methods} for w in bank]
        curves=[{m:candidate.LossCurve(w[m],w['loss'],w['neutral'],w['gain'],w['common']) for m in methods} for w in bank]
        thresholds_checked=0;endpoints_checked=0;effects_checked=0
        rng=np.random.default_rng(candidate.SEED);inspect=set(range(start,min(start+64,stop)))|{stop-1}
        for k in range(stop):
            ri=int(rng.integers(1,9));cal=rng.choice(ci,24,replace=True);chosen=rng.choice(ei,24,replace=True)
            if k<start:continue
            row=k-start;assert np.array_equal(indices[row],np.r_[ri,cal,chosen])
            if k not in inspect:continue
            w=bank[ri];scored={}
            for col,m in enumerate(methods):
                threshold=ranks.threshold_reference(w[m],w['neutral']&w['common'],cal)
                assert threshold==thresholds[row,col]==cached[ri][m].threshold(cal);thresholds_checked+=1
                expected=scalar_scores(w[m],w['loss'],w['neutral'],w['gain'],w['common'],threshold)
                assert np.array_equal(expected,curves[ri][m].at(threshold),equal_nan=True)
                endpoints_checked+=expected.size
                scored[m]=np.column_stack([corr[ri][m],expected])
            for col,r in enumerate(contrasts):
                j=candidate.ENDPOINTS.index(r['endpoint'])
                vector=[float(scored[r['method']][i,j]-scored[r['control']][i,j])*r['direction'] for i in chosen]
                vector=[v for v in vector if np.isfinite(v)];assert len(vector)==paired[row,col]
                if not vector:assert np.isnan(draws[row,col]);continue
                expected=math.fsum(vector)/len(vector);u=np.finfo(float).eps/2
                gamma=(len(vector)+3)*u/(1-(len(vector)+3)*u)
                guard=gamma*math.fsum(abs(v) for v in vector)/len(vector)
                assert abs(expected-draws[row,col])<=guard,(k,col,expected,draws[row,col],guard)
                effects_checked+=1
        # Anchor point effects to the preceding independent paired audit.
        previous=a.task.parents[1]/'dev/bulk_amplitude_reference_20261109/audit_000_v1/paired_effects.json'
        original={(r['method'],r['control'],r['endpoint']):r for r in json.loads(previous.read_text())}
        for i,r in enumerate(contrasts):
            row=original[r['method'],r['control'],r['endpoint']]
            assert row['paired_patients']==definition['point_paired_patients'][i]
            expected=r['direction']*row['mean_paired_effect'];n=row['paired_patients'];u=np.finfo(float).eps/2
            assert abs(expected-definition['point_effects'][i])<=(n+3)*u*(1+abs(expected))
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        (a.out/'report.json').write_text(json.dumps(dict(all_checks_passed=True,fixture_checks=fixture_checks,
            indexed_draws_verified=stop-start,exact_fraction_threshold_checks=thresholds_checked,
            independent_patient_endpoints=endpoints_checked,scalar_paired_effects=effects_checked,
            original_point_effects_anchored=80,coverage_validated=False),indent=2)+'\n')
        finish_run(a.out,rec,backend='CPU Fraction thresholds,scalar truth counts and compensated paired sums')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
