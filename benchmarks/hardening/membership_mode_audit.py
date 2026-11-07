"""Audit the thinning comparison and its paired even/odd control without relaxing identity checks."""
import argparse
from dataclasses import asdict
import importlib.util
import json
from pathlib import Path
import numpy as np
from scipy.stats import binomtest
from haplodose.clones import CloneParams, Rows, fit_specimen


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--generator',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    audit=module('independent_bounds',Path(__file__).with_name('membership_evidence_audit.py'))
    generator=module('fixed_generator',a.generator)
    read,key,bound=audit.read,audit.key,audit.bound
    even=a.task/'predictive_rank_v4';thin=a.task/'thin_rank_v1'
    previous=json.loads((a.task/'evidence_v2/report.json').read_text())
    assert previous['status']=='completed_independent_recount'
    old={key(r):r for f in ['calibration.jsonl','evaluation.jsonl'] for r in read(even/f)}
    new={key(r):r for f in ['calibration.jsonl','evaluation.jsonl'] for r in read(thin/f)}
    assert set(new)==set(old) and len(new)==998
    if any(r['input_sha256']!=old[k]['input_sha256'] for k,r in new.items()):
        raise ValueError('Thinning experiment changed input observations or latent truth')
    params=CloneParams(n_programmes=0,spatial_neighbors=0,seed=20261107,split_mode='thin')
    actual=json.loads((thin/'manifest.json').read_text())['actual_params']
    assert actual==json.loads(json.dumps(asdict(params)))
    replay=[]
    for index,condition in enumerate(['matched_nb','correlated_state','two_clones']):
        row=new[(condition,20,0)]
        rng=np.random.default_rng(np.random.SeedSequence(row['seed_key']))
        Y,base,phi,truth,latent=generator.generate(rng,condition,600,240,return_latent=True)
        reset=np.zeros(240,bool);reset[::60]=True
        result=fit_specimen(Rows(Y,np.arange(len(Y))),Y.sum(axis=1,dtype=np.float64),base[:,None],reset,phi,params)
        saved=np.load(row['assignments']['path'])
        if not np.array_equal(result['R'],saved['raw_R']):
            raise ValueError('Thinning observer changed posterior probabilities in a prespecified control')
        replay.append(dict(condition=condition,replicate=0,posterior_identical=True))
    records=read(thin/'comparison_evaluation.jsonl');cal=read(thin/'calibration.jsonl')
    for row in records:
        for prefix,score in [('', 'score'),('hard_','hard_score')]:
            for name,c in [('NB_calibrated_gate','matched_nb'),('state_calibrated_gate','correlated_state')]:
                prob=(1+sum(r[score]>=row[score] for r in cal if r['condition']==c))/200
                method=prefix+name
                if prob!=row['methods'][method]['rank_p']:raise ValueError('Rank mismatch')
                expected=row['conditional_proposal_count']>0 and row[score]>0 and prob<=.05
                if expected!=row['methods'][method]['accepted_multigroup']:raise ValueError('Gate mismatch')
    intervals=[];paired=[]
    for mode,root in [('evenodd',even),('thin',thin)]:
        cmp=json.loads((root/'membership_comparison.json').read_text())
        rr=read(root/'comparison_evaluation.jsonl')
        for row in cmp['summary']:
            selected=[r for r in rr if r['condition']==row['condition']]
            n=len(selected);k=sum(r['methods'][row['method']]['accepted_multigroup'] for r in selected)
            if (k,n)!=(row['events'],row['n']):raise ValueError('Count mismatch')
            independent=bound(k,n,.05/15,row['direction'])
            intervals.append(dict(mode=mode,**row,independent_family15_bound=independent,
                saved_minus_reference=row['simultaneous_family15_bound']-independent,
                simultaneous_family30_bound=bound(k,n,.05/30,row['direction'])))
        for c in ['matched_nb','correlated_state','two_clones']:
            selected=[r for r in rr if r['condition']==c]
            for name in ['NB_calibrated_gate','state_calibrated_gate']:
                soft=np.array([r['methods'][name]['accepted_multigroup'] for r in selected])
                hard=np.array([r['methods']['hard_'+name]['accepted_multigroup'] for r in selected])
                a_only=int((soft & ~hard).sum());b_only=int((hard & ~soft).sum());n=a_only+b_only
                paired.append(dict(mode=mode,condition=c,method=name,soft_only=a_only,hard_only=b_only,
                    difference=(a_only-b_only)/len(selected),nominal_p=1. if n==0 else float(binomtest(a_only,n,.5).pvalue)))
    running=0.
    for rank,i in enumerate(sorted(range(len(paired)),key=lambda i:paired[i]['nominal_p'])):
        running=max(running,min(1.,(len(paired)-rank)*paired[i]['nominal_p']))
        paired[i]['holm_family12_p']=running
    # Preserve the first heterogeneous-CPU replay and report changed decisions.
    mixed={key(r):r for r in read(a.task/'predictive_rank_v3/comparison_evaluation.jsonl')}
    fixed={key(r):r for r in read(even/'comparison_evaluation.jsonl')}
    hardware=[]
    for k,r in fixed.items():
        for name in r['methods']:
            if r['methods'][name]['accepted_multigroup']!=mixed[k]['methods'][name]['accepted_multigroup']:
                hardware.append(dict(identity=k,method=name,first=mixed[k]['methods'][name],same_host=r['methods'][name],
                    first_scores={s:mixed[k][s] for s in ['score','hard_score','conditional_proposal_count']},
                    same_host_scores={s:r[s] for s in ['score','hard_score','conditional_proposal_count']},
                    scope='CPU-sensitive experimental statistic; no cross-hardware numerical-certification claim'))
    report=dict(status='completed_independent_mode_audit',identical_input_specimens=998,
        observer_controls=replay,intervals=intervals,paired=paired,hardware_sensitive_decisions=hardware,
        sources=[audit.record(root/f) for root in [even,thin] for f in ['membership_comparison.json','comparison_evaluation.jsonl','manifest.json']],
        scope='Known-reference binned simulations, all paired regression. Thirty simultaneous one-sided endpoints; no biological guarantee.')
    (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],observer_controls=replay,hardware_sensitive_decisions=hardware)))


if __name__=='__main__':main()
