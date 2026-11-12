"""Frozen153-endpoint fresh comparison and independent numerical/data audit."""
import argparse
from decimal import Decimal,localcontext
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq,linear_sum_assignment
from scipy.stats import beta,binomtest
from haplodose.provenance import start_run,finish_run,file_record

ARMS=('joint_proposal','joint_full')
ROLES={'cal_nb':('matched_nb',110,199),'cal_state':('correlated_state',111,199),
       'eval_nb':('matched_nb',120,200),'eval_state':('correlated_state',120,200),'eval_two':('two_clones',120,200)}
CONDITIONS=('matched_nb','correlated_state','two_clones')
RULES=[('raw',None,None)]+[(method,bank,alpha) for method,bank in
    [('NB_calibrated_gate','cal_nb'),('state_calibrated_gate','cal_state')] for alpha in (.05,.005,.01,.025)]
FAMILY=153


def verify(r):
    if file_record(r['path'])['sha256']!=r['sha256']:raise AssertionError('Recorded artifact changed')


def collect(a):
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',previous_summary=a.previous/'summary_v1/summary.json',
        previous_audit=a.previous/'audit_v1/report.json');rows={}
    for arm in ARMS:
        rows[arm]={}
        for role,(condition,seed_role,n) in ROLES.items():
            rows[arm][role]=[]
            for start in (0,50,100,150):
                directory=a.task/'chunks_v1'/arm/f'{role}_{start}'
                mp=directory/'run_manifest.json';m=json.loads(mp.read_text());params=m['execution']['parameters']
                if m['status']!='complete':raise AssertionError('Incomplete chunk')
                if (params['arm'],params['role'],params['seed_role'],params['start'],params['stop'],params['confidence_family'])!=(arm,role,seed_role,start,min(start+50,n),FAMILY):raise AssertionError('Changed protocol')
                for r in m['input']['files'].values():verify(r)
                verify(m['output']['files']['records.jsonl'])
                path=directory/'records.jsonl';part=[]
                # Traces contain full cell matrices. Verify their operative
                # mechanism one record at a time, then retain only sufficient
                # metadata for the prespecified summaries. Raw files stay intact.
                with path.open() as stream:
                    for line in stream:
                        row=json.loads(line);trace=row.pop('stage_trace');row.pop('initialization_trace')
                        counts=dict(joint_proposals=0,joint_scoring_calls=0,merge_pairs=0,full_assignment_calls=0)
                        for entry in trace:
                            stage=entry['stage']
                            if stage=='joint_proposal':counts['joint_proposals']+=1
                            elif stage=='joint_scoring':counts['joint_scoring_calls']+=len(entry['calls'])
                            elif stage=='full_merge_assignment':counts['full_assignment_calls']+=1
                            elif stage=='merge_membership':
                                counts['merge_pairs']+=len(entry['pairs'])
                                expected='joint_full' if arm=='joint_full' else 'proposal'
                                if entry['assignment_source']!=expected:raise AssertionError('Wrong merge mechanism')
                                if arm=='joint_full' and entry['merge_responsibility']!=entry['full_responsibility']:
                                    raise AssertionError('Full merge did not use supplied memberships')
                        row['mechanism_counts']=counts;part.append(row)
                if [r['replicate'] for r in part]!=list(range(start,min(start+50,n))):raise AssertionError('Missing/reordered draws')
                for row in part:
                    if row['seed_key']!=[20261107,seed_role,CONDITIONS.index(condition),row['replicate']]:raise AssertionError('Changed/familiar seed')
                rows[arm][role].extend(part);inputs[f'{arm}:{role}:{start}:records']=path;inputs[f'{arm}:{role}:{start}:manifest']=mp
    return rows,inputs


def scipy_bound(k,n,power):
    tail=.05/FAMILY
    return (0. if k==0 else float(beta.ppf(tail,k,n-k+1))) if power else (1. if k==n else float(beta.isf(tail,k+1,n-k)))


def bound(k,n,power,precision=80):
    """Invert the exact binomial tail by high-precision recurrence/bisection.

    Fixed100 bisections enclose the root to2^-100 before float conversion. The
    preflight checks60/80 digit convergence against a separately written binary
    polynomial enclosure for every possible event count; no empirical epsilon.
    """
    if power and k==0:return 0.
    if not power and k==n:return 1.
    with localcontext() as context:
        context.prec=precision;one=Decimal(1);tail=one/(Decimal(20)*FAMILY)
        lo=Decimal(0);hi=one
        for _ in range(100):
            p=(lo+hi)/2;odds=p/(one-p);term=(one-p)**n;probability=Decimal(0)
            for j in range(n+1):
                if (j>=k if power else j<=k):probability+=term
                if j<n:term=term*(n-j)/(j+1)*odds
            if (probability<tail if power else probability>tail):lo=p
            else:hi=p
        return float((lo+hi)/2)


def member(row,independent=False):
    verify(row['assignments'])
    with np.load(row['assignments']['path']) as data:
        truth=data['simulated_CNA_labels'];labels=data['raw_MAP'];eligible=data['eligible_tumour_groups']
        if set(np.unique(truth))!={0,1}:raise AssertionError('Alternative truth changed')
        table=np.array([[int(((truth==t)&(labels==k)).sum()) for k in eligible] for t in (0,1)])
        if independent:
            best=max((sum(int(((truth==t)&(labels==k)).sum()) for t,k in enumerate(pair))
                for pair in itertools.permutations(list(eligible)+[-2,-3],2)),default=0)
        elif table.size:
            rr,cc=linear_sum_assignment(-table);best=int(table[rr,cc].sum())
        else:best=0
        return dict(correct=best/len(truth),coverage=float(np.isin(labels,eligible).mean()),exactly_two=len(eligible)==2)


def summarize(a,rows,inputs):
    previous=json.loads(inputs['previous_summary'].read_text())
    if previous['family']!=99 or len(previous['summary'])+len(previous['retained_historical'])!=99:raise AssertionError('Historical family incomplete')
    retained=[]
    for row in previous['retained_historical']+previous['summary']:
        old={k:v for k,v in row.items() if k!='independent_bound'};power=row['direction']=='lower';value=bound(row['events'],row['n'],power)
        old.update(previous_bound=dict(family=row['family'],value=row['simultaneous_bound']),family=FAMILY,
            simultaneous_bound=value,error_criterion_met=None if power else value<=.05,evidence_role='historical_exposed_regression')
        retained.append(old)
    summary=[];decisions=[];coverage=[]
    for arm in ARMS:
        for role in ROLES:
            coverage.append(dict(arm=arm,role=role,**{key:sum(r['mechanism_counts'][key] for r in rows[arm][role])
                for key in ('joint_proposals','joint_scoring_calls','merge_pairs','full_assignment_calls')}))
        membership={r['replicate']:member(r) for r in rows[arm]['eval_two']}
        for role in ('eval_nb','eval_state','eval_two'):
            condition=ROLES[role][0];power=role=='eval_two'
            for method,bank,alpha in RULES:
                scores=None if bank is None else np.array([r['score'] for r in rows[arm][bank]])
                accepted=[]
                for row in rows[arm][role]:
                    p=None if scores is None else float((1+np.count_nonzero(scores>=row['score']))/200)
                    event=bool(row['raw_candidate_count']>1 and (bank is None or (row['score']>0 and p<=alpha)))
                    accepted.append(event);decisions.append(dict(arm=arm,condition=condition,method=method,alpha=alpha,
                        replicate=row['replicate'],event=event,rank_p=p,score=row['score'],raw_groups=row['raw_candidate_count']))
                k,n=sum(accepted),len(accepted);value=bound(k,n,power)
                record=dict(arm=arm,condition=condition,method=method,alpha=alpha,events=k,n=n,family=FAMILY,
                    simultaneous_bound=value,direction='lower' if power else 'upper',error_criterion_met=None if power else value<=.05,
                    evidence_role='fresh_conditional_simulation')
                if power:
                    record.update(exactly_two_and_accepted=sum(bool(e and membership[i]['exactly_two']) for i,e in enumerate(accepted)),
                        mean_correct_fraction_all_cells_if_gate_applied=sum(membership[i]['correct']*e for i,e in enumerate(accepted))/n,
                        mean_coverage_if_gate_applied=sum(membership[i]['coverage']*e for i,e in enumerate(accepted))/n)
                summary.append(record)
    paired=[]
    for condition in CONDITIONS:
        for method,_,alpha in RULES:
            before,after=[[r['event'] for r in decisions if (r['arm'],r['condition'],r['method'],r['alpha'])==(arm,condition,method,alpha)] for arm in ARMS]
            lost=sum(x and not y for x,y in zip(before,after));gained=sum(y and not x for x,y in zip(before,after))
            paired.append(dict(condition=condition,method=method,alpha=alpha,lost=lost,gained=gained,
                nominal_p=1. if lost+gained==0 else float(binomtest(gained,lost+gained,.5).pvalue)))
    largest=0.
    for rank,i in enumerate(sorted(range(27),key=lambda i:paired[i]['nominal_p'])):
        largest=max(largest,min(1.,(27-rank)*paired[i]['nominal_p']));paired[i]['holm_p']=largest
    (a.out/'summary.json').write_text(json.dumps(dict(summary=summary,retained_historical=retained,paired=paired,family=FAMILY,mechanism_coverage=coverage,
        scope='Fresh paired known-reference bin-only challenge; joint data, not biological validation or public raw-input error control'),indent=2))
    (a.out/'decisions.json').write_text(json.dumps(decisions,indent=2))


def independent_enclosure(k,n,power):
    if (power and k==0) or (not power and k==n):return (0.,0.) if power else (1.,1.)
    # Integer binomial coefficients and scalar powers/sums, independently of
    # inverse beta. Bound rounding of at most2n+16 operations per weighted term;
    # bracket the root under both signs of that relative error.
    js=range(k,n+1) if power else range(k+1);coefficients=[(j,math.comb(n,j)) for j in js]
    eps=np.finfo(float).eps;ops=2*n+16;gamma=ops*eps/(1-ops*eps);tail=.05/FAMILY
    def probability(p):return math.fsum(c*p**j*(1-p)**(n-j) for j,c in coefficients)
    roots=[brentq(lambda p:factor*probability(p)-tail,0.,1.,xtol=np.nextafter(0.,1.),rtol=4*eps)
           for factor in (1-gamma,1+gamma)]
    error=4*eps*max(map(abs,roots))+np.nextafter(0.,1.)
    return min(roots)-error,max(roots)+error


def audit(a,rows,inputs):
    summary=json.loads((a.task/'summary_v1/summary.json').read_text());decisions=json.loads((a.task/'summary_v1/decisions.json').read_text())
    if summary['family']!=FAMILY or len(summary['summary'])!=54 or len(summary['retained_historical'])!=99 or len(decisions)!=10800:raise AssertionError('Endpoint population changed')
    inputs_paired=0
    for role,(_,_,n) in ROLES.items():
        for before,after in zip(rows[ARMS[0]][role],rows[ARMS[1]][role]):
            if before['input_sha256']!=after['input_sha256'] or before['seed_key']!=after['seed_key']:raise AssertionError('Unpaired fresh inputs')
            for row in (before,after):
                verify(row['assignments'])
                with np.load(row['assignments']['path']) as z:
                    if not np.array_equal(z['raw_R'].argmax(axis=1),z['raw_MAP']):raise AssertionError('MAP differs from probabilities')
            inputs_paired+=1
    index={}
    for r in decisions:
        arm=r['arm'];role=next(k for k in ('eval_nb','eval_state','eval_two') if ROLES[k][0]==r['condition']);row=rows[arm][role][r['replicate']]
        bank=None if r['method']=='raw' else rows[arm]['cal_nb' if r['method']=='NB_calibrated_gate' else 'cal_state']
        integer_rank=None if bank is None else 1+sum(v['score']>=row['score'] for v in bank)
        cutoff=None if bank is None else {0.05:10,0.005:1,0.01:2,0.025:5}[r['alpha']]
        event=bool(row['raw_candidate_count']>1 and (bank is None or (row['score']>0 and integer_rank<=cutoff)))
        if event!=r['event'] or (None if bank is None else integer_rank/200)!=r['rank_p']:raise AssertionError('Independent rank decision differs')
        key=(arm,r['condition'],r['method'],r['alpha'],r['replicate'])
        if key in index:raise AssertionError('Duplicate decision')
        index[key]=event
    enclosures=[]
    prior=json.loads(inputs['previous_summary'].read_text());history=prior['retained_historical']+prior['summary']
    for pos,r in enumerate(summary['summary']+summary['retained_historical']):
        if r['family']!=FAMILY:raise AssertionError('Incorrect per-endpoint family')
        if pos<54:
            k=sum(index[(r['arm'],r['condition'],r['method'],r['alpha'],i)] for i in range(200))
            if k!=r['events'] or r['n']!=200:raise AssertionError('Endpoint count differs')
        else:
            old=history[pos-54]
            if (r['events'],r['n'],r['direction'])!=(old['events'],old['n'],old['direction']):raise AssertionError('Historical counts changed')
            if r['previous_bound']!={'family':99,'value':old['simultaneous_bound']}:raise AssertionError('History overwritten')
        power=r['direction']=='lower';lo,hi=independent_enclosure(r['events'],r['n'],power)
        if not lo<=r['simultaneous_bound']<=hi:raise AssertionError('Bound outside independent rounding enclosure')
        if r['error_criterion_met']!=(None if power else r['simultaneous_bound']<=.05):raise AssertionError('Stale pass flag')
        enclosures.append(dict(events=r['events'],n=r['n'],lo=lo,hi=hi))
    for arm in ARMS:
        metrics={r['replicate']:member(r,independent=True) for r in rows[arm]['eval_two']}
        for r in [x for x in summary['summary'] if x['arm']==arm and x['condition']=='two_clones']:
            accepted=[index[(arm,'two_clones',r['method'],r['alpha'],i)] for i in range(200)]
            for field,key in [('mean_correct_fraction_all_cells_if_gate_applied','correct'),('mean_coverage_if_gate_applied','coverage')]:
                value=sum(metrics[i][key]*e for i,e in enumerate(accepted))/200
                if value!=r[field]:raise AssertionError('Independent membership metric differs')
            if sum(e and metrics[i]['exactly_two'] for i,e in enumerate(accepted))!=r['exactly_two_and_accepted']:raise AssertionError('Wrong exact recovery count')
    if len(summary['paired'])!=27:raise AssertionError('Paired family changed')
    for r in summary['paired']:
        keys=[(r['condition'],r['method'],r['alpha'],i) for i in range(200)]
        lost=sum(index[(ARMS[0],)+k] and not index[(ARMS[1],)+k] for k in keys)
        gained=sum(index[(ARMS[1],)+k] and not index[(ARMS[0],)+k] for k in keys)
        if (lost,gained)!=(r['lost'],r['gained']):raise AssertionError('Paired contrast differs')
    baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());core=0
    for name,digest in baseline.items():
        if file_record(a.task/'baseline'/name)['sha256']!=digest:raise AssertionError('Baseline changed')
        if name.startswith(('python/','tests/','crates/')):
            if file_record(a.task/'candidate_v1'/name)['sha256']!=digest:raise AssertionError('Executed core changed')
            core+=1
    artifact=json.loads((a.task.parents[1]/'dev/terminal_wheel_20261108/validation_22991908/installed_artifact.json').read_text())
    verify(artifact['wheel']);verify(artifact['native']);jobs=json.loads((a.task/'jobs.json').read_text())
    if len(jobs)!=40:raise AssertionError('Wrong inference job population')
    for job in jobs:
        jobid=job['stdout'].strip().split()[-1];execution=json.loads((a.task/f'analysis_{jobid}/execution.json').read_text())
        if Path(execution['source']).resolve()!=(a.task/'candidate_v1').resolve():raise AssertionError('Wrong executed source')
        if execution['wheel']['sha256']!=artifact['wheel']['sha256'] or execution['native_sha256']!=artifact['native']['sha256']:raise AssertionError('Wrong execution artifact')
    (a.out/'report.json').write_text(json.dumps(dict(paired_fresh_inputs=inputs_paired,independent_decisions=len(index),
        independent_bound_enclosures=enclosures,paired_contrasts=27,membership_endpoints=18,baseline_files=len(baseline),
        executed_core_files=core,verified_inference_jobs=len(jobs),scope='Conditional simulation, not biological validation'),indent=2))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);p.add_argument('--mode',choices=['summary','audit'],required=True);a=p.parse_args()
    rows,inputs=collect(a)
    if a.mode=='audit':inputs.update(summary=a.task/'summary_v1/summary.json',decisions=a.task/'summary_v1/decisions.json',jobs=a.task/'jobs.json')
    a.out.mkdir(exist_ok=False);rec=start_run(a.out,inputs,dict(mode=a.mode,family=FAMILY),
        dict(scope='Fresh known-reference binned experiment with retained exposed history'),mode='joint_fresh_'+a.mode)
    try:
        (summarize if a.mode=='summary' else audit)(a,rows,inputs);finish_run(a.out,rec,backend='cpu')
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
