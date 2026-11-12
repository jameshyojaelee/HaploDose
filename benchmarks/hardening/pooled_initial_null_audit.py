"""Independent input/assignment/decision/count/bound audit of initial pooling."""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import optimize, stats


def verify(record):
    if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()!=record['sha256']:
        raise AssertionError('Frozen file changed')


def read(path):return [json.loads(line) for line in path.read_text().splitlines()]


def preservation(task,root,validated_source=None,validated_artifact=None):
    baseline=json.loads((task/'baseline/source_manifest.json').read_text())
    frozen_count=live_unchanged=0;live_changes=[]
    if (validated_source is None)!=(validated_artifact is None):raise ValueError('Both correction identities are required')
    expected=None
    if validated_source is not None:
        expected=json.loads((validated_source/'source_manifest.json').read_text())
        artifact=json.loads(validated_artifact.read_text())
        verify(artifact['wheel']);verify(artifact['native']);verify(artifact['package'])
        test_report=validated_artifact.parent/'pytest.txt'
        if not test_report.read_text().strip().splitlines()[-1].startswith('321 passed, 1 skipped'):
            raise AssertionError('Specific precision candidate has not passed its full test suite')
        if 'test result: ok. 4 passed; 0 failed;' not in (validated_artifact.parent/'cargo_test.txt').read_text():
            raise AssertionError('Specific precision candidate Rust checks did not pass')
        for path,h in expected.items():
            if path.startswith(('python/','tests/','crates/')):
                if hashlib.sha256((validated_source/path).read_bytes()).hexdigest()!=h:
                    raise AssertionError('Validated correction source changed')
                if hashlib.sha256((root/path).read_bytes()).hexdigest()!=h:
                    raise AssertionError('Live correction differs from validated source')
    for path,h in baseline.items():
        if hashlib.sha256((task/'baseline'/path).read_bytes()).hexdigest()!=h:raise AssertionError('Baseline changed')
        if path.startswith(('python/','tests/','crates/')):
            # The old experiment executes this frozen source and its matching
            # immutable wheel, even when a separately validated live version exists.
            if hashlib.sha256((task/'candidate_v1'/path).read_bytes()).hexdigest()!=h:
                raise AssertionError('Executed old package/test/native source changed')
            frozen_count+=1
            actual=hashlib.sha256((root/path).read_bytes()).hexdigest()
            if actual==h:live_unchanged+=1
            elif expected is not None and path=='python/haplodose/clones.py' and actual==expected[path]:
                live_changes.append(dict(path=path,old_sha256=h,validated_sha256=actual))
            else:raise AssertionError('Unrecorded production/test/native change')
    # Every old shard launcher records the actual source, imported native object
    # and wheel, independently of the live checkout. Verify those identities too.
    executions=[]
    for job in (22964075,22964077,22964078,22964079,22964080):
        record=json.loads((task/f'analysis_{job}/execution.json').read_text())
        if Path(record['source']).resolve()!=(task/'candidate_v1').resolve():raise AssertionError('Unexpected executed source')
        verify(record['wheel'])
        installed=json.loads(Path(record['validated_artifact']).read_text())
        verify(installed['native']);verify(installed['package'])
        if record['native_sha256']!=installed['native']['sha256']:raise AssertionError('Old inference native identity changed')
        executions.append(job)
    return dict(baseline_files_preserved=len(baseline),executed_core_test_native_preserved=frozen_count,
                live_original_core_test_native_unchanged=live_unchanged,validated_live_changes=live_changes,
                validated_new_test=None if expected is None else dict(path='tests/test_proposal_graph_precision.py',
                    sha256=expected['tests/test_proposal_graph_precision.py']),old_execution_jobs_verified=executions,
                validated_correction_source=None if validated_source is None else str(validated_source),
                validated_correction_artifact=None if validated_artifact is None else str(validated_artifact))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--previous-task',type=Path,required=True)
    p.add_argument('--validated-source',type=Path);p.add_argument('--validated-artifact',type=Path)
    p.add_argument('--preservation-only',action='store_true')
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    preserved=preservation(a.task,a.task.parents[1],a.validated_source,a.validated_artifact)
    if a.preservation_only:
        (a.out/'report.json').write_text(json.dumps(dict(preserved,scope='Source/version preservation only; scientific audit remains pending'),indent=2))
        return
    summary=json.loads((a.task/'summary_v1/summary.json').read_text())
    evaluations=json.loads((a.task/'summary_v1/evaluation.json').read_text())
    root=a.task.parents[1]
    if summary['confidence_family']!=54 or len(summary['summary'])!=9 or len(summary['retained_historical'])!=45:
        raise AssertionError('Confidence family changed')
    previous=json.loads((a.previous_task/'null_summary_v1/evaluation.json').read_text())
    older=json.loads((root/'dev/confirmation_20261107/null_summary_v1/evaluation.json').read_text())
    if len(evaluations)!=1800:raise AssertionError('Missing evaluation decisions')
    index={(r['condition'],r['method'],r['replicate']):r for r in evaluations}
    if len(index)!=1800:raise AssertionError('Duplicated evaluation decisions')
    banks={role:[r['score'] for r in read(a.task/'shards_v1'/role/'records.jsonl')] for role in ('cal_nb','cal_state')}
    if any(len(bank)!=199 or not np.isfinite(bank).all() for bank in banks.values()):raise AssertionError('Bad calibration bank')
    input_count=decision_count=assignment_count=hashes=0;mechanism={}
    for role,n,condition,role_id,ci in [('cal_nb',199,'matched_nb',10,0),('cal_state',199,'correlated_state',11,1),
                                      ('eval_nb',200,'matched_nb',20,0),('eval_state',200,'correlated_state',20,1),
                                      ('eval_two',200,'two_clones',20,2)]:
        directory=a.task/'shards_v1'/role
        manifest=json.loads((directory/'run_manifest.json').read_text())
        if manifest['status']!='complete':raise AssertionError('Incomplete shard')
        for record in manifest['input']['files'].values():verify(record);hashes+=1
        current=read(directory/'records.jsonl');old=read(a.previous_task/'null_shards_v1/exact_matched'/role/'records.jsonl')
        if len(current)!=n or len(old)!=n:raise AssertionError('Wrong paired population size')
        changed=0;exact_calls=0
        for i,(r,before) in enumerate(zip(current,old)):
            if r['seed_key']!=[20261107,role_id,ci,i] or r['input_sha256']!=before['input_sha256']:
                raise AssertionError('Seed or counts/nuisance/reset changed')
            input_count+=1
            if len(r['initialization_trace'])!=1:raise AssertionError('Missing/extra initialization')
            init=r['initialization_trace'][0]
            if init['changed']:
                if init['input_groups']<=2 or init['output_groups']!=2:raise AssertionError('Invalid pooling intervention')
                changed+=1
            elif init['input_groups']!=init['output_groups']:raise AssertionError('Unrecorded group change')
            verify(r['assignments'])
            with np.load(r['assignments']['path']) as z:
                labels=z['raw_MAP'];R=z['raw_R']
                if R.shape[0]!=600 or not np.isfinite(R).all():raise AssertionError('Invalid membership array')
                np.testing.assert_array_equal(labels,R.argmax(axis=1))
                sizes=np.bincount(labels,minlength=R.shape[1]);eligible=np.flatnonzero(sizes[1:]>=30)+1
                np.testing.assert_array_equal(eligible,z['eligible_tumour_groups'])
                if len(eligible)!=r['raw_candidate_count']:raise AssertionError('Eligible group count differs')
            assignment_count+=1
            pair_scores=[];allowances=[]
            for d in r['final_pair_diagnostics']:
                values=[abs(s['z']) for s in d['segments']]
                value=math.fsum(values);pair_scores.append(value)
                # Positive-term NumPy sum versus independently rounded fsum.
                n_terms=len(values);eps_sum=np.finfo(float).eps
                allowances.append((n_terms*eps_sum/(1-n_terms*eps_sum)+2*eps_sum)*value)
            if abs(max(pair_scores,default=0.)-r['score'])>max(allowances,default=0.):
                raise AssertionError('Final statistic differs from retained segments')
            if not np.isfinite(r['score']):raise AssertionError('Nonfinite statistic')
            exact_calls+=sum(len(e.get('calls',[])) for e in r['stage_trace'] if e['stage']=='exact_confirmation')
            if role.startswith('eval_'):
                for method,bank in [('raw',None),('NB_calibrated_gate',banks['cal_nb']),('state_calibrated_gate',banks['cal_state'])]:
                    probability=None if bank is None else (1+sum(v>=r['score'] for v in bank))/200
                    event=bool(r['raw_candidate_count']>1 and (bank is None or r['score']>0 and probability<=.05))
                    saved=index[(condition,method,i)]
                    if saved['event']!=event or saved['rank_p']!=probability or saved['score']!=r['score']:
                        raise AssertionError('Independent acceptance mismatch')
                    decision_count+=1
        mechanism[role]=dict(pooled_initializations=changed,exact_calls=exact_calls)
    checks=[];alpha=.05/54;eps=np.finfo(float).eps
    for historical,group in [(False,summary['summary']),(True,summary['retained_historical'])]:
        population=previous+older if historical else evaluations
        for r in group:
            selected=[x for x in population if all(x[k]==r[k] for k in ('stage','condition','method'))]
            k,n=sum(x['event'] for x in selected),len(selected)
            if n!=r['n'] or k!=r['events']:raise AssertionError('Endpoint count mismatch')
            if r['direction']=='upper':
                ref=1. if k==n else optimize.brentq(lambda x:stats.binom.cdf(k,n,x)-alpha,0.,1.,xtol=eps,rtol=4*eps)
            else:
                ref=0. if not k else optimize.brentq(lambda x:stats.binom.sf(k-1,n,x)-alpha,0.,1.,xtol=eps,rtol=4*eps)
            discrepancy=abs(ref-r['simultaneous_bound'])
            if discrepancy>64*eps:raise ArithmeticError('Independent binomial-bound disagreement')
            checks.append(dict(stage=r['stage'],condition=r['condition'],method=r['method'],events=k,n=n,
                               bound=ref,absolute_discrepancy=discrepancy))
    if len(summary['paired'])!=18:raise AssertionError('Paired family changed')
    for r in summary['paired']:
        before={x['replicate']:x['event'] for x in previous if x['stage']==r['comparator'] and
                x['condition']==r['condition'] and x['method']==r['method']}
        after={x['replicate']:x['event'] for x in evaluations if x['condition']==r['condition'] and x['method']==r['method']}
        if len(before)!=200 or before.keys()!=after.keys():raise AssertionError('Unpaired decision population')
        lost=sum(before[i] and not after[i] for i in before);gained=sum(after[i] and not before[i] for i in before)
        if (lost,gained)!=(r['lost'],r['gained']):raise AssertionError('Paired discordance differs')
    report=dict(input_hashes_verified=hashes,paired_generated_inputs=input_count,assignments_verified=assignment_count,
                independent_decisions=decision_count,endpoint_counts_and_bounds=checks,paired_discordances=18,
                mechanism=mechanism,preservation=preserved,
                scope='Conditional simulation audit; does not supply biological replication or suitable truth')
    (a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)


if __name__=='__main__':main()
