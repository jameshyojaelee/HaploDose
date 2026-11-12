"""Label-invariant genomic coassignment diagnostic on frozen RNA fits."""
import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.provenance import start_run, finish_run, file_record


def coassignment(raw):
    raw = np.asarray(raw)
    if raw.ndim != 2 or raw.shape[1] == 0 or len(raw) < 2:
        raise ValueError('At least two cells and one component required')
    if not np.isfinite(raw).all() or np.any(raw < 0):
        raise ValueError('Finite nonnegative responsibilities required')
    eps = np.finfo(raw.dtype).eps
    rows = raw.astype(float).sum(axis=1)
    k = raw.shape[1]
    budget = (k+2)*eps/(1-(k+2)*eps)
    if (rows <= 0).any() or np.any(np.abs(rows-1) > budget):
        raise ValueError('Responsibilities are not probability vectors')
    R = raw.astype(float)/rows[:, None]
    i, j = np.triu_indices(len(R), 1)
    q = np.einsum('ij,ij->i', R[i], R[j])
    # Roundoff outside [0,1] is a failed numerical check, not silently clipped.
    if not ((q >= 0) & (q <= 1)).all():
        raise ArithmeticError('Unresolved coassignment probability range')
    return i, j, q, float(np.max(np.abs(R-raw)))


def metrics(q, same):
    confident = (q >= .95) | (q <= .05)
    correct = (q >= .5) == same
    n = len(q); selected = int(confident.sum())
    return dict(pairs=n, brier=float(np.mean((q-same)**2)),
        one_population_brier=float(np.mean((1-same)**2)), uninformative_brier=.25,
        confident_pairs=selected, confident_errors=int(np.sum(confident & ~correct)),
        confident_coverage=selected/n,
        confident_conditional_error=float(np.mean(~correct[confident])) if selected else np.nan,
        confident_unconditional_error=float(np.mean(confident & ~correct)),
        confident_unconditional_correct=float(np.mean(confident & correct)),
        high_confidence_separations=int(np.sum(q <= .05)),
        high_confidence_same=int(np.sum(q >= .95)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True); p.add_argument('--task',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(exist_ok=False)
    source=a.task.parents[1]/'dev/joint_fresh_null_20261108'
    inputs=dict(script=Path(__file__),design=a.task/'DESIGN.md',prior_audit=source/'audit_v1/run_manifest.json')
    assert json.loads(inputs['prior_audit'].read_text())['status']=='complete'
    paths=[]
    for arm,role,start in itertools.product(('joint_proposal','joint_full'),
            ('eval_nb','eval_state','eval_two'),(0,50,100,150)):
        directory=source/'chunks_v1'/arm/f'{role}_{start}'
        manifest=json.loads((directory/'run_manifest.json').read_text()); assert manifest['status']=='complete'
        record=manifest['output']['files']['records.jsonl']
        assert file_record(record['path'])['sha256']==record['sha256']
        inputs[f'{arm}:{role}:{start}:manifest']=directory/'run_manifest.json'
        inputs[f'{arm}:{role}:{start}:records']=directory/'records.jsonl'
        paths.append((arm,role,start,Path(record['path'])))
    run=start_run(a.out,inputs,dict(cases=1200,refits=0),
        dict(scope='Descriptive pairwise membership confidence; no biological calibration or promotion'),
        mode='clone_membership_confidence')
    try:
        rows=[]; reliability=[]; inventory=[]
        for arm,role,start,path in paths:
            records=[json.loads(s) for s in path.read_text().splitlines()]
            assert [r['replicate'] for r in records]==list(range(start,start+50))
            for record in records:
                assert record['seed_key']==[20261107,120,('eval_nb','eval_state','eval_two').index(role),record['replicate']]
                saved=record['assignments']; assert file_record(saved['path'])['sha256']==saved['sha256']
                with np.load(saved['path']) as data:
                    truth=data['simulated_CNA_labels']; raw=data['raw_R']; labels=data['raw_MAP']
                assert raw.shape[0]==len(truth)==600
                assert np.array_equal(raw.argmax(axis=1),labels)
                if role!='eval_two': assert len(np.unique(truth))==1
                else: assert np.array_equal(truth,np.repeat([0,1],300))
                i,j,q,change=coassignment(raw); same=truth[i]==truth[j]
                key=dict(arm=arm,role=role,replicate=record['replicate'])
                rows.append(dict(**key,**metrics(q,same),normalization_change=change))
                decile=np.minimum((q*10).astype(int),9)
                for d in range(10):
                    mask=decile==d; n=int(mask.sum())
                    reliability.append(dict(**key,decile=d,pairs=n,
                        probability_sum=float(q[mask].sum()), same_pairs=int(same[mask].sum())))
                inventory.append(dict(**key,**saved))
        frame=pd.DataFrame(rows); frame.to_csv(a.out/'cases.tsv',sep='\t',index=False)
        rel=pd.DataFrame(reliability); rel.to_csv(a.out/'reliability.tsv',sep='\t',index=False)
        summary=[]
        for (arm,role),f in frame.groupby(['arm','role']):
            assert len(f)==200
            summary.append(dict(arm=arm,role=role,specimens=len(f),
                **{c:float(f[c].mean()) for c in ('brier','one_population_brier','uninformative_brier',
                    'confident_coverage','confident_unconditional_error','confident_unconditional_correct')},
                specimens_with_confident_pairs=int(f.confident_pairs.gt(0).sum()),
                macro_conditional_error=float(f.confident_conditional_error.mean()),
                specimens_with_confident_separation=int(f.high_confidence_separations.gt(0).sum()),
                total_confident_pairs=int(f.confident_pairs.sum()), total_confident_errors=int(f.confident_errors.sum())))
        pd.DataFrame(summary).to_csv(a.out/'summary.tsv',sep='\t',index=False)
        (a.out/'assignments.json').write_text(json.dumps(inventory,indent=2)+'\n')
        (a.out/'report.json').write_text(json.dumps(dict(cases=1200,biological_specimens=0,
            pairs_per_specimen=179700,summary=summary,max_normalization_change=float(frame.normalization_change.max()),
            no_new_calibration_claim=True,original_error_power_and_recovery_criteria_unchanged=True),indent=2)+'\n')
        for record in run['input']['files'].values(): assert file_record(record['path'])['sha256']==record['sha256']
    except BaseException as error:
        finish_run(a.out,run,error=error); raise
    finish_run(a.out,run,backend='CPU NumPy, no inference')


if __name__=='__main__': main()
