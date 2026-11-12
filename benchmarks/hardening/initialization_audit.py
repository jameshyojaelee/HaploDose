"""Independent saved-output checks for initialization tracing and its ablation."""
import argparse
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


def verify_file(record):
    if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest() != record['sha256']:
        raise AssertionError('Recorded input/output changed')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--rank-source', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    files = 0
    for directory in ('trace_v2','pool_mechanism_v1','pool_recovery_v1','pool_recovery_v1/panel'):
        manifest = json.loads((a.task/directory/'run_manifest.json').read_text())
        if manifest['status'] != 'complete': raise AssertionError('Incomplete required experiment')
        for record in manifest['input']['files'].values():
            verify_file(record)
            files += 1
    traces = json.loads((a.task/'trace_v2/traces.json').read_text())
    if len(traces) != 6: raise AssertionError('Missing replay')
    phases = []
    for row in traces:
        for key in ('frozen_assignment','replay_assignment'): verify_file(row[key])
        with np.load(row['frozen_assignment']['path']) as old, np.load(row['replay_assignment']['path']) as new:
            for name in old.files: np.testing.assert_array_equal(old[name],new[name])
        em = [r for r in row['phases'] if r['phase']=='em']
        splits = [r for r in row['phases'] if r['phase']=='try_splits']
        phases.append(dict(stage=row['stage'],role=row['role'],initial_groups=em[0]['input_groups']-1,
                           raw_groups=row['raw_candidate_count'],accepted_splits=sum(r['accepted'] for r in splits),
                           matched_merges=sum(r['action']=='merge_matched' for r in row['model_log'])))
    spec = importlib.util.spec_from_file_location('independent_generator',a.rank_source/'global_rank_generator.py')
    gen = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gen)
    records = json.loads((a.task/'pool_mechanism_v1/records.json').read_text())
    pools = json.loads((a.task/'pool_mechanism_v1/initialization_trace.json').read_text())
    if len(records)!=3 or len(pools)!=3: raise AssertionError('Incomplete mechanism panel')
    conservation = []
    for row,pool in zip(records,pools):
        Y,_,_,_=gen.generate(np.random.default_rng(np.random.SeedSequence(row['seed_key'])),row['condition'],600,240)
        verify_file(row['assignments'])
        if pool['changed']:
            # For these preselected cases every initial group is tumour-like.
            # Check that premise independently against the exact replay trace.
            reference = next(r for r in traces if r['stage']=='exact_matched' and r['condition']==row['condition'])
            first = reference['phases'][0]
            members = first['initial_agglomerated_members']
            if first['initial_tumour_indices'] != list(range(len(members))):
                raise AssertionError('Cannot apply all-cell conservation reference')
            np.testing.assert_array_equal(Y.sum(axis=0,dtype=np.float64),pool['pooled_counts'])
            if float(Y.sum(dtype=np.float64)) != pool['pooled_depth']:
                raise AssertionError('Pooled depth not conserved')
            L=Y.sum(axis=1,dtype=np.float64); x=np.log(L)-np.log(L).mean()
            expected=np.array([sum(float(l*t) for l,t in zip(L,x)),sum(float(l*t*t) for l,t in zip(L,x))])
            allowance=4*len(L)*np.finfo(float).eps*np.array([sum(abs(L*x)),sum(abs(L*x*x))])
            if np.any(abs(expected-pool['pooled_depth_moments'])>allowance):
                raise AssertionError('Pooled depth moments disagree with scalar reference')
        conservation.append(dict(condition=row['condition'],changed=pool['changed'],groups=row['raw_candidate_count']))
    panel=a.task/'pool_recovery_v1/panel'
    summary=json.loads((panel/'summary.json').read_text())
    if [r['seed'] for r in summary['cases']] != [20261031,20261107,20261108,20261109,20261110]:
        raise AssertionError('Recovery panel changed')
    recovery=[]
    for row in summary['cases']:
        seed=panel/f"seed_{row['seed']}"
        cells=pd.read_csv(seed/'assertion_view/cells.tsv.gz',sep='\t').to_dict('records')
        differences=pd.read_csv(seed/'fit/differences.tsv',sep='\t').to_dict('records')
        normal=[r for r in cells if r['specimen']=='S1' and r['truth_clone']=='normal']
        tumour=[r for r in cells if r['specimen']=='S1' and r['truth_clone']!='normal']
        sub=[r for r in tumour if r['truth_clone'] in ('clone1','clone2') and r['malignant_call']]
        truth=sorted({r['truth_clone'] for r in sub}); labels=sorted({r['clone'] for r in sub})
        counts=Counter((r['truth_clone'],r['clone']) for r in sub)
        diagonal=lambda names:sum(counts[(t,c)] for t,c in zip(names,labels))
        accuracy=max(diagonal(truth),diagonal(truth[::-1]))/len(sub) if sub else None
        chrom=sorted({r['chromosome'] for r in differences if r['specimen']=='S1'
                      and not r['clone_a'].endswith('_c0') and not r['clone_b'].endswith('_c0')})
        values=dict(normal_fpr=sum(r['malignant_call'] for r in normal)/len(normal),
                    tumour_tpr=sum(r['malignant_call'] for r in tumour)/len(tumour),
                    s1_two_groups=len(labels),s2_one_group=len({r['clone'] for r in cells
                    if r['specimen']=='S2' and r['malignant_call']}),assignment_accuracy=accuracy,
                    difference_chromosomes=chrom)
        expected=dict(normal_fpr=values['normal_fpr']<.1,tumour_tpr=values['tumour_tpr']>.9,
                      s1_two_groups=len(labels)==2,s2_one_group=values['s2_one_group']==1,
                      assignment_accuracy=accuracy is not None and accuracy>.65,
                      difference_chromosomes=chrom==['chr12'])
        for key,value in values.items():
            if value != row['criteria'][key]['value'] or expected[key] != row['criteria'][key]['passed']:
                raise AssertionError(f'Independent recovery recount differs: {key}')
        if all(expected.values()) != row['original_assertion_passed']:
            raise AssertionError('Original assertion and exhaustive criteria disagree')
        recovery.append(dict(seed=row['seed'],passed=all(expected.values()),values=values))
    base=json.loads((a.task/'baseline/source_manifest.json').read_text()); unchanged=0
    for relative,digest in base.items():
        if hashlib.sha256((a.task/'baseline'/relative).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Initial baseline changed')
        if relative.startswith(('python/','tests/','crates/')):
            if hashlib.sha256((a.task.parents[1]/relative).read_bytes()).hexdigest()!=digest:
                raise AssertionError('Production/test/native source changed')
            unchanged+=1
    report=dict(recorded_input_hashes_verified=files,identity_replays=6,phases=phases,
                independent_pooling_conservation=conservation,recovery=recovery,
                baseline_files_preserved=len(base),production_test_native_unchanged=unchanged,
                scope='Mechanism and original recovery checks; no population error or biological validation')
    (a.out/'report.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report),flush=True)


if __name__=='__main__':main()
