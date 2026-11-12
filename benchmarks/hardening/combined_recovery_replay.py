"""Test the combined support/refinement fixes on four exposed phase8 cases."""
import argparse
from contextlib import contextmanager
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


CASES = (('two_clone', 1), ('two_clone', 2), ('two_clone', 4), ('state', 3))


def load(name):
    spec = importlib.util.spec_from_file_location('combined_replay_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_control(before, after):
    """Check complete production output sets, generated inputs and gene folds."""
    counts = dict(arrays=0, tables=0, folds=0, input_arrays=0)
    eligible = lambda p: p.suffix == '.npy' or p.name.endswith(('.tsv', '.tsv.gz'))
    a, b = before/'production/fit', after/'production/fit'
    names = {p.name for p in a.iterdir() if eligible(p)}
    assert names == {p.name for p in b.iterdir() if eligible(p)}
    for name in sorted(names):
        if name.endswith('.npy'):
            np.testing.assert_array_equal(np.load(a/name), np.load(b/name))
            counts['arrays'] += 1
        else:
            pd.testing.assert_frame_equal(pd.read_csv(a/name, sep='\t'),
                pd.read_csv(b/name, sep='\t'), check_exact=True)
            counts['tables'] += 1
    for relative in ('generation.npz', 'pooled_prevalence/folds'):
        if relative.endswith('.npz'):
            pairs = [(before/relative, after/relative)]
        else:
            a, b = before/relative, after/relative
            names = {p.name for p in a.glob('*.npz')}
            assert names == {p.name for p in b.glob('*.npz')}
            pairs = [(a/name, b/name) for name in sorted(names)]
        for a, b in pairs:
            with np.load(a) as x, np.load(b) as y:
                assert x.files == y.files
                for key in x.files:
                    np.testing.assert_array_equal(x[key], y[key])
            counts['folds' if relative.endswith('folds') else 'input_arrays'] += 1
    with h5py.File(before/'input/counts.h5') as x, h5py.File(after/'input/counts.h5') as y:
        for key in ('X/data', 'X/indices', 'X/indptr', 'genes', 'cell_ids'):
            np.testing.assert_array_equal(x[key][:], y[key][:])
            counts['input_arrays'] += 1
    for name in ('model_cells.tsv', 'truth_cells.tsv', 'positions.tsv'):
        pd.testing.assert_frame_equal(pd.read_csv(before/'input'/name, sep='\t'),
            pd.read_csv(after/'input'/name, sep='\t'), check_exact=True)
    assert counts == dict(arrays=2, tables=14, folds=1, input_arrays=6)
    return dict(counts, exact=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--condition', choices=('two_clone', 'state'), required=True)
    parser.add_argument('--replicate', type=int, required=True)
    args = parser.parse_args()
    assert (args.condition, args.replicate) in CASES
    # Check scientific reports before calling the original generator or pipeline.
    integrated = args.task/'combined_full_v1'
    report = json.loads((integrated/'report.json').read_text())
    assert json.loads((integrated/'run_manifest.json').read_text())['status'] == 'complete'
    assert report['all_original_criteria_passed'] and report['recovered'] == 5
    assert len(report['cases']) == 5
    assert all(row['original_assertion_passed'] and len(row['criteria']) == 6
        and all(value['passed'] for value in row['criteria'].values())
        for row in report['cases'])
    assert report['returned_support_checks'] and all(
        r['all_supported'] and r['shapes_aligned'] for r in report['returned_support_checks'])
    expected = json.loads((args.task/'source_after/source_manifest.json').read_text())
    assert hashlib.sha256(Path(clones.__file__).read_bytes()).hexdigest() == expected['python/haplodose/clones.py']
    seed = 202611098000+args.replicate
    old = args.task.parent/'pooled_prevalence_pilot_20261110'/(
        'pilot_'+args.condition+'_00_v1')/('seed_'+str(seed))
    assert json.loads((old/'production/fit/run_manifest.json').read_text())['status'] == 'complete'
    assert json.loads((old/'pooled_prevalence/fit/run_manifest.json').read_text())['status'] == 'complete'
    method, observer = load('refinement_checked_split'), load('final_membership_full_recovery')
    sham = args.task.parent/'refinement_checked_split_20261110/source_v1/benchmarks/hardening/refinement_checked_split.py'
    assert Path(method.__file__).read_bytes() == sham.read_bytes()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'PHASE8_COMBINED_DESIGN.md',
        combined_gate=integrated/'report.json', combined_manifest=integrated/'run_manifest.json',
        wrapper=Path(method.__file__), support_observer=Path(observer.__file__),
        installed_artifact=args.task/'validation_23169639/installed_artifact.json',
        original_results=old.parent/'results.json',
        pilot_manifest=old.parent/'run_manifest.json')
    run = start_run(args.out, inputs, dict(seed=seed, condition=args.condition,
        replicate=args.replicate, full_workflow=True, minimum_unchanged=30),
        dict(scope='Four locked exposed regressions; no independent error/power or biological claim'),
        mode='combined_recovery_replay')
    pilot = load('pooled_prevalence_pilot')
    guard = pilot.load('pooled_prevalence_workflow')
    original, argv = guard.safeguard, sys.argv
    previews, support = [], []

    @contextmanager
    def context(*a, **kw):
        with original(*a, **kw):
            with observer.observe_support(support):
                with method.checked_splits(previews, apply=True):
                    yield

    guard.safeguard = context
    try:
        sys.argv = [str(pilot.__file__), '--task', str(old.parents[1]),
            '--out', str(args.out/'replay'), '--condition', args.condition,
            '--start', str(args.replicate), '--stop', str(args.replicate+1)]
        pilot.main()
        after = args.out/'replay'/('seed_'+str(seed))
        parity = exact_control(old, after)
        assert previews and support
        assert all(r['all_supported'] and r['shapes_aligned'] for r in support)
        rows = json.loads((args.out/'replay/results.json').read_text())
        assert len(rows) == 2
        candidate = next(r for r in rows if r['arm'] == 'pooled_prevalence')
        old_rows = json.loads(inputs['original_results'].read_text())
        baseline = [r for r in old_rows if r['replicate'] == args.replicate]
        assert len(baseline) == 2
        passed = (not candidate['null_error'] if args.condition == 'state'
                  else candidate['exact_chr12_recovery'])
        retained = None
        if args.condition == 'two_clone' and args.replicate in (1, 4):
            required = {1: 445, 4: 438}[args.replicate]
            retained = bool(passed and candidate['assignment_n'] == 600
                            and candidate['assignment_correct'] >= required)
        write_json(dict(results=rows, original_results=baseline,
            original_recovery_or_null_criterion_passed=bool(passed),
            constituent_benefit_retained=retained, production_and_input_parity=parity,
            preview_records=previews, returned_support_checks=support,
            scientific_criteria_changed=False, biological_specimens=0,
            independent_evaluation=False, false_split_control_established=False,
            promotion=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU corrected wheel plus frozen refinement constraint')
    except BaseException as error:
        write_json(dict(previews=previews, support=support), args.out/'records_failed.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        guard.safeguard, sys.argv = original, argv


if __name__ == '__main__':
    main()
