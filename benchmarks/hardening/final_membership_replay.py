"""Validate returned-support repair on its exact exposed complete-workflow case."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import numpy as np
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name):
    spec = importlib.util.spec_from_file_location('final_support_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    old = args.task.parent/'phase8_rep4_trace_20261110'
    assert json.loads((old/'summary_v1/run_manifest.json').read_text())['status'] == 'complete'
    assert json.loads((old/'summary_v1/report.json').read_text())['exact_replay_gate_passed']
    tests = (args.validation/'pytest.txt').read_text()
    assert re.search(r'\b\d+ passed\b', tests)
    assert not re.search(r'\b\d+ (failed|errors?)\b', tests)
    jobs = json.loads((args.task/'jobs.json').read_text())
    baseline = next(r['job_id'] for r in jobs if r['name'] == 'baseline')
    counterexample = args.task/f'baseline_{baseline}/report.json'
    assert json.loads(counterexample.read_text())['expected_support_failures'] == 2
    manifest = json.loads((args.task/'source_after/source_manifest.json').read_text())
    installed = Path(clones.__file__)
    assert hashlib.sha256(installed.read_bytes()).hexdigest() == manifest['python/haplodose/clones.py']
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'REPLAY_DESIGN.md',
        baseline=counterexample, pytest=args.validation/'pytest.txt',
        installed_artifact=args.validation/'installed_artifact.json',
        original_trace=old/'trace_v1/report.json', original_summary=old/'summary_v1/report.json',
        original_results=old/'trace_v1/replay/results.json')
    run = start_run(args.out, inputs, dict(seed=202611098004, minimum_unchanged=30),
        dict(scope='Exposed complete-workflow correctness regression; no error-rate or biological claim'),
        mode='final_membership_support_replay')
    pilot = load('pooled_prevalence_pilot')
    argv = sys.argv
    try:
        sys.argv = [str(pilot.__file__), '--out', str(args.out/'replay'),
            '--task', str(args.task.parent/'pooled_prevalence_pilot_20261110'),
            '--condition', 'two_clone', '--start', '4', '--stop', '5']
        pilot.main()
        before = old/'trace_v1/replay/seed_202611098004'
        after = args.out/'replay/seed_202611098004'
        production = load('refinement_checked_replay').exact_outputs(before, after, arms=('production',))
        folds = 0
        for file in (before/'pooled_prevalence/folds').glob('*.npz'):
            with np.load(file) as a, np.load(after/'pooled_prevalence/folds'/file.name) as b:
                assert a.files == b.files
                for key in a.files:
                    np.testing.assert_array_equal(a[key], b[key])
            folds += 1
        assert folds == 1
        rows = json.loads((args.out/'replay/results.json').read_text())
        original = json.loads(inputs['original_results'].read_text())
        write_json(dict(results=rows, original_results=original,
            production_parity=production, count_fold_archives_exact=folds,
            scientific_criteria_changed=False, biological_specimens=0,
            independent_evaluation=False, false_split_control_established=False,
            promotion=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU installed corrected wheel with frozen complete phase8 workflow')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise
    finally:
        sys.argv = argv


if __name__ == '__main__':
    main()
