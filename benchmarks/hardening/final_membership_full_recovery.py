"""Check the combined refinement/support correction on all original cases."""
import argparse
from contextlib import contextmanager
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
    spec = importlib.util.spec_from_file_location('final_support_full_'+name,
        Path(__file__).with_name(name+'.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def observe_support(records):
    """Check returned support without replacing a numerical method or its caller."""
    if sys.getprofile() is not None:
        raise RuntimeError('Support observer requires an otherwise free profile hook')
    code = clones.SpecimenModel.em.__code__

    def observe(frame, event, result):
        if event != 'return' or frame.f_code is not code:
            return
        if result is None:
            raise RuntimeError('Observed refinement failed')
        model = frame.f_locals['self']
        doses, pi, R, LL = result[:4]
        mass = np.asarray(R).sum(axis=0)
        record = dict(source=str(model.rows.Y.filename), components=len(doses),
            minimum=model.params.min_cells, mass=mass.tolist(),
            shapes_aligned=bool(R.shape[1] == len(doses) == len(pi) == LL.shape[1]),
            all_supported=bool(np.isfinite(mass).all()
                and np.all(mass[1:] >= model.params.min_cells)))
        records.append(record)
        assert record['shapes_aligned'] and record['all_supported']

    sys.setprofile(observe)
    try:
        yield
    finally:
        sys.setprofile(None)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prior = args.task.parent/'refinement_checked_split_20261110'
    inputs = dict(script=Path(__file__), design=args.task/'FULL_COMBINED_DESIGN.md',
        pytest=args.validation/'pytest.txt',
        corrected_case=args.task/'replay_v1/report.json',
        corrected_case_manifest=args.task/'replay_v1/run_manifest.json',
        original_full=prior/'full_recovery_v1/report.json',
        original_full_manifest=prior/'full_recovery_v1/run_manifest.json',
        wrapper=Path(__file__).with_name('refinement_checked_split.py'),
        runner=Path(__file__).with_name('refinement_checked_full_recovery.py'))
    tests = inputs['pytest'].read_text()
    assert re.search(r'\b\d+ passed\b', tests)
    assert not re.search(r'\b\d+ (failed|errors?)\b', tests)
    assert json.loads(inputs['corrected_case_manifest'].read_text())['status'] == 'complete'
    assert json.loads(inputs['original_full_manifest'].read_text())['status'] == 'complete'
    old = json.loads(inputs['original_full'].read_text())
    assert old['all_original_criteria_passed'] and old['recovered'] == 5
    corrected = json.loads(inputs['corrected_case'].read_text())
    assert corrected['production_parity']['exact']
    rows = [r for r in corrected['results'] if r['arm'] == 'pooled_prevalence']
    assert len(rows) == 1
    assert rows[0]['reported_groups'] == 2 and rows[0]['exact_chr12_recovery']
    assert rows[0]['assignment_n'] == 600 and rows[0]['assignment_correct'] >= 438
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(original_cases=5, original_criteria=30),
        dict(scope='Combined corrected-wheel/refinement original regression; biological n=0'),
        mode='final_membership_full_recovery')
    runner = load('refinement_checked_full_recovery')
    source = runner.load('pooled_prevalence_workflow')
    original = source.safeguard
    records = []

    @contextmanager
    def context(*a, **kw):
        with original(*a, **kw):
            with observe_support(records):
                yield

    source.safeguard = context
    original_loader = runner.load

    def select(name):
        # The original full runner's loader is intentionally not cached. Route
        # its workflow request to the SAME observed module prepared above.
        return source if name == 'pooled_prevalence_workflow' else original_loader(name)

    runner.load = select
    argv = sys.argv
    try:
        # Preserve every original gate and criterion. The runner itself installs
        # the byte-checked refinement constraint around this observed context.
        sys.argv = [str(runner.__file__), '--task', str(prior),
                    '--out', str(args.out/'full')]
        runner.main()
        report = json.loads((args.out/'full/report.json').read_text())
        assert records and all(r['all_supported'] and r['shapes_aligned'] for r in records)
        write_json(dict(recovered=report['recovered'], cases=report['cases'],
            all_original_criteria_passed=report['all_original_criteria_passed'],
            returned_support_checks=records, full_workflow_including_merge=True,
            scientific_criteria_changed=False, biological_specimens=0,
            independent_evaluation=False, false_split_control_established=False,
            promotion=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU combined corrected wheel/refinement full original recovery')
    except BaseException as error:
        write_json(records, args.out/'support_records_failed.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        source.safeguard, sys.argv = original, argv
        runner.load = original_loader


if __name__ == '__main__':
    main()
