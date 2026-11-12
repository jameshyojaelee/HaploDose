"""Test a consistent profile estimator after merging, on locked development inputs.

Research only. All merge decisions and original recovery assertions stay intact.
An observed improvement is neither genomic identification nor error calibration.
"""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@contextmanager
def consistent_merge(records, directory, *, apply):
    original = clones.SpecimenModel.profiles
    merge_code = clones.SpecimenModel.merge_unsupported.__code__

    def profiles(self, Y, L, XL):
        caller = sys._getframe(1)
        if caller.f_code is not merge_code:
            return original(self, Y, L, XL)
        if self.alleles is not None or self.graph is not None:
            raise ValueError('This research probe is expression-only and nonspatial')
        assert caller.f_locals['statistical_merge']
        R = caller.f_locals['R']
        assert len(Y) == R.shape[1] == len(caller.f_locals['doses'])
        if apply:
            # shared_profiles calls profiles internally. Those calls originate
            # in shared_profiles, so they take the unchanged branch above.
            result = self.shared_profiles(Y, L, XL, labels=R.argmax(axis=1))
        else:
            result = original(self, Y, L, XL)
        path = directory/f'merge_{len(records):04d}.npz'
        np.savez_compressed(path, Y=Y, L=L, XL=XL, R=R,
            dose=result['dose'], own_dose=result.get('dose_own', result['dose']),
            rows=self.rows.rows, reset=self.reset,
            supplied_doses=np.asarray(caller.f_locals['doses']))
        records.append(dict(arrays=str(path), apply=apply, groups=len(Y)-1,
            counts_source=str(self.rows.Y.filename),
            estimator='shared_profiles' if apply else 'profiles',
            merge_criteria_changed=False))
        return result

    clones.SpecimenModel.profiles = profiles
    try:
        yield
    finally:
        clones.SpecimenModel.profiles = original


def exact_outputs(before, after):
    records = []
    eligible = lambda p: p.suffix == '.npy' or p.name.endswith(('.tsv', '.tsv.gz'))
    for case in sorted(before.glob('seed_*')):
        old, new = case/'fit', after/case.name/'fit'
        names = sorted(p.name for p in old.iterdir() if eligible(p))
        assert names == sorted(p.name for p in new.iterdir() if eligible(p))
        assert names
        for name in names:
            if name.endswith('.npy'):
                np.testing.assert_array_equal(np.load(old/name), np.load(new/name))
            else:
                pd.testing.assert_frame_equal(pd.read_csv(old/name, sep='\t'),
                    pd.read_csv(new/name, sep='\t'), check_exact=True)
        records.append(dict(case=case.name, outputs_exact=names))
    assert len(records) == 5
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.task.parent.parent
    recovery = root/'dev/recovery_20261107'
    assertions = recovery/'baseline/tests/test_split_mode.py'
    panel = recovery/'recovery_panel_v1'
    runner_path = Path(__file__).with_name('membership_stage_panel.py')
    runner = load('consistent_merge_fixed_panel', runner_path)
    # Resource placement is allowed only after BOTH executed fixture controls.
    control_root = root/'dev/production_support_fresh_20261110'
    for version in ('historical', 'corrected'):
        control = control_root/f'node003_{version}_v1'
        assert json.loads((control/'run_manifest.json').read_text())['status'] == 'complete'
        assert json.loads((control/'report.json').read_text())['passed']
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, dict(script=Path(__file__), runner=runner_path,
        design=args.task/'DESIGN.md', assertions=assertions,
        panel_manifest=panel/'manifest.json'),
        dict(seeds=[20261031, 20261107, 20261108, 20261109, 20261110],
             arms=['baseline', 'sham', 'candidate'], original_criteria=30),
        dict(scope='Exposed causal recovery experiment; biological n=0; no calibration claim'),
        mode='shared_profile_merge_probe')
    reports, traces, argv = {}, {}, sys.argv
    original_method = clones.SpecimenModel.profiles
    try:
        for arm in ('baseline', 'sham', 'candidate'):
            directory = args.out/f'{arm}_arrays'; directory.mkdir()
            records = []
            sys.argv = [str(runner_path), '--out', str(args.out/arm),
                '--panel-root', str(panel), '--assertion-source', str(assertions), '--stage', 'hard']
            if arm == 'baseline':
                runner.main()
            else:
                with consistent_merge(records, directory, apply=arm == 'candidate'):
                    runner.main()
                assert records, 'An unexercised merge intervention proves nothing'
            assert clones.SpecimenModel.profiles is original_method
            assert json.loads((args.out/arm/'run_manifest.json').read_text())['status'] == 'complete'
            reports[arm] = json.loads((args.out/arm/'summary.json').read_text())
            traces[arm] = records
            write_json(records, args.out/f'{arm}_merge_records.json')
            # Require exact sham output before running the intervention.
            if arm == 'sham':
                sham = exact_outputs(args.out/'baseline', args.out/'sham')
                assert reports['baseline']['cases'] and reports['sham']['cases']
                for old, new in zip(reports['baseline']['cases'], reports['sham']['cases']):
                    for key in ('seed', 'original_assertion_passed', 'criteria', 's3_groups'):
                        assert old[key] == new[key], key
                write_json(dict(passed=True, cases=sham), args.out/'SHAM_LOCKED.json')
        old, new = reports['baseline']['cases'], reports['candidate']['cases']
        assert len(old) == len(new) == 5
        comparisons = []
        for before, after in zip(old, new):
            assert before['seed'] == after['seed']
            comparisons.append(dict(seed=before['seed'], before=before, after=after,
                recovered=not before['original_assertion_passed'] and after['original_assertion_passed'],
                lost=before['original_assertion_passed'] and not after['original_assertion_passed'],
                criterion_losses=[k for k, v in before['criteria'].items()
                    if v['passed'] and not after['criteria'][k]['passed']],
                additional_s3_groups=after['s3_groups'] > before['s3_groups']))
        advance = (any(r['recovered'] for r in comparisons)
            and not any(r['lost'] or r['criterion_losses'] or r['additional_s3_groups']
                        for r in comparisons))
        write_json(dict(sham_exact=True, comparisons=comparisons,
            baseline_recovered=sum(r['original_assertion_passed'] for r in old),
            candidate_recovered=sum(r['original_assertion_passed'] for r in new),
            prespecified_continuation_passed=advance, merge_calls=traces,
            scientific_criteria_changed=False, independent_evaluation=False,
            biological_specimens=0, false_split_control_established=False, promotion=False),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed wheel; shared-profile merge research only')
    except BaseException as error:
        write_json(traces, args.out/'merge_records_failed.json')
        finish_run(args.out, run, error=error)
        raise
    finally:
        sys.argv = argv
        assert clones.SpecimenModel.profiles is original_method


if __name__ == '__main__':
    main()
