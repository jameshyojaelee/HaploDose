"""Regression-only application of the previously frozen production specimen gate.

No new calibration, inference, cutoff selection or prospective validation.
"""
import argparse
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import file_record, start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def checked_gate(gate, folder, banks):
    assert json.loads((folder/'run_manifest.json').read_text())['status'] == 'complete'
    value = gate.statistic(folder)
    decision = gate.decision(value, banks)
    independent = [Fraction(1+sum(float(s) >= value for s in bank), 87) for bank in banks]
    assert decision['passed'] == (max(independent) <= Fraction(1, 87))
    assert decision['passed'] == (value > max(max(bank) for bank in banks))
    assert decision['inclusive_rank_p'] == [float(p) for p in independent]
    return dict(score=value, accepted=decision['passed'], ranks=decision['inclusive_rank_p'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = args.task.parent.parent
    prior = root/'dev/specimen_workflow_gate_20261109'
    fresh = root/'dev/production_support_fresh_20261110'
    development = root/'dev/final_membership_support_20261110'
    original = root/'dev/shared_profile_merge_20261110/probe_v1/baseline'
    gate_path = prior/'source_v1/benchmarks/hardening/specimen_workflow_gate.py'
    gate = load('unchanged_historical_gate', gate_path)
    scorer = load('original_six_criteria', Path(__file__).with_name('recovery_baseline.py'))
    assertion_path = root/'dev/recovery_20261107/baseline/tests/test_split_mode.py'
    tests = load('unchanged_original_assertion', assertion_path)
    bank_manifest = json.loads((prior/'bank_v1/run_manifest.json').read_text())
    assert bank_manifest['status'] == 'complete'
    bank_record = bank_manifest['output']['files']['banks.json']
    assert file_record(bank_record['path'])['sha256'] == bank_record['sha256']
    assert file_record(gate_path)['sha256'] == bank_manifest['input']['files']['script']['sha256']
    banks = json.loads(Path(bank_record['path']).read_text())['production']
    banks = [banks['matched'], banks['state']]
    assert all(len(bank) == 86 for bank in banks)
    assert gate.decision(3., [np.full(86, 2.), np.full(86, 1.)])['passed']
    assert not gate.decision(2., [np.full(86, 2.), np.full(86, 1.)])['passed']
    assert not gate.decision(0., [np.zeros(86), np.zeros(86)])['passed']
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', gate=gate_path,
        bank=Path(bank_record['path']), bank_manifest=prior/'bank_v1/run_manifest.json',
        frozen_primary=fresh/'summary_v1/report.json', raw_scores=fresh/'summary_v1/scores.tsv',
        raw_manifest=fresh/'summary_v1/run_manifest.json',
        original_manifest=original/'run_manifest.json', original_summary=original/'summary.json',
        scorer=Path(scorer.__file__), assertions=assertion_path)
    for condition in ('matched', 'state', 'two_clone'):
        for start in (0, 100):
            inputs[f'{condition}{start}_lock'] = fresh/f'corrected_{condition}_{start:03d}_v1/FITS_LOCKED.json'
        for start in (0, 5):
            inputs[f'development_{condition}{start}'] = development/f'production_regression_{condition}_{start:02d}_v1/report.json'
    assert json.loads(inputs['raw_manifest'].read_text())['status'] == 'complete'
    assert json.loads(inputs['original_manifest'].read_text())['status'] == 'complete'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(bank_size=86, rank_cutoff='1/87',
        statistic='UNCHANGED maximum final pair sum absolute z', arms=['raw', 'frozen_production_gate']),
        dict(scope='Previously exposed regression; not new primary testing or biological calibration'),
        mode='production_frozen_gate_replay')
    try:
        raw = pd.read_csv(inputs['raw_scores'], sep='\t')
        raw = raw[raw.version.eq('corrected')].set_index(['condition', 'replicate'])
        assert len(raw) == 600 and raw.index.is_unique
        rows = []
        for condition in ('matched', 'state', 'two_clone'):
            for start in (0, 100):
                lock = json.loads(inputs[f'{condition}{start}_lock'].read_text())
                assert lock['count'] == 100 and not lock['truth_scored']
                for record in lock['cases']:
                    assert file_record(record['fit_manifest']['path'])['sha256'] == record['fit_manifest']['sha256']
                    before = raw.loc[(condition, record['replicate'])]
                    result = checked_gate(gate, Path(record['fit']), banks)
                    rows.append(dict(condition=condition, replicate=record['replicate'],
                        seed=record['seed'], **result, raw_groups=int(before.reported_groups),
                        raw_null_error=bool(before.null_error),
                        gated_null_error=bool(before.null_error and result['accepted']),
                        raw_recovery=bool(before.exact_chr12_recovery),
                        gated_recovery=bool(before.exact_chr12_recovery and result['accepted']),
                        raw_assignment_correct=int(before.assignment_correct),
                        supported_assignment_correct=int(before.assignment_correct) if result['accepted'] else 0,
                        assignment_n=600))
        assert len(rows) == 600
        frame = pd.DataFrame(rows)
        frame.to_csv(args.out/'all600_regression_cases.tsv', sep='\t', index=False)
        descriptions = []
        for condition in ('matched', 'state', 'two_clone'):
            selected = frame[frame.condition.eq(condition)]
            descriptions.append(dict(condition=condition, cases=len(selected),
                raw_false_splits=int(selected.raw_null_error.sum()) if condition != 'two_clone' else None,
                gated_false_splits=int(selected.gated_null_error.sum()) if condition != 'two_clone' else None,
                raw_recovery=int(selected.raw_recovery.sum()) if condition == 'two_clone' else None,
                gated_recovery=int(selected.gated_recovery.sum()) if condition == 'two_clone' else None,
                raw_assignment_correct=int(selected.raw_assignment_correct.sum()),
                supported_assignment_correct=int(selected.supported_assignment_correct.sum()),
                total_cells=600*len(selected)))
        development_rows = []
        for condition in ('matched', 'state', 'two_clone'):
            for start in (0, 5):
                report = json.loads(inputs[f'development_{condition}{start}'].read_text())
                assert report['count'] == 5
                for row in report['cases']:
                    manifest = row['after_manifest']
                    assert file_record(manifest['path'])['sha256'] == manifest['sha256']
                    decision = checked_gate(gate, Path(manifest['path']).parent, banks)
                    development_rows.append(dict(condition=condition, replicate=row['replicate'],
                        **decision, raw=row['after'],
                        gated_recovery=bool(row['after']['exact_chr12_recovery'] and decision['accepted']),
                        gated_null_error=bool(row['after']['null_error'] and decision['accepted'])))
        assert len(development_rows) == 30
        original_rows = []
        for row in json.loads(inputs['original_summary'].read_text())['cases']:
            seed = row['seed'];folder = original/f'seed_{seed}'
            result = checked_gate(gate, folder/'fit', banks)
            cells = pd.read_csv(folder/'assertion_view/cells.tsv.gz', sep='\t')
            differences = pd.read_csv(folder/'fit/differences.tsv', sep='\t')
            profiles = pd.read_csv(folder/'fit/clones.tsv', sep='\t')
            names = set(profiles.loc[profiles.specimen.eq('S1') & ~profiles.normal_like, 'clone'])
            before = scorer.criteria(cells, differences)
            assert before == row['criteria']
            if len(names) > 1 and not result['accepted']:
                selected = cells.specimen.eq('S1') & cells.malignant_call
                cells.loc[selected, 'clone'] = 'S1_unresolved'
                cells.loc[selected, 'status'] = 'unresolved_workflow_evidence'
                differences = differences.loc[~(differences.specimen.eq('S1')
                    & differences.clone_a.isin(names) & differences.clone_b.isin(names))]
            view = args.out/f'original_{seed}';view.mkdir()
            cells.to_csv(view/'cells.tsv.gz', sep='\t', index=False)
            differences.to_csv(view/'differences.tsv', sep='\t', index=False)
            shutil.copyfile(folder/'fit/qc.json', view/'qc.json')
            passed = True
            try:
                tests.test_thin_records_mode_and_recovers_clones({'thin':view})
            except AssertionError:
                passed = False
            original_rows.append(dict(seed=seed, **result,
                raw_passed=row['original_assertion_passed'], gated_passed=passed,
                raw_criteria=before, gated_criteria=scorer.criteria(cells, differences)))
        write_json(dict(descriptions=descriptions, development=development_rows,
            original_recovery=original_rows,
            original_raw_passed=sum(r['raw_passed'] for r in original_rows),
            original_gated_passed=sum(r['gated_passed'] for r in original_rows),
            thresholds=[max(b) for b in banks],
            historical_candidate_failure_reclassified=False, old_primary_result_overridden=False,
            calibration_source_matches_current_production_established=False,
            prospective_evaluation=False, scientific_criteria_changed=False,
            biological_specimens=0, promotion=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU saved outputs; unchanged frozen production rank gate')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
