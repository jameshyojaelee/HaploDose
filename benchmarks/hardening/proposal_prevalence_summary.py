"""Check completed saved inference after its final helper import failed."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose import clones
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    diagnostic = args.task/'diagnostic_v1'
    original = args.task.parent/'phase8_rep2_trace_20261110'
    failed = json.loads((diagnostic/'run_manifest.json').read_text())
    assert failed['status'] == 'failed'
    assert failed['error']['type'] == 'FileNotFoundError'
    assert failed['error']['message'].endswith("/refinement_checked_replay.py'")
    assert json.loads((diagnostic/'replay/run_manifest.json').read_text())['status'] == 'complete'
    assert json.loads((original/'summary_v1/report.json').read_text())['exact_replay_gate_passed']
    before = original/'trace_v1/replay/seed_202611098002'
    after = diagnostic/'replay/seed_202611098002'
    records = json.loads((diagnostic/'prevalence/records.json').read_text())
    assert len(records) == 2 and all('fit' in row for row in records)
    inputs = dict(script=Path(__file__), design=args.task/'READER_DESIGN.md',
        failed_manifest=diagnostic/'run_manifest.json',
        completed_inference=diagnostic/'replay/run_manifest.json',
        numerics=diagnostic/'numerics.json', records=diagnostic/'prevalence/records.json',
        original_trace=before/'pooled_prevalence/trace.json',
        original_results=original/'trace_v1/replay/results.json',
        completed_results=diagnostic/'replay/results.json',
        profiled_method=Path(__file__).with_name('proposal_prevalence_merge.py'))
    for arm in ('production', 'pooled_prevalence'):
        for side, folder in (('old', before/arm), ('new', after/arm)):
            for file in sorted((folder/'fit').iterdir()):
                inputs[f'{side}_{arm}_{file.name}'] = file
            for file in sorted((folder/'folds').glob('*.npz')):
                inputs[f'{side}_{arm}_fold_{file.name}'] = file
    for n, row in enumerate(records):
        inputs[f'prevalence_{n}'] = Path(row['arrays'])
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(inference_reruns=0, criteria_changed=False),
        dict(scope='Independent full saved-output gate; original failed checker retained'),
        mode='proposal_prevalence_summary')
    try:
        arrays = tables = folds = 0
        for arm in ('production', 'pooled_prevalence'):
            a, b = before/arm, after/arm
            eligible = lambda p: p.suffix == '.npy' or p.name.endswith(('.tsv', '.tsv.gz'))
            old_names = {p.name for p in (a/'fit').iterdir() if eligible(p)}
            new_names = {p.name for p in (b/'fit').iterdir() if eligible(p)}
            assert old_names == new_names
            for name in sorted(old_names):
                old_file, new_file = a/'fit'/name, b/'fit'/name
                if old_file.suffix == '.npy':
                    np.testing.assert_array_equal(np.load(old_file), np.load(new_file))
                    arrays += 1
                else:
                    pd.testing.assert_frame_equal(pd.read_csv(old_file, sep='\t'),
                        pd.read_csv(new_file, sep='\t'), check_exact=True)
                    tables += 1
            old_folds = {p.name for p in (a/'folds').glob('*.npz')}
            assert old_folds == {p.name for p in (b/'folds').glob('*.npz')}
            for name in sorted(old_folds):
                with np.load(a/'folds'/name) as x, np.load(b/'folds'/name) as y:
                    assert x.files == y.files
                    for key in x.files:
                        np.testing.assert_array_equal(x[key], y[key])
                folds += 1
        assert (arrays, tables, folds) == (4, 28, 1)
        old_trace = json.loads(inputs['original_trace'].read_text())
        old_merge = [r for r in old_trace['stages'] if r.get('stage') == 'merge']
        assert len(old_merge) == len(records)
        spec = importlib.util.spec_from_file_location('saved_proposal_profile', inputs['profiled_method'])
        method = importlib.util.module_from_spec(spec); spec.loader.exec_module(method)
        checks = json.loads(inputs['numerics'].read_text())
        assert checks['exact_rational_maximizer'] == '1/8'
        assert checks['exact_error'] <= checks['arithmetic_allowance']
        assert checks['independent_brent_error'] <= checks['arithmetic_allowance']
        for n, (old, row) in enumerate(zip(old_merge, records)):
            assert row['minimum'] == 30 and not row['inference_changed']
            np.testing.assert_array_equal(old['effective_cells'], row['original']['effective_cells'])
            assert old['testable'] == row['original']['testable']
            assert old['segments'] == len(row['original']['pre_veto_segments'])
            with np.load(inputs[f'prevalence_{n}']) as z:
                updated, fitted = method.fit_prevalence(z['proposal_loglik'], z['original_log_pi'])
                np.testing.assert_array_equal(updated, z['profiled_log_pi'])
                assert fitted == row['fit']
                np.testing.assert_array_equal(clones.posteriors(z['proposal_loglik'], z['original_log_pi']),
                                              z['original_R'])
                np.testing.assert_array_equal(clones.posteriors(z['proposal_loglik'], updated),
                                              z['profiled_R'])
        restored = [bool(r['proposal_profiled']['testable']
                        and r['proposal_profiled']['pre_veto_segments']) for r in records]
        write_json(dict(exact_output_parity=dict(arrays=arrays,tables=tables,folds=folds,exact=True),
            numerics=checks, records=records, support_and_pre_veto_restored=restored,
            sufficient_correction=bool(all(restored)),
            original_checker_failed=True, original_failure_preserved=True,
            inference_reruns=0, inference_changed=False, scientific_criteria_changed=False,
            biological_specimens=0, independent_evaluation=False, promotion=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU saved-output parity and scalar prevalence reference checks')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
