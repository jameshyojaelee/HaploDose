"""Exact installed/default reproduction of all qualified profile predictions."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.average_segments import _fit
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record
from haplodose.stacking import stacked_rna_profile


CHUNKS = ((0, 46), (46, 92), (92, 138), (138, 181))


def gate(args):
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/partition_averaged_stack_20261110'
    report = json.loads((prior/'summary_v1/report.json').read_text())
    assert json.loads((prior/'summary_v1/run_manifest.json').read_text())['status'] == 'complete'
    assert report['conditional_gate'] and report['family_size'] == 8
    assert report['old_profiles_exact'] == 181 and report['previous_scores_exact'] == 1448
    assert all(r['superior'] for r in report['contrasts']) and len(report['contrasts']) == 8
    assert '\nexit\t0\n' in ('\n'+(args.validation/'job.tsv').read_text())
    text = (args.validation/'pytest.txt').read_text()
    assert 'passed' in text and 'failed' not in text.splitlines()[-1]
    assert 'test result: ok' in (args.validation/'cargo_test.txt').read_text()
    return root, prior


def worker(args):
    assert (args.start, args.stop) in CHUNKS
    root, prior = gate(args)
    genes_prior = root/'dev/gene_laplace_profile_20261110'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        qualification=prior/'summary_v1/report.json', qualification_manifest=prior/'summary_v1/run_manifest.json',
        validated_wheel=args.validation/'installed_artifact.json', full_tests=args.validation/'pytest.txt',
        prior_decisions=genes_prior/f'regression_{args.start:03d}_v2/RNA_LOCKED.json')
    records = json.loads(inputs['prior_decisions'].read_text())['cases']
    assert len(records) == args.stop-args.start
    for row in records:
        i = row['prediction_index']
        folder = prior/f'regression_{args.start:03d}_v1'
        for name, filename in (('RNA', f'predictions_case{i}.tsv.gz'),
                               ('fit', f'fit_case{i}.npz'), ('decision', f'decision_case{i}.json')):
            inputs[f'{name}{i}'] = folder/filename
        inputs[f'bins{i}'] = root/f'dev/reference_weighted_pseudobulk_20261109/regression_v1/bins_case{i}.tsv'
    run = start_run(args.out, inputs, dict(start=args.start, stop=args.stop, equality='exact'),
        dict(scope='Installed descriptive RNA reproduction; no DNA parsed'), mode='partition_average_native_parity')
    try:
        output = []
        for frozen in records:
            i = frozen['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            blocks = pd.read_csv(inputs[f'bins{i}'], sep='\t').smoothing_block.to_numpy()
            values, bins = gene.rna_log_ratio.to_numpy(), gene.bin.to_numpy()
            default, previous = stacked_rna_profile(values, bins, blocks)
            np.testing.assert_array_equal(default[bins], gene.stacked_smoothed)
            assert previous == frozen['previous']
            old, old_decision = stacked_rna_profile(values, bins, blocks, gene_segments=True)
            np.testing.assert_array_equal(old[bins], gene.gene_laplace_stack)
            for key, value in frozen.items():
                if key not in ('case', 'prediction_index'):
                    assert old_decision[key] == value, (frozen['case'], key)
            candidate, standalone, decision, arrays = _fit(values, bins, blocks)
            expected = json.loads(inputs[f'decision{i}'].read_text())
            assert decision == expected, frozen['case']
            with np.load(inputs[f'fit{i}']) as fit:
                assert set(arrays) == set(fit.files)
                for key in arrays:
                    np.testing.assert_array_equal(arrays[key], fit[key])
            np.testing.assert_array_equal(candidate[bins], gene.partition_average_stack)
            np.testing.assert_array_equal(standalone[bins], gene.partition_average_alone)
            public, metadata = stacked_rna_profile(values, bins, blocks,
                gene_segments=True, average_segments=True)
            np.testing.assert_array_equal(public, candidate)
            for key, value in expected.items():
                if key != 'units':
                    assert metadata[key] == value, (frozen['case'], key)
            assert metadata['units'] == 'natural_log_relative_rna_rate'
            output.append(dict(case=frozen['case'], prediction_index=i,
                profile_exact=True, arrays_exact=True, decision_exact=True,
                previous_profile_exact=True, previous_decision_exact=True, default_exact=True))
            print(json.dumps(dict(case=frozen['case'], completed=len(output), exact=True)), flush=True)
        write_json(dict(cases=output, passed=True, DNA_read=False,
            frozen_outputs_changed=False), args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed wheel, local exact interval captures and RNA APIs')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def summary(args):
    root, _ = gate(args)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        identities=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1/predictions_locked.json',
        counts=args.task/'count_parity_v1/report.json', counts_manifest=args.task/'count_parity_v1/run_manifest.json')
    for start, _ in CHUNKS:
        inputs[f'part{start}'] = args.task/f'parity_{start:03d}_v1/report.json'
        inputs[f'manifest{start}'] = args.task/f'parity_{start:03d}_v1/run_manifest.json'
    run = start_run(args.out, inputs, dict(cases=181, equality='exact'),
        dict(scope='Exhaustive installed RNA reproduction and two source-count cases; no DNA parsed'),
        mode='partition_average_native_summary')
    try:
        cases = []
        for start, _ in CHUNKS:
            report = json.loads(inputs[f'part{start}'].read_text())
            assert report['passed'] and report['DNA_read'] is False
            assert json.loads(inputs[f'manifest{start}'].read_text())['status'] == 'complete'
            cases.extend(report['cases'])
        expected = json.loads(inputs['identities'].read_text())['cases']
        assert [r['case'] for r in cases] == [r['case'] for r in expected]
        assert [r['prediction_index'] for r in cases] == list(range(181))
        assert all(all(row[k] for k in ('profile_exact', 'arrays_exact', 'decision_exact',
            'previous_profile_exact', 'previous_decision_exact', 'default_exact')) for row in cases)
        counts = json.loads(inputs['counts'].read_text())
        assert counts['passed'] and len(counts['cases']) == 2 and counts['DNA_read'] is False
        assert json.loads(inputs['counts_manifest'].read_text())['status'] == 'complete'
        write_json(dict(cases=181, count_cases=2, passed=True, all_predictions_arrays_diagnostics_exact=True,
            previous_profiles_and_defaults_exact=True, DNA_read=False, frozen_results_changed=False,
            independent_evaluation=False, biological_calibration=False, default_promoted=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU exact installed/count-chain exhaustive coverage gate')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    parser.add_argument('--mode', choices=('worker', 'summary'), required=True)
    parser.add_argument('--start', type=int)
    parser.add_argument('--stop', type=int)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    (worker if args.mode == 'worker' else summary)(args)


if __name__ == '__main__':
    main()
