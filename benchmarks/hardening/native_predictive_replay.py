"""Installed native profile selector must reproduce the retained RNA-only candidate."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.smoothing import predictive_rna_profile
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = Path('/gpfs/commons/home/jameslee/Haplodose')
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    prior = root/'dev/rna_predictive_smoothing_20261109/regression_v1'
    source = root/'dev/rna_predictive_smoothing_20261109/source_v3/benchmarks/hardening/rna_predictive_smoothing.py'
    records = json.loads((base/'predictions_locked.json').read_text())['cases']
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', reference_source=source,
                  identities=base/'predictions_locked.json', decisions=prior/'RNA_decisions_locked.json')
    for row in records:
        i = row['prediction_index']
        inputs[f'RNA{i}'] = base/f'predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
        inputs[f'prior{i}'] = prior/f'predictions_case{i}.tsv.gz'
    run = start_run(args.out, inputs, {}, dict(scope='Exact installed-selector regression; no prospective evaluation'),
                    mode='native_predictive_replay')
    try:
        spec = importlib.util.spec_from_file_location('frozen_rna_selector', source)
        reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        decisions = {r['case']: r for r in json.loads(inputs['decisions'].read_text())['cases']}
        comparisons = []
        for row in records:
            i = row['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'], sep='\t')
            old = pd.read_csv(inputs[f'prior{i}'], sep='\t', float_precision='round_trip')
            values = gene.rna_log_ratio.to_numpy()
            bin_ids, blocks = gene.bin.to_numpy(), bins.smoothing_block.to_numpy()
            result, decision = predictive_rna_profile(values, bin_ids, blocks)
            expected, expected_decision, _ = reference.select(values, bin_ids, blocks)
            np.testing.assert_array_equal(result, expected)
            np.testing.assert_array_equal(result[bin_ids], old.predictive_smoothed)
            for key, value in expected_decision.items():
                assert decision[key] == value, (row['case'], key)
            assert decision['selected_width'] == decisions[row['case']]['selected_width']
            comparisons.append(dict(case=row['case'], selected_width=decision['selected_width'],
                identical_to_retained_candidate=True, profile_sha256=hashlib.sha256(result.tobytes()).hexdigest()))
        assert len(comparisons) == 181
        write_json(dict(cases=comparisons, exact_identity=True, prospective_evaluation=False,
                        new_accuracy_evaluation=False, DNA_values_used=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU installed-wheel conditional RNA selector')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
