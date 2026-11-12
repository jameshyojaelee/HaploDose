"""Installed optional gene segment stack versus every frozen RNA regression."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run
from haplodose.stacking import stacked_rna_profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int, default=181)
    args = parser.parse_args()
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    bins_root = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    inputs = dict(script=Path(__file__), design=args.task/'REGRESSION_DESIGN.md')
    records = []
    for start in (0, 46, 92, 138):
        folder = args.task/f'regression_{start:03d}_v2'
        inputs[f'decisions{start}'] = folder/'RNA_LOCKED.json'
        for record in json.loads(inputs[f'decisions{start}'].read_text())['cases']:
            i = record['prediction_index']
            inputs[f'RNA{i}'] = folder/f'predictions_case{i}.tsv.gz'
            inputs[f'bins{i}'] = bins_root/f'bins_case{i}.tsv'
            records.append(record)
    assert len(records) == 181
    assert 0 <= args.start < args.stop <= 181
    records = records[args.start:args.stop]
    args.out.mkdir(exist_ok=False)
    rec = start_run(args.out, inputs, dict(equality='exact', start=args.start, stop=args.stop),
        dict(scope='Installed descriptive RNA parity; no DNA read'), mode='gene_laplace_native_parity')
    try:
        output = []
        for frozen in records:
            i = frozen['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            blocks = pd.read_csv(inputs[f'bins{i}'], sep='\t').smoothing_block.to_numpy()
            values, bins = gene.rna_log_ratio.to_numpy(), gene.bin.to_numpy()
            old, previous = stacked_rna_profile(values, bins, blocks)
            np.testing.assert_array_equal(old[bins], gene.stacked_smoothed)
            assert previous == frozen['previous']
            new, decision = stacked_rna_profile(values, bins, blocks, gene_segments=True)
            np.testing.assert_array_equal(new[bins], gene.gene_laplace_stack)
            for key, value in frozen.items():
                if key not in ('case', 'prediction_index'):
                    assert decision[key] == value, (frozen['case'], key)
            output.append(dict(case=frozen['case'], prediction_index=i, profile_exact=True, decision_exact=True,
                               previous_profile_exact=True, previous_decision_exact=True))
            print(f"Exact parity: {i} {frozen['case']}", flush=True)
        write_json(dict(cases=output, passed=True, DNA_read=False,
            frozen_outputs_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU installed wheel; exact gene interval costs and sparse L1 stack')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
