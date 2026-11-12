"""Replay the installed stacking API against all frozen RNA-only profiles."""
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
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    research = root/'dev/rna_stacked_smoothing_20261110'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    sources, cases = {}, []
    for directory in sorted(research.glob('regression_*_v3')):
        path = directory/'RNA_DECISIONS_LOCKED.json'
        decisions = json.loads(path.read_text())['cases']
        for decision in decisions:
            i = decision['prediction_index']
            sources[f'RNA{i}'] = directory/f'predictions_case{i}.tsv.gz'
            sources[f'bins{i}'] = base/f'bins_case{i}.tsv'
            cases.append(decision)
        sources[directory.name] = path
    assert len(cases) == 181 and len({r['case'] for r in cases}) == 181
    sources['script'] = Path(__file__)
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, sources, dict(exact_equality=True, cases=181),
        dict(scope='RNA-only installed API parity; no DNA or score-driven fitting'), mode='stacked_package_parity')
    try:
        rows = []
        for previous in sorted(cases, key=lambda r: r['prediction_index']):
            i = previous['prediction_index']
            gene = pd.read_csv(sources[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(sources[f'bins{i}'], sep='\t')
            output, decision = stacked_rna_profile(gene.rna_log_ratio.to_numpy(), gene.bin.to_numpy(),
                                                  bins.smoothing_block.to_numpy())
            np.testing.assert_array_equal(output[gene.bin.to_numpy()], gene.stacked_smoothed.to_numpy())
            # Ignore added public contract fields; require every frozen numerical
            # decision and diagnostic, including weights, to match exactly.
            for key, value in previous.items():
                if key not in ('case', 'prediction_index'):
                    assert decision[key] == value, (previous['case'], key)
            rows.append(dict(case=previous['case'], genes=len(gene), bins=len(output),
                             exact_profile=True, exact_diagnostics=True))
        pd.DataFrame(rows).to_csv(args.out/'parity.tsv', sep='\t', index=False)
        write_json(dict(cases=len(rows), all_profiles_exact=True, all_decisions_exact=True,
                        DNA_read=False, changed_frozen_outputs=False), args.out/'report.json')
        finish_run(args.out, run, backend='CPU installed wheel; sparse HiGHS')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
