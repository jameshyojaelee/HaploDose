"""Reproduce the qualified averaged profiles from both original count matrices."""
import argparse
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.pipeline import gene_loader
from haplodose.provenance import start_run, finish_run
from haplodose.stacking import stacked_count_profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    args = parser.parse_args()
    assert '\nexit\t0\n' in ('\n'+(args.validation/'job.tsv').read_text())
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    prior = root/'dev/gene_laplace_profile_20261110'
    qualified = root/'dev/partition_averaged_stack_20261110'
    assert json.loads((qualified/'summary_v1/report.json').read_text())['conditional_gate']
    assert json.loads((qualified/'summary_v1/run_manifest.json').read_text())['status'] == 'complete'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    identities = json.loads((base/'predictions_locked.json').read_text())['cases']
    locked = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')
    selected = [r for r in identities if r['case'] in locked]
    assert len(selected) == 2
    paths, decisions = {}, {}
    for start in (0, 46, 92, 138):
        directory = prior/f'regression_{start:03d}_v2'
        for r in json.loads((directory/'RNA_LOCKED.json').read_text())['cases']:
            if r['case'] in locked:
                paths[r['case']] = qualified/f"regression_{start:03d}_v1/predictions_case{r['prediction_index']}.tsv.gz"
                decisions[r['case']] = json.loads((qualified/f"regression_{start:03d}_v1/decision_case{r['prediction_index']}.json").read_text())
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        counts=root/'results/kinker_input_v2_20261101_1551/counts.h5',
        identities=base/'predictions_locked.json')
    for r in selected:
        rotation = r['rotation']
        run = root/f'results/kinker_field_20261102_0529/haplodose_r{rotation}'
        inputs[f'cells{rotation}'] = run/'cells.tsv.gz'
        inputs[f'genes{rotation}'] = run/'genes.tsv.gz'
        inputs['RNA_'+r['case']] = paths[r['case']]
        inputs['bins_'+r['case']] = base/f"bins_case{r['prediction_index']}.tsv"
    args.out.mkdir(exist_ok=False)
    rec = start_run(args.out, inputs, dict(cases=locked, equality='exact'),
        dict(scope='Source-count to installed API regression; reconstructed CPM, no count calibration; no DNA read'),
        mode='partition_average_count_parity')
    try:
        with h5py.File(inputs['counts']) as h:
            cell_ids, gene_ids = h['cell_ids'][:].astype(str), h['genes'][:].astype(str)
        assert pd.Index(cell_ids).is_unique and pd.Index(gene_ids).is_unique
        outputs = []
        for r in selected:
            meta = pd.read_csv(inputs[f"cells{r['rotation']}"], sep='\t', keep_default_na=False)
            used = pd.read_csv(inputs[f"genes{r['rotation']}"], sep='\t')
            assert meta.cell_id.is_unique and used.gene.is_unique
            query = meta.specimen.eq(r['case']) & meta.role.eq('query') & meta.design.eq('single_line')
            assert meta.loc[query, 'dataset'].nunique() == 1
            pool = meta.loc[query, 'dataset'].iloc[0]
            ref = meta.dataset.eq(pool) & meta.role.eq('reference') & meta.ref_split.eq('fit')
            local = meta.loc[query | ref].reset_index(drop=True)
            assert local.cell_id.tolist() == r['ordered_cells']
            rows = pd.Index(cell_ids).get_indexer(local.cell_id)
            cols = pd.Index(gene_ids).get_indexer(used.gene)
            assert (rows >= 0).all() and (cols >= 0).all()
            np.testing.assert_array_equal(rows, r['input_rows'])
            np.testing.assert_array_equal(cols, used.column)
            X, bins = gene_loader(inputs['counts'], used, rows)(np.arange(len(rows)))
            frozen = pd.read_csv(paths[r['case']], sep='\t', float_precision='round_trip')
            block = pd.read_csv(inputs['bins_'+r['case']], sep='\t').smoothing_block.to_numpy()
            np.testing.assert_array_equal(used.gene, frozen.gene)
            np.testing.assert_array_equal(bins, frozen.bin)
            profile, effects, result = stacked_count_profile(X, local.cell_id, used.gene, used[['gene','bin']],
                block, local.loc[local.role.eq('query'), 'cell_id'],
                local.loc[local.role.eq('reference'), 'cell_id'], counts_kind='reconstructed', gene_segments=True, average_segments=True)
            np.testing.assert_array_equal(effects, frozen.rna_log_ratio)
            np.testing.assert_array_equal(profile[bins], frozen.partition_average_stack)
            for key, value in decisions[r['case']].items():
                if key != 'units':
                    assert result[key] == value, (r['case'], key)
            assert result['units'] == 'natural_log_relative_rna_rate'
            write_json(result, args.out/f"decision_case{r['prediction_index']}.json")
            outputs.append(dict(case=r['case'], query_cells=r['query_cells'], reference_cells=r['reference_cells'],
                genes=len(effects), matrix_rows=rows.tolist(), cell_ids=local.cell_id.tolist(),
                gene_effects_exact=True, profile_exact=True, selection_exact=True))
        write_json(dict(cases=outputs, passed=True, counts_kind='reconstructed', DNA_read=False,
            frozen_outputs_changed=False, count_calibration_claimed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU identified source CSR counts and installed sparse RNA stack')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()

