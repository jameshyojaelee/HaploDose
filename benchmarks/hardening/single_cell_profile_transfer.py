"""Bounded transfer of the delivered pooled RNA smoother to individual cells."""
import argparse
import json
import math
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose.io import write_json
from haplodose.pipeline import gene_loader
from haplodose.provenance import file_record, start_run, finish_run
from haplodose.smoothing import _running_median
from haplodose.stacking import stacked_count_profile, stacked_rna_profile


LOCKED = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')
METHODS = ('bin_median', 'naive_sm5', 'gene_laplace_stack', 'partition_average_stack')


def correlation(x, y):
    """Independent compensated Pearson calculation with a reduction budget."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if (x.shape != y.shape or len(x) < 2 or not np.isfinite(x).all()
            or not np.isfinite(y).all()):
        return dict(value=None, reason='nonfinite_or_insufficient_support')
    a = x-math.fsum(map(float, x))/len(x)
    b = y-math.fsum(map(float, y))/len(y)
    aa = math.fsum(float(v*v) for v in a)
    bb = math.fsum(float(v*v) for v in b)
    if aa == 0 or bb == 0:
        return dict(value=None, reason='constant_prediction_or_truth')
    ab = math.fsum(float(v*w) for v, w in zip(a, b))
    direct = ab/math.sqrt(aa*bb)
    scipy = float(pearsonr(x, y).statistic)
    # Centring cancellation is explicitly amplified by the uncentred magnitude
    # relative to each centred norm; this is arithmetic, not sampling error.
    operations = 64*len(x)+256
    epsilon = np.finfo(float).eps
    gamma = operations*epsilon/(1-operations*epsilon)
    condition = (1+np.linalg.norm(x)/math.sqrt(aa))*(1+np.linalg.norm(y)/math.sqrt(bb))
    budget = float(gamma*condition)
    if abs(direct-scipy) > budget:
        raise ArithmeticError('Independent Pearson reference exceeds reduction budget')
    return dict(value=scipy, independent_value=direct,
        independent_discrepancy=abs(direct-scipy), rounding_budget=budget)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    delivered = root/'dev/native_partition_average_20261110'
    assert json.loads((delivered/'summary_v1/report.json').read_text())['passed']
    assert json.loads((delivered/'summary_v1/run_manifest.json').read_text())['status'] == 'complete'
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    records = json.loads((base/'predictions_locked.json').read_text())['cases']
    chosen = [next(r for r in records if r['case'] == name) for name in LOCKED]
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        counts=root/'results/kinker_input_v2_20261101_1551/counts.h5',
        identities=base/'predictions_locked.json', delivery=delivered/'summary_v1/report.json',
        delivery_manifest=delivered/'summary_v1/run_manifest.json',
        truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz')
    for r in chosen:
        i, rotation = r['prediction_index'], r['rotation']
        chunk = max(v for v in (0, 46, 92, 138) if v <= i)
        inputs[f'cells{i}'] = root/f'results/kinker_field_20261102_0529/haplodose_r{rotation}/cells.tsv.gz'
        inputs[f'genes{i}'] = root/f'results/kinker_field_20261102_0529/haplodose_r{rotation}/genes.tsv.gz'
        inputs[f'RNA{i}'] = root/f'dev/partition_averaged_stack_20261110/regression_{chunk:03d}_v1/predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
    args.out.mkdir(exist_ok=False)
    run = start_run(args.out, inputs, dict(cases=LOCKED, methods=METHODS,
        ranks='0,(n-1)//2,n-1 in selected-feature depth then lexical cell ID'),
        dict(scope='Six dependent exposed cells; reconstructed measurements; bulk DNA compatibility only'),
        mode='single_cell_profile_transfer')
    try:
        with h5py.File(inputs['counts']) as h:
            cell_ids = h['cell_ids'][:].astype(str)
            gene_ids = h['genes'][:].astype(str)
        assert pd.Index(cell_ids).is_unique and pd.Index(gene_ids).is_unique
        outputs, selections, pooled_controls = [], [], []
        for r in chosen:
            i = r['prediction_index']
            meta = pd.read_csv(inputs[f'cells{i}'], sep='\t', keep_default_na=False)
            used = pd.read_csv(inputs[f'genes{i}'], sep='\t')
            query = meta.specimen.eq(r['case']) & meta.role.eq('query') & meta.design.eq('single_line')
            assert meta.loc[query, 'dataset'].nunique() == 1
            pool = meta.loc[query, 'dataset'].iloc[0]
            reference = meta.dataset.eq(pool) & meta.role.eq('reference') & meta.ref_split.eq('fit')
            local = meta.loc[query | reference].reset_index(drop=True)
            assert local.cell_id.tolist() == r['ordered_cells']
            rows = pd.Index(cell_ids).get_indexer(local.cell_id)
            columns = pd.Index(gene_ids).get_indexer(used.gene)
            assert (rows >= 0).all() and (columns >= 0).all()
            np.testing.assert_array_equal(rows, r['input_rows'])
            np.testing.assert_array_equal(columns, used.column)
            X, bins = gene_loader(inputs['counts'], used, rows)(np.arange(len(rows)))
            frozen = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            np.testing.assert_array_equal(used.gene, frozen.gene)
            np.testing.assert_array_equal(bins, frozen.bin)
            block = pd.read_csv(inputs[f'bins{i}'], sep='\t').smoothing_block.to_numpy()
            query_ids = local.loc[local.role.eq('query'), 'cell_id']
            ref_ids = local.loc[local.role.eq('reference'), 'cell_id']
            pooled, effects, decision = stacked_count_profile(X, local.cell_id, used.gene,
                used[['gene', 'bin']], block, query_ids, ref_ids,
                counts_kind='reconstructed', gene_segments=True, average_segments=True)
            np.testing.assert_array_equal(effects, frozen.rna_log_ratio)
            np.testing.assert_array_equal(pooled[bins], frozen.partition_average_stack)
            pooled_controls.append(dict(case=r['case'], effects_exact=True, profile_exact=True))
            qrows = pd.Index(local.cell_id).get_indexer(query_ids)
            totals = np.asarray(X[qrows].sum(axis=1)).ravel()
            ranked = sorted(range(len(qrows)), key=lambda k: (int(totals[k]), str(query_ids.iloc[k])))
            ranks = (0, (len(ranked)-1)//2, len(ranked)-1)
            assert len(set(ranks)) == 3
            for label, rank in zip(('low', 'median', 'high'), ranks):
                index = qrows[ranked[rank]]
                identity = str(local.cell_id.iloc[index])
                profile, effects, decision = stacked_count_profile(X, local.cell_id, used.gene,
                    used[['gene', 'bin']], block, [identity], ref_ids,
                    counts_kind='reconstructed', gene_segments=True, average_segments=True)
                old, old_decision = stacked_rna_profile(effects, bins, block, gene_segments=True)
                naive = np.array([np.median(effects[bins == b]) for b in range(len(block))])
                name = f'case{i}_{label}'
                frame = frozen[['gene', 'chromosome', 'start', 'end', 'bin']].copy()
                frame['rna_log_ratio'] = effects
                frame['bin_median'] = naive[bins]
                frame['naive_sm5'] = _running_median(naive, block, 5)[bins]
                frame['gene_laplace_stack'] = old[bins]
                frame['partition_average_stack'] = profile[bins]
                frame['pooled_reference_only'] = pooled[bins]
                frame.to_csv(args.out/f'{name}.tsv.gz', sep='\t', index=False)
                write_json(dict(new=decision, previous=old_decision), args.out/f'{name}_decision.json')
                selections.append(dict(case=r['case'], label=label, rank=rank,
                    query_cells=len(qrows), cell_id=identity, input_row=int(rows[index]),
                    selected_feature_total=int(X[index].sum()),
                    nonzero_features=int(X[index].count_nonzero()),
                    ordered_references=ref_ids.tolist(), ordered_genes=used.gene.tolist(),
                    prediction=file_record(args.out/f'{name}.tsv.gz'),
                    diagnostics=file_record(args.out/f'{name}_decision.json')))
                outputs.append((r, label, frame, decision))
                print(json.dumps(dict(case=r['case'], cell=identity, depth=label,
                    selected=decision['selected'], completed=len(outputs))), flush=True)
        write_json(dict(selections=selections, pooled_controls=pooled_controls,
            DNA_parsed=False), args.out/'RNA_LOCKED.json')
        # Only this scoring step sees corresponding bulk DNA measurements.
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0, float_precision='round_trip')
        assert cn.index.is_unique and cn.columns.is_unique
        scores = []
        for r, label, frame, decision in outputs:
            weights = pd.Series(r['reference_DNA_weights'])
            gene = frame[frame.gene.isin(cn.index)].copy()
            gene['truth'] = (cn.loc[gene.gene, r['query_DNA_line'][0]].to_numpy()
                            -cn.loc[gene.gene, weights.index].to_numpy()@weights.to_numpy())
            gene['mb'] = gene.chromosome+':'+((gene.start+gene.end)//2//1000000).astype(str)
            projected = gene.groupby('mb')[['truth', *METHODS, 'pooled_reference_only']].mean()
            for method in (*METHODS, 'pooled_reference_only'):
                result = correlation(projected[method], projected.truth)
                value = projected[method].to_numpy()
                norm = float(np.linalg.norm(value-value.mean()))
                export = (decision.get('candidate_profile_rounding', 0.)
                          if method == 'partition_average_stack' else 0.)
                perturbation = math.sqrt(len(projected))*export
                bound = 2*perturbation/(norm-perturbation) if norm > perturbation else None
                scores.append(dict(case=r['case'], depth=label, method=method,
                    windows=len(projected), common_genes=len(gene), **result,
                    conditional_export_bound=bound))
        pd.DataFrame(scores).to_csv(args.out/'scores.tsv', sep='\t', index=False)
        comparisons = []
        for name in LOCKED:
            subset = [v for v in scores if v['case'] == name]
            new = [v for v in subset if v['method'] == 'partition_average_stack']
            assert len(new) == 3 and {v['depth'] for v in new} == {'low', 'median', 'high'}
            for comparator in ('gene_laplace_stack', 'naive_sm5'):
                old = [v for v in subset if v['method'] == comparator]
                finite = all(v['value'] is not None and v['conditional_export_bound'] is not None
                             for v in new+old)
                delta = (math.fsum(v['value'] for v in new)/3
                         -math.fsum(v['value'] for v in old)/3) if finite else None
                budget = math.fsum(v['rounding_budget']+v['conditional_export_bound']
                                  for v in new+old)/3 if finite else None
                comparisons.append(dict(case=name, comparator=comparator, mean_difference=delta,
                    conditional_numerical_budget=budget, survived=bool(finite and delta > budget)))
        write_json(dict(cells=len(outputs), lines=2, pooled_controls=pooled_controls,
            comparisons=comparisons, small_panel_transfer_survived=all(v['survived'] for v in comparisons),
            interpretation='Bulk-line DNA compatibility; not same-cell truth or independent validation',
            biological_replication_for_error_control=0, independent_evaluation=False,
            calibrated_uncertainty=False, default_promoted=False, frozen_results_changed=False),
            args.out/'report.json')
        for record in run['input']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        finish_run(args.out, run, backend='CPU installed wheel; six identified reconstructed single-cell profiles')
    except BaseException as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
