"""RNA-only masked-gene selection of a robust regional smoothing width."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


LOCKED = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')


def running_median(values, blocks, width):
    out = np.empty(len(values))
    for block in pd.unique(blocks):
        ix = np.flatnonzero(blocks == block)
        out[ix] = pd.Series(values[ix]).rolling(width, center=True, min_periods=1).median().to_numpy()
    return out


def masked_prediction(values, gene_bins, blocks, fold, excluded, widths):
    train = fold != excluded
    medians = np.array([np.median(values[train & (gene_bins == b)])
                       if np.any(train & (gene_bins == b)) else np.nan for b in range(len(blocks))])
    return np.array([running_median(medians, blocks, w) for w in widths])


def select(values, gene_bins, blocks):
    values, gene_bins, blocks = np.asarray(values), np.asarray(gene_bins), np.asarray(blocks)
    assert values.shape == gene_bins.shape and np.isfinite(values).all()
    assert np.array_equal(np.unique(gene_bins), np.arange(len(blocks)))
    largest = max(int(np.sum(blocks == block)) for block in pd.unique(blocks))
    widths = [1]+[2**k+1 for k in range(1, largest.bit_length()) if 2**k+1 <= largest]
    folds = pd.Series(gene_bins).groupby(gene_bins, sort=False).cumcount().to_numpy()%10
    sizes = np.bincount(gene_bins, minlength=len(blocks))
    eligible = sizes[gene_bins] >= 2
    assert eligible.any()
    losses = np.zeros((10, len(widths)))
    held_out_n = np.zeros(10, int)
    for f in range(10):
        holdout = (folds == f) & eligible
        if not holdout.any():
            continue
        predicted = masked_prediction(values, gene_bins, blocks, folds, f, widths)
        target_predictions = predicted[:, gene_bins[holdout]]
        assert np.isfinite(target_predictions).all()
        losses[f] = np.sum(np.abs(target_predictions-values[holdout][None]), axis=1)
        held_out_n[f] = int(holdout.sum())
    assert held_out_n.sum() == eligible.sum()
    total = losses.sum(axis=0)
    eps = np.finfo(float).eps
    operations = 4*len(values)+32
    allowance = operations*eps/(1-operations*eps)*max(1., float(total.max()))
    tied = np.flatnonzero(total <= total.min()+2*allowance)
    index = int(tied[0])  # Smaller width when the numerical risk ordering is unresolved.
    full = np.array([np.median(values[gene_bins == b]) for b in range(len(blocks))])
    fitted = running_median(full, blocks, widths[index])
    return fitted, dict(widths=widths, selected_width=widths[index], loss='absolute RNA log-rate error',
        mean_prediction_error=(total/eligible.sum()).tolist(), fold_error_sums=losses.tolist(),
        held_out_genes=held_out_n.tolist(), ineligible_singleton_genes=int((~eligible).sum()),
        risk_sum_rounding_allowance=float(allowance), numerically_tied_widths=[widths[i] for i in tied]), folds


def checks():
    # Noiseless alternating regional effects have exactly known optimal width1.
    gbin = np.repeat(np.arange(7), 4)
    values = np.repeat([0., 2., 0., 2., 0., 2., 0.], 4)
    blocks = np.repeat('chr1p', 7)
    fitted, record, folds = select(values, gbin, blocks)
    assert record['selected_width'] == 1 and np.array_equal(fitted, values[::4])
    # Changing ONLY excluded observations must leave every training prediction identical.
    altered = values.copy()
    altered[folds == 0] += 100.
    before = masked_prediction(values, gbin, blocks, folds, 0, record['widths'])
    after = masked_prediction(altered, gbin, blocks, folds, 0, record['widths'])
    assert np.array_equal(before, after)
    # Independent Python-list reference for a masked width3 rolling median.
    import statistics
    train = [statistics.median([float(v) for v, b, f in zip(values, gbin, folds) if b == j and f != 0]) for j in range(7)]
    expected = [statistics.median(train[max(0, j-1):min(7, j+2)]) for j in range(7)]
    assert np.array_equal(before[record['widths'].index(3)], expected)
    constant, constant_record, _ = select(np.ones(len(values)), gbin, blocks)
    assert constant_record['selected_width'] == 1 and np.array_equal(constant, np.ones(7))
    # Arms cannot influence each other even with a large width.
    reset = running_median(np.r_[np.zeros(4), np.ones(4)*100], np.repeat(['p', 'q'], 4), 5)
    assert np.array_equal(reset, np.r_[np.zeros(4), np.ones(4)*100])
    return dict(noiseless_known_optimum=True, masked_value_perturbation_invariant=True,
        independent_scalar_medians_agree=True, constant_risk_tie_selects_smallest=True,
        genomic_resets_respected=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--all-lines', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    base = root/'dev/reference_weighted_pseudobulk_20261109/regression_v1'
    identities = json.loads((base/'predictions_locked.json').read_text())['cases']
    selected = identities if args.all_lines else [r for r in identities if r['case'] in LOCKED]
    assert len(selected) == (181 if args.all_lines else 2)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md',
        identities=base/'predictions_locked.json', truth=root/'results/kinker_input_v2_20261101_1551/depmap_cn.tsv.gz')
    if args.all_lines:
        inputs['regression_design'] = args.task/'REGRESSION_DESIGN.md'
    for r in selected:
        i = r['prediction_index']
        inputs[f'RNA{i}'] = base/f'predictions_case{i}.tsv.gz'
        inputs[f'bins{i}'] = base/f'bins_case{i}.tsv'
    rec = start_run(args.out, inputs, dict(all_lines=args.all_lines, folds=10,
        width_family='1 and 2^k+1 within largest block, k>=1', criterion='minimum masked-gene absolute RNA prediction error'),
        dict(scope='Conditional RNA model selection on preselected genes; no independent confirmation or count calibration',
             measurement='CPM-derived integer reconstruction', evaluation='historical regression only'),
        mode='rna_predictive_smoothing')
    try:
        numerical = checks()
        predictions, selection_records = [], []
        for r in selected:
            i = r['prediction_index']
            gene = pd.read_csv(inputs[f'RNA{i}'], sep='\t', float_precision='round_trip')
            bins = pd.read_csv(inputs[f'bins{i}'], sep='\t')
            assert gene.case.eq(r['case']).all() and gene.gene.is_unique
            assert np.array_equal(bins.bin, np.arange(len(bins)))
            fitted, decision, folds = select(gene.rna_log_ratio.to_numpy(), gene.bin.to_numpy(), bins.smoothing_block.to_numpy())
            gene['predictive_smoothed'] = fitted[gene.bin.to_numpy()]
            gene['masked_gene_fold'] = folds
            gene.to_csv(args.out/f'predictions_case{i}.tsv.gz', sep='\t', index=False)
            selection_records.append(dict(case=r['case'], prediction_index=i, **decision))
            predictions.append(gene)
        write_json(dict(cases=selection_records, files=[file_record(args.out/f'predictions_case{r["prediction_index"]}.tsv.gz') for r in selected]),
                   args.out/'RNA_decisions_locked.json')
        # DNA values only become available to the scorer after every RNA decision is saved.
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0)
        assert cn.index.is_unique and cn.columns.is_unique
        methods = ['pseudobulk', 'pseudobulk_sm5', 'pseudobulk_mdl', 'predictive_smoothed']
        scores = []
        for gene, r in zip(predictions, selected):
            assert len(r['query_DNA_line']) == 1
            weights = pd.Series(r['reference_DNA_weights'])
            common = gene.gene.isin(cn.index)
            g = gene.loc[common].copy()
            g['truth'] = cn.loc[g.gene, r['query_DNA_line'][0]].to_numpy()-cn.loc[g.gene, weights.index].to_numpy()@weights.to_numpy()
            assert np.isfinite(g[['truth']+methods].to_numpy()).all()
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1_000_000).astype(str)
            mb = g.groupby('mb')[['truth']+methods].mean()
            mb.to_csv(args.out/f'windows_case{r["prediction_index"]}.tsv', sep='\t')
            for method in methods:
                scores.append(dict(case=r['case'], method=method, common_genes=len(g), windows=len(mb),
                                   pearson=float(pearsonr(mb[method], mb.truth).statistic)))
        score_frame = pd.DataFrame(scores)
        score_frame.to_csv(args.out/'scores.tsv', sep='\t', index=False)
        matrix = score_frame.pivot(index='case', columns='method', values='pearson')
        delta = matrix.predictive_smoothed-matrix.pseudobulk_sm5
        write_json(dict(cases=len(selected), computational_checks=numerical,
            means=score_frame.groupby('method').pearson.mean().to_dict(),
            comparison_to_sm5=dict(wins=int((delta > 0).sum()), losses=int((delta < 0).sum()), ties=int((delta == 0).sum()),
                mean_paired_change=float(delta.mean()), median_paired_change=float(delta.median())),
            selected_widths={r['case']: r['selected_width'] for r in selection_records},
            independent_evaluation=False, biological_calibration=False, production_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU masked-gene RNA prediction and matched-support DNA regression')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
