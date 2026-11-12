"""Locked Kinker prediction probe: reference-line variability weighted pseudobulk."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from haplodose import pipeline
from haplodose.bins import arm_labels
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


CASES = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')


def weighted_median(values, weights):
    values, weights = np.asarray(values), np.asarray(weights)
    if len(values) == 0 or not np.isfinite(values).all() or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError('Weighted median requires finite observations and positive weights')
    order = np.argsort(values, kind='stable')
    v, w = values[order], weights[order]
    cumulative = np.cumsum(w)
    midpoint = cumulative[-1]/2
    i = np.searchsorted(cumulative, midpoint)
    return (v[i]+v[i+1])/2 if cumulative[i] == midpoint else v[i]


def smooth(values, blocks):
    out = np.empty(len(values))
    for block in pd.unique(blocks):
        ix = np.flatnonzero(blocks == block)
        out[ix] = pd.Series(values[ix]).rolling(5, center=True, min_periods=1).median().to_numpy()
    return out


def reference_residuals(shares, gbin):
    """Remove each reference line's local shift using the other bin genes."""
    residuals = shares.copy()
    for b in np.unique(gbin):
        ix = np.flatnonzero(gbin == b)
        if len(ix) == 1:
            continue
        for g in ix:
            residuals[:, g] -= np.median(shares[:, ix[ix != g]], axis=1)
    return residuals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--variance', choices=('total', 'within-bin'), default='total')
    parser.add_argument('--partition', action='store_true')
    parser.add_argument('--all-lines', action='store_true')
    parser.add_argument('--unweighted-only', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inp = root/'results/kinker_input_v2_20261101_1551'
    scorer_path = Path(__file__).parents[1]/'kinker/evaluate.py'
    spec = importlib.util.spec_from_file_location('weighted_pb_scorer', scorer_path)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    inputs = dict(script=Path(__file__), scorer=scorer_path, counts=inp/'counts.h5',
                  truth=inp/'depmap_cn.tsv.gz', design=args.task/'DESIGN.md',
                  identities=root/'dev/kinker_measurement_chain_20261109/chain_v1/identities.json')
    runs = [root/f'results/kinker_field_20261102_0529/haplodose_r{r}' for r in range(2)]
    for r, run in enumerate(runs):
        for name in ('cells.tsv.gz', 'genes.tsv.gz', 'bins.tsv'):
            inputs[f'r{r}_{name}'] = run/name
    if args.variance == 'within-bin':
        inputs['within_bin_design'] = args.task/'WITHIN_BIN_DESIGN.md'
    if args.partition:
        inputs['partition_design'] = args.task/'PARTITION_DESIGN.md'
        inputs['partition_source'] = Path(__file__).with_name('performance_candidates.py')
    if args.all_lines:
        inputs['regression_design'] = args.task/'REGRESSION_DESIGN.md'
    rec = start_run(args.out, inputs, dict(cases='all historical single-line cases' if args.all_lines else CASES,
        all_lines=args.all_lines, unweighted_only=args.unweighted_only,
        reference_variance_ddof=1, variance=args.variance, partition=args.partition,
        bin_estimator='inverse reference-line log-share variance weighted median', smoothing_bins=5),
        dict(scope='Historical regression replay, no prospective inference' if args.all_lines else
                   'Locked retrospective prediction probe; n=2 cell lines, no population inference',
             measurement='CPM-derived integer reconstruction; no raw-count uncertainty claims'),
        mode='reference_weighted_pseudobulk')
    try:
        # Independent algebraic fixtures, including even-sized median ties.
        for v in ([1., 4.], [1., 2., 4.], [4., 1., 2., 3.]):
            assert weighted_median(v, np.ones(len(v))) == np.median(v)
        assert weighted_median([0., 1., 2.], [1., 4., 1.]) == 1.
        assert weighted_median([0., 1., 2.], [7., 28., 7.]) == 1.
        fixture = np.array([[1., 2., 3., 5.], [3., 5., 7., 2.]])
        fixture_bins = np.array([0, 0, 1, 1])
        shift = np.array([[2., 2., 4., 4.], [8., 8., 16., 16.]])
        assert np.array_equal(reference_residuals(fixture, fixture_bins), reference_residuals(fixture+shift, fixture_bins))
        source_counts = None
        if args.all_lines:
            source_counts, source_genes, source_cells = scorer.load_counts(inputs['counts'])
        else:
            with h5py.File(inputs['counts']) as h:
                source_cells = h['cell_ids'][:].astype(str)
                source_genes = h['genes'][:].astype(str)
        scorer.unique_ids(source_cells, 'source cells')
        scorer.unique_ids(source_genes, 'source genes')
        locked = {r['case']: r for r in json.loads(inputs['identities'].read_text())}
        predictions, records = [], []
        tables = [(pd.read_csv(run/'cells.tsv.gz', sep='\t', keep_default_na=False),
                   pd.read_csv(run/'genes.tsv.gz', sep='\t'), pd.read_csv(run/'bins.tsv', sep='\t')) for run in runs]
        case_runs = [(r, case) for r, (meta, _, _) in enumerate(tables)
                     for case in sorted(meta.loc[meta.role.eq('query') & meta.design.eq('single_line'), 'specimen'].unique())]
        if not args.all_lines:
            case_runs = list(enumerate(CASES))
        for number, (rotation, case) in enumerate(case_runs):
            meta, used, bins = tables[rotation]
            scorer.unique_ids(meta.cell_id, 'metadata cells')
            scorer.unique_ids(used.gene, 'used genes')
            q = meta.specimen.eq(case) & meta.role.eq('query') & meta.design.eq('single_line')
            assert q.any() and meta.loc[q, 'dataset'].nunique() == 1
            dataset = meta.loc[q, 'dataset'].iloc[0]
            ref = meta.dataset.eq(dataset) & meta.role.eq('reference') & meta.ref_split.eq('fit')
            selected = np.flatnonzero(q | ref)
            local = meta.iloc[selected].reset_index(drop=True)
            rows = scorer.checked_indexer(source_cells, local.cell_id, 'selected cells')
            gi = scorer.checked_indexer(source_genes, used.gene, 'used genes')
            assert np.array_equal(gi, used.column)
            if case in locked:
                assert rows.tolist() == locked[case]['input_rows']
                assert local.cell_id.tolist() == locked[case]['ordered_cells']
                assert used.gene.tolist() == locked[case]['ordered_genes']
                assert used.bin.tolist() == locked[case]['gene_bins']
            if source_counts is None:
                X, gbin = pipeline.gene_loader(inputs['counts'], used, rows)(np.arange(len(rows)))
            else:
                X, gbin = source_counts[rows][:, gi], used.bin.to_numpy()
            assert X.shape == (len(rows), len(used)) and np.isfinite(X.data).all() and (X.data >= 0).all()
            is_query = local.role.eq('query').to_numpy()
            y = np.asarray(X[is_query].sum(axis=0)).ravel().astype(float)
            r = np.asarray(X[~is_query].sum(axis=0)).ravel().astype(float)
            assert y.sum() > 0 and r.sum() > 0
            effect = np.log((y+.5)/(r/r.sum()*y.sum()+.5))
            effect -= np.median(effect)
            lines = local.loc[~is_query, 'line'].to_numpy()
            reference_lines = sorted(set(lines))
            assert len(reference_lines) >= 2
            if args.unweighted_only:
                variability, precision = np.full(len(used), np.nan), np.ones(len(used))
            else:
                reference_counts = X[~is_query]
                pooled = np.array([np.asarray(reference_counts[lines == line].sum(axis=0)).ravel() for line in reference_lines], float)
                shares = np.log((pooled+.5)/(pooled.sum(axis=1, keepdims=True)+.5*len(used)))
                if args.variance == 'within-bin':
                    shares = reference_residuals(shares, gbin)
                variability = shares.var(axis=0, ddof=1)
                if args.variance == 'within-bin':
                    # Singleton bins have only one possible estimate, regardless of weight.
                    # A unit weight represents that identity, not an estimated precision.
                    singleton = np.bincount(gbin, minlength=len(bins))[gbin] == 1
                    variability[singleton] = 1.
                if not np.isfinite(variability).all() or np.any(variability <= 0):
                    raise ValueError('Reference-line variability not identified; no epsilon or fallback')
                precision = 1/variability
            assert np.array_equal(bins.bin, np.arange(len(bins)))
            naive = np.array([np.median(effect[gbin == b]) for b in bins.bin])
            candidate = np.array([weighted_median(effect[gbin == b], precision[gbin == b]) for b in bins.bin])
            labels = arm_labels(bins, genome_build='hg38')
            if 'genomic_block' in bins:
                labels = np.array([f'{a}|{b}' for a, b in zip(labels, bins.genomic_block)])
            naive_sm, candidate_sm = smooth(naive, labels), smooth(candidate, labels)
            gene = used[['gene', 'chromosome', 'start', 'end', 'bin']].copy()
            gene['case'] = case
            gene['rna_log_ratio'] = effect
            gene['reference_log_share_variance'] = variability
            gene['weight'] = precision
            for name, values in (('pseudobulk', naive), ('candidate', candidate),
                                 ('pseudobulk_sm5', naive_sm), ('candidate_sm5', candidate_sm)):
                if args.unweighted_only and name.startswith('candidate'):
                    continue
                gene[name] = values[gbin]
            if args.partition:
                partition_spec = importlib.util.spec_from_file_location('profile_partition', inputs['partition_source'])
                partition_module = importlib.util.module_from_spec(partition_spec)
                partition_spec.loader.exec_module(partition_module)
                reset = np.r_[True, labels[1:] != labels[:-1]]
                partition_records = []
                for name, track in (('pseudobulk', naive), ('candidate', candidate)):
                    if args.unweighted_only and name == 'candidate':
                        continue
                    # Robust Gaussian scale: MAD consistency constant and Var(diff)=2*Var(bin).
                    # This is a smoothing working model, not a coverage claim.
                    differences = np.diff(track)[~reset[1:]]
                    scale = np.median(np.abs(differences-np.median(differences)))/.6744897501960817/np.sqrt(2)
                    if not np.isfinite(scale) or scale <= 0:
                        raise ValueError('Positive smoothing scale not identified; no floor')
                    fitted, parts = partition_module.partition(track, np.full(len(track), scale**2), reset)
                    gene[name+'_mdl'] = fitted[gbin]
                    partition_records.append(dict(method=name, working_sd=float(scale), segments=parts))
                write_json(partition_records, args.out/f'partitions_case{number}.json')
            gene.to_csv(args.out/f'predictions_case{number}.tsv.gz', sep='\t', index=False)
            bins.assign(smoothing_block=labels).to_csv(args.out/f'bins_case{number}.tsv', sep='\t', index=False)
            rw = local.loc[~is_query].groupby('line').umi.sum()
            rw /= rw.sum()
            records.append(dict(case=case, rotation=rotation, prediction_index=number, query_cells=int(is_query.sum()),
                reference_cells=int((~is_query).sum()), reference_lines=reference_lines,
                reference_DNA_weights=rw.to_dict(), query_DNA_line=local.loc[is_query, 'truth_line'].unique().tolist(),
                input_rows=rows.tolist(), ordered_cells=local.cell_id.tolist(), genes=len(used)))
            predictions.append(gene)
        # Commit RNA-only predictions before loading any DNA values.
        write_json(dict(cases=records, files=[file_record(args.out/f'predictions_case{r}.tsv.gz') for r in range(len(records))]), args.out/'predictions_locked.json')
        cn = pd.read_csv(inputs['truth'], sep='\t', index_col=0)
        scorer.unique_ids(cn.index, 'DNA genes')
        scorer.unique_ids(cn.columns, 'DNA lines')
        scores = []
        for gene, record in zip(predictions, records):
            assert len(record['query_DNA_line']) == 1
            line = record['query_DNA_line'][0]
            rw = pd.Series(record['reference_DNA_weights'])
            assert line in cn and all(c in cn for c in rw.index)
            common = gene.gene.isin(cn.index)
            g = gene.loc[common].copy()
            truth = cn.loc[g.gene, line].to_numpy()-cn.loc[g.gene, rw.index].to_numpy()@rw.to_numpy()
            g['truth'] = truth
            g['mb'] = g.chromosome+':'+((g.start.to_numpy()+g.end.to_numpy())//2//1_000_000).astype(str)
            methods = ['pseudobulk', 'pseudobulk_sm5'] if args.unweighted_only else ['pseudobulk', 'candidate', 'pseudobulk_sm5', 'candidate_sm5']
            if args.partition:
                methods += ['pseudobulk_mdl'] if args.unweighted_only else ['pseudobulk_mdl', 'candidate_mdl']
            finite = np.isfinite(g[['truth']+methods].to_numpy()).all(axis=1)
            g = g.loc[finite]
            mb = g.groupby('mb')[['truth']+methods].mean()
            mb.to_csv(args.out/f'scored_windows_case{record["prediction_index"]}.tsv', sep='\t')
            for method in methods:
                scores.append(dict(case=record['case'], method=method, common_genes=len(g),
                    windows=len(mb), pearson=float(pearsonr(mb[method], mb.truth).statistic)))
        pd.DataFrame(scores).to_csv(args.out/'scores.tsv', sep='\t', index=False)
        write_json(dict(cases=len(records), biological_unit='cell line', scores=scores,
            population_superiority_established=False, production_changed=False,
            raw_count_calibration_claimed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU RNA-only reference variability weighting and matched-support DNA prediction')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
