"""Locked development-stage diagnostics from frozen artifacts; no model fitting.

Kinker: 42MGBA and 2313287_STOMACH. TCGA: TCGA-05-4389-01, whose
normalized-RSEM pseudo-counts are a processing diagnostic only. Optionally add
the completed locked U-CAN run. Reference membership, clone weights and HD
segment boundaries are taken from the recorded runs, never reselected.

Naive-on-HD-boundaries is a conditional replay control: boundaries were chosen
using HD's data. It cannot establish independent HMM or causal component benefit.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys

import h5py
import numpy as np
import pandas as pd

CODE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CODE_ROOT / 'benchmarks'))
from hardening.ucan import file_record, metrics, smooth

spec = importlib.util.spec_from_file_location('artifact_subarm', CODE_ROOT / 'benchmarks/tcga_bulk/subarm.py')
subarm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(subarm)

KINKER_CASES = ((0, '2313287_STOMACH'), (1, '42MGBA_CENTRAL_NERVOUS_SYSTEM'))
LUAD_CASE = 'TCGA-05-4389-01'


def segment_projection(values, segments):
    """Mean over available bins per frozen boundary; require one cover per bin.

    RNA methods being contrasted must have matched finite support before replay.
    Segments without any supported bins remain undefined.
    """
    values = np.asarray(values, float)
    out, coverage = np.full_like(values, np.nan), np.zeros(len(values), int)
    for r in segments.itertuples():
        a, b = int(r.first_bin), int(r.last_bin)
        if not 0 <= a <= b < len(values):
            raise ValueError('frozen segment lies outside its profile')
        coverage[a:b + 1] += 1
        v = values[a:b + 1]
        out[a:b + 1] = float(np.mean(v[np.isfinite(v)])) if np.isfinite(v).any() else np.nan
    if not np.all(coverage == 1):
        raise ValueError('recorded segments must cover every bin once without overlap')
    return out


def read_specimen(path, specimen):
    selected = []
    for chunk in pd.read_csv(path, sep='\t', chunksize=100000):
        selected.append(chunk[chunk.specimen.eq(specimen)])
    result = pd.concat(selected, ignore_index=True)
    if result.empty:
        raise ValueError(f'locked specimen {specimen} absent from {path}')
    return result


def sum_h5_rows(handle, indptr, rows, n_genes):
    """Read only the locked query/reference rows; do not materialize the full CSR."""
    result = np.zeros(n_genes, np.float64)
    for row in np.asarray(rows, int):
        lo, hi = int(indptr[row]), int(indptr[row + 1])
        np.add.at(result, handle['X/indices'][lo:hi], handle['X/data'][lo:hi])
    return result


def score_case(out, case, bins, truth, methods, width, note, sources, metadata):
    """Same bin/arm support for all listed methods, with no biological-unit p-values."""
    bins = bins.reset_index(drop=True).copy()
    truth = np.asarray(truth, float)
    if len(bins) != len(truth) or any(np.shape(v) != truth.shape for v in methods.values()):
        raise ValueError('diagnostic methods/truth must share one genomic grid')
    common = np.isfinite(truth) & np.logical_and.reduce([np.isfinite(v) for v in methods.values()])
    table = bins.copy()
    table['case'], table['truth'], table['common_callable'] = case, truth, common
    records = []
    for method, values in methods.items():
        table[method] = values
        records.append(dict(case=case, method=method, resolution='common_bin', **metrics(values[common], truth[common])))
        if 'focal_loss' in bins and 'neutral' in bins:
            loss = common & bins.focal_loss.to_numpy(bool)
            neutral = common & bins.neutral.to_numpy(bool)
            records[-1].update(n_focal_loss=int(loss.sum()), n_neutral=int(neutral.sum()),
                               focal_loss_signal_median=float(np.median(values[loss])) if loss.any() else np.nan,
                               neutral_signal_median=float(np.median(values[neutral])) if neutral.any() else np.nan)
        frame = pd.DataFrame(dict(arm=bins.arm.to_numpy()[common], value=values[common], truth=truth[common]))
        for aggregation in ('median', 'mean'):
            arms = frame.groupby('arm')[['value', 'truth']].agg(aggregation)
            records.append(dict(case=case, method=method, resolution='common_arm_' + aggregation,
                                **metrics(arms.value.to_numpy(), arms.truth.to_numpy())))
    table.to_csv(out / f'{case}_profiles.tsv.gz', sep='\t', index=False)
    result = dict(case=case, evidence_role='development_artifact_diagnostic', independent_evaluation=False,
                  biological_n=1, interpretation=note, n_common_bins=int(common.sum()), n_bins=len(bins),
                  fixed_smoothing_width_bins=width, metadata=metadata,
                  source_files=[file_record(p) for p in sources], scores=records)
    (out / f'{case}_summary.json').write_text(json.dumps(result, indent=2))
    return result


def kinker_cases(data_root, out, width):
    inp = data_root / 'results/kinker_input_v2_20261101_1551'
    run_root = data_root / 'results/kinker_field_20261102_0529'
    cn_path = inp / 'depmap_cn.tsv.gz'
    cn = pd.read_csv(cn_path, sep='\t', index_col=0)
    if cn.index.has_duplicates:
        raise ValueError('ambiguous gene symbols in recorded DepMap truth')
    cen = pd.read_csv(subarm.CEN38, sep='\t', comment='#').set_index('chromosome').p_arm_end
    results = []
    with h5py.File(inp / 'counts.h5', 'r') as handle:
        genes, cell_ids = handle['genes'][:].astype(str), handle['cell_ids'][:].astype(str)
        subarm.unique_ids(genes, 'frozen count gene IDs')
        cell_index = subarm.unique_ids(cell_ids, 'frozen count cell IDs')
        indptr = handle['X/indptr'][:]
        for rotation, case in KINKER_CASES:
            run = run_root / f'haplodose_r{rotation}'
            cells = pd.read_csv(run / 'cells.tsv.gz', sep='\t')
            query = cells[cells.specimen.eq(case) & cells.design.eq('single_line')]
            if query.empty or query.dataset.nunique() != 1 or query.truth_line.nunique() != 1:
                raise ValueError('locked Kinker identity is ambiguous')
            reference = cells[cells.dataset.eq(query.dataset.iloc[0]) & cells.role.eq('reference') & cells.ref_split.eq('fit')]
            qr, rr = cell_index.get_indexer(query.cell_id), cell_index.get_indexer(reference.cell_id)
            if not len(rr) or np.any(qr < 0) or np.any(rr < 0):
                raise ValueError('locked Kinker cells missing from original counts')
            used = pd.read_csv(run / 'genes.tsv.gz', sep='\t')
            used = used[used.gene.isin(genes)].copy()
            gidx = pd.Index(genes).get_indexer(used.gene)
            y = sum_h5_rows(handle, indptr, qr, len(genes))[gidx]
            r = sum_h5_rows(handle, indptr, rr, len(genes))[gidx]
            expected = r / r.sum() * y.sum()
            pseudobulk = np.log((y + .5) / (expected + .5))
            pseudobulk -= np.median(pseudobulk)
            weights = reference.groupby('line').umi.sum()
            truth_cn = cn.reindex(used.gene)
            truth_g = truth_cn[query.truth_line.iloc[0]].to_numpy() - (
                truth_cn[weights.index].to_numpy() @ (weights / weights.sum()).to_numpy())
            bins = pd.read_csv(run / 'bins.tsv', sep='\t')
            bins['arm'] = bins.chromosome + np.where(bins.start < bins.chromosome.map(cen), 'p', 'q')
            gene_bin = used.bin.to_numpy(int)
            pb = subarm.window_median(pseudobulk[:, None], gene_bin, len(bins))[:, 0]
            truth = subarm.window_median(truth_g[:, None], gene_bin, len(bins))[:, 0]
            profiles = read_specimen(run / 'profiles.tsv.gz', case)
            segment = read_specimen(run / 'segments.tsv', case)
            clone_weights = query.clone.value_counts()
            clone_weights = clone_weights[clone_weights.index.isin(profiles.clone)]
            if clone_weights.empty:
                raise ValueError('locked query has no recorded clone profile')
            clone_weights = clone_weights / clone_weights.sum()
            methods = dict(pseudobulk=pb)
            for name, column in (('hd_gene_ratio', 'log_ratio_gene'), ('hd_bin_ratio', 'log_ratio'), ('hd_posterior', 'log_dose')):
                p = profiles.pivot(index='bin', columns='clone', values=column).reindex(range(len(bins)))
                methods[name] = (p[clone_weights.index] * clone_weights.to_numpy()).sum(axis=1, min_count=len(clone_weights)).to_numpy()
            exported, pb_boundary, hd_boundary = (np.zeros(len(bins)) for _ in range(3))
            for clone, weight in clone_weights.items():
                intervals = segment[segment.clone.eq(clone)]
                p = profiles[profiles.clone.eq(clone)].set_index('bin').reindex(range(len(bins)))
                pb_boundary += weight * segment_projection(pb, intervals)
                hd_boundary += weight * segment_projection(p.log_ratio_gene.to_numpy(), intervals)
                value = np.full(len(bins), np.nan)
                for row in intervals.itertuples():
                    value[row.first_bin:row.last_bin + 1] = row.log_dose
                exported += weight * value
            methods.update(hd_segment_export=exported,
                           pseudobulk_on_hd_boundaries=pb_boundary, hd_gene_ratio_on_hd_boundaries=hd_boundary)
            for name, value in list(methods.items()):
                if not name.endswith('boundaries'):
                    methods[name + '_smoothed'] = smooth(value, bins.arm.to_numpy(), width)
            results.append(score_case(out, case, bins, truth, methods, width,
                'Frozen development run. Boundary replay is conditional on HD-selected clone/segment structure; '
                'no independent segmentation or causal HMM benefit is identified. DepMap target retains the '
                'historical line-minus-UMI-weighted-reference-line definition in its native transformed units.',
                [inp / 'counts.h5', cn_path, run / 'cells.tsv.gz', run / 'profiles.tsv.gz', run / 'segments.tsv',
                 run / 'genes.tsv.gz', run / 'bins.tsv'],
                dict(rotation=rotation, n_query_cells=len(query), n_fit_reference_cells=len(reference),
                     reference_dataset=query.dataset.iloc[0], clone_cell_weights=clone_weights.to_dict(),
                     no_model_fitting=True, source_coordinate_system='legacy frozen genomic positions')))
    return results


def luad_case(data_root, out, width, gene_width):
    inp = data_root / 'results/tcga_bulk_20261101_1405/input/LUAD'
    hd = inp / 'haplodose'
    subarm.GTF37 = data_root / 'data/external/ensembl_grch37_r75/Homo_sapiens.GRCh37.75.gtf.gz'
    g37 = subarm.hg19_genes()
    pend, size = subarm.hg19_arms()
    tab, segments = subarm.absolute()
    d = subarm.type_windows('LUAD_diagnostic', inp, None, g37, tab, segments, pend, size)
    found = np.flatnonzero(d['tum'] == LUAD_CASE)
    if len(found) != 1:
        raise ValueError('locked LUAD specimen absent or ambiguous')
    j = int(found[0])
    raw = {name: values[:, j] for name, values in d['gene_values'].items()}
    methods = {name: values[:, j] for name, values in d['W'].items()}
    gb = pd.read_csv(hd / 'genes.tsv.gz', sep='\t').set_index('gene').bin
    gb.index = subarm.canonical_symbols(gb.index)
    mapping = gb.reindex(d['genes']).to_numpy(float)
    n_bins = int(np.load(hd / 'log_ratio.npy', mmap_mode='r').shape[1])
    intervals = pd.read_csv(hd / 'segments.tsv', sep='\t')
    intervals = intervals[intervals['sample'].eq(LUAD_CASE)]
    for name in ('naive', 'hd_lr'):
        # Replay the exact frozen HD segmentation on both methods on its native
        # bin grid, then map to the same gene-anchored evaluation windows.
        base = np.full(n_bins, np.nan)
        for k in range(n_bins):
            v = raw[name][mapping == k]
            if np.isfinite(v).any():
                base[k] = np.nanmedian(v)
        projected = segment_projection(base, intervals)
        values = np.full(len(mapping), np.nan)
        has = np.isfinite(mapping)
        values[has] = projected[mapping[has].astype(int)]
        methods[name + '_on_hd_boundaries'] = subarm.window_median(values[:, None], d['window_id'], len(d['warm']))[:, 0]
        filtered = smooth(raw[name], d['gene_positions'].chromosome.to_numpy(), gene_width)
        methods[name + '_matched_gene_smooth'] = subarm.window_median(filtered[:, None], d['window_id'], len(d['warm']))[:, 0]
    for name, values in list(methods.items()):
        if name in ('naive', 'hd_lr', 'hd_dose'):
            methods[name + '_bin_smoothed'] = smooth(values, d['warm'], width)
    bins = pd.DataFrame(dict(bin=np.arange(len(d['warm'])), arm=d['warm'],
                             focal_loss=d['FL'][:, j], neutral=d['NEU'][:, j]))
    return score_case(out, LUAD_CASE, bins, d['D'][:, j], methods, width,
        'Normalized, batch-adjusted RSEM pseudo-count processing diagnostic only. No count likelihood is fit. '
        'Per-bin posterior, expected counts and programme shifts were not saved, so the HMM/posterior/export '
        'stages cannot be separated. Naive boundary replay uses HD-selected boundaries and is conditional.',
        [inp / 'counts.tsv', inp / 'samples.tsv', inp / 'positions.tsv', hd / 'log_ratio.npy', hd / 'segments.tsv',
         hd / 'genes.tsv.gz', subarm.GTF37, subarm.CN / 'TCGA_mastercalls.abs_segtabs.fixed.txt'],
        dict(input_scale='normalized_RSEM_pseudocounts', no_model_fitting=True,
             matched_gene_smoothing_width=gene_width, threshold_estimation=False,
             historical_focal_recall=dict(hd_dose=.05, hd_lr=.50, naive=.25, naive_smooth=.30),
             historical_recall_source=str(data_root / 'results/tcga_bulk_20261101_1405/subarm_20261101_1431/score_fprmatched_random_20261101_1442/per_patient.tsv.gz')))


def ucan_boundaries(run, out, width):
    trace = pd.read_csv(run / 'bin_stage_trace.tsv.gz', sep='\t')
    intervals = pd.read_csv(run / 'haplodose/segments.tsv', sep='\t')
    results = []
    for sample, frame in trace.groupby('sample', sort=False):
        frame = frame.sort_values('bin').reset_index(drop=True)
        native = intervals[intervals['sample'].eq(sample)]
        names = ['gene_ratio', 'naive', 'pooled_bin_ratio', 'hd_ratio', 'hd_posterior', 'hd_segment_export']
        methods = {name: frame[name].to_numpy() for name in names}
        for name in ('gene_ratio', 'naive', 'pooled_bin_ratio', 'hd_ratio'):
            methods[name + '_on_hd_boundaries'] = segment_projection(methods[name], native)
        blocks = frame.genomic_block.to_numpy() if 'genomic_block' in frame else frame.chromosome.to_numpy()
        for name in names:
            methods[name + '_smoothed'] = smooth(methods[name], blocks, width)
        results.append(score_case(out, sample, frame[['bin', 'arm']], frame.dna_ln_tcn_over_ploidy.to_numpy(),
                                 methods, width, 'Two locked raw-count development examples; boundary replay is '
                                 'conditional on HD-selected intervals. It is not independent naive segmentation.',
                                 [run / 'bin_stage_trace.tsv.gz', run / 'haplodose/segments.tsv', run / 'provenance.json'],
                                 dict(no_model_fitting=True, prior_current_wheel_run=str(run))))
    return results


def summarize_profile(run, out):
    """Describe saved external-threshold transfer; no recalibration or new testing."""
    endpoints = ['fpr_loss', 'recall_loss', 'precision_loss',
                 'classified_truth_coverage', 'recall_loss_failure_aware',
                 'rho_window', 'pearson_window']
    raw = pd.read_csv(run / 'raw_patient_scores.tsv', sep='\t')
    xena = pd.read_csv(run / 'xena_processing_diagnostic.tsv', sep='\t')
    manifest = raw[['sample', 'split']].drop_duplicates()
    xena = xena.merge(manifest, on='sample', how='left', validate='many_to_one')
    records = []
    for processing, frame in [('raw_STAR', raw), ('Xena_processing_diagnostic', xena)]:
        for (split, method), group in frame.groupby(['split', 'method']):
            for endpoint in endpoints:
                values = group[endpoint]
                records.append(dict(processing=processing, split=split, method=method,
                                    endpoint=endpoint, patient_weighted_mean=values.mean(),
                                    n_patients=len(group), n_defined=int(values.notna().sum()),
                                    n_undefined=int(values.isna().sum())))
    pd.DataFrame(records).to_csv(out / 'threshold_transfer_summary.tsv', sep='\t', index=False)
    pairs = pd.read_csv(run / 'matched_processing_baseline_pairs.tsv', sep='\t')
    records = []
    for (split, method), group in pairs.groupby(['split', 'method']):
        for endpoint in endpoints:
            valid = group[[endpoint + '_raw', endpoint + '_xena']].dropna()
            difference = valid[endpoint + '_raw'] - valid[endpoint + '_xena']
            records.append(dict(split=split, method=method, endpoint=endpoint,
                                mean_raw_minus_xena=difference.mean(), n_paired=len(valid),
                                n_patients=len(group), n_undefined_pairs=len(group) - len(valid)))
    pd.DataFrame(records).to_csv(out / 'processing_paired_summary.tsv', sep='\t', index=False)
    (out / 'threshold_transfer_provenance.json').write_text(json.dumps(dict(
        source_files=[file_record(run / name) for name in
                      ['raw_patient_scores.tsv', 'xena_processing_diagnostic.tsv',
                       'matched_processing_baseline_pairs.tsv', 'summary.json']],
        replication_unit='patient', independent_evaluation=False, new_testing=False,
        interpretation='Descriptive patient means at previously frozen neutral-only calibration thresholds. '
        'Reserved development patients assess transfer; calibration rows are resubstitution. '
        'Precision excludes unclassified truth and is undefined without classified calls. '
        'Raw-Xena paired baselines share normal donors and features; frozen Xena HD has a different reference panel.'), indent=2))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--smooth-bins', type=int, required=True)
    parser.add_argument('--smooth-genes', type=int, required=True)
    parser.add_argument('--ucan-run', type=Path)
    parser.add_argument('--profile-run', type=Path)
    args = parser.parse_args(argv)
    if any(x < 1 or x % 2 == 0 for x in (args.smooth_bins, args.smooth_genes)):
        raise ValueError('explicit smoothing widths must be positive odd integers')
    args.out.mkdir(parents=True, exist_ok=False)
    subarm.CEN38 = args.data_root / 'python/haplodose/data/centromeres_hg38.tsv'
    results = kinker_cases(args.data_root, args.out, args.smooth_bins)
    print('locked Kinker artifacts complete', flush=True)
    results.append(luad_case(args.data_root, args.out, args.smooth_bins, args.smooth_genes))
    print('locked LUAD processing diagnostic complete', flush=True)
    if args.ucan_run:
        results.extend(ucan_boundaries(args.ucan_run, args.out, args.smooth_bins))
    if args.profile_run:
        summarize_profile(args.profile_run, args.out)
    pd.DataFrame([row for result in results for row in result['scores']]).to_csv(args.out / 'scores.tsv', sep='\t', index=False)
    (args.out / 'summary.json').write_text(json.dumps(dict(
        executed_script=file_record(__file__), no_model_fitting=True, independent_evaluation=False,
        evidence_role='locked_development_artifact_diagnostics', no_new_thresholds=True,
        operator_control='HD-selected fixed-boundary projection; conditional, not independent naive HMM',
        cases=results), indent=2))


if __name__ == '__main__':
    main()
