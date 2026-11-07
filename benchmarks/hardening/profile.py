"""Small, versioned development falsifiers for bulk/profile inference.

Only existing raw STAR COAD/READ is refitted. Its paired normalized-Xena output
is a frozen processing diagnostic; it is never fed to the count likelihood.
Former evaluation cohorts are neither opened nor tuned. All writes require a
new output root. Run on a compute node, e.g.:

  python benchmarks/hardening/profile.py --out NEWROOT --reference-refits 8

The fixed method contrasts test segmentation, matched smoothing, common support
and transfer of independently calibrated loss thresholds. This is development
evidence. Few reference bootstrap refits cannot establish nominal CI coverage.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.stats import pearsonr

ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get('HAPLODOSE_DATA_ROOT', str(ROOT)))
sys.path.insert(0, str(ROOT / 'benchmarks'))
sys.path.insert(0, str(ROOT / 'python'))
from hardening.statistics import SCORE_VERSION, adjust_family

spec = importlib.util.spec_from_file_location('hardening_subarm', ROOT / 'benchmarks/tcga_bulk/subarm.py')
subarm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(subarm)

METHODS = ('hd_dose', 'hd_lr', 'hd_lr_smooth', 'naive', 'naive_smooth')
CONTRASTS = (('hd_dose', 'hd_lr'), ('hd_dose', 'naive'), ('hd_dose', 'naive_smooth'),
             ('hd_lr_smooth', 'naive_smooth'), ('hd_lr', 'naive'), ('hd_lr_smooth', 'hd_lr'))
ENDPOINTS = ('recall_loss', 'rho_window', 'pearson_window')
SEED = 20261106
FPR = 0.05
DEFAULT_INPUT = DATA_ROOT / 'results/tcga_bulk_20261101_1405/input'


def deterministic_order(ids, salt):
    return sorted(ids, key=lambda x: (hashlib.sha256(f'{salt}:{x}'.encode()).hexdigest(), x))


def weighted_quantile(values, weights, q):
    values, weights = np.asarray(values, float), np.asarray(weights, float)
    if values.shape != weights.shape or values.ndim != 1 or not 0 <= q <= 1:
        raise ValueError('quantile requires equal-length vectors and q in [0,1]')
    ok = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not ok.any():
        return np.nan
    order = np.argsort(values[ok], kind='stable')
    v, w = values[ok][order], weights[ok][order]
    return float(v[min(np.searchsorted(np.cumsum(w), q * w.sum(), side='left'), len(v) - 1)])


def loss_threshold(scores, neutral, calibration_columns, fpr=FPR):
    """Each calibration patient contributes equal neutral-window probability mass.

    Strict '< threshold' keeps exact ties conservative. The achieved FPR is
    reported; no jitter is added to manufacture a nominal 5% error rate.
    """
    scores, neutral = np.asarray(scores, float), np.asarray(neutral, bool)
    if scores.shape != neutral.shape or scores.ndim != 2 or not 0 < fpr < 1:
        raise ValueError('scores/neutral must be matched 2D arrays and fpr in (0,1)')
    idx = np.asarray(calibration_columns, int)
    if np.any((idx < 0) | (idx >= scores.shape[1])):
        raise ValueError('calibration column outside matrix')
    values, weights = [], []
    for j in idx:
        v = scores[:, j][neutral[:, j] & np.isfinite(scores[:, j])]
        if len(v):
            values.append(v)
            weights.append(np.full(len(v), 1.0 / len(v)))
    return weighted_quantile(np.concatenate(values), np.concatenate(weights), fpr) if values else np.nan


def score_loss(scores, loss, neutral, common, threshold, negative=None):
    """Patient-level loss metrics with explicit callability and precision truth.

    ``negative`` supplies confident nonloss truth (focal gains plus neutral
    windows in the profile runner). Its default is neutral, preserving the
    narrower fixture/API contract. Precision excludes unclassified windows;
    their calls are counted separately. An empty call denominator is undefined.
    Calibration still uses neutral windows exclusively in ``loss_threshold``.
    """
    scores = np.asarray(scores, float)
    loss, neutral, common = (np.asarray(x, bool) for x in (loss, neutral, common))
    negative = neutral if negative is None else np.asarray(negative, bool)
    if any(x.shape != scores.shape for x in (loss, neutral, common, negative)):
        raise ValueError('loss scoring masks must match scores')
    if np.any(loss & negative):
        raise ValueError('positive and confident negative loss truth must not overlap')
    out = []
    for j in range(scores.shape[1]):
        finite = np.isfinite(scores[:, j])
        l, n = loss[:, j] & common[:, j], neutral[:, j] & common[:, j]
        calibrated = np.isfinite(threshold)
        hit = scores[:, j] < threshold
        classified = loss[:, j] | negative[:, j]
        callable_classified = classified & common[:, j] & finite
        calls = hit & callable_classified
        n_calls = int(calls.sum()) if calibrated else None
        true_calls = int((calls & loss[:, j]).sum()) if calibrated else None
        out.append(dict(recall_loss=float(hit[l].mean()) if l.any() and calibrated else np.nan,
                        fpr_loss=float(hit[n].mean()) if n.any() and calibrated else np.nan,
                        recall_loss_failure_aware=float(hit[loss[:, j]].mean()) if loss[:, j].any() and calibrated else np.nan,
                        precision_loss=true_calls / n_calls if calibrated and n_calls else np.nan,
                        n_calls=n_calls, n_true_loss_calls=true_calls,
                        n_false_loss_calls=int((calls & negative[:, j]).sum()) if calibrated else None,
                        n_classified_truth=int(classified.sum()), n_classified_callable=int(callable_classified.sum()),
                        classified_truth_coverage=float(callable_classified.sum() / classified.sum()) if classified.any() else np.nan,
                        n_unclassified_truth=int((~classified).sum()),
                        n_unclassified_callable=int((~classified & common[:, j] & finite).sum()),
                        n_unclassified_call_exclusions=int((hit & ~classified & common[:, j] & finite).sum()) if calibrated else None,
                        n_loss=int(loss[:, j].sum()), n_neutral=int(neutral[:, j].sum()),
                        n_common_loss=int(l.sum()), n_common_neutral=int(n.sum()),
                        n_callable_loss=int((loss[:, j] & finite).sum()),
                        n_callable_neutral=int((neutral[:, j] & finite).sum()), threshold=threshold))
    return out


def profiles(inp, hd, truth, genes):
    d = subarm.type_windows('COADREAD', inp, None, *truth, hd_path=hd, fixed_genes=genes)
    chrom = d['gene_positions'].chromosome.to_numpy()
    sm = np.full_like(d['gene_values']['hd_lr'], np.nan)
    for c in pd.unique(chrom):
        m = chrom == c
        sm[m] = median_filter(d['gene_values']['hd_lr'][m], size=(subarm.SMOOTH, 1), mode='nearest')
    d['W']['hd_lr_smooth'] = subarm.window_median(sm, d['window_id'], len(d['warm']))
    return d


def correlations(d, common):
    out = {}
    for m in METHODS:
        r, p = [], []
        for j in range(len(d['tum'])):
            ok = common[:, j]
            x, y = d['W'][m][ok, j], d['D'][ok, j]
            r.append(subarm.rho(x, y))
            p.append(float(pearsonr(x, y).statistic) if len(x) > 5 and np.ptp(y) > 0 and np.ptp(x) > 0
                     else 0.0 if len(x) > 5 and np.ptp(y) > 0 else np.nan)
        out[m] = dict(rho_window=np.asarray(r), pearson_window=np.asarray(p))
    return out


def score_profiles(d, calibration_columns):
    common = np.isfinite(d['D']) & np.logical_and.reduce([np.isfinite(d['W'][m]) for m in METHODS])
    neutral = d['NEU'] & common
    corr = correlations(d, common)
    scores, thresholds, records = {}, {}, []
    for m in METHODS:
        threshold = loss_threshold(d['W'][m], neutral, calibration_columns)
        thresholds[m] = threshold
        rows = score_loss(d['W'][m], d['FL'], d['NEU'], common, threshold,
                          negative=d['FG'] | d['NEU'])
        scores[m] = dict(corr[m], recall_loss=np.asarray([r['recall_loss'] for r in rows]))
        for j, row in enumerate(rows):
            records.append(dict(sample=str(d['tum'][j]), method=m, **row,
                                rho_window=float(corr[m]['rho_window'][j]),
                                pearson_window=float(corr[m]['pearson_window'][j])))
    return scores, thresholds, records, common


def _paired_mean(a, b, idx):
    d = (np.asarray(a) - np.asarray(b))[idx]
    return float(d[np.isfinite(d)].mean()) if np.isfinite(d).any() else np.nan


def nested_comparisons(point, refits, calibration_columns, evaluation_columns, n_boot, seed):
    """Resample reference refits, calibration patients and evaluation patients.

    Reference refits are biological normal-donor bootstrap draws, not new
    independent patients. Each inner draw re-estimates neutral thresholds.
    """
    rng = np.random.default_rng(seed)
    point_scores = score_profiles(point, calibration_columns)[0]
    family_size = len(CONTRASTS) * len(ENDPOINTS)
    records = []
    keys = [(e, a, b) for e in ENDPOINTS for a, b in CONTRASTS]
    estimates = {key: _paired_mean(point_scores[key[1]][key[0]], point_scores[key[2]][key[0]],
                                  evaluation_columns) for key in keys}
    draws = {key: [] for key in keys}
    # Correlations do not depend on calibration thresholds; cache them per normal refit.
    caches = []
    for d in refits:
        common = np.isfinite(d['D']) & np.logical_and.reduce([np.isfinite(d['W'][m]) for m in METHODS])
        caches.append((common, correlations(d, common)))
    for _ in range(n_boot):
        ri = int(rng.integers(len(refits)))
        d = refits[ri]
        common, corr = caches[ri]
        ci = rng.choice(calibration_columns, len(calibration_columns), replace=True)
        ei = rng.choice(evaluation_columns, len(evaluation_columns), replace=True)
        scored = {}
        for m in METHODS:
            threshold = loss_threshold(d['W'][m], d['NEU'] & common, ci)
            loss = score_loss(d['W'][m], d['FL'], d['NEU'], common, threshold,
                              negative=d['FG'] | d['NEU'])
            scored[m] = dict(corr[m], recall_loss=np.asarray([r['recall_loss'] for r in loss]))
        for e, a, b in keys:
            draws[(e, a, b)].append(_paired_mean(scored[a][e], scored[b][e], ei))
    for key in keys:
        e, a, b = key
        boot = np.asarray(draws[key])
        boot = boot[np.isfinite(boot)]
        delta = estimates[key]
        paired = np.isfinite(point_scores[a][e][evaluation_columns]) & np.isfinite(point_scores[b][e][evaluation_columns])
        testable = np.isfinite(delta) and len(boot) >= max(99, int(0.9 * n_boot)) and np.ptp(boot) > 0
        r = dict(score_version=SCORE_VERSION, endpoint=e, a=a, b=b, mean_diff=delta,
                 n=int(paired.sum()), n_missing=int((~paired).sum()), family_size=family_size,
                 p=1.0, ci_low=np.nan, ci_high=np.nan, family_ci_low=np.nan, family_ci_high=np.nan,
                 testable=bool(testable), valid_bootstrap_replicates=len(boot),
                 uncertainty='normal donor, calibration patient, evaluation patient resampling',
                 evidence_role='internal_development', independent_evaluation=False)
        if testable:
            error = boot - delta
            r['p'] = float((1 + np.count_nonzero(np.abs(error) >= abs(delta))) / (len(boot) + 1))
            for prefix, tail in (('', .025), ('family_', .05 / (2 * family_size))):
                lo, hi = delta - np.quantile(error, [1 - tail, tail])
                r[prefix + 'ci_low'], r[prefix + 'ci_high'] = float(lo), float(hi)
        records.append(r)
    return adjust_family(records)


def write_locked_diagnostics(out, data_root=DATA_ROOT):
    """Only previously exposed development scores; no selection from former evaluations."""
    definitions = [
        ('ucan', data_root / 'results/ucan_bulk_frozen1424_20261102_1612/score_20261102_162020/per_sample_rho.tsv',
         ['CRC.SW.U0008.T', 'CRC.SW.U0001.T'], 'raw_counts'),
        ('tcga', data_root / 'results/tcga_bulk_20261101_1405/score_centred/per_sample_rho.tsv',
         ['TCGA-05-4244-01', 'TCGA-05-4249-01'], 'normalized_pseudocounts_processing_diagnostic'),
        ('kinker', data_root / 'results/kinker_field_20261102_0529/eval_haplodose/line_scores.tsv',
         ['42MGBA_CENTRAL_NERVOUS_SYSTEM', '2313287_STOMACH'], 'raw_counts')]
    records = []
    for name, path, ids, scale in definitions:
        table = pd.read_csv(path, sep='\t')
        idcol = 'sample' if 'sample' in table else 'line'
        table = table[table[idcol].isin(ids)].copy()
        if set(table[idcol]) != set(ids):
            raise ValueError(f'locked diagnostic case missing from {path}')
        table['source'], table['source_table'], table['input_scale'] = name, str(path), scale
        records.append(table)
    pd.concat(records, ignore_index=True).to_csv(out / 'locked_diagnostic_scores.tsv', sep='\t', index=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True)
    parser.add_argument('--data-root', type=Path, default=DATA_ROOT,
                        help='live repository data/results root, separate from frozen code')
    parser.add_argument('--input-root', type=Path,
                        help='explicit existing paired raw/Xena input directory; overrides --data-root')
    parser.add_argument('--device', choices=['cpu', 'gpu'], default='cpu')
    parser.add_argument('--calibration-size', type=int, default=24)
    parser.add_argument('--evaluation-size', type=int, default=24)
    parser.add_argument('--reference-size', type=int, default=20)
    parser.add_argument('--reference-refits', type=int, default=8)
    parser.add_argument('--bootstrap', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=SEED)
    args = parser.parse_args(argv)
    args.data_root = args.data_root.resolve()
    if args.input_root is None:
        args.input_root = args.data_root / 'results/tcga_bulk_20261101_1405/input'
    subarm.GTF37 = args.data_root / 'data/external/ensembl_grch37_r75/Homo_sapiens.GRCh37.75.gtf.gz'
    if min(args.calibration_size, args.evaluation_size) < 2 or args.reference_size < 3 or args.reference_refits < 0 or args.bootstrap < 99:
        raise ValueError('require >=2 patients per split, >=3 normals, >=0 refits and >=99 bootstrap draws')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    raw, xena = args.input_root / 'coadread_raw', args.input_root / 'coadread_xena'
    meta = pd.read_csv(raw / 'samples.tsv', sep='\t', dtype=str)
    xm = pd.read_csv(xena / 'samples.tsv', sep='\t', dtype=str)
    for table in (meta, xm):
        subarm.unique_ids(table['sample'], 'source sample metadata')
        if table.patient.isna().any() or table.groupby('patient').role.apply(lambda x: x.eq('tumour').sum()).gt(1).any():
            raise ValueError('one tumour per biological patient is required')
    if set(meta['sample']) != set(xm['sample']):
        raise ValueError('raw and Xena diagnostic must contain identical source samples')
    if not meta.set_index('sample')[['role', 'patient']].sort_index().equals(
            xm.set_index('sample')[['role', 'patient']].sort_index()):
        raise ValueError('raw and Xena source roles/patient identities disagree')
    pos37 = subarm.hg19_genes()
    position_exclusions = pos37.attrs.get('symbol_position_exclusions', [])
    pd.DataFrame(position_exclusions, columns=['source_symbol', 'canonical_symbol', 'chromosome', 'start', 'reason']).to_csv(
        out / 'DNA_position_exclusions.tsv', sep='\t', index=False)
    arms = subarm.hg19_arms()
    tab, seg = subarm.absolute()
    truth = (pos37, tab, seg, *arms)
    # Select without outcome/likelihood access. All derivatives stay in this new root.
    eligible_tumour = meta.role.eq('tumour') & meta['sample'].isin(tab.index)
    patients = deterministic_order(meta.patient[eligible_tumour].unique(), f'{args.seed}:patients')
    need = args.calibration_size + args.evaluation_size
    if len(patients) < need:
        raise ValueError('not enough independent existing patients for the requested development split')
    cal = set(patients[:args.calibration_size])
    evaluation = set(patients[args.calibration_size:need])
    normal_meta = meta[meta.role.eq('normal') & ~meta.patient.isin(evaluation)].copy()
    if normal_meta.patient.duplicated().any():
        raise ValueError('normal reference bootstrap requires one sample per biological donor')
    normal_ids = deterministic_order(normal_meta['sample'], f'{args.seed}:normals')[:args.reference_size]
    if len(normal_ids) < args.reference_size:
        raise ValueError('not enough normal donors independent of evaluation patients')
    chosen = meta[(meta.role.eq('tumour') & meta.patient.isin(cal | evaluation)) | meta['sample'].isin(normal_ids)].copy()
    chosen['split'] = np.where(chosen.role.eq('normal'), 'reference',
                               np.where(chosen.patient.isin(cal), 'calibration', 'evaluation'))
    chosen.to_csv(out / 'split_manifest.tsv', sep='\t', index=False)
    selected = chosen['sample'].to_numpy()
    X = pd.read_csv(raw / 'counts.tsv', sep='\t', index_col=0, usecols=['gene', *selected])
    X = X.loc[:, selected]
    position = pd.read_csv(raw / 'positions.tsv', sep='\t')
    raw_gene_names = X.index.to_numpy(dtype=str)
    X.index = subarm.canonical_symbols(raw_gene_names, 'raw count genes')
    X.index.name = 'gene'
    position['gene'] = subarm.canonical_symbols(position.gene, 'raw position symbols')
    pd.DataFrame(dict(source_gene=raw_gene_names, canonical_gene=X.index)).to_csv(
        out / 'raw_gene_identity_map.tsv', sep='\t', index=False)
    values = X.to_numpy(float)
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values != np.floor(values)):
        raise ValueError('raw STAR development arm must contain nonnegative integer counts')
    if np.any(values.sum(axis=0) <= 0):
        raise ValueError('raw STAR development samples require positive count totals')
    normal = chosen.role.eq('normal').to_numpy()
    cpm = values / values.sum(axis=0, keepdims=True) * 1e6
    keep = (cpm[:, normal].mean(axis=1) >= subarm.MIN_CPM) & (values[:, normal] > 0).all(axis=1)
    keep &= X.index.isin(position.gene)
    keep &= X.index.isin(pos37.index)
    frozen_genes = pd.read_csv(xena / 'haplodose/genes.tsv.gz', sep='\t', usecols=['gene']).gene
    keep &= X.index.isin(subarm.canonical_symbols(frozen_genes, 'frozen Xena output symbols'))
    from haplodose.bins import RIBOSOMAL
    keep &= ~X.index.to_series().str.match(RIBOSOMAL).to_numpy()
    cen = pd.read_csv(subarm.CEN38, sep='\t', comment='#').set_index('chromosome').p_arm_end
    gene_position = position.set_index('gene').reindex(X.index)
    gene_chromosome = 'chr' + gene_position.chromosome.astype(str).str.removeprefix('chr')
    boundary = gene_chromosome.map(cen)
    keep &= ~((gene_position.start < boundary) & (gene_position.end > boundary)).to_numpy()
    X = X.loc[keep]
    if len(X) < subarm.MIN_GENES:
        raise ValueError('too few genes with fixed normal-panel support')
    genes = X.index.to_numpy()
    position = position[position.gene.isin(genes)]
    # This fixed universe is reused across refits; it cannot expand after outcomes are seen.
    (out / 'gene_universe.txt').write_text('\n'.join(genes) + '\n')
    from haplodose.bulk import run_bulk
    refits = []
    rng = np.random.default_rng(args.seed)
    nref = args.reference_refits
    for refit in range(nref + 1):
        directory = out / f'raw_refit_{refit:03d}'
        directory.mkdir()
        counts, samples = X.copy(), chosen.copy()
        if refit:
            draw = rng.choice(normal_ids, len(normal_ids), replace=True)
            tum = chosen[chosen.role.eq('tumour')].copy()
            reference = normal_meta.set_index('sample').loc[draw].reset_index()
            reference['source_sample'] = reference['sample']
            reference['sample'] = [f'reference_bootstrap_{j}' for j in range(len(draw))]
            reference['split'] = 'reference'
            counts = pd.concat([X.loc[:, tum['sample']],
                                pd.DataFrame(X.loc[:, draw].to_numpy(), index=X.index, columns=reference['sample'])], axis=1)
            samples = pd.concat([tum, reference], ignore_index=True)
        counts.to_csv(directory / 'counts.tsv', sep='\t')
        samples.to_csv(directory / 'samples.tsv', sep='\t', index=False)
        position.to_csv(directory / 'positions.tsv', sep='\t', index=False)
        run_bulk(directory / 'counts.tsv', directory / 'samples.tsv', directory / 'positions.tsv',
                 directory / 'haplodose', min_cpm=0.0, device=args.device, counts_kind='raw',
                 genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')
        d = profiles(directory, directory / 'haplodose', truth, genes)
        if refits and not np.array_equal(d['tum'], refits[0]['tum']):
            raise ValueError('refit changed target specimen order')
        refits.append(d)
        print(f'raw reference refit {refit}/{nref} complete', flush=True)
    patient = chosen.set_index('sample').patient.reindex(refits[0]['tum'])
    ci = np.flatnonzero(patient.isin(cal).to_numpy())
    ei = np.flatnonzero(patient.isin(evaluation).to_numpy())
    if len(ci) != len(cal) or len(ei) != len(evaluation):
        raise ValueError('selected patients lack matched DNA truth; freeze a new metadata eligibility manifest')
    score, thresholds, records, common = score_profiles(refits[0], ci)
    result = pd.DataFrame(records)
    result['split'] = result['sample'].map(chosen.set_index('sample').split)
    result.to_csv(out / 'raw_patient_scores.tsv', sep='\t', index=False)
    np.savez_compressed(out / 'raw_support_and_truth.npz', common=common, loss=refits[0]['FL'],
                        neutral=refits[0]['NEU'], dna=refits[0]['D'], samples=refits[0]['tum'])
    comparisons = nested_comparisons(refits[0], refits[1:] or refits, ci, ei, args.bootstrap, args.seed)
    pd.DataFrame(comparisons).to_csv(out / 'raw_comparisons.tsv', sep='\t', index=False)
    # Processed-expression diagnostic: score existing outputs only, with the same patients and genes.
    xdiag = out / 'xena_diagnostic_input'
    xdiag.mkdir()
    xv = pd.read_csv(xena / 'counts.tsv', sep='\t', index_col=0, usecols=['gene', *selected])
    xv.index = subarm.canonical_symbols(xv.index, 'Xena genes')
    xv.index.name = 'gene'
    if not set(genes) <= set(xv.index):
        raise ValueError('paired processing diagnostic is missing a fixed-universe gene')
    xv.loc[genes, selected].to_csv(xdiag / 'counts.tsv', sep='\t')
    chosen.to_csv(xdiag / 'samples.tsv', sep='\t', index=False)
    position.to_csv(xdiag / 'positions.tsv', sep='\t', index=False)
    xd = profiles(xdiag, xena / 'haplodose', truth, genes)
    xpatient = chosen.set_index('sample').patient.reindex(xd['tum'])
    keep_sample = xpatient.isin(cal | evaluation).to_numpy()
    xd['tum'] = xd['tum'][keep_sample]
    for key in ('D', 'FG', 'FL', 'NEU'):
        xd[key] = xd[key][:, keep_sample]
    xd['W'] = {m: values[:, keep_sample] for m, values in xd['W'].items()}
    xci = np.flatnonzero(xpatient[keep_sample].isin(cal).to_numpy())
    _, xthresholds, xrecords, _ = score_profiles(xd, xci)
    xr = pd.DataFrame(xrecords)
    xr['evidence_role'] = 'frozen_normalized_Xena_processing_diagnostic_only'
    xr['reference_note'] = np.where(xr.method.str.startswith('naive'),
                                    'same selected normal donors and gene universe as raw STAR',
                                    'frozen HaploDose output used its original larger normal panel; confounded diagnostic')
    xr.to_csv(out / 'xena_processing_diagnostic.tsv', sep='\t', index=False)
    pairing = result[result.method.str.startswith('naive')].merge(
        xr[xr.method.str.startswith('naive')], on=['sample', 'method'], suffixes=('_raw', '_xena'), validate='one_to_one')
    pairing.to_csv(out / 'matched_processing_baseline_pairs.tsv', sep='\t', index=False)
    write_locked_diagnostics(out, args.data_root)
    summary = dict(score_version=SCORE_VERSION, evidence_role='internal_development',
                   independent_evaluation=False, confirmatory_claim=False,
                   code_root=str(ROOT), data_root=str(args.data_root), input_root=str(args.input_root),
                   seed=args.seed, methods=METHODS, contrasts=CONTRASTS, endpoints=ENDPOINTS,
                   n_calibration_patients=len(ci), n_evaluation_patients=len(ei),
                   n_reference_donors=len(normal_ids), normal_reference_refits=nref,
                   reference_uncertainty='included by finite normal-donor refits' if nref else 'not included',
                   limitation='small development falsifier; finite reference-refit ensemble is not validated CI coverage',
                   threshold_scope='separate calibration patients; patient-weighted neutral windows',
                   fpr_loss_target=FPR, tie_policy='conservative strict inequality, no jitter',
                   precision_scope='descriptive; common callable focal-loss/focal-gain/neutral truth; unclassified excluded and counted',
                   family='Holm over 18 fixed mean-effect contrasts; Bonferroni basic-bootstrap intervals',
                   input_scale=dict(raw='raw STAR counts', xena='frozen normalized-pseudocount processing diagnostic'),
                   raw_thresholds=thresholds, xena_diagnostic_thresholds=xthresholds,
                   common_callable_windows=int(common.sum()), comparisons=comparisons)
    summary['ambiguous_DNA_position_loci_excluded'] = len(position_exclusions)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({k: v for k, v in summary.items() if k != 'comparisons'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
