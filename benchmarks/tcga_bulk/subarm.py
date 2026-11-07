"""Does Haplodose add anything below arm level? TCGA bulk RNA vs ABSOLUTE at ~10 Mb gene-anchored windows.

Version 2 uses validated specimen alignment, complete mapped DNA truth, common
callable windows, separate failure-aware recall and paired-mean bootstrap inference.
The historical protocol below describes version 1, preserved in frozen result
copies and the hardening baseline. New executions are retrospective diagnostics;
same-cohort DNA-oracle thresholds are not deployable false-positive calibration.

Same run as score.py (results/tcga_bulk_20261101_1405), no re-profiling except infercnvpy windows
(subarm_infercnvpy.py, same settings as the arm run).

Windows (gene-anchored, so no liftover): the naive gene set of each type (mean normal CPM >= 1, protein coding,
GRCh38 positions), sorted by position within each hg38 arm; a window closes when the next gene starts > 10 Mb
after the window's first gene. Windows with < MIN_GENES genes are dropped.
DNA per window: each gene's hg19 position from Ensembl GRCh37 release 75 (by gene symbol, same chromosome;
data/external/ensembl_grch37_r75), the ABSOLUTE segment covering that position, log2(max(TCN, 0.25) / ploidy);
window value = median over its genes. Focal: a gene is in a focal event when its segment is shorter than 25%
of its hg19 arm and |TCN - ploidy| >= 1; a window is a focal gain (loss) when >= half its genes are in focal
gains (losses). A window is DNA-neutral when all its genes have |TCN - ploidy| < 0.5.
RNA methods, each reduced to gene values (a gene takes its bin's or window's value), window = median over genes:
  hd_dose        Haplodose HMM log_dose per bin (segments.tsv)
  hd_lr          Haplodose log ratio per bin, plain-median centred (as diag_centre.py)
  naive          gene log2(CPM + 1) - mean normal, minus sample median (naive.py)
  naive_smooth   naive, then running median of 101 genes within each chromosome
  infercnvpy     infercnvpy windows (gene -> its infercnvpy window)
Metrics, unit = patient (one primary tumour each):
  rho_window      Spearman across windows, RNA vs DNA
  rho_within_arm  same after subtracting each arm's median window value from RNA and from DNA (sub-arm signal only)
  focal recall    thresholds set on DNA-neutral windows: per type and method, the 95th (gain) / 5th (loss)
                  percentile of window values over all DNA-neutral windows of that type's tumours (FPR 5% on
                  DNA-neutral windows, not on normals). Reported pooled over windows and per patient (patients
                  with >= 1 focal window of that direction). The HMM dose is piecewise constant with a spike at 0,
                  so for recall its ties are broken by adding TIEBREAK x the Haplodose log-ratio window value; every
                  method then sits at the same realised FPR per side on DNA-neutral windows (realised rates are
                  reported). Recall-vs-FPR curve over FPR_GRID per side (curve_*.tsv).
Paired Wilcoxon signed-rank per metric, Holm over the 10 method pairs within each type (and within pooled).
Usage: subarm.py OUT TYPE=INPUT_DIR:ICNV_DIR ...
"""
import argparse
import gzip
import itertools
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.stats import spearmanr

H = Path(__file__).resolve().parents[2]
DATA_ROOT = Path(os.environ.get('HAPLODOSE_DATA_ROOT', str(H)))
sys.path.insert(0, str(H / 'benchmarks'))
from hardening.statistics import SCORE_VERSION, adjust_family, paired_mean_comparison
CN = Path('/gpfs/commons/home/jameslee/CRC/data/external/tcga_pancan_cn')
CYTO_HG19 = Path('/gpfs/commons/home/jameslee/CRC/data/external/bian_cnv_annotations/cytoBand.txt.gz')
GTF37 = DATA_ROOT / 'data/external/ensembl_grch37_r75/Homo_sapiens.GRCh37.75.gtf.gz'
CEN38 = H / 'python/haplodose/data/centromeres_hg38.tsv'
WINDOW_BP = 10_000_000
MIN_GENES = 10
MIN_CPM = 1.0
SMOOTH = 101
FPR = 0.05
FPR_GRID = (0.01, 0.02, 0.05, 0.10, 0.20)
SEED = 20261101
TIEBREAK_MODE = 'none'   # legacy 'lr' or 'random' remain explicit --tiebreak options
TIEBREAK = 1e-6   # hd_dose for focal recall: dose + 1e-6 x Haplodose log-ratio window value (orders tied doses)
METHODS = ['hd_dose', 'hd_lr', 'naive', 'naive_smooth', 'infercnvpy']


def unique_ids(values, label):
    ids = pd.Index(values)
    if ids.has_duplicates or ids.isna().any() or any(not str(x).strip() for x in ids):
        raise ValueError(f'{label} must have nonempty, unique IDs')
    return ids


def _symbol_labels(values, label):
    raw = pd.Series(values, dtype=object)
    if raw.isna().any() or raw.map(lambda x: not isinstance(x, str) or not x.strip()).any():
        raise ValueError(f'{label} must contain nonempty strings')
    canonical = raw.astype(str).str.strip().str.upper()
    if canonical.str.lower().isin(['nan', 'none', 'null', 'na']).any():
        raise ValueError(f'{label} must not contain missing-value placeholders')
    return canonical.to_numpy(dtype=str)


def canonical_symbols(values, label='gene symbols'):
    """Declared symbol namespace: match package strip/uppercase identity without aggregation."""
    canonical = _symbol_labels(values, label)
    unique_ids(canonical, f'canonical {label}')
    return canonical


def unambiguous_position_symbols(table):
    """Extend the existing exclude-all ambiguous annotation policy to canonical symbols.

    Position annotations are not count observations: a name with several possible
    loci is unavailable as truth. Exclude every candidate locus, report it, and
    never aggregate counts or select a locus from the RNA/DNA outcome.
    """
    canonical = _symbol_labels(table.index, 'DNA position symbols')
    ambiguous = pd.Index(canonical).duplicated(keep=False)
    excluded = table.loc[ambiguous].copy()
    excluded.insert(0, 'source_symbol', table.index[ambiguous])
    excluded.insert(1, 'canonical_symbol', canonical[ambiguous])
    excluded['reason'] = 'ambiguous_canonical_symbol_position_all_loci_excluded'
    out = table.loc[~ambiguous].copy()
    out.index = canonical[~ambiguous]
    unique_ids(out.index, 'unambiguous canonical DNA position symbols')
    out.attrs['symbol_position_exclusions'] = excluded.to_dict('records')
    return out


def align_rows(matrix, source_ids, target_ids):
    """Missing output IDs become NaN rows, never Python's -1 (last) row."""
    source = unique_ids(source_ids, 'method samples')
    target = unique_ids(target_ids, 'target samples')
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != len(source):
        raise ValueError('method matrix rows must match its sample IDs')
    row = source.get_indexer(target)
    out = np.full((len(target), matrix.shape[1]), np.nan)
    present = row >= 0
    out[present] = matrix[row[present]]
    return out


def gene_values(matrix, mapping, genes):
    """Map columns without silently coercing invalid or duplicated gene/bin IDs."""
    unique_ids(mapping.index, 'gene mapping')
    bins = mapping.reindex(genes)
    has = bins.notna().to_numpy()
    bi = bins.fillna(0).to_numpy(float)
    if np.any(has & ((bi < 0) | (bi >= matrix.shape[1]) | (bi != np.floor(bi)))):
        raise ValueError('gene mapping contains an out-of-range or noninteger bin')
    out = np.full((len(genes), matrix.shape[0]), np.nan)
    if has.any():
        out[has] = matrix[:, bi[has].astype(int)].T
    return out


def hg19_genes():
    rows = []
    with gzip.open(GTF37, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 9)
            if f[2] != 'gene':
                continue
            m = re.search(r'gene_name "([^"]+)"', f[8])
            if m:
                rows.append((m.group(1), 'chr' + f[0], int(f[3])))
    g = pd.DataFrame(rows, columns=['gene', 'chromosome', 'start'])
    g = g[g.chromosome.isin([f'chr{i}' for i in range(1, 23)])]
    return unambiguous_position_symbols(g.set_index('gene'))


def hg19_arms():
    with gzip.open(CYTO_HG19, 'rt') as fh:
        c = pd.read_csv(fh, sep='\t', header=None, names=['chromosome', 'start', 'end', 'band', 'stain'])
    pend = c[c.stain.eq('acen') & c.band.str.startswith('p')].set_index('chromosome').end
    size = c.groupby('chromosome').end.max()
    return pend, size


def absolute():
    tab = pd.read_csv(CN / 'TCGA_mastercalls.abs_tables_JSedit.fixed.txt', sep='\t')
    tab = tab[tab['call status'].eq('called')].set_index('array')
    seg = pd.read_csv(CN / 'TCGA_mastercalls.abs_segtabs.fixed.txt', sep='\t',
                      usecols=['Sample', 'Chromosome', 'Start', 'End', 'Modal_Total_CN'])
    seg = seg[seg.Sample.isin(tab.index) & seg.Modal_Total_CN.notna()]
    seg['chromosome'] = 'chr' + seg.Chromosome.astype(float).astype(int).astype(str)
    return tab, seg.sort_values(['Sample', 'chromosome', 'Start'])


def windows(pos38):
    """pos38: genes (index) with chromosome, start, arm, sorted. Returns window id per gene (-1 = dropped)."""
    wid = np.full(len(pos38), -1)
    k = 0
    for arm, g in pos38.groupby('arm', sort=False):
        idx = np.flatnonzero(pos38.arm.to_numpy() == arm)
        starts = g.start.to_numpy()
        first = 0
        groups = []
        for i in range(1, len(idx) + 1):
            if i == len(idx) or starts[i] - starts[first] > WINDOW_BP:
                groups.append(idx[first:i])
                first = i
        for gi in groups:
            if len(gi) >= MIN_GENES:
                wid[gi] = k
                k += 1
    return wid


def window_median(G, wid, nw):
    """G genes x samples (NaN allowed) -> windows x samples, median over genes."""
    out = np.full((nw, G.shape[1]), np.nan, np.float32)
    for w in range(nw):
        out[w] = np.nanmedian(G[wid == w], axis=0)
    return out


def dna_genes(samples, genes, g37, tab, seg, pend, size):
    """Per tumour and gene: rel log2(TCN/ploidy), focal gain/loss, neutral. Arrays genes x samples."""
    unique_ids(tab.index, 'DNA sample table')
    unique_ids(g37.index, 'DNA gene positions')
    chrom = g37.chromosome.reindex(genes).to_numpy()
    p19 = g37.start.reindex(genes).to_numpy()
    shape = (len(genes), len(samples))
    rel = np.full(shape, np.nan, np.float32)
    fg, fl, neu = (np.zeros(shape, bool) for _ in range(3))
    armlen = np.where(p19 < pd.Series(chrom).map(pend).to_numpy(), pd.Series(chrom).map(pend).to_numpy(),
                      pd.Series(chrom).map(size).to_numpy() - pd.Series(chrom).map(pend).to_numpy())
    by = dict(tuple(seg[seg.Sample.isin(samples)].groupby('Sample')))
    for j, s in enumerate(samples):
        if s not in by:
            continue
        ploidy = tab.loc[s, 'ploidy']
        if not np.isfinite(ploidy) or ploidy <= 0:
            raise ValueError(f'invalid DNA ploidy for {s}')
        for c, sc in by[s].groupby('chromosome'):
            m = np.flatnonzero(chrom == c)
            if not len(m):
                continue
            sc = sc.sort_values('Start')
            st, en = sc.Start.to_numpy(), sc.End.to_numpy()
            if not np.isfinite(sc[['Start', 'End', 'Modal_Total_CN']].to_numpy(float)).all():
                raise ValueError(f'nonfinite DNA segment for {s} {c}')
            if np.any(sc.Modal_Total_CN.to_numpy() < 0):
                raise ValueError(f'negative DNA copy number for {s} {c}')
            if np.any(en < st) or np.any(st[1:] < en[:-1]):
                raise ValueError(f'inverted or overlapping DNA segments for {s} {c}')
            k = np.searchsorted(st, p19[m], side='right') - 1
            ok = (k >= 0) & (p19[m] <= en[np.maximum(k, 0)])
            m, k = m[ok], k[ok]
            tcn = sc.Modal_Total_CN.to_numpy()[k]
            length = (en - st)[k]
            rel[m, j] = np.log2(np.maximum(tcn, 0.25) / ploidy)
            focal = length < 0.25 * armlen[m]
            fg[m, j] = focal & (tcn - ploidy >= 1)
            fl[m, j] = focal & (tcn - ploidy <= -1)
            neu[m, j] = np.abs(tcn - ploidy) < 0.5
    return rel, fg, fl, neu


def window_flags(rel, fg, fl, neu, wid, nw):
    has = ~np.isnan(rel)
    W = window_median(rel, wid, nw)
    FG, FL, NEU = (np.zeros((nw, rel.shape[1]), bool) for _ in range(3))
    for w in range(nw):
        m = wid == w
        n = has[m].sum(axis=0)
        complete = (n == int(m.sum())) & (n > 0)
        gain = (fg[m].sum(axis=0) >= 0.5 * n) & complete
        loss = (fl[m].sum(axis=0) >= 0.5 * n) & complete
        FG[w], FL[w] = gain & ~loss, loss & ~gain
        NEU[w] = (neu[m].sum(axis=0) == n) & complete
    return W, FG, FL, NEU


def type_windows(t, inp, icnv, g37, tab, seg, pend, size, hd_path=None, fixed_genes=None):
    t0 = time.time()
    inp = Path(inp)
    X = pd.read_csv(inp / 'counts.tsv', sep='\t', index_col=0)
    X.index = canonical_symbols(X.index, 'count genes')
    unique_ids(X.columns, 'count samples')
    smeta = pd.read_csv(inp / 'samples.tsv', sep='\t')
    unique_ids(smeta['sample'], 'sample metadata')
    s = smeta.set_index('sample').loc[X.columns]
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').set_index('gene')
    pos.index = canonical_symbols(pos.index, 'RNA position symbols')
    g37 = g37.copy()
    g37.index = canonical_symbols(g37.index, 'DNA position symbols')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str).str.removeprefix('chr')
    cen = pd.read_csv(CEN38, sep='\t', comment='#').set_index('chromosome').p_arm_end
    values = X.to_numpy(float)
    if not np.isfinite(values).all() or np.any(values < 0) or np.any(values.sum(axis=0) <= 0):
        raise ValueError('expression inputs must be finite, nonnegative with positive library totals')
    cpm = X / X.sum(axis=0) * 1e6
    normal = s.role.eq('normal').to_numpy()
    if not normal.any():
        raise ValueError('a normal reference panel is required')
    keep = (cpm.loc[:, normal].mean(axis=1) >= MIN_CPM) & cpm.index.isin(pos.index)
    if fixed_genes is not None:
        keep = cpm.index.isin(canonical_symbols(fixed_genes, 'fixed scoring symbols')) & cpm.index.isin(pos.index)
    n_expression_genes = int(np.count_nonzero(keep))
    # Fixed positional availability, not observed CN values, defines the scoring
    # universe. Do not let unmappable genes create mostly incomplete windows.
    mapped = g37.reindex(X.index)
    same_chromosome = mapped.chromosome.eq(pos.chromosome.reindex(X.index))
    keep &= X.index.isin(g37.index) & same_chromosome.to_numpy()
    L = np.log2(cpm[keep] + 1)
    R = L.sub(L.loc[:, normal].mean(axis=1), axis=0)
    R = R - R.median(axis=0)
    p = pos.loc[R.index].copy()
    p['chrom_n'] = p.chromosome.str[3:].astype(int)
    p = p.sort_values(['chrom_n', 'start'])
    p['arm'] = p.chromosome + np.where(p.start.to_numpy() >= p.chromosome.map(cen).to_numpy(), 'q', 'p')
    genes = p.index.to_numpy()
    tum = s.index[s.role.eq('tumour')].to_numpy()
    tum = np.array([x for x in tum if x in tab.index])
    G = {}
    G['naive'] = R.loc[genes, tum].to_numpy(np.float32)
    sm = np.empty_like(G['naive'])
    for c in pd.unique(p.chromosome):
        m = (p.chromosome == c).to_numpy()
        sm[m] = median_filter(G['naive'][m], size=(SMOOTH, 1), mode='nearest')
    G['naive_smooth'] = sm
    # Haplodose: per-bin values mapped to genes.
    hd = Path(hd_path) if hd_path is not None else inp / 'haplodose'
    hs = pd.read_csv(hd / 'samples.tsv', sep='\t')
    unique_ids(hs['sample'], 'HaploDose samples')
    lr = np.load(hd / 'log_ratio.npy').astype(np.float64) + hs.centre.to_numpy()[:, None]
    if lr.ndim != 2 or lr.shape[0] != len(hs):
        raise ValueError('HaploDose log-ratio matrix does not match sample metadata')
    lr -= np.median(lr, axis=1, keepdims=True)
    nb = lr.shape[1]
    dose = np.zeros((len(hs), nb), np.float32)
    sg = pd.read_csv(hd / 'segments.tsv', sep='\t', usecols=['sample', 'first_bin', 'last_bin', 'log_dose'])
    srow = pd.Index(hs['sample']).get_indexer(sg['sample'])
    if np.any(srow < 0):
        raise ValueError('HaploDose segment refers to an unknown sample')
    for r_, a, b, v in zip(srow, sg.first_bin, sg.last_bin, sg.log_dose):
        if not np.isfinite([a, b, v]).all() or a != int(a) or b != int(b) or not 0 <= a <= b < nb:
            raise ValueError('invalid HaploDose segment coordinates or dosage')
        dose[r_, a:b + 1] = v
    gb = pd.read_csv(hd / 'genes.tsv.gz', sep='\t').set_index('gene').bin
    gb.index = canonical_symbols(gb.index, 'HaploDose mapped symbols')
    shared_gene_mapping = gb.reindex(genes).notna().to_numpy()
    for name, M in (('hd_lr', lr), ('hd_dose', dose)):
        G[name] = gene_values(align_rows(M, hs['sample'], tum), gb, genes)
    # infercnvpy windows mapped to genes.
    if icnv is not None:
        with np.load(Path(icnv) / 'infercnvpy_windows.npz') as z:
            aligned = align_rows(z['Z'], z['samples'], tum)
        gw = pd.read_csv(Path(icnv) / 'infercnvpy_gene_window.tsv', sep='\t').set_index('gene').bin
        gw.index = canonical_symbols(gw.index, 'infercnvpy mapped symbols')
        G['infercnvpy'] = gene_values(aligned, gw, genes)
        shared_gene_mapping &= gw.reindex(genes).notna().to_numpy()

    wid = windows(p[['chromosome', 'start', 'arm']])
    nw = int(wid.max() + 1) if len(wid) else 0
    warm = pd.Series(p.arm.to_numpy()).groupby(wid).first().reindex(range(nw)).to_numpy()
    rel, fg, fl, neu = dna_genes(tum, genes, g37, tab, seg, pend, size)
    D, FG, FL, NEU = window_flags(rel, fg, fl, neu, wid, nw)
    failure_aware = {m: window_median(G[m], wid, nw) for m in G}
    shared_G = {m: np.where(shared_gene_mapping[:, None], values, np.nan) for m, values in G.items()}
    Wm = {m: window_median(shared_G[m], wid, nw) for m in G}
    D = window_median(np.where(shared_gene_mapping[:, None], rel, np.nan), wid, nw)
    info = dict(type=t, n_tumours=len(tum), n_genes=len(genes), n_windows=int(nw),
                n_expression_genes=n_expression_genes, n_genes_without_concordant_DNA_position=n_expression_genes - len(genes),
                n_common_mapped_genes=int(shared_gene_mapping.sum()),
                genes_with_hg19=int(g37.index.isin(genes).sum()), seconds_windows=round(time.time() - t0, 1))
    return dict(tum=tum, warm=warm, D=D, FG=FG, FL=FL, NEU=NEU, W=Wm, info=info,
                W_failure_aware=Wm, W_own_feature_diagnostic=failure_aware,
                gene_values=shared_G, genes=genes,
                shared_gene_mapping=shared_gene_mapping, gene_positions=p, window_id=wid)


def within_arm(V, warm):
    out = V.copy()
    for a in pd.unique(warm):
        m = warm == a
        out[m] -= np.nanmedian(V[m], axis=0)
    return out


def rho(x, y):
    """Spearman; NaN if DNA is constant or < 6 windows; 0 if DNA varies but the RNA method is constant
    (e.g. an HMM dose with no change inside any arm: it detected nothing, which counts, not drops out)."""
    if len(x) <= 5 or np.ptp(y) == 0:
        return np.nan
    return 0.0 if np.ptp(x) == 0 else float(spearmanr(x, y)[0])


def per_patient(t, d):
    """Version 2: common callable domain plus separate failure-aware denominators.

    Thresholds here deliberately remain a same-cohort DNA oracle. Independent
    calibration lives in hardening/profile.py; this endpoint cannot support a
    deployable FPR claim. Historical results use the archived version-1 scorer.
    """
    rows = []
    Dr = within_arm(d['D'], d['warm'])
    V = dict(d['W'])
    if TIEBREAK_MODE == 'lr':
        V['hd_dose'] = d['W']['hd_dose'] + TIEBREAK * np.nan_to_num(d['W']['hd_lr'])
    elif TIEBREAK_MODE == 'random':
        rng = np.random.default_rng(SEED)
        for m in ('hd_dose', 'infercnvpy'):
            V[m] = d['W'][m] + TIEBREAK * rng.random(d['W'][m].shape)
    common = np.isfinite(d['D']) & np.logical_and.reduce([np.isfinite(V[m]) for m in METHODS])
    gain, loss, neutral = (d[k] & common for k in ('FG', 'FL', 'NEU'))
    thr = {}
    for m in METHODS:
        v = V[m][neutral]
        thr[m] = (np.quantile(v, 1 - FPR), np.quantile(v, FPR)) if len(v) else (np.nan, np.nan)
    Wr = {m: within_arm(d['W'][m], d['warm']) for m in METHODS}
    for j, s in enumerate(d['tum']):
        r = dict(type=t, sample=s, n_focal_gain=int(d['FG'][:, j].sum()), n_focal_loss=int(d['FL'][:, j].sum()),
                 n_neutral=int(d['NEU'][:, j].sum()),
                 n_common_gain=int(gain[:, j].sum()), n_common_loss=int(loss[:, j].sum()),
                 n_common_neutral=int(neutral[:, j].sum()), n_common_callable=int(common[:, j].sum()))
        for m in METHODS:
            x, y = d['W'][m][:, j], d['D'][:, j]
            ok = common[:, j]
            r[f'rho_window_{m}'] = rho(x[ok], y[ok])
            xr, yr = Wr[m][:, j], Dr[:, j]
            ok = common[:, j] & np.isfinite(xr) & np.isfinite(yr)
            r[f'rho_within_arm_{m}'] = rho(xr[ok], yr[ok])
            g, l, n = gain[:, j], loss[:, j], neutral[:, j]
            x = V[m][:, j]
            calibrated = np.isfinite(thr[m]).all()
            r[f'recall_gain_{m}'] = float((x[g] > thr[m][0]).mean()) if g.any() and calibrated else np.nan
            r[f'recall_loss_{m}'] = float((x[l] < thr[m][1]).mean()) if l.any() and calibrated else np.nan
            r[f'fpr_gain_{m}'] = float((x[n] > thr[m][0]).mean()) if n.any() and calibrated else np.nan
            r[f'fpr_loss_{m}'] = float((x[n] < thr[m][1]).mean()) if n.any() and calibrated else np.nan
            r[f'fpr_neutral_{m}'] = float(((x[n] > thr[m][0]) | (x[n] < thr[m][1])).mean()) if n.any() and calibrated else np.nan
            full_x = d.get('W_failure_aware', d['W'])[m][:, j]
            for direction, truth, hit in (('gain', d['FG'][:, j], full_x > thr[m][0]),
                                           ('loss', d['FL'][:, j], full_x < thr[m][1])):
                r[f'recall_{direction}_failure_aware_{m}'] = float(hit[truth].mean()) if truth.any() and calibrated else np.nan
                r[f'callable_{direction}_{m}'] = int((truth & np.isfinite(full_x)).sum())
            r[f'callable_neutral_{m}'] = int((d['NEU'][:, j] & np.isfinite(full_x)).sum())
        rows.append(r)
    pooled = dict(type=t)
    curve_win, curve_pat = [], []
    for m in METHODS:
        x = V[m]
        nv = x[neutral]
        pooled[f'realised_fpr_gain_{m}'] = float((nv > thr[m][0]).mean()) if len(nv) else np.nan
        pooled[f'realised_fpr_loss_{m}'] = float((nv < thr[m][1]).mean()) if len(nv) else np.nan
        for f in FPR_GRID:
            hi, lo = (np.quantile(nv, 1 - f), np.quantile(nv, f)) if len(nv) else (np.nan, np.nan)
            curve_win.append(dict(type=t, method=m, fpr=f, realised_fpr_gain=int((nv > hi).sum()), realised_fpr_loss=int((nv < lo).sum()),
                                  neutral=len(nv), gain_hits=int((x[gain] > hi).sum()), gain_windows=int(gain.sum()),
                                  loss_hits=int((x[loss] < lo).sum()), loss_windows=int(loss.sum())))
            for j, s_ in enumerate(d['tum']):
                g, l = gain[:, j], loss[:, j]
                if g.any() or l.any():
                    curve_pat.append(dict(type=t, sample=s_, method=m, fpr=f,
                                          recall_gain=float((x[g, j] > hi).mean()) if g.any() and len(nv) else np.nan,
                                          recall_loss=float((x[l, j] < lo).mean()) if l.any() and len(nv) else np.nan))
        pooled[f'pooled_recall_gain_{m}'] = float((x[gain] > thr[m][0]).mean()) if gain.any() and len(nv) else np.nan
        pooled[f'pooled_recall_loss_{m}'] = float((x[loss] < thr[m][1]).mean()) if loss.any() and len(nv) else np.nan
    pooled['focal_gain_windows'] = int(d['FG'].sum())
    pooled['focal_loss_windows'] = int(d['FL'].sum())
    pooled['neutral_windows'] = int(d['NEU'].sum())
    columns = ['type', 'sample', 'n_focal_gain', 'n_focal_loss', 'n_neutral', 'n_common_gain',
               'n_common_loss', 'n_common_neutral', 'n_common_callable']
    for m in METHODS:
        columns += [f'{k}_{m}' for k in ('rho_window', 'rho_within_arm', 'recall_gain', 'recall_loss',
                    'fpr_gain', 'fpr_loss', 'fpr_neutral', 'recall_gain_failure_aware',
                    'recall_loss_failure_aware', 'callable_gain', 'callable_loss', 'callable_neutral')]
    curve_columns = ['type', 'sample', 'method', 'fpr', 'recall_gain', 'recall_loss']
    return (pd.DataFrame(rows, columns=columns), pooled, pd.DataFrame(curve_win),
            pd.DataFrame(curve_pat, columns=curve_columns))


def tests(per, t):
    out = []
    for metric in ['rho_window', 'rho_within_arm', 'recall_gain', 'recall_loss']:
        rows = []
        for a, b in itertools.combinations(METHODS, 2):
            r = paired_mean_comparison(per[f'{metric}_{a}'], per[f'{metric}_{b}'], seed=SEED,
                                       family_size=len(METHODS) * (len(METHODS) - 1) // 2)
            rows.append(dict(type=t, metric=metric, a=a, b=b, **r))
        out.append(pd.DataFrame(adjust_family(rows)))
    return pd.concat(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('types', nargs='+', help='TYPE=INPUT_DIR:ICNV_DIR')
    ap.add_argument('--tiebreak', choices=['none', 'lr', 'random'], default='none',
                    help='lr: dose + 1e-6 x log-ratio window value; random: seeded 1e-6 jitter on hd_dose and infercnvpy')
    ap.add_argument('--input-scale', choices=['raw_counts', 'normalized_pseudocounts', 'unspecified'],
                    default='unspecified', help='Required scientific provenance; unspecified is diagnostic only')
    args = ap.parse_args()
    global TIEBREAK_MODE
    TIEBREAK_MODE = args.tiebreak
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    g37 = hg19_genes()
    pend, size = hg19_arms()
    tab, seg = absolute()
    pers, pooled, infos, tst, cws, cps = [], [], [], [], [], []
    for spec in args.types:
        t, rest = spec.split('=', 1)
        inp, icnv = rest.split(':', 1)
        d = type_windows(t, inp, icnv, g37, tab, seg, pend, size)
        per, pl, cw, cp = per_patient(t, d)
        cws.append(cw)
        cps.append(cp)
        pers.append(per)
        pooled.append(pl)
        seconds = Path(icnv) / 'infercnvpy_seconds.txt'
        infos.append(d['info'] | dict(seconds_infercnvpy=float(seconds.read_text()) if seconds.exists() else np.nan))
        tst.append(tests(per, t))
        print(d['info'], flush=True)
    per = pd.concat(pers, ignore_index=True)
    per.to_csv(out / 'per_patient.tsv.gz', sep='\t', index=False)
    tst.append(tests(per, 'pooled'))
    pd.concat(tst).to_csv(out / 'paired_tests.tsv', sep='\t', index=False)
    pd.DataFrame(infos).to_csv(out / 'type_info.tsv', sep='\t', index=False)
    summ = []
    for t, g in list(per.groupby('type')) + [('pooled', per)]:
        r = dict(type=t, n_tumours=len(g), n_with_focal_gain=int((g.n_focal_gain > 0).sum()),
                 n_with_focal_loss=int((g.n_focal_loss > 0).sum()))
        for m in METHODS:
            r[f'median_rho_window_{m}'] = g[f'rho_window_{m}'].median()
            r[f'median_rho_within_arm_{m}'] = g[f'rho_within_arm_{m}'].median()
            r[f'mean_recall_gain_{m}'] = g[f'recall_gain_{m}'].mean()
            r[f'mean_recall_loss_{m}'] = g[f'recall_loss_{m}'].mean()
            r[f'mean_fpr_neutral_{m}'] = g[f'fpr_neutral_{m}'].mean()
        summ.append(r)
    summ = pd.DataFrame(summ).merge(pd.DataFrame(pooled), on='type', how='left')
    cw, cp = pd.concat(cws), pd.concat(cps)
    cw.to_csv(out / 'curve_windows_by_type.tsv', sep='\t', index=False)
    cp.to_csv(out / 'curve_per_patient.tsv.gz', sep='\t', index=False)
    agg = cw.groupby(['method', 'fpr'])[['realised_fpr_gain', 'realised_fpr_loss', 'neutral', 'gain_hits', 'gain_windows', 'loss_hits', 'loss_windows']].sum()
    curve = pd.DataFrame({'realised_fpr_gain': agg.realised_fpr_gain / agg.neutral, 'realised_fpr_loss': agg.realised_fpr_loss / agg.neutral,
                          'pooled_window_recall_gain': agg.gain_hits / agg.gain_windows, 'pooled_window_recall_loss': agg.loss_hits / agg.loss_windows})
    pm = cp.groupby(['method', 'fpr']).agg(patient_mean_recall_gain=('recall_gain', 'mean'), n_gain=('recall_gain', 'count'),
                                           patient_mean_recall_loss=('recall_loss', 'mean'), n_loss=('recall_loss', 'count'))
    curve.join(pm).reset_index().to_csv(out / 'curve_pooled.tsv', sep='\t', index=False)
    summ.to_csv(out / 'summary.tsv', sep='\t', index=False)
    (out / 'params.json').write_text(json.dumps(dict(score_version=SCORE_VERSION,
                                                     input_scale=args.input_scale, diagnostic_only=True,
                                                     threshold_scope='same-cohort DNA oracle, not external calibration',
                                                     callable_scope='common finite windows; failure-aware recall separate',
                                                     truth_coverage='all window genes require DNA coverage; mixed gain/loss ties excluded',
                                                     uncertainty='conditional on reference panel and oracle threshold',
                                                     window_bp=WINDOW_BP, min_genes=MIN_GENES, smooth_genes=SMOOTH,
                                                     fpr=FPR, fpr_grid=FPR_GRID, dose_tiebreak=TIEBREAK, tiebreak_mode=TIEBREAK_MODE, seed=SEED, fpr_set_on='DNA-neutral windows of the type', methods=METHODS,
                                                     hg19_genes=str(GTF37)), indent=1))
    pd.set_option('display.width', 250)
    print(summ[['type', 'n_tumours'] + [f'median_rho_within_arm_{m}' for m in METHODS]].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
