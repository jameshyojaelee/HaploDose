"""Does Haplodose add anything below arm level? TCGA bulk RNA vs ABSOLUTE at ~10 Mb gene-anchored windows.

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
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

H = Path('/gpfs/commons/home/jameslee/Haplodose')
CN = Path('/gpfs/commons/home/jameslee/CRC/data/external/tcga_pancan_cn')
CYTO_HG19 = Path('/gpfs/commons/home/jameslee/CRC/data/external/bian_cnv_annotations/cytoBand.txt.gz')
GTF37 = H / 'data/external/ensembl_grch37_r75/Homo_sapiens.GRCh37.75.gtf.gz'
CEN38 = H / 'python/haplodose/data/centromeres_hg38.tsv'
WINDOW_BP = 10_000_000
MIN_GENES = 10
MIN_CPM = 1.0
SMOOTH = 101
FPR = 0.05
FPR_GRID = (0.01, 0.02, 0.05, 0.10, 0.20)
SEED = 20261101
TIEBREAK_MODE = 'lr'   # 'lr' or 'random' (--tiebreak)
TIEBREAK = 1e-6   # hd_dose for focal recall: dose + 1e-6 x Haplodose log-ratio window value (orders tied doses)
METHODS = ['hd_dose', 'hd_lr', 'naive', 'naive_smooth', 'infercnvpy']


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
    return g[~g.gene.duplicated(keep=False)].set_index('gene')


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
        for c, sc in by[s].groupby('chromosome'):
            m = np.flatnonzero(chrom == c)
            if not len(m):
                continue
            st, en = sc.Start.to_numpy(), sc.End.to_numpy()
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
        FG[w] = (fg[m].sum(axis=0) >= 0.5 * n) & (n > 0)
        FL[w] = (fl[m].sum(axis=0) >= 0.5 * n) & (n > 0)
        NEU[w] = (neu[m].sum(axis=0) == n) & (n > 0)
    return W, FG, FL, NEU


def type_windows(t, inp, icnv, g37, tab, seg, pend, size):
    t0 = time.time()
    inp = Path(inp)
    X = pd.read_csv(inp / 'counts.tsv', sep='\t', index_col=0)
    s = pd.read_csv(inp / 'samples.tsv', sep='\t').set_index('sample').loc[X.columns]
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').set_index('gene')
    pos['chromosome'] = 'chr' + pos.chromosome.astype(str)
    cen = pd.read_csv(CEN38, sep='\t', comment='#').set_index('chromosome').p_arm_end
    cpm = X / X.sum(axis=0) * 1e6
    normal = s.role.eq('normal').to_numpy()
    keep = (cpm.loc[:, normal].mean(axis=1) >= MIN_CPM) & cpm.index.isin(pos.index)
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
    hd = inp / 'haplodose'
    hs = pd.read_csv(hd / 'samples.tsv', sep='\t')
    lr = np.load(hd / 'log_ratio.npy').astype(np.float64) + hs.centre.to_numpy()[:, None]
    lr -= np.median(lr, axis=1, keepdims=True)
    row = pd.Index(hs['sample']).get_indexer(tum)
    nb = lr.shape[1]
    dose = np.zeros((len(hs), nb), np.float32)
    sg = pd.read_csv(hd / 'segments.tsv', sep='\t', usecols=['sample', 'first_bin', 'last_bin', 'log_dose'])
    srow = pd.Index(hs['sample']).get_indexer(sg['sample'])
    for r_, a, b, v in zip(srow, sg.first_bin, sg.last_bin, sg.log_dose):
        dose[r_, a:b + 1] = v
    gb = pd.read_csv(hd / 'genes.tsv.gz', sep='\t').set_index('gene').bin.reindex(genes)
    has = gb.notna().to_numpy()
    bi = gb.fillna(0).astype(int).to_numpy()
    for name, M in (('hd_lr', lr), ('hd_dose', dose)):
        A = M[row][:, bi].T.astype(np.float32)
        A[~has] = np.nan
        G[name] = A
    # infercnvpy windows mapped to genes.
    z = np.load(Path(icnv) / 'infercnvpy_windows.npz')
    zrow = pd.Index(z['samples']).get_indexer(tum)
    gw = pd.read_csv(Path(icnv) / 'infercnvpy_gene_window.tsv', sep='\t').set_index('gene').bin.reindex(genes)
    has = gw.notna().to_numpy()
    A = z['Z'][zrow][:, gw.fillna(0).astype(int).to_numpy()].T.astype(np.float32)
    A[~has] = np.nan
    G['infercnvpy'] = A

    wid = windows(p[['chromosome', 'start', 'arm']])
    nw = wid.max() + 1
    warm = pd.Series(p.arm.to_numpy()).groupby(wid).first().reindex(range(nw)).to_numpy()
    rel, fg, fl, neu = dna_genes(tum, genes, g37, tab, seg, pend, size)
    D, FG, FL, NEU = window_flags(rel, fg, fl, neu, wid, nw)
    Wm = {m: window_median(G[m], wid, nw) for m in METHODS}
    info = dict(type=t, n_tumours=len(tum), n_genes=len(genes), n_windows=int(nw),
                genes_with_hg19=int(g37.index.isin(genes).sum()), seconds_windows=round(time.time() - t0, 1))
    return dict(tum=tum, warm=warm, D=D, FG=FG, FL=FL, NEU=NEU, W=Wm, info=info)


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
    rows = []
    Dr = within_arm(d['D'], d['warm'])
    V = dict(d['W'])
    if TIEBREAK_MODE == 'lr':
        V['hd_dose'] = d['W']['hd_dose'] + TIEBREAK * np.nan_to_num(d['W']['hd_lr'])
    else:   # seeded uniform jitter far below any real difference, for the methods with exact ties
        rng = np.random.default_rng(SEED)
        for m in ('hd_dose', 'infercnvpy'):
            V[m] = d['W'][m] + TIEBREAK * rng.random(d['W'][m].shape)
    thr = {}
    for m in METHODS:
        v = V[m][d['NEU']]
        v = v[~np.isnan(v)]
        thr[m] = (np.quantile(v, 1 - FPR), np.quantile(v, FPR)) if len(v) else (np.nan, np.nan)
    Wr = {m: within_arm(d['W'][m], d['warm']) for m in METHODS}
    for j, s in enumerate(d['tum']):
        r = dict(type=t, sample=s, n_focal_gain=int(d['FG'][:, j].sum()), n_focal_loss=int(d['FL'][:, j].sum()),
                 n_neutral=int(d['NEU'][:, j].sum()))
        for m in METHODS:
            x, y = d['W'][m][:, j], d['D'][:, j]
            ok = ~np.isnan(x) & ~np.isnan(y)
            r[f'rho_window_{m}'] = rho(x[ok], y[ok])
            xr, yr = Wr[m][:, j], Dr[:, j]
            ok = ~np.isnan(xr) & ~np.isnan(yr)
            r[f'rho_within_arm_{m}'] = rho(xr[ok], yr[ok])
            g, l = d['FG'][:, j], d['FL'][:, j]
            x = V[m][:, j]
            r[f'recall_gain_{m}'] = float((x[g] > thr[m][0]).mean()) if g.any() else np.nan
            r[f'recall_loss_{m}'] = float((x[l] < thr[m][1]).mean()) if l.any() else np.nan
            r[f'fpr_neutral_{m}'] = float(((x[d['NEU'][:, j]] > thr[m][0]) | (x[d['NEU'][:, j]] < thr[m][1])).mean()) if d['NEU'][:, j].any() else np.nan
        rows.append(r)
    pooled = dict(type=t)
    curve_win, curve_pat = [], []
    for m in METHODS:
        x = V[m]
        nv = x[d['NEU']]
        nv = nv[~np.isnan(nv)]
        pooled[f'realised_fpr_gain_{m}'] = float((nv > thr[m][0]).mean())
        pooled[f'realised_fpr_loss_{m}'] = float((nv < thr[m][1]).mean())
        for f in FPR_GRID:
            hi, lo = np.quantile(nv, 1 - f), np.quantile(nv, f)
            curve_win.append(dict(type=t, method=m, fpr=f, realised_fpr_gain=int((nv > hi).sum()), realised_fpr_loss=int((nv < lo).sum()),
                                  neutral=len(nv), gain_hits=int((x[d['FG']] > hi).sum()), gain_windows=int(d['FG'].sum()),
                                  loss_hits=int((x[d['FL']] < lo).sum()), loss_windows=int(d['FL'].sum())))
            for j, s_ in enumerate(d['tum']):
                g, l = d['FG'][:, j], d['FL'][:, j]
                if g.any() or l.any():
                    curve_pat.append(dict(type=t, sample=s_, method=m, fpr=f,
                                          recall_gain=float((x[g, j] > hi).mean()) if g.any() else np.nan,
                                          recall_loss=float((x[l, j] < lo).mean()) if l.any() else np.nan))
        pooled[f'pooled_recall_gain_{m}'] = float((x[d['FG']] > thr[m][0]).mean()) if d['FG'].any() else np.nan
        pooled[f'pooled_recall_loss_{m}'] = float((x[d['FL']] < thr[m][1]).mean()) if d['FL'].any() else np.nan
    pooled['focal_gain_windows'] = int(d['FG'].sum())
    pooled['focal_loss_windows'] = int(d['FL'].sum())
    pooled['neutral_windows'] = int(d['NEU'].sum())
    return pd.DataFrame(rows), pooled, pd.DataFrame(curve_win), pd.DataFrame(curve_pat)


def tests(per, t):
    out = []
    for metric in ['rho_window', 'rho_within_arm', 'recall_gain', 'recall_loss']:
        rows = []
        for a, b in itertools.combinations(METHODS, 2):
            x = per[[f'{metric}_{a}', f'{metric}_{b}']].dropna()
            diff = x.iloc[:, 0] - x.iloc[:, 1]
            p = wilcoxon(x.iloc[:, 0], x.iloc[:, 1]).pvalue if len(x) >= 5 and (diff != 0).sum() >= 5 else np.nan
            rows.append(dict(type=t, metric=metric, a=a, b=b, n=len(x), median_a=x.iloc[:, 0].median(),
                             median_b=x.iloc[:, 1].median(), mean_a=x.iloc[:, 0].mean(), mean_b=x.iloc[:, 1].mean(),
                             a_higher=int((diff > 0).sum()), b_higher=int((diff < 0).sum()), p=p))
        r = pd.DataFrame(rows)
        ok = r.p.notna()
        r['p_holm'] = np.nan
        if ok.any():
            r.loc[ok, 'p_holm'] = multipletests(r.p[ok], method='holm')[1]
        out.append(r)
    return pd.concat(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('types', nargs='+', help='TYPE=INPUT_DIR:ICNV_DIR')
    ap.add_argument('--tiebreak', choices=['lr', 'random'], default='lr',
                    help='lr: dose + 1e-6 x log-ratio window value; random: seeded 1e-6 jitter on hd_dose and infercnvpy')
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
        infos.append(d['info'] | dict(seconds_infercnvpy=float((Path(icnv) / 'infercnvpy_seconds.txt').read_text())))
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
    (out / 'params.json').write_text(json.dumps(dict(window_bp=WINDOW_BP, min_genes=MIN_GENES, smooth_genes=SMOOTH,
                                                     fpr=FPR, fpr_grid=FPR_GRID, dose_tiebreak=TIEBREAK, tiebreak_mode=TIEBREAK_MODE, seed=SEED, fpr_set_on='DNA-neutral windows of the type', methods=METHODS,
                                                     hg19_genes=str(GTF37)), indent=1))
    pd.set_option('display.width', 250)
    print(summ[['type', 'n_tumours'] + [f'median_rho_within_arm_{m}' for m in METHODS]].round(3).to_string(index=False))


if __name__ == '__main__':
    main()
