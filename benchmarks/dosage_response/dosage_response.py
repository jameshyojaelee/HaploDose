"""Per-gene RNA dosage response from bulk tumour RNA with matched DNA copy number.

Bulk RNA is a mixture of tumour and normal cells. With tumour fraction p, a gene whose tumour-cell
expression scales as (CN / ploidy)^b gives, to first order in log space,
    log2 E_bulk ~= a + p * b * log2(CN / ploidy)
so the regressor is the purity-scaled relative copy number x = p * log2(TCN / ploidy), and its slope b
is the dosage response in tumour cells (b = 1 full dosage, b = 0 fully compensated). Purity also
enters as a main effect (normal-cell expression differs from tumour expression).

Per gene: log2(CPM + 0.5) ~ 1 + purity + covariates + x, fitted by OLS and by Huber M-estimation
(statsmodels RLM, HuberT t=1.345, MAD scale). The Huber slope is the primary slope; it is shrunk by
normal-normal empirical Bayes with prior mean and variance estimated within mean-expression deciles.

Cohorts:
  U-CAN (training): counts from E-MTAB-12862 processed matrix, FACETS segments (purity, ploidy from VCF header).
  TCGA COAD/READ (replication): STAR counts, ASCAT2 allele-specific segments (purity-corrected integer CN).
    No purity is available locally for TCGA, so purity is predicted from expression by a ridge model
    trained on U-CAN FACETS purity (training patients only). Sex is inferred from Y-gene expression.
The 20 Visium HD validation patients are removed before anything is fitted.

Usage: dosage_response.py OUTDIR [--threads N] [--seed S]
"""
import argparse
import gzip
import json
import os
import re
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

UCAN = Path('/gpfs/commons/groups/innovation/Bulk_WGS_EGAD50000000169/processed')
UCAN_COUNTS = UCAN / 'CRC.SW.mRNA.count.txt.gz'
UCAN_META = UCAN / 'Supplementary_Table_01.tsv'
VCF_DIR = Path('/gpfs/commons/home/jameslee/CRC/data/external/eva_cnv/vcf')
GTF = Path('/gpfs/commons/home/jameslee/CRC/data/external/ensembl112/Homo_sapiens.GRCh38.112.gtf.gz')
TCGA = Path('/gpfs/commons/home/jameslee/CRC/data/external/tcga_coadread')
HOLDOUT = ('U0007 U0046 U0178 U0179 U0209 U0245 U0247 U0285 U0354 U0383 U0396 U0438 U0440 U0469 U0991 '
           'U2761 U2955 U2964 U2972 U2975').split()
AUTOSOMES = [str(i) for i in range(1, 23)]
PSEUDO_CPM = 0.5        # log2(CPM + 0.5)
MIN_MEDIAN_CPM = 0.5    # gene filter, U-CAN training patients
TCN_FLOOR = 0.25        # as benchmarks/ucan_spatial/dna_truth.py
REL_CLIP = 2.0          # |log2(TCN / ploidy)| capped at 2 to limit leverage of amplicons and homozygous deletions
Y_GENES = ['RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY']


def log(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


# ---------------------------------------------------------------- annotation
def read_gtf_genes():
    rows = []
    pat = {k: re.compile(k + r' "([^"]+)"') for k in ('gene_id', 'gene_name', 'gene_biotype')}
    with gzip.open(GTF, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 8)
            if f[2] != 'gene':
                continue
            g = {k: (m.group(1) if (m := p.search(f[8])) else '') for k, p in pat.items()}
            rows.append((g['gene_id'], g['gene_name'], f[0], int(f[3]), int(f[4]), g['gene_biotype']))
    return pd.DataFrame(rows, columns=['ensembl_id', 'gtf_symbol', 'chromosome', 'start', 'end', 'biotype']
                        ).set_index('ensembl_id')


# ---------------------------------------------------------------- copy number
def read_facets(patient):
    vcfs = list(VCF_DIR.glob(f'*CRC-SW-{patient}-T*.vcf.gz'))
    if len(vcfs) != 1:
        return None
    head, rows = {}, []
    with gzip.open(vcfs[0], 'rt') as fh:
        for line in fh:
            if line.startswith('##purity=') or line.startswith('##ploidy='):
                k, v = line[2:].strip().split('=')
                head[k] = float(v) if v not in ('NA', '.', '') else np.nan
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            info = dict(x.split('=', 1) for x in f[7].split(';') if '=' in x)
            tcn = info.get('TCN_EM', '.')
            rows.append((f[0].removeprefix('chr'), int(f[1]), int(info['END']), float(tcn) if tcn != '.' else np.nan))
    seg = pd.DataFrame(rows, columns=['chromosome', 'start', 'end', 'tcn'])
    return seg, head.get('purity', np.nan), head.get('ploidy', np.nan)


def gene_rel_cn(seg, ploidy, genes):
    """log2(TCN / ploidy) at each gene midpoint (NaN where no segment covers it)."""
    out = np.full(len(genes), np.nan, dtype=np.float32)
    mid = ((genes.start.to_numpy() + genes.end.to_numpy()) // 2)
    chrom = genes.chromosome.to_numpy()
    for c, s in seg.groupby('chromosome'):
        s = s.sort_values('start')
        gi = np.flatnonzero(chrom == c)
        if not len(gi):
            continue
        k = np.searchsorted(s.start.to_numpy(), mid[gi], side='right') - 1
        ok = k >= 0
        ok[ok] &= mid[gi][ok] <= s.end.to_numpy()[k[ok]]
        tcn = s.tcn.to_numpy()[np.where(ok, k, 0)]
        out[gi[ok]] = tcn[ok]
    rel = np.log2(np.maximum(out, TCN_FLOOR) / ploidy)
    return np.clip(rel, -REL_CLIP, REL_CLIP).astype(np.float32)


# ---------------------------------------------------------------- cohorts
def load_ucan(genes_all):
    log('U-CAN counts')
    cnt = pd.read_csv(UCAN_COUNTS, sep='\t', index_col=0)
    cnt = cnt[[c for c in cnt.columns if c.endswith('.T')]]
    cnt.columns = [c.split('.')[2] for c in cnt.columns]
    n_all = cnt.shape[1]
    held = [p for p in HOLDOUT if p in cnt.columns]
    cnt = cnt.drop(columns=held)
    meta = pd.read_csv(UCAN_META, sep='\t').set_index('Sample ID')
    log(f'U-CAN tumour RNA samples {n_all}; holdout present and removed {len(held)}/{len(HOLDOUT)}')
    pur, plo, segs = {}, {}, {}
    for p in cnt.columns:
        r = read_facets(p)
        if r is None:
            continue
        segs[p], pur[p], plo[p] = r
    keep = [p for p in cnt.columns if p in segs and np.isfinite(pur[p]) and np.isfinite(plo[p]) and p in meta.index]
    log(f'U-CAN with FACETS purity and ploidy: {len(keep)}')
    cnt = cnt[keep]
    m = meta.loc[keep]
    cov = pd.DataFrame(index=keep)
    cov['purity'] = [pur[p] for p in keep]
    cov['ploidy'] = [plo[p] for p in keep]
    cov['male'] = m['Sex'].eq('Male').astype(float).to_numpy()
    cov['msi'] = m['MSI Status'].eq('MSI').astype(float).to_numpy()
    cov['biopsy'] = m['Tumour Specimen Type'].eq('Biopsy').astype(float).to_numpy()
    cov['pretreated'] = m['Pre-Treated'].eq('Treated').astype(float).to_numpy()
    return cnt, cov, segs


def load_tcga(universe):
    log('TCGA expression and ASCAT segments')
    fmap = pd.read_csv(TCGA / 'manifests/file_case_map.tsv', sep='\t')
    fmap = fmap[fmap.sample_type.eq('Primary Tumor')]
    fmap = fmap.assign(vial=fmap.sample_id.str[-3:])
    ex = fmap[fmap.data_type.eq('Gene Expression Quantification')].sort_values(['patient', 'vial', 'file_id'])
    sg = fmap[fmap.data_type.eq('Allele-specific Copy Number Segment')].sort_values(['patient', 'vial', 'file_id'])
    ex = ex.drop_duplicates('patient')
    sg = sg.drop_duplicates('patient')
    pts = sorted(set(ex.patient) & set(sg.patient))
    ex = ex.set_index('patient').loc[pts]
    sg = sg.set_index('patient').loc[pts]
    cols, keep = [], []
    for p in pts:
        d = TCGA / 'raw/gene_expression_quantification' / ex.loc[p, 'file_id']
        fs = list(d.glob('*.rna_seq.augmented_star_gene_counts.tsv'))
        if len(fs) != 1:
            continue
        t = pd.read_csv(fs[0], sep='\t', comment='#', usecols=['gene_id', 'unstranded'])
        t = t[t.gene_id.str.startswith('ENSG')]
        t.index = t.gene_id.str.split('.').str[0]
        t = t[~t.index.duplicated()]   # PAR_Y copies share the stable ID
        cols.append(t.unstranded.reindex(universe).astype(np.float32).rename(p))
        keep.append(p)
    cnt = pd.concat(cols, axis=1)
    segs, plo = {}, {}
    for p in keep:
        d = TCGA / 'raw/allele-specific_copy_number_segment' / sg.loc[p, 'file_id']
        fs = list(d.glob('*.seg.txt'))
        if len(fs) != 1:
            continue
        s = pd.read_csv(fs[0], sep='\t')
        s = pd.DataFrame(dict(chromosome=s.Chromosome.astype(str).str.removeprefix('chr'), start=s.Start,
                              end=s.End, tcn=s.Copy_Number.astype(float)))
        a = s[s.chromosome.isin(AUTOSOMES)]
        w = (a.end - a.start).to_numpy()
        plo[p] = float(np.average(a.tcn, weights=w))
        segs[p] = s
    keep = [p for p in keep if p in segs]
    cnt = cnt[keep]
    cov = pd.DataFrame(index=keep)
    cov['ploidy'] = [plo[p] for p in keep]
    cov['read'] = (ex.loc[keep, 'project'] == 'TCGA-READ').astype(float).to_numpy()
    log(f'TCGA patients with primary-tumour RNA and ASCAT segments: {len(keep)}')
    return cnt, cov, segs


def log_cpm(cnt):
    lib = cnt.sum(axis=0).to_numpy(dtype=np.float64)
    return np.log2(cnt.to_numpy(dtype=np.float64) / lib * 1e6 + PSEUDO_CPM).astype(np.float32), lib


# ---------------------------------------------------------------- per-gene fits
_G = {}


def _fit_chunk(idx):
    import statsmodels.api as sm
    Y, R, Z, pur = _G['Y'], _G['R'], _G['Z'], _G['S']
    out = np.full((len(idx), 9), np.nan)
    for j, g in enumerate(idx):
        y, r = Y[g].astype(np.float64), R[g].astype(np.float64)
        ok = np.isfinite(r)
        n = int(ok.sum())
        x = pur * r
        if n < Z.shape[1] + 10 or np.std(x[ok]) < 1e-6:
            out[j, 0] = n
            continue
        X = np.column_stack([Z[ok], x[ok]])
        yy = y[ok]
        XtX_inv = np.linalg.pinv(X.T @ X)
        beta = XtX_inv @ X.T @ yy
        res = yy - X @ beta
        df = n - X.shape[1]
        s2 = res @ res / df
        se = np.sqrt(s2 * XtX_inv[-1, -1])
        t = beta[-1] / se
        try:
            rl = sm.RLM(yy, X, M=sm.robust.norms.HuberT(t=1.345)).fit(maxiter=50, tol=1e-6)
            hb, hse = rl.params[-1], rl.bse[-1]
        except Exception:
            hb = hse = np.nan
        out[j] = (n, beta[-1], se, t, t * t / (t * t + df), stats.t.sf(abs(t), df) * 2, hb, hse,
                  np.std(r[ok]))
    return idx, out


def fit_cohort(Y, R, Z, threads, scale=None):
    """Regressor is scale * rel_cn; scale defaults to purity (column 1 of Z)."""
    _G.update(Y=Y, R=R, Z=Z, S=Z[:, 1] if scale is None else scale)
    chunks = np.array_split(np.arange(Y.shape[0]), max(1, threads * 8))
    res = np.full((Y.shape[0], 9), np.nan)
    with Pool(threads) as pool:
        for idx, out in pool.imap_unordered(_fit_chunk, chunks):
            res[idx] = out
    cols = ['n', 'ols_slope', 'ols_se', 'ols_t', 'r2', 'ols_p', 'slope', 'se', 'rel_cn_sd']
    return pd.DataFrame(res, columns=cols)


def eb_shrink(slope, se, mean_expr, n_bins=10):
    """Normal-normal EB per mean-expression decile: prior N(m, tau2), m and tau2 by method of moments."""
    shr = np.full(len(slope), np.nan)
    prior = []
    ok = np.isfinite(slope) & np.isfinite(se) & (se > 0)
    bins = pd.qcut(mean_expr[ok], n_bins, labels=False, duplicates='drop')
    idx = np.flatnonzero(ok)
    for b in np.unique(bins):
        k = idx[bins == b]
        w = 1 / se[k] ** 2
        m = np.sum(w * slope[k]) / np.sum(w)
        tau2 = max(np.var(slope[k]) - np.mean(se[k] ** 2), 1e-4)
        shr[k] = (slope[k] / se[k] ** 2 + m / tau2) / (1 / se[k] ** 2 + 1 / tau2)
        prior.append(dict(expr_bin=int(b), n_genes=len(k), mean_expr_lo=float(mean_expr[k].min()),
                          mean_expr_hi=float(mean_expr[k].max()), prior_mean=float(m), prior_sd=float(np.sqrt(tau2)),
                          median_se=float(np.median(se[k]))))
    return shr, pd.DataFrame(prior)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    o = np.argsort(pv)
    r = pv[o] * len(pv) / np.arange(1, len(pv) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty_like(pv)
    qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q


def cor(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    return dict(n=int(ok.sum()), pearson=float(np.corrcoef(a[ok], b[ok])[0, 1]),
                spearman=float(stats.spearmanr(a[ok], b[ok]).statistic))


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--threads', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', 4)))
    ap.add_argument('--seed', type=int, default=20261101)
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    summary = dict(seed=a.seed, holdout=HOLDOUT, pseudo_cpm=PSEUDO_CPM, min_median_cpm=MIN_MEDIAN_CPM,
                   tcn_floor=TCN_FLOOR, rel_clip=REL_CLIP)

    gtf = read_gtf_genes()
    log(f'GTF genes {len(gtf)}')

    # ---- U-CAN
    ucnt, ucov, usegs = load_ucan(gtf)
    # Gene universe for CPM: U-CAN matrix genes with a GTF record (TCGA STAR covers all Ensembl genes).
    universe = ucnt.index[ucnt.index.isin(gtf.index)]
    ucnt = ucnt.loc[universe]
    uY, ulib = log_cpm(ucnt)
    genes = gtf.loc[universe].copy()
    cpm = 2 ** uY - PSEUDO_CPM
    med = np.median(cpm, axis=1)
    fit_mask = genes.chromosome.isin(AUTOSOMES).to_numpy() & (med >= MIN_MEDIAN_CPM)
    fg = genes[fit_mask].copy()
    log(f'universe {len(universe)}; autosomal genes with median CPM >= {MIN_MEDIAN_CPM}: {len(fg)}')
    uR = np.vstack([gene_rel_cn(usegs[p], ucov.loc[p, 'ploidy'], fg) for p in ucnt.columns]).T
    uYf = uY[fit_mask]
    # Trans control: each gene paired with the gene half the fitted genome away (other chromosome).
    G = len(fg)
    shift = (np.arange(G) + G // 2) % G
    same_chr = fg.chromosome.to_numpy() == fg.chromosome.to_numpy()[shift]
    # CIN score: mean |log2(TCN/ploidy)| over fitted genes; adjusts for expression programmes that differ
    # between chromosomally unstable and stable tumours (confounders shared by all altered arms).
    ucov['cin'] = np.nanmean(np.abs(uR), axis=0)
    Zu_base = np.column_stack([np.ones(len(ucov)), ucov.purity, ucov.male, ucov.msi, ucov.biopsy, ucov.pretreated])
    Zu = np.column_stack([Zu_base, ucov.purity ** 2, np.log2(ucov.ploidy), ucov.cin])
    summary['ucan'] = dict(n_patients=len(ucov), purity_median=float(ucov.purity.median()),
                           purity_iqr=[float(x) for x in ucov.purity.quantile([.25, .75])],
                           ploidy_median=float(ucov.ploidy.median()), n_male=int(ucov.male.sum()),
                           n_msi=int(ucov.msi.sum()), n_biopsy=int(ucov.biopsy.sum()),
                           n_pretreated=int(ucov.pretreated.sum()), n_genes_fitted=int(len(fg)),
                           covariates=['purity', 'purity^2', 'log2 ploidy', 'CIN (mean |rel_cn|)', 'male', 'msi',
                                       'biopsy', 'pretreated'],
                           base_covariates=['purity', 'male', 'msi', 'biopsy', 'pretreated'])
    ucov.to_csv(out / 'ucan_patients.tsv', sep='\t')
    log('fitting U-CAN')
    fu = fit_cohort(uYf, uR, Zu, a.threads)
    fu_base = fit_cohort(uYf, uR, Zu_base, a.threads)
    uRt = uR[shift].copy()
    uRt[same_chr] = np.nan
    fu_trans = fit_cohort(uYf, uRt, Zu, a.threads)
    fu['mean_expr'] = uYf.mean(axis=1)
    fu['shrunk_slope'], prior_u = eb_shrink(fu.slope.to_numpy(), fu.se.to_numpy(), fu.mean_expr.to_numpy())
    fu['ols_q_bh'] = bh(fu.ols_p)
    prior_u.to_csv(out / 'eb_prior_ucan.tsv', sep='\t', index=False)

    # Within-cohort split-half reliability (random halves of U-CAN patients, OLS slope).
    perm = rng.permutation(len(ucov))
    h1, h2 = perm[: len(perm) // 2], perm[len(perm) // 2:]
    fa = fit_cohort(uYf[:, h1], uR[:, h1], Zu[h1], a.threads)
    fb = fit_cohort(uYf[:, h2], uR[:, h2], Zu[h2], a.threads)
    summary['ucan_split_half'] = dict(huber=cor(fa.slope.to_numpy(), fb.slope.to_numpy()),
                                      ols=cor(fa.ols_slope.to_numpy(), fb.ols_slope.to_numpy()))
    fu['split_half_a'] = fa.slope.to_numpy()
    fu['split_half_b'] = fb.slope.to_numpy()

    # Sensitivity: unscaled CN (x = rel_cn; purity main effect and covariates kept).
    fu_unscaled = fit_cohort(uYf, uR, Zu, a.threads, scale=np.ones(len(ucov)))

    # ---- purity predictor for TCGA (ridge on U-CAN training patients)
    log('TCGA load')
    tcnt, tcov, tsegs = load_tcga(universe)
    tY, tlib = log_cpm(tcnt.fillna(0))
    from sklearn.linear_model import RidgeCV
    from sklearn.model_selection import cross_val_predict
    hv = np.argsort(-uY.var(axis=1) * (med >= 1))[:3000]           # 3000 most variable expressed genes
    mu_u, sd_u = uY[hv].mean(1, keepdims=True), uY[hv].std(1, keepdims=True) + 1e-6
    mu_t, sd_t = tY[hv].mean(1, keepdims=True), tY[hv].std(1, keepdims=True) + 1e-6
    Xu = ((uY[hv] - mu_u) / sd_u).T
    Xt = ((tY[hv] - mu_t) / sd_t).T
    ridge = RidgeCV(alphas=np.logspace(0, 5, 21))
    cvp = cross_val_predict(ridge, Xu, ucov.purity.to_numpy(), cv=5)
    ridge.fit(Xu, ucov.purity.to_numpy())
    tcov['purity'] = np.clip(ridge.predict(Xt), 0.05, 1.0)
    summary['tcga_purity_model'] = dict(
        genes=3000, standardisation='per gene z-score within cohort', alpha=float(ridge.alpha_),
        ucan_cv5_pearson=float(np.corrcoef(cvp, ucov.purity)[0, 1]),
        ucan_cv5_mae=float(np.mean(np.abs(cvp - ucov.purity))),
        tcga_pred_median=float(tcov.purity.median()))
    # Sex from Y-chromosome genes.
    sym = gtf.loc[universe, 'gtf_symbol'].to_numpy()
    yi = np.flatnonzero(np.isin(sym, Y_GENES))
    yscore = tY[yi].mean(axis=0)
    tcov['male'] = (yscore > 2.0).astype(float)
    uy = uY[yi].mean(axis=0)
    summary['sex_from_y_check_ucan'] = dict(
        agreement=float(np.mean((uy > 2.0) == ucov.male.astype(bool))), y_genes=list(sym[yi]))
    tR = np.vstack([gene_rel_cn(tsegs[p], tcov.loc[p, 'ploidy'], fg) for p in tcnt.columns]).T
    tYf = tY[fit_mask]
    tcov['cin'] = np.nanmean(np.abs(tR), axis=0)
    Zt_base = np.column_stack([np.ones(len(tcov)), tcov.purity, tcov.male, tcov.read])
    Zt = np.column_stack([Zt_base, tcov.purity ** 2, np.log2(tcov.ploidy), tcov.cin])
    summary['tcga'] = dict(n_patients=len(tcov), n_male_inferred=int(tcov.male.sum()),
                           n_read=int(tcov.read.sum()), ploidy_median=float(tcov.ploidy.median()),
                           covariates=['predicted purity', 'predicted purity^2', 'log2 ploidy', 'CIN',
                                       'male (Y genes)', 'READ vs COAD'],
                           base_covariates=['predicted purity', 'male (Y genes)', 'READ vs COAD'])
    tcov.to_csv(out / 'tcga_patients.tsv', sep='\t')
    log('fitting TCGA')
    ft = fit_cohort(tYf, tR, Zt, a.threads)
    ft_base = fit_cohort(tYf, tR, Zt_base, a.threads)
    tRt = tR[shift].copy()
    tRt[same_chr] = np.nan
    ft_trans = fit_cohort(tYf, tRt, Zt, a.threads)
    ft['mean_expr'] = tYf.mean(axis=1)
    ft['shrunk_slope'], prior_t = eb_shrink(ft.slope.to_numpy(), ft.se.to_numpy(), ft.mean_expr.to_numpy())
    prior_t.to_csv(out / 'eb_prior_tcga.tsv', sep='\t', index=False)
    # TCGA sensitivity: no purity at all (x = rel_cn, purity column dropped).
    ft0 = fit_cohort(tYf, tR, np.column_stack([np.ones(len(tcov)), tcov.male, tcov.read]), a.threads,
                     scale=np.ones(len(tcov)))

    # ---- table
    tab = pd.DataFrame(dict(
        gene=fg.gtf_symbol.where(fg.gtf_symbol.ne(''), fg.index.to_series()).to_numpy(),
        ensembl_id=fg.index.to_numpy(), chromosome=fg.chromosome.to_numpy(), start=fg.start.to_numpy(),
        end=fg.end.to_numpy(), biotype=fg.biotype.to_numpy(),
        slope=fu.slope, se=fu.se, shrunk_slope=fu.shrunk_slope, r2=fu.r2, n=fu.n.astype('Int64'),
        mean_expr=fu.mean_expr, ols_slope=fu.ols_slope, ols_se=fu.ols_se, ols_p=fu.ols_p, ols_q_bh=fu.ols_q_bh,
        rel_cn_sd=fu.rel_cn_sd, split_half_a=fu.split_half_a, split_half_b=fu.split_half_b,
        tcga_slope=ft.slope, tcga_se=ft.se, tcga_shrunk_slope=ft.shrunk_slope, tcga_r2=ft.r2,
        tcga_n=ft.n.astype('Int64'), tcga_mean_expr=ft.mean_expr, tcga_ols_slope=ft.ols_slope,
        slope_base=fu_base.slope, tcga_slope_base=ft_base.slope,
        trans_slope=fu_trans.slope, tcga_trans_slope=ft_trans.slope,
        slope_unscaled=fu_unscaled.slope, tcga_slope_nopurity=ft0.slope))
    tab = (tab.assign(_c=tab.chromosome.astype(int)).sort_values(['_c', 'start']).drop(columns='_c')
           .reset_index(drop=True))
    num = tab.select_dtypes('float').columns
    tab[num] = tab[num].round(5)
    tab.to_csv(out / 'gene_dosage.tsv', sep='\t', index=False)

    # ---- replication and distribution
    s, t = tab.slope.to_numpy(float), tab.tcga_slope.to_numpy(float)
    prec = (tab.se < 0.15).to_numpy() & (tab.tcga_se < 0.15).to_numpy()
    summary['replication'] = dict(
        huber_all=cor(s, t), shrunk_all=cor(tab.shrunk_slope.to_numpy(float), tab.tcga_shrunk_slope.to_numpy(float)),
        ols_all=cor(tab.ols_slope.to_numpy(float), tab.tcga_ols_slope.to_numpy(float)),
        huber_both_se_lt_0p15=cor(s[prec], t[prec]),
        ucan_vs_tcga_nopurity=cor(s, tab.tcga_slope_nopurity.to_numpy(float)),
        ucan_scaled_vs_unscaled=cor(s, tab.slope_unscaled.to_numpy(float)),
        ucan_unscaled_vs_tcga_nopurity=cor(tab.slope_unscaled.to_numpy(float), tab.tcga_slope_nopurity.to_numpy(float)),
        # attenuation-corrected correlation: r_xy / sqrt(rel_x rel_y), rel = 1 - mean(se^2)/var(slope)
        reliability_ucan=float(1 - np.nanmean(tab.se ** 2) / np.nanvar(s)),
        reliability_tcga=float(1 - np.nanmean(tab.tcga_se ** 2) / np.nanvar(t)),
        slope_median_ucan=float(np.nanmedian(s)), slope_median_tcga=float(np.nanmedian(t)))
    rep = summary['replication']
    rep['huber_all_disattenuated'] = rep['huber_all']['pearson'] / np.sqrt(rep['reliability_ucan'] * rep['reliability_tcga'])
    ts, tt = tab.trans_slope.to_numpy(float), tab.tcga_trans_slope.to_numpy(float)
    summary['confounding_checks'] = dict(
        primary_vs_base_ucan=cor(s, tab.slope_base.to_numpy(float)),
        base_replication=cor(tab.slope_base.to_numpy(float), tab.tcga_slope_base.to_numpy(float)),
        base_median_ucan=float(np.nanmedian(tab.slope_base)),
        trans_pairs=int(np.isfinite(ts).sum()), trans_same_chr_dropped=int(same_chr.sum()),
        trans_median_ucan=float(np.nanmedian(ts)), trans_median_abs_ucan=float(np.nanmedian(np.abs(ts))),
        cis_median_abs_ucan=float(np.nanmedian(np.abs(s))),
        trans_sd_ucan=float(np.nanstd(ts)), cis_sd_ucan=float(np.nanstd(s)),
        trans_replication=cor(ts, tt),
        chromosome_median_shrunk=tab.groupby('chromosome').shrunk_slope.median().round(3).to_dict())
    sh = tab.shrunk_slope.to_numpy(float)
    sh = sh[np.isfinite(sh)]
    summary['distribution_shrunk_ucan'] = dict(
        n=int(len(sh)), quantiles={str(q): float(np.quantile(sh, q)) for q in (.05, .1, .25, .5, .75, .9, .95)},
        frac_lt_0p25=float(np.mean(sh < 0.25)), frac_0p25_0p5=float(np.mean((sh >= 0.25) & (sh < 0.5))),
        frac_0p5_0p75=float(np.mean((sh >= 0.5) & (sh < 0.75))), frac_ge_0p75=float(np.mean(sh >= 0.75)),
        frac_lt_0=float(np.mean(sh < 0)))
    summary['ols_tests'] = dict(test='two-sided t on OLS slope = 0', correction='Benjamini-Hochberg',
                                family='all fitted autosomal genes in U-CAN',
                                n_q_lt_0p05=int((tab.ols_q_bh < 0.05).sum()), n_tested=int(tab.ols_p.notna().sum()))
    with open(out / 'summary.json', 'w') as fh:
        json.dump(summary, fh, indent=2, default=float)
    log(json.dumps(summary, indent=2, default=float))


if __name__ == '__main__':
    main()
