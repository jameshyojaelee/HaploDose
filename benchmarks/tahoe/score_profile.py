"""Claim 8 scorer (docs/STATUS.md "Claim 8 prespecified 2026-11-06"): 1 Mb relative dosage profile vs DNA truth
per query cell line, pooled-line design (Kinker / MIX-seq / Tahoe input format). Written from
dev/profile_gap_20261106_0007/profile_gap.py (same readers, truth and 1 Mb rule) and checked on Kinker before use.

Per query line (single-line specimen, role query):
  truth     DepMap CN of the line minus the UMI-weighted mean of the pool's fit reference lines (evaluate.py rule)
  hd_sm5    HaploDose log_ratio_gene, clone values weighted by assigned query cells, then a running median over
            5 HaploDose bins within each chromosome (PRIMARY HaploDose export, fixed in advance)
  copykat   mean over the line's query cells of CopyKAT CNA_raw_results_gene_by_cell
  infercnvpy  infercnvpy window log_ratio per clone (wrapper layout), weighted by assigned query cells
  pb_sm5    pseudobulk per-gene log((y+0.5)/(E+0.5)), median-centred, median per HaploDose bin, same 5-bin running
            median (smoothed pseudobulk baseline)
  descriptive: hd_lrg (unsmoothed export), pb_bin10 (unsmoothed bin pseudobulk)
1 Mb rule: each method's genes averaged per 1 Mb bin; truth averaged over the same genes' bins; bins kept only
where every method in --methods-common has a value (the prespecified set: hd_sm5, copykat, infercnvpy, pb_sm5).
Endpoint: Pearson per line. Primary: hd_sm5 vs each comparator, paired Wilcoxon over lines (unit = cell line),
Holm over the comparators; claim holds if every Holm p < 0.05 and every mean difference > 0; bootstrap 95% CI of
the mean difference (10,000 resamples of lines, seed 20261106). A method missing on a line drops that line from that
pair only (its 1 Mb bins are then taken over the methods present). Power floor: < 20 primary lines -> not evaluable.
Usage: score_profile.py --input INP --hd R0 R1 --tool-root ROOT [--infercnvpy [--infercnvpy-root DIR]] [--lines LINES_TSV] OUT
  ROOT holds r{r}_pool*/copykat_*/native and r{r}_pool*/infercnvpy_* plus cells/r{r}_pool*.tsv (held-out layout).
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import pearsonr, wilcoxon
from statsmodels.stats.multitest import multipletests

np.random.seed(20261106)
BIN = 1_000_000
SEED, B, MIN_LINES = 20261106, 10_000, 20
PRIMARY = 'hd_sm5'


def load_counts(path):
    with h5py.File(path) as h:
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]))
        return X, h['genes'][:].astype(str), h['cell_ids'][:].astype(str)


def running_median(v, chrom, k=5):
    return pd.Series(v).groupby(chrom).transform(lambda x: x.rolling(k, center=True, min_periods=1).median()).to_numpy()


def query_cells(root, pool):
    meta = pd.read_csv(root / 'cells' / f'{pool}.tsv', sep='\t')
    return meta[meta.role.eq('query') & meta.design.eq('single_line')]


def read_copykat(root, rot):
    out = {}
    for d in sorted(root.glob(f'r{rot}_pool*/copykat_*')):
        f = d / 'native/run_copykat_raw_results_gene_by_cell.txt'
        if not (d / 'qc.json').exists() or not f.exists():
            continue
        q = query_cells(root, d.parent.name)
        head = pd.read_csv(f, sep='\t', nrows=0).columns
        keep = [c for c in head if c in set(q.cell_id)]
        m = pd.read_csv(f, sep='\t', usecols=['hgnc_symbol'] + keep, dtype={c: np.float32 for c in keep})
        m = m.drop_duplicates('hgnc_symbol').set_index('hgnc_symbol')
        for spec, c in q[q.cell_id.isin(keep)].groupby('specimen'):
            out[spec] = m[c.cell_id].mean(axis=1)
        print('copykat', d.parent.name, len(keep), 'cells', flush=True)
    return out


def clone_weighted(prof, cells, spec, col):
    """Per-bin value of a specimen: clone profiles weighted by its assigned cells per (MAP) clone."""
    a = cells[cells.specimen.eq(spec) & cells.status.eq('assigned')] if 'status' in cells else cells[cells.specimen.eq(spec)]
    n = a.clone.astype(str).value_counts()
    p = prof[prof.specimen.eq(spec) & prof.clone.astype(str).isin(n.index)]
    if p.empty:
        return None
    pv = p.pivot(index='bin', columns='clone', values=col)
    pv.columns = pv.columns.astype(str)
    nn = n[n.index.isin(pv.columns)]
    return (pv[nn.index] * nn.to_numpy()).sum(axis=1) / nn.sum()


def read_infercnvpy(root, rot):
    """{specimen: per-gene value} from the wrapper layout (genes -> window, clone-weighted window log_ratio)."""
    out = {}
    for d in sorted(root.glob(f'r{rot}_pool*/infercnvpy_*')):
        if not (d / 'qc.json').exists():
            continue
        cells = pd.read_csv(d / 'cells.tsv.gz', sep='\t', dtype={'clone': str}, keep_default_na=False)
        prof = pd.read_csv(d / 'profiles.tsv.gz', sep='\t', dtype={'clone': str})
        genes = pd.read_csv(d / 'genes.tsv.gz', sep='\t').drop_duplicates('gene').set_index('gene').bin
        q = query_cells(root, d.parent.name)
        for spec in q.specimen.unique():
            b = clone_weighted(prof, cells, spec, 'log_ratio')
            if b is not None:
                out[spec] = genes.map(b).dropna()
        print('infercnvpy', d.parent.name, flush=True)
    return out


def primary_tests(res, comps):
    rows, rng = [], np.random.default_rng(SEED)
    for m in comps:
        x = res[[f'mbP_{PRIMARY}', f'mbP_{m}']].dropna()
        d = (x.iloc[:, 0] - x.iloc[:, 1]).to_numpy()
        p = wilcoxon(d).pvalue if (d != 0).sum() >= 5 else np.nan
        boot = rng.choice(d, (B, len(d))).mean(axis=1) if len(d) else np.array([np.nan])
        rows.append(dict(comparator=m, n_lines=len(d), median_hd=float(x.iloc[:, 0].median()),
                         median_comparator=float(x.iloc[:, 1].median()), mean_diff=float(d.mean()) if len(d) else np.nan,
                         ci_low=float(np.quantile(boot, 0.025)), ci_high=float(np.quantile(boot, 0.975)),
                         hd_higher=int((d > 0).sum()), comparator_higher=int((d < 0).sum()), p=float(p)))
    r = pd.DataFrame(rows)
    r['p_holm'] = multipletests(r.p.fillna(1.0), method='holm')[1]
    r['superior'] = (r.p_holm < 0.05) & (r.mean_diff > 0)
    n = int(r.n_lines.min()) if len(r) else 0
    return dict(comparisons=r.to_dict('records'), n_min=n, evaluable=n >= MIN_LINES,
                claim_holds=bool(n >= MIN_LINES and r.superior.all()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', required=True, type=Path)
    ap.add_argument('--hd', nargs=2, required=True, type=Path)
    ap.add_argument('--tool-root', required=True, type=Path)
    ap.add_argument('--infercnvpy', action='store_true')
    ap.add_argument('--infercnvpy-root', type=Path, help='default: --tool-root')
    ap.add_argument('--lines', type=Path)
    ap.add_argument('out', type=Path)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    inp = a.input
    common_methods = ['hd_sm5', 'copykat', 'pb_sm5'] + (['infercnvpy'] if a.infercnvpy else [])
    methods = common_methods + ['hd_lrg', 'pb_bin10']
    comps = [m for m in ['copykat', 'infercnvpy', 'pb_sm5'] if m in common_methods]
    cn = pd.read_csv(inp / 'depmap_cn.tsv.gz', sep='\t', index_col=0)
    cn = cn[~cn.index.duplicated()]
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').drop_duplicates('gene').set_index('gene')
    X, genes, cells = load_counts(inp / 'counts.h5')
    row = pd.Series(np.arange(len(cells)), index=cells)
    rows = []
    for rot, run in enumerate(a.hd):
        ck = read_copykat(a.tool_root, rot)
        ic = read_infercnvpy(a.infercnvpy_root or a.tool_root, rot) if a.infercnvpy else {}
        meta = pd.read_csv(run / 'cells.tsv.gz', sep='\t', dtype={'clone': str}, keep_default_na=False)
        used = pd.read_csv(run / 'genes.tsv.gz', sep='\t')
        bins = pd.read_csv(run / 'bins.tsv', sep='\t')
        prof = pd.read_csv(run / 'profiles.tsv.gz', sep='\t', dtype={'clone': str})
        used = used[used.gene.isin(pos.index) & used.gene.isin(cn.index)].copy()
        used['chrom'] = 'chr' + pos.loc[used.gene, 'chromosome'].astype(str).str.removeprefix('chr').to_numpy()
        used['start'] = pos.loc[used.gene, 'start'].to_numpy()
        used['mb'] = used.chrom + ':' + (used.start // BIN).astype(str)
        gidx = pd.Index(genes).get_indexer(used.gene)
        bchrom = bins.chromosome.to_numpy()
        allbins = np.arange(len(bins))
        for spec, c in meta[meta.design.eq('single_line') & meta.role.eq('query')].groupby('specimen'):
            line, ds = c.truth_line.iloc[0], c.dataset.iloc[0]
            ref = meta[meta.dataset.eq(ds) & meta.role.eq('reference') & meta.ref_split.eq('fit')]
            w = ref.groupby('line').umi.sum()
            g = used[['gene', 'chrom', 'start', 'mb', 'bin']].copy()
            g['truth'] = cn.loc[used.gene, line].to_numpy() - cn.loc[used.gene, w.index].to_numpy() @ (w / w.sum()).to_numpy()
            y = np.asarray(X[row[c.cell_id].to_numpy()][:, gidx].sum(axis=0)).ravel()
            r = np.asarray(X[row[ref.cell_id].to_numpy()][:, gidx].sum(axis=0)).ravel()
            lr = np.log((y + 0.5) / (r / r.sum() * y.sum() + 0.5))
            g['pb'] = lr - np.median(lr)
            pbb = g.groupby('bin').pb.median().reindex(allbins)
            g['pb_bin10'] = g.bin.map(pbb).to_numpy()
            g['pb_sm5'] = g.bin.map(pd.Series(running_median(pbb.to_numpy(), bchrom), index=allbins)).to_numpy()
            b = clone_weighted(prof, meta, spec, 'log_ratio_gene')
            if b is not None:
                bv = b.reindex(allbins).to_numpy()
                g['hd_lrg'] = g.bin.map(b).to_numpy()
                g['hd_sm5'] = g.bin.map(pd.Series(running_median(bv, bchrom), index=allbins)).to_numpy()
            else:
                g['hd_lrg'] = g['hd_sm5'] = np.nan
            g['copykat'] = g.gene.map(ck[spec]).to_numpy() if spec in ck else np.nan
            g['infercnvpy'] = g.gene.map(ic[spec]).to_numpy() if spec in ic else np.nan
            g = g.dropna(subset=['truth'])
            present = [m for m in common_methods if g[m].notna().any()]
            d = dict(specimen=spec, line=line, rotation=rot, n_cells=len(c),
                     missing=','.join(m for m in common_methods if m not in present))
            mb = pd.concat({m: g.dropna(subset=[m]).groupby('mb')[m].mean() for m in methods if g[m].notna().any()}, axis=1)
            mb['truth'] = g.groupby('mb').truth.mean()
            mb = mb.dropna(subset=present + ['truth'])
            d['n_mb_bins'] = len(mb)
            for m in methods:
                ok = m in mb and mb[m].notna().sum() > 2 and mb[m].std() > 0
                d[f'mbP_{m}'] = pearsonr(mb[m], mb.truth)[0] if ok else np.nan
            rows.append(d)
    res = pd.DataFrame(rows)
    if a.lines is not None:
        lines = pd.read_csv(a.lines, sep='\t')
        prim = set(lines.loc[lines.primary.astype(bool), 'line'])
        res['primary'] = res.line.isin(prim)
    else:
        res['primary'] = True
    res.to_csv(a.out / 'line_scores.tsv', sep='\t', index=False)
    out = dict(primary=primary_tests(res[res.primary], comps),
               all_query_lines=primary_tests(res, comps),
               medians={m: float(res.loc[res.primary, f'mbP_{m}'].median()) for m in methods},
               n_query_lines=len(res), n_primary=int(res.primary.sum()),
               family=f'Holm over {len(comps)}: {PRIMARY} vs ' + ', '.join(comps) + ' (1 Mb Pearson per line)',
               unit='query cell line', descriptive='all_query_lines, hd_lrg and pb_bin10 are descriptive')
    (a.out / 'primary.json').write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps(out, indent=1, default=float))
    print('CLAIM8', out['primary']['claim_holds'], 'evaluable', out['primary']['evaluable'])


if __name__ == '__main__':
    main()
