#!/usr/bin/env python3
"""Per-cell copy number and clones from Andor 2020 10x single-cell CNV bin counts (fetch_bins.sbatch).

Modality: single-cell whole-genome DNA (low-pass, ~0.01-0.05x per cell). Counts per 500 kb bin are
read-depth; the model is the Ginkgo/HMMcopy family (GC-corrected depth ratio, piecewise-constant
segments, integer copy number by a ploidy scale), built for this modality.

Steps (per line):
1. Bins: 500 kb on chr1-22; drop bins with > 10% N, > 10% ENCODE blacklist v2 overlap, or < 50%
   100-mer unique mappability (bin_annotation.tsv, built once from the reference).
2. Cells: read groups (Cell Ranger DNA cell calls). Keep cells with >= --min-fragments fragments.
3. Per cell: GC correction (median count in 25 GC-quantile strata, linear interpolation), then
   mappability correction (divide by bin mappability), ratio to the cell mean.
4. Noise: MAPD = median |r_{i+1} - r_i| / median r over adjacent bins of a chromosome. Cells above
   --max-mapd are excluded (noisy and S-phase cells; S-phase replication timing raises MAPD).
5. Segmentation per cell: recursive binary segmentation within each chromosome arm, split kept if
   |t| > --t-split and both parts >= 4 bins.
6. Integer CN: one ploidy scale per line, s in [1.5, 6] minimising sum of squared distances of
   segment means x s to integers over the line's pooled segments (Ginkgo). CN = round(mean x s),
   capped at 12. A scale is relative: whole-genome doublings are not resolved.
7. Clones: Leiden (resolution 1.0) on PCA of segment ratios; clusters then merged greedily while
   any pair of consensus profiles (median CN) differs in < --min-diff-mb Mb, counting only runs of
   >= 6 adjacent bins (3 Mb) that differ by >= 1 copy. Clusters < --min-clone cells -> 'unassigned'.
Outputs in OUT/<line>/: cn.npz (cells x bins CN int8, segment ratio float32), cells.tsv,
clone_bins.tsv, clone_arms.tsv, clone_differences.tsv, qc.json.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261101
BIN = 500_000
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]
CEN = Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv'
FASTA = '/gpfs/commons/groups/vickovic_lab/refs/human/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa'
MAPPABILITY = '/gpfs/commons/home/jameslee/CRC/data/external/ase_artifacts/unique_mappability.100bp.nochr.bed'
BLACKLIST = '/gpfs/commons/home/jameslee/CRC/data/external/encode_blacklist/hg38-blacklist.v2.bed'


def overlap_fraction(bed, chrom_prefix=''):
    """Fraction of each 500 kb bin covered by BED intervals; dict chrom -> array."""
    frac = {}
    for chunk in pd.read_csv(bed, sep='\t', header=None, usecols=[0, 1, 2], names=['c', 's', 'e'],
                             dtype={'c': str}, chunksize=2_000_000):
        chunk['c'] = chrom_prefix + chunk.c
        for c, g in chunk[chunk.c.isin(AUTOSOMES)].groupby('c'):
            a = frac.setdefault(c, np.zeros(260_000_000 // BIN + 1))
            s, e = g.s.to_numpy(), g.e.to_numpy()
            b0, b1 = s // BIN, (e - 1) // BIN
            same = b0 == b1
            np.add.at(a, b0[same], (e - s)[same])
            for x, y, bs, be in zip(s[~same], e[~same], b0[~same], b1[~same]):
                for b in range(bs, be + 1):
                    a[b] += min(y, (b + 1) * BIN) - max(x, b * BIN)
    return {c: a / BIN for c, a in frac.items()}


def build_annotation(dest):
    """GC and N fraction per bin from the reference FASTA, mappability and blacklist overlap."""
    rows, chrom, seq = [], None, []

    def flush():
        if chrom in AUTOSOMES:
            s = ''.join(seq).upper()
            for b in range(0, len(s), BIN):
                w = s[b:b + BIN]
                n = w.count('N')
                gc = (w.count('G') + w.count('C')) / max(len(w) - n, 1)
                rows.append(dict(chromosome=chrom, bin=b // BIN, start=b, end=b + len(w), gc=gc, n_frac=n / len(w)))
    with open(FASTA) as fh:
        for line in fh:
            if line.startswith('>'):
                flush()
                chrom, seq = line[1:].split()[0], []
            elif chrom in AUTOSOMES:
                seq.append(line.strip())
        flush()
    ann = pd.DataFrame(rows)
    mp, bl = overlap_fraction(MAPPABILITY, 'chr'), overlap_fraction(BLACKLIST)
    ann['mappability'] = [mp.get(c, np.zeros(1))[b] if b < len(mp.get(c, [])) else 0 for c, b in zip(ann.chromosome, ann.bin)]
    ann['blacklist'] = [bl[c][b] if c in bl and b < len(bl[c]) else 0 for c, b in zip(ann.chromosome, ann.bin)]
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    ann['arm'] = ann.chromosome + np.where(ann.start >= ann.chromosome.map(cen), 'q', 'p')
    ann.to_csv(dest, sep='\t', index=False)
    return ann


def load_counts(path, ann):
    df = pd.read_csv(path, sep='\t', header=None, names=['rg', 'chromosome', 'bin', 'n'], dtype={'rg': str})
    df = df[df.chromosome.isin(AUTOSOMES)]
    key = pd.MultiIndex.from_frame(ann[['chromosome', 'bin']])
    col = key.get_indexer(pd.MultiIndex.from_frame(df[['chromosome', 'bin']]))
    df = df[col >= 0].assign(col=col[col >= 0])
    cells = np.sort(df.rg.unique())
    row = pd.Index(cells).get_indexer(df.rg)
    M = np.zeros((len(cells), len(ann)), dtype=np.float32)
    np.add.at(M, (row, df.col.to_numpy()), df.n.to_numpy())
    return cells, M


def gc_correct(M, gc, strata=25):
    q = np.quantile(gc, np.linspace(0, 1, strata + 1))
    idx = np.clip(np.searchsorted(q, gc, side='right') - 1, 0, strata - 1)
    centres = np.array([np.median(gc[idx == k]) for k in range(strata)])
    out = np.empty_like(M)
    for i in range(M.shape[0]):
        med = np.array([np.median(M[i, idx == k]) for k in range(strata)])
        med = np.maximum(med, 1e-3)
        out[i] = M[i] / np.interp(gc, centres, med)
    return out


def segment(x, t_split, min_len=4):
    """Recursive binary segmentation of x; returns segment-mean vector."""
    out = np.empty_like(x)
    stack = [(0, len(x))]
    while stack:
        a, b = stack.pop()
        y = x[a:b]
        n = len(y)
        best = None
        if n >= 2 * min_len:
            cs = np.cumsum(y)
            k = np.arange(min_len, n - min_len + 1)
            m1, m2 = cs[k - 1] / k, (cs[-1] - cs[k - 1]) / (n - k)
            sd = max(np.std(np.diff(y)) / np.sqrt(2), 1e-6)
            t = np.abs(m1 - m2) / (sd * np.sqrt(1 / k + 1 / (n - k)))
            j = int(np.argmax(t))
            if t[j] > t_split:
                best = a + int(k[j])
        if best is None:
            out[a:b] = y.mean()
        else:
            stack += [(a, best), (best, b)]
    return out


def clone_consensus(cn, lab, k):
    return np.median(cn[lab == k], axis=0)


def diff_mb(p, q, min_run=6):
    d = np.abs(p - q) >= 1
    tot, run = 0, 0
    for v in np.append(d, False):
        if v:
            run += 1
        else:
            if run >= min_run:
                tot += run
            run = 0
    return tot * BIN / 1e6


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bins-dir', default='data/external/andor2020_scdna/bins')
    ap.add_argument('--lines', nargs='+', required=True)
    ap.add_argument('--out', required=True, help='new results/andor_<stamp>/scdna directory')
    ap.add_argument('--min-fragments', type=int, default=50_000)
    ap.add_argument('--max-mapd', type=float, default=0.45)
    ap.add_argument('--t-split', type=float, default=5.0)
    ap.add_argument('--min-diff-mb', type=float, default=10.0)
    ap.add_argument('--min-clone', type=int, default=10)
    args = ap.parse_args(argv)
    import anndata as ad
    import scanpy as sc
    np.random.seed(SEED)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ann_path = Path(args.bins_dir).parent / 'bin_annotation_500kb.tsv'
    ann_all = pd.read_csv(ann_path, sep='\t') if ann_path.exists() else build_annotation(ann_path)
    keep_bin = (ann_all.n_frac <= 0.1) & (ann_all.blacklist <= 0.1) & (ann_all.mappability >= 0.5)
    ann = ann_all[keep_bin].reset_index(drop=True)
    arm_id = ann.arm.to_numpy()
    for line in args.lines:
        d = out / line
        d.mkdir(exist_ok=False)
        rg = pd.read_csv(Path(args.bins_dir) / f'{line}.rg.txt', sep='\t', header=None)
        tags = [dict(f.split(':', 1) for f in r[1:] if ':' in str(f)) for r in rg.itertuples(index=False)]
        rg_map = {t['ID']: t.get('SM', '') for t in tags}
        cells, M = load_counts(Path(args.bins_dir) / f'{line}.tsv.gz', ann)
        frag = M.sum(axis=1)
        ok = frag >= args.min_fragments
        G = gc_correct(M[ok], ann.gc.to_numpy()) / np.maximum(ann.mappability.to_numpy(), 0.5)
        R = G / G.mean(axis=1, keepdims=True)
        same = np.append(arm_id[1:] == arm_id[:-1], False)
        mapd = np.array([np.median(np.abs(np.diff(r)[same[:-1]])) / max(np.median(r), 1e-6) for r in R])
        good = mapd <= args.max_mapd
        S = np.vstack([np.concatenate([segment(r[arm_id == a], args.t_split) for a in pd.unique(arm_id)]) for r in R[good]])
        segs = np.unique(np.round(S, 4))
        grid = np.arange(1.5, 6.0001, 0.01)
        sos = [np.mean((np.unique(S.round(3)) * s - np.round(np.unique(S.round(3)) * s)) ** 2) for s in grid]
        scale = float(grid[int(np.argmin(sos))])
        CN = np.clip(np.round(S * scale), 0, 12).astype(np.int8)
        # Clones
        a = ad.AnnData(X=S.astype(np.float32))
        sc.pp.pca(a, n_comps=min(30, S.shape[0] - 1, S.shape[1] - 1), random_state=SEED)
        sc.pp.neighbors(a, n_neighbors=15, random_state=SEED)
        sc.tl.leiden(a, resolution=1.0, flavor='igraph', n_iterations=2, random_state=SEED)
        lab = a.obs.leiden.astype(int).to_numpy()
        while True:
            ks = [k for k in np.unique(lab) if (lab == k).sum() >= args.min_clone]
            cons = {k: clone_consensus(CN, lab, k) for k in ks}
            pairs = [(diff_mb(cons[i], cons[j]), i, j) for x, i in enumerate(ks) for j in ks[x + 1:]]
            if not pairs:
                break
            dmin, i, j = min(pairs)
            if dmin >= args.min_diff_mb:
                break
            lab[lab == j] = i
        sizes = pd.Series(lab).value_counts()
        big = sizes[sizes >= args.min_clone].index
        order = {k: n for n, k in enumerate(sizes[big].index, start=1)}
        clone = np.array([f'{line}_d{order[k]}' if k in order else 'unassigned' for k in lab], dtype=object)
        idx_ok = np.flatnonzero(ok)
        cells_tab = pd.DataFrame({'rg': cells, 'barcode': [rg_map.get(c, '') for c in cells], 'fragments': frag,
                                  'pass_depth': ok})
        cells_tab['mapd'] = np.nan
        cells_tab.loc[idx_ok, 'mapd'] = mapd
        cells_tab['pass_mapd'] = False
        cells_tab.loc[idx_ok[good], 'pass_mapd'] = True
        cells_tab['clone'] = 'excluded'
        cells_tab.loc[idx_ok[good], 'clone'] = clone
        cells_tab.to_csv(d / 'cells.tsv', sep='\t', index=False)
        np.savez_compressed(d / 'cn.npz', cn=CN, seg_ratio=S.astype(np.float32), cells=cells[idx_ok[good]],
                            bins=ann[['chromosome', 'start', 'end', 'arm']].to_numpy(dtype=str))
        names = [c for c in pd.unique(clone) if c != 'unassigned']
        cb = pd.DataFrame({c: np.median(CN[clone == c], axis=0) for c in names})
        cb = pd.concat([ann[['chromosome', 'start', 'end', 'arm']], cb], axis=1)
        cb.to_csv(d / 'clone_bins.tsv', sep='\t', index=False)
        ca = cb.groupby('arm', sort=False)[names].mean().reset_index()
        ca.to_csv(d / 'clone_arms.tsv', sep='\t', index=False)
        drows = []
        for x, ca_ in enumerate(names):
            for cb_ in names[x + 1:]:
                dd = (cb[ca_] - cb[cb_]).to_numpy()
                for arm, g in pd.Series(dd).groupby(arm_id, sort=False):
                    big_run = np.abs(g.to_numpy()) >= 1
                    frac = float(big_run.mean())
                    if frac >= 0.3:
                        drows.append(dict(clone_a=ca_, clone_b=cb_, arm=arm, frac_bins_diff=frac,
                                          mean_diff=float(g[big_run].mean())))
        pd.DataFrame(drows, columns=['clone_a', 'clone_b', 'arm', 'frac_bins_diff', 'mean_diff']).to_csv(
            d / 'clone_differences.tsv', sep='\t', index=False)
        qc = dict(line=line, n_read_groups=len(cells), n_pass_depth=int(ok.sum()), n_pass_mapd=int(good.sum()),
                  median_fragments=float(np.median(frag)), median_mapd=float(np.median(mapd)), ploidy_scale=scale,
                  n_bins=len(ann), clones={c: int((clone == c).sum()) for c in names},
                  n_unassigned=int((clone == 'unassigned').sum()), params=vars(args), seed=SEED)
        (d / 'qc.json').write_text(json.dumps(qc, indent=2) + '\n')
        print(json.dumps({k: qc[k] for k in ['line', 'n_pass_mapd', 'median_mapd', 'ploidy_scale', 'clones']}))


if __name__ == '__main__':
    main()
