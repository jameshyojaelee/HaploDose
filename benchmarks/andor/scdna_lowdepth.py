#!/usr/bin/env python3
"""Low-depth clone truth from the cell-tagged subset of the SRA Andor 2020 scDNA runs.

Why: in the SRA alignment load, only a few thousand fragments per cell keep the cell read group
(ID like 50941-0001); the rest carry the 8-mer sample index as read group. The cell-tagged subset is
biased by genome position (chr1 about 30% under, chr22 about 30% over, in SNU-638 and SNU-668), the
same for every cell. So each cell is normalised to the line's pooled cell-tagged profile, which cancels
the shared bias and keeps only how a cell differs from the line average. Absolute line copy number
comes from the index read groups (bulk of the line), GC corrected, with a Ginkgo ploidy scale.

Modality: single-cell low-pass WGS read depth, as scdna_cn.py, at ~3k fragments per cell.
Steps (per line):
1. Same 500 kb bin filter as scdna_cn.py, then coarse bins of --coarse consecutive kept bins within an
   arm (a short tail joins the previous coarse bin).
2. Cells: cell read groups with >= --min-fragments fragments.
3. Per cell: share per coarse bin / pooled cell share, then divided by its median within 10 GC strata
   (per-cell GC bias), giving ratio to the line average.
4. Clustering on per-cell arm log2 ratios (pooled counts per arm, so ~40 deep features): PCA, Leiden
   (resolution 1.0). Clusters < --min-clone cells are unassigned.
5. Clone consensus from pooled counts of the clone's cells (deep), on coarse bins:
   CN = bulk ratio x clone relative ratio x line ploidy scale (scale from the bulk profile, as Ginkgo).
   Clusters merged greedily while no arm has >= 30% of its coarse bins differing by >= 1 copy.
Outputs in OUT/<line>/ as scdna_cn.py (cells.tsv, clone_bins.tsv on coarse bins, clone_arms.tsv,
clone_differences.tsv, qc.json, cn.npz with per-cell relative ratios) plus bulk_bins.tsv.
Absolute clone CN here rests on the bulk profile; the tested quantity is the between-clone difference.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from scdna_cn import SEED, gc_correct, load_counts  # noqa: E402


def coarse_groups(arm_id, size):
    g = np.empty(len(arm_id), dtype=int)
    k = 0
    for a in pd.unique(arm_id):
        idx = np.flatnonzero(arm_id == a)
        n = max(1, len(idx) // size)
        lab = np.minimum(np.arange(len(idx)) // size, n - 1)
        g[idx] = k + lab
        k += n
    return g


def ploidy_scale(x):
    grid = np.arange(1.5, 6.0001, 0.01)
    sos = [np.mean((x * s - np.round(x * s)) ** 2) for s in grid]
    return float(grid[int(np.argmin(sos))])


def arm_diffs(cn, arm_c, a, b):
    d = cn[a] - cn[b]
    out = []
    for arm in pd.unique(arm_c):
        m = arm_c == arm
        big = np.abs(d[m]) >= 1
        if big.mean() >= 0.3:
            out.append((arm, float(big.mean()), float(d[m][big].mean())))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bins-dir', default='data/external/andor2020_scdna/bins')
    ap.add_argument('--lines', nargs='+', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--coarse', type=int, default=10, help='500 kb bins per coarse bin (10 = 5 Mb)')
    ap.add_argument('--min-fragments', type=int, default=1000)
    ap.add_argument('--min-clone', type=int, default=10)
    args = ap.parse_args(argv)
    import anndata as ad
    import scanpy as sc
    np.random.seed(SEED)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    ann_all = pd.read_csv(Path(args.bins_dir).parent / 'bin_annotation_500kb.tsv', sep='\t')
    keep = (ann_all.n_frac <= 0.1) & (ann_all.blacklist <= 0.1) & (ann_all.mappability >= 0.5)
    ann = ann_all[keep].reset_index(drop=True)
    grp = coarse_groups(ann.arm.to_numpy(), args.coarse)
    K = grp.max() + 1
    onehot = np.zeros((len(ann), K), np.float32)
    onehot[np.arange(len(ann)), grp] = 1
    cbins = ann.groupby(grp).agg(chromosome=('chromosome', 'first'), start=('start', 'min'), end=('end', 'max'),
                                 arm=('arm', 'first'), gc=('gc', 'mean'), mappability=('mappability', 'mean'))
    arm_c = cbins.arm.to_numpy()
    arms = pd.unique(arm_c)
    arm_onehot = np.stack([(arm_c == a).astype(np.float32) for a in arms], axis=1)
    for line in args.lines:
        d = out / line
        d.mkdir(exist_ok=False)
        rg = pd.read_csv(Path(args.bins_dir) / f'{line}.rg.txt', sep='\t', header=None)
        tags = [dict(f.split(':', 1) for f in r[1:] if ':' in str(f)) for r in rg.itertuples(index=False)]
        rg_map = {t['ID']: t.get('SM', '') for t in tags}
        names, M = load_counts(Path(args.bins_dir) / f'{line}.tsv.gz', ann)
        is_cell = np.array(['-' in n for n in names])
        # Bulk line profile from the index read groups: GC + mappability corrected, ploidy scale.
        bulk = gc_correct(M[~is_cell].sum(axis=0, keepdims=True), ann.gc.to_numpy())[0] / np.maximum(ann.mappability.to_numpy(), 0.5)
        bulk_c = (bulk @ onehot) / onehot.sum(axis=0)
        bulk_c = bulk_c / np.mean(bulk_c)
        scale = ploidy_scale(bulk_c)
        C = M[is_cell] @ onehot
        cells = names[is_cell]
        frag = C.sum(axis=1)
        ok = frag >= args.min_fragments
        C, cells_ok = C[ok], cells[ok]
        pool = C.sum(axis=0)
        R = (C / C.sum(axis=1, keepdims=True)) / (pool / pool.sum())
        R = gc_correct(R, cbins.gc.to_numpy(), strata=10)
        # Per-cell arm features from pooled counts per arm, relative to the pool.
        A = (C @ arm_onehot) / C.sum(axis=1, keepdims=True) / ((pool @ arm_onehot) / pool.sum())
        X = np.log2(np.maximum(A, 0.05)).astype(np.float32)
        a = ad.AnnData(X=X - X.mean(axis=0))
        sc.pp.pca(a, n_comps=min(20, X.shape[0] - 1, X.shape[1] - 1), random_state=SEED)
        sc.pp.neighbors(a, n_neighbors=15, random_state=SEED)
        sc.tl.leiden(a, resolution=1.0, flavor='igraph', n_iterations=2, random_state=SEED)
        lab = a.obs.leiden.astype(int).to_numpy()

        def consensus(mask):
            p = C[mask].sum(axis=0)
            rel = (p / p.sum()) / (pool / pool.sum())
            return np.clip(bulk_c * rel * scale, 0, 12)

        while True:
            ks = [k for k in np.unique(lab) if (lab == k).sum() >= args.min_clone]
            cons = {k: consensus(lab == k) for k in ks}
            pairs = [(len(arm_diffs(cons, arm_c, i, j)), i, j) for x, i in enumerate(ks) for j in ks[x + 1:]]
            if not pairs:
                break
            n, i, j = min(pairs)
            if n >= 1:
                break
            lab[lab == j] = i
        sizes = pd.Series(lab).value_counts()
        big = sizes[sizes >= args.min_clone].index
        order = {k: n for n, k in enumerate(sizes[big].index, start=1)}
        clone = np.array([f'{line}_d{order[k]}' if k in order else 'unassigned' for k in lab], dtype=object)
        cons = {f'{line}_d{order[k]}': consensus(lab == k) for k in order}
        cn_names = list(cons)

        cells_tab = pd.DataFrame({'rg': cells, 'barcode': [rg_map.get(c, '') for c in cells], 'fragments': frag,
                                  'pass_depth': ok})
        cells_tab['pass_mapd'] = ok          # no MAPD filter at this depth; column kept for score.py
        cells_tab['clone'] = 'excluded'
        cells_tab.loc[np.flatnonzero(ok), 'clone'] = clone
        cells_tab.to_csv(d / 'cells.tsv', sep='\t', index=False)
        np.savez_compressed(d / 'cn.npz', rel_ratio=R.astype(np.float32), arm_log2=X, cells=cells_ok,
                            bins=cbins[['chromosome', 'start', 'end', 'arm']].to_numpy(dtype=str), arms=arms.astype(str))
        cb = pd.concat([cbins[['chromosome', 'start', 'end', 'arm']].reset_index(drop=True),
                        pd.DataFrame(cons)], axis=1)
        cb.to_csv(d / 'clone_bins.tsv', sep='\t', index=False)
        cb.groupby('arm', sort=False)[cn_names].mean().reset_index().to_csv(d / 'clone_arms.tsv', sep='\t', index=False)
        cbins.assign(bulk_ratio=bulk_c, bulk_cn=bulk_c * scale).to_csv(d / 'bulk_bins.tsv', sep='\t', index=False)
        drows = [dict(clone_a=x, clone_b=y, arm=arm, frac_bins_diff=f, mean_diff=m)
                 for i, x in enumerate(cn_names) for y in cn_names[i + 1:] for arm, f, m in arm_diffs(cons, arm_c, x, y)]
        pd.DataFrame(drows, columns=['clone_a', 'clone_b', 'arm', 'frac_bins_diff', 'mean_diff']).to_csv(
            d / 'clone_differences.tsv', sep='\t', index=False)
        qc = dict(line=line, truth='low_depth_cell_tagged_subset', n_read_groups=int(is_cell.sum()),
                  n_index_read_groups=int((~is_cell).sum()), index_fragments=float(M[~is_cell].sum()),
                  n_pass_depth=int(ok.sum()), n_pass_mapd=int(ok.sum()), median_fragments=float(np.median(frag)),
                  ploidy_scale=scale, n_coarse_bins=int(K), clones={c: int((clone == c).sum()) for c in cn_names},
                  n_unassigned=int((clone == 'unassigned').sum()), params=vars(args), seed=SEED)
        (d / 'qc.json').write_text(json.dumps(qc, indent=2) + '\n')
        print(json.dumps({k: qc[k] for k in ['line', 'n_pass_depth', 'median_fragments', 'ploidy_scale', 'clones']}), flush=True)


if __name__ == '__main__':
    main()
