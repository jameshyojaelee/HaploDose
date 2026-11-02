#!/usr/bin/env python3
"""infercnvpy 0.6.1 comparator (micromamba env scrna_cnv) on the Haplodose input contract.

Built for scRNA-seq UMI counts. On Visium HD (probe-based, segmented cells) it runs unchanged:
the reference-ratio model needs only counts and gene positions, but its smoothing window
assumes droplet-like per-gene noise [Inference].
Steps:
- normalize_total(1e4) + log1p on all cells; autosomal genes with mean count >= --min-gene-mean.
- infercnv(window 100, step 10) against fit reference cells, one category per ref_group.
- malignant calls and clones: cnv_post.call (see its docstring).
- Output in the Haplodose layout (cmp_common.write_layout); units are infercnvpy windows.
"""
from __future__ import annotations

import argparse
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmp_common as C  # noqa: E402
import cnv_post  # noqa: E402


def window_table(var, chr_pos, n_win, window, step):
    """Genomic span of each infercnvpy window and each gene's window, reproducing its gene order
    (per chromosome, sorted by start; window k covers genes k*step .. k*step+window-1)."""
    order = sorted(chr_pos.items(), key=lambda kv: kv[1])
    bins, genes = [], []
    for i, (chrom, lo) in enumerate(order):
        hi = order[i + 1][1] if i + 1 < len(order) else n_win
        g = var.loc[var.chromosome == chrom].sort_values('start')
        n = hi - lo
        for k in range(n):
            a, b = (k * step, min(k * step + window, len(g))) if len(g) > window else (0, len(g))
            blk = g.iloc[a:b]
            bins.append(dict(bin=lo + k, chromosome=chrom, start=int(blk.start.min()), end=int(blk.end.max()),
                             n_genes=len(blk)))
        centre = np.array([k * step + window / 2 for k in range(n)]) if len(g) > window else np.zeros(1)
        idx = np.abs(np.arange(len(g))[:, None] - centre[None, :]).argmin(axis=1)
        genes.append(pd.DataFrame(dict(gene=g.index, chromosome=chrom, start=g.start.to_numpy(),
                                       end=g.end.to_numpy(), bin=lo + idx)))
    return pd.DataFrame(bins), pd.concat(genes, ignore_index=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--cells', required=True)
    ap.add_argument('--positions', required=True)
    ap.add_argument('--out-parent', required=True, help='a new <out-parent>/infercnvpy_<timestamp> is created')
    ap.add_argument('--exclude-genes', help='one gene per line, never used')
    ap.add_argument('--query-roles', default='tumour_epithelium,query')
    ap.add_argument('--min-gene-mean', type=float, default=0.02)
    ap.add_argument('--window', type=int, default=100)
    ap.add_argument('--step', type=int, default=10)
    ap.add_argument('--rule', choices=['cluster', 'cell', 'all'], default='cluster')
    args = ap.parse_args(argv)

    import anndata as ad
    import infercnvpy as cnv
    import scanpy as sc

    np.random.seed(C.SEED)
    t0 = time.monotonic()
    out = C.new_run_dir(args.out_parent, 'infercnvpy')
    C.provenance(out, 'infercnvpy', args, infercnvpy=cnv.__version__, scanpy=sc.__version__)
    C.pip_freeze(out)

    cells = C.read_cells(args.cells)
    X, genes, _ = C.read_counts(args.counts, cells.cell_id)
    X, genes = C.collapse_genes(X, genes)
    pos = C.read_positions(args.positions).set_index('gene')
    keep = pd.Index(genes).isin(pos.index[pos.chromosome.isin(C.AUTOSOMES)])
    keep &= np.asarray(X.mean(axis=0)).ravel() >= args.min_gene_mean
    if args.exclude_genes:
        keep &= ~pd.Index(genes).isin(Path(args.exclude_genes).read_text().split())
    X, genes = X[:, keep], genes[keep]
    var = pos.loc[genes, ['chromosome', 'start', 'end']]
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=cells.cell_id.to_numpy()), var=var)
    del X

    fit = (cells.role.eq('reference') & cells.ref_split.eq('fit')).to_numpy()
    hold = (cells.role.eq('reference') & cells.ref_split.eq('holdout')).to_numpy()
    query = C.query_mask(cells, args.query_roles.split(','))
    if fit.sum() < 20 or hold.sum() < 20:
        raise ValueError('Need at least 20 fit and 20 held-out reference cells')
    groups = sorted(cells.ref_group[fit].unique())
    adata.obs['ref_key'] = np.where(fit, cells.ref_group, 'not_reference')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    t_load = time.monotonic()
    cnv.tl.infercnv(adata, reference_key='ref_key', reference_cat=groups, window_size=args.window,
                    step=args.step, n_jobs=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    t_infer = time.monotonic()
    Z = adata.obsm['X_cnv']
    Z = (Z.toarray() if hasattr(Z, 'toarray') else np.asarray(Z)).astype(np.float32)
    bins, genes_df = window_table(adata.var, adata.uns['cnv']['chr_pos'], Z.shape[1], args.window, args.step)
    del adata

    score, malignant, clone, clusters = cnv_post.call(Z, cells, query, hold, rule=args.rule, seed=C.SEED)
    profiles = C.clone_profiles(Z, cells, clone, hold)
    normal_like = [c for c in pd.unique(clone) if str(c).endswith('_normal')]
    C.write_layout(out, cells, clone, malignant, bins, genes_df, profiles, normal_like, hold)
    C.write_table(pd.DataFrame({'cell_id': cells.cell_id, 'cnv_score': score}), out / 'scores.tsv')
    C.write_table(clusters, out / 'call_clusters.tsv')
    C.write_json(dict(tool='infercnvpy', n_cells=len(cells), n_genes=int(keep.sum()), n_windows=Z.shape[1],
                      n_fit=int(fit.sum()), n_holdout=int(hold.sum()), rule=args.rule,
                      n_malignant_query=int((malignant & query).sum()),
                      holdout_malignant_rate=float(malignant[hold].mean()),
                      load_seconds=t_load - t0, inference_seconds=t_infer - t_load,
                      wall_seconds=time.monotonic() - t0,
                      max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2),
                 out / 'qc.json')
    print(out)


if __name__ == '__main__':
    main()
