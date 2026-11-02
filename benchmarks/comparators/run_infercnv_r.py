#!/usr/bin/env python3
"""R inferCNV comparator on the Haplodose input contract. Run with the scrna_cnv python;
the R step needs RSCRIPT and R_LD_LIBRARY_PATH from `module load R/4.4.1 JAGS/4.3.2-foss-2023a`
(module use /nfs/sw/modules), exported by the caller after it purges the modules (see smoke2.sbatch).

Built for scRNA-seq UMI counts (Smart-seq2 and 10x). On Visium HD it runs unchanged [Inference:
its gene-level residual model needs only counts and positions].
Cells sent to inferCNV: fit reference cells (one group per ref_group) and query plus held-out
reference cells as observations. inferCNV drops cells under 300 UMI; they get clone NA.
Post-processing: log2 of the final residual matrix, averaged over consecutive 50-gene bins, then
cnv_post.call (malignant calls and clones) and cmp_common.write_layout.
"""
from __future__ import annotations

import argparse
import gzip
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmp_common as C  # noqa: E402
import cnv_post  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--cells', required=True)
    ap.add_argument('--positions', required=True)
    ap.add_argument('--out-parent', required=True, help='a new <out-parent>/infercnv_r_<timestamp> is created')
    ap.add_argument('--exclude-genes', help='one gene per line, never used')
    ap.add_argument('--query-roles', default='tumour_epithelium,query')
    ap.add_argument('--bin-genes', type=int, default=50)
    ap.add_argument('--rule', choices=['cluster', 'cell', 'all'], default='cluster')
    ap.add_argument('--ncores', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', '8')))
    args = ap.parse_args(argv)

    np.random.seed(C.SEED)
    t0 = time.monotonic()
    out = C.new_run_dir(args.out_parent, 'infercnv_r')
    C.provenance(out, 'infercnv_r', args)
    C.pip_freeze(out)
    prep = C.guard(out / 'prep')
    prep.mkdir()

    cells = C.read_cells(args.cells)
    X, genes, _ = C.read_counts(args.counts, cells.cell_id)
    X, genes = C.collapse_genes(X, genes)
    pos = C.read_positions(args.positions)
    auto = pos[pos.chromosome.isin(C.AUTOSOMES)].set_index('gene')
    keep = pd.Index(genes).isin(auto.index)
    if args.exclude_genes:
        keep &= ~pd.Index(genes).isin(Path(args.exclude_genes).read_text().split())
    X, genes = X[:, keep], genes[keep]
    fit = (cells.role.eq('reference') & cells.ref_split.eq('fit')).to_numpy()
    hold = (cells.role.eq('reference') & cells.ref_split.eq('holdout')).to_numpy()
    query = C.query_mask(cells, args.query_roles.split(','))
    send = np.flatnonzero(fit | query | hold)
    sent = cells.iloc[send]
    pd.DataFrame({'gene': genes}).to_csv(prep / 'genes.tsv', sep='\t', index=False)
    sent[['cell_id']].to_csv(prep / 'cells.tsv', sep='\t', index=False)
    group = np.where(fit[send], 'ref_' + sent.ref_group, 'observation')
    pd.DataFrame({'c': sent.cell_id, 'g': group}).to_csv(prep / 'annotations.tsv', sep='\t', index=False, header=False)
    auto.loc[genes, ['chromosome', 'start', 'end']].reset_index().to_csv(
        prep / 'gene_order.tsv', sep='\t', index=False, header=False)
    with gzip.open(prep / 'counts.mtx.gz', 'wb', compresslevel=1) as fh:
        sio.mmwrite(fh, X[send].T.tocoo())
    del X

    native = C.guard(out / 'native')
    t_prep = time.monotonic()
    # The R and JAGS modules set PYTHONPATH to a system Python, which breaks numpy here, so the
    # calling script loads them, saves RSCRIPT and R_LD_LIBRARY_PATH, then purges them.
    env = dict(os.environ, LD_LIBRARY_PATH=os.environ.get('R_LD_LIBRARY_PATH', os.environ.get('LD_LIBRARY_PATH', '')))
    env.pop('PYTHONPATH', None)
    subprocess.run([os.environ.get('RSCRIPT', 'Rscript'), '--vanilla', str(Path(__file__).with_name('infercnv_run.R')),
                    str(prep), str(native), str(args.ncores)], check=True, env=env)
    t_r = time.monotonic()

    eg = (native / 'expr_genes.txt').read_text().split('\n')[:-1]
    ec = (native / 'expr_cells.txt').read_text().split('\n')[:-1]
    E = np.fromfile(native / 'expr_cells_by_genes.f32', dtype=np.float32).reshape(len(ec), len(eg))
    bins, genes_df = C.gene_bins(np.asarray(eg), pos, args.bin_genes)
    gidx = pd.Index(eg).get_indexer(genes_df.gene)
    E = E[:, gidx]
    n_nan = int(np.isnan(E).sum())
    L = np.log2(np.maximum(np.nan_to_num(E, nan=1.0), 1e-3))   # NaN residual -> 1 (neutral); count in qc.json
    Zb = np.zeros((len(ec), len(bins)), dtype=np.float32)
    nb = np.bincount(genes_df.bin, minlength=len(bins))
    for j, b in enumerate(genes_df.bin.to_numpy()):
        Zb[:, b] += L[:, j]
    Zb /= nb
    del E, L
    row = pd.Index(ec).get_indexer(cells.cell_id)
    present = row >= 0
    Z = np.zeros((len(cells), len(bins)), dtype=np.float32)
    Z[present] = Zb[row[present]]
    del Zb
    q, h = query & present, hold & present
    if args.rule == 'all':
        q = query                                     # every query cell is analysed; absent cells have Z = 0
    score, malignant, clone, clusters = cnv_post.call(Z, cells, q, h, rule=args.rule, seed=C.SEED)
    profiles = C.clone_profiles(Z, cells, clone, h)
    normal_like = [c for c in pd.unique(clone) if str(c).endswith('_normal')]
    C.write_layout(out, cells, clone, malignant, bins, genes_df, profiles, normal_like, h)
    C.write_table(clusters, out / 'call_clusters.tsv')
    (native / 'expr_cells_by_genes.f32').unlink()          # large; bins are kept in profiles
    C.write_json(dict(tool='infercnv_r', n_cells=len(cells), n_sent=len(send), n_returned=int(present.sum()),
                      n_genes=len(eg), n_bins=len(bins), n_nan_values_set_neutral=n_nan, rule=args.rule,
                      n_malignant_query=int((malignant & query).sum()),
                      holdout_malignant_rate=float(malignant[h].mean()) if h.any() else None,
                      prep_seconds=t_prep - t0, r_seconds=t_r - t_prep, wall_seconds=time.monotonic() - t0,
                      python_max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2,
                      r_max_rss_gib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 ** 2),
                 out / 'qc.json')
    print(out)


if __name__ == '__main__':
    main()
