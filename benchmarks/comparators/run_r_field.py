#!/usr/bin/env python3
"""CopyKAT or SCEVAN comparator on the Haplodose input contract. Run with the scrna_cnv python; the R step
uses $FIELD_RSCRIPT (default: the cnv_field env Rscript, CopyKAT @01924c0 and SCEVAN @5a49b88).

Both tools were built for 10x scRNA-seq UMI counts with diploid normal cells as the baseline; here the
baseline is the fit reference cells (in Kinker, other cancer lines), as for inferCNV.
Cells sent: fit reference cells (as norm cells) and query plus held-out reference cells.
Post-processing: the tool's CNA matrix (CopyKAT: 220 kb windows; SCEVAN: genes) is mapped to each of our
genes by the nearest row on the same chromosome, averaged over consecutive 50-gene bins, then
cnv_post.call (malignant calls and Leiden clones) as for the inferCNV runners. With --clones native, query
cells take the tool's own subclone labels where it gives them (SCEVAN); CopyKAT has no subclone output, so
its clones stay Leiden.
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

RSCRIPT = '/gpfs/commons/home/jameslee/micromamba/envs/cnv_field/bin/Rscript'


def rows_to_genes(rows, genes_df):
    """Index of the nearest CNA-matrix row on the same chromosome for each gene (-1 if none)."""
    idx = np.full(len(genes_df), -1)
    for chrom, g in genes_df.groupby('chromosome'):
        r = rows.index[rows.chromosome.eq(chrom)].to_numpy()
        if not len(r):
            continue
        st = rows.start.to_numpy()[r]
        o = np.argsort(st)
        r, st = r[o], st[o]
        j = np.clip(np.searchsorted(st, g.start.to_numpy()), 1, len(st) - 1) if len(st) > 1 else np.zeros(len(g), int)
        if len(st) > 1:
            left = np.abs(g.start.to_numpy() - st[j - 1]) <= np.abs(st[j] - g.start.to_numpy())
            j = np.where(left, j - 1, j)
        idx[g.index.to_numpy()] = r[j]
    return idx


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--tool', required=True, choices=['copykat', 'scevan'])
    ap.add_argument('--counts', required=True)
    ap.add_argument('--cells', required=True)
    ap.add_argument('--positions', required=True)
    ap.add_argument('--out-parent', required=True, help='a new <out-parent>/<tool>_<timestamp> is created')
    ap.add_argument('--query-roles', default='tumour_epithelium,query')
    ap.add_argument('--bin-genes', type=int, default=50)
    ap.add_argument('--rule', choices=['cluster', 'cell', 'all'], default='cluster')
    ap.add_argument('--clones', choices=['native', 'leiden'], default='native')
    ap.add_argument('--ncores', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', '8')))
    args = ap.parse_args(argv)

    np.random.seed(C.SEED)
    t0 = time.monotonic()
    out = C.new_run_dir(args.out_parent, args.tool)
    C.provenance(out, args.tool, args)
    C.pip_freeze(out)
    prep = C.guard(out / 'prep')
    prep.mkdir()

    cells = C.read_cells(args.cells)
    X, genes, _ = C.read_counts(args.counts, cells.cell_id)
    X, genes = C.collapse_genes(X, genes)
    pos = C.read_positions(args.positions)
    auto = pos[pos.chromosome.isin(C.AUTOSOMES)].set_index('gene')
    keep = pd.Index(genes).isin(auto.index)
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
    with gzip.open(prep / 'counts.mtx.gz', 'wb', compresslevel=1) as fh:
        sio.mmwrite(fh, X[send].T.tocoo())
    del X

    native = C.guard(out / 'native')
    t_prep = time.monotonic()
    env = {k: v for k, v in os.environ.items() if k not in ('PYTHONPATH', 'R_LIBS', 'R_LIBS_USER', 'LD_LIBRARY_PATH')}
    script = Path(__file__).with_name(f'{args.tool}_run.R')
    subprocess.run([os.environ.get('FIELD_RSCRIPT', RSCRIPT), '--vanilla', str(script), str(prep), str(native),
                    str(args.ncores)], check=True, env=env)
    t_r = time.monotonic()

    rows = pd.read_csv(native / 'rows.tsv', sep='\t')
    mc = (native / 'cells.txt').read_text().split('\n')[:-1]
    M = np.fromfile(native / 'mat.f32', dtype=np.float32).reshape(len(mc), len(rows))
    bins, genes_df = C.gene_bins(np.asarray(genes), pos, args.bin_genes)
    gi = rows_to_genes(rows, genes_df)
    n_unmapped = int((gi < 0).sum())
    G = np.where(gi[None] >= 0, M[:, np.maximum(gi, 0)], 0.0)
    G = np.nan_to_num(G, nan=0.0)
    Zb = np.zeros((len(mc), len(bins)), dtype=np.float32)
    for j, b in enumerate(genes_df.bin.to_numpy()):
        Zb[:, b] += G[:, j]
    Zb /= np.bincount(genes_df.bin, minlength=len(bins))
    del M, G
    row = pd.Index(mc).get_indexer(cells.cell_id)
    present = row >= 0
    Z = np.zeros((len(cells), len(bins)), dtype=np.float32)
    Z[present] = Zb[row[present]]
    del Zb
    q, h = query & present, hold & present
    if args.rule == 'all':
        q = query                                     # every query cell is analysed; absent cells have Z = 0
    score, malignant, clone, clusters = cnv_post.call(Z, cells, q, h, rule=args.rule, seed=C.SEED)
    calls = pd.read_csv(native / 'calls.tsv', sep='\t').set_index('cell_id')
    sub = calls.subclone.reindex(cells.cell_id).to_numpy()
    n_native = 0
    if args.clones == 'native' and pd.notna(sub).any():
        spec = cells.specimen.to_numpy()
        m = malignant & query
        lab = np.where(pd.notna(sub), pd.Series(sub).astype(str).str.replace(r'\.0$', '', regex=True), 'none')
        clone[m] = [f'{s}_c{x}' for s, x in zip(spec[m], lab[m])]
        n_native = int(m.sum())
    profiles = C.clone_profiles(Z, cells, clone, h)
    normal_like = [c for c in pd.unique(clone) if str(c).endswith('_normal')]
    C.write_layout(out, cells, clone, malignant, bins, genes_df, profiles, normal_like, h)
    C.write_table(clusters, out / 'call_clusters.tsv')
    (native / 'mat.f32').unlink()                     # large; bins are kept in profiles
    C.write_json(dict(tool=args.tool, n_cells=len(cells), n_sent=len(send), n_in_matrix=int(present.sum()),
                      n_rows=len(rows), n_bins=len(bins), n_genes_unmapped=n_unmapped, rule=args.rule,
                      clones=args.clones, n_query_native_clone=n_native,
                      native_malignant_query=float(calls.malignant.reindex(cells.cell_id[query]).mean()),
                      n_malignant_query=int((malignant & query).sum()),
                      holdout_malignant_rate=float(malignant[h].mean()) if h.any() else None,
                      prep_seconds=t_prep - t0, r_seconds=t_r - t_prep, wall_seconds=time.monotonic() - t0,
                      python_max_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2,
                      r_max_rss_gib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024 ** 2),
                 out / 'qc.json')
    print(out)


if __name__ == '__main__':
    main()
