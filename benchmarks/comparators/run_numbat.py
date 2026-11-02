#!/usr/bin/env python3
"""Numbat comparator on one specimen with phased allele counts (~/CRC/data/interim/numbat_env).

Numbat is built for scRNA-seq UMI counts plus phased SNP pileups, like Haplodose with
--alleles. It needs real hg38 gene symbols and a Numbat allele table (pileup_and_phase.R
output, *_allele_counts.tsv.gz), so it does not run on `haplodose simulate` output.

Cells passed to Numbat: the specimen's query cells plus its held-out reference cells,
with library size >= --min-umi (300 in the CRC call pattern). The expression reference
is aggregate_counts over fit reference cells by ref_group. Genes in --exclude-genes are
dropped from counts and from the allele table.

Common output: clone = <specimen>_c<clone_opt> from the last clone_post_<i>.tsv,
malignant_call = compartment_opt == 'tumor'. arms.tsv per clone from
bulk_clones_final.tsv.gz in the Haplodose layout (cmp_common.write_layout), units 50-gene bins
(see layout_profiles). Clones with < 50% tumour cells are normal_like.
"""
from __future__ import annotations

import argparse
import gzip
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmp_common as C  # noqa: E402

ENV = Path('/gpfs/commons/home/jameslee/CRC/data/interim/numbat_env')


def last_iteration(numbat_out, stem):
    files = sorted(numbat_out.glob(f'{stem}_*.tsv'), key=lambda p: int(p.stem.split('_')[-1]))
    if not files:
        raise FileNotFoundError(f'No {stem}_*.tsv in {numbat_out}')
    return files[-1]


NUMBAT_CN = {'amp': 3, 'bamp': 3, 'del': 1, 'bdel': 1, 'neu': 2, 'loh': 2}


def layout_profiles(numbat_out, specimen, pos, bin_genes):
    """Per clone and 50-gene bin: log_ratio = mean Numbat gene lnFC / ln 2 (log2), cn_state = modal
    cnv_state_post mapped to 1 loss / 2 neutral (incl. cnLOH) / 3 gain."""
    bulk = pd.read_csv(numbat_out / 'bulk_clones_final.tsv.gz', sep='\t', low_memory=False)
    g = bulk[bulk.gene.notna() & (bulk.gene.astype(str) != '') & bulk.lnFC.notna()].drop_duplicates(['sample', 'gene'])
    bins, genes_df = C.gene_bins(g.gene.unique(), pos, bin_genes)
    g = g.merge(genes_df[['gene', 'bin']], on='gene')
    g['cn'] = g.cnv_state_post.map(NUMBAT_CN).fillna(2).astype(int)
    prof = (g.groupby(['sample', 'bin']).agg(log_ratio=('lnFC', lambda v: float(v.mean() / np.log(2))),
                                             cn_state=('cn', lambda v: int(v.mode().iloc[0]))).reset_index())
    prof['clone'] = specimen + '_c' + prof['sample'].astype(str)
    prof['specimen'] = specimen
    return bins, genes_df, prof[['specimen', 'clone', 'bin', 'log_ratio', 'cn_state']]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--counts', required=True)
    ap.add_argument('--cells', required=True)
    ap.add_argument('--positions', required=True, help='contract input; Numbat uses its own hg38 gene table')
    ap.add_argument('--alleles', required=True, help='Numbat allele table, e.g. <sample>_allele_counts.tsv.gz')
    ap.add_argument('--specimen', required=True)
    ap.add_argument('--out-parent', required=True, help='a new <out-parent>/numbat_<timestamp> is created')
    ap.add_argument('--exclude-genes', help='one gene per line, dropped from counts and alleles')
    ap.add_argument('--query-roles', default='tumour_epithelium,query')
    ap.add_argument('--min-umi', type=float, default=300)
    ap.add_argument('--bin-genes', type=int, default=50)
    ap.add_argument('--ncores', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', '8')))
    ap.add_argument('--genome', default='hg38', choices=['hg38', 'hg19'], help='Numbat genome (alleles and genes)')
    ap.add_argument('--reference-scope', default='patient', choices=['patient', 'dataset'],
                    help="fit reference cells of the specimen's patient (CRC, in-sample normals) or of its dataset "
                         '(Kinker: the pool\'s other lines)')
    args = ap.parse_args(argv)

    t0 = time.monotonic()
    out = C.new_run_dir(args.out_parent, 'numbat')
    C.provenance(out, 'numbat', args, env=str(ENV))
    C.pip_freeze(out)
    prep = out / 'prep'
    C.guard(prep).mkdir()

    cells = C.read_cells(args.cells)
    X, genes, _ = C.read_counts(args.counts, cells.cell_id)
    X, genes = C.collapse_genes(X, genes)
    exclude = set(Path(args.exclude_genes).read_text().split()) if args.exclude_genes else set()
    if exclude:
        keep = ~pd.Index(genes).isin(exclude)
        X, genes = X[:, keep], genes[keep]
    lib = np.asarray(X.sum(axis=1)).ravel()
    spec = cells.specimen.eq(args.specimen).to_numpy()
    query = C.query_mask(cells, args.query_roles.split(','))
    is_ref = cells.role.eq('reference').to_numpy()
    fit = is_ref & cells.ref_split.eq('fit').to_numpy()
    hold = is_ref & cells.ref_split.eq('holdout').to_numpy()
    passed = spec & (query | hold) & (lib >= args.min_umi)
    if args.reference_scope == 'dataset':
        scope = (cells.dataset == cells.dataset[spec].iloc[0]).to_numpy()
    else:
        scope = spec | (cells.patient == cells.patient[spec].iloc[0]).to_numpy()
    ref_rows = np.flatnonzero(fit & (lib >= args.min_umi) & scope)
    q_rows = np.flatnonzero(passed)
    if len(ref_rows) < 40 or len(q_rows) < 50:
        raise ValueError(f'Too few cells: {len(ref_rows)} reference, {len(q_rows)} query')
    pd.DataFrame({'gene': genes}).to_csv(C.guard(prep / 'genes.tsv'), sep='\t', index=False)
    pd.DataFrame({'cell': cells.cell_id.to_numpy()[q_rows]}).to_csv(C.guard(prep / 'query_cells.tsv'), sep='\t', index=False)
    pd.DataFrame({'cell': cells.cell_id.to_numpy()[ref_rows], 'group': cells.ref_group.to_numpy()[ref_rows]}).to_csv(
        C.guard(prep / 'ref_cells.tsv'), sep='\t', index=False)
    for name, rows in [('query.mtx.gz', q_rows), ('ref.mtx.gz', ref_rows)]:
        with gzip.open(C.guard(prep / name), 'wb') as fh:
            sio.mmwrite(fh, X[rows].T.tocoo())
    alleles = pd.read_csv(args.alleles, sep='\t', dtype={'CHROM': str})
    alleles = alleles[alleles.cell.isin(cells.cell_id.to_numpy()[q_rows]) & ~alleles.gene.isin(exclude)]
    alleles.to_csv(C.guard(prep / 'alleles.tsv.gz'), sep='\t', index=False, compression='gzip')
    del X

    numbat_out = out / 'native'
    env = dict(os.environ, PATH=f'{ENV}/bin:' + os.environ['PATH'], R_LIBS_USER=f'{ENV}/lib/R/library',
               OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    env.pop('LD_LIBRARY_PATH', None)
    subprocess.run([str(ENV / 'bin/Rscript'), '--vanilla', str(Path(__file__).with_name('numbat_run.R')),
                    str(prep), str(C.guard(numbat_out)), str(args.ncores), args.genome], check=True, env=env)

    post = pd.read_csv(last_iteration(numbat_out, 'clone_post'), sep='\t', dtype={'clone_opt': str})
    post = post.set_index('cell').reindex(cells.cell_id)
    called = post.clone_opt.notna().to_numpy()
    clone = np.where(called, args.specimen + '_c' + post.clone_opt.fillna('').astype(str), 'NA').astype(object)
    malignant = (post.compartment_opt == 'tumor').fillna(False).to_numpy()
    bins, genes_df, profiles = layout_profiles(numbat_out, args.specimen, C.read_positions(args.positions),
                                               args.bin_genes)
    comp = pd.DataFrame({'clone': clone[called], 'tumour': malignant[called]}).groupby('clone').tumour.mean()
    normal_like = comp.index[comp < 0.5].tolist()
    C.write_layout(out, cells, clone, malignant, bins, genes_df, profiles, normal_like, hold & called)
    C.write_json(dict(tool='numbat', n_query=len(q_rows), n_reference=len(ref_rows), n_called=int(called.sum()),
                      n_malignant=int(malignant.sum()),
                      holdout_malignant_rate=float(malignant[hold & called].mean()) if (hold & called).any() else None,
                      wall_seconds=time.monotonic() - t0), out / 'qc.json')
    print(out)


if __name__ == '__main__':
    main()
