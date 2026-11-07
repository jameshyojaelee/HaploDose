"""HaploDose as one more caller in the Schmid et al. 2025 benchmark pipeline.

Pipeline: github.com/colomemaria/benchmark_scrnaseq_cnv_callers, pinned in external/ (see
benchmarks/schmid/PLAN.md). Every caller there reads the same three inputs and writes its own files,
which a per-tool reader in their scripts/lib.R turns into one pseudobulk score per gene (2 = neutral).

Inputs (their format, as written by their preprocessing_datasets/*.R):
  count_matrix.txt       genes x cells, tab separated, gene symbols as row names; R write.table layout,
                         so the header has one field fewer than the data rows
  sample_annotation.txt  no header: cell, group
  ref_groups.txt         header 'ref_groups': the groups that are normal reference cells
  gene positions         their inferCNV gene order file (data/annotations/hg38_gencode_v27.txt):
                         no header: gene, chromosome, start, end
Reference cells are the cells whose group is in ref_groups, exactly as for inferCNV, copyKat, SCEVAN and
Numbat in their scripts. Each reference group is one HaploDose reference type (ref_group); within each
group a seeded share of cells is held out (ref_split holdout, default 0.3, as in the Kinker and MIX-seq
runs) for HaploDose's own reference check. Every other annotated cell is a query cell of one specimen.

Outputs in OUT_DIR (read by benchmarks/schmid/lib_haplodose.R, read_haplodose):
  <dataset>_haplodose_gene_by_clone.tsv  gene, chromosome, start, end, then one column per clone: the
                                         clone's profile value (COLUMN, natural log ratio, 0 = neutral)
                                         at the bin that holds the gene; only genes HaploDose used
  <dataset>_haplodose_cell_clone.tsv     cell, clone: query cells with an assigned clone
  haplodose_run/                         the full HaploDose output directory
  run_haplodose.json                     arguments, code path and SHA-256 of every code file, timings
The reader averages the clone values over the dataset's non-reference cells (each cell takes its MAP
clone's profile) and adds 2, as read_copykat and read_numbat_expr do.

The HaploDose package is taken from --code (put first on PYTHONPATH), never from the environment's
install; the run stops if the imported package is not under --code. CPU by default (JAX_PLATFORMS=cpu).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

PY = '/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python'
PROTECTED = '/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets'


def read_count_matrix(path, chunk_rows=2000):
    """count_matrix.txt -> (cells x genes CSR int32, genes, cells), read in row chunks."""
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        first = fh.readline().rstrip('\n').split('\t')
    n_data = len(first)
    if len(header) == n_data - 1:             # R write.table(row.names=TRUE): no field above the row names
        cells = header
    elif len(header) == n_data:               # a leading field (empty, 'V1' or 'gene') above the row names
        cells = header[1:]
    else:
        raise ValueError(f'header has {len(header)} fields, data rows {n_data}')
    genes, blocks = [], []
    reader = pd.read_csv(path, sep='\t', header=None, skiprows=1, chunksize=chunk_rows,
                         dtype={0: str}, engine='c')
    for part in reader:
        if part.shape[1] != len(cells) + 1:
            raise ValueError('ragged count matrix')
        genes.extend(part.iloc[:, 0].astype(str))
        v = part.iloc[:, 1:].to_numpy(dtype=np.float64)
        if not np.isfinite(v).all() or (v < 0).any():
            raise ValueError('counts must be finite and nonnegative')
        if np.abs(v - np.round(v)).max(initial=0.0) > 1e-6:
            raise ValueError('counts must be integers (raw UMI counts)')
        blocks.append(sparse.csr_matrix(v.astype(np.int32)))
    X = sparse.vstack(blocks).T.tocsr()       # cells x genes
    return X, np.array(genes), np.array(cells)


def write_csr_h5(path, X, genes, cells):
    with h5py.File(path, 'w') as h:
        g = h.create_group('X')
        g.create_dataset('data', data=X.data.astype(np.int32), compression='gzip')
        g.create_dataset('indices', data=X.indices.astype(np.int32), compression='gzip')
        g.create_dataset('indptr', data=X.indptr.astype(np.int64))
        h.create_dataset('genes', data=np.array(genes, dtype='S'))
        h.create_dataset('cell_ids', data=np.array(cells, dtype='S'))


def cell_table(cells, annot_path, ref_path, dataset, holdout, seed):
    annot = pd.read_csv(annot_path, sep='\t', header=None, dtype=str, keep_default_na=False).iloc[:, :2]
    annot.columns = ['cell', 'group']
    if annot.cell.duplicated().any():
        raise ValueError('duplicated cells in sample_annotation.txt')
    refs = set(pd.read_csv(ref_path, sep='\t', dtype=str).ref_groups)
    group = annot.set_index('cell').group.reindex(cells)
    meta = pd.DataFrame({'cell_id': cells, 'specimen': dataset, 'patient': dataset, 'dataset': dataset,
                         'role': 'other', 'ref_group': '', 'ref_split': '', 'group': group.fillna('').to_numpy()})
    is_ref = group.isin(refs).to_numpy()
    is_query = group.notna().to_numpy() & ~is_ref
    meta.loc[is_query, 'role'] = 'query'
    meta.loc[is_ref, 'role'] = 'reference'
    meta.loc[is_ref, 'ref_group'] = meta.group[is_ref]
    rng = np.random.default_rng(seed)
    for g, idx in meta[is_ref].groupby('group').groups.items():
        idx = np.asarray(idx)
        n_hold = int(round(holdout * len(idx)))
        hold = rng.choice(idx, size=n_hold, replace=False) if n_hold else np.array([], dtype=int)
        meta.loc[idx, 'ref_split'] = 'fit'
        meta.loc[hold, 'ref_split'] = 'holdout'
    counts = dict(n_annotated=int(group.notna().sum()), n_not_annotated=int(group.isna().sum()),
                  n_annotated_missing_from_matrix=int((~annot.cell.isin(cells)).sum()),
                  n_query=int(is_query.sum()), n_reference=int(is_ref.sum()),
                  reference_groups=sorted(refs), groups_found=sorted(meta.group[is_ref].unique()))
    missing_groups = refs - set(meta.group[is_ref])
    if missing_groups:
        raise ValueError(f'ref_groups with no cells in the matrix: {sorted(missing_groups)}')
    return meta.drop(columns='group'), counts


def positions_table(gene_pos_path):
    # names has a spare fifth column: hg38_gencode_v27.txt has one row (TNFRSF18) with a trailing tab
    pos = pd.read_csv(gene_pos_path, sep='\t', header=None, names=['gene', 'chromosome', 'start', 'end', 'extra'],
                      dtype={'gene': str, 'chromosome': str})[['gene', 'chromosome', 'start', 'end']]
    pos['chromosome'] = np.where(pos.chromosome.str.startswith('chr'), pos.chromosome, 'chr' + pos.chromosome)
    return pos.drop_duplicates('gene', keep='first')


def code_hashes(code):
    rows = []
    for f in sorted(Path(code).rglob('*')):
        if f.is_file() and '__pycache__' not in f.parts:
            rows.append(dict(file=str(f.relative_to(code)), sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
    return rows


def convert_outputs(run, out, dataset, column):
    genes = pd.read_csv(run / 'genes.tsv.gz', sep='\t')
    prof = pd.read_csv(run / 'profiles.tsv.gz', sep='\t')
    cells = pd.read_csv(run / 'cells.tsv.gz', sep='\t', dtype={'cell_id': str, 'clone': str},
                        keep_default_na=False)
    if column not in prof.columns:
        raise ValueError(f'profiles.tsv.gz has no column {column}')
    wide = prof.pivot(index='bin', columns='clone', values=column)
    clones = list(wide.columns)
    table = genes[['gene', 'chromosome', 'start', 'end']].copy()
    table[['start', 'end']] = table[['start', 'end']].astype(np.int64)
    vals = wide.reindex(genes.bin.to_numpy()).to_numpy()
    for i, c in enumerate(clones):
        table[c] = vals[:, i]
    table.to_csv(out / f'{dataset}_haplodose_gene_by_clone.tsv', sep='\t', index=False)
    assigned = cells[(cells.role == 'query') & (cells.status == 'assigned') & cells.clone.isin(clones)]
    assigned[['cell_id', 'clone']].rename(columns={'cell_id': 'cell'}).to_csv(
        out / f'{dataset}_haplodose_cell_clone.tsv', sep='\t', index=False)
    return dict(n_genes_out=len(table), n_clones=len(clones), n_query_assigned=len(assigned),
                n_query=int((cells.role == 'query').sum()),
                status_counts=cells[cells.role == 'query'].status.value_counts().to_dict())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--matrix', required=True)
    ap.add_argument('--annot', required=True)
    ap.add_argument('--ref-groups', required=True)
    ap.add_argument('--gene-pos', required=True)
    ap.add_argument('--out-dir', required=True, help='new or empty; e.g. results/output_<dataset>/haplodose')
    ap.add_argument('--dataset', required=True)
    ap.add_argument('--code', required=True, help='directory that contains the haplodose package (frozen copy)')
    ap.add_argument('--column', default='log_ratio_gene', choices=['log_ratio_gene', 'log_ratio', 'log_dose',
                                                                   'assignment_dose'])
    ap.add_argument('--ref-holdout', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=20261103)
    ap.add_argument('--jax-platform', default='cpu', choices=['cpu', 'cuda'])
    ap.add_argument('--hd-args', default='', help='extra haplodose run arguments (default: none, CLI defaults)')
    a = ap.parse_args()

    out = Path(a.out_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f'output directory is not empty: {out}')
    out.mkdir(parents=True, exist_ok=True)
    code = Path(a.code).resolve()
    env = dict(os.environ, PYTHONPATH=str(code), PYTHONNOUSERSITE='1', PYTHONDONTWRITEBYTECODE='1',
               JAX_PLATFORMS=a.jax_platform, XLA_PYTHON_CLIENT_PREALLOCATE='false',
               HAPLODOSE_PROTECTED_PATHS=PROTECTED)
    imported = subprocess.run([PY, '-s', '-c', 'import haplodose, os; print(os.path.realpath(haplodose.__file__))'],
                              env=env, capture_output=True, text=True, check=True).stdout.strip()
    if not Path(imported).is_relative_to(code):    # code is resolve()d, so symlinked code paths compare by target
        raise RuntimeError(f'haplodose imported from {imported}, not from --code {code}')

    t0 = time.monotonic()
    log = dict(args=vars(a), code=str(code), code_realpath=os.path.realpath(code), imported=imported,
               started=dt.datetime.now().isoformat(timespec='seconds'), host=os.uname().nodename)
    inp = out / 'input'
    inp.mkdir()
    X, genes, cells = read_count_matrix(a.matrix)
    if pd.Index(cells).duplicated().any() or pd.Index(genes).duplicated().any():
        raise ValueError('duplicated cell or gene names in count_matrix.txt')
    write_csr_h5(inp / 'counts.h5', X, genes, cells)
    meta, cell_counts = cell_table(cells, a.annot, a.ref_groups, a.dataset, a.ref_holdout, a.seed)
    meta.to_csv(inp / 'cells.tsv', sep='\t', index=False)
    pos = positions_table(a.gene_pos)
    pos.to_csv(inp / 'positions.tsv', sep='\t', index=False)
    log.update(n_cells=len(cells), n_genes=len(genes), n_genes_with_position=int(pd.Index(genes).isin(pos.gene).sum()),
               cells=cell_counts, prepare_seconds=time.monotonic() - t0)
    del X

    run = out / 'haplodose_run'
    cmd = [PY, '-s', '-m', 'haplodose.cli', 'run', '--counts', str(inp / 'counts.h5'), '--cells',
           str(inp / 'cells.tsv'), '--positions', str(inp / 'positions.tsv'), '--out', str(run)] + a.hd_args.split()
    log['command'] = cmd
    t1 = time.monotonic()
    subprocess.run(cmd, env=env, check=True)
    log['haplodose_seconds'] = time.monotonic() - t1
    log.update(convert_outputs(run, out, a.dataset, a.column))
    log['code_sha256'] = code_hashes(os.path.realpath(code))
    log['finished'] = dt.datetime.now().isoformat(timespec='seconds')
    log['wall_seconds'] = time.monotonic() - t0
    (out / 'run_haplodose.json').write_text(json.dumps(log, indent=2, default=str) + '\n')
    print(json.dumps({k: log[k] for k in ('n_cells', 'n_genes', 'n_query_assigned', 'n_clones', 'wall_seconds')}))


if __name__ == '__main__':
    sys.exit(main())
