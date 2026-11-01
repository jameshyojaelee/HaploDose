"""Inputs, protected-path guard and atomic outputs.

Counts stay sparse on disk and are read in bounded cell chunks. Every write goes
through `guard`, which refuses paths under HAPLODOSE_PROTECTED_PATHS after
resolving symlinks, including for paths that do not exist yet.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse

REQUIRED_CELL_COLUMNS = ['cell_id', 'specimen', 'patient', 'role', 'ref_group', 'ref_split']
AUTOSOMES = [f'chr{i}' for i in range(1, 23)]


def protected_roots():
    raw = os.environ.get('HAPLODOSE_PROTECTED_PATHS', '')
    return [Path(os.path.realpath(p)) for p in raw.split(os.pathsep) if p]


def guard(path) -> Path:
    path = Path(os.path.realpath(os.path.abspath(path)))
    for root in protected_roots():
        if path == root or root in path.parents:
            raise PermissionError(f'Refusing to write under protected path {root}: {path}')
    return path


def new_output_dir(path) -> Path:
    path = guard(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f'Output directory is not empty: {path}')
    path.mkdir(parents=True, exist_ok=True)
    return path


def write_table(frame: pd.DataFrame, path) -> Path:
    path = guard(path)
    tmp = guard(path.with_name('.' + path.name + '.tmp'))
    frame.to_csv(tmp, sep='\t', index=False,
                 compression='gzip' if path.suffix == '.gz' else None)
    os.replace(tmp, path)
    return path


def write_json(obj, path) -> Path:
    path = guard(path)
    tmp = guard(path.with_name('.' + path.name + '.tmp'))
    tmp.write_text(json.dumps(obj, indent=2, default=_json_default) + '\n')
    os.replace(tmp, path)
    return path


def _json_default(x):
    if isinstance(x, (np.integer, np.floating, np.bool_)):
        return x.item()
    if isinstance(x, np.ndarray):
        return x.tolist()
    raise TypeError(type(x))


def open_memmap(path, dtype, shape, fill=None):
    arr = np.lib.format.open_memmap(guard(path), mode='w+', dtype=dtype, shape=shape)
    if fill is not None:
        arr[:] = fill
    return arr


class CountSource:
    """Cells x genes counts from 10x H5, h5ad (CSR X) or CSR H5 (X/*, genes, cell_ids).

    10x stores genes x cells CSC, which is the same layout as cells x genes CSR.
    """

    def __init__(self, path):
        self.path = Path(path)
        self.h = h5py.File(self.path, 'r')
        h = self.h
        if 'matrix' in h:
            self.g = h['matrix']
            self.genes = self.g['features/name'].asstr()[:]
            self.cell_ids = self.g['barcodes'].asstr()[:]
        elif 'obs' in h and 'var' in h:
            self.g = h['X']
            if self.g.attrs.get('encoding-type') != 'csr_matrix':
                raise ValueError('h5ad X must be stored as csr_matrix')
            self.genes = _h5ad_index(h['var'])
            self.cell_ids = _h5ad_index(h['obs'])
        else:
            self.g = h['X']
            self.genes = h['genes'].asstr()[:]
            self.cell_ids = h['cell_ids'].asstr()[:]
        self.n_cells, self.n_genes = len(self.cell_ids), len(self.genes)
        self.indptr = self.g['indptr'][:].astype(np.int64)
        if len(self.indptr) != self.n_cells + 1:
            raise ValueError('Expected cell-major sparse storage (indptr length = cells + 1)')

    def read(self, start, stop):
        lo, hi = int(self.indptr[start]), int(self.indptr[stop])
        X = sparse.csr_matrix((self.g['data'][lo:hi].astype(np.float32),
                               self.g['indices'][lo:hi].astype(np.int32),
                               self.indptr[start:stop + 1] - lo),
                              shape=(stop - start, self.n_genes))
        if X.data.size and (not np.isfinite(X.data).all() or X.data.min() < 0):
            raise ValueError('Counts must be finite and nonnegative')
        return X

    def chunks(self, size):
        for start in range(0, self.n_cells, size):
            stop = min(start + size, self.n_cells)
            yield start, stop, self.read(start, stop)

    def close(self):
        self.h.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def _h5ad_index(group):
    key = group.attrs.get('_index', '_index')
    return group[key].asstr()[:]


def read_cells(path, cell_ids) -> pd.DataFrame:
    meta = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    missing = [c for c in REQUIRED_CELL_COLUMNS if c not in meta.columns]
    if missing:
        raise ValueError(f'cells table lacks columns: {missing}')
    if meta.cell_id.duplicated().any():
        raise ValueError('Duplicated cell IDs')
    if not np.array_equal(meta.cell_id.to_numpy(), np.asarray(cell_ids, dtype=str)):
        raise ValueError('cells table must list cells in exact count-matrix order')
    if 'dataset' not in meta.columns:
        meta['dataset'] = 'dataset'
    return meta


def read_positions(path) -> pd.DataFrame:
    pos = pd.read_csv(path, sep='\t')
    pos = pos.rename(columns={pos.columns[0]: 'gene'})
    for c in ['chromosome', 'start', 'end']:
        if c not in pos.columns:
            raise ValueError(f'positions table lacks {c}')
    pos['chromosome'] = pos.chromosome.astype(str)
    pos.loc[~pos.chromosome.str.startswith('chr'), 'chromosome'] = 'chr' + pos.chromosome
    return pos.drop_duplicates('gene', keep='first')


def read_site_list(path):
    """Genomic sites, one per line as chrN:pos (or N:pos), normalized to chrN:pos."""
    if path is None:
        return set()
    sites = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)[0].str.strip()
    return {x if x.startswith('chr') else 'chr' + x for x in sites}


def read_gene_list(path):
    if path is None:
        return set()
    names = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)[0]
    return set(names.str.strip().str.upper())
