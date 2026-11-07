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
GENE_NAMESPACES = ('symbol', 'ensembl', 'custom')
COORDINATE_SYSTEMS = ('one-based-closed', 'zero-based-half-open')


def validate_feature_ids(values, gene_namespace='custom'):
    """Canonical, unique feature identities; ambiguous loci must be resolved upstream."""
    if gene_namespace not in GENE_NAMESPACES:
        raise ValueError(f'gene_namespace must be one of {GENE_NAMESPACES}')
    raw = pd.Series(values, dtype=object)
    if raw.isna().any() or raw.map(lambda x: not isinstance(x, str) or not x.strip()).any():
        raise ValueError('Feature IDs must be nonempty strings')
    names = raw.astype(str).str.strip()
    if names.str.lower().isin(['nan', 'none', 'null', 'na']).any():
        raise ValueError('Feature IDs must not be missing-value placeholders')
    if gene_namespace == 'symbol':
        names = names.str.upper()
    if names.duplicated().any():
        raise ValueError(f'Duplicated canonical feature IDs: {names[names.duplicated()].iloc[0]}; '
                         'aggregate raw counts or resolve feature identities upstream')
    if gene_namespace == 'ensembl' and not names.str.match(r'^ENS[A-Z]*G\d+(\.\d+)?$').all():
        raise ValueError('Ensembl gene IDs must use ENS...G<number> (optional version suffix)')
    return names.to_numpy(dtype=str)


def validate_input_contract(counts_kind, genome_build, coordinate_system, gene_namespace, centromeres=None):
    if counts_kind != 'raw':
        raise ValueError('counts_kind must explicitly be raw; normalized or estimated counts are unsupported')
    if not isinstance(genome_build, str) or not genome_build.strip():
        raise ValueError('genome_build must be explicitly declared')
    if coordinate_system not in COORDINATE_SYSTEMS:
        raise ValueError(f'coordinate_system must explicitly be one of {COORDINATE_SYSTEMS}')
    if gene_namespace not in GENE_NAMESPACES:
        raise ValueError(f'gene_namespace must explicitly be one of {GENE_NAMESPACES}')
    if genome_build != 'hg38' and centromeres is None:
        raise ValueError('Only hg38 has bundled centromeres; other builds require an explicit centromeres table')
    return dict(counts_kind=counts_kind, genome_build=genome_build, coordinate_system=coordinate_system,
                gene_namespace=gene_namespace, centromeres=None if centromeres is None else str(centromeres))


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
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)
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

    def __init__(self, path, *, gene_namespace='custom'):
        try:
            self._open(path, gene_namespace=gene_namespace)
        except Exception:
            if hasattr(self, 'h'):
                self.h.close()
            raise

    def _open(self, path, *, gene_namespace):
        self.path = Path(path)
        self.h = h5py.File(self.path, 'r')
        h = self.h
        if 'matrix' in h:
            self.g = h['matrix']
            feature_key = 'features/id' if gene_namespace == 'ensembl' else 'features/name'
            self.genes = self.g[feature_key].asstr()[:]
            if 'feature_type' in self.g['features']:
                feature_types = self.g['features/feature_type'].asstr()[:]
                if np.any(feature_types != 'Gene Expression'):
                    self.h.close()
                    raise ValueError('10x input includes non-Gene Expression features; subset RNA features upstream')
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
        self.genes = validate_feature_ids(self.genes, gene_namespace)
        ids = pd.Series(self.cell_ids, dtype=str)
        if (ids.str.strip().eq('').any() or ids.duplicated().any()
                or ids.str.lower().isin(['nan', 'none', 'null', 'na']).any()):
            raise ValueError('Count matrix cell IDs must be nonempty and unique')
        self.n_cells, self.n_genes = len(self.cell_ids), len(self.genes)
        if not self.n_cells or not self.n_genes:
            raise ValueError('Count matrix must have cells and genes')
        if any(self.g[key].ndim != 1 for key in ['indptr', 'indices', 'data']):
            raise ValueError('Sparse indptr, indices and data must be one-dimensional')
        pointers = self.g['indptr'][:]
        if pointers.dtype.kind not in 'iu':
            raise ValueError('Sparse indptr must have integer storage')
        self.indptr = pointers.astype(np.int64)
        if len(self.indptr) != self.n_cells + 1:
            raise ValueError('Expected cell-major sparse storage (indptr length = cells + 1)')
        if (self.indptr[0] != 0 or np.any(self.indptr < 0) or np.any(np.diff(self.indptr) < 0)
                or self.indptr[-1] != len(self.g['data']) or len(self.g['indices']) != len(self.g['data'])):
            raise ValueError('Invalid CSR pointers or sparse data/indices lengths')
        if self.g['indices'].dtype.kind not in 'iu':
            raise ValueError('Sparse indices must have integer storage')
        if 'shape' in self.g:
            expected = (self.n_genes, self.n_cells) if 'matrix' in h else (self.n_cells, self.n_genes)
            if tuple(self.g['shape'][:]) != expected:
                raise ValueError('Sparse matrix shape disagrees with cell/feature IDs')
        elif 'shape' in self.g.attrs and tuple(self.g.attrs['shape']) != (self.n_cells, self.n_genes):
            raise ValueError('Sparse matrix shape disagrees with cell/feature IDs')

    def read(self, start, stop):
        if not 0 <= start <= stop <= self.n_cells:
            raise ValueError('Count chunk bounds are outside the matrix')
        lo, hi = int(self.indptr[start]), int(self.indptr[stop])
        data, indices = self.g['data'][lo:hi], self.g['indices'][lo:hi]
        if indices.size and (indices.min() < 0 or indices.max() >= self.n_genes):
            raise ValueError('Sparse gene index out of range')
        if data.size and (not np.isfinite(data).all() or data.min() < 0):
            raise ValueError('Counts must be finite and nonnegative')
        if data.size and np.any(data != np.floor(data)):
            raise ValueError('Counts must be raw integer UMI/read counts; got fractional values')
        # float32 is exact only through 2**24; do not silently round before splitting.
        if data.size and np.any(data.astype(np.float32).astype(np.float64) != data.astype(np.float64)):
            raise ValueError('Counts cannot be represented exactly in float32; input exceeds supported precision')
        X = sparse.csr_matrix((data.astype(np.float32), indices.astype(np.int32),
                               self.indptr[start:stop + 1] - lo),
                              shape=(stop - start, self.n_genes))
        X.sort_indices()
        if not X.has_canonical_format:
            raise ValueError('Sparse rows contain duplicate feature entries; aggregate integers upstream')
        if X.data.size and (not np.isfinite(X.data).all() or X.data.min() < 0):
            raise ValueError('Counts must be finite and nonnegative')
        if X.data.size and np.any(X.data != np.floor(X.data)):
            raise ValueError('Counts must be raw integer UMI/read counts; got fractional values (normalised or '
                             'estimated counts are not supported: binning and count splitting assume integers)')
        return X

    def chunks(self, size):
        if not isinstance(size, (int, np.integer)) or size <= 0:
            raise ValueError('Count chunk size must be a positive integer')
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
    for col in ['cell_id', 'specimen', 'patient', 'role']:
        if meta[col].str.strip().eq('').any() or meta[col].str.lower().isin(['nan', 'none', 'null', 'na']).any():
            raise ValueError(f'cells table has missing {col}')
    ref = meta.role.eq('reference')
    if (ref & ~meta.ref_split.isin(['fit', 'holdout'])).any():
        raise ValueError('Reference cells require ref_split fit or holdout')
    if (ref & meta.ref_group.str.strip().isin(['', 'NA'])).any():
        raise ValueError('Reference cells require a nonempty ref_group')
    if not np.array_equal(meta.cell_id.to_numpy(), np.asarray(cell_ids, dtype=str)):
        raise ValueError('cells table must list cells in exact count-matrix order')
    if 'dataset' not in meta.columns:
        meta['dataset'] = 'dataset'
    return meta


def read_positions(path, *, genome_build=None, coordinate_system=None, gene_namespace=None,
                   centromeres=None) -> pd.DataFrame:
    validate_input_contract('raw', genome_build, coordinate_system, gene_namespace, centromeres)
    pos = pd.read_csv(path, sep='\t', keep_default_na=False, dtype={'gene': str, 'chromosome': str})
    for c in ['gene', 'chromosome', 'start', 'end']:
        if c not in pos.columns:
            raise ValueError(f'positions table lacks {c}')
    pos['gene'] = validate_feature_ids(pos.gene, gene_namespace)
    if (pos.chromosome.astype(str).str.strip().eq('').any()
            or pos.chromosome.astype(str).str.lower().isin(['nan', 'none', 'null', 'na']).any()):
        raise ValueError('positions table has missing chromosome')
    pos['chromosome'] = pos.chromosome.astype(str).str.strip()
    pos.loc[~pos.chromosome.str.startswith('chr'), 'chromosome'] = 'chr' + pos.chromosome
    for c in ['start', 'end']:
        values = pd.to_numeric(pos[c], errors='coerce').to_numpy(dtype=float)
        if not np.isfinite(values).all() or np.any(values != np.floor(values)):
            raise ValueError(f'positions {c} must be finite integer coordinates')
        if np.any(values < 0) or np.any(values >= np.iinfo(np.int64).max):
            raise ValueError(f'positions {c} is outside supported coordinate bounds')
        pos[c] = values.astype(np.int64)
    if coordinate_system == 'one-based-closed':
        if (pos.start < 1).any() or (pos.end < pos.start).any():
            raise ValueError('One-based closed intervals require 1 <= start <= end')
        pos['start'] -= 1
    if (pos.end <= pos.start).any():
        raise ValueError('Zero-based half-open intervals require end > start')
    pos.attrs.update(genome_build=genome_build, coordinate_system='zero-based-half-open', gene_namespace=gene_namespace)
    return pos


def read_site_list(path):
    """Genomic sites, one per line as chrN:pos (or N:pos), normalized to chrN:pos."""
    if path is None:
        return set()
    sites = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str)[0].str.strip()
    return {x if x.startswith('chr') else 'chr' + x for x in sites}


def read_gene_list(path, gene_namespace='symbol'):
    if path is None:
        return set()
    names = pd.read_csv(path, sep='\t', header=None, comment='#', dtype=str, keep_default_na=False)[0]
    names = names.str.strip()
    if gene_namespace == 'symbol':
        names = names.str.upper()
    return set(validate_feature_ids(names.drop_duplicates(), gene_namespace))
