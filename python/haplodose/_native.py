"""Rust kernels with NumPy/SciPy fallbacks; the fallbacks define the expected results."""
from __future__ import annotations

import numpy as np
from scipy import sparse

try:
    from . import _core
except ImportError:  # pure-Python install or unbuilt checkout
    _core = None


def native_available() -> bool:
    return _core is not None


def _csr(X):
    X = sparse.csr_matrix(X)
    X.check_format(full_check=True)
    data = X.data
    if (not np.isfinite(data).all() or (data < 0).any()
            or not np.equal(data, np.floor(data)).all()):
        raise ValueError('counts must be finite nonnegative integers')
    # Integer addition in float64 is exact in this bounded range. Conservatively
    # reject chunks beyond it rather than promise exactness after double rounding.
    if data.size and int(data.max()) * len(data) > 2 ** 53:
        raise ValueError('count chunk exceeds supported exact float64 accumulation range')
    if (data.size and np.any(data.astype(np.float32).astype(np.float64) != data.astype(np.float64))):
        raise ValueError('input counts cannot be represented exactly in float32')
    if X.shape[1] > np.iinfo(np.int32).max:
        raise ValueError('column count exceeds supported int32 identity range')
    return (np.ascontiguousarray(X.indptr, dtype=np.int64),
            np.ascontiguousarray(X.indices, dtype=np.int32),
            np.ascontiguousarray(X.data, dtype=np.float32))


def _groups(values, length, count, name):
    if not isinstance(count, (int, np.integer)) or count < 0 or count > np.iinfo(np.int32).max:
        raise ValueError(f'{name} count must be a nonnegative supported integer')
    values = np.asarray(values)
    if (values.shape != (length,) or not np.issubdtype(values.dtype, np.integer)
            or np.any((values < -1) | (values >= count))):
        raise ValueError(f'{name} must have the aligned length and contain -1 or valid indices')
    return np.ascontiguousarray(values, dtype=np.int32)


def bin_sums(X, gene_bin, n_bins, native=True):
    """Exact integer bin sums in float32, or fail; row totals are float64."""
    X = sparse.csr_matrix(X)
    gene_bin = _groups(gene_bin, X.shape[1], n_bins, 'gene_bin')
    csr = _csr(X)
    if native and _core is not None:
        return _core.bin_sums(*csr, gene_bin, int(n_bins))
    X = sparse.csr_matrix(X, dtype=np.float64)
    keep = np.flatnonzero(gene_bin >= 0)
    W = sparse.csr_matrix((np.ones(len(keep), np.float64), (keep, gene_bin[keep])),
                          shape=(len(gene_bin), int(n_bins)))
    exact = np.asarray((X @ W).toarray())
    out = exact.astype(np.float32)
    if not np.isfinite(out).all() or np.any(out.astype(np.float64) != exact):
        raise ValueError('bin sums cannot be represented exactly in float32')
    totals = np.asarray(X.sum(axis=1)).ravel()
    return out, totals


def grouped_column_sums(X, row_group, n_groups, native=True):
    """Groups x genes column sums; rows with group -1 are skipped."""
    X = sparse.csr_matrix(X)
    row_group = _groups(row_group, X.shape[0], n_groups, 'row_group')
    csr = _csr(X)
    if native and _core is not None:
        return _core.grouped_column_sums(*csr, row_group, int(n_groups), X.shape[1])
    keep = np.flatnonzero(row_group >= 0)
    G = sparse.csr_matrix((np.ones(len(keep)), (row_group[keep], keep)),
                          shape=(int(n_groups), X.shape[0]))
    return np.asarray((G @ X.astype(np.float64)).toarray())
