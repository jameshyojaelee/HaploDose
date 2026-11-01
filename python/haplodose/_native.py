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
    return (np.ascontiguousarray(X.indptr, dtype=np.int64),
            np.ascontiguousarray(X.indices, dtype=np.int32),
            np.ascontiguousarray(X.data, dtype=np.float32))


def bin_sums(X, gene_bin, n_bins, native=True):
    """Cells x bins sums over genes with gene_bin >= 0, plus each row's total."""
    gene_bin = np.ascontiguousarray(gene_bin, dtype=np.int32)
    if native and _core is not None:
        return _core.bin_sums(*_csr(X), gene_bin, int(n_bins))
    X = sparse.csr_matrix(X, dtype=np.float32)
    keep = np.flatnonzero(gene_bin >= 0)
    W = sparse.csr_matrix((np.ones(len(keep), np.float32), (keep, gene_bin[keep])),
                          shape=(len(gene_bin), int(n_bins)))
    out = np.asarray((X @ W).toarray(), dtype=np.float32)
    totals = np.asarray(X.astype(np.float64).sum(axis=1)).ravel()
    return out, totals


def grouped_column_sums(X, row_group, n_groups, native=True):
    """Groups x genes column sums; rows with group -1 are skipped."""
    row_group = np.ascontiguousarray(row_group, dtype=np.int32)
    X = sparse.csr_matrix(X)
    if native and _core is not None:
        return _core.grouped_column_sums(*_csr(X), row_group, int(n_groups), X.shape[1])
    keep = np.flatnonzero(row_group >= 0)
    G = sparse.csr_matrix((np.ones(len(keep)), (row_group[keep], keep)),
                          shape=(int(n_groups), X.shape[0]))
    return np.asarray((G @ X.astype(np.float64)).toarray())
