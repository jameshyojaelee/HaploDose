"""Rust aggregation must equal the SciPy reference for any chunking."""
import numpy as np
import pytest
from scipy import sparse

from haplodose._native import bin_sums, grouped_column_sums, native_available


def random_counts(seed, n=257, g=301):
    rng = np.random.default_rng(seed)
    X = sparse.random(n, g, density=0.07, random_state=rng, format='csr')
    X.data = rng.poisson(5, X.nnz).astype(np.float32) + 1
    return X


@pytest.mark.skipif(not native_available(), reason='Rust extension not built')
def test_native_matches_scipy():
    X = random_counts(1)
    gene_bin = np.random.default_rng(2).integers(-1, 17, X.shape[1]).astype(np.int32)
    a, ta = bin_sums(X, gene_bin, 17, native=True)
    b, tb = bin_sums(X, gene_bin, 17, native=False)
    np.testing.assert_array_equal(a, b)
    np.testing.assert_array_equal(ta, tb)
    group = np.random.default_rng(3).integers(-1, 4, X.shape[0]).astype(np.int32)
    np.testing.assert_array_equal(grouped_column_sums(X, group, 4, native=True),
                                  grouped_column_sums(X, group, 4, native=False))


@pytest.mark.parametrize('size', [1, 7, 64, 1000])
def test_chunked_bin_sums_equal_whole(size):
    X = random_counts(4)
    gene_bin = (np.arange(X.shape[1]) // 10).astype(np.int32)
    whole, total = bin_sums(X, gene_bin, gene_bin.max() + 1)
    parts = [bin_sums(X[a:a + size], gene_bin, gene_bin.max() + 1) for a in range(0, X.shape[0], size)]
    np.testing.assert_array_equal(whole, np.vstack([p[0] for p in parts]))
    np.testing.assert_array_equal(total, np.concatenate([p[1] for p in parts]))


@pytest.mark.skipif(not native_available(), reason='Rust extension not built')
def test_native_rejects_bad_bins():
    X = random_counts(5, 4, 5)
    with pytest.raises(ValueError):
        bin_sums(X, np.array([0, 1, 9, 0, 0], np.int32), 2)


@pytest.mark.parametrize('native', [False, True])
def test_bin_sums_reject_unrepresentable_aggregate(native):
    if native and not native_available():
        pytest.skip('Rust extension not built')
    X = sparse.csr_matrix(np.array([[2 ** 24, 1]], dtype=np.float64))
    with pytest.raises(ValueError, match='bin sums'):
        bin_sums(X, np.array([0, 0]), 1, native=native)


@pytest.mark.parametrize('native', [False, True])
def test_bin_sums_accumulate_before_casting(native):
    if native and not native_available():
        pytest.skip('Rust extension not built')
    # The intermediate float32 sum 2**24 + 1 would round down; the final
    # sum 2**24 + 2 is representable and must survive aggregation exactly.
    X = sparse.csr_matrix(np.array([[2 ** 24, 1, 1]], dtype=np.float64))
    bins, total = bin_sums(X, np.array([0, 0, 0]), 1, native=native)
    assert float(bins[0, 0]) == 2 ** 24 + 2
    assert total[0] == 2 ** 24 + 2
    X = sparse.csr_matrix(np.array([[2 ** 24, 2]], dtype=np.float64))
    bins, _ = bin_sums(X, np.array([0, 0]), 1, native=native)
    assert float(bins[0, 0]) == 2 ** 24 + 2


@pytest.mark.parametrize('native', [False, True])
@pytest.mark.parametrize('mapping', [np.array([-2, 0]), np.array([0, 2]),
                                     np.array([0]), np.array([0., 0.]), np.array([0, 2 ** 32])])
def test_bin_identity_validation_matches_backends(native, mapping):
    if native and not native_available():
        pytest.skip('Rust extension not built')
    with pytest.raises(ValueError):
        bin_sums(sparse.csr_matrix([[1, 2]]), mapping, 1, native=native)


@pytest.mark.parametrize('native', [False, True])
@pytest.mark.parametrize('groups', [np.array([-2, 0]), np.array([0, 2]), np.array([0]), np.array([0., 0.])])
def test_row_group_identity_validation_matches_backends(native, groups):
    if native and not native_available():
        pytest.skip('Rust extension not built')
    with pytest.raises(ValueError):
        grouped_column_sums(sparse.csr_matrix([[1], [2]]), groups, 1, native=native)


@pytest.mark.parametrize('native', [False, True])
def test_input_precision_is_checked_before_float32_conversion(native):
    if native and not native_available():
        pytest.skip('Rust extension not built')
    X = sparse.csr_matrix(np.array([[2 ** 24 + 1]], dtype=np.int64))
    with pytest.raises(ValueError, match='input counts'):
        bin_sums(X, np.array([0]), 1, native=native)
