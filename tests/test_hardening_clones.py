"""Scientific safety boundaries and an independently evaluated allele observation model."""
import numpy as np
import pandas as pd
import pytest
from scipy import sparse, stats

from haplodose.alleles import AlleleParams, SpecimenAlleles, cell_loglik, imbalance_runs, segment_alleles
from haplodose.clones import CloneParams, Rows, SpecimenModel


def test_query_derived_normal_reference_requires_opt_in():
    assert not CloneParams().in_sample_reference


def test_joint_channel_null_is_rejected_before_fitting():
    Y = np.ones((60, 40), np.float32)
    reset = np.zeros(40, bool)
    reset[0] = True
    with pytest.raises(ValueError, match='does not reproduce allele proposals'):
        SpecimenModel(Rows(Y, np.arange(60)), Y.sum(axis=1), np.full((40, 1), 1 / 40), reset,
                      np.full(40, .1), CloneParams(split_null_draws=19), alleles=object())


def test_sparse_beta_binomial_matches_scipy_probability_contrast():
    A = sparse.csr_matrix([[0., 2., 0.], [1., 0., 8.], [0., 0., 0.]])
    D = sparse.csr_matrix([[0., 5., 1.], [2., 0., 10.], [0., 0., 0.]])
    B = np.array([[.5, .5, .5], [.8, .2, .7]])
    concentration, shrink = 7., .9
    got = cell_loglik(A, D, B, shrink, concentration)
    b = .5 + shrink * (B - .5)
    exact = np.array([[stats.betabinom.logpmf(A.toarray()[c], D.toarray()[c],
                                             concentration * bk, concentration * (1 - bk)).sum()
                       for bk in b] for c in range(A.shape[0])])
    np.testing.assert_allclose(got[:, 1] - got[:, 0], exact[:, 1] - exact[:, 0], atol=2e-6)
    assert np.array_equal(got[2], [0., 0.])


def test_beta_binomial_binomial_limit_and_correlated_depth():
    B = np.array([[.5], [.8]])
    A, D = sparse.csr_matrix([[8.]]), sparse.csr_matrix([[10.]])
    exact = cell_loglik(A, D, B, kappa=1e6)
    legacy = cell_loglik(A, D, B)
    np.testing.assert_allclose(exact[:, 1] - exact[:, 0], legacy[:, 1] - legacy[:, 0], atol=1e-4)
    # Repeating the same imbalanced fraction at greater depth is correlated evidence under finite kappa.
    deep = cell_loglik(1000 * A, 1000 * D, B, kappa=10.)
    deep_binomial = cell_loglik(1000 * A, 1000 * D, B)
    assert abs(deep[0, 1] - deep[0, 0]) < abs(deep_binomial[0, 1] - deep_binomial[0, 0]) / 100


@pytest.mark.parametrize('kappa', [0., -1., np.inf, np.nan])
def test_invalid_cell_concentration_is_rejected(kappa):
    with pytest.raises(ValueError, match='finite and positive'):
        cell_loglik(sparse.csr_matrix([[1.]]), sparse.csr_matrix([[2.]]), np.array([[.5]]), kappa=kappa)


def test_beta_binomial_does_not_ignore_orphan_alt_counts():
    with pytest.raises(ValueError, match='0 <= A <= D'):
        cell_loglik(sparse.csr_matrix([[1., 0.]]), sparse.csr_matrix([[0., 2.]]),
                    np.array([[.5, .5]]), kappa=10.)


def test_one_based_snps_map_once_at_half_open_bin_boundary():
    snps = pd.DataFrame(dict(chromosome='chr1', pos=[100, 101, 200, 201, 300, 301], cM=np.arange(6) * .01))
    bins = pd.DataFrame(dict(chromosome='chr1', start=[100, 200], end=[200, 300]))
    alleles = SpecimenAlleles(np.zeros((1, 6)), np.ones((1, 6)), snps, 50.)
    np.testing.assert_array_equal(alleles.snps_in(bins, 0, 0), [1, 2])
    np.testing.assert_array_equal(alleles.snps_in(bins, 1, 1), [3, 4])
    post = np.zeros((6, 9))
    post[:, 0] = 1.
    for k in range(2):
        summary = segment_alleles(post, snps, np.ones(6), bins, k, k)
        assert summary['n_snps'] == 2
    post[:, 0], post[:, 1] = 0., 1.
    runs = imbalance_runs(post, snps, np.ones(6) * 10, AlleleParams(min_snps=1))
    assert runs[0]['start'] == 99 and runs[0]['end'] == 301
