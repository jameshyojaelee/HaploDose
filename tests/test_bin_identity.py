"""Named observations retain bin membership across equivalent matrix layouts."""
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from haplodose.bins import chromosome_starts, make_bins


@pytest.mark.parametrize('order', [list(range(6)), [0,1,3,2,4,5], list(reversed(range(6))),
                                  [3,0,5,2,1,4], [2,3,4,5,0,1]])
def test_tied_start_membership_and_raw_aggregation_follow_gene_identity(order):
    names = np.array(list('ABCDEF')); order = np.array(order)
    positions = pd.DataFrame(dict(gene=names, chromosome='chr1', start=[1,2,3,3,4,5],
                                  end=[2,3,9,4,5,6]))
    counts = np.array([[2,3,5,101,7,11],[13,17,19,23,29,31]])
    bins, mapping, used = make_bins(names[order], positions.iloc[::-1], np.ones(6),
                                   bin_genes=3, exclude_pattern=None)
    assert dict(zip(names[order], mapping)) == dict(A=0,B=0,C=0,D=1,E=1,F=1)
    assert used.gene.tolist() == list('ABCDEF')
    assert bins.first_gene.tolist() == ['A','D'] and bins.last_gene.tolist() == ['C','F']
    onehot = sparse.csr_matrix((np.ones(6, dtype=int), (np.arange(6), mapping)), shape=(6,2))
    actual = (sparse.csr_matrix(counts[:, order]) @ onehot).toarray()
    assert np.array_equal(actual, [[10,119],[49,83]])
    assert np.array_equal(actual.sum(axis=1), counts.sum(axis=1))


def test_gene_means_remain_identity_aligned_when_tied_features_are_filtered():
    names = np.array(list('ABCDEFG')); order = np.array([3,6,2,1,4,0,5])
    positions = pd.DataFrame(dict(gene=names, chromosome='chr1', start=[1,2,3,3,3,4,5],
                                  end=[2,3,4,4,4,5,6]))
    means = np.array([1,1,1,0,1,1,1])
    _, mapping, used = make_bins(names[order], positions, means[order], bin_genes=3, exclude_pattern=None)
    assert dict(zip(names[order], mapping)) == dict(A=0,B=0,C=0,D=-1,E=1,F=1,G=1)
    assert used.gene.tolist() == list('ABCEFG')


def test_tied_gene_identity_keeps_explicit_gap_resets():
    names = np.array(list('ABCDEFGHIJKL')); order = np.arange(12)[::-1]
    starts = np.array([1,2,3,3,4,5,100,101,102,102,103,104])
    positions = pd.DataFrame(dict(gene=names, chromosome='chr1', start=starts, end=starts+1))
    bins, mapping, used = make_bins(names[order], positions, np.ones(12), bin_genes=3,
                                   max_gap_bp=10, exclude_pattern=None)
    assert dict(zip(names[order], mapping)) == dict(zip(names, np.repeat(np.arange(4),3)))
    assert used.gene.tolist() == names.tolist()
    assert chromosome_starts(bins).tolist() == [True,False,True,False]
