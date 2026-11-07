"""Malformed allele measurements cannot become truncated or overflowing counts."""
import numpy as np
import pandas as pd
import pytest

from haplodose.alleles import read_alleles


@pytest.mark.parametrize('build,message', [(None, 'explicitly declared'), ('hg19', 'must match')])
def test_joint_inference_requires_matching_declared_assemblies(tmp_path, build, message):
    from haplodose.pipeline import run
    with pytest.raises(ValueError, match=message):
        run('unread-counts', 'unread-cells', 'unread-positions', tmp_path / 'out',
            alleles='unread-alleles', allele_genome_build=build, counts_kind='raw',
            genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')


@pytest.mark.parametrize('column,value', [('AD', np.inf), ('DP', np.inf), ('AD', np.nan), ('DP', np.nan),
                                         ('AD', -1), ('DP', -1), ('AD', 0.9999999), ('DP', 1.9999999),
                                         ('AD', 2 ** 24 + 1), ('DP', 2 ** 31), ('POS', np.inf),
                                         ('POS', np.nan), ('POS', 100.9999999), ('POS', 0),
                                         ('POS', 2 ** 31)])
def test_malformed_allele_measurement_is_rejected_before_cast(tmp_path, column, value):
    row = dict(cell='c', CHROM='1', POS=101, cM=.1, AD=1, DP=2, GT='0|1', gene='G')
    row[column] = value
    path = tmp_path / 'alleles.tsv'
    pd.DataFrame([row]).to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match=column):
        read_alleles(path, ['c'])


def test_alt_count_cannot_exceed_total_count(tmp_path):
    path = tmp_path / 'alleles.tsv'
    pd.DataFrame([dict(cell='c', CHROM='1', POS=101, cM=.1, AD=3, DP=2, GT='0|1', gene='G')]).to_csv(
        path, sep='\t', index=False)
    with pytest.raises(ValueError, match='AD cannot exceed DP'):
        read_alleles(path, ['c'])


def test_duplicate_cell_site_counts_are_summed_exactly_after_orientation(tmp_path):
    path = tmp_path / 'alleles.tsv'
    pd.DataFrame([dict(cell='c', CHROM='1', POS=101, cM=.1, AD=1, DP=3, GT='0|1', gene='G'),
                  dict(cell='c', CHROM='1', POS=101, cM=.1, AD=1, DP=4, GT='1|0', gene='G')]).to_csv(
        path, sep='\t', index=False)
    data = read_alleles(path, ['c'], chunksize=1)
    assert data.A[0, 0] == 4 and data.D[0, 0] == 7
    assert data.snps.pos.tolist() == [101]


def test_duplicate_depth_cannot_overflow_exact_count_channel(tmp_path):
    path = tmp_path / 'alleles.tsv'
    row = dict(cell='c', CHROM='1', POS=101, cM=.1, AD=1, DP=2 ** 24, GT='0|1', gene='G')
    pd.DataFrame([row, row]).to_csv(path, sep='\t', index=False)
    with pytest.raises(ValueError, match='Aggregated allele counts'):
        read_alleles(path, ['c'], chunksize=1)


def test_valid_unknown_cells_remain_explicitly_excluded(tmp_path):
    path = tmp_path / 'alleles.tsv'
    pd.DataFrame([dict(cell='c', CHROM='1', POS=101, cM=.1, AD=1, DP=2, GT='0|1', gene='G'),
                  dict(cell='other', CHROM='1', POS=101, cM=.1, AD=1, DP=2, GT='0|1', gene='G')]).to_csv(
        path, sep='\t', index=False)
    data = read_alleles(path, ['c'])
    assert data.excluded['unknown_cell'] == 1 and data.D[0, 0] == 2


def test_integer_valued_float_position_obeys_site_mask(tmp_path):
    path = tmp_path / 'alleles.tsv'
    pd.DataFrame([dict(cell='c', CHROM='1', POS=101., cM=.1, AD=1, DP=2, GT='0|1', gene='G'),
                  dict(cell='c', CHROM='1', POS=102., cM=.1, AD=1, DP=2, GT='0|1', gene='G')]).to_csv(
        path, sep='\t', index=False)
    data = read_alleles(path, ['c'], exclude_sites={'chr1:101'})
    assert data.excluded['site'] == 1 and data.snps.pos.tolist() == [102]
