"""Independent small fixtures for allele-table source reconciliation."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy import sparse


@pytest.fixture
def trace():
    path = Path(__file__).resolve().parents[1] / 'benchmarks/hardening/allele_trace.py'
    spec = importlib.util.spec_from_file_location('independent_allele_trace', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sources():
    sites = pd.DataFrame(dict(snp_id=['1_10_A_G', '1_20_C_T', '1_30_G_A', '1_40_A_C'],
                              DP_all=[3, 1, 4, 5], OTH_all=[0, 0, 1, 0]))
    phase = pd.DataFrame(dict(snp_id=sites.snp_id, GT=['0|1', '1|0', '0|1', '0|0']))
    counts = dict(DP=sparse.csr_matrix([[2, 1], [1, 0], [2, 2], [3, 2]]),
                  AD=sparse.csr_matrix([[1, 0], [1, 0], [1, 1], [2, 1]]))
    return sites, phase, counts


def test_expected_final_reconciles_filters_without_using_final_table(trace):
    sites, phased, counts = sources()
    expected, reconciliation = trace.expected_final('channel', sites, counts, ['C1', 'C2'], phased)
    assert expected[['cell', 'snp_id', 'AD', 'DP', 'GT']].to_dict('records') == [
        dict(cell='C1', snp_id='1_10_A_G', AD=1, DP=2, GT='0|1'),
        dict(cell='C2', snp_id='1_10_A_G', AD=0, DP=1, GT='0|1')]
    assert reconciliation == dict(stored_DP_cell_sites=7, nonzero_DP_cell_sites=7,
                                  explicit_zero_DP_cell_sites=0, excluded_pooled_DP_or_OTH=3,
                                  excluded_no_phased_het=2, expected_final_rows=2)


def test_final_count_mismatch_is_found_despite_matching_global_totals(trace):
    phased = pd.DataFrame(dict(snp_id=['1_10_A_G'], CHROM=['1'], POS=[10], REF=['A'], ALT=['G'], GT=['0|1']))
    expected = pd.DataFrame(dict(cell=['C1', 'C2'], snp_id=['1_10_A_G'] * 2,
                                 AD=[1, 0], DP=[2, 1], GT=['0|1'] * 2, channel=['Pool7_16'] * 2))
    # Swapping cell-specific counts preserves every global sum but violates identity.
    final = pd.DataFrame(dict(cell=['C1', 'C2'], snp_id=['1_10_A_G'] * 2, CHROM=['1'] * 2,
                              POS=[10] * 2, REF=['A'] * 2, ALT=['G'] * 2,
                              AD=[0, 1], DP=[1, 2], GT=['0|1'] * 2))
    mismatches = []
    result = trace.compare_final(final, expected, phased, mismatches)
    assert result['count_or_GT_mismatches'] == 2
    assert len(mismatches) == 2


def test_final_orientation_mismatch_is_not_hidden_by_equal_counts(trace):
    phased = pd.DataFrame(dict(snp_id=['1_10_A_G'], CHROM=['1'], POS=[10], REF=['A'], ALT=['G'], GT=['1|0']))
    expected = pd.DataFrame(dict(cell=['C'], snp_id=['1_10_A_G'], AD=[1], DP=[3], GT=['1|0'], channel=['Pool7_16']))
    final = pd.DataFrame(dict(cell=['C'], snp_id=['1_10_A_G'], CHROM=['1'], POS=[10], REF=['A'], ALT=['G'],
                              AD=[1], DP=[3], GT=['0|1']))
    mismatches = []
    result = trace.compare_final(final, expected, phased, mismatches)
    assert result['count_or_GT_mismatches'] == 1 and mismatches


def test_matrix_parser_does_not_round_fractional_counts(tmp_path, trace):
    path = tmp_path / 'fractional.mtx'
    path.write_text('%%MatrixMarket matrix coordinate real general\n1 1 1\n1 1 0.9999999\n')
    with pytest.raises(ValueError, match='Malformed'):
        trace.matrix(path)
