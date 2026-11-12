"""Count-to-profile identity, conservation and scalar normalization references."""
import math

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from haplodose.stacking import stacked_count_profile, stacked_rna_profile


def fixture():
    counts = np.array([[4, 1, 0, 2], [2, 3, 1, 0], [1, 2, 3, 4], [9, 7, 6, 8], [1, 0, 0, 0]])
    cells, genes = ['q1', 'q2', 'r1', 'r2', 'unused'], ['a', 'b', 'c', 'd']
    features = pd.DataFrame(dict(gene=['a', 'b', 'c', 'd'], bin=[0, 0, 1, 1]))
    return counts, cells, genes, features


def run(counts, cells, genes, features, **kwargs):
    return stacked_count_profile(counts, cells, genes, features, ['p', 'p'],
        ['q1', 'q2'], ['r1', 'r2'], counts_kind=kwargs.pop('counts_kind', 'raw'), **kwargs)


def test_scalar_reference_and_sparse_dense_conservation():
    counts, cells, genes, features = fixture()
    profile, effects, record = run(counts, cells, genes, features)
    y, r = [6, 4, 1, 2], [10, 9, 9, 12]
    scalar = np.array([math.log((a+.5)/(b/sum(r)*sum(y)+.5)) for a, b in zip(y, r)])
    scalar -= np.median(scalar)
    # Scalar libm and vector log may differ by a few rounded operations.
    budget = 64*len(scalar)*np.finfo(float).eps*max(1., max(abs(scalar)))
    assert max(abs(effects-scalar)) <= budget
    native, _ = stacked_rna_profile(effects, features.bin.to_numpy(), ['p', 'p'])
    np.testing.assert_array_equal(profile, native)
    sp, se, sm = run(sparse.csr_matrix(counts), cells, genes, features)
    np.testing.assert_array_equal(sp, profile)
    np.testing.assert_array_equal(se, effects)
    assert sm == record
    assert record['query_feature_total'] == 13 and record['reference_feature_total'] == 40
    assert record['calibration_status'] == 'unvalidated'


def test_count_entry_point_forwards_dual_solver():
    counts, cells, genes, features = fixture()
    profile, effects, record = run(counts, cells, genes, features, solver='dual')
    native, diagnostics = stacked_rna_profile(effects, features.bin.to_numpy(), ['p', 'p'], solver='dual')
    np.testing.assert_array_equal(profile, native)
    assert record['stacking'] == diagnostics['stacking']
    assert record['stacking']['solver'] == 'dual'
    with pytest.raises(ValueError, match='solver'):
        run(counts, cells, genes, features, solver='unknown')


def test_axis_permutations_preserve_canonical_gene_folds_and_ignore_unused_cells():
    counts, cells, genes, features = fixture()
    expected = run(counts, cells, genes, features)
    rows, cols = np.array([3, 0, 4, 1, 2]), np.array([3, 0, 2, 1])
    changed = counts.copy()
    changed[-1] += 100000
    actual = run(changed[rows][:, cols], np.array(cells)[rows], np.array(genes)[cols], features)
    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])
    assert actual[2] == expected[2]


@pytest.mark.parametrize('problem', ['missing_cell', 'missing_gene', 'duplicate_cell', 'duplicate_gene',
    'duplicate_feature', 'fractional_bins', 'overlap', 'fractional_counts', 'negative_counts',
    'nonfinite_counts', 'zero_reference', 'ambiguous_scale', 'inexact_total'])
def test_ambiguous_or_invalid_measurements_reject(problem):
    counts, cells, genes, features = fixture()
    counts = counts.astype(float)
    query, reference, kind = ['q1', 'q2'], ['r1', 'r2'], 'raw'
    if problem == 'missing_cell':query = ['absent']
    if problem == 'missing_gene':features.loc[0, 'gene'] = 'absent'
    if problem == 'duplicate_cell':cells[1] = cells[0]
    if problem == 'duplicate_gene':genes[1] = genes[0]
    if problem == 'duplicate_feature':features.loc[1, 'gene'] = features.loc[0, 'gene']
    if problem == 'fractional_bins':features['bin'] = features.bin.astype(float)
    if problem == 'overlap':reference = ['q1']
    if problem == 'fractional_counts':counts[0, 0] = .25
    if problem == 'negative_counts':counts[0, 0] = -1
    if problem == 'nonfinite_counts':counts[0, 0] = np.nan
    if problem == 'zero_reference':counts[2:4] = 0
    if problem == 'ambiguous_scale':kind = 'normalized'
    if problem == 'inexact_total':counts[0, 0] = 2**53-1
    with pytest.raises(ValueError):
        stacked_count_profile(counts, cells, genes, features, ['p', 'p'], query, reference, counts_kind=kind)


def test_reconstructed_label_does_not_gain_count_calibration():
    counts, cells, genes, features = fixture()
    p, e, metadata = run(counts, cells, genes, features, counts_kind='reconstructed')
    assert metadata['counts_kind'] == 'reconstructed'
    assert metadata['input_scale_interpretation'] == 'descriptive_only_no_count_model_calibration'
    assert metadata['calibration_status'] == 'unvalidated'


def test_excluded_observations_do_not_enter_conversion_or_count_validation():
    counts, cells, genes, features = fixture()
    expected = run(counts, cells, genes, features)
    matrix = counts.astype(float)
    matrix[-1] = [np.nan, np.inf, -.5, .25]
    actual = run(sparse.csr_matrix(matrix), cells, genes, features)
    np.testing.assert_array_equal(actual[0], expected[0])
    np.testing.assert_array_equal(actual[1], expected[1])
    assert actual[2] == expected[2]
    assert actual[2]['count_validation_scope'] == 'selected_cells_and_features'
