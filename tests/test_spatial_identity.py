"""A graph must retain cell identity under permutations, including distance ties."""
import numpy as np
import pandas as pd
import pytest

from haplodose.clones import neighbour_graph
from haplodose.pipeline import _specimen_graph


@pytest.mark.parametrize('kind', ['duplicates', 'grid'])
@pytest.mark.parametrize('k', [1, 2, 4, 8])
def test_spatial_graph_ties_preserve_cell_identity(kind, k):
    xy = (np.repeat(np.arange(20, dtype=float)[:, None] * [10., 0.], 3, axis=0)
          if kind == 'duplicates' else np.column_stack(np.unravel_index(np.arange(64), (8, 8))).astype(float))
    ids = np.array([f'cell_{i:03d}' for i in range(len(xy))])
    permutation = np.random.default_rng(20261109).permutation(len(xy))
    inverse = np.argsort(permutation)
    graph = neighbour_graph(xy, k, cell_ids=ids)
    permuted = neighbour_graph(xy[permutation], k, cell_ids=ids[permutation])
    assert (graph != permuted[inverse][:, inverse]).nnz == 0
    assert np.all(graph.diagonal() == 0)
    assert np.all(np.diff(graph.indptr) == k)
    np.testing.assert_allclose(np.asarray(graph.sum(axis=1)).ravel(), 1.)
    distances = ((xy[:, None] - xy[None]) ** 2).sum(axis=-1)
    np.fill_diagonal(distances, np.inf)
    cutoff = np.sort(distances, axis=1)[:, k - 1]
    for i in range(len(xy)):
        selected = graph.indices[graph.indptr[i]:graph.indptr[i + 1]]
        assert np.all(distances[i, selected] <= cutoff[i])


@pytest.mark.parametrize('ids', [['a', 'a'], ['a'], ['a', ' '], ['a', None], [['a'], ['b']]])
def test_spatial_graph_rejects_ambiguous_cell_ids(ids):
    with pytest.raises(ValueError, match='spatial cell IDs'):
        neighbour_graph(np.array([[0., 0.], [1., 0.]]), 1, cell_ids=ids)


def test_identity_order_does_not_change_untied_geometry():
    xy = np.column_stack(([0., 1., 4., 13., 35., 80.], np.zeros(6)))
    ids = np.array(['010', '001', 'NA', 'z', 'b', 'a'])
    before = neighbour_graph(xy, 2)
    after = neighbour_graph(xy, 2, cell_ids=ids)
    assert (before != after).nnz == 0


def test_spatial_graph_separates_overlapping_section_coordinates():
    meta = pd.DataFrame(dict(cell_id=list('abcde'), section=['A', 'B', 'A', 'B', 'C'],
                             x=[0, 0, 2, 2, 0], y=[0, 0, 0, 0, 0]))
    expected = np.zeros((5, 5))
    expected[[0, 2, 1, 3], [2, 0, 3, 1]] = 1
    graph, status = _specimen_graph(meta, 8)
    assert status == 'section_graph'
    np.testing.assert_array_equal(graph.toarray(), expected)
    order = [3, 0, 4, 2, 1]
    permuted, _ = _specimen_graph(meta.iloc[order], 8)
    np.testing.assert_array_equal(permuted.toarray(), expected[np.ix_(order, order)])


@pytest.mark.parametrize('columns,status', [
    ({'x': [0, 1]}, 'incomplete_coordinate_columns'),
    ({'x': [0, 1], 'y': ['', 1]}, 'nonfinite_coordinates'),
    ({'x': [0, 1], 'y': [0, 1], 'section': ['A', '']}, 'missing_section_ids'),
    ({'x': [0, 1], 'y': [0, 1], 'section': ['A', 'B']}, 'insufficient_spots'),
])
def test_spatial_graph_reports_fallback_reason(columns, status):
    meta = pd.DataFrame(dict(cell_id=['a', 'b'], **columns))
    assert _specimen_graph(meta, 1) == (None, status)
