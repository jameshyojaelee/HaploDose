"""Neighbour identity invariants, with an independent small distance reference."""
import numpy as np
import igraph as ig
from scipy.spatial.distance import cdist
from haplodose.clones import CloneParams, leiden_groups


def capture_edges(monkeypatch):
    original = ig.Graph
    records = []

    def graph(*args, **kwargs):
        records.append(np.asarray(kwargs['edges']).copy())
        return original(*args, **kwargs)

    monkeypatch.setattr(ig, 'Graph', graph)
    return records


def test_duplicate_feature_points_retain_all_nonself_neighbours(monkeypatch):
    # Two triples at exactly equal locations: each point's two closest other
    # vertices are independently known. Skipping an arbitrary first zero-distance
    # neighbour can retain self and lose a real neighbour after graph simplification.
    features = np.repeat([[0., 0.], [1., 1.]], 3, axis=0).astype(np.float32)
    records = capture_edges(monkeypatch)
    leiden_groups(features, CloneParams(min_cells=1, n_neighbors=2, n_pcs=1))
    assert len(records) == 1
    for i in range(6):
        targets = records[0][records[0][:, 0] == i, 1]
        assert set(targets) == set(range(3 * (i // 3), 3 * (i // 3) + 3)) - {i}


def test_neighbours_agree_with_independent_rank_one_geometry(monkeypatch):
    # Centred rank-one geometry has analytically known distances regardless of
    # PCA sign: multiplying by [1,2,3] only scales every squared distance by14.
    coordinate = np.array([-8., -3., -1., 0., 4., 10., 18.])
    features = coordinate[:, None] * np.array([[1., 2., 3.]])
    records = capture_edges(monkeypatch)
    leiden_groups(features.astype(np.float32), CloneParams(min_cells=1, n_neighbors=2, n_pcs=1))
    distances = cdist(coordinate[:, None], coordinate[:, None], metric='sqeuclidean')
    np.fill_diagonal(distances, np.inf)
    expected = np.argsort(distances, axis=1)[:, :2]
    for i in range(len(coordinate)):
        assert set(records[0][records[0][:, 0] == i, 1]) == set(expected[i])
