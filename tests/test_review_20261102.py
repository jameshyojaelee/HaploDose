"""Regression tests for review 2026-11-02 findings 1 (shared-dose boundaries) and 7 (spatial graph self-edges)."""
import numpy as np

from haplodose.clones import CloneParams, Rows, SpecimenModel, neighbour_graph


def _probe(order=(1, 2, 3)):
    J = 240
    reset = np.zeros(J, bool)
    reset[0] = True
    Y = np.full((4, J), 100000.)
    Y[1, 40:80], Y[1, 80:120] = 180000, 120000
    Y[2, 40:80], Y[2, 80:120] = 120000, 180000
    Y[3, 40:120] = 60000
    Y = Y[[0, *order]]
    L = Y.sum(axis=1)
    XL = np.zeros((4, 2))
    m = SpecimenModel(Rows(Y, np.arange(4)), L, np.full((J, 1), 1 / J), reset, np.full(J, .1),
                      CloneParams(n_programmes=0), gamma=np.zeros(J))
    return m.shared_profiles(Y, L, XL)


def test_shared_doses_keep_opposite_events_a_third_clone_spans():
    """Tumour 1 and 2 differ by opposite events in bins 40-79 and 80-119; tumour 3's loss spans both. The shared
    dose contrast between 1 and 2 must keep the independently fitted contrast (it fell 0.099 -> 0.005)."""
    p = _probe()
    own = np.abs(p['dose_own'][1] - p['dose_own'][2]).max()
    shared = p['dose'][1] - p['dose'][2]
    assert np.abs(shared).max() >= 0.9 * own
    assert np.sign(shared[40:80].mean()) == -np.sign(shared[80:120].mean())


def test_shared_doses_do_not_depend_on_clone_order():
    a = _probe((1, 2, 3))['dose']
    b = _probe((3, 1, 2))['dose']          # group k of b holds original tumour (3, 1, 2)[k - 1]
    np.testing.assert_allclose(b[[2, 3, 1]], a[1:], atol=1e-6)


def test_neighbour_graph_has_no_self_edges_at_tied_coordinates():
    xy = np.repeat(np.arange(20, dtype=float)[:, None] * [10.0, 0.0], 3, axis=0)   # 20 points, 3 cells each
    for k in (1, 2, 4):
        G = neighbour_graph(xy, k)
        assert G.diagonal().sum() == 0
        assert np.allclose(np.asarray(G.sum(axis=1)).ravel(), 1.0)
    G1, G2 = neighbour_graph(xy, 2), neighbour_graph(xy, 2)
    assert (G1 != G2).nnz == 0
