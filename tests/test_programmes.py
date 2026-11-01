import numpy as np

import pandas as pd

from haplodose.bins import arm_starts
from haplodose.kernels.estep import programme_scores
from haplodose.programmes import fit_amplitudes, high_pass


def test_high_pass_keeps_short_structure_only():
    reset = np.zeros(200, bool)
    reset[[0, 100]] = True
    arm = np.where(np.arange(200) < 100, 0.4, -0.2)            # chromosome-wide levels
    spike = np.zeros(200)
    spike[140:145] = 1.0                                        # 5-bin cluster
    W = high_pass(np.stack([arm, spike]), reset, window=25)
    assert np.abs(W[0]).max() < 1e-9
    assert W[1, 140:145].min() > 0.7
    # A chromosome-wide gain is not explained by high-passed programmes.
    shift = fit_amplitudes(arm, np.ones(200), W)[0]
    assert np.abs(shift).max() < 1e-6


def test_programme_scores_bounded_in_sparse_cells():
    rng = np.random.default_rng(0)
    J, F = 300, 4
    W = rng.normal(0, 3.0, (F, J))                              # large loadings
    Y = np.zeros((5, J), np.float32)
    Y[:, :3] = 50                                               # extreme cells
    u = programme_scores(Y, Y.sum(axis=1) + 1, np.ones((5, 1), np.float32), np.full((1, J), -np.log(J)),
                         np.full(J, 0.1), W, np.full(F, 1e-6))
    assert np.isfinite(u).all()
    assert np.abs(u @ W).max() <= 1.5 + 1e-4                   # three steps of at most 0.5


def test_arm_mean_removal_keeps_sub_arm_structure():
    bins = pd.DataFrame({'chromosome': ['chr4'] * 100, 'start': np.arange(100) * 1_000_000 + 1})
    reset = arm_starts(bins)                                    # hg38 chr4 p-arm ends at 50 Mb
    assert np.flatnonzero(reset).tolist() == [0, 50]
    arm = np.where(np.arange(100) < 50, 0.4, -0.2)              # p-arm gain, q-arm loss
    block = np.zeros(100)
    block[60:75] = 1.0                                          # 15-bin cluster
    W = high_pass(np.stack([arm, block]), reset, window=0)
    assert np.abs(W[0]).max() < 1e-9
    assert W[1, 60:75].min() > 0.65
    assert np.abs(fit_amplitudes(arm, np.ones(100), W)[0]).max() < 1e-6
