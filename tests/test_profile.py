"""Pooled segmentation and the clone-difference test on synthetic profiles."""
import numpy as np

from haplodose.profile import ProfileParams, difference_segments, fit_pooled, segment


def chromosomes(n_chrom=10, per=60):
    reset = np.zeros(n_chrom * per, bool)
    reset[::per] = True
    return reset


def test_pooled_gain_and_loss_recovered():
    rng = np.random.default_rng(0)
    reset = chromosomes()
    J = len(reset)
    E = np.full(J, 400.0)
    dose = np.zeros(J)
    dose[60:120] = 0.6 * np.log(1.5)
    dose[300:330] = 0.6 * np.log(0.5)
    Y = rng.poisson(E * np.exp(dose + rng.normal(0, 0.1, J))).astype(float)
    p = ProfileParams()
    centre, beta, tau = fit_pooled(Y, E, reset, p)
    post, path, mean, _ = segment(Y[None], E[None], reset, [centre], beta, [tau], p)
    assert (path[0, 62:118] == 2).mean() > 0.95
    assert (path[0, 302:328] == 0).mean() > 0.95
    assert (path[0, 130:290] == 1).mean() > 0.97
    assert 0.08 < tau < 0.14


def test_difference_ignores_shared_baseline_and_finds_subclonal_step():
    rng = np.random.default_rng(1)
    reset = chromosomes()
    J = len(reset)
    shared = rng.lognormal(0, 0.8, J)              # baseline mismatch shared by both groups
    step = np.zeros(J)
    step[425:475] = 0.6 * np.log(1.5)          # inside one chromosome (bins 420-479)
    Ya = rng.poisson(5e4 * shared / shared.sum() * 10).astype(float)
    rate_b = shared * np.exp(step)
    Yb = rng.poisson(5e4 * rate_b / rate_b.sum() * 10).astype(float)
    p = ProfileParams()
    mask = np.arange(J) % 2 == 1
    same = difference_segments(Ya, Ya.sum(), rng.poisson(Ya).astype(float), Ya.sum(), reset, 0.6, p, mask)[0]
    diff = difference_segments(Ya, Ya.sum(), Yb, Yb.sum(), reset, 0.6, p, mask)[0]
    assert same['segments'] == []
    assert len(diff['segments']) == 1
    s = diff['segments'][0]
    assert s["direction"] == "lower" and 421 <= s["first_bin"] <= 429 and 471 <= s["last_bin"] <= 479
