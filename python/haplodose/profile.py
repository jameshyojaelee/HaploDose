"""Segmentation of pooled (clone or group) profiles and tests for clone differences.

Pooled bin counts Y_j are compared with expected counts E_j from the reference
baseline. log E[Y_j] = log E_j + centre + beta * log(cn / 2) + o_j, with o_j the
bin-level baseline mismatch (sd tau). States are relative copy number 1-4 on a
diploid scale; beta is the expression-dosage slope (dosage compensation and
purity). Not absolute DNA copy number.

Two clones of one tumour share baseline and mismatch, so their log-rate
difference does not depend on the reference at all. `difference_segments`
tests that difference for contiguous copy-number steps.
"""
from __future__ import annotations

from dataclasses import dataclass

import jax.numpy as jnp
import numpy as np

from .kernels.hmm import forward_backward_batch, sticky_transition, viterbi_batch
from .kernels.likelihood import gaussian_logpdf, pln_logpmf

CN = np.array([1, 2, 3, 4])
LOG_DOSE = np.log(CN / 2.0).astype(np.float32)
STATE_LABELS = np.array(['loss', 'neutral', 'gain', 'amplification'])
NEUTRAL = 1
STATE_PRIOR = np.array([0.2, 0.5, 0.2, 0.1], np.float32)
# Difference states: lower by one copy step, none, higher. (Fractional steps, tried 2026-11-01,
# raised sensitivity to non-genetic differences more than to diluted real ones.)
DIFF_LEVELS = np.array([-1.0, 0.0, 1.0], np.float32)
DIFF_PRIOR = np.array([0.15, 0.7, 0.15], np.float32)


@dataclass
class ProfileParams:
    beta: float = 0.6
    beta_bounds: tuple = (0.3, 1.0)
    p_change: float = 2e-3
    min_bins: int = 3
    posterior: float = 0.95
    tau_floor: float = 0.02


def log_ratio(Y, E):
    return np.log((Y + 0.5) / (E + 0.5))


def weighted_median(x, w):
    order = np.argsort(x)
    cw = np.cumsum(w[order])
    return float(x[order][np.searchsorted(cw, 0.5 * cw[-1])])


def estimate_tau(r, Y, reset, floor=0.02, extra_var=None):
    """Extra-Poisson sd from first differences within chromosomes (robust to breakpoints).

    Var(r_j - r_{j-1}) = 2 tau^2 + sampling variances; the MAD ignores the few
    differences that span a breakpoint.
    """
    pv = 1.0 / (Y + 0.5) if extra_var is None else extra_var
    inside = ~reset[1:]
    d = np.diff(r)[inside]
    s = pv[1:][inside] + pv[:-1][inside]
    if len(d) < 10:
        return floor
    mad = 1.4826 * np.median(np.abs(d - np.median(d)))
    return float(np.sqrt(max(floor ** 2, (mad ** 2 - np.median(s)) / 2)))


def runs(path, reset):
    """Contiguous runs of equal state that do not cross chromosome starts."""
    out, start = [], 0
    for j in range(1, len(path) + 1):
        if j == len(path) or reset[j] or path[j] != path[start]:
            out.append((start, j - 1, int(path[start])))
            start = j
    return out


def segment(Y, E, reset, centre, beta, tau, params: ProfileParams):
    """Batch HMM over groups. Y, E [G, J]; centre [G]; tau [G] or per bin [G, J]. Returns posteriors, paths, mean log dose."""
    Y = np.atleast_2d(Y).astype(np.float32)
    E = np.atleast_2d(E).astype(np.float32)
    means = np.log(E)[..., None] + np.asarray(centre, np.float32)[:, None, None] + beta * LOG_DOSE
    tau = np.asarray(tau, np.float32)
    emit = pln_logpmf(Y[..., None], means, tau[:, None, None] if tau.ndim == 1 else tau[:, :, None])
    trans = sticky_transition(len(CN), params.p_change, STATE_PRIOR)
    init = jnp.log(jnp.asarray(STATE_PRIOR))
    post, ll = forward_backward_batch(emit, trans, init, jnp.asarray(reset))
    path = viterbi_batch(emit, trans, init, jnp.asarray(reset))
    post = np.asarray(post)
    return post, np.asarray(path), post @ (beta * LOG_DOSE), np.asarray(ll)


def fit_pooled(Y, E, reset, params: ProfileParams, rounds=4, v=None):
    """Centre, dosage slope beta and tau from one pooled profile (the whole specimen).

    v [J]: optional per-bin mismatch variance added to tau^2 (reference heterogeneity)."""
    v = np.zeros(len(Y)) if v is None else np.asarray(v)
    r = log_ratio(Y, E)
    w = 1.0 / (0.04 + 1.0 / (Y + 0.5))
    centre = weighted_median(r, w)
    tau = estimate_tau(r, Y, reset, params.tau_floor)
    beta = params.beta
    for _ in range(rounds):
        post, path, dose, _ = segment(Y[None], E[None], reset, [centre], beta, np.sqrt(tau ** 2 + v)[None], params)
        post, dose = post[0], dose[0]
        resid = r - dose
        w = 1.0 / (tau ** 2 + v + 1.0 / (Y + 0.5))
        centre = weighted_median(resid, w)
        aberrant = 1.0 - post[:, NEUTRAL]
        if aberrant.sum() >= 3 * params.min_bins:
            x = post @ LOG_DOSE                        # expected log(cn/2) per bin
            keep = aberrant > 0.5
            num = np.sum(w[keep] * x[keep] * (r[keep] - centre))
            den = np.sum(w[keep] * x[keep] ** 2)
            if den > 0:
                beta = float(np.clip(num / den, *params.beta_bounds))
    return centre, beta, tau


def group_tau(Y, E, centre, reset, floor):
    return np.array([estimate_tau(log_ratio(y, e) - c, y, reset, floor)
                     for y, e, c in zip(np.atleast_2d(Y), np.atleast_2d(E), np.atleast_1d(centre))])


def group_centres(Y, E, tau):
    out = []
    for y, e, t in zip(np.atleast_2d(Y), np.atleast_2d(E), np.atleast_1d(tau)):
        out.append(weighted_median(log_ratio(y, e), 1.0 / (t ** 2 + 1.0 / (y + 0.5))))
    return np.array(out)


def segments_table(post, path, dose, Y, E, reset, beta, tau, params: ProfileParams, z=3.0):
    """One row per HMM run, with posterior support and the smallest detectable log dosage.

    The run mean has se = sqrt(tau^2 / n_bins + 1 / U), U = pooled counts; a shift
    must exceed z * se to be seen, i.e. a log dosage of z * se / beta.
    """
    rows = []
    for a, b, s in runs(path, reset):
        sl = slice(a, b + 1)
        n = b - a + 1
        U = float(Y[sl].sum())
        se = float(np.sqrt(tau ** 2 / n + 1.0 / max(U, 1.0)))
        rows.append(dict(first_bin=a, last_bin=b, n_bins=n, cn_state=int(CN[s]),
                         state=STATE_LABELS[s], log_dose=float(dose[sl].mean()),
                         posterior=float(post[sl, s].mean()), pooled_UMI=U,
                         expected_UMI=float(E[sl].sum()), se_log_ratio=se,
                         min_detectable_log_dose=z * se / beta,
                         resolved=bool(s != NEUTRAL and n >= params.min_bins
                                       and post[sl, s].mean() >= params.posterior)))
    return rows


def difference_segments(Ya, La, Yb, Lb, reset, beta, params: ProfileParams, mask=None, shift=None):
    """Contiguous copy-number differences between pooled groups a and b.

    Inputs are [P, J] for P pairs. The test track uses bins where `mask` is True
    (cross-fitting: proposals are made on the other bins). `shift` [P, J] is the
    expected log-rate difference from known covariates (sequencing depth), removed
    before testing. Returns, per pair, resolved segments (bin indices in the full
    track) and a distance.
    """
    Ya, Yb = np.atleast_2d(Ya).astype(np.float64), np.atleast_2d(Yb).astype(np.float64)
    La, Lb = np.atleast_1d(La).astype(np.float64), np.atleast_1d(Lb).astype(np.float64)
    J = Ya.shape[1]
    idx = np.arange(J) if mask is None else np.flatnonzero(mask)
    sub_reset = _subset_reset(reset, idx)
    d = np.log((Ya[:, idx] + 0.5) / La[:, None]) - np.log((Yb[:, idx] + 0.5) / Lb[:, None])
    if shift is not None:
        d = d - np.atleast_2d(shift)[:, idx]
    v = 1.0 / (Ya[:, idx] + 0.5) + 1.0 / (Yb[:, idx] + 0.5)
    tau = np.empty(len(d))
    for p in range(len(d)):
        d[p] -= weighted_median(d[p], 1.0 / (0.04 + v[p]))
        tau[p] = estimate_tau(d[p], None, sub_reset, params.tau_floor, extra_var=v[p])
    var = tau[:, None] ** 2 + v
    step = beta * np.log(1.5)
    emit = gaussian_logpdf(d[..., None], jnp.asarray(DIFF_LEVELS * step), var[..., None])
    trans = sticky_transition(len(DIFF_LEVELS), params.p_change, DIFF_PRIOR)
    init = jnp.log(jnp.asarray(DIFF_PRIOR))
    post, _ = forward_backward_batch(jnp.asarray(emit, jnp.float32), trans, init, jnp.asarray(sub_reset))
    path = np.asarray(viterbi_batch(jnp.asarray(emit, jnp.float32), trans, init, jnp.asarray(sub_reset)))
    post = np.asarray(post)
    sign = np.sign(DIFF_LEVELS).astype(int)              # -1, 0, +1 per state
    results = []
    for p in range(len(d)):
        segs = []
        signed = sign[path[p]]
        for a, b, s in runs(signed + 1, sub_reset):
            n = b - a + 1
            s -= 1
            support = float(post[p, a:b + 1][:, sign == s].sum(axis=1).mean())
            if s == 0 or n < params.min_bins or support < params.posterior:
                continue
            w = 1.0 / var[p, a:b + 1]
            mean = float(np.sum(w * d[p, a:b + 1]) / w.sum())
            segs.append(dict(first_bin=int(idx[a]), last_bin=int(idx[b]), n_test_bins=n,
                             direction='higher' if s > 0 else 'lower', mean_diff=mean,
                             z=mean * np.sqrt(w.sum()), posterior=support))
        dist = float(np.mean(d[p] ** 2 / var[p]))
        results.append(dict(segments=segs, distance=dist, tau=float(tau[p])))
    return results


def _subset_reset(reset, idx):
    chrom_id = np.cumsum(reset)[idx]
    out = np.ones(len(idx), dtype=bool)
    out[1:] = chrom_id[1:] != chrom_id[:-1]
    return out
