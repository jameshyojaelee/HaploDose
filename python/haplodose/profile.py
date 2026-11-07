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


def credible_difference_runs(path, post, reset, params):
    """Maximal contiguous difference intervals supported at every retained bin.

    A Viterbi run can include uncertain boundaries. Those boundaries must not
    erase a sufficiently long interior that already meets the same posterior
    criterion. A low-support interior bin breaks the interval; it is never
    bridged or retained by averaging stronger neighbours.
    """
    path, post = np.asarray(path), np.asarray(post)
    support = post[np.arange(len(path)), path]
    confident = np.where(support >= params.posterior, path, 1)
    return [(a, b, s) for a, b, s in runs(confident, reset)
            if s != 1 and b - a + 1 >= params.min_bins]


def segment(Y, E, reset, centre, beta, tau, params: ProfileParams, decode=True):
    """Batch HMM over groups. Y, E [G, J]; centre [G]; tau [G] or per bin [G, J]. Returns posteriors, Viterbi
    paths (None when decode is False: fits that use only posteriors skip the decoding pass), mean log dose and
    log-likelihood."""
    Y = np.atleast_2d(Y).astype(np.float32)
    E = np.atleast_2d(E).astype(np.float32)
    means = np.log(E)[..., None] + np.asarray(centre, np.float32)[:, None, None] + beta * LOG_DOSE
    tau = np.asarray(tau, np.float32)
    emit = pln_logpmf(Y[..., None], means, tau[:, None, None] if tau.ndim == 1 else tau[:, :, None])
    trans = sticky_transition(len(CN), params.p_change, STATE_PRIOR)
    init = jnp.log(jnp.asarray(STATE_PRIOR))
    post, ll = forward_backward_batch(emit, trans, init, jnp.asarray(reset))
    path = np.asarray(viterbi_batch(emit, trans, init, jnp.asarray(reset))) if decode else None
    post = np.asarray(post)
    return post, path, post @ (beta * LOG_DOSE), np.asarray(ll)


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
        post, _, dose, _ = segment(Y[None], E[None], reset, [centre], beta, np.sqrt(tau ** 2 + v)[None], params,
                                   decode=False)
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


def segments_table(post, path, dose, Y, E, reset, beta, tau, params: ProfileParams, z=3.0,
                   reference_var=None):
    """One row per HMM run, with explicit RNA-effect units and conditional uncertainty.

    The independent-bin approximation is the variance of the unweighted mean
    log ratio: sum_j(tau_j^2 + reference_var_j + 1/(Y_j+0.5)) / n_bins^2.
    It omits parameter estimation, genomic dependence and interval selection;
    z times this SE is a scale diagnostic, not a power/detection guarantee.
    log_dose and se_log_ratio are deprecated aliases in RNA log-rate units;
    The historical detection-limit field is deliberately absent from schema 2.
    """
    Y, E = np.asarray(Y, np.float64), np.asarray(E, np.float64)
    tau = np.broadcast_to(np.asarray(tau, np.float64), Y.shape)
    reference_var = np.zeros_like(Y) if reference_var is None else np.broadcast_to(
        np.asarray(reference_var, np.float64), Y.shape)
    if (not np.isfinite(tau).all() or (tau < 0).any() or not np.isfinite(reference_var).all()
            or (reference_var < 0).any() or not np.isfinite(beta) or beta <= 0):
        raise ValueError('Segment uncertainty requires finite nonnegative variances and positive beta')
    rows = []
    for a, b, s in runs(path, reset):
        sl = slice(a, b + 1)
        n = b - a + 1
        U = float(Y[sl].sum())
        se = float(np.sqrt(np.sum(tau[sl] ** 2 + reference_var[sl] + 1.0 / (Y[sl] + 0.5))) / n)
        rna_effect = float(dose[sl].mean())
        rows.append(dict(first_bin=a, last_bin=b, n_bins=n, cn_state=int(CN[s]),
                         state=STATE_LABELS[s], log_dose=rna_effect,
                         rna_log_effect=rna_effect, relative_state_log_ratio=rna_effect / beta,
                         observed_mean_rna_log_ratio=float(log_ratio(Y[sl], E[sl]).mean()),
                         rna_effect_estimator='HMM_posterior_mean_state_effect',
                         posterior=float(post[sl, s].mean()), pooled_UMI=U,
                         pooled_count=U, expected_count=float(E[sl].sum()),
                         count_units='input_count_units',
                         expected_UMI=float(E[sl].sum()), se_log_ratio=se,
                         se_rna_log_ratio_approx=se, three_sigma_rna_log_ratio_approx=3.0 * se,
                         uncertainty_method='independent_bin_log_ratio_approximation',
                         uncertainty_statistic='unweighted_mean_pseudocount_log_count_ratio',
                         uncertainty_center_field='observed_mean_rna_log_ratio',
                         uncertainty_estimand=None,
                         uncertainty_status='approximate_scale_diagnostic_without_validated_parameter_interval',
                         posterior_effect_interval_validated=False,
                         uncertainty_scope='conditional_on_fit_and_selected_interval',
                         uncertainty_calibrated=False,
                         uncertainty_scale_z=float(z), scaled_se_rna_log_ratio_approx=z * se,
                         resolved=bool(s != NEUTRAL and n >= params.min_bins
                                       and post[sl, s].mean() >= params.posterior)))
    return rows


def pooled_lograte(Y, L, count_variance=None, weight_scale=None):
    """Conditional Poisson sampling term for weighted counts, with weight-unit smoothing.

    S=sum(wY), Q=sum(w²Y), c=sum(w²)/sum(w). The pseudocount .5*c
    and variance regularizer .5*c² preserve information under weight rescaling.
    This is not a correction for NB noise, fitted weights or interval selection.
    Zero-exposure rows get finite placeholders; callers must mark them untestable.
    """
    Y = np.atleast_2d(np.asarray(Y, np.float64))
    L = np.broadcast_to(np.atleast_1d(np.asarray(L, np.float64)), (len(Y),))
    Q = Y if count_variance is None else np.broadcast_to(np.asarray(count_variance, np.float64), Y.shape)
    c = np.ones(len(Y)) if weight_scale is None else np.broadcast_to(
        np.atleast_1d(np.asarray(weight_scale, np.float64)), (len(Y),))
    if any(not np.isfinite(x).all() or np.any(x < 0) for x in (Y, L, Q, c)):
        raise ValueError('Pooled counts, exposure, count variance and weight scale must be finite and nonnegative')
    if np.any((L > 0) & (c <= 0)):
        raise ValueError('Positive pooled exposure requires positive weight scale')
    safe_c = np.where(c > 0, c, 1.)[:, None]
    adjusted = Y + .5 * safe_c
    lr = np.log(adjusted / np.where(L > 0, L, 1.)[:, None])
    v = ((Q + .5 * safe_c ** 2) / adjusted ** 2
         if not (np.array_equal(Q, Y) and np.all(safe_c == 1)) else 1. / adjusted)
    lr[L <= 0] = 0.
    return lr, v


def difference_segments(Ya, La, Yb, Lb, reset, beta, params: ProfileParams, mask=None, shift=None,
                        *, count_variance=None, weight_scale=None):
    """Contiguous copy-number differences between pooled groups a and b.

    Inputs are [P, J] for P pairs. The test track uses bins where `mask` is True
    (cross-fitting: proposals are made on the other bins). `shift` [P, J] is the
    expected log-rate difference from known covariates (sequencing depth), removed
    before testing. Returns, per pair, resolved segments (bin indices in the full
    track), a distance and `testable`.

    Retained intervals are contiguous credible cores of the Viterbi differences:
    every retained bin must meet the existing posterior criterion, and each core
    must still meet min_bins. Uncertain flanks and interior gaps are excluded.

    A pair is untestable when either group has no UMIs (e.g. depth-matched pooling of
    groups with no depth overlap gives both zero weight) or its log ratios are not
    finite; it gets no segments, distance inf and testable False, so callers treat it
    as showing no difference (merged last, split rejected) instead of reading NaN.
    """
    Ya, Yb = np.atleast_2d(Ya).astype(np.float64), np.atleast_2d(Yb).astype(np.float64)
    La, Lb = np.atleast_1d(La).astype(np.float64), np.atleast_1d(Lb).astype(np.float64)
    J = Ya.shape[1]
    idx = np.arange(J) if mask is None else np.flatnonzero(mask)
    sub_reset = _subset_reset(reset, idx)
    Qa, Qb = (None, None) if count_variance is None else count_variance
    ca, cb = (None, None) if weight_scale is None else weight_scale
    ra, va = pooled_lograte(Ya, La, Qa, ca)
    rb, vb = pooled_lograte(Yb, Lb, Qb, cb)
    d = (ra - rb)[:, idx]
    if shift is not None:
        d = d - np.atleast_2d(shift)[:, idx]
    testable = (La > 0) & (Lb > 0) & np.isfinite(d).all(axis=1)
    d[~testable] = 0.0                                   # placeholder so the batched HMM stays finite
    v = (va + vb)[:, idx]
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
        if not testable[p]:
            results.append(dict(segments=[], distance=float('inf'), tau=float('nan'), testable=False))
            continue
        segs = []
        for a, b, s in credible_difference_runs(path[p], post[p], sub_reset, params):
            n = b - a + 1
            s -= 1
            support = float(post[p, a:b + 1][:, sign == s].sum(axis=1).mean())
            if s == 0 or n < params.min_bins or support < params.posterior:
                continue
            w = 1.0 / var[p, a:b + 1]
            mean = float(np.sum(w * d[p, a:b + 1]) / w.sum())
            segs.append(dict(first_bin=int(idx[a]), last_bin=int(idx[b]), n_test_bins=n,
                             direction='higher' if s > 0 else 'lower', mean_diff=mean,
                             z=mean * np.sqrt(w.sum()), posterior=support,
                             min_bin_posterior=float(post[p, a:b + 1, s + 1].min()),
                             posterior_rule='every_retained_bin'))
        dist = float(np.mean(d[p] ** 2 / var[p]))
        results.append(dict(segments=segs, distance=dist, tau=float(tau[p]), testable=True))
    return results


def _subset_reset(reset, idx):
    chrom_id = np.cumsum(reset)[idx]
    out = np.ones(len(idx), dtype=bool)
    out[1:] = chrom_id[1:] != chrom_id[:-1]
    return out
