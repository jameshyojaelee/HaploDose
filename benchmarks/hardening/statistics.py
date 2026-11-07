"""Versioned, biological-unit inference for development benchmarks.

The bootstrap here conditions on supplied profiles and thresholds. When those are
estimated, callers must refit them within an outer resampling loop; this module
does not turn a conditional interval into an unconditional one.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import beta, binom

SCORE_VERSION = '2.0'


def holm(pvalues):
    """Holm adjustment over the entire supplied family, including untestable p=1."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1:
        raise ValueError('p-values must be one dimensional')
    p = np.where(np.isfinite(p), p, 1.0)
    if np.any((p < 0) | (p > 1)):
        raise ValueError('p-values must lie in [0, 1]')
    order = np.argsort(p, kind='stable')
    out = np.empty(len(p))
    out[order] = np.minimum(1.0, np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return out


def paired_mean_comparison(a, b, *, seed=20261106, n_boot=10000, alpha=0.05,
                           family_size=1, clusters=None):
    """Mean paired effect, basic-bootstrap intervals and centered-bootstrap p-value.

    Missing pairs are reported, never independently realigned. Optional clusters
    resample whole biological/reference groups, retaining all their observations.
    This is asymptotic bootstrap inference, not an exact randomization test.
    A constant observed difference cannot estimate sampling variance and is
    conservatively untestable. Family intervals use Bonferroni; p-values use Holm
    in ``adjust_family``. Both refer to the mean, unlike a signed-rank/mean hybrid.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    if a.ndim != 1 or a.shape != b.shape:
        raise ValueError('paired inputs must be equal-length vectors')
    if not 0 < alpha < 1 or family_size < 1 or n_boot < 99:
        raise ValueError('require alpha in (0,1), family_size >= 1 and n_boot >= 99')
    valid = np.isfinite(a) & np.isfinite(b)
    d = (a - b)[valid]
    groups = None
    if clusters is not None:
        clusters = np.asarray(clusters)
        if clusters.shape != a.shape:
            raise ValueError('cluster labels must match paired inputs')
        labels = clusters[valid]
        if any(x is None or str(x) in ('', 'nan', 'None') for x in labels):
            raise ValueError('cluster labels must be present for every usable pair')
        groups = [np.flatnonzero(labels == x) for x in np.unique(labels)]
    units = len(d) if groups is None else len(groups)
    mean = float(d.mean()) if len(d) else np.nan
    testable = units >= 2 and len(d) >= 2 and np.ptp(d) > 0
    result = dict(score_version=SCORE_VERSION, n=len(d), n_missing=int((~valid).sum()),
                  missing_a=int((~np.isfinite(a)).sum()), missing_b=int((~np.isfinite(b)).sum()),
                  n_resampling_units=units, mean_a=float(a[valid].mean()) if len(d) else np.nan,
                  mean_b=float(b[valid].mean()) if len(d) else np.nan, mean_diff=mean,
                  median_a=float(np.median(a[valid])) if len(d) else np.nan,
                  median_b=float(np.median(b[valid])) if len(d) else np.nan,
                  median_diff=float(np.median(d)) if len(d) else np.nan,
                  a_higher=int((d > 0).sum()), b_higher=int((d < 0).sum()), ties=int((d == 0).sum()),
                  testable=bool(testable), p=1.0, ci_low=np.nan, ci_high=np.nan,
                  family_ci_low=np.nan, family_ci_high=np.nan, family_size=int(family_size),
                  alpha=alpha, bootstrap_replicates=n_boot, seed=seed,
                  inference='paired mean, centered bootstrap; basic intervals; conditional on supplied inputs',
                  resampling='independent supplied units' if groups is None else 'whole supplied clusters')
    if not testable:
        result['reason'] = 'fewer than two units or no estimable paired-difference variance'
        return result
    rng = np.random.default_rng(seed)
    boot = np.empty(n_boot)
    # Bound allocations even for large cohorts; no B x N matrix is retained.
    if groups is None:
        batch = max(1, min(256, 1_000_000 // len(d)))
        for lo in range(0, n_boot, batch):
            hi = min(lo + batch, n_boot)
            boot[lo:hi] = d[rng.integers(len(d), size=(hi - lo, len(d)))].mean(axis=1)
    else:
        sums = np.asarray([d[g].sum() for g in groups])
        sizes = np.asarray([len(g) for g in groups])
        for lo in range(0, n_boot, 256):
            idx = rng.integers(len(groups), size=(min(256, n_boot - lo), len(groups)))
            boot[lo:lo + len(idx)] = sums[idx].sum(axis=1) / sizes[idx].sum(axis=1)
    error = boot - mean
    if np.ptp(boot) == 0:
        result['testable'] = False
        result['reason'] = 'resampling units have no estimable mean-effect variance'
        return result
    result['p'] = float((1 + np.count_nonzero(np.abs(error) >= abs(mean))) / (n_boot + 1))
    for prefix, tail in (('', alpha / 2), ('family_', alpha / (2 * family_size))):
        low, high = mean - np.quantile(error, [1 - tail, tail])
        result[prefix + 'ci_low'], result[prefix + 'ci_high'] = float(low), float(high)
    return result


def adjust_family(records, *, alpha=0.05):
    """Return records with a fixed Holm family and a simultaneous interval gate."""
    adjusted = holm([r.get('p', 1.0) if r.get('testable', False) else 1.0 for r in records])
    for r, p in zip(records, adjusted):
        if r.get('family_size') != len(records):
            raise ValueError('compute intervals with the full fixed family_size')
        r['p_holm'] = float(p)
        r['superior'] = bool(r.get('testable') and p < alpha and r['family_ci_low'] > 0)
    return records


def specimen_false_split(errors, total, *, target=0.05, alpha=0.05):
    """Exact one-sided binomial upper bound; each trial must be a biological specimen."""
    if int(errors) != errors or int(total) != total or not 0 <= errors <= total:
        raise ValueError('require integer 0 <= errors <= total')
    if not 0 < target < 1 or not 0 < alpha < 1:
        raise ValueError('target and alpha must lie in (0,1)')
    errors, total = int(errors), int(total)
    upper = 1.0 if errors == total else float(beta.ppf(1 - alpha, errors + 1, total - errors))
    return dict(errors=errors, n_biological_specimens=total,
                rate=errors / total if total else np.nan, upper=upper, confidence=1 - alpha,
                target=target, p_control=float(binom.cdf(errors, total, target)) if total else 1.0,
                controlled=bool(total and upper <= target),
                interval='one-sided Clopper-Pearson; assumes independent biological specimens')
