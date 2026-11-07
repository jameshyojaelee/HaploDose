"""Optional model-conditional intervals for systematic relative RNA effects.

No pipeline calls this helper automatically. It does not estimate a covariance
from fitted HMM residuals, nor provide DNA or clone-identity confidence intervals.
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def systematic_rna_intervals(counts, expectation, covariance, reset, intervals, *,
                             counts_kind, count_provenance,
                             expectation_provenance, covariance_provenance,
                             bin_universe_provenance, confidence=.95):
    """Simultaneous systematic RNA intervals, including data-selected intervals.

    Contract: Y_j | Z ~ Poisson(E_j exp(mu_j + Z_j)), Z ~ N(0, Sigma).
    E is positive, independently frozen, and Sigma is independently known. The
    provenance strings identify those external artifacts; strings cannot verify
    their scientific independence. Fitted tau/v do not satisfy this contract by
    declaration alone. Sampling dispersion beyond conditional Poisson and error
    in E or Sigma are outside the guarantee.
    Counts must explicitly be declared raw, with their source and aggregation
    recorded. Integrality and provenance declarations cannot establish raw-count
    provenance or conditional Poisson sampling by themselves.

    The target for inclusive bin interval [a,b] is mean(mu[a:b+1]). Garwood
    intervals control all bin intensities with half the total error budget.
    Gaussian bounds control all nonempty intervals within reset-delimited blocks
    with the other half. Their union therefore covers arbitrary subsequently
    selected intervals under this model, without requiring independent bins.
    A reset must also mark any genomic gap across which intervals are ineligible.
    The bin universe and reset boundaries must be frozen independently of these
    observations; only the requested intervals may be selected from the data.

    Zero counts give an infinite lower endpoint. Bounds are not centered on an
    HMM posterior effect; the pseudocount log ratio is reported only descriptively.
    """
    if counts_kind != 'raw':
        raise ValueError('counts_kind must explicitly be raw; transformed or reconstructed counts are ineligible')
    for name, value in (('count_provenance', count_provenance),
                        ('expectation_provenance', expectation_provenance),
                        ('covariance_provenance', covariance_provenance),
                        ('bin_universe_provenance', bin_universe_provenance)):
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f'{name} must identify the independently frozen external input')
    if not np.isfinite(confidence) or not 0 < confidence < 1:
        raise ValueError('confidence must be finite and between zero and one')
    y, E = np.asarray(counts, np.float64), np.asarray(expectation, np.float64)
    if (y.ndim != 1 or y.size == 0 or E.shape != y.shape or not np.isfinite(y).all()
            or (y < 0).any() or (y != np.floor(y)).any() or (y > 2 ** 53 - 1).any()):
        raise ValueError('counts must be a nonempty vector of finite exactly representable nonnegative integers')
    if not np.isfinite(E).all() or (E <= 0).any():
        raise ValueError('expectation must be positive and finite with one entry per bin')
    J = len(y)
    starts = np.asarray(reset)
    if starts.shape != y.shape or starts.dtype.kind != 'b' or not starts[0]:
        raise ValueError('reset must be a boolean vector with its first entry True')
    block = np.cumsum(starts)
    boundaries = np.r_[np.flatnonzero(starts), J]
    lengths = np.diff(boundaries)
    family = sum(int(n) * (int(n) + 1) // 2 for n in lengths)
    requested = np.asarray(intervals)
    if requested.size == 0:
        requested = np.empty((0, 2), dtype=np.int64)
    if requested.ndim != 2 or requested.shape[1] != 2 or requested.dtype.kind not in 'iu':
        raise ValueError('intervals must contain inclusive integer [first_bin,last_bin] pairs')
    if ((requested < 0).any() or (requested >= J).any()
            or (requested[:, 0] > requested[:, 1]).any()):
        raise ValueError('interval bounds are outside the eligible bin vector')
    if np.any(block[requested[:, 0]] != block[requested[:, 1]]):
        raise ValueError('intervals cannot cross a chromosome or gap reset')
    S = np.asarray(covariance, np.float64)
    if S.shape != (J, J) or not np.isfinite(S).all():
        raise ValueError('covariance must be a finite square matrix with one row per bin')
    eps = np.finfo(np.float64).eps
    # This is an input diagnostic, not a certified eigensolver error bound. Do
    # not repair a negative eigenvalue using an arbitrary numerical tolerance.
    # Scientifically known PSD covariance remains part of the external contract.
    if not np.array_equal(S, S.T):
        raise ValueError('covariance must be exactly symmetric; no automatic repair is performed')
    minimum = float(np.linalg.eigvalsh(S)[0])
    if minimum < 0:
        raise ValueError('covariance must be positive semidefinite; a negative computed eigenvalue '
                         'requires independent resolution, not automatic repair')
    alpha = 1 - confidence
    alpha_count = alpha_latent = alpha / 2
    count_tail = alpha_count / (2 * J)
    latent_tail = alpha_latent / (2 * family)
    if count_tail == 0 or latent_tail == 0:
        raise ValueError('requested error allocation is not representable in float64')
    low_rate = np.zeros(J)
    positive = y > 0
    low_rate[positive] = stats.gamma.ppf(count_tail, y[positive])
    high_rate = stats.gamma.isf(count_tail, y + 1)
    if (not np.isfinite(low_rate).all() or (low_rate < 0).any()
            or not np.isfinite(high_rate).all() or (high_rate <= 0).any()):
        raise ArithmeticError('Poisson rate quantiles were not finite positive values')
    with np.errstate(divide='ignore'):
        low_log = np.log(low_rate) - np.log(E)
    high_log = np.log(high_rate) - np.log(E)
    critical = float(stats.norm.isf(latent_tail))
    rows = []
    for first, last in requested:
        a, b = int(first), int(last)
        n = b - a + 1
        region = S[a:b + 1, a:b + 1]
        # A summation-error upper bound prevents roundoff from reducing variance.
        terms = n * n
        sum_gamma = terms * eps / (1 - terms * eps)
        variance = max(0., float(region.sum()) + sum_gamma * float(np.abs(region).sum())) / n ** 2
        if not np.isfinite(variance):
            raise ArithmeticError('Covariance summation overflowed; interval uncertainty is unresolved')
        latent_halfwidth = critical * np.sqrt(variance)
        count_low, count_high = float(low_log[a:b + 1].mean()), float(high_log[a:b + 1].mean())
        low, high = count_low - latent_halfwidth, count_high + latent_halfwidth
        rows.append(dict(first_bin=a, last_bin=b, n_bins=n, lower=float(low), upper=float(high),
                         width=float(high - low), count_log_lower=count_low, count_log_upper=count_high,
                         count_log_width=count_high - count_low,
                         latent_mean_variance_upper=float(variance), latent_halfwidth=float(latent_halfwidth),
                         observed_mean_pseudocount_log_ratio=float(
                             (np.log(y[a:b + 1] + .5) - np.log(E[a:b + 1] + .5)).mean())))
    return dict(intervals=rows, confidence=float(confidence), alpha_count=alpha_count,
                alpha_latent=alpha_latent, count_family_size=J, latent_family_size=family,
                count_two_sided_tail=count_tail, latent_two_sided_tail=latent_tail,
                latent_critical_value=critical, covariance_diagonal_padding=0.,
                covariance_minimum_computed_eigenvalue=minimum,
                covariance_psd_check='float64_eigenvalue_diagnostic_without_repair',
                counts_kind=counts_kind, count_provenance=count_provenance,
                expectation_provenance=expectation_provenance, covariance_provenance=covariance_provenance,
                bin_universe_provenance=bin_universe_provenance,
                method='simultaneous_Garwood_and_Gaussian_union_bound',
                estimand='population_mean_systematic_rna_log_intensity_ratio_relative_to_frozen_expectation',
                scope='model_conditional_simultaneous_over_all_nonempty_chromosome_or_gap_contained_intervals',
                error_control='Bonferroni_within_each_family_then_union_bound_between_families',
                posterior_effect_interval=False, biological_validation=False)
