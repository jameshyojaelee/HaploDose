"""Optional gene-level Laplace working smoother for descriptive RNA profiles.

The penalty and scale are working losses, not calibrated likelihood criteria.
Exact integer costs avoid cancellation and arbitrary decision tolerances.
"""
from fractions import Fraction
import heapq
import math

import numpy as np
import pandas as pd

from .smoothing import _running_median
from .stacking import stacked_rna_profile, _training_predictions, _segmented, _convex_l1


def _block_partition(groups, scale, *, capture=None):
    """Exact median interval costs and DP on one common integer objective grid."""
    ratios = [[float(v).as_integer_ratio() for v in group] for group in groups]
    denominator = max(d for group in ratios for _, d in group)
    integer = [[n*(denominator//d) for n, d in group] for group in ratios]
    n = len(groups)
    samples = sum(map(len, groups))
    penalty = 3*math.log(samples) if samples > 1 else 0.
    sn, sd = float(scale).as_integer_ratio()
    pn, pd = penalty.as_integer_ratio()
    cost_denominator = denominator*sn
    unit = math.lcm(cost_denominator, pd)
    multiplier = sd*(unit//cost_denominator)
    penalty_units = pn*(unit//pd)
    costs, medians, counts = ([[None]*n for _ in range(n)] for _ in range(3))
    for lo in range(n):
        left, right, lower_sum, upper_sum = [], [], 0, 0
        for hi in range(lo, n):
            for value in integer[hi]:
                if not left or value <= -left[0]:
                    heapq.heappush(left, -value)
                    lower_sum += value
                else:
                    heapq.heappush(right, value)
                    upper_sum += value
                if len(left) > len(right)+1:
                    moved = -heapq.heappop(left)
                    heapq.heappush(right, moved)
                    lower_sum -= moved
                    upper_sum += moved
                if len(right) > len(left):
                    moved = heapq.heappop(right)
                    heapq.heappush(left, -moved)
                    upper_sum -= moved
                    lower_sum += moved
            median_twice = -2*left[0] if len(left) > len(right) else -left[0]+right[0]
            twice_cost = median_twice*(len(left)-len(right))-2*lower_sum+2*upper_sum
            assert twice_cost >= 0 and (not right or -left[0] <= right[0])
            costs[lo][hi] = twice_cost*multiplier
            medians[lo][hi] = median_twice
            counts[lo][hi] = len(left)+len(right)
    dp, back, ties = [-penalty_units]+[None]*n, [0]*(n+1), []
    for hi in range(1, n+1):
        candidates = [dp[lo]+costs[lo][hi-1]+penalty_units for lo in range(hi)]
        minimum = min(candidates)
        winners = [lo for lo, value in enumerate(candidates) if value == minimum]
        dp[hi], back[hi] = minimum, winners[0]
        if len(winners) > 1:
            ties.append(dict(end_bin=hi-1, predecessors=winners))
    output, parts, hi = np.empty(n), [], n
    while hi:
        lo = back[hi]
        value = float(Fraction(medians[lo][hi-1], 2*denominator))
        output[lo:hi] = value
        parts.append(dict(first_bin=lo, last_bin=hi-1, effect=value, genes=counts[lo][hi-1]))
        hi = lo
    if capture is not None:
        capture(costs, medians, unit, denominator, penalty_units)
    return output, dict(segments=parts, penalty=penalty, observations=samples,
        objective=float(Fraction(dp[-1], unit)), exact_objective_numerator=str(dp[-1]),
        exact_objective_denominator=str(unit), exact_tied_predecessors=ties,
        arithmetic='Exact integer interval costs and DP for observed float effects, scale and penalty')


def _partition_profile(values, gene_bins, blocks, *, capture=None):
    values, gene_bins, blocks = np.asarray(values, float), np.asarray(gene_bins), np.asarray(blocks)
    groups = [values[gene_bins == j] for j in range(len(blocks))]
    if any(not len(g) for g in groups) or not np.isfinite(values).all():
        raise ValueError('Finite effects and at least one training gene per bin required')
    residual = np.concatenate([abs(g-np.median(g)) for g in groups])
    scale = float(np.median(residual)/math.log(2.))
    if not np.isfinite(scale) or scale <= 0:
        return None, dict(status='positive_working_scale_not_identified')
    output, records = np.empty(len(blocks)), []
    starts = np.flatnonzero(np.r_[True, blocks[1:] != blocks[:-1]])
    for lo, hi in zip(starts, np.r_[starts[1:], len(blocks)]):
        fit, record = _block_partition(groups[lo:hi], scale,
            capture=None if capture is None else lambda *args: capture(lo, hi, *args))
        output[lo:hi] = fit
        for segment in record['segments']:
            segment['first_bin'] += int(lo)
            segment['last_bin'] += int(lo)
        records.append(dict(first_bin=int(lo), last_bin=int(hi-1), **record))
    return output, dict(status='estimated', working_Laplace_scale=scale, blocks=records,
        scale='Training gene within-bin absolute residual median divided by log2',
        calibration_status='unvalidated', interpretation='Descriptive loss on prepared RNA effects; not a count model')


def _select(values, gene_bins, blocks, *, capture=None, solver='primal'):
    original, previous = stacked_rna_profile(values, gene_bins, blocks, solver=solver)
    values, gene_bins, blocks = np.asarray(values, float), np.asarray(gene_bins), np.asarray(blocks)
    base = previous['baseline']
    full, full_record = _partition_profile(values, gene_bins, blocks,
        capture=None if capture is None else capture['partition']('full'))
    if base['selected_width'] is None or full is None:
        return original, dict(selected=False, status='insufficient_new_expert_support', previous=previous,
                              full_expert=full_record), full
    widths, names = base['widths'], previous['experts']
    use_mdl = names[-1] == 'piecewise_mean_mdl'
    folds = pd.Series(gene_bins).groupby(gene_bins, sort=False).cumcount().to_numpy()%10
    sizes = np.bincount(gene_bins, minlength=len(blocks))
    eligible = sizes[gene_bins] >= 2
    P, fold_records = np.full((len(values), len(names)+1), np.nan), []
    for fold in range(10):
        medians, mdl, hold, _ = _training_predictions(values, gene_bins, blocks, widths, folds, fold)
        if not hold.any():
            continue
        candidate, record = _partition_profile(values[~hold], gene_bins[~hold], blocks,
            capture=None if capture is None else capture['partition'](fold))
        fold_records.append(dict(fold=fold, targets=int(hold.sum()), **record))
        if candidate is None:
            return original, dict(selected=False, status='masked_new_expert_unidentified', previous=previous,
                                  full_expert=full_record, folds=fold_records), full
        P[np.ix_(hold, np.arange(len(widths)))] = medians[:, gene_bins[hold]].T
        if use_mdl:
            assert mdl is not None
            P[hold, len(widths)] = mdl[gene_bins[hold]]
        P[hold, -1] = candidate[gene_bins[hold]]
    train = np.array([np.median(values[gene_bins == j]) for j in range(len(blocks))])
    F = np.array([_running_median(train, blocks, width) for width in widths])
    if use_mdl:
        mdl, _ = _segmented(train, blocks)
        assert mdl is not None
        F = np.vstack([F, mdl])
    assert P[eligible].shape[1] == len(F)+1 and np.isfinite(P[eligible]).all()
    if capture is not None:
        capture.update(P=P.copy(), eligible=eligible.copy(), bank=np.vstack([F, full]),
                       previous=previous, widths=widths)
    if previous['selected']:
        old_weights = np.asarray(previous['stacking']['weights'])
        np.testing.assert_array_equal(old_weights@F, original)
        old_risk = previous['stacking']['mean_absolute_error']
        old_uncertainty = previous['stacking']['optimization_gap_upper']+2*previous['stacking']['rounding_budget']
    else:
        index = widths.index(base['selected_width'])
        old_risk = float(math.fsum(map(float, abs(P[eligible, index]-values[eligible])))/eligible.sum())
        old_uncertainty = base['risk_sum_rounding_allowance']/eligible.sum()
    weights, fit = _convex_l1(P[eligible], values[eligible], solver=solver)
    chosen = bool(weights is not None and fit['mean_absolute_error']+fit['optimization_gap_upper']+
                  2*fit['rounding_budget'] < old_risk-old_uncertainty)
    output = weights@np.vstack([F, full]) if chosen else original
    return output, dict(selected=chosen, status=fit['status'], previous=previous,
        previous_prediction_risk=old_risk, previous_numerical_budget=old_uncertainty,
        stacking=fit, experts=names+['gene_Laplace_median'], full_expert=full_record,
        folds=fold_records, interpretation='Conditional RNA fitting, not independent confirmation'), full

