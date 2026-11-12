"""Descriptive RNA profile stacking; no DNA-state or calibration inference."""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import linprog

from .smoothing import _running_median, predictive_rna_profile


def _l1_primal(P, target):
    n, k = P.shape
    identity = sparse.eye(n, format='csr')
    A = sparse.vstack([sparse.hstack([sparse.csr_matrix(P), -identity]),
                      sparse.hstack([-sparse.csr_matrix(P), -identity])], format='csr')
    eq = sparse.csr_matrix((np.ones(k), (np.zeros(k, int), np.arange(k))), shape=(1, k+n))
    result = linprog(np.r_[np.zeros(k), np.full(n, 1/n)], A_ub=A,
                     b_ub=np.r_[target, -target], A_eq=eq, b_eq=[1.],
                     bounds=(0., None), method='highs', options={'threads': 1})
    if not result.success:
        return result, None, None
    return result, result.x[:k], n*(result.ineqlin.marginals[n:]-result.ineqlin.marginals[:n])


def _l1_dual(P, target):
    n, k = P.shape
    # The dual objective maximizes t - target @ s, where t is at most
    # every expert's inner product with s. Inequality prices recover weights.
    constraints = np.column_stack((-P.T, np.ones(k)))
    result = linprog(np.r_[target, -1.], A_ub=constraints, b_ub=np.zeros(k),
                     bounds=[(-1., 1.)]*n + [(None, None)],
                     method='highs', options={'threads': 1})
    if not result.success:
        return result, None, None
    return result, -result.ineqlin.marginals, result.x[:-1]


def _convex_l1(predictions, target, *, solver='primal'):
    """Simplex L1 fit, checked against an independently feasible dual bound."""
    P, target = np.asarray(predictions, float), np.asarray(target, float)
    if (P.ndim != 2 or target.shape != (len(P),) or not P.size
            or not np.isfinite(P).all() or not np.isfinite(target).all()):
        raise ValueError('Nonempty finite aligned expert predictions and targets required')
    n, k = P.shape
    diagnostics = {}
    if solver == 'dual':
        result, weights, dual = _l1_dual(P, target)
        diagnostics['solver'] = 'dual'
        if (not result.success or not np.isfinite(weights).all() or not np.isfinite(dual).all()
                or not np.any(weights > 0)):
            diagnostics.update(solver='primal_fallback', dual_solver_status=int(result.status))
            result, weights, dual = _l1_primal(P, target)
    elif solver == 'primal':
        result, weights, dual = _l1_primal(P, target)
    else:
        raise ValueError("solver must be 'primal' or 'dual'")
    if not result.success:
        return None, dict(status='unresolved_solver', solver_status=int(result.status),
                          message=result.message, **diagnostics)
    positive = np.maximum(weights, 0.)
    total = math.fsum(map(float, positive))
    if not total > 0:
        raise ArithmeticError('No positive simplex mass')
    weights = positive/total
    prediction = P@weights
    primal = math.fsum(map(float, abs(prediction-target)))/n
    # min_w max_|s|<=1 s'(Pw-y)/n. Clipping supplies a feasible s even
    # when the solver's residuals are imperfect. This is not an interval.
    dual = np.clip(dual, -1., 1.)
    products = [math.fsum(float(a*b) for a, b in zip(P[:, j], dual)) for j in range(k)]
    lower = (min(products)-math.fsum(float(a*b) for a, b in zip(target, dual)))/n
    operations = 32*n*(k+1)+256
    gamma = operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    magnitude = max(1., float(np.max(abs(P))), float(np.max(abs(target))))
    rounding = gamma*magnitude
    if primal+rounding < lower-rounding:
        raise ArithmeticError('Primal-dual ordering contradicts derived rounding bounds')
    return weights, dict(status='estimated', weights=weights.tolist(), mean_absolute_error=primal,
        feasible_dual_lower=lower, optimization_gap_upper=max(0., primal-lower)+2*rounding,
        rounding_budget=rounding, simplex_sum=float(weights.sum()), solver_message=result.message,
        certificate='Feasible clipped L1 dual; not an inferential confidence interval', **diagnostics)


def _partition(track, variance, reset):
    """Exact dynamic programming for the fixed Gaussian working segment score."""
    output, parts = np.empty(len(track)), []
    starts = np.flatnonzero(reset)
    for begin, end in zip(starts, np.r_[starts[1:], len(track)]):
        y, w = track[begin:end], 1/variance[begin:end]
        n = len(y)
        penalty = 3*np.log(n) if n > 1 else 0.
        A, B, C = [np.r_[0., np.cumsum(v)] for v in (w, w*y, w*y*y)]
        dp, back = np.full(n+1, np.inf), np.zeros(n+1, int)
        dp[0] = -penalty
        for hi in range(1, n+1):
            cost = (C[hi]-C[:hi])-(B[hi]-B[:hi])**2/(A[hi]-A[:hi])
            values = dp[:hi]+np.maximum(cost, 0)+penalty
            lo = int(np.argmin(values))
            dp[hi], back[hi] = values[lo], lo
        hi = n
        while hi:
            lo = int(back[hi])
            value = float(np.sum(w[lo:hi]*y[lo:hi])/w[lo:hi].sum())
            output[begin+lo:begin+hi] = value
            parts.append(dict(first_bin=int(begin+lo), last_bin=int(begin+hi-1), effect=value))
            hi = lo
    return output, parts


def _segmented(track, blocks):
    reset = np.r_[True, blocks[1:] != blocks[:-1]]
    differences = np.diff(track)[~reset[1:]]
    if not len(differences):
        return None, dict(status='no_within_block_differences')
    scale = float(np.median(np.abs(differences-np.median(differences)))/.6744897501960817/np.sqrt(2))
    if not np.isfinite(scale) or scale <= 0:
        return None, dict(status='positive_working_scale_not_identified')
    fitted, parts = _partition(track, np.full(len(track), scale**2), reset)
    return fitted, dict(status='estimated', working_sd=scale, segments=parts,
        penalty='3 log(block bins) per added segment', calibration_status='unvalidated')


def _training_predictions(values, gene_bins, blocks, widths, folds, excluded):
    sizes = np.bincount(gene_bins, minlength=len(blocks))
    hold = (folds == excluded) & (sizes[gene_bins] >= 2)
    track = np.array([np.median(values[(~hold) & (gene_bins == b)]) for b in range(len(blocks))])
    medians = np.array([_running_median(track, blocks, w) for w in widths])
    segmented, decision = _segmented(track, blocks)
    return medians, segmented, hold, decision


def stacked_rna_profile(gene_log_ratios, gene_bins, blocks, *, gene_segments=False,
                        average_segments=False, solver='primal'):
    """Combine smoothing scales by conditional masked-gene RNA prediction.

    Inputs have the same contract as ``predictive_rna_profile``: finite prepared
    natural-log RNA contrasts, complete zero-based integer bin IDs, and one
    contiguous arm/gap label per bin. Preserve gene identities and order: within
    each bin, rank modulo ten defines prediction folds. Singleton genes train
    neighbouring predictions but are never prediction targets.

    Each target is excluded from its base predictors. Simplex weights minimize
    absolute error on these predictions. The weights themselves use the target
    values; this is not independent confirmation or an unbiased held-out score.
    Upstream normalization, feature selection and reference fitting are fixed.
    The Gaussian segment expert is only a working smoother, not a count model.

    Return a descriptive bin profile and diagnostics. Retain the existing median
    selector unless the stack beats every expert beyond its feasible-dual gap
    and operation/magnitude rounding budget. No eligible targets means explicit
    abstention. This opt-in function supplies no CNA states, clone memberships,
    calibrated intervals or detection limits. Its benefit is dataset-specific;
    the fixed stack failed the bulk regression replacement gate.

    Set gene_segments=True to include an additional gene-level median segment
    expert with a Laplace working loss. It preserves the original expert bank
    and falls back unless RNA fit improves beyond numerical bounds. This mode
    improved the single-cell regression but failed bulk transfer. Its scale
    and 3 log(training genes) segment penalty are not calibrated error control.

    With gene_segments=True, average_segments=True also adds the average over
    all contiguous partitions under that same fitted loss. These are descriptive
    Gibbs weights, not integrated Bayesian evidence or calibrated boundaries.
    The added predictor uses the same masked genes, loss and numerical gates.
    Its benefit is currently supported only by exposed single-cell regression.

    Set solver='dual' to solve the same stacking loss with one constraint per
    expert. Both solvers check the objective with a feasible dual bound. If
    multiple weight vectors optimize the masked loss, their full profiles can
    differ. The default primal formulation preserves the historical choice.
    """
    if not isinstance(gene_segments, (bool, np.bool_)):
        raise ValueError('gene_segments must be boolean')
    if not isinstance(average_segments, (bool, np.bool_)):
        raise ValueError('average_segments must be boolean')
    if solver not in ('primal', 'dual'):
        raise ValueError("solver must be 'primal' or 'dual'")
    if average_segments:
        if not gene_segments:
            raise ValueError('average_segments requires gene_segments=True')
        from .average_segments import _fit
        profile, _, decision, _ = _fit(gene_log_ratios, gene_bins, blocks, solver=solver)
        return profile, dict(decision, method='conditional_partition_average_L1_stack',
            units='natural_log_relative_rna_rate', calibration_status='unvalidated')
    if gene_segments:
        from .gene_segments import _select
        profile, decision, _ = _select(gene_log_ratios, gene_bins, blocks, solver=solver)
        return profile, dict(decision, method='conditional_gene_Laplace_L1_stack',
            units='natural_log_relative_rna_rate', calibration_status='unvalidated')
    # Validate before integer conversion so fractional IDs cannot be truncated.
    baseline, record = predictive_rna_profile(gene_log_ratios, gene_bins, blocks)
    values = np.asarray(gene_log_ratios, float)
    gene_bins, blocks = np.asarray(gene_bins, np.int64), np.asarray(blocks)
    metadata = dict(method='conditional_masked_gene_L1_stack',
        units='natural_log_relative_rna_rate', calibration_status='unvalidated',
        interpretation='Conditional descriptive RNA stack; no independent confirmation or count calibration')
    if record['selected_width'] is None:
        return baseline, dict(metadata, status='insufficient_masked_gene_support',
                              selected=False, baseline=record)
    widths = record['widths']
    folds = pd.Series(gene_bins).groupby(gene_bins, sort=False).cumcount().to_numpy()%10
    eligible = np.bincount(gene_bins, minlength=len(blocks))[gene_bins] >= 2
    P = np.full((len(values), len(widths)+1), np.nan)
    partition_ok, fold_records = True, []
    for f in range(10):
        medians, segmented, hold, decision = _training_predictions(
            values, gene_bins, blocks, widths, folds, f)
        if not hold.any():
            continue
        P[np.ix_(hold, np.arange(len(widths)))] = medians[:, gene_bins[hold]].T
        if segmented is None:
            partition_ok = False
        else:
            P[hold, -1] = segmented[gene_bins[hold]]
        fold_records.append(dict(fold=f, targets=int(hold.sum()), partition=decision))
    full = np.array([np.median(values[gene_bins == b]) for b in range(len(blocks))])
    F = np.array([_running_median(full, blocks, w) for w in widths])
    segmented, full_partition = _segmented(full, blocks)
    partition_ok = partition_ok and segmented is not None
    names = [f'median_width_{w}' for w in widths]
    if partition_ok:
        F = np.vstack([F, segmented])
        names.append('piecewise_mean_mdl')
    else:
        P = P[:, :-1]
    if not np.isfinite(P[eligible]).all():
        raise ArithmeticError('Incomplete out-of-training prediction matrix')
    weights, fit = _convex_l1(P[eligible], values[eligible], solver=solver)
    single_risk = float(np.min(np.mean(abs(P[eligible]-values[eligible, None]), axis=0)))
    selected = bool(weights is not None and fit['mean_absolute_error']+fit['optimization_gap_upper']+
                    2*fit['rounding_budget'] < single_risk)
    return (weights@F if selected else baseline), dict(metadata, status=fit['status'], selected=selected,
        baseline=record, experts=names, best_single_expert_risk=single_risk, stacking=fit,
        eligible_genes=int(eligible.sum()), full_partition=full_partition, fold_partitions=fold_records)


def stacked_count_profile(counts, cell_ids, gene_ids, features, blocks,
                          query_cells, reference_cells, *, counts_kind, gene_segments=False,
                          average_segments=False, solver='primal'):
    """Descriptive specimen profile directly from an identified count matrix.

    ``counts`` is cells by genes, dense or sparse. ``features`` has unique
    ``gene`` and integer ``bin`` columns; its row order defines the canonical
    gene order used for prediction folds, independently of matrix column order.
    Query and reference IDs must be nonempty, disjoint and present in cell_ids.
    The reference is supplied by the caller, never chosen using DNA agreement.

    Declare counts_kind='raw' or 'reconstructed'. Selected observations require
    finite nonnegative integer values. Unselected observations are unused.
    This function is descriptive: reconstructed measurements do
    not acquire raw-count calibration, and the inferential run/bulk raw-only
    contract is unchanged. The fixed half-count convention matches the profile
    regression; it is not a confidence threshold. Library scaling uses only the
    selected features. Reference composition and transcription remain potential
    explanations for regional RNA effects.

    Returns (bin_profile, ordered_gene_log_ratios, diagnostics). No cell-level
    clone assignments, DNA states, calibrated intervals or detection limits are
    inferred. Single-cell profile regression gains do not establish a bulk or
    spatial benefit. No defaults or existing smoothing behavior are changed.
    gene_segments=True forwards the optional gene-level working smoother; the
    same RNA-only interpretation and failed bulk-transfer limitation apply.
    The solver argument is passed to ``stacked_rna_profile``.
    """
    if counts_kind not in ('raw', 'reconstructed'):
        raise ValueError('Declare counts_kind as raw or reconstructed; normalized inputs are unsupported')

    def ids(values, name):
        result = pd.Index(values)
        if (not len(result) or result.hasnans or not result.is_unique
                or any(not str(v).strip() for v in result)):
            raise ValueError(f'{name} must be nonempty, unique and nonmissing')
        return result

    cells, genes = ids(cell_ids, 'Matrix cells'), ids(gene_ids, 'Matrix genes')
    query, reference = ids(query_cells, 'Query cells'), ids(reference_cells, 'Reference cells')
    if len(query.intersection(reference)):
        raise ValueError('Query and reference populations must be disjoint')
    features = pd.DataFrame(features)
    if not features.columns.is_unique or not {'gene', 'bin'}.issubset(features.columns):
        raise ValueError('Ordered features require unique gene/bin columns')
    selected_genes = ids(features.gene, 'Feature genes')
    qr, rr, gc = cells.get_indexer(query), cells.get_indexer(reference), genes.get_indexer(selected_genes)
    if any(np.any(v < 0) for v in (qr, rr, gc)):
        raise ValueError('Every query/reference cell and feature gene must occur in the matrix axes')
    if np.ndim(counts) != 2:
        raise ValueError('A cells-by-genes matrix is required')
    matrix = sparse.csr_matrix(counts)
    if matrix.shape != (len(cells), len(genes)):
        raise ValueError('Count shape disagrees with cell/gene identities')
    if matrix.dtype.kind not in 'iuf':
        raise ValueError('A numeric count matrix is required')
    # Restrict first: converting a large, otherwise unused specimen collection
    # to float64 would allocate memory proportional to the full input matrix.
    Q, R = matrix[qr][:, gc], matrix[rr][:, gc]
    for population in (Q, R):
        if (not np.isfinite(population.data).all() or np.any(population.data < 0)
                or np.any(population.data != np.floor(population.data))
                or np.any(population.data >= 2**53)):
            raise ValueError('Finite nonnegative integer counts below the exact float64 integer limit required')
    # Positive integer subtotals below2**53 are exactly representable. Check
    # correctly rounded population totals before any sparse floating summation.
    Q, R = Q.astype(np.float64, copy=False), R.astype(np.float64, copy=False)
    nq, nr = math.fsum(map(float, Q.data)), math.fsum(map(float, R.data))
    if not (0 < nq < 2**53 and 0 < nr < 2**53):
        raise ValueError('Both selected-feature totals must be positive and exactly representable below2**53')
    y, r = np.asarray(Q.sum(axis=0)).ravel(), np.asarray(R.sum(axis=0)).ravel()
    if y.sum() != nq or r.sum() != nr:
        raise ArithmeticError('Population pooling did not conserve counts')
    expected = r/r.sum()*y.sum()
    effects = np.log((y+.5)/(expected+.5))
    effects -= np.median(effects)
    profile, selection = stacked_rna_profile(effects, features.bin.to_numpy(), blocks,
                                           gene_segments=gene_segments,
                                           average_segments=average_segments, solver=solver)
    return profile, effects, dict(selection, counts_kind=counts_kind,
        input_scale_interpretation='descriptive_only_no_count_model_calibration',
        query_cells=query.tolist(), reference_cells=reference.tolist(),
        ordered_genes=selected_genes.tolist(), query_feature_total=nq, reference_feature_total=nr,
        count_validation_scope='selected_cells_and_features',
        normalization='reference_gene_share_times_query_selected_feature_total',
        pseudocount=.5, gene_centering='median', canonical_order='features_rows')
