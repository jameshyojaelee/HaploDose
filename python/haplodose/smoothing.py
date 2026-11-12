"""Conditional RNA profile smoothing, without DNA-state or calibration claims."""
from __future__ import annotations

import numpy as np
import pandas as pd


def rna_profile_on_regions(gene_log_ratios, gene_bins, state_index, blocks):
    """Gene medians on contiguous native RNA-state regions, split at arms/gaps.

    The path defines boundaries only. Its labels do not supply DNA copy numbers
    or independent evidence for the resulting descriptive RNA effects.
    """
    values, bins = np.asarray(gene_log_ratios, dtype=float), np.asarray(gene_bins)
    path, blocks = np.asarray(state_index), np.asarray(blocks)
    if (values.ndim != 1 or bins.shape != values.shape or not len(values)
            or path.ndim != 1 or blocks.shape != path.shape or not len(path)):
        raise ValueError('Aligned gene contrasts/bin IDs and nonempty state/block vectors required')
    if (not np.isfinite(values).all() or bins.dtype.kind not in 'iu'
            or path.dtype.kind not in 'iu' or np.any(path < 0)
            or not np.array_equal(np.unique(bins), np.arange(len(path)))):
        raise ValueError('Finite RNA contrasts, complete integer bin IDs and nonnegative integer states required')
    if (pd.isna(blocks).any() or any(not str(v).strip() for v in blocks)
            or len(pd.unique(blocks)) != 1+np.sum(blocks[1:] != blocks[:-1])):
        raise ValueError('Genomic blocks must be nonmissing, nonempty and contiguous')
    cut = np.r_[True, (path[1:] != path[:-1]) | (blocks[1:] != blocks[:-1])]
    region = np.cumsum(cut)-1
    output = np.empty(len(path))
    for label in range(int(region[-1])+1):
        ix = np.flatnonzero(region == label)
        output[ix] = np.median(values[np.isin(bins, ix)])
    return output, dict(method='gene_median_on_native_RNA_regions', n_regions=int(region[-1])+1,
        units='natural_log_relative_rna_rate', calibration_status='unvalidated',
        boundary_selection='native_full_RNA_HMM_with_arm_and_gap_resets',
        interpretation='descriptive_RNA_not_DNA_copy_number_or_independent_confirmation')


def _running_median(values, blocks, width):
    values, blocks = np.asarray(values), np.asarray(blocks)
    out = np.empty(len(values))
    for block in pd.unique(blocks):
        ix = np.flatnonzero(blocks == block)
        out[ix] = pd.Series(values[ix]).rolling(width, center=True, min_periods=1).median().to_numpy()
    return out


def predictive_rna_profile(gene_log_ratios, gene_bins, blocks):
    """Choose a regional median width by conditional masked-gene prediction.

    Inputs are finite, prepared natural-log RNA contrasts, their zero-based bin
    IDs, and one contiguous arm/gap block label per bin. Every bin needs at least
    one gene. Gene order within bins defines rank-modulo10 prediction folds and
    must be preserved with the gene identities. Upstream feature selection,
    normalization and reference estimation are conditioned on, not held out.

    Returns a descriptive bin profile and selection diagnostics. It supplies no
    DNA copy states, clone assignments, calibrated intervals or detection limit.
    All-singleton inputs retain their observed bin values and explicitly abstain
    from choosing a predictive width.
    """
    values = np.asarray(gene_log_ratios, dtype=np.float64)
    raw_bins, blocks = np.asarray(gene_bins), np.asarray(blocks)
    if values.ndim != 1 or raw_bins.shape != values.shape or blocks.ndim != 1 or not len(blocks):
        raise ValueError('Expected gene contrasts/bin IDs and a nonempty vector of genomic blocks')
    if (not len(values) or not np.isfinite(values).all() or raw_bins.dtype.kind not in 'iu'
            or raw_bins.dtype.kind == 'b'):
        raise ValueError('Finite RNA contrasts and integer bin IDs required')
    if pd.isna(blocks).any() or any(not str(v).strip() for v in blocks):
        raise ValueError('Genomic block labels must be nonmissing and nonempty')
    if len(pd.unique(blocks)) != 1+np.sum(blocks[1:] != blocks[:-1]):
        raise ValueError('Each genomic block must be contiguous; use distinct labels across gaps')
    gene_bins = raw_bins.astype(np.int64)
    if not np.array_equal(np.unique(gene_bins), np.arange(len(blocks))):
        raise ValueError('Bin IDs must cover every genomic bin exactly, starting at zero')
    full = np.array([np.median(values[gene_bins == b]) for b in range(len(blocks))])
    largest = max(int(np.sum(blocks == block)) for block in pd.unique(blocks))
    widths = [1]+[2**k+1 for k in range(1, largest.bit_length()) if 2**k+1 <= largest]
    folds = pd.Series(gene_bins).groupby(gene_bins, sort=False).cumcount().to_numpy()%10
    sizes = np.bincount(gene_bins, minlength=len(blocks))
    eligible = sizes[gene_bins] >= 2
    metadata = dict(method='conditional_masked_gene_median', units='natural_log_relative_rna_rate',
        calibration_status='unvalidated', interpretation='descriptive_RNA_not_DNA_copy_number',
        widths=widths, selected_width=None, loss='absolute RNA log-rate error',
        ineligible_singleton_genes=int((~eligible).sum()))
    if not eligible.any():
        return full, dict(metadata, status='insufficient_masked_gene_support',
                          held_out_genes=[0]*10, mean_prediction_error=None)
    losses = np.zeros((10, len(widths)))
    held_out_n = np.zeros(10, int)
    for f in range(10):
        holdout = (folds == f) & eligible
        if not holdout.any():
            continue
        # Singleton genes are never prediction targets, so do not remove their
        # only observation from neighbouring-bin training profiles.
        train = ~holdout
        medians = np.array([np.median(values[train & (gene_bins == b)])
                           if np.any(train & (gene_bins == b)) else np.nan for b in range(len(blocks))])
        predicted = np.array([_running_median(medians, blocks, w) for w in widths])
        target_predictions = predicted[:, gene_bins[holdout]]
        if not np.isfinite(target_predictions).all():
            raise ArithmeticError('An eligible masked gene has no training prediction')
        losses[f] = np.sum(np.abs(target_predictions-values[holdout][None]), axis=1)
        held_out_n[f] = int(holdout.sum())
    if held_out_n.sum() != eligible.sum():
        raise ArithmeticError('Masked gene coverage is not conserved')
    total = losses.sum(axis=0)
    eps = np.finfo(float).eps
    operations = 4*len(values)+32
    allowance = operations*eps/(1-operations*eps)*max(1., float(total.max()))
    tied = np.flatnonzero(total <= total.min()+2*allowance)
    index = int(tied[0])
    fitted = _running_median(full, blocks, widths[index])
    return fitted, dict(metadata, status='estimated', selected_width=widths[index],
        mean_prediction_error=(total/eligible.sum()).tolist(), fold_error_sums=losses.tolist(),
        held_out_genes=held_out_n.tolist(), risk_sum_rounding_allowance=float(allowance),
        numerically_tied_widths=[widths[i] for i in tied])
