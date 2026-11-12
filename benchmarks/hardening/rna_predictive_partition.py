"""Conditional RNA selection between the retained smoother and a segment model."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.smoothing import predictive_rna_profile

_path = Path(__file__).with_name('performance_candidates.py')
_spec = importlib.util.spec_from_file_location('predictive_partition_model', _path)
_model = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_model)


def medians(values, gene_bins, train, n_bins):
    return np.array([np.median(values[train & (gene_bins == b)]) for b in range(n_bins)])


def segmented(track, blocks):
    reset = np.r_[True, blocks[1:] != blocks[:-1]]
    differences = np.diff(track)[~reset[1:]]
    if not len(differences):
        return None, dict(status='no_within_block_differences')
    scale = float(np.median(np.abs(differences-np.median(differences)))/.6744897501960817/np.sqrt(2))
    if not np.isfinite(scale) or scale <= 0:
        return None, dict(status='positive_working_scale_not_identified')
    fitted, parts = _model.partition(track,np.full(len(track),scale**2),reset)
    return fitted, dict(status='estimated',working_sd=scale,segments=parts,
        penalty='3 log(block bins) per added segment',calibration_status='unvalidated')


def masked_segment(values, gene_bins, blocks, folds, f):
    sizes = np.bincount(gene_bins,minlength=len(blocks))
    holdout = (folds == f) & (sizes[gene_bins] >= 2)
    track = medians(values,gene_bins,~holdout,len(blocks))
    fitted, metadata = segmented(track,blocks)
    return fitted, metadata, holdout


def select(values, gene_bins, blocks):
    values, gene_bins, blocks = np.asarray(values), np.asarray(gene_bins), np.asarray(blocks)
    base, baseline = predictive_rna_profile(values,gene_bins,blocks)
    record = dict(baseline=baseline,selected_model='running_median',
        estimand='descriptive_relative_RNA_profile',calibration_status='unvalidated',
        selection='conditional masked-gene absolute RNA prediction error')
    if baseline['selected_width'] is None:
        return base, dict(record,status='insufficient_masked_gene_support'), base
    full = medians(values,gene_bins,np.ones(len(values),bool),len(blocks))
    candidate, full_metadata = segmented(full,blocks)
    record['full_partition'] = full_metadata
    if candidate is None:
        return base, dict(record,status='segmentation_ineligible'), base
    folds = pd.Series(gene_bins).groupby(gene_bins,sort=False).cumcount().to_numpy()%10
    losses, diagnostics = [], []
    for f in range(10):
        fitted, metadata, holdout = masked_segment(values,gene_bins,blocks,folds,f)
        diagnostics.append(metadata)
        if not holdout.any():
            losses.append(0.)
        elif fitted is None:
            return base, dict(record,status='masked_segmentation_ineligible',folds=diagnostics), base
        else:
            losses.append(float(np.sum(np.abs(fitted[gene_bins[holdout]]-values[holdout]))))
    total = float(np.sum(losses))
    median_risk = np.sum(baseline['fold_error_sums'],axis=0)
    operations = 4*len(values)+32
    unit = np.finfo(float).eps
    allowance = operations*unit/(1-operations*unit)*max(1.,float(median_risk.max()),abs(total))
    choose = bool(total+2*allowance < median_risk.min())
    record.update(status='estimated',selected_model='piecewise_mean_mdl' if choose else 'running_median',
        partition_fold_error_sums=losses,partition_risk_sum=total,
        best_median_risk_sum=float(median_risk.min()),risk_sum_rounding_allowance=float(allowance),
        folds=diagnostics,added_model_selected=choose)
    return candidate if choose else base, record, base
