"""Bulk RNA mode: one copy-number profile per sample against a panel of normal samples.

Each bulk sample is a single pooled profile, so there is no clone model. The steps are those of a pooled
clone profile, with the noise structure of bulk RNA: counts are deep (Poisson noise is small), and the
main error is gene-level expression that differs from the reference.

1. Bins of consecutive expressed genes (`bins.make_bins`), bin counts per sample.
2. Reference rate mu_j = pooled normal counts. Per-bin reference variance v_j = robust variance of the
   normals' log ratios to mu (how much a bin varies between normal samples).
3. Expression programmes: top singular vectors of the normals' log ratios, high-passed within arms
   (`programmes.high_pass`, window 0), so they can never explain a whole-arm change. Each tumour gets
   free programme amplitudes, fitted around its current dosage profile and removed before segmentation.
4. `profile.fit_pooled` (centre, dosage slope beta, extra sd tau with v_j added) and the four-state HMM
   (`profile.segment`). beta absorbs tumour purity times dosage response, so its lower bound is lower
   than for single-cell clones.
5. Arm calls calibrated on the normal panel: each normal is profiled leave-one-out (reference rate, bin
   variance and programmes all refitted without it, so its own noise never enters its null); a sample's arm median log ratio becomes z = (value - normal median) / (normal MAD sd for that
   arm x the sample's tau over the normals' median tau, at least 1), called gain or loss at |z| > ARM_Z. The HMM segments stay as the finer-scale output.

Outputs relative RNA dosage, not DNA copy number: a uniform genome doubling is invisible to expression.
"""
from __future__ import annotations

from dataclasses import asdict, replace

import warnings

import numpy as np
import pandas as pd
from scipy import sparse

from .bins import arm_labels, arm_starts, chromosome_starts, make_bins
from .io import new_output_dir, read_gene_list, read_positions, write_json, write_table
from .profile import ProfileParams, fit_pooled, log_ratio, segment, segments_table
from .programmes import fit_amplitudes, high_pass

BULK_PROFILE = ProfileParams(beta_bounds=(0.05, 1.0))
BETA_GRID = (0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0)
ARM_Z = 3.0
# Below this many normals the arm z scale (median/MAD over the normals, treated as known) is not a calibrated
# null: for independent Gaussian panels, the fraction of a null sample beyond |z| = 3 was 30.6% with 3 normals,
# 5.3% with 10, 1.3% with 30 and 0.50% with 100 (docs/efficiency_review_2026-11-02.md, 50,000 repetitions;
# not a measured pipeline false-call rate). run.json records it; calls are not suppressed.
SMALL_PANEL = 10


def reference_model(Yn, n_programmes, arm_reset):
    """Pooled normal rate mu [J], per-bin variance v [J] and programme loadings W [F, J] from normals Yn [N, J]."""
    Ln = Yn.sum(axis=1)
    mu = Yn.sum(axis=0) / Ln.sum()
    R = log_ratio(Yn, Ln[:, None] * mu[None])
    R -= np.median(R, axis=1, keepdims=True)
    med = np.median(R, axis=0)
    mad = 1.4826 * np.median(np.abs(R - med), axis=0)
    v = np.maximum(mad ** 2 - np.median(1.0 / (Yn + 0.5), axis=0), 0.0)
    W = np.zeros((0, Yn.shape[1]))
    if n_programmes > 0 and len(Yn) > 2:
        _, sv, Vt = np.linalg.svd(R - R.mean(axis=0), full_matrices=False)
        W = high_pass(Vt[:min(n_programmes, len(sv))], arm_reset, 0)
    return mu, v, W


def start_beta(y, E, v, reset, params: ProfileParams, grid=BETA_GRID):
    """Dosage slope with the highest HMM likelihood on a grid. fit_pooled only updates beta from bins it
    already calls aberrant, so in low-purity bulk a start far above the true slope calls nothing and stays."""
    _, _, tau = fit_pooled(y, E, reset, params, rounds=0, v=v)
    best, best_ll = params.beta, -np.inf
    for b in grid:
        if not params.beta_bounds[0] <= b <= params.beta_bounds[1]:
            continue
        c, _, t = fit_pooled(y, E, reset, replace(params, beta=b, beta_bounds=(b, b)), rounds=2, v=v)
        ll = float(segment(y[None], E[None], reset, [c], b, np.sqrt(t ** 2 + v)[None], params, decode=False)[3].sum())
        if ll > best_ll:
            best, best_ll = b, ll
    return replace(params, beta=float(best))


def fit_sample(y, mu, v, W, reset, params: ProfileParams, rounds=2):
    """Programme shift, pooled fit and HMM for one sample y [J] against reference rate mu [J]."""
    E0 = y.sum() * mu
    E, shift = E0, np.zeros_like(E0)
    params = start_beta(y, E, v, reset, params)
    centre, beta, tau = fit_pooled(y, E, reset, params, v=v)
    for _ in range(rounds if len(W) else 0):
        sd = np.sqrt(tau ** 2 + v)
        _, _, dose, _ = segment(y[None], E[None], reset, [centre], beta, sd[None], params, decode=False)
        w = 1.0 / (sd ** 2 + 1.0 / (y + 0.5))
        shift = fit_amplitudes(log_ratio(y, E0) - centre - dose[0], w, W)[0]
        E = E0 * np.exp(shift)
        E *= E0.sum() / E.sum()
        centre, beta, tau = fit_pooled(y, E, reset, params, v=v)
    sd = np.sqrt(tau ** 2 + v)
    post, path, dose, _ = segment(y[None], E[None], reset, [centre], beta, sd[None], params)
    return dict(E=E, shift=shift, centre=centre, beta=beta, tau=tau, sd=sd, post=post[0], path=path[0], dose=dose[0])


def median_pseudo_counts(G, gbin, onehot, normal):
    """bin_stat 'median': bin value = median over its genes of per-gene log ratios to the pooled normals (columns
    `normal` of G, genes x samples), returned as pseudo counts with the same bin totals: E_bin * exp(median).
    Each gene counts once, so a few highly expressed or tumour-specific genes cannot dominate a bin."""
    L = G.sum(axis=0)
    mu_g = G[:, normal].sum(axis=1) / L[normal].sum()
    lr = np.log((G + 0.5) / (L[None] * mu_g[:, None] + 0.5))
    med = np.column_stack([np.median(lr[gbin == j], axis=0) for j in range(onehot.shape[1])])
    return (L[:, None] * np.asarray(onehot.T @ mu_g).ravel()[None]) * np.exp(med)


def run_bulk(counts, samples, positions, out, exclude_genes=None, bin_genes=10, min_cpm=1.0, n_programmes=10,
             params: ProfileParams = BULK_PROFILE, bin_stat='sum'):
    """counts: TSV genes x samples (first column gene); samples: TSV with sample, role (tumour | normal).
    bin_stat: 'sum' (pooled counts per bin) or 'median' (median per-gene log ratio, as pseudo counts)."""
    out = new_output_dir(out)
    meta = pd.read_csv(samples, sep='\t', dtype=str)
    X = pd.read_csv(counts, sep='\t', index_col=0)
    X = X.loc[:, meta['sample']]
    normal = meta.role.eq('normal').to_numpy()
    if normal.sum() < 3:
        raise ValueError('need at least 3 normal samples')
    if normal.sum() < SMALL_PANEL:
        warnings.warn(f'{int(normal.sum())} normal samples: arm calls (|z| > {ARM_Z}) use a median/MAD scale from '
                      f'the normals as if known; below {SMALL_PANEL} normals they are not calibrated (null '
                      'exceedance about 31% at 3 normals, 5% at 10 in a Gaussian simulation). Use the z values '
                      'as ranks, not calls.')
    cpm = X.to_numpy(np.float64) / X.to_numpy(np.float64).sum(axis=0, keepdims=True) * 1e6
    bins, gene_bin, used = make_bins(X.index, read_positions(positions), cpm[:, normal].mean(axis=1),
                                     exclude=read_gene_list(exclude_genes), bin_genes=bin_genes, min_mean=min_cpm)
    J = len(bins)
    keep = gene_bin >= 0
    n_keep = int(keep.sum())
    onehot = sparse.csr_matrix((np.ones(n_keep), (np.arange(n_keep), gene_bin[keep])), shape=(n_keep, J))
    G = X.to_numpy(np.float64)[keep]                      # genes x samples
    Y = np.asarray((onehot.T @ G).T)                      # samples x bins
    if bin_stat == 'median':
        Y = median_pseudo_counts(G, gene_bin[keep], onehot, normal)
    reset, arm_reset, arm = chromosome_starts(bins), arm_starts(bins), arm_labels(bins)
    mu, v, W = reference_model(Y[normal], n_programmes, arm_reset)
    normal_rows = np.flatnonzero(normal)

    seg_rows, arm_rows, sample_rows = [], [], []
    lr_out = np.full((len(meta), J), np.nan, np.float32)
    for i, s in enumerate(meta['sample']):
        # Normals are profiled against a reference fitted without them (their own null distribution). In median
        # mode the pseudo counts themselves depend on the normals, so they are rebuilt without the held-out one.
        Yi = Y
        if normal[i] and bin_stat == 'median':
            Yi = median_pseudo_counts(G, gene_bin[keep], onehot, normal & (np.arange(len(meta)) != i))
        ref = reference_model(Yi[normal_rows[normal_rows != i]], n_programmes, arm_reset) if normal[i] else (mu, v, W)
        f = fit_sample(Yi[i], *ref, reset, params)
        # Reported log ratios are centred by the plain median over bins, not the HMM fit's centre: the fit
        # centre moves with the sample's own calls, which made arm values less comparable across samples
        # (U-CAN, 1,034 tumours: arm AUC loss 0.854 -> 0.904, gain 0.871 -> 0.887; per-sample rho unchanged).
        r = log_ratio(Yi[i], f['E'])
        r -= np.median(r)
        lr_out[i] = r
        for row in segments_table(f['post'], f['path'], f['dose'], Yi[i], f['E'],
                                  reset, f['beta'], float(np.median(f['sd'])), params):
            a, b = row['first_bin'], row['last_bin']
            seg_rows.append(dict(sample=s, chromosome=bins.chromosome.iloc[a], start=int(bins.start.iloc[a]),
                                 end=int(bins.end.iloc[b]), **row))
        frame = pd.DataFrame({'arm': arm, 'r': r, 'dose': f['dose'], 'state': f['path'] + 1})
        for a_name, g in frame.groupby('arm', sort=False):
            arm_rows.append(dict(sample=s, arm=a_name, n_bins=len(g), median_log_ratio=float(g.r.median()),
                                 mean_log_dose=float(g.dose.mean()),
                                 gain_fraction=float(g.state.ge(3).mean()), loss_fraction=float(g.state.eq(1).mean())))
        sample_rows.append(dict(sample=s, role=meta.role.iloc[i], binned_counts=float(Yi[i].sum()), centre=f['centre'],
                                median_centre=float(np.median(log_ratio(Yi[i], f['E']))),
                                beta=f['beta'], tau=f['tau'], programme_shift_rms=float(np.sqrt((f['shift'] ** 2).mean()))))
    arms_t = pd.DataFrame(arm_rows)
    null = arms_t[arms_t['sample'].isin(meta['sample'][normal])].groupby('arm').median_log_ratio
    centre_n = null.median()
    sd_n = null.apply(lambda x: 1.4826 * np.median(np.abs(x - x.median()))).clip(lower=1e-3)
    # Tumours deviate from normals more than normals do from each other: scale the normal-panel sd by the
    # sample's own bin-level noise (tau) relative to the normals' median tau, never below 1.
    tau = pd.DataFrame(sample_rows).set_index('sample').tau
    inflate = (tau / tau[meta['sample'][normal]].median()).clip(lower=1.0)
    arms_t['z'] = ((arms_t.median_log_ratio - arms_t.arm.map(centre_n))
                   / (arms_t.arm.map(sd_n) * arms_t['sample'].map(inflate)))
    arms_t['call'] = np.select([arms_t.z > ARM_Z, arms_t.z < -ARM_Z], ['gain', 'loss'], 'neutral')
    write_table(bins.assign(arm=arm, reference_sd=np.sqrt(v)), out / 'bins.tsv')
    write_table(used, out / 'genes.tsv.gz')
    write_table(pd.DataFrame(seg_rows), out / 'segments.tsv')
    write_table(arms_t, out / 'arms.tsv')
    write_table(pd.DataFrame(sample_rows), out / 'samples.tsv')
    np.save(out / 'log_ratio.npy', lr_out)
    write_json(dict(mode='bulk', n_samples=len(meta), n_normals=int(normal.sum()), n_bins=J,
                    n_genes_used=int(keep.sum()), n_programmes=len(W), bin_stat=bin_stat, min_cpm=min_cpm, bin_genes=bin_genes,
                    params=asdict(params), arm_z=ARM_Z, beta_grid=BETA_GRID,
                    arm_calls_calibrated=bool(normal.sum() >= SMALL_PANEL), exclude_genes=str(exclude_genes)),
               out / 'run.json')
    return out

