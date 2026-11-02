"""Expression programmes of normal cells, used as per-cell nuisance terms.

Normal cells carry co-expression programmes (immune subtypes, cell state, gene
clusters such as HLA) whose genes can sit in neighbouring bins, so in single
cells they look like copy-number change. Their loadings W are learned on fit
reference cells only, which carry no copy-number change; held-out normals test
the result. A clone profile that lies largely in the programme span is reported
(`programme_overlap`), since it may be expression rather than dosage.

Loadings are high-passed along the genome: within each segment between `reset`
boundaries, either a running mean over `window` bins is removed, or (window 0)
the segment mean. With arm boundaries and window 0, programmes cannot explain a
whole-arm change, so arm-level copy-number change is always left to dosage.
Shorter events compete with programmes that match their shape; that is the
stated resolution limit. The converse limit: the arm-mean part of a programme
(its width over the arm width times its log fold) stays unexplained, so on short
arms (e.g. chr21q, about 20 bins) a strong programme can still look like dosage.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .kernels.estep import programme_scores


@dataclass
class Programmes:
    W: np.ndarray            # [F, J] log-rate loadings, high-passed
    prior_prec: np.ndarray   # [F] precision of the per-cell scores
    phi: np.ndarray          # [J]
    depth: float             # median binned UMIs of the fit sample
    variance_share: float    # share of fit-normal residual variance in the top F components

    def overlap(self, dose, base):
        """Share of a log-dose profile's Fisher information (one cell at median fit depth) the programmes explain."""
        m = self.depth * base
        w = m / (1.0 + self.phi * m)
        shift = fit_amplitudes(dose, w, self.W)[0]
        total = float((w * dose ** 2).sum())
        return float((w * shift ** 2).sum() / total) if total > 0 else np.nan


def high_pass(W, reset, window):
    """Remove a centred running mean (window bins, truncated at segment ends) from each row; window 0 removes the segment mean."""
    out = np.array(W, np.float64)
    starts = np.flatnonzero(reset)
    ends = np.append(starts[1:], len(reset))
    half = window // 2
    for a, b in zip(starts, ends):
        seg = out[:, a:b]
        if window <= 0:
            seg -= seg.mean(axis=1, keepdims=True)
            continue
        c = np.concatenate([np.zeros((len(seg), 1)), np.cumsum(seg, axis=1)], axis=1)
        idx = np.arange(b - a)
        lo, hi = np.maximum(idx - half, 0), np.minimum(idx + half + 1, b - a)
        out[:, a:b] = seg - (c[:, hi] - c[:, lo]) / (hi - lo)
    return out


def learn_programmes(Y, L, mu, phi, rows_by_type, n_programmes, reset, seed, window=0, per_type=1500):
    """Top residual programmes of fit normal cells, each scored on its own type baseline.

    mu [J, T] reference rates; rows_by_type: T arrays of fit-cell rows (aligned with mu).
    The score prior variance is the mean squared flat-prior score on the same cells
    (includes sampling noise, so shrinkage is mild)."""
    rng = np.random.default_rng(seed)
    samples, resid = [], []
    for t, rows in enumerate(rows_by_type):
        sub = np.sort(rng.choice(rows, min(len(rows), per_type), replace=False))
        Ys = np.asarray(Y[sub], np.float64)
        m = L[sub, None] * mu[None, :, t]
        resid.append((Ys - m) / np.sqrt(m + phi[None] * m ** 2))
        samples.append((t, sub, Ys))
    Rz = np.vstack(resid)
    Rz -= Rz.mean(axis=0)
    # Only the top n_programmes components are kept, so a seeded randomized SVD replaces the full one.
    from sklearn.utils.extmath import randomized_svd
    F = min(n_programmes, *Rz.shape)
    _, sv, Vt = randomized_svd(Rz, F, n_oversamples=20, n_iter=7, random_state=seed)
    basis = Vt[:F]
    weights = mu.mean(axis=1)
    depth = float(np.median(np.concatenate([L[sub] for _, sub, _ in samples])))
    W = high_pass(basis * np.sqrt(1.0 / (depth * weights) + phi)[None], reset, window)
    U = [programme_scores(Ys.astype(np.float32), L[sub], np.ones((len(sub), 1), np.float32),
                          np.log(mu[:, t])[None], phi, W, np.full(F, 1e-6))
         for t, sub, Ys in samples]
    var = np.maximum(np.mean(np.vstack(U) ** 2, axis=0), 1e-6)
    return Programmes(W=W, prior_prec=1.0 / var, phi=np.asarray(phi, np.float64), depth=depth,
                      variance_share=float((sv[:F] ** 2).sum() / (Rz ** 2).sum()))


def fit_amplitudes(r, w, W, ridge=1e-6):
    """Programme shift W^T beta fitted to log-ratio profiles r [G, J] by weighted least squares (w [G, J])."""
    r, w = np.atleast_2d(r), np.atleast_2d(w)
    out = np.empty_like(r, dtype=np.float64)
    for g in range(len(r)):
        A = (W * w[g][None]) @ W.T
        A[np.diag_indices_from(A)] += ridge * np.trace(A) / len(A)
        out[g] = W.T @ np.linalg.solve(A, (W * w[g][None]) @ r[g])
    return out
