"""Frozen candidate evaluation with explicit cell and feature roles.

This module enforces computational role isolation. Disjoint counts are not a proof of statistical
independence: composition, biological state, fitted discovery quantities and correlated genes remain
load-bearing assumptions. Returned contrasts are nominal diagnostics, not validated DNA clone calls.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib

import numpy as np
from scipy.special import gammaln, logsumexp, ndtr, xlogy


def _ids(cell_ids):
    ids = np.asarray(cell_ids, dtype=str)
    if ids.ndim != 1 or any(not value.strip() for value in ids) or len(set(ids)) != len(ids):
        raise ValueError('cell IDs must be one-dimensional, unique and nonempty')
    return ids


def assign_roles(cell_ids, *, discovery_fraction, seed, min_discovery, min_confirmation):
    """Hash each ID to a role; existing roles do not change with row order or additional IDs.

    The fraction and minimum counts are required scientific design choices. A hash allocation can fail
    the minimums; it is never silently rebalanced using outcomes or counts.
    """
    ids = _ids(cell_ids)
    if not np.isfinite(discovery_fraction) or not 0 < discovery_fraction < 1:
        raise ValueError('discovery_fraction must be strictly between 0 and 1')
    if any(not isinstance(v, (int, np.integer)) or v < 1 for v in (min_discovery, min_confirmation)):
        raise ValueError('minimum discovery and confirmation counts must be positive integers')
    threshold = int(discovery_fraction * 2 ** 64)
    discovery = np.array([int.from_bytes(hashlib.blake2b(
        f'{int(seed)}\0{value}'.encode(), digest_size=8).digest(), 'big') < threshold for value in ids])
    if discovery.sum() < min_discovery or (~discovery).sum() < min_confirmation:
        raise ValueError('ID-based role allocation does not meet the requested minimum cell counts')
    return dict(cell_ids=ids.copy(), discovery=discovery, confirmation=~discovery,
                assignment='blake2b_id_threshold', seed=int(seed), discovery_fraction=float(discovery_fraction))


@dataclass(frozen=True)
class Contrast:
    """A discovery-frozen candidate pair, feature region and expected effect direction (+1 or -1)."""
    candidate_a: int
    candidate_b: int
    features: tuple[int, ...]
    direction: int
    name: str

    def __post_init__(self):
        object.__setattr__(self, 'features', tuple(self.features))


@dataclass(frozen=True)
class FrozenDiscovery:
    log_rates: np.ndarray                 # candidates x features; arbitrary additive row centres allowed
    phi: np.ndarray                       # NB dispersion per feature
    proposal_mask: np.ndarray
    test_mask: np.ndarray
    contrasts: tuple[Contrast, ...]
    discovery_cell_ids: tuple[str, ...]
    depth_gamma: np.ndarray | None = None
    proposal_log_depth_centre: float = 0.
    programme_loadings: np.ndarray | None = None   # programmes x features
    programme_precision: np.ndarray | None = None
    candidate_prior: np.ndarray | None = None

    def __post_init__(self):
        def frozen(name, value, dtype=float):
            array = np.array(value, dtype=dtype, copy=True)
            array.setflags(write=False)
            object.__setattr__(self, name, array)
            return array
        lr = frozen('log_rates', self.log_rates)
        phi = frozen('phi', self.phi)
        prop = frozen('proposal_mask', self.proposal_mask, bool)
        test = frozen('test_mask', self.test_mask, bool)
        if lr.ndim != 2 or not lr.shape[0] or not lr.shape[1] or not np.isfinite(lr).all():
            raise ValueError('log_rates must be a finite nonempty candidate x feature matrix')
        K, J = lr.shape
        if phi.shape != (J,) or not np.isfinite(phi).all() or (phi <= 0).any():
            raise ValueError('phi must contain a finite positive dispersion per feature')
        if prop.shape != (J,) or test.shape != (J,) or not prop.any() or not test.any() or (prop & test).any():
            raise ValueError('proposal and test masks must be nonempty, correctly shaped and disjoint')
        ids = tuple(_ids(self.discovery_cell_ids))
        object.__setattr__(self, 'discovery_cell_ids', ids)
        gamma = frozen('depth_gamma', np.zeros(J) if self.depth_gamma is None else self.depth_gamma)
        W = frozen('programme_loadings', np.zeros((0, J)) if self.programme_loadings is None
                   else self.programme_loadings)
        precision = frozen('programme_precision', np.ones(len(W)) if self.programme_precision is None
                           else self.programme_precision)
        prior = frozen('candidate_prior', np.ones(K) / K if self.candidate_prior is None else self.candidate_prior)
        if (gamma.shape != (J,) or not np.isfinite(gamma).all() or not np.isfinite(self.proposal_log_depth_centre)
                or W.ndim != 2 or W.shape[1] != J or not np.isfinite(W).all()
                or precision.shape != (len(W),) or not np.isfinite(precision).all() or (precision <= 0).any()
                or prior.shape != (K,) or not np.isfinite(prior).all() or (prior <= 0).any()):
            raise ValueError('invalid frozen depth, programme or candidate-prior parameters')
        frozen('candidate_prior', prior / prior.sum())
        contrasts = tuple(self.contrasts)
        for contrast in contrasts:
            features = np.asarray(contrast.features)
            if (contrast.direction not in (-1, 1) or not 0 <= contrast.candidate_a < K
                    or not 0 <= contrast.candidate_b < K or contrast.candidate_a == contrast.candidate_b
                    or features.ndim != 1 or not len(features) or not np.issubdtype(features.dtype, np.integer)
                    or (features < 0).any() or (features >= J).any() or len(set(features)) != len(features)
                    or not test[features].any() or not (test & ~np.isin(np.arange(J), features)).any()):
                raise ValueError('each contrast needs distinct valid candidates and a test region plus test background')
        object.__setattr__(self, 'contrasts', contrasts)


def contrast_diagnostic(contrast, sums, variance, n_cells):
    """Summarize additive sufficient diagnostics; permits bounded confirmation chunks."""
    sums, variance = np.asarray(sums, float), np.asarray(variance, float)
    testable = bool(all(n_cells) and (sums > 0).all())
    effect = float(np.log((sums[0, 0] + .5) / (sums[0, 1] + .5))
                   - np.log((sums[1, 0] + .5) / (sums[1, 1] + .5))) if testable else None
    se = float(np.sqrt((variance / (sums + .5) ** 2).sum())) if testable else None
    z = contrast.direction * effect / se if testable and se > 0 else None
    return dict(name=contrast.name, candidate_a=contrast.candidate_a, candidate_b=contrast.candidate_b,
                direction=contrast.direction, n_cells=list(n_cells), test_counts=sums.tolist(),
                test_variance=variance.tolist(), testable=testable, log_odds_contrast=effect,
                approximate_se=se, signed_z=z, nominal_one_sided_p=float(ndtr(-z)) if z is not None else None,
                multiple_testing='none; nominal diagnostics')


def score_confirmation(Y, cell_ids, artifact: FrozenDiscovery):
    """Score fixed candidates using proposal counts only, then summarize frozen test contrasts.

    Proposal-normalized NB scores are a composite likelihood, not an exact conditional count model.
    Programme scores are fixed-ridge approximations fitted separately for each frozen candidate on
    proposal features only. No parameters are written back and no candidate is added or selected from
    test evidence. Test contrasts use a delta-method NB variance without cross-feature covariance;
    their one-sided p-values are nominal, unadjusted and not an independent-confirmation guarantee.
    """
    ids = _ids(cell_ids)
    if set(ids) & set(artifact.discovery_cell_ids):
        raise ValueError('discovery cells cannot also serve as confirmation cells')
    Y = np.asarray(Y, dtype=float)
    if (Y.shape != (len(ids), artifact.log_rates.shape[1]) or not np.isfinite(Y).all()
            or (Y < 0).any() or (np.abs(Y - np.rint(Y)) > 1e-6).any()):
        raise ValueError('confirmation counts must be finite nonnegative integers with aligned IDs and features')
    prop, test = artifact.proposal_mask, artifact.test_mask
    P = Y[:, prop]
    L = P.sum(axis=1)
    x = np.log(np.maximum(L, 1.)) - artifact.proposal_log_depth_centre
    K = len(artifact.log_rates)
    scores = np.full((len(Y), K), -np.inf)
    U = np.zeros((len(Y), len(artifact.programme_loadings)))
    W = artifact.programme_loadings[:, prop]
    size = 1 / artifact.phi[prop]
    for c in np.flatnonzero((L > 0) & np.isfinite(L)):
        best_score = -np.inf
        for k in range(K):
            u = np.zeros(len(W))
            eta = artifact.log_rates[k, prop] + artifact.depth_gamma[prop] * x[c]
            if len(W):
                mean0 = L[c] * np.exp(eta - logsumexp(eta))
                residual = np.log((P[c] + .5) / (mean0 + .5))
                weight = 1 / (artifact.phi[prop] + 1 / (P[c] + .5))
                u = np.linalg.solve((W * weight) @ W.T + np.diag(artifact.programme_precision),
                                    W @ (weight * residual))
                eta = eta + u @ W
            mean = L[c] * np.exp(eta - logsumexp(eta))
            logp = (gammaln(P[c] + size) - gammaln(size) - gammaln(P[c] + 1)
                    + size * (np.log(size) - np.log(size + mean)) + xlogy(P[c], mean / (size + mean)))
            value = (logp.sum() + np.log(artifact.candidate_prior[k])
                     - .5 * np.sum(artifact.programme_precision * u ** 2))
            scores[c, k] = value if np.isfinite(value) else -np.inf
            if scores[c, k] > best_score:
                best_score = scores[c, k]
                U[c] = u
    supported = (L > 0) & np.isfinite(L) & np.isfinite(scores).any(axis=1)
    assigned = np.where(supported, scores.argmax(axis=1), -1)
    unresolved_reason = np.where(supported, '', np.where(L <= 0, 'no_proposal_depth', 'no_finite_candidate_score'))
    diagnostics = []
    for contrast in artifact.contrasts:
        region = test & np.isin(np.arange(Y.shape[1]), contrast.features)
        background = test & ~region
        groups = [Y[assigned == k] for k in (contrast.candidate_a, contrast.candidate_b)]
        sums = np.array([[g[:, region].sum(), g[:, background].sum()] for g in groups])
        variance = np.array([[(g[:, mask] + artifact.phi[mask] * g[:, mask] ** 2).sum()
                              for mask in (region, background)] for g in groups])
        diagnostics.append(contrast_diagnostic(contrast, sums, variance, [len(g) for g in groups]))
    return dict(cell_ids=ids, candidate=assigned, unresolved_reason=unresolved_reason,
                proposal_log_scores=scores, programme_scores=U,
                contrasts=diagnostics, evidence=dict(role_isolation='discovery_ids_excluded',
                    membership_features='proposal_only', nuisance_features='proposal_only',
                    contrast_features='test_only', statistical_independence='not_established',
                    specimen_error_control=False, biological_clone_validation=False))
