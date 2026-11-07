"""Cell x clone log-likelihoods from per-cell bin counts.

Clone k predicts cell c's bin counts as L_c * rate_kj with rates summing to 1,
so only the compositional (dosage) differences between clones matter. Negative
binomial with bin dispersion phi_j; terms constant across clones are dropped.

Optional per-cell nuisance terms enter the log rate before renormalization:
depth, x_c * gamma_j, and expression programmes, sum_f u_cf * W_fj, with W
learned on normal cells (see `haplodose.programmes`).
"""
from __future__ import annotations

import warnings

import jax
import jax.numpy as jnp
import numpy as np


@jax.jit
def _chunk_loglik(Y, L, x, U, rate0, delta, gamma, W, inv_phi):
    # NB log-likelihood of each clone relative to clone 0, per bin and then summed. Summing each clone's large
    # absolute terms in float32 and subtracting later lost small contrasts (2,000 bins, a +0.01 log-rate
    # difference in 20 bins: float64 difference -0.005, float32 0.0). rate0 [J] is clone 0's log rate and
    # delta [K, J] = log_rate_k - log_rate_0 (formed in float64 by the caller, so small contrasts stay exact).
    # With s0 clone 0's normalised composition, clone k's normalised log-rate difference is
    # d = delta_k - log1p(sum_j s0_j expm1(delta_kj)); per bin, with r = 1 / phi,
    # Y d - (Y + r) log((r + m_k) / (r + m_0)) and log((r + m_k) / (r + m_0)) = log1p(m_0 / (r + m_0) expm1(d)).
    off = x[:, None] * gamma[None] + U @ W                                # n, J
    lr0 = rate0[None] + off
    s0 = jnp.exp(lr0 - jax.scipy.special.logsumexp(lr0, axis=-1, keepdims=True))   # n, J
    shift = jnp.log1p(s0 @ jnp.expm1(delta).T)                            # n, K
    d = delta[None] - shift[:, :, None]                                   # n, K, J
    m0 = (L[:, None] * s0)[:, None, :]
    return jnp.sum(Y[:, None, :] * d - (Y[:, None, :] + inv_phi) * jnp.log1p(m0 / (inv_phi + m0) * jnp.expm1(d)), axis=-1)


@jax.jit
def _chunk_scores_fit(Y, L, x, R, log_rate, gamma, W, inv_phi, prior_prec):
    """Fisher scoring of sum_k R_ck log NB(Y_c | clone k, u_c) plus the score prior.

    Clone accumulation uses n x J temporaries, not n x K x J. Backtracking
    compares stable NB log-likelihood ratios, so large shared count terms do not
    determine acceptance. Iteration limits are resource safeguards; reaching one
    never reports convergence. The step tolerance is tied to float32 precision.
    """
    n, J = Y.shape
    F = W.shape[0]
    phi = 1.0 / inv_phi
    eps = np.finfo(np.float32).eps
    step_tol = np.sqrt(eps)
    offset = x[:, None] * gamma[None]

    def composition(u, k):
        lr = log_rate[k][None] + offset + u @ W
        return jax.nn.softmax(lr, axis=-1)

    def derivatives(u):
        def add(k, accum):
            g, H = accum
            s = composition(u, k)
            m = L[:, None] * s
            a = (Y - m) / (1.0 + phi[None] * m)
            w = m / (1.0 + phi[None] * m)
            Wbar, Ww = s @ W.T, w @ W.T
            grad = a @ W.T - Wbar * a.sum(axis=1, keepdims=True)
            fisher = ((W[None] * w[:, None, :]) @ W.T
                      - Wbar[:, :, None] * Ww[:, None, :] - Ww[:, :, None] * Wbar[:, None, :]
                      + Wbar[:, :, None] * Wbar[:, None, :] * w.sum(axis=1)[:, None, None])
            return g + R[:, k, None] * grad, H + R[:, k, None, None] * fisher
        return jax.lax.fori_loop(0, log_rate.shape[0], add,
                                (-u * prior_prec[None], jnp.broadcast_to(jnp.diag(prior_prec), (n, F, F))))

    def improvement(u, du):
        delta = du @ W
        def add(k, accum):
            value, magnitude = accum
            s = composition(u, k)
            d = delta - jnp.log1p(jnp.sum(s * jnp.expm1(delta), axis=1))[:, None]
            m = L[:, None] * s
            a = Y * d
            b = (Y + inv_phi) * jnp.log1p(m / (inv_phi + m) * jnp.expm1(d))
            return (value + R[:, k] * jnp.sum(a - b, axis=1),
                    magnitude + R[:, k] * jnp.sum(jnp.abs(a) + jnp.abs(b), axis=1))
        value, magnitude = jax.lax.fori_loop(0, log_rate.shape[0], add, (jnp.zeros(n), jnp.zeros(n)))
        penalty = jnp.sum(prior_prec[None] * du * (u + 0.5 * du), axis=1)
        # Standard gamma_N reduction bound, with elementary per-bin arithmetic.
        gamma_n = (J + 16) * eps / jnp.maximum(1.0 - (J + 16) * eps, eps)
        error = gamma_n * (magnitude + jnp.sum(jnp.abs(prior_prec[None] * du * (u + 0.5 * du)), axis=1))
        return value - penalty, error

    def iterate(state):
        u, active, converged, failed, iterations, gain = state
        g, H = derivatives(u)
        du = jnp.linalg.solve(H, g[..., None])[..., 0]
        finite = jnp.isfinite(du).all(axis=1)
        # Bound the initial step for the local approximation; acceptance is still
        # decided by the actual responsibility-weighted objective, not this cap.
        size = jnp.max(jnp.abs(du @ W), axis=1)
        small = size <= step_tol * (1.0 + jnp.max(jnp.abs(u @ W), axis=1))
        converged = converged | (active & finite & small)
        working = active & finite & ~small
        direction = jnp.where(working[:, None], du, 0.0)
        direction *= jnp.minimum(1.0, 0.5 / jnp.maximum(size, eps))[:, None]

        def searching(search):
            count, scale, accepted, _, _ = search
            return (count < np.finfo(np.float32).nmant + 2) & jnp.any(working & ~accepted)

        def shrink(search):
            count, scale, accepted, accepted_du, accepted_gain = search
            trial = direction * scale[:, None]
            delta_q, error = improvement(u, trial)
            accept = working & ~accepted & jnp.isfinite(delta_q) & (delta_q >= -error)
            accepted_du = jnp.where(accept[:, None], trial, accepted_du)
            accepted_gain = jnp.where(accept, delta_q, accepted_gain)
            return count + 1, jnp.where(accepted | accept, scale, 0.5 * scale), accepted | accept, accepted_du, accepted_gain

        _, _, accepted, move, delta_q = jax.lax.while_loop(
            searching, shrink, (0, jnp.ones(n), jnp.zeros(n, bool), jnp.zeros_like(u), jnp.zeros(n)))
        failed = failed | (active & (~finite | (working & ~accepted)))
        iterations += active.astype(jnp.int32)
        return u + move, working & accepted, converged, failed, iterations, gain + delta_q

    initial = (jnp.zeros((n, F), Y.dtype), jnp.ones(n, bool), jnp.zeros(n, bool),
               jnp.zeros(n, bool), jnp.zeros(n, jnp.int32), jnp.zeros(n))
    # The cap bounds compute only; the returned convergence flags are explicit.
    state = jax.lax.while_loop(lambda s: jnp.any(s[1]) & (jnp.max(s[4]) < 100), iterate, initial)
    u, _, converged, failed, iterations, gain = state
    return u, converged, failed, iterations, gain


def _chunk_scores(Y, L, x, R, log_rate, gamma, W, inv_phi, prior_prec):
    """Array-only compatibility entry point; programme_scores exposes diagnostics."""
    if W.shape[0] == 0:
        return jnp.zeros((Y.shape[0], 0), Y.dtype)
    return _chunk_scores_fit(Y, L, x, R, log_rate, gamma, W, inv_phi, prior_prec)[0]


def _chunks(n, per_row, max_elements):
    # Chunk rows are padded to a power of two (at most the memory cap), so specimens of different sizes reuse
    # a few compiled shapes instead of compiling one per cell count.
    cap = int(max(16, max_elements // max(1, per_row)))
    size = min(cap, max(16, 1 << (int(n) - 1).bit_length()))
    for a in range(0, n, size):
        yield a, min(a + size, n), size


def _pad(arr, size, fill=0.0):
    # Pad to one compiled shape.
    out = np.full((size,) + np.shape(arr)[1:], fill, np.float32)
    out[:len(arr)] = arr
    return jnp.asarray(out)


def cell_loglik(Y, L, log_rate, phi, x=None, gamma=None, U=None, W=None, max_elements=2 ** 24):
    """Y [n, J] counts, L [n] binned totals, log_rate [K, J], phi [J] -> [n, K] NB log-likelihoods relative to
    clone 0 (row 0 of log_rate): each cell's values are shifted by a per-cell constant (clone 0's log-likelihood
    without the gammaln terms), which posteriors, assignment margins and the Potts lambda search ignore.

    Optional x [n] (centred log depth) and gamma [J] (bin depth slopes);
    U [n, F] programme scores and W [F, J] programme loadings."""
    n, J = Y.shape
    K = log_rate.shape[0]
    log_rate = np.asarray(log_rate, np.float64)
    rate0 = jnp.asarray(log_rate[0], jnp.float32)
    delta = jnp.asarray(log_rate - log_rate[:1], jnp.float32)
    inv_phi = jnp.asarray(1.0 / np.asarray(phi), jnp.float32)
    x = np.zeros(n, np.float32) if x is None else np.asarray(x, np.float32)
    gamma = jnp.asarray(np.zeros(J) if gamma is None else gamma, jnp.float32)
    W = jnp.zeros((1, J), jnp.float32) if W is None else jnp.asarray(W, jnp.float32)
    U = np.zeros((n, W.shape[0]), np.float32) if U is None else np.asarray(U, np.float32)
    out = np.empty((n, K), np.float32)
    for a, b, size in _chunks(n, K * J, max_elements):
        out[a:b] = np.asarray(_chunk_loglik(_pad(Y[a:b], size), _pad(L[a:b], size, 1.0), _pad(x[a:b], size),
                                            _pad(U[a:b], size), rate0, delta, gamma, W, inv_phi))[:b - a]
    return out


def programme_scores(Y, L, R, log_rate, phi, W, prior_prec, x=None, gamma=None, max_elements=2 ** 24,
                     return_diagnostics=False):
    """Per-cell programme scores u [n, F], shrunk towards 0 with precision prior_prec [F].

    Optimizes the responsibility-weighted sum of clone-specific NB log-likelihoods,
    not the likelihood at an averaged log rate. Diagnostics describe numerical
    convergence of this conditional working objective, not biological validity.
    """
    n, J = Y.shape
    F = W.shape[0]
    R = np.asarray(R)
    if R.ndim != 2 or R.shape != (n, np.shape(log_rate)[0]):
        raise ValueError('Programme responsibilities must have one row per cell and one column per clone')
    if (not np.isfinite(R).all() or (R < 0).any()
            or not np.allclose(R.sum(axis=1), 1.0, rtol=0,
                               atol=8 * np.finfo(np.float32).eps * max(1, R.shape[1]))):
        raise ValueError('Programme responsibilities must be finite, nonnegative and sum to one')
    if (not np.isfinite(phi).all() or (np.asarray(phi) <= 0).any()
            or not np.isfinite(prior_prec).all() or (np.asarray(prior_prec) <= 0).any()):
        raise ValueError('Programme dispersions and prior precisions must be finite and positive')
    if (not np.isfinite(L).all() or (np.asarray(L) < 0).any() or not np.isfinite(log_rate).all()
            or not np.isfinite(W).all()):
        raise ValueError('Programme depths, clone rates and loadings must be finite, with nonnegative depths')
    diagnostics = dict(converged=np.ones(n, bool), failed=np.zeros(n, bool),
                       iterations=np.zeros(n, np.int32), objective_gain=np.zeros(n, np.float64))
    if F == 0:
        out = np.zeros((n, 0), np.float32)
        return (out, diagnostics) if return_diagnostics else out
    x = np.zeros(n, np.float32) if x is None else np.asarray(x, np.float32)
    gamma = jnp.asarray(np.zeros(J) if gamma is None else gamma, jnp.float32)
    if not np.isfinite(x).all() or not np.isfinite(gamma).all():
        raise ValueError('Programme depth offsets must be finite')
    args = (jnp.asarray(log_rate, jnp.float32), gamma, jnp.asarray(W, jnp.float32),
            jnp.asarray(1.0 / np.asarray(phi), jnp.float32), jnp.asarray(prior_prec, jnp.float32))
    out = np.empty((n, F), np.float32)
    for a, b, size in _chunks(n, max(F, R.shape[1]) * J, max_elements):
        if not np.isfinite(Y[a:b]).all() or (Y[a:b] < 0).any():
            raise ValueError('Programme counts must be finite and nonnegative')
        result = _chunk_scores_fit(_pad(Y[a:b], size), _pad(L[a:b], size, 1.0), _pad(x[a:b], size),
                                   _pad(R[a:b], size), *args)
        out[a:b] = np.asarray(result[0])[:b - a]
        for key, value in zip(diagnostics, result[1:]):
            diagnostics[key][a:b] = np.asarray(value)[:b - a]
    if not np.isfinite(out).all():
        raise FloatingPointError('Nonfinite programme scores')
    if not diagnostics['converged'].all():
        warnings.warn(f"Programme optimizer did not converge for {int((~diagnostics['converged']).sum())} cells; "
                      'returning the last finite objective-checked scores', RuntimeWarning)
    return (out, diagnostics) if return_diagnostics else out


def posteriors(loglik, log_pi):
    z = loglik + log_pi[None, :]
    z = z - z.max(axis=1, keepdims=True)
    p = np.exp(z)
    return p / p.sum(axis=1, keepdims=True)
