"""Cell x clone log-likelihoods from per-cell bin counts.

Clone k predicts cell c's bin counts as L_c * rate_kj with rates summing to 1,
so only the compositional (dosage) differences between clones matter. Negative
binomial with bin dispersion phi_j; terms constant across clones are dropped.

Optional per-cell nuisance terms enter the log rate before renormalization:
depth, x_c * gamma_j, and expression programmes, sum_f u_cf * W_fj, with W
learned on normal cells (see `haplodose.programmes`).
"""
from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np


@jax.jit
def _chunk_loglik(Y, L, x, U, log_rate, gamma, W, inv_phi):
    off = x[:, None] * gamma[None] + U @ W                                # n, J
    lr = log_rate[None] + off[:, None, :]                                 # n, K, J
    lr = lr - jax.scipy.special.logsumexp(lr, axis=-1, keepdims=True)
    m = L[:, None, None] * jnp.exp(lr)
    return jnp.sum(Y[:, None, :] * lr - (Y[:, None, :] + inv_phi) * jnp.log1p(m / inv_phi), axis=-1)


@jax.jit
def _chunk_scores(Y, L, x, R, log_rate, gamma, W, inv_phi, prior_prec):
    # Fisher scoring for u_c on the normalised NB objective around the responsibility-weighted clone log rate
    # (R @ log_rate, a surrogate for the responsibility-weighted likelihood). The rate is a softmax of
    # lr0 + u W, so d log m_j / d u_f = W_fj - sum_k s_k W_fk (loadings centred on the cell's composition s).
    lr0 = R @ log_rate + x[:, None] * gamma[None]                         # n, J
    phi = 1.0 / inv_phi

    def step(u, _):
        lr = lr0 + u @ W
        s = jnp.exp(lr - jax.scipy.special.logsumexp(lr, axis=-1, keepdims=True))
        m = L[:, None] * s
        a = (Y - m) / (1.0 + phi[None] * m)                               # d loglik / d log m
        w = m / (1.0 + phi[None] * m)                                     # Fisher weight per log m
        Wbar = s @ W.T                                                    # n, F
        g = a @ W.T - Wbar * a.sum(axis=1, keepdims=True) - u * prior_prec[None]
        Ww = w @ W.T                                                      # n, F
        H = ((W[None] * w[:, None, :]) @ W.T - Wbar[:, :, None] * Ww[:, None, :] - Ww[:, :, None] * Wbar[:, None, :]
             + Wbar[:, :, None] * Wbar[:, None, :] * w.sum(axis=1)[:, None, None] + jnp.diag(prior_prec)[None])
        du = jnp.linalg.solve(H, g[..., None])[..., 0]
        # Damped: no bin's log rate moves more than 0.5 per step (Newton can overshoot in sparse cells).
        size = jnp.max(jnp.abs(du @ W), axis=1, keepdims=True)
        return u + du * jnp.minimum(1.0, 0.5 / jnp.maximum(size, 1e-12)), None

    u, _ = jax.lax.scan(step, jnp.zeros((Y.shape[0], W.shape[0]), Y.dtype), None, length=5)
    return u


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
    """Y [n, J] counts, L [n] binned totals, log_rate [K, J], phi [J] -> [n, K].

    Optional x [n] (centred log depth) and gamma [J] (bin depth slopes);
    U [n, F] programme scores and W [F, J] programme loadings."""
    n, J = Y.shape
    K = log_rate.shape[0]
    log_rate = jnp.asarray(log_rate, jnp.float32)
    inv_phi = jnp.asarray(1.0 / np.asarray(phi), jnp.float32)
    x = np.zeros(n, np.float32) if x is None else np.asarray(x, np.float32)
    gamma = jnp.asarray(np.zeros(J) if gamma is None else gamma, jnp.float32)
    W = jnp.zeros((1, J), jnp.float32) if W is None else jnp.asarray(W, jnp.float32)
    U = np.zeros((n, W.shape[0]), np.float32) if U is None else np.asarray(U, np.float32)
    out = np.empty((n, K), np.float32)
    for a, b, size in _chunks(n, K * J, max_elements):
        out[a:b] = np.asarray(_chunk_loglik(_pad(Y[a:b], size), _pad(L[a:b], size, 1.0), _pad(x[a:b], size),
                                            _pad(U[a:b], size), log_rate, gamma, W, inv_phi))[:b - a]
    return out


def programme_scores(Y, L, R, log_rate, phi, W, prior_prec, x=None, gamma=None, max_elements=2 ** 24):
    """Per-cell programme scores u [n, F], shrunk towards 0 with precision prior_prec [F].

    The expected composition is the responsibility-weighted clone log rate
    (R [n, K]), so the cell's copy-number profile is not absorbed into u."""
    n, J = Y.shape
    F = W.shape[0]
    x = np.zeros(n, np.float32) if x is None else np.asarray(x, np.float32)
    gamma = jnp.asarray(np.zeros(J) if gamma is None else gamma, jnp.float32)
    args = (jnp.asarray(log_rate, jnp.float32), gamma, jnp.asarray(W, jnp.float32),
            jnp.asarray(1.0 / np.asarray(phi), jnp.float32), jnp.asarray(prior_prec, jnp.float32))
    out = np.empty((n, F), np.float32)
    for a, b, size in _chunks(n, max(F, R.shape[1]) * J, max_elements):
        out[a:b] = np.asarray(_chunk_scores(_pad(Y[a:b], size), _pad(L[a:b], size, 1.0), _pad(x[a:b], size),
                                            _pad(R[a:b], size), *args))[:b - a]
    return out


def posteriors(loglik, log_pi):
    z = loglik + log_pi[None, :]
    z = z - z.max(axis=1, keepdims=True)
    p = np.exp(z)
    return p / p.sum(axis=1, keepdims=True)
