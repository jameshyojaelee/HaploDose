"""Log-space HMM kernels with chain restarts (one chain per chromosome).

Each position t is an element M_t[i, j] = A_t[i, j] + e_t[j], where A_t is the
transition matrix, or the initial distribution broadcast over rows at a restart.
Forward and backward variables are prefix and suffix products of the M_t in the
(logsumexp, +) semiring. `lax.scan` computes them sequentially (best on CPU);
`associative_scan` computes them in O(log T) depth (best on GPU). Both carry the
mass of earlier chains as a per-position constant, which cancels in posteriors.
"""
from __future__ import annotations

from functools import partial

import jax
import jax.numpy as jnp
from jax import lax
from jax.scipy.special import logsumexp


def _elements(log_emit, log_trans, log_init, reset):
    # log_trans is [S, S] (shared) or [T, S, S] (row t: transition into position t).
    reset = reset.at[0].set(True)
    trans = log_trans if log_trans.ndim == 3 else log_trans[None]
    step = jnp.where(reset[:, None, None], log_init[None, None, :], trans)
    return step + log_emit[:, None, :]


def _lmm(a, b):
    """Batched log-semiring matrix product a (x) b."""
    return logsumexp(a[..., :, :, None] + b[..., None, :, :], axis=-2)


def _forward_scan(M):
    def step(alpha, m):
        alpha = logsumexp(alpha[:, None] + m, axis=0)
        return alpha, alpha
    alpha0 = M[0, 0]
    _, rest = lax.scan(step, alpha0, M[1:])
    return jnp.concatenate([alpha0[None], rest])


def _backward_scan(M):
    def step(beta, m):
        beta = logsumexp(m + beta[None, :], axis=1)
        return beta, beta
    last = jnp.zeros(M.shape[-1], M.dtype)
    _, rest = lax.scan(step, last, M[1:], reverse=True)
    return jnp.concatenate([rest, last[None]])


def _forward_assoc(M):
    return lax.associative_scan(_lmm, M)[:, 0, :]


def _backward_assoc(M):
    # Suffix products Q_t = M_t (x) ... (x) M_{T-1}: prefix scan of the reversed
    # sequence with the operands swapped. beta_t = rowwise logsumexp of Q_{t+1}.
    Q = lax.associative_scan(lambda a, b: _lmm(b, a), M[::-1])[::-1]
    last = jnp.zeros((1, M.shape[-1]), M.dtype)
    return jnp.concatenate([logsumexp(Q[1:], axis=-1), last])


@partial(jax.jit, static_argnames=('method',))
def forward_backward(log_emit, log_trans, log_init, reset, method='scan'):
    """Posterior state probabilities [T, S] and total log-likelihood."""
    M = _elements(log_emit, log_trans, log_init, reset)
    if method == 'scan':
        alpha, beta = _forward_scan(M), _backward_scan(M)
    elif method == 'assoc':
        alpha, beta = _forward_assoc(M), _backward_assoc(M)
    else:
        raise ValueError(method)
    gamma = alpha + beta
    gamma = gamma - logsumexp(gamma, axis=1, keepdims=True)
    return jnp.exp(gamma), logsumexp(alpha[-1])


@jax.jit
def viterbi(log_emit, log_trans, log_init, reset):
    """Most probable state path [T] (int32); restarts give independent chains."""
    M = _elements(log_emit, log_trans, log_init, reset)

    def step(delta, m):
        cand = delta[:, None] + m
        return jnp.max(cand, axis=0), jnp.argmax(cand, axis=0).astype(jnp.int32)
    last, back = lax.scan(step, M[0, 0], M[1:])
    s_last = jnp.argmax(last).astype(jnp.int32)

    def trace(s, b):
        prev = b[s]
        return prev, prev
    _, path = lax.scan(trace, s_last, back, reverse=True)
    return jnp.concatenate([path, s_last[None]])


def default_method():
    return 'assoc' if jax.default_backend() == 'gpu' else 'scan'


def forward_backward_batch(log_emit, log_trans, log_init, reset, method=None):
    """vmap over a leading batch axis of log_emit; transitions and resets are shared."""
    method = method or default_method()
    f = jax.vmap(lambda e: forward_backward(e, log_trans, log_init, reset, method))
    return f(log_emit)


viterbi_batch = jax.jit(jax.vmap(viterbi, in_axes=(0, None, None, None)))


def phase_transitions(p_switch, p_change):
    """[T, 3, 3] log transitions over (balanced, imbalanced +, imbalanced -).

    Balanced and imbalanced swap with probability p_change (entering imbalance in
    either orientation equally). Within imbalance, p_switch [T] is the probability
    that the input phase flips between SNP t-1 and t, which swaps the orientation."""
    ps = jnp.clip(jnp.asarray(p_switch, jnp.float32), 1e-12, 0.5)
    pc = jnp.float32(p_change)
    one = jnp.ones_like(ps)
    rows = [jnp.stack([(1 - pc) * one, pc / 2 * one, pc / 2 * one], -1),
            jnp.stack([pc * one, (1 - pc) * (1 - ps), (1 - pc) * ps], -1),
            jnp.stack([pc * one, (1 - pc) * ps, (1 - pc) * (1 - ps)], -1)]
    return jnp.log(jnp.stack(rows, -2))


def sticky_transition(n_states, p_change, prior):
    """Stay with probability 1 - p_change; otherwise move in proportion to prior."""
    prior = jnp.asarray(prior, jnp.float32)
    off = p_change * prior[None, :] / (1 - prior[:, None])
    T = jnp.where(jnp.eye(n_states, dtype=bool), 1 - p_change, off)
    return jnp.log(T)
