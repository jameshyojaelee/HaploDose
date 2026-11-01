"""Count likelihoods.

Pooled counts: Poisson-lognormal. y ~ Poisson(exp(log_mu + o)), o ~ N(0, tau^2).
Sampling noise shrinks with depth; the gene/bin-level offset o does not. The
integral over o uses adaptive Gauss-Hermite quadrature centred at the posterior
mode (Newton steps), which stays accurate at high depth where fixed nodes fail.

Single cells: negative binomial with bin-specific dispersion.
"""
from __future__ import annotations

import jax
import jax.numpy as jnp
import numpy as np
from jax.scipy.special import gammaln, logsumexp

_GH_X, _GH_W = np.polynomial.hermite.hermgauss(12)
GH_X = jnp.asarray(_GH_X, jnp.float32)
GH_LOGW = jnp.asarray(np.log(_GH_W) + _GH_X ** 2, jnp.float32)


@jax.jit
def pln_logpmf(y, log_mu, tau):
    y, log_mu, tau = jnp.broadcast_arrays(jnp.asarray(y, jnp.float32),
                                          jnp.asarray(log_mu, jnp.float32),
                                          jnp.asarray(tau, jnp.float32))
    log_mu = jnp.maximum(log_mu, -30.0)
    inv = 1.0 / tau ** 2
    # Start from the shrunken log ratio, then Newton on h(o) = y o - mu e^o - o^2/(2 tau^2).
    shrink = tau ** 2 / (tau ** 2 + 1.0 / (y + 0.5))
    o = shrink * (jnp.log(y + 0.5) - jnp.logaddexp(log_mu, jnp.log(0.5)))
    for _ in range(8):
        e = jnp.exp(log_mu + o)
        step = (y - e - o * inv) / (e + inv)
        o = o + jnp.clip(step, -1.0, 1.0)
    sd = 1.0 / jnp.sqrt(jnp.exp(log_mu + o) + inv)
    nodes = o[..., None] + jnp.sqrt(2.0) * sd[..., None] * GH_X
    h = (y[..., None] * (log_mu[..., None] + nodes) - jnp.exp(log_mu[..., None] + nodes)
         - 0.5 * nodes ** 2 * inv[..., None])
    log_int = jnp.log(jnp.sqrt(2.0) * sd) + logsumexp(h + GH_LOGW, axis=-1)
    return log_int - 0.5 * jnp.log(2 * jnp.pi * tau ** 2) - gammaln(y + 1.0)


@jax.jit
def nb_logpmf(y, mu, phi):
    """Negative binomial with mean mu and variance mu + phi mu^2."""
    r = 1.0 / phi
    return (gammaln(y + r) - gammaln(r) - gammaln(y + 1.0)
            + y * jnp.log(mu / (mu + r)) + r * jnp.log(r / (mu + r)))


@jax.jit
def gaussian_logpdf(x, mean, var):
    return -0.5 * (jnp.log(2 * jnp.pi * var) + (x - mean) ** 2 / var)
