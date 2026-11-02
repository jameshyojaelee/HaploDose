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


def _log_factorial_rest(y):
    """y log y - y - log y! (0 at y = 0), without the cancellation of the direct form at large y."""
    ys = jnp.maximum(y, 1.0)
    small = jnp.where(y > 0, y * jnp.log(ys) - y, 0.0) - gammaln(y + 1.0)
    big = -0.5 * jnp.log(2 * jnp.pi * ys) - (1 / 12 - (1 / 360 - 1 / (1260 * ys ** 2)) / ys ** 2) / ys
    return jnp.where(y < 10, small, big)


def _excess(u):
    """expm1(u) - u, by series near 0 (the direct form loses relative precision there)."""
    series = u * u * (0.5 + u * (1 / 6 + u * (1 / 24 + u * (1 / 120 + u / 720))))
    return jnp.where(jnp.abs(u) < 0.1, series, jnp.expm1(u) - u)


@jax.jit
def pln_logpmf(y, log_mu, tau):
    """log PLN(y). With lambda = e^(log_mu + o) and u = log(lambda / y), the Poisson part is written as
    y log lambda - lambda - log y! = -y (e^u - 1 - u) + (y log y - y - log y!), so no large, nearly equal
    terms are subtracted in float32 (the direct form gave log pmf +16 at y = mu = 1e7)."""
    y, log_mu, tau = jnp.broadcast_arrays(jnp.asarray(y, jnp.float32),
                                          jnp.asarray(log_mu, jnp.float32),
                                          jnp.asarray(tau, jnp.float32))
    log_mu = jnp.maximum(log_mu, -30.0)
    inv = 1.0 / tau ** 2
    pos = y > 0
    a = log_mu - jnp.log(jnp.maximum(y, 1.0))              # u = a + o for y > 0

    def resid(o):
        # y - lambda, stable at large y.
        return jnp.where(pos, -y * jnp.expm1(a + o), -jnp.exp(log_mu + o))

    def lam(o):
        return jnp.exp(log_mu + o)

    # Start from the shrunken log ratio, then Newton on h(o) = y o - lambda(o) - o^2/(2 tau^2).
    shrink = tau ** 2 / (tau ** 2 + 1.0 / (y + 0.5))
    o = shrink * (jnp.log(y + 0.5) - jnp.logaddexp(log_mu, jnp.log(0.5)))
    for _ in range(8):
        step = (resid(o) - o * inv) / (lam(o) + inv)
        o = o + jnp.clip(step, -1.0, 1.0)
    sd = 1.0 / jnp.sqrt(lam(o) + inv)
    nodes = o[..., None] + jnp.sqrt(2.0) * sd[..., None] * GH_X
    core = jnp.where(pos[..., None], -y[..., None] * _excess(a[..., None] + nodes),
                     -jnp.exp(log_mu[..., None] + nodes))
    h = core - 0.5 * nodes ** 2 * inv[..., None]
    log_int = jnp.log(jnp.sqrt(2.0) * sd) + logsumexp(h + GH_LOGW, axis=-1)
    return log_int - 0.5 * jnp.log(2 * jnp.pi * tau ** 2) + _log_factorial_rest(y)


@jax.jit
def nb_logpmf(y, mu, phi):
    """Negative binomial with mean mu and variance mu + phi mu^2."""
    r = 1.0 / phi
    return (gammaln(y + r) - gammaln(r) - gammaln(y + 1.0)
            + y * jnp.log(mu / (mu + r)) + r * jnp.log(r / (mu + r)))


@jax.jit
def gaussian_logpdf(x, mean, var):
    return -0.5 * (jnp.log(2 * jnp.pi * var) + (x - mean) ** 2 / var)
