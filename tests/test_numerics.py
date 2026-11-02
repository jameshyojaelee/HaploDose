"""Numerical invariants from the 2026-11-01 efficiency review (docs/efficiency_review_2026-11-01.md)."""
import jax.numpy as jnp
import numpy as np
import pytest
from scipy import integrate, optimize, stats
from scipy.special import gammaln, logsumexp

from haplodose.kernels.estep import _chunk_scores
from haplodose.kernels.hmm import forward_backward, sticky_transition
from haplodose.kernels.likelihood import pln_logpmf


def pln_quad_scaled(y, mu, tau):
    """log PLN(y) by quadrature around the mode with the Poisson log pmf in float64 (no underflow at high y)."""
    lp = lambda o: stats.poisson.logpmf(y, mu * np.exp(o)) + stats.norm.logpdf(o, 0, tau)
    mode = optimize.minimize_scalar(lambda o: -lp(o), bounds=(-10 * tau, 10 * tau), method='bounded').x
    top = lp(mode)
    val = integrate.quad(lambda o: np.exp(lp(o) - top), mode - 12 * tau, mode + 12 * tau, points=[mode],
                         limit=400, epsabs=0)[0]
    return top + np.log(val)


@pytest.mark.parametrize('y', [1e5, 1e6, 1e7])
def test_pln_high_counts(y):
    """Review finding 1: float32 cancellation gave +16 at y = mu = 1e7 (a log pmf must be negative)."""
    got = float(pln_logpmf(y, np.log(y), 0.15))
    assert got < 0
    assert abs(got - pln_quad_scaled(y, y, 0.15)) < 0.02


@pytest.mark.parametrize('y,mu', [(0, 0.3), (1, 0.05), (7, 3.0), (300, 450.0), (5e4, 4e4)])
def test_pln_low_and_mid_counts(y, mu):
    assert abs(float(pln_logpmf(y, np.log(mu), 0.2)) - pln_quad_scaled(y, mu, 0.2)) < 5e-3


@pytest.mark.parametrize('method', ['scan', 'assoc'])
def test_hmm_posterior_invariant_to_emission_offset(method):
    """Review finding 2: adding the same constant to every state's emission must not change posteriors."""
    rng = np.random.default_rng(0)
    T, S = 10000, 4
    e = rng.normal(0, 1, (T, S)).astype(np.float32)
    trans = sticky_transition(S, 2e-3, [0.2, 0.5, 0.2, 0.1])
    init = jnp.log(jnp.asarray([0.2, 0.5, 0.2, 0.1]))
    reset = np.zeros(T, bool)
    reset[::250] = True
    a, la = forward_backward(jnp.asarray(e), trans, init, jnp.asarray(reset), method=method)
    b, lb = forward_backward(jnp.asarray(e - 100.0), trans, init, jnp.asarray(reset), method=method)
    a, b = np.asarray(a), np.asarray(b)
    assert np.abs(a - b).max() < 1e-4
    np.testing.assert_allclose(b.sum(axis=1), 1, atol=1e-5)
    assert abs((float(la) - float(lb)) - 100.0 * T) / (100.0 * T) < 1e-5


def _nb_objective(u, y, L, lr0, W, phi, prec):
    lr = lr0 + u @ W
    lr = lr - logsumexp(lr)
    m = L * np.exp(lr)
    r = 1 / phi
    ll = np.sum(gammaln(y + r) - gammaln(r) + y * np.log(m / (m + r)) + r * np.log(r / (m + r)))
    return ll - 0.5 * np.sum(prec * u ** 2)


def test_programme_scores_maximise_normalised_nb():
    """Review finding 3: the score update must use the gradient of the softmax-normalised NB objective."""
    rng = np.random.default_rng(1)
    J, F = 60, 2
    W = rng.normal(0, 0.4, (F, J))
    lr0 = np.log(rng.dirichlet(np.ones(J) * 2))
    phi = rng.uniform(0.05, 2.0, J)                    # unequal dispersions: the softmax term no longer cancels
    u_true = np.array([0.8, -0.6])
    p = np.exp(lr0 + u_true @ W)
    p /= p.sum()
    y = rng.negative_binomial(1 / phi, 1 / (1 + phi * 4000 * p)).astype(np.float64)
    L = y.sum()
    prec = np.full(F, 1e-3)
    u = np.asarray(_chunk_scores(jnp.asarray(y[None], jnp.float32), jnp.asarray([L], jnp.float32),
                                 jnp.zeros(1, jnp.float32), jnp.ones((1, 1), jnp.float32),
                                 jnp.asarray(lr0[None], jnp.float32), jnp.zeros(J, jnp.float32),
                                 jnp.asarray(W, jnp.float32), jnp.asarray(1 / phi, jnp.float32),
                                 jnp.asarray(prec, jnp.float32)))[0]
    best = optimize.minimize(lambda v: -_nb_objective(v, y, L, lr0, W, phi, prec), np.zeros(F)).x
    np.testing.assert_allclose(u, best, atol=0.01)


def test_nb_rate_score_equation():
    """Review finding 4: the weighted NB rate update must solve the NB score equation, not the Poisson one."""
    from haplodose.clones import nb_rate
    rng = np.random.default_rng(2)
    E = rng.gamma(1.0, 50.0, 400)
    r = 3.0
    O = rng.negative_binomial(r, r / (r + 0.5 * E)).astype(np.float64)
    w = rng.uniform(0.2, 1.0, 400)
    nll = lambda t: -np.sum(w * (gammaln(O + r) - gammaln(r) + O * np.log(t * E / (t * E + r))
                                 + r * np.log(r / (t * E + r))))
    best = optimize.minimize_scalar(nll, bounds=(1e-4, 10), method='bounded', options=dict(xatol=1e-8)).x
    assert abs(nb_rate(O, E, r, w) - best) < 1e-4


def test_pooled_expectation_includes_depth_variance():
    """Review finding 5: the pooled expectation averages softmax(log base + gamma x_c) over cells, which the
    pooled mean depth alone misses for bins with steep depth slopes."""
    from types import SimpleNamespace
    from haplodose.clones import SpecimenModel, depth_moments
    rng = np.random.default_rng(3)
    J, n = 200, 3000
    base = rng.dirichlet(np.ones(J) * 3)
    gamma = rng.uniform(-1.5, 1.5, J)
    gamma -= np.average(gamma, weights=base)
    L = rng.lognormal(7, 0.6, n)
    x = np.log(L) - np.log(L).mean()
    A = np.log(base)[None] + gamma[None] * x[:, None]
    exact = (L[:, None] * np.exp(A - logsumexp(A, axis=1, keepdims=True))).sum(axis=0)
    m = SimpleNamespace(base=base, gamma=gamma)
    m.depth_log_rate = lambda Lk, XLk: SpecimenModel.depth_log_rate(m, Lk, XLk)
    got = SpecimenModel.expected(m, L.sum(), depth_moments(L, x).sum(axis=0))[0]
    mean_only = SpecimenModel.expected(m, L.sum(), np.array([(L * x).sum(), (L * x).sum() ** 2 / L.sum()]))[0]
    assert np.abs(got / exact - 1).max() < 0.03
    assert np.abs(mean_only / exact - 1).max() > 0.1       # the test can fail for the reason it tests


def test_nb_rate_per_bin_sizes():
    """Per-bin NB sizes (cells x bins, per-cell weights [n, 1]) must solve the full NB score equation."""
    from haplodose.clones import nb_rate
    rng = np.random.default_rng(4)
    n, J = 300, 6
    E = rng.gamma(2.0, 5.0, (n, J))
    r = np.array([0.5, 1.0, 2.0, 5.0, 20.0, 100.0])[None]
    O = rng.negative_binomial(r, r / (r + 1.3 * E)).astype(np.float64)
    w = rng.uniform(0.2, 1.0, (n, 1))
    nll = lambda t: -np.sum(w * (O * np.log(t * E / (t * E + r)) + r * np.log(r / (t * E + r))))
    best = optimize.minimize_scalar(nll, bounds=(1e-3, 10), method='bounded', options=dict(xatol=1e-9)).x
    assert abs(nb_rate(O, E, r, w) - best) < 1e-4


def test_event_mixtures_ignores_one_overdispersed_gene():
    """A window event in 30% of cells must outrank a chromosome holding one very overdispersed, highly
    expressed bin (FTH1-like); region totals with a median dispersion ranked the latter first."""
    from types import SimpleNamespace
    from haplodose.clones import CloneParams, SpecimenModel
    rng = np.random.default_rng(5)
    n, per_chrom, n_chrom = 400, 120, 3
    J = per_chrom * n_chrom
    reset = np.zeros(J, bool)
    reset[::per_chrom] = True
    base = rng.dirichlet(np.ones(J) * 5)
    phi = np.full(J, 0.1)
    base[10] *= 40                                        # chrom 1, one highly expressed bin ...
    phi[10] = 3.0                                         # ... with large cell-to-cell variance
    base /= base.sum()
    L = rng.lognormal(np.log(4000), 0.3, n)
    shift = np.zeros((n, J))
    sub = rng.random(n) < 0.3
    ev = np.arange(per_chrom + 40, per_chrom + 64)        # chrom 2: 24 bins (12 proposal bins) at 0.6x
    shift[np.ix_(sub, ev)] = np.log(0.6)
    mu = L[:, None] * base[None] * np.exp(shift)
    Y = rng.negative_binomial(1 / phi, 1 / (1 + phi * mu)).astype(np.float64)
    m = SimpleNamespace(params=CloneParams(event_windows=(8, 24)), propose_mask=np.arange(J) % 2 == 0,
                        gamma=np.zeros(J), phi=phi,
                        graph=None, reset=reset, beta=1.0)
    out = SpecimenModel.event_mixtures(m, Y, L, np.zeros(n), np.log(base), np.zeros((n, J)))
    g, src = out[0]
    assert src.startswith('event_chrom2_')
    agree = max(np.mean(g.astype(bool) == sub), np.mean(g.astype(bool) != sub))
    assert agree > 0.9


def _nb_loglik64(Y, L, log_rate, phi, off):
    """float64 reference: absolute NB log-likelihood without the gammaln terms, [n, K]."""
    lr = log_rate[None] + off[:, None, :]
    lr = lr - np.log(np.exp(lr).sum(axis=-1, keepdims=True))
    m = L[:, None, None] * np.exp(lr)
    r = 1.0 / phi
    return (Y[:, None, :] * lr - (Y[:, None, :] + r) * np.log1p(m / r)).sum(axis=-1)


def test_cell_loglik_keeps_small_clone_contrasts():
    """Review 2026-11-02 finding 9: a +0.01 log-rate difference in 20 of 2,000 bins, counts at clone 0's mean.
    The float32 absolute sums returned a difference of exactly 0; the relative form must match float64."""
    from haplodose.kernels.estep import cell_loglik
    J, phi = 2000, np.full(2000, 0.1)
    base = np.full(J, np.log(1.0 / J))
    log_rate = np.stack([base, base + np.r_[np.full(20, 0.01), np.zeros(J - 20)]])
    for L in (2e4, 2e5):
        Y = (L * np.exp(base))[None]
        ref = _nb_loglik64(Y, np.array([L]), log_rate, phi, np.zeros((1, J)))
        got = cell_loglik(Y, np.array([L]), log_rate, phi)
        want = ref[0, 1] - ref[0, 0]
        assert want < -0.004 and abs((got[0, 1] - got[0, 0]) - want) < 0.02 * abs(want)


def test_cell_loglik_matches_float64_with_depth_and_programmes():
    from haplodose.kernels.estep import cell_loglik
    rng = np.random.default_rng(7)
    n, K, J, F = 50, 4, 300, 2
    log_rate = np.log(rng.dirichlet(np.ones(J), K)) + rng.normal(0, 0.3, (K, J))
    phi = rng.uniform(0.05, 0.5, J)
    x, gamma = rng.normal(0, 0.5, n), rng.normal(0, 0.2, J)
    U, W = rng.normal(0, 0.5, (n, F)), rng.normal(0, 0.3, (F, J))
    L = rng.uniform(1e3, 5e4, n)
    Y = rng.poisson(L[:, None] / J * 2.0, (n, J)).astype(float)
    ref = _nb_loglik64(Y, L, log_rate, phi, x[:, None] * gamma[None] + U @ W)
    got = cell_loglik(Y, L, log_rate, phi, x=x, gamma=gamma, U=U, W=W)
    assert np.allclose(got - got[:, :1], ref - ref[:, :1], rtol=1e-4, atol=1e-2)
