"""Likelihoods against SciPy (NB) and numerical integration (Poisson-lognormal)."""
import numpy as np
import pytest
from scipy import integrate, stats

from haplodose.kernels.likelihood import nb_logpmf, pln_logpmf


def pln_quad(y, mu, tau):
    f = lambda o: stats.poisson.pmf(y, mu * np.exp(o)) * stats.norm.pdf(o, 0, tau)
    lo, hi = -10 * tau, 10 * tau
    pts = [np.log((y + 0.5) / mu)] if lo < np.log((y + 0.5) / mu) < hi else None
    return np.log(integrate.quad(f, lo, hi, points=pts, limit=400, epsabs=0)[0])


@pytest.mark.parametrize('y,mu,tau', [(0, 0.3, 0.2), (3, 2.0, 0.1), (40, 25.0, 0.3), (2000, 1500.0, 0.15),
                                      (120, 400.0, 0.5), (5, 5.0, 0.02)])
def test_pln_matches_quadrature(y, mu, tau):
    assert abs(float(pln_logpmf(y, np.log(mu), tau)) - pln_quad(y, mu, tau)) < 2e-3


def test_nb_matches_scipy():
    y = np.array([0, 1, 5, 40, 300])
    mu, phi = 30.0, 0.2
    r = 1 / phi
    expected = stats.nbinom.logpmf(y, r, r / (r + mu))
    np.testing.assert_allclose(np.asarray(nb_logpmf(y, mu, phi)), expected, rtol=1e-4, atol=1e-4)
