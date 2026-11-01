"""HMM kernels against brute-force enumeration, with chromosome restarts."""
import itertools

import jax.numpy as jnp
import numpy as np
import pytest
from scipy.special import logsumexp

from haplodose.kernels.hmm import forward_backward, sticky_transition, viterbi


def brute(log_emit, log_trans, log_init, reset):
    T, S = log_emit.shape
    paths = list(itertools.product(range(S), repeat=T))
    scores = []
    for p in paths:
        s = 0.0
        for t in range(T):
            s += log_emit[t, p[t]] + (log_init[p[t]] if reset[t] else log_trans[p[t - 1], p[t]])
        scores.append(s)
    scores = np.array(scores)
    post = np.zeros((T, S))
    w = np.exp(scores - logsumexp(scores))
    for p, wi in zip(paths, w):
        post[np.arange(T), p] += wi
    return post, logsumexp(scores), np.array(paths[int(np.argmax(scores))])


@pytest.mark.parametrize('seed', [0, 1, 2])
def test_forward_backward_and_viterbi_match_enumeration(seed):
    rng = np.random.default_rng(seed)
    T, S = 7, 3
    log_emit = rng.normal(0, 2, (T, S)).astype(np.float32)
    log_trans = np.log(rng.dirichlet(np.ones(S), S)).astype(np.float32)
    log_init = np.log(rng.dirichlet(np.ones(S))).astype(np.float32)
    reset = np.array([1, 0, 0, 1, 0, 0, 0], bool)
    post, ll, best = brute(log_emit, log_trans, log_init, reset)
    for method in ['scan', 'assoc']:
        p, l = forward_backward(jnp.asarray(log_emit), jnp.asarray(log_trans), jnp.asarray(log_init),
                                jnp.asarray(reset), method=method)
        np.testing.assert_allclose(np.asarray(p), post, atol=2e-5)
        assert abs(float(l) - ll) < 1e-3
    path = viterbi(jnp.asarray(log_emit), jnp.asarray(log_trans), jnp.asarray(log_init), jnp.asarray(reset))
    np.testing.assert_array_equal(np.asarray(path), best)


def test_scan_and_assoc_agree_on_long_chain():
    rng = np.random.default_rng(3)
    T, S = 3000, 4
    log_emit = rng.normal(0, 1, (T, S)).astype(np.float32)
    trans = sticky_transition(S, 1e-3, [0.2, 0.5, 0.2, 0.1])
    init = jnp.log(jnp.asarray([0.2, 0.5, 0.2, 0.1]))
    reset = np.zeros(T, bool)
    reset[::250] = True
    a, la = forward_backward(jnp.asarray(log_emit), trans, init, jnp.asarray(reset), method='scan')
    b, lb = forward_backward(jnp.asarray(log_emit), trans, init, jnp.asarray(reset), method='assoc')
    np.testing.assert_allclose(np.asarray(a), np.asarray(b), atol=2e-4)
    assert abs(float(la) - float(lb)) / abs(float(la)) < 1e-4


def test_sticky_transition_rows_sum_to_one():
    T = np.exp(np.asarray(sticky_transition(4, 0.01, [0.2, 0.5, 0.2, 0.1])))
    np.testing.assert_allclose(T.sum(axis=1), 1, atol=1e-6)
    np.testing.assert_allclose(np.diag(T), 0.99, atol=1e-6)
