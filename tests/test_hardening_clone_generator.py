"""Independent latent-effect checks of the prespecified same-nuisance two-clone alternative."""
import importlib.util
from pathlib import Path

import numpy as np


def runner():
    path = Path(__file__).resolve().parents[1] / 'benchmarks/hardening/clones.py'
    spec = importlib.util.spec_from_file_location('hardening_clone_generator', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_alternative_shares_the_exact_same_nuisance_law_and_realization_at_equal_seed():
    module = runner()
    n, J = 100, 40
    null = module.generate(np.random.default_rng(13), 'correlated_state', n, J, return_latent=True)
    alternative = module.generate(np.random.default_rng(13), 'two_clones', n, J, return_latent=True)
    np.testing.assert_array_equal(null[4]['depth'], alternative[4]['depth'])
    np.testing.assert_array_equal(null[4]['state_effect'], alternative[4]['state_effect'])
    assert np.std(alternative[4]['state_effect']) > 0
    expected_difference = np.zeros((n, J))
    expected_difference[n // 2:, J // 2:3 * J // 4] = .6 * np.log(1.5)
    np.testing.assert_allclose(alternative[4]['log_rate_effect'] - null[4]['log_rate_effect'], expected_difference,
                               rtol=1e-14, atol=1e-14)
    np.testing.assert_allclose(alternative[4]['mean'].sum(axis=1), alternative[4]['depth'])
    assert null[3].max() == 0 and alternative[3].sum() == n // 2


def test_matched_nb_latent_state_is_absent_and_fixture_is_deterministic():
    module = runner()
    a = module.generate(np.random.default_rng(11), 'matched_nb', 60, 40, return_latent=True)
    b = module.generate(np.random.default_rng(11), 'matched_nb', 60, 40, return_latent=True)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[4]['state_effect'], np.zeros(60))
