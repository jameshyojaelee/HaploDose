"""Independent mathematical checks for the schema-2 model corrections."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import optimize, special, stats

from haplodose.bulk import reference_model, run_bulk
from haplodose.kernels.estep import programme_scores
from haplodose.profile import ProfileParams, segments_table
from haplodose.programmes import Programmes, high_pass, nonzero_row_basis


def test_normal_programmes_respect_centered_rank_and_zero_variation():
    reset = np.zeros(40, bool)
    reset[[0, 20]] = True
    constant = np.tile(np.arange(1, 41, dtype=float), (4, 1))
    _, _, W = reference_model(constant, 10, reset)
    assert W.shape == (0, 40)
    rng = np.random.default_rng(81)
    normals = 1 + rng.poisson(constant * rng.lognormal(0, 0.3, constant.shape))
    _, _, W = reference_model(normals, 10, reset)
    assert len(W) <= len(normals) - 1
    _, _, reverse = reference_model(normals[::-1], 10, reset)
    np.testing.assert_allclose(W.T @ W, reverse.T @ reverse, atol=200 * np.finfo(float).eps)


def test_high_pass_null_directions_are_removed():
    reset = np.r_[True, np.zeros(19, bool)]
    arm_constant = np.ones((2, 20))
    filtered = high_pass(arm_constant, reset, 0)
    assert nonzero_row_basis(filtered, reference_scale=np.linalg.norm(arm_constant)).shape == (0, 20)
    wave = np.tile(np.arange(20, dtype=float), (2, 1))
    assert len(nonzero_row_basis(high_pass(wave, reset, 0), reference_scale=np.linalg.norm(wave))) == 1


def _weighted_nb_objective(u, y, L, R, log_rate, phi, W, prior):
    means = L * special.softmax(log_rate + u @ W, axis=1)
    ll = stats.nbinom.logpmf(y[None], 1.0 / phi, 1.0 / (1.0 + phi * means)).sum(axis=1)
    return float(R @ ll - 0.5 * np.sum(prior * u ** 2))


def test_programme_objective_uses_uncertain_clone_responsibilities():
    rng = np.random.default_rng(82)
    J, F = 30, 2
    W = rng.normal(0, 0.3, (F, J))
    log_rate = rng.normal(-np.log(J), 0.5, (3, J))
    phi = rng.uniform(0.05, 0.6, J)
    R = np.array([0.15, 0.55, 0.30])
    prior = np.array([0.8, 1.1])
    y = rng.negative_binomial(1 / phi, 1 / (1 + phi * 1000 * special.softmax(log_rate[1] + np.array([0.4, -0.7]) @ W)))
    objective = lambda u: _weighted_nb_objective(u, y, y.sum(), R, log_rate, phi, W, prior)
    independent = optimize.minimize(lambda u: -objective(u), np.zeros(F), method='BFGS')
    fitted, diagnostic = programme_scores(y[None], np.array([y.sum()]), R[None], log_rate, phi, W, prior,
                                          return_diagnostics=True)
    # The numerical acceptance budget follows float32 reduction precision,
    # not a scientific performance tolerance.
    error_budget = (J + 16) * np.finfo(np.float32).eps * (1 + abs(objective(independent.x)))
    assert objective(fitted[0]) >= objective(np.zeros(F)) - error_budget
    assert objective(independent.x) - objective(fitted[0]) <= error_budget
    assert diagnostic['converged'][0]
    assert not diagnostic['failed'][0]


def test_programme_zero_rank_and_invalid_responsibilities():
    y = np.ones((2, 5))
    args = (y, y.sum(axis=1), np.ones((2, 1)), np.zeros((1, 5)), np.ones(5), np.zeros((0, 5)), np.zeros(0))
    fitted, diagnostic = programme_scores(*args, return_diagnostics=True)
    assert fitted.shape == (2, 0)
    assert diagnostic['converged'].all()
    with pytest.raises(ValueError, match='responsibilities'):
        programme_scores(y, y.sum(axis=1), np.zeros((2, 1)), *args[3:])


def test_compositional_overlap_is_invariant_to_global_dose_shift():
    rng = np.random.default_rng(83)
    J = 15
    base = rng.dirichlet(np.ones(J))
    W = rng.normal(size=(3, J))
    programme = Programmes(W, np.ones(3), np.full(J, 0.2), 2000, 0.5)
    dose = W[0] + 0.2 * rng.normal(size=J)
    np.testing.assert_allclose(programme.overlap(dose, base), programme.overlap(dose + 0.7, base),
                               atol=J * np.finfo(float).eps)
    np.testing.assert_allclose(programme.overlap(W[0], base), 1.0, atol=J * np.finfo(float).eps)
    assert 0 <= programme.overlap(dose, base) <= 1 + J * np.finfo(float).eps


def test_segment_uncertainty_uses_the_segments_bin_variances_and_rna_units():
    Y = np.array([1e6, 1e6, 1e6, 10.0, 20.0, 30.0])
    path = np.array([2, 2, 2, 0, 0, 0])
    post = np.eye(4)[path]
    dose = np.array([0.2, 0.2, 0.2, -0.3, -0.3, -0.3])
    variance = np.array([0.0, 1.0, 2.0, 0.0, 0.0, 0.0])
    rows = segments_table(post, path, dose, Y, Y, np.r_[True, np.zeros(5, bool)], 0.5, 0.02,
                          ProfileParams(), reference_var=variance)
    expected = np.sqrt(np.sum(0.02 ** 2 + variance[:3] + 1.0 / (Y[:3] + 0.5))) / 3
    assert rows[0]['se_rna_log_ratio_approx'] == pytest.approx(expected)
    assert rows[0]['rna_log_effect'] == pytest.approx(0.2)
    assert rows[0]['observed_mean_rna_log_ratio'] == pytest.approx(0.)
    assert rows[0]['uncertainty_estimand'] is None
    assert rows[0]['uncertainty_statistic'] == 'unweighted_mean_pseudocount_log_count_ratio'
    assert rows[0]['uncertainty_center_field'] == 'observed_mean_rna_log_ratio'
    assert rows[0]['uncertainty_status'] == 'approximate_scale_diagnostic_without_validated_parameter_interval'
    assert not rows[0]['posterior_effect_interval_validated']
    assert rows[0]['relative_state_log_ratio'] == pytest.approx(0.4)
    assert rows[0]['three_sigma_rna_log_ratio_approx'] == pytest.approx(3 * expected)
    assert not rows[0]['uncertainty_calibrated']
    assert rows[0]['uncertainty_scope'] == 'conditional_on_fit_and_selected_interval'
    assert 'min_detectable_log_dose' not in rows[0]


def test_bulk_requires_declared_raw_integer_counts(tmp_path):
    with pytest.raises(ValueError, match='explicitly be raw'):
        run_bulk('absent', 'absent', 'absent', tmp_path / 'absent')
    genes = ['G0', 'G1']
    pd.DataFrame({'N0': [1.25, 2], 'N1': [3, 4], 'N2': [4, 5]}, index=genes).to_csv(
        tmp_path / 'counts.tsv', sep='\t')
    pd.DataFrame({'sample': ['N0', 'N1', 'N2'], 'role': ['normal'] * 3}).to_csv(
        tmp_path / 'samples.tsv', sep='\t', index=False)
    pd.DataFrame({'gene': genes, 'chromosome': ['chr1'] * 2, 'start': [0, 1000], 'end': [500, 1500]}).to_csv(
        tmp_path / 'positions.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='integer-valued'):
        run_bulk(tmp_path / 'counts.tsv', tmp_path / 'samples.tsv', tmp_path / 'positions.tsv', tmp_path / 'out',
                 counts_kind='raw', genome_build='hg38', coordinate_system='zero-based-half-open',
                 gene_namespace='custom')
    record = json.loads((tmp_path / 'out' / 'run_manifest.json').read_text())
    assert record['status'] == 'failed'
    assert record['error']['type'] == 'ValueError'
    assert record['execution']['parameters']['variance_mode'] == 'conservative_roughness_plus_reference'
    assert record['input']['files']['counts']['sha256']


def test_bulk_never_reuses_existing_outputs(tmp_path):
    (tmp_path / 'keep.txt').write_text('preserved')
    with pytest.raises(FileExistsError):
        run_bulk('absent', 'absent', 'absent', tmp_path, counts_kind='raw', genome_build='hg38',
                 coordinate_system='zero-based-half-open', gene_namespace='custom')
    assert (tmp_path / 'keep.txt').read_text() == 'preserved'
