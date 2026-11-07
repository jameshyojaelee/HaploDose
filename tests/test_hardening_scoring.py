"""Scorer invariants: specimen identity, callable domains, calibration and simulation truth."""
import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'benchmarks'))
from hardening.profile import loss_threshold, score_loss, nested_comparisons
from hardening.statistics import (adjust_family, holm, paired_mean_comparison,
                                  specimen_false_split)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


scorer = load('hardening_scoring_fixture', 'benchmarks/tcga_bulk/subarm.py')
implant = load('hardening_implant_fixture', 'benchmarks/kinker/implant.py')
ccle = load('hardening_ccle_fixture', 'benchmarks/ccle_bulk/score.py')
ucan = load('hardening_ucan_fixture', 'benchmarks/hardening/ucan.py')


def test_raw_verified_observations_align_to_canonical_package_symbols():
    counts = pd.DataFrame([[7, 11], [13, 17]], index=['C1orf1', 'TP53'], columns=['normal', 'tumour'])
    positions = pd.DataFrame(dict(gene=['C1orf1', 'TP53'], chromosome=['1', '17'], start=[10, 20], end=[19, 29]))
    canonical, position, mapping = ucan.canonicalize_observations(counts, positions)
    np.testing.assert_array_equal(canonical.loc[['C1ORF1', 'TP53']].to_numpy(), counts.to_numpy())
    assert position.gene.tolist() == ['C1ORF1', 'TP53']
    assert mapping.source_gene.tolist() == ['C1orf1', 'TP53']
    assert mapping.canonical_gene.tolist() == ['C1ORF1', 'TP53']
    np.testing.assert_array_equal(scorer.canonical_symbols(['C1orf1', 'TP53']), ['C1ORF1', 'TP53'])


def test_canonical_gene_ambiguity_is_rejected_without_count_collapse():
    counts = pd.DataFrame([[7], [11]], index=['C1orf1', 'C1ORF1'], columns=['tumour'])
    positions = pd.DataFrame(dict(gene=['C1orf1'], chromosome=['1'], start=[10], end=[19]))
    with pytest.raises(ValueError, match='canonical'):
        ucan.canonicalize_observations(counts, positions)
    with pytest.raises(ValueError, match='canonical'):
        scorer.canonical_symbols(counts.index)


def test_ambiguous_position_loci_are_all_excluded_not_selected_or_count_collapsed():
    positions = pd.DataFrame(dict(chromosome=['chr1', 'chr2', 'chr17'], start=[10, 20, 30]),
                             index=['C1orf1', 'C1ORF1', 'TP53'])
    resolved = scorer.unambiguous_position_symbols(positions)
    assert resolved.index.tolist() == ['TP53']
    excluded = resolved.attrs['symbol_position_exclusions']
    assert len(excluded) == 2
    assert {r['source_symbol'] for r in excluded} == {'C1orf1', 'C1ORF1'}
    assert {r['canonical_symbol'] for r in excluded} == {'C1ORF1'}


def test_ucan_source_selection_resets_nonconsecutive_metadata_index(tmp_path):
    inp = tmp_path / 'input'
    inp.mkdir()
    samples = [ucan.LOCKED_SAMPLES[0], 'CRC.SW.U9999.T', ucan.LOCKED_SAMPLES[1],
               'CRC.SW.U0001.N', 'CRC.SW.U0002.N', 'CRC.SW.U0003.N']
    meta = pd.DataFrame(dict(sample=samples, role=['tumour'] * 3 + ['normal'] * 3,
                             patient=['U0001', 'U9999', 'U0008', 'U0001', 'U0002', 'U0003']))
    meta.to_csv(inp / 'samples.tsv', sep='\t', index=False)
    X = pd.DataFrame(np.arange(12).reshape(2, 6) + 1, index=['C1orf1', 'TP53'], columns=samples)
    X.index.name = 'gene'
    X.to_csv(inp / 'counts.tsv', sep='\t')
    pd.DataFrame(dict(gene=X.index, chromosome=['1', '17'], start=[10, 20], end=[19, 29])).to_csv(
        inp / 'positions.tsv', sep='\t', index=False)
    annotation = tmp_path / 'annotation.tsv'
    pd.DataFrame(dict(gene_id=['ENSG000001', 'ENSG000002'], gene_name=X.index,
                      biotype=['protein_coding'] * 2, contig=['1', '17'], start=[10, 20], end=[19, 29])).to_csv(
        annotation, sep='\t', index=False)
    source = X.copy()
    source.index = ['ENSG000001', 'ENSG000002']
    source.index.name = 'gene_id'
    raw = tmp_path / 'exon_counts.tsv.gz'
    source.to_csv(raw, sep='\t')
    (tmp_path / 'per_sample').mkdir()
    for sample in ucan.LOCKED_SAMPLES:
        source[[sample]].to_csv(tmp_path / 'per_sample' / f'{sample}.exon.tsv.gz', sep='\t')
    verified_counts, selected, _, evidence = ucan.verify_inputs(inp, raw, annotation)
    assert selected.index.tolist() == list(range(5))
    assert selected.iloc[1]['patient'] == 'U0008'
    assert verified_counts.columns.tolist() == [s for s in samples if s != 'CRC.SW.U9999.T']
    assert evidence['n_values_verified'] == 10


def test_missing_specimen_is_never_last_specimen():
    values = np.array([[11, 12], [21, 22]])
    got = scorer.align_rows(values, ['first', 'last'], ['missing', 'last', 'first'])
    assert np.isnan(got[0]).all()
    np.testing.assert_array_equal(got[1:], values[[1, 0]])


@pytest.mark.parametrize('source,target', [(['a', 'a'], ['a']), (['a', 'b'], ['a', 'a']),
                                          (['a', ''], ['a'])])
def test_ambiguous_ids_fail(source, target):
    with pytest.raises(ValueError, match='unique IDs'):
        scorer.align_rows(np.ones((2, 3)), source, target)


@pytest.mark.parametrize('bad_bin', [-1, 2, .5])
def test_invalid_gene_mapping_fails(bad_bin):
    with pytest.raises(ValueError, match='out-of-range or noninteger'):
        scorer.gene_values(np.ones((2, 2)), pd.Series([bad_bin], index=['g']), ['g'])


def test_window_truth_requires_complete_coverage_and_excludes_mixed_tie():
    rel = np.array([[-1.], [-1.], [1.], [1.]])
    fg, fl = rel > 0, rel < 0
    _, gain, loss, neutral = scorer.window_flags(rel, fg, fl, rel == 0, np.zeros(4, int), 1)
    assert not gain.any() and not loss.any() and not neutral.any()
    rel[0] = np.nan
    _, gain, loss, neutral = scorer.window_flags(rel, rel > 0, rel < 0, rel == 0, np.zeros(4, int), 1)
    assert not gain.any() and not loss.any() and not neutral.any()


def dataset(n=2):
    shape = (4, n)
    loss = np.zeros(shape, bool)
    loss[0] = True
    return dict(tum=np.array([f's{i}' for i in range(n)]), warm=np.array(['a'] * 4),
                D=np.tile(np.array([-1., 0, 0, 0])[:, None], (1, n)), FG=np.zeros(shape, bool),
                FL=loss, NEU=~loss,
                W={m: np.tile(np.array([-2., 0, 1, 2])[:, None], (1, n)) for m in scorer.METHODS})


def test_common_support_and_failure_aware_denominators_are_separate():
    d = dataset()
    d['W']['infercnvpy'][0, 0] = np.nan
    per, _, _, _ = scorer.per_patient('test', d)
    first = per.iloc[0]
    assert first.n_focal_loss == 1 and first.n_common_loss == 0
    assert np.isnan(first.recall_loss_hd_dose)
    assert first.recall_loss_failure_aware_hd_dose == 1
    assert first.recall_loss_failure_aware_infercnvpy == 0
    assert first.callable_loss_infercnvpy == 0


def test_no_focal_events_produces_typed_empty_curve():
    d = dataset(0)
    per, _, curve_windows, curve_patients = scorer.per_patient('empty', d)
    assert per.empty and 'recall_loss_hd_dose' in per
    assert not curve_windows.empty
    assert curve_patients.empty and {'method', 'fpr', 'recall_gain', 'recall_loss'} <= set(curve_patients)
    tests = scorer.tests(per, 'empty')
    assert tests.p_holm.eq(1).all() and not tests.superior.any()


def test_mean_effect_and_ci_have_same_estimand_and_paired_population():
    a, b = np.array([.2, .3, .4, 100]), np.array([.1, .1, .1, np.nan])
    result = paired_mean_comparison(a, b, n_boot=499)
    assert result['n'] == 3 and result['n_missing'] == 1
    assert result['mean_a'] == pytest.approx(.3)
    assert result['mean_diff'] == pytest.approx(.2)
    assert result['median_a'] == pytest.approx(.3)
    assert result['ci_low'] <= result['mean_diff'] <= result['ci_high']


def test_untestable_comparators_stay_in_full_holm_family():
    np.testing.assert_allclose(holm([.01, .04, np.nan]), [.03, .08, 1.])
    records = [paired_mean_comparison([1, 2, 3], [0, 0, 0], n_boot=499, family_size=2),
               paired_mean_comparison([], [], n_boot=499, family_size=2)]
    adjusted = adjust_family(records)
    assert adjusted[1]['p_holm'] == 1 and not adjusted[1]['superior']
    assert adjusted[0]['family_ci_low'] <= adjusted[0]['ci_low']


def test_ties_and_constant_bootstrap_are_not_spurious_evidence():
    result = paired_mean_comparison([1, 1, 1], [0, 0, 0], n_boot=499)
    assert not result['testable'] and result['p'] == 1
    result = paired_mean_comparison([1, 2, 2, 3], [0, 0, 0, 0],
                                    clusters=['a', 'a', 'b', 'b'], n_boot=499)
    assert result['n_resampling_units'] == 2


def test_old_bulk_claim_is_not_reissued_by_retrospective_scorer():
    empty = pd.DataFrame(columns=['sample', 'recall_loss_hd_dose', 'recall_loss_naive', 'recall_loss_infercnvpy'])
    result = ccle.primary(empty, set())['all']
    assert not result['claim_holds'] and not result['descriptive_rule_met']
    assert len(result['comparisons']) == 2


def test_external_threshold_does_not_consult_evaluation_patients():
    scores = np.column_stack([np.arange(20.), np.arange(20.)])
    neutral = np.ones_like(scores, bool)
    threshold = loss_threshold(scores, neutral, [0])
    scores[:, 1] = -1e6
    assert loss_threshold(scores, neutral, [0]) == threshold
    # Exact ties remain conservative, not jittered into false discoveries.
    assert loss_threshold(np.zeros((20, 2)), neutral, [0]) == 0


def test_threshold_weights_patients_equally_not_by_window_count():
    scores = np.ones((100, 2))
    scores[:, 0] = np.nan
    scores[0, 0] = 0
    neutral = np.isfinite(scores)
    assert loss_threshold(scores, neutral, [0, 1], fpr=.25) == 0


def test_uncalibrated_loss_is_missing_not_a_false_zero():
    values = np.array([[-2.], [0.]])
    loss, neutral = values < 0, values == 0
    row = score_loss(values, loss, neutral, np.ones_like(loss), np.nan)[0]
    assert np.isnan(row['recall_loss']) and np.isnan(row['fpr_loss'])


@pytest.mark.parametrize('missing_positive,expected_precision,expected_calls,expected_coverage',
                         [(False, 2 / 3, 3, 1.0), (True, .5, 2, .75)])
def test_loss_precision_uses_classified_callable_truth(missing_positive, expected_precision,
                                                      expected_calls, expected_coverage):
    values = np.array([-2., -1., -3., 0., -4.])[:, None]
    loss = np.array([True, True, False, False, False])[:, None]
    neutral = np.array([False, False, False, True, False])[:, None]
    gain = np.array([False, False, True, False, False])[:, None]
    common = np.ones_like(loss)
    if missing_positive:
        common[1] = False
    row = score_loss(values, loss, neutral, common, -.5, negative=gain | neutral)[0]
    assert row['precision_loss'] == pytest.approx(expected_precision)
    assert row['n_calls'] == expected_calls
    assert row['n_false_loss_calls'] == 1
    assert row['n_classified_truth'] == 4
    assert row['classified_truth_coverage'] == pytest.approx(expected_coverage)
    assert row['n_unclassified_truth'] == 1 and row['n_unclassified_call_exclusions'] == 1


def test_loss_precision_empty_call_denominator_is_undefined():
    values = np.zeros((2, 1))
    loss = np.array([[True], [False]])
    neutral = ~loss
    row = score_loss(values, loss, neutral, np.ones_like(loss), -.5)[0]
    assert row['n_calls'] == 0 and np.isnan(row['precision_loss'])


def test_conflicting_precision_truth_is_rejected():
    with pytest.raises(ValueError, match='must not overlap'):
        score_loss(np.zeros((2, 1)), np.ones((2, 1), bool), np.zeros((2, 1), bool),
                   np.ones((2, 1), bool), -.5, negative=np.ones((2, 1), bool))


def test_exact_specimen_control_has_correct_boundary():
    result = specimen_false_split(0, 59)
    assert result['upper'] == pytest.approx(1 - .05 ** (1 / 59))
    assert result['controlled']
    assert not specimen_false_split(0, 58)['controlled']
    assert not specimen_false_split(0, 0)['controlled']


def test_fixed_library_implants_preserve_exact_total_for_both_groups():
    counts = np.array([10, 20, 30, 40])
    region = np.array([True, True, False, False])
    for factor in (.5, 1., 1.5):
        got = implant.resample_event(counts, region, factor, np.random.default_rng(1))
        assert got.sum() == counts.sum() and (got >= 0).all()
    zero = implant.resample_event(np.zeros(4, int), region, 1.5, np.random.default_rng(1))
    assert not zero.any()


def test_compensated_null_has_no_simulated_clone_signal():
    counts, region = np.array([4, 8, 16]), np.array([True, False, False])
    control = implant.resample_event(counts, region, 1., np.random.default_rng(19))
    compensated = implant.resample_event(counts, region, 1.5 ** 0., np.random.default_rng(19))
    np.testing.assert_array_equal(control, compensated)
