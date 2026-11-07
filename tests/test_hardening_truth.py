"""Scientific invariants for corrected truth and scoring. Run on a compute node."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

H = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, H / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


dna = load('truth_dna_test', 'benchmarks/andor/scdna_cn.py')
low = load('truth_low_test', 'benchmarks/andor/scdna_lowdepth.py')
prep = load('truth_prep_test', 'benchmarks/spectrum/prepare_inputs.py')
clones = load('truth_clones_test', 'benchmarks/spectrum/dna_clones.py')
corr = load('truth_corr_test', 'benchmarks/spectrum/correspondence.py')
spatial = load('truth_spatial_test', 'benchmarks/ucan_spatial/dna_truth_v2.py')
tahoe = load('truth_tahoe_test', 'benchmarks/tahoe/score_profile_v2.py')
andor_prep = load('truth_andor_prep_test', 'benchmarks/andor/prepare_rna.py')
andor_rna = load('truth_andor_rna_test', 'benchmarks/hardening/andor_rna.py')


def test_andor_raw_feature_identity_is_unambiguous():
    with pytest.raises(ValueError, match='ambiguous.*A'):
        andor_prep.retained_gene_indices(['A', 'A', 'B'], ['A'])
    with pytest.raises(ValueError, match='missing.*C'):
        andor_prep.retained_gene_indices(['A', 'B'], ['C'])
    np.testing.assert_array_equal(andor_prep.retained_gene_indices(['A', 'A', 'B'], ['B']), [2])


def test_andor_duplicate_identity_requires_unique_historical_first_row(tmp_path):
    import gzip
    import io
    import tarfile
    path = tmp_path / 'raw.tar'
    payload = gzip.compress(b'ENSG1\tA\nENSG2\tA\n')
    with tarfile.open(path, 'w') as archive:
        entry = tarfile.TarInfo('GSM_L_genes.tsv.gz')
        entry.size = len(payload)
        archive.addfile(entry, io.BytesIO(payload))
    inputs = {'query': (pd.DataFrame({'library': ['L']}), pd.DataFrame({'gene': ['A']}), None)}
    params = {'libraries': {'GSM_L': {'line': 'L'}}}
    diagnostics = {}
    assert andor_rna.raw_gate(inputs, path, params, {'A': {'ENSG1'}}, diagnostics, True) == {}
    assert diagnostics['raw_duplicate_mapping'][0]['historical_first_row_certified']
    for targets in [{'ENSG2'}, {'ENSG1', 'ENSG2'}, set()]:
        with pytest.raises(andor_rna.IdentityGap):
            andor_rna.raw_gate(inputs, path, params, {'A': targets}, diagnostics, True)
        assert not diagnostics['raw_duplicate_mapping'][0]['historical_first_row_certified']


def test_andor_derivative_changes_only_independent_canonical_reference_count(tmp_path):
    import gzip
    import io
    import tarfile
    from scipy import sparse
    path = tmp_path / 'raw.tar'
    payloads = {'genes.tsv.gz': b'ENSG00000124593\tPRICKLE4\nENSG00000278224\tPRICKLE4\n',
                'barcodes.tsv.gz': b'cell\n',
                'matrix.mtx.gz': b'%%MatrixMarket matrix coordinate integer general\n2 1 2\n1 1 7\n2 1 2\n'}
    with tarfile.open(path, 'w') as archive:
        for suffix, raw in payloads.items():
            payload = gzip.compress(raw)
            entry = tarfile.TarInfo('GSM_HGC-27_1_' + suffix)
            entry.size = len(payload)
            archive.addfile(entry, io.BytesIO(payload))
    cells = pd.DataFrame(dict(library=['HGC-27_1'], cell_id=['HGC-27_1_cell'], patient=['HGC-27'],
                              role=['reference'], ref_group=['HGC-27'], ref_split=['fit'], specimen=['ref_HGC-27']))
    inputs = {'SNU-601': (cells, pd.DataFrame({'gene': ['PRICKLE4']}), sparse.csr_matrix([[7]]))}
    params = {'libraries': {'GSM_HGC-27_1': {'line': 'HGC-27'}}, 'min_umi': 0}
    diagnostics = {}
    identities = {'PRICKLE4': {'ENSG00000278224'}}
    assert andor_rna.raw_gate(inputs, path, params, identities, diagnostics, derive=True) == {'SNU-601': 1}
    assert inputs['SNU-601'][2].toarray().tolist() == [[2]]
    assert diagnostics['count_deltas'][0]['old'] == 7
    assert diagnostics['count_deltas'][0]['new'] == 2
    pd.testing.assert_frame_equal(inputs['SNU-601'][0], cells)
    assert andor_rna.raw_gate(inputs, path, params, identities, canonical=True) == {'SNU-601': 1}


def bins(chrom, start, end, arm=None):
    return pd.DataFrame(dict(chromosome=chrom, start=start, end=end,
                             arm=arm if arm is not None else chrom))


def test_dna_run_cannot_bridge_chromosomes_or_missing_bins():
    b = bins(['chr1'] * 3 + ['chr2'] * 3, [0, 1, 2, 0, 1, 2], [1, 2, 3, 1, 2, 3])
    assert dna.diff_mb(np.ones(6), np.zeros(6), b) == 0
    b = bins(['chr1'] * 6, [0, 1, 2, 10, 11, 12], [1, 2, 3, 11, 12, 13])
    assert dna.diff_mb(np.ones(6), np.zeros(6), b) == 0


def test_dna_supported_bp_and_arm_boundaries():
    b = bins(['chr1'] * 6, [0, 10, 20, 30, 40, 50], [10, 20, 30, 40, 50, 52])
    assert dna.diff_mb(np.ones(6), np.zeros(6), b) == pytest.approx(52 / 1e6)
    b['arm'] = ['p'] * 3 + ['q'] * 3
    assert dna.diff_mb(np.ones(6), np.zeros(6), b) == 0


def test_segmentation_cannot_smooth_across_filtered_gaps():
    b = bins(['chr1'] * 4, [0, 1, 20, 21], [1, 2, 21, 22])
    x = np.array([1., 1., 4., 4.])
    np.testing.assert_array_equal(dna.segment_genome(x, b, 1e9), x)
    with pytest.raises(ValueError, match='non-overlapping'):
        dna.genomic_runs(bins(['chr1', 'chr1'], [0, 1], [2, 3]))


def test_dna_truth_fraction_weights_supported_bp():
    b = bins(['chr1'] * 3, [0, 10, 20], [10, 20, 1020])
    b['a'], b['b'] = [1., 1., 0.], 0.
    assert dna.clone_differences(b, ['a', 'b'], b.arm.to_numpy()).empty
    assert low.arm_diffs({'a': b.a.to_numpy(), 'b': b.b.to_numpy()}, b.arm.to_numpy(),
                         'a', 'b', b.end - b.start) == []


def test_lowdepth_coarsening_retains_gaps():
    b = bins(['chr1'] * 4, [0, 1, 20, 21], [1, 2, 21, 22])
    group = low.coarse_groups(b.arm.to_numpy(), 10, bins=b)
    assert group[0] == group[1] and group[2] == group[3] and group[1] != group[2]


def test_gc_strata_with_ties_remain_finite():
    result = dna.gc_correct(np.array([[1., 2., 3.]]), np.array([.5, .5, .5]))
    assert np.isfinite(result).all()


def test_missing_gene_rejected_before_negative_index_conversion():
    with pytest.raises(ValueError, match='gene sets differ'):
        prep.align_gene_indices(['A', 'B'], ['A', 'MISSING'])
    np.testing.assert_array_equal(prep.align_gene_indices(['A', 'B', 'A'], ['B', 'A']), [1, 0])
    with pytest.raises(ValueError, match='unique'):
        prep.align_gene_indices(['A'], ['A', 'A'])


def test_missing_categorical_code_cannot_become_last_chromosome():
    with pytest.raises(ValueError, match='categorical'):
        clones.categorical({'x/categories': np.array(['chr1', 'chr2']), 'x/codes': np.array([-1])}, 'x')


def test_site_label_agreement_does_not_establish_same_biopsy(tmp_path):
    pd.DataFrame(dict(cell_id=['a', 'b'], leaf_id=['0', '1'], tumor_site=['Omentum', 'Omentum'])).to_csv(
        tmp_path / 'P_cell_info.csv.gz', index=False)
    sample = pd.DataFrame(dict(patient=['P'], sample_title=['P_S1_UNSORTED_OMENTUM']))
    result = corr.correspondence(sample, tmp_path)
    assert result.match_tier.iloc[0] == 'same_anatomical_label_only'
    assert not result.same_biopsy_verified.any()
    assert result.unseen_clone_prevalence_upper95.iloc[0] == pytest.approx(1 - .05 ** .5)


def test_spatial_overlap_uses_callable_union_and_unknown_direction():
    calls = pd.DataFrame(dict(start=[0, 0, 20], end=[10, 10, 30], mean_diff=[1., 1., -1.]))
    support = [(0, 10), (20, 30)]
    metric = spatial.event_overlap(0, 30, calls, support)
    assert metric['callable_event_bp'] == 20 and metric['overlap_bp'] == 20
    assert metric['overlap_fraction'] == 1 and not metric['direction_evaluable']
    assert np.isnan(metric['signed_overlap_fraction'])
    signed = spatial.event_overlap(0, 30, calls, support, expected_direction=1)
    assert signed['signed_overlap_fraction'] == .5
    assert spatial.event_overlap(10, 20, calls, [(0, 30)])['overlap_bp'] == 0


def test_callable_shift_null_does_not_include_unassayed_gaps():
    calls = pd.DataFrame(dict(start=[0, 100], end=[10, 110]))
    null = spatial.callable_shift_null(0, 110, calls, [(0, 10), (100, 110)],
                                       np.random.default_rng(3), draws=20)
    np.testing.assert_allclose(null, 1.)
    empty = calls.iloc[:0]
    null = spatial.callable_shift_null(0, 110, empty, [(0, 10), (100, 110)],
                                       np.random.default_rng(3), draws=20)
    np.testing.assert_allclose(null, 0.)


def test_missing_tahoe_comparator_cannot_shrink_family_or_pass():
    res = pd.DataFrame({'mbP_hd_sm5': np.ones(25), 'mbP_copykat': np.zeros(25), 'mbP_pb_sm5': np.zeros(25)})
    result = tahoe.primary_tests(res)
    assert len(result['comparisons']) == 3 and not result['claim_holds']
    missing = next(x for x in result['comparisons'] if x['comparator'] == 'infercnvpy')
    assert missing['n_lines'] == 0 and missing['p_holm'] == 1
    with pytest.raises(ValueError, match='family'):
        tahoe.primary_tests(res, ['copykat', 'pb_sm5'])


def test_profile_truth_and_methods_have_identical_gene_support():
    g = pd.DataFrame({'mb': ['chr1:0'] * 2, 'truth': [0., 100.],
                      'hd': [2., 4.], 'other': [3., np.nan]})
    result, present = tahoe.common_gene_bins(g, ['hd', 'other'])
    assert present == ['hd', 'other']
    assert result.truth.iloc[0] == 0 and result.hd.iloc[0] == 2 and result.other.iloc[0] == 3
