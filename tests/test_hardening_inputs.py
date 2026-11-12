"""Input identities, coordinate equivalence, sparse structure and genomic boundaries."""
import h5py
import numpy as np
import pandas as pd
import pytest

from haplodose.bins import arm_labels, chromosome_starts, make_bins
from haplodose.io import CountSource, read_cells, read_positions, validate_feature_ids, validate_input_contract


CONTRACT = dict(genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')


def counts_file(path, *, genes=('G0', 'G1'), cells=('c0', 'c1'), data=(2, 3), indices=(0, 1),
                indptr=(0, 1, 2)):
    with h5py.File(path, 'w') as h:
        h.create_dataset('genes', data=np.asarray(genes, object), dtype=h5py.string_dtype())
        h.create_dataset('cell_ids', data=np.asarray(cells, object), dtype=h5py.string_dtype())
        h['X/data'], h['X/indices'], h['X/indptr'] = data, indices, indptr
    return path


@pytest.mark.parametrize('change', [dict(counts_kind=None), dict(genome_build=None),
                                  dict(coordinate_system=None), dict(gene_namespace=None),
                                  dict(genome_build='hg19')])
def test_contract_declarations_cannot_be_inferred(change):
    args = dict(counts_kind='raw', **CONTRACT)
    args.update(change)
    with pytest.raises(ValueError):
        validate_input_contract(**args)


@pytest.mark.parametrize('names', [['G', 'g'], ['G', ''], ['G', None], ['G', 'NA']])
def test_ambiguous_or_missing_feature_identities_are_rejected(names):
    with pytest.raises(ValueError):
        validate_feature_ids(names, 'symbol')


def test_feature_namespace_does_not_discard_ensembl_versions():
    assert list(validate_feature_ids(['ENSG000001.2', 'ENSG000001.3'], 'ensembl')) == [
        'ENSG000001.2', 'ENSG000001.3']


def test_coordinate_conventions_produce_identical_internal_intervals(tmp_path):
    one = pd.DataFrame(dict(gene=['g0', 'g1'], chromosome=['1', '1'], start=[1, 11], end=[5, 20]))
    zero = one.assign(start=one.start - 1)
    one.to_csv(tmp_path / 'one.tsv', sep='\t', index=False)
    zero.to_csv(tmp_path / 'zero.tsv', sep='\t', index=False)
    a = read_positions(tmp_path / 'one.tsv', **CONTRACT)
    b = read_positions(tmp_path / 'zero.tsv', **dict(CONTRACT, coordinate_system='zero-based-half-open'))
    pd.testing.assert_frame_equal(a, b)
    assert list(a.start) == [0, 10]


@pytest.mark.parametrize('start,end', [(0, 2), (2, 1), (1.5, 3), (1, np.inf)])
def test_positions_reject_invalid_closed_intervals(tmp_path, start, end):
    pd.DataFrame(dict(gene=['G'], chromosome=['chr1'], start=[start], end=[end])).to_csv(
        tmp_path / 'pos.tsv', sep='\t', index=False)
    with pytest.raises(ValueError):
        read_positions(tmp_path / 'pos.tsv', **CONTRACT)


@pytest.mark.parametrize('kwargs', [dict(indptr=(1, 1, 2)), dict(indptr=(0, 2, 1)),
                                   dict(indices=(0, 2)), dict(indices=(0., 1.)),
                                   dict(data=(1.25, 3)), dict(data=(1.0000000001, 3)),
                                   dict(indptr=(0, 2, 2), indices=(0, 0)), dict(genes=('G', 'g')),
                                   dict(cells=('c0', 'c0')), dict(data=(2 ** 24 + 1, 3))])
def test_invalid_sparse_counts_fail_before_aggregation(tmp_path, kwargs):
    path = counts_file(tmp_path / 'counts.h5', **kwargs)
    with pytest.raises(ValueError):
        with CountSource(path, gene_namespace='symbol') as src:
            src.read(0, 2)


def test_count_chunks_conserve_valid_integer_counts(tmp_path):
    with CountSource(counts_file(tmp_path / 'counts.h5'), gene_namespace='symbol') as src:
        assert sum(X.sum() for _, _, X in src.chunks(1)) == 5
        with pytest.raises(ValueError, match='positive integer'):
            list(src.chunks(0))


def test_10x_uses_declared_feature_namespace_and_rejects_mixed_assays(tmp_path):
    path = tmp_path / 'tenx.h5'
    with h5py.File(path, 'w') as h:
        g = h.create_group('matrix')
        for name, values in [('features/name', ['SAME', 'SAME']),
                             ('features/id', ['ENSG000001', 'ENSG000002']),
                             ('features/feature_type', ['Gene Expression', 'Antibody Capture']),
                             ('barcodes', ['c0', 'c1'])]:
            g.create_dataset(name, data=np.asarray(values, object), dtype=h5py.string_dtype())
        g['data'], g['indices'], g['indptr'], g['shape'] = [2, 3], [0, 1], [0, 1, 2], [2, 2]
    with pytest.raises(ValueError, match='non-Gene Expression'):
        CountSource(path, gene_namespace='ensembl')
    with h5py.File(path, 'r+') as h:
        h['matrix/features/feature_type'][1] = 'Gene Expression'
    with CountSource(path, gene_namespace='ensembl') as src:
        assert list(src.genes) == ['ENSG000001', 'ENSG000002']
        assert src.read(0, 2).sum() == 5
    with pytest.raises(ValueError, match='Duplicated canonical'):
        CountSource(path, gene_namespace='symbol')


def test_bins_and_hmm_restart_at_arm_and_declared_gap(tmp_path):
    pd.DataFrame(dict(chromosome=['chr1'], p_arm_end=[100])).to_csv(tmp_path / 'cen.tsv', sep='\t', index=False)
    starts = np.r_[np.arange(0, 50, 10), np.arange(110, 160, 10), np.arange(1000, 1050, 10)]
    genes = [f'G{i}' for i in range(len(starts))]
    positions = pd.DataFrame(dict(gene=genes, chromosome='chr1', start=starts, end=starts + 5))
    bins, _, used = make_bins(genes, positions, np.ones(len(genes)), bin_genes=5,
                             genome_build='custom', centromeres=tmp_path / 'cen.tsv', max_gap_bp=100)
    assert list(bins.n_genes) == [5, 5, 5]
    assert list(arm_labels(bins, tmp_path / 'cen.tsv', 'custom')) == ['chr1p', 'chr1q', 'chr1q']
    assert chromosome_starts(bins).all()
    assert len(used) == len(genes)
    assert not ((bins.start < 100) & (bins.end > 100)).any()


def test_centromere_spanning_gene_does_not_join_both_arms(tmp_path):
    pd.DataFrame(dict(chromosome=['chr1'], p_arm_end=[100])).to_csv(tmp_path / 'cen.tsv', sep='\t', index=False)
    starts = np.r_[np.arange(0, 50, 10), 95, np.arange(110, 160, 10)]
    ends = starts + 5
    ends[5] = 105
    genes = [f'G{i}' for i in range(len(starts))]
    pos = pd.DataFrame(dict(gene=genes, chromosome='chr1', start=starts, end=ends))
    _, mapping, used = make_bins(genes, pos, np.ones(len(genes)), bin_genes=5, centromeres=tmp_path / 'cen.tsv')
    assert mapping[5] == -1 and genes[5] not in set(used.gene)


def test_missing_patient_identifier_is_rejected(tmp_path):
    pd.DataFrame(dict(cell_id=['c'], specimen=['S'], patient=[''], role=['query'], ref_group=['NA'],
                      ref_split=['NA'])).to_csv(tmp_path / 'cells.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='patient'):
        read_cells(tmp_path / 'cells.tsv', ['c'])


@pytest.fixture
def kinker_identity_scorer():
    import importlib.util
    from pathlib import Path
    path = Path(__file__).resolve().parents[1] / 'benchmarks/kinker/evaluate.py'
    spec = importlib.util.spec_from_file_location('identity_checked_kinker', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_kinker_missing_scored_gene_never_reuses_last_count_column(tmp_path, kinker_identity_scorer):
    from scipy import sparse
    inp, run = tmp_path / 'input', tmp_path / 'run'
    inp.mkdir()
    run.mkdir()
    pd.DataFrame(dict(gene=['MISSING'], chromosome=['chr1'], start=[1], end=[10])).to_csv(
        inp / 'positions.tsv', sep='\t', index=False)
    pd.DataFrame(dict(cell_id=['c'], design=['single_line'], specimen=['S'], truth_line=['L'],
                      dataset=['pool1'], role=['query'], ref_split=[''])).to_csv(
        run / 'cells.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(gene=['MISSING'], bin=[0])).to_csv(run / 'genes.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(bin=[0])).to_csv(run / 'bins.tsv', sep='\t', index=False)
    pd.DataFrame(dict(specimen=['S'], clone=['S_c1'], bin=[0], log_ratio=[.2], log_dose=[.2])).to_csv(
        run / 'profiles.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(specimen=['S'], clone=['S_c1'], normal_like=[False])).to_csv(
        run / 'clones.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='scored genes absent'):
        kinker_identity_scorer.score_run(inp, run, sparse.csr_matrix([[1, 999]]),
                                         ['G1', 'LAST'], ['c'], pd.DataFrame({'L': [0.]}, index=['MISSING']),
                                         pd.Series({'chr1': 100}), 'haplodose')


@pytest.mark.parametrize('source,target', [(['G', 'G'], ['G']), (['G'], ['G', 'G']),
                                         (['G', None], ['G']), (['G'], [''])])
def test_kinker_ambiguous_identities_fail_instead_of_changing_weight(source, target, kinker_identity_scorer):
    with pytest.raises(ValueError, match='unique, nonmissing'):
        kinker_identity_scorer.checked_indexer(source, target, 'genes')


def test_kinker_permuted_identity_alignment_preserves_named_counts(kinker_identity_scorer):
    source = np.array([[7, 2, 11]])
    indexer = kinker_identity_scorer.checked_indexer(['B', 'A', 'C'], ['A', 'C'], 'genes')
    np.testing.assert_array_equal(source[:, indexer], [[2, 11]])


def test_kinker_profile_mean_preserves_cell_population_and_bin_identity(kinker_identity_scorer):
    profile = pd.DataFrame(dict(clone=['A', 'B', 'A', 'B'], bin=[9, 9, 2, 2],
                                log_ratio=[2., -2., 4., 0.]))
    actual, fraction, reason = kinker_identity_scorer.specimen_profile(
        profile.iloc[[3, 0, 2, 1]], pd.Series(['A', 'B', 'A']), [2, 9], 'log_ratio')
    np.testing.assert_allclose(actual, [8/3, 2/3], rtol=0, atol=0)
    assert fraction == 1 and reason == 'complete'


@pytest.mark.parametrize('defect', ['missing_clone', 'missing_bin', 'nan', 'infinite', 'unassigned'])
def test_kinker_incomplete_profile_cannot_become_zero_or_smaller_population(kinker_identity_scorer, defect):
    profile = pd.DataFrame(dict(clone=['A', 'B', 'A', 'B'], bin=[0, 0, 1, 1],
                                log_ratio=[2., -2., 4., 0.]))
    cells = pd.Series(['A', 'A', 'B'])
    if defect == 'missing_clone':
        profile = profile[profile.clone.ne('B')]
    elif defect == 'missing_bin':
        profile = profile.drop(index=3)
    elif defect in ('nan', 'infinite'):
        profile.loc[3, 'log_ratio'] = np.nan if defect == 'nan' else np.inf
    else:
        cells.iloc[2] = np.nan
    result, fraction, reason = kinker_identity_scorer.specimen_profile(profile, cells, [0, 1], 'log_ratio')
    assert result.isna().all()
    assert fraction == 2/3 and reason == 'incomplete_population_or_profile'


def test_kinker_duplicate_clone_bin_is_ambiguous(kinker_identity_scorer):
    profile = pd.DataFrame(dict(clone=['A', 'A'], bin=[0, 0], log_ratio=[1., 2.]))
    with pytest.raises(ValueError, match='duplicate clone/bin'):
        kinker_identity_scorer.specimen_profile(profile, pd.Series(['A']), [0], 'log_ratio')


def test_kinker_comparator_family_uses_same_population(kinker_identity_scorer):
    lines = pd.DataFrame(dict(rho_A=[.2, .3, .4], rho_B=[.1, .3, .2], rho_C=[np.nan, .3, .2]))
    population, tests = kinker_identity_scorer.paired_method_tests(lines, ['A', 'B', 'C'])
    assert population.index.tolist() == [1, 2]
    assert tests.n.tolist() == [2, 2, 2]
    assert set(zip(tests.a, tests.b)) == {('A', 'B'), ('A', 'C'), ('B', 'C')}
    assert tests.loc[tests.a.eq('B') & tests.b.eq('C'), 'p'].iloc[0] == 1


def test_kinker_unavailable_comparator_keeps_family_without_a_verdict(kinker_identity_scorer):
    lines = pd.DataFrame(dict(rho_A=[.2, .3], rho_B=[.1, .4], rho_C=[np.nan, np.nan]))
    population, tests = kinker_identity_scorer.paired_method_tests(lines, ['A', 'B', 'C'])
    assert len(population) == 0 and len(tests) == 3
    assert tests.n.eq(0).all() and tests.p.isna().all() and tests.p_holm.isna().all()


def test_kinker_unavailable_population_has_no_sign_score(tmp_path, kinker_identity_scorer):
    from scipy import sparse
    inp, run = tmp_path/'input', tmp_path/'run'
    inp.mkdir(); run.mkdir()
    genes = [f'G{i}' for i in range(5)]
    pd.DataFrame(dict(gene=genes, chromosome=['1']*5, start=range(5))).to_csv(
        inp/'positions.tsv', sep='\t', index=False)
    pd.DataFrame(dict(cell_id=['q1', 'q2', 'r'], design=['single_line', 'single_line', 'reference'],
        specimen=['S', 'S', 'R'], truth_line=['L', 'L', 'REF'], dataset=['D']*3,
        role=['query', 'query', 'reference'], ref_split=['', '', 'fit'], line=['L', 'L', 'REF'],
        umi=[5]*3, clone=['A', 'B', ''], malignant_call=[True, True, False])).to_csv(
            run/'cells.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(gene=genes, bin=[0]*5)).to_csv(run/'genes.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(bin=[0])).to_csv(run/'bins.tsv', sep='\t', index=False)
    pd.DataFrame(dict(specimen=['S'], clone=['A'], bin=[0], log_ratio=[1.], log_dose=[1.])).to_csv(
        run/'profiles.tsv.gz', sep='\t', index=False)
    pd.DataFrame(dict(specimen=['S'], clone=['A'], normal_like=[False])).to_csv(
        run/'clones.tsv', sep='\t', index=False)
    scores, _ = kinker_identity_scorer.score_run(inp, run, sparse.csr_matrix(np.ones((3, 5))),
        genes, ['q1', 'q2', 'r'], pd.DataFrame(dict(L=np.ones(5), REF=np.zeros(5)), index=genes),
        pd.Series({'chr1': 100}), 'haplodose')
    assert np.isnan(scores.sign_haplodose_log_ratio.iloc[0])
    assert scores.profile_fraction_haplodose_log_ratio.iloc[0] == .5
