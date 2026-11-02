"""End-to-end on simulated specimens with known clones."""
import os

import numpy as np
import pandas as pd
import pytest

from haplodose.io import guard, new_output_dir
from haplodose.pipeline import run
from haplodose.simulate import write_dataset


@pytest.fixture(scope='module')
def result(tmp_path_factory):
    root = tmp_path_factory.mktemp('hd')
    data = write_dataset(root / 'data')
    (root / 'mask.txt').write_text('IGF2\n')
    qc = run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', root / 'out',
             exclude_genes=root / 'mask.txt', chunk=97)
    cells = pd.read_csv(root / 'out/cells.tsv.gz', sep='\t')
    return root, qc, cells


def test_recovers_clones_and_normals(result):
    root, qc, cells = result
    s1 = cells[cells.specimen.eq('S1')]
    normal = s1.truth_clone.eq('normal')
    assert s1.loc[normal, 'malignant_call'].mean() < 0.1
    assert s1.loc[~normal, 'malignant_call'].mean() > 0.9
    s2 = cells[cells.specimen.eq('S2') & cells.malignant_call]
    assert s2.clone.nunique() == 1


def test_subclone_assignment_above_chance(result):
    """clone2 differs from clone1 by one 100-gene gain. At ~3000 UMIs per cell the
    region holds ~70 expected UMIs, so per-cell separation is limited (d' ~ 1.6,
    ~80% accuracy at best). Require the split and assignment clearly above chance."""
    root, qc, cells = result
    s1 = cells[cells.specimen.eq('S1') & cells.truth_clone.isin(['clone1', 'clone2']) & cells.malignant_call]
    assert s1.clone.nunique() == 2
    table = pd.crosstab(s1.truth_clone, s1.clone)
    accuracy = max(table.values.diagonal().sum(), table.values[::-1].diagonal().sum()) / table.values.sum()
    assert accuracy > 0.65


def test_subclonal_difference_reported(result):
    root, _, _ = result
    diff = pd.read_csv(root / 'out/differences.tsv', sep='\t')
    s1 = diff[diff.specimen.eq('S1') & ~diff.clone_a.str.endswith('_c0') & ~diff.clone_b.str.endswith('_c0')]
    assert set(s1.chromosome) == {'chr12'}


def test_masked_gene_and_outputs(result):
    root, qc, cells = result
    genes = pd.read_csv(root / 'out/genes.tsv.gz', sep='\t')
    assert 'IGF2' not in set(genes.gene)
    seg = pd.read_csv(root / 'out/segments.tsv', sep='\t')
    assert set(seg.state) <= {'loss', 'neutral', 'gain', 'amplification'}
    P = np.load(root / 'out/posteriors.npy')
    assert P.shape[0] == len(cells)
    ref = pd.read_csv(root / 'out/reference_check.tsv', sep='\t')
    # Flat prior; expected false rate for these profiles is about Phi(-d'/2) ~ 0.05.
    assert (ref.fraction_called_aneuploid < 0.10).all()


def test_guard_and_nonempty_output(tmp_path, monkeypatch):
    protected = tmp_path / 'raw'
    protected.mkdir()
    (tmp_path / 'link').symlink_to(protected)
    monkeypatch.setenv('HAPLODOSE_PROTECTED_PATHS', str(protected))
    with pytest.raises(PermissionError):
        guard(tmp_path / 'link' / 'new' / 'file.tsv')
    out = tmp_path / 'out'
    out.mkdir()
    (out / 'x').write_text('1')
    with pytest.raises(FileExistsError):
        new_output_dir(out)


def test_em_drops_unsupported_clone():
    """A clone that attracts no cells is dropped without breaking later updates."""
    from haplodose._native import bin_sums
    from haplodose.bins import chromosome_starts, make_bins
    from haplodose.clones import CloneParams, Rows, SpecimenModel
    from haplodose.simulate import simulate
    X, genes, cells, positions, _, _ = simulate(
        ref_cells_per_type=100, specimens=[dict(name='S', n_cells=300, normal_frac=0.0,
                                                clones=[dict(frac=1.0, events=[('chr7', 0.0, 1.0, 3)])])])
    bins, gene_bin, _ = make_bins(genes, positions, np.asarray(X.mean(axis=0)).ravel())
    Y, _ = bin_sums(X, gene_bin, len(bins))
    L = Y.sum(axis=1).astype(np.float64)
    ref = cells.ref_split.eq('fit').to_numpy()
    mu = np.stack([Y[ref & cells.ref_group.eq(t).to_numpy()].sum(0) for t in ['type0', 'type1', 'type2']], 1)
    mu = (mu + 0.1) / (mu + 0.1).sum(0)
    rows = np.flatnonzero(cells.specimen.eq('S').to_numpy())
    reset = chromosome_starts(bins)
    model = SpecimenModel(Rows(Y, rows), L[rows], mu, reset, np.full(len(bins), 0.05), CloneParams())
    chr7 = (bins.chromosome == 'chr7').to_numpy()
    real = np.where(chr7, model.beta * np.log(1.5), 0.0)
    fake = np.where((bins.chromosome == 'chr2').to_numpy(), -2.0, 0.0)   # nobody looks like this
    doses, pi, R, *_ = model.em([np.zeros(len(bins)), real, fake], np.array([0.1, 0.6, 0.3]))
    assert len(doses) == 2 and R.shape[1] == 2
    assert R[:, 1].mean() > 0.9


def test_depth_does_not_split_clones(result):
    """S3: one clone, half the cells at a quarter of the depth. Must stay one clone."""
    root, _, cells = result
    s3 = cells[cells.specimen.eq('S3') & cells.malignant_call]
    assert s3.clone.nunique() == 1
    assert len(s3) > 0.9 * cells.specimen.eq('S3').sum()


def test_reference_check_near_expected(result):
    root, _, _ = result
    ref = pd.read_csv(root / 'out/reference_check.tsv', sep='\t')
    assert (ref.fraction_called_aneuploid <= ref.expected_fraction + 0.05).all()


def test_expression_programme_is_not_a_clone(tmp_path):
    """Gene clusters switched on in 25-30% of all cells (normals too) must not become clones.

    Chromosomes of 1000 genes and 500 fit cells per type, near real arm sizes and SMC
    references. On short arms a programme's arm mean is left to dosage by design, so a
    20-bin chromosome would not test the programme model (see haplodose.programmes)."""
    spec = [dict(name='P', n_cells=500, normal_frac=0.1,
                 clones=[dict(frac=1.0, events=[('chr1', 0.0, 1.0, 3), ('chr4', 0.0, 0.5, 1)])])]
    prog = [dict(chromosome='chr2', start=0.40, end=0.45, fold=2.5, frac=0.3),
            dict(chromosome='chr3', start=0.60, end=0.75, fold=1.8, frac=0.25)]
    data = write_dataset(tmp_path / 'data', n_chrom=6, genes_per_chrom=1000, ref_cells_per_type=1000,
                         specimens=spec, programmes=prog)
    run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', tmp_path / 'out')
    clones = pd.read_csv(tmp_path / 'out/clones.tsv', sep='\t')
    assert (~clones.normal_like).sum() == 1
    ref = pd.read_csv(tmp_path / 'out/reference_check.tsv', sep='\t')
    assert (ref.fraction_called_aneuploid <= ref.expected_fraction + 0.05).all()
    assert clones.programme_overlap.dropna().max() < 0.5


def test_arms_table(result):
    """arms.tsv covers every bin once per clone and follows S2's chr7 gain and chr18 loss."""
    root, _, _ = result
    arms = pd.read_csv(root / 'out/arms.tsv', sep='\t')
    n_bins = len(pd.read_csv(root / 'out/bins.tsv', sep='\t'))
    assert (arms.groupby(['specimen', 'clone']).n_bins.sum() == n_bins).all()
    s2 = arms[arms.specimen.eq('S2') & ~arms.clone.str.endswith('_c0')]
    assert s2.clone.nunique() == 1
    chrom = s2.arm.str.extract(r'^(chr\w+?)[pq]$', expand=False)
    r = s2.groupby(chrom).median_log_ratio.median()
    assert r['chr7'] > 0.1 and r['chr18'] < -0.1 and abs(r['chr1']) < 0.1


def test_reference_depth_slopes(result):
    """Depth slopes from reference cells, held fixed: S3's depth mix still gives one clone and S1's
    normal cells stay normal."""
    from haplodose.clones import CloneParams
    root, _, _ = result
    data = root / 'data'
    run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', root / 'out_refdepth',
        exclude_genes=root / 'mask.txt', chunk=97, params=CloneParams(depth_source='reference'))
    cells = pd.read_csv(root / 'out_refdepth/cells.tsv.gz', sep='\t')
    s3 = cells[cells.specimen.eq('S3') & cells.malignant_call]
    assert s3.clone.nunique() == 1
    assert len(s3) > 0.9 * cells.specimen.eq('S3').sum()
    s1 = cells[cells.specimen.eq('S1')]
    normal = s1.truth_clone.eq('normal')
    assert s1.loc[normal, 'malignant_call'].mean() < 0.1
    assert s1.loc[~normal, 'malignant_call'].mean() > 0.9


def test_gene_level_log_ratio_ignores_one_dominant_gene():
    """A bin of 10 genes where one highly expressed gene is 4x off: the summed ratio follows it, the
    gene median does not; genes below the expected-count floor are left out."""
    from haplodose.pipeline import gene_level_log_ratio, gene_shares
    gbin = np.repeat([0, 1], 10)
    ref = np.r_[1000.0, np.full(9, 20.0), np.full(10, 50.0)]
    share = gene_shares(ref[None], [1.0], gbin, 2)
    E_bin = np.array([1180.0, 500.0])
    y = E_bin[gbin] * share
    y[0] *= 4.0                                           # the dominant gene, not copy number
    summed = np.log((np.bincount(gbin, weights=y) + 0.5) / (E_bin + 0.5))
    out = gene_level_log_ratio(y, E_bin, share, gbin, summed - summed[1])
    assert summed[0] - summed[1] > 0.9 and abs(out[0] - out[1]) < 0.05
    low = gene_level_log_ratio(y * 1e-3, E_bin * 1e-3, share, gbin, np.array([0.3, 0.0]))
    assert np.allclose(low, [0.3, 0.0])                   # every gene below the floor: summed value kept


def test_gene_level_columns_reported(result):
    """log_ratio_gene sits next to log_ratio and follows S2's chr7 gain and chr18 loss (as test_arms_table)."""
    root, _, _ = result
    prof = pd.read_csv(root / 'out/profiles.tsv.gz', sep='\t')
    arms = pd.read_csv(root / 'out/arms.tsv', sep='\t')
    assert {'log_ratio', 'log_ratio_gene'} <= set(prof.columns) and 'median_log_ratio_gene' in arms.columns
    assert prof.log_ratio_gene.notna().all()
    s2 = arms[arms.specimen.eq('S2') & ~arms.clone.str.endswith('_c0')]
    r = s2.groupby(s2.arm.str.extract(r'^(chr\w+?)[pq]$', expand=False)).median_log_ratio_gene.median()
    assert r['chr7'] > 0.1 and r['chr18'] < -0.1 and abs(r['chr1']) < 0.1


def test_gene_shares_follow_fitted_reference_weights():
    """Review 2026-11-02 finding 3: the query is exactly reference type A; type B has the same bin totals but a
    different gene composition in bin 2. Shares from A's fitted weight give a flat gene profile; shares from
    raw pooling of A and B report a large false loss (about -1.0)."""
    from haplodose.pipeline import gene_level_log_ratio, gene_shares
    gbin = np.repeat(np.arange(5), 3)
    A = np.full(15, 100.0)
    A[6:9] = [240, 30, 30]
    B = A.copy()
    B[6:9] = [30, 135, 135]
    E_bin = np.bincount(gbin, weights=A)                  # query = type A, so bin expectation = A's bin totals
    r_bin = np.zeros(5)
    fitted = gene_level_log_ratio(A, E_bin, gene_shares(np.vstack([A, B]), [1.0, 0.0], gbin, 5), gbin, r_bin)
    pooled = gene_level_log_ratio(A, E_bin, gene_shares(np.vstack([A, B]), [0.5, 0.5], gbin, 5), gbin, r_bin)
    assert np.abs(fitted).max() < 0.02
    assert pooled[2] < -0.9


def test_assignment_dose_reported(result):
    """profiles.tsv carries the dose used for cell assignment next to the clone's own fitted log_dose; the
    normal-like clone's assignment dose is zero (review 2026-11-02 #2)."""
    root, _, _ = result
    prof = pd.read_csv(root / 'out/profiles.tsv.gz', sep='\t')
    assert prof.assignment_dose.notna().all()
    assert (prof[prof.clone.str.endswith('_c0')].assignment_dose == 0).all()
    s1 = prof[prof.specimen.eq('S1') & ~prof.clone.str.endswith('_c0')]
    assert s1.groupby('clone').assignment_dose.apply(lambda d: d.abs().max()).gt(0.05).all()
