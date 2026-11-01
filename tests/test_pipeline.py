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
