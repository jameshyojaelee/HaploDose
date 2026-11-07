"""Regressions for scientific code review 2026-11-03 findings F1, F7 and F8.

Built from the counterexamples in dev/scientific_review_20261103/probes.py."""
import h5py
import numpy as np
import pandas as pd
import pytest

from haplodose.alleles import read_alleles
from haplodose.bulk import run_bulk
from haplodose.clones import CloneParams, Rows, SpecimenModel, thin_counts
from haplodose.io import CountSource
from haplodose.pipeline import run
from haplodose.programmes import Programmes
from haplodose.simulate import write_dataset


INPUT_CONTRACT = dict(counts_kind='raw', genome_build='hg38',
                      coordinate_system='one-based-closed', gene_namespace='symbol')


def test_f1_dropping_confident_small_clone_keeps_responsibilities_finite():
    """Ten high-depth cells confidently in a clone below min_cells: after the drop their responsibilities over
    the surviving clone underflow to 0 and plain renormalisation gave 0/0."""
    J, n = 40, 40
    Y = np.full((n, J), 10.0, np.float32)
    Y[:10, :20], Y[:10, 20:] = 1000., 100.
    L = Y.sum(axis=1, dtype=np.float64)
    reset = np.zeros(J, bool)
    reset[[0, 20]] = True
    W = np.tile([0.1, -0.1], J // 2)[None]
    prog = Programmes(W=W, prior_prec=np.ones(1), phi=np.full(J, .001), depth=400., variance_share=.1)
    params = CloneParams(min_cells=30, outer_rounds=1, n_programmes=1)
    model = SpecimenModel(Rows(Y, np.arange(n)), L, np.full((J, 1), 1 / J), reset,
                          np.full(J, .001), params, programmes=prog, gamma=np.zeros(J))
    doses = [np.zeros(J), np.r_[np.full(20, .6 * np.log(1.5)), np.full(20, .6 * np.log(.5))]]
    with np.errstate(invalid='ignore', divide='ignore'):
        _, _, R, LL, *_ = model.em(doses, np.array([.75, .25]))
    assert any(e['action'] == 'drop_small' for e in model.log)
    assert np.isfinite(model.U).all()
    assert np.isfinite(LL).all()
    assert np.isfinite(R).all()
    np.testing.assert_allclose(R.sum(axis=1), 1.0, rtol=1e-5)


def test_f7_reader_rejects_fractional_counts(tmp_path):
    with h5py.File(tmp_path / 'fractional.h5', 'w') as h:
        h['genes'] = np.array(['g'], dtype=h5py.string_dtype())
        h['cell_ids'] = np.array(['c'], dtype=h5py.string_dtype())
        h['X/data'], h['X/indices'], h['X/indptr'] = [0.6], [0], [0, 1]
    with CountSource(tmp_path / 'fractional.h5') as src, pytest.raises(ValueError, match='integer'):
        src.read(0, 1)
    with pytest.raises(ValueError, match='integer'):
        thin_counts(np.full((100, 1), 0.6, np.float32), np.arange(100), np.array([5.0]), seed=1)


def test_f7_bulk_and_allele_inputs_are_validated(tmp_path):
    genes = [f'G{i}' for i in range(5)]
    pd.DataFrame({'gene': genes, 'chromosome': 'chr1', 'start': np.arange(5) * 1000,
                  'end': np.arange(5) * 1000 + 500}).to_csv(tmp_path / 'pos.tsv', sep='\t', index=False)
    pd.DataFrame({'sample': ['N0', 'N1', 'N2', 'T0'], 'role': ['normal'] * 3 + ['tumour']}).to_csv(
        tmp_path / 'samples.tsv', sep='\t', index=False)
    good = pd.DataFrame({s: np.full(5, 10.0) for s in ['N0', 'N1', 'N2', 'T0']}, index=pd.Index(genes, name='gene'))
    for name, edit, msg in [('neg', ('T0', 0, -3.0), 'finite and nonnegative'),
                            ('zero', ('T0', slice(None), 0.0), 'zero total')]:
        X = good.copy()
        X.iloc[edit[1], X.columns.get_loc(edit[0])] = edit[2]
        X.to_csv(tmp_path / f'{name}.tsv', sep='\t')
        with pytest.raises(ValueError, match=msg):
            run_bulk(tmp_path / f'{name}.tsv', tmp_path / 'samples.tsv', tmp_path / 'pos.tsv', tmp_path / f'out_{name}', **dict(INPUT_CONTRACT, coordinate_system='zero-based-half-open'))
    pd.DataFrame({'cell': ['c'], 'CHROM': ['1'], 'POS': [100], 'cM': [0.1], 'AD': [0.5], 'DP': [2.5],
                  'GT': ['0|1'], 'gene': ['G0']}).to_csv(tmp_path / 'alleles.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='integer'):
        read_alleles(tmp_path / 'alleles.tsv', ['c'])


def test_f8_specimen_spanning_two_patients_is_rejected(tmp_path):
    data = write_dataset(tmp_path / 'data')
    cells = pd.read_csv(data / 'cells.tsv', sep='\t', dtype=str, keep_default_na=False)
    s1 = np.flatnonzero(cells.specimen.eq('S1').to_numpy())
    cells.loc[s1[::2], 'patient'] = 'OTHER_PATIENT'
    cells.to_csv(data / 'cells.tsv', sep='\t', index=False)
    with pytest.raises(ValueError, match='more than one patient or dataset'):
        run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', tmp_path / 'out', chunk=97, **INPUT_CONTRACT)
