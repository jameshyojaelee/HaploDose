"""split_mode 'thin' (NB count splitting) against the default even/odd bin split."""
import json

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from haplodose.clones import CloneParams, Rows, thin_counts, thin_gene_counts
from haplodose.pipeline import run
from haplodose.simulate import write_dataset


INPUT_CONTRACT = dict(counts_kind='raw', genome_build='hg38',
                      coordinate_system='one-based-closed', gene_namespace='symbol')


def _counts(seed=0, n=300, J=40):
    rng = np.random.default_rng(seed)
    m = rng.gamma(2.0, 2.0, J)
    return rng.negative_binomial(4.0, 4.0 / (4.0 + m[None] * rng.uniform(0.5, 2.0, (n, 1))), (n, J)).astype(
        np.float32)


def test_thinning_halves_sum_and_are_deterministic_and_chunk_free():
    Y = _counts()
    rows = np.arange(10, 10 + 2 * len(Y), 2)                  # a specimen's rows inside a larger matrix
    full = np.zeros((rows.max() + 1, Y.shape[1]), np.float32)
    full[rows] = Y
    size = np.full(Y.shape[1], 3.0)
    a, b = Rows(full, rows, chunk=7), Rows(full, rows, chunk=4096)
    for r in (a, b):
        r.set_thinning(11, size)
    prop = np.vstack([Yc for _, _, Yc in a.chunks('prop')])
    test = np.vstack([Yc for _, _, Yc in a.chunks('test')])
    assert np.array_equal(prop + test, Y)
    assert (prop >= 0).all() and (test >= 0).all()
    assert np.array_equal(prop, np.vstack([Yc for _, _, Yc in b.chunks('prop')]))
    idx = np.array([250, 3, 77])
    assert np.array_equal(a.take(idx, 'prop'), prop[idx])
    assert np.array_equal(b.take(idx, 'test'), test[idx])
    R = sparse.random(len(Y), 3, density=0.3, random_state=1, format='csr')
    assert np.allclose(a.weighted_sum(R, 'test'), R.T @ test.astype(np.float64))
    assert np.allclose(b.weighted_sum(R.toarray(), 'prop'), R.T @ prop.astype(np.float64))
    assert np.array_equal(a.take(np.arange(len(Y))), Y)            # no fold: full counts
    other = Rows(full, rows)
    other.set_thinning(12, size)
    assert not np.array_equal(other.take(np.arange(len(Y)), 'prop'), prop)


def test_thinned_halves_have_nb_moments():
    """NB(m, b) splits into NB(m/2, b/2) halves: means m/2 and variance m/2 + m^2 / (2 b), nearly uncorrelated."""
    rng = np.random.default_rng(3)
    m, b, n = 8.0, 2.0, 40000
    Y = rng.negative_binomial(b, b / (b + m), (n, 1)).astype(np.float32)
    prop = thin_counts(Y, np.arange(n), np.array([b]), seed=5)[:, 0]
    test = Y[:, 0] - prop
    for h in (prop, test):
        assert abs(h.mean() - m / 2) < 0.1
        assert abs(h.var() / (m / 2 + m * m / (2 * b)) - 1) < 0.05
    assert abs(np.corrcoef(prop, test)[0, 1]) < 0.02


def test_gene_test_half_nested_in_bin_split():
    rng = np.random.default_rng(4)
    n, G = 50, 30
    gbin = np.repeat(np.arange(6), 5)
    X = rng.poisson(2.0, (n, G)).astype(np.float64)
    Ybin = np.zeros((n, 6))
    np.add.at(Ybin.T, gbin, X.T)
    Ytest = Ybin - thin_counts(Ybin.astype(np.float32), np.arange(n), np.full(6, 5.0), seed=1)
    Xt = thin_gene_counts(sparse.csr_matrix(X), gbin, Ybin, Ytest, seed=1).toarray()
    assert (Xt >= 0).all() and (Xt <= X).all()
    got = np.zeros((n, 6))
    np.add.at(got.T, gbin, Xt.T)
    assert np.array_equal(got, Ytest)
    assert np.array_equal(Xt, thin_gene_counts(sparse.csr_matrix(X), gbin, Ybin, Ytest, seed=1).toarray())


@pytest.fixture(scope='module')
def runs(tmp_path_factory):
    root = tmp_path_factory.mktemp('split')
    data = write_dataset(root / 'data')
    out = {}
    for name, params in [('default', None), ('evenodd', CloneParams(split_mode='evenodd')),
                         ('thin', CloneParams(split_mode='thin'))]:
        run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', root / name, chunk=97, params=params, **INPUT_CONTRACT)
        out[name] = root / name
    return out


def test_default_is_evenodd_and_unchanged(runs):
    a, b = runs['default'], runs['evenodd']
    assert np.array_equal(np.load(a / 'posteriors.npy'), np.load(b / 'posteriors.npy'), equal_nan=True)
    for f in ('cells.tsv.gz', 'differences.tsv', 'segments.tsv'):
        pd.testing.assert_frame_equal(pd.read_csv(a / f, sep='\t'), pd.read_csv(b / f, sep='\t'))
    assert json.loads((a / 'qc.json').read_text())['params']['split_mode'] == 'evenodd'


def test_thin_records_mode_and_recovers_clones(runs):
    out = runs['thin']
    assert json.loads((out / 'qc.json').read_text())['params']['split_mode'] == 'thin'
    cells = pd.read_csv(out / 'cells.tsv.gz', sep='\t')
    s1 = cells[cells.specimen.eq('S1')]
    normal = s1.truth_clone.eq('normal')
    assert s1.loc[normal, 'malignant_call'].mean() < 0.1
    assert s1.loc[~normal, 'malignant_call'].mean() > 0.9
    s2 = cells[cells.specimen.eq('S2') & cells.malignant_call]
    assert s2.clone.nunique() == 1
    # Same bar as test_pipeline.test_subclone_assignment_above_chance (one 100-gene chr12 gain).
    sub = s1[s1.truth_clone.isin(['clone1', 'clone2']) & s1.malignant_call]
    assert sub.clone.nunique() == 2
    table = pd.crosstab(sub.truth_clone, sub.clone)
    accuracy = max(table.values.diagonal().sum(), table.values[::-1].diagonal().sum()) / table.values.sum()
    assert accuracy > 0.65
    diff = pd.read_csv(out / 'differences.tsv', sep='\t')
    d1 = diff[diff.specimen.eq('S1') & ~diff.clone_a.str.endswith('_c0') & ~diff.clone_b.str.endswith('_c0')]
    assert set(d1.chromosome) == {'chr12'}
