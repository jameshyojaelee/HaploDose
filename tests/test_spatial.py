"""Potts prior: helps when clones are spatially coherent, harmless when they are not."""
import numpy as np
import pandas as pd
import pytest

from haplodose.clones import CloneParams
from haplodose.io import write_table
from haplodose.pipeline import run
from haplodose.simulate import write_dataset


def accuracy(out):
    cells = pd.read_csv(out / 'cells.tsv.gz', sep='\t')
    t = cells[cells.truth_clone.isin(['clone1', 'clone2'])]
    table = pd.crosstab(t.truth_clone, t.clone)
    if table.shape[1] < 2:
        return 0.5
    best = 0
    for a in table.columns:
        for b in table.columns:
            if a != b:
                best = max(best, table.loc['clone1', a] + table.loc['clone2', b])
    return best / len(t)


@pytest.fixture(scope='module')
def spatial(tmp_path_factory):
    root = tmp_path_factory.mktemp('sp')
    spec = [dict(name='V', n_cells=800, normal_frac=0.1, clones=[
        dict(frac=0.5, events=[('chr3', 0.0, 1.0, 3)]),
        dict(frac=0.5, events=[('chr3', 0.0, 1.0, 3), ('chr12', 0.0, 1.0, 3)])])]
    data = write_dataset(root / 'data', specimens=spec, depth=800.0, spatial=True)
    cells = pd.read_csv(data / 'cells.tsv', sep='\t', keep_default_na=False)
    shuffled = cells.copy()
    spec_rows = np.flatnonzero(cells.specimen.eq('V'))
    perm = np.random.default_rng(1).permutation(spec_rows)
    shuffled.loc[spec_rows, ['x', 'y']] = cells.loc[perm, ['x', 'y']].to_numpy()
    write_table(shuffled, root / 'shuffled_cells.tsv')
    args = (data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv')
    run(*args, root / 'potts')
    run(*args, root / 'plain', params=CloneParams(spatial_neighbors=0))
    run(data / 'counts.h5', root / 'shuffled_cells.tsv', data / 'positions.tsv', root / 'shuffled')
    return root


def test_potts_improves_coherent_clones(spatial):
    a_potts, a_plain = accuracy(spatial / 'potts'), accuracy(spatial / 'plain')
    lam = pd.read_csv(spatial / 'potts' / 'specimens.tsv', sep='\t').potts_lambda.iloc[0]
    assert lam > 0
    assert a_potts >= a_plain + 0.05 and a_potts > 0.85, (a_potts, a_plain)


def test_potts_harmless_without_spatial_structure(spatial):
    lam = pd.read_csv(spatial / 'shuffled' / 'specimens.tsv', sep='\t').potts_lambda.iloc[0]
    assert lam <= 0.5
    assert accuracy(spatial / 'shuffled') >= accuracy(spatial / 'plain') - 0.03
