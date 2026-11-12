"""Actual pipeline isolation across confirmation cells and independent datasets."""
import json
from types import SimpleNamespace

import h5py
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from haplodose.clones import CloneParams
from haplodose.pipeline import confirm_specimen, run
from haplodose.simulate import write_dataset


CONTRACT = dict(counts_kind='raw', genome_build='hg38', coordinate_system='one-based-closed', gene_namespace='symbol')


def write_counts(path, X, genes, ids):
    X = sparse.csr_matrix(X)
    with h5py.File(path, 'x') as handle:
        handle['genes'] = np.asarray(genes, dtype='S')
        handle['cell_ids'] = np.asarray(ids, dtype='S')
        handle['X/data'], handle['X/indices'], handle['X/indptr'] = X.data, X.indices, X.indptr


@pytest.fixture(scope='module')
def split_run(tmp_path_factory):
    root = tmp_path_factory.mktemp('pipeline_roles')
    data = write_dataset(root / 'data', n_chrom=4, genes_per_chrom=60, extra_genes=(),
                         n_types=2, ref_cells_per_type=60, depth=1200.,
                         specimens=[dict(name='A', n_cells=160, normal_frac=.1,
                                         clones=[dict(frac=1., events=[('chr3', 0., 1., 3)])])])
    meta = pd.read_csv(data / 'cells.tsv', sep='\t', dtype=str, keep_default_na=False)
    meta['dataset'] = 'A'
    query = meta.role.ne('reference')
    meta['analysis_split'] = 'discovery'
    meta.loc[np.flatnonzero(query)[::2], 'analysis_split'] = 'confirmation'
    meta.to_csv(root / 'cells.tsv', sep='\t', index=False)
    params = CloneParams(min_cells=10, n_programmes=0, n_pcs=5, n_neighbors=5,
                         max_groups=5, outer_rounds=2, cycles=1, spatial_neighbors=0)
    run(data / 'counts.h5', root / 'cells.tsv', data / 'positions.tsv', root / 'original',
        params=params, bin_genes=5, **CONTRACT)
    return root, data, meta, params


def test_confirmation_changes_do_not_change_discovery_fit(split_run):
    root, data, meta, params = split_run
    with h5py.File(data / 'counts.h5') as handle:
        X = sparse.csr_matrix((handle['X/data'][:], handle['X/indices'][:], handle['X/indptr'][:]),
                              shape=(len(meta), len(handle['genes'])))
        genes = handle['genes'].asstr()[:]
    factors = np.where(meta.analysis_split.eq('confirmation'), 7, 1)
    changed = X.multiply(factors[:, None]).tocsr()
    write_counts(root / 'changed.h5', changed, genes, meta.cell_id)
    run(root / 'changed.h5', root / 'cells.tsv', data / 'positions.tsv', root / 'changed',
        params=params, bin_genes=5, **CONTRACT)
    for name in ['bins.tsv', 'clones.tsv', 'profiles.tsv.gz', 'differences.tsv', 'specimens.tsv',
                 'type_mismatch_genes.tsv']:
        pd.testing.assert_frame_equal(pd.read_csv(root / 'original' / name, sep='\t'),
                                      pd.read_csv(root / 'changed' / name, sep='\t'))
    a = pd.read_csv(root / 'original/cells.tsv.gz', sep='\t')
    b = pd.read_csv(root / 'changed/cells.tsv.gz', sep='\t')
    discovery = a.analysis_split.eq('discovery')
    for column in ['clone', 'clone_posterior', 'p_normal', 'evidence_status']:
        pd.testing.assert_series_equal(a.loc[discovery, column], b.loc[discovery, column])
    manifests = [json.loads((root / out / 'confirmation.json').read_text()) for out in ('original', 'changed')]
    assert manifests[0]['specimens'][0]['artifact']['sha256'] == manifests[1]['specimens'][0]['artifact']['sha256']
    assert a.loc[a.analysis_split.eq('confirmation'), 'clone_posterior'].isna().all()
    assert not a.loc[a.analysis_split.eq('confirmation'), 'malignant_call'].any()
    confirmation = a.analysis_split.eq('confirmation')
    assert a.loc[~confirmation, 'confirmation_unresolved_reason'].isna().all()
    assigned = confirmation & a.status.eq('confirmation_candidate')
    assert a.loc[assigned, 'confirmation_unresolved_reason'].isna().all()
    unresolved = confirmation & a.status.eq('confirmation_unresolved')
    assert a.loc[unresolved, 'confirmation_unresolved_reason'].notna().all()


def test_unrelated_dataset_does_not_change_discovery_or_confirmation_artifact(split_run):
    root, data, meta, params = split_run
    other = meta.copy()
    for column in ('cell_id', 'specimen', 'patient'):
        other[column] = 'B_' + other[column]
    other['dataset'] = 'B'
    combined = pd.concat([meta, other], ignore_index=True)
    combined.to_csv(root / 'combined_cells.tsv', sep='\t', index=False)
    with h5py.File(data / 'counts.h5') as handle:
        genes = handle['genes'].asstr()[:]
        X = sparse.csr_matrix((handle['X/data'][:], handle['X/indices'][:], handle['X/indptr'][:]),
                              shape=(len(meta), len(genes)))
    write_counts(root / 'combined.h5', sparse.vstack([X, X * 3]), genes, combined.cell_id)
    run(root / 'combined.h5', root / 'combined_cells.tsv', data / 'positions.tsv', root / 'combined',
        params=params, bin_genes=5, **CONTRACT)
    datasets = json.loads((root / 'combined/datasets.json').read_text())
    child = next(item for item in datasets if item['dataset'] == 'A')['output']
    for name in ['bins.tsv', 'clones.tsv', 'profiles.tsv.gz', 'differences.tsv']:
        pd.testing.assert_frame_equal(pd.read_csv(root / 'original' / name, sep='\t'),
                                      pd.read_csv(f'{child}/{name}', sep='\t'))
    a = json.loads((root / 'original/confirmation.json').read_text())
    b = json.loads(open(f'{child}/confirmation.json').read())
    assert a['specimens'][0]['artifact']['sha256'] == b['specimens'][0]['artifact']['sha256']


def test_fraction_and_explicit_metadata_are_not_silently_combined(split_run):
    root, data, _, params = split_run
    with pytest.raises(ValueError, match='not both'):
        run(data / 'counts.h5', root / 'cells.tsv', data / 'positions.tsv', root / 'invalid',
            params=params, discovery_fraction=.5, **CONTRACT)


def test_confirmation_reasons_survive_chunk_boundaries(tmp_path):
    # A frozen fit isolates the export path without fitting to confirmation data.
    rates = np.log([[.4, .4, .1, .1], [.1, .1, .4, .4]])
    model = SimpleNamespace(phi=np.full(4, .1), gamma=np.zeros(4), programmes=None,
                            log_rates=lambda doses: rates)
    res = dict(model=model, R=np.array([[.9, .1], [.1, .9]]), pi=np.array([.5, .5]), doses=None)
    Y = np.tile([80, 8, 20, 2], (4099, 1))
    Y[3::2] = [10, 100, 10, 1]
    Y[-1] = [0, 1000, 0, 1]
    meta = pd.DataFrame(dict(cell_id=[f'c{i}' for i in range(len(Y))]))
    (tmp_path/'bins.tsv').write_text('bin\n0\n1\n2\n3\n')
    result = confirm_specimen('A', np.array([0, 1]), np.arange(2, len(Y)),
                              meta, Y, res, [], tmp_path)
    scored = result['scored']
    expected = np.full(4097, '', dtype=object)
    expected[1::2] = 'tied_candidate_scores'
    expected[-1] = 'no_proposal_depth'
    np.testing.assert_array_equal(scored['unresolved_reason'], expected)
    np.testing.assert_array_equal(scored['candidate'], np.where(expected == '', 0, -1))
