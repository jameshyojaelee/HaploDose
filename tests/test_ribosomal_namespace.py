"""Category filtering requires biological symbols independently of count IDs."""
import numpy as np
import pandas as pd
import pytest

from haplodose.bins import make_bins


def positions(namespace):
    symbols = ['RPL13A', 'RPS19', 'RPLP0', 'RPSA', 'RPS6KA1', 'TP53', 'ACTB', 'GAPDH']
    genes = ([f'ENSG{i:011d}' for i in range(8)] if namespace == 'ensembl'
             else [f'id{i}' for i in range(8)])
    table = pd.DataFrame(dict(gene=genes, gene_symbol=symbols, chromosome='chr19',
        start=np.arange(8)*1000, end=np.arange(8)*1000+500))
    table.attrs['gene_namespace'] = namespace
    return table


@pytest.mark.parametrize('namespace', ['ensembl', 'custom'])
def test_default_exclusion_tracks_biological_symbols(namespace):
    pos = positions(namespace)
    _, mapping, used = make_bins(pos.gene, pos, np.ones(8), bin_genes=2)
    # Expect identical biological retention to the literal symbol namespace.
    assert set(used.gene) == set(pos.gene.iloc[4:])
    assert set(used.gene_symbol) == set(pos.gene_symbol.iloc[4:])
    assert np.array_equal(mapping < 0, np.arange(8) < 4)
    symbolic = pos.assign(gene=pos.gene_symbol).copy()
    symbolic.attrs['gene_namespace'] = 'symbol'
    _, other, kept = make_bins(symbolic.gene, symbolic, np.ones(8), bin_genes=2)
    np.testing.assert_array_equal(mapping, other)
    assert used.gene_symbol.tolist() == kept.gene.tolist()


@pytest.mark.parametrize('namespace', ['ensembl', 'custom'])
def test_missing_annotation_rejects_but_explicit_no_category_filter_is_allowed(namespace):
    pos = positions(namespace).drop(columns='gene_symbol')
    with pytest.raises(ValueError, match='gene_symbol'):
        make_bins(pos.gene, pos, np.ones(8), bin_genes=2)
    _, _, used = make_bins(pos.gene, pos, np.ones(8), bin_genes=2, exclude_pattern=None)
    assert set(used.gene) == set(pos.gene)


@pytest.mark.parametrize('bad', ['', None, 'NA', 'nan', 12])
def test_retained_features_require_unambiguous_symbol_annotation(bad):
    pos = positions('ensembl')
    pos['gene_symbol'] = pos.gene_symbol.astype(object)
    pos.loc[0, 'gene_symbol'] = bad
    with pytest.raises(ValueError, match='gene_symbol'):
        make_bins(pos.gene, pos, np.ones(8), bin_genes=2)


def test_category_mapping_does_not_guess_from_custom_ids_or_validate_unused_rows():
    pos = positions('custom')
    pos.loc[4, 'gene'] = 'RPL13A'  # ID resembles a ribosomal symbol, annotation does not.
    pos.loc[5, 'gene_symbol'] = None
    _, _, used = make_bins(pos.gene, pos, np.ones(8), bin_genes=2, exclude={pos.gene[5]})
    assert 'RPL13A' in set(used.gene)
    assert pos.gene[0] not in set(used.gene)
