import numpy as np
import pandas as pd

from haplodose.bins import make_bins


def test_ribosomal_genes_are_not_features():
    genes = ['RPL13A', 'RPS19', 'RPLP0', 'RPSA', 'RPS6KA1'] + [f'G{i}' for i in range(20)]
    pos = pd.DataFrame({'gene': genes, 'chromosome': 'chr19', 'start': np.arange(25) * 1000 + 1,
                        'end': np.arange(25) * 1000 + 500})
    bins, gene_bin, used = make_bins(genes, pos, np.ones(25), bin_genes=5)
    assert set(used.gene) == {'RPS6KA1'} | {f'G{i}' for i in range(20)}
    _, _, kept = make_bins(genes, pos, np.ones(25), bin_genes=5, exclude_pattern=None)
    assert len(kept) == 25
