"""Fine genomic bins of consecutive expressed genes.

Bins are the unit for per-cell sufficient statistics (UMI sums). Segments are
unions of bins, so one pass over the counts serves every later iteration.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .io import AUTOSOMES

# Ribosomal protein genes (and their pseudogenes): the most highly expressed genes,
# they dominate their bins' counts and shift together with cell state and quality,
# which made whole-chromosome (chr19) clone splits on SMC and Moorman. RPS6K*
# kinases do not match.
RIBOSOMAL = r'^(RPL\d+[A-Z]*\d*|RPS\d+[A-Z]?\d*|RPLP\d|RPSA(P\d+)?)$'


def chromosome_rank(chrom: pd.Series) -> pd.Series:
    order = {c: i for i, c in enumerate(AUTOSOMES + ['chrX'])}
    return chrom.map(order)


def make_bins(genes, positions: pd.DataFrame, gene_mean, exclude=frozenset(),
              bin_genes=10, min_mean=0.02, chromosomes=tuple(AUTOSOMES), exclude_pattern=RIBOSOMAL):
    """Returns (bins table, gene_bin array with -1 for unused genes, used-gene table).

    Genes named in `exclude` or matching `exclude_pattern` (None: none) are never features."""
    genes = pd.Index(np.asarray(genes, dtype=str))
    pos = positions.set_index(positions.gene.str.upper())
    table = pd.DataFrame({'gene': genes, 'column': np.arange(len(genes)),
                          'mean': np.asarray(gene_mean, dtype=float)})
    key = table.gene.str.upper()
    table['chromosome'] = key.map(pos.chromosome)
    table['start'] = key.map(pos.start)
    table['end'] = key.map(pos.end)
    keep = (table.chromosome.isin(chromosomes) & table.start.notna()
            & (table['mean'] >= min_mean) & ~key.isin(exclude))
    if exclude_pattern is not None:
        keep &= ~key.str.match(exclude_pattern)
    table = table[keep].copy()
    table['rank'] = chromosome_rank(table.chromosome)
    table = table.sort_values(['rank', 'start', 'column'], kind='stable')
    gene_bin = np.full(len(genes), -1, dtype=np.int32)
    rows, used = [], []
    for chrom, frame in table.groupby('rank', sort=True):
        n = len(frame)
        if n < max(3, bin_genes // 2):
            continue
        n_bins = max(1, n // bin_genes)  # remainder spread over bins, none smaller than bin_genes
        edges = np.linspace(0, n, n_bins + 1).round().astype(int)
        for a, b in zip(edges[:-1], edges[1:]):
            block = frame.iloc[a:b]
            k = len(rows)
            gene_bin[block.column.to_numpy()] = k
            rows.append(dict(bin=k, chromosome=block.chromosome.iloc[0],
                             start=int(block.start.min()), end=int(block.end.max()),
                             n_genes=len(block), first_gene=block.gene.iloc[0],
                             last_gene=block.gene.iloc[-1]))
            used.append(block.assign(bin=k))
    if not rows:
        raise ValueError('No genomic bins: check gene names, positions and expression filter')
    bins = pd.DataFrame(rows)
    used = pd.concat(used)[['gene', 'column', 'chromosome', 'start', 'end', 'mean', 'bin']]
    return bins, gene_bin, used


def chromosome_starts(bins: pd.DataFrame) -> np.ndarray:
    """Boolean per bin: True where an HMM chain restarts (first bin of a chromosome)."""
    chrom = bins.chromosome.to_numpy()
    reset = np.ones(len(chrom), dtype=bool)
    reset[1:] = chrom[1:] != chrom[:-1]
    return reset


def arm_labels(bins: pd.DataFrame, centromeres=None) -> np.ndarray:
    """'chrNp' / 'chrNq' per bin (bin start before / after the p-arm end; default hg38 table)."""
    if centromeres is None:
        centromeres = Path(__file__).with_name('data') / 'centromeres_hg38.tsv'
    cen = pd.read_csv(centromeres, sep='\t', comment='#').set_index('chromosome').p_arm_end
    q = bins.start.to_numpy() >= bins.chromosome.map(cen).fillna(np.inf).to_numpy()
    return (bins.chromosome + np.where(q, 'q', 'p')).to_numpy()


def arm_starts(bins: pd.DataFrame, centromeres=None) -> np.ndarray:
    """chromosome_starts plus the first q-arm bin of each chromosome (bin start past the p-arm end).

    centromeres: table with chromosome, p_arm_end; default is the bundled hg38 table.
    Chromosomes missing from the table get no arm boundary."""
    if centromeres is None:
        centromeres = Path(__file__).with_name('data') / 'centromeres_hg38.tsv'
    cen = pd.read_csv(centromeres, sep='\t', comment='#').set_index('chromosome').p_arm_end
    reset = chromosome_starts(bins)
    q = bins.start.to_numpy() >= bins.chromosome.map(cen).fillna(np.inf).to_numpy()
    return reset | (q & ~np.r_[False, q[:-1]])
