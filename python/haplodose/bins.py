"""Fine genomic bins of consecutive expressed genes.

Bins are the unit for per-cell sufficient statistics (UMI sums). Segments are
unions of bins, so one pass over the counts serves every later iteration.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .io import AUTOSOMES, validate_feature_ids

# Ribosomal protein genes (and their pseudogenes): the most highly expressed genes,
# they dominate their bins' counts and shift together with cell state and quality,
# which made whole-chromosome (chr19) clone splits on SMC and Moorman. RPS6K*
# kinases do not match.
RIBOSOMAL = r'^(RPL\d+[A-Z]*\d*|RPS\d+[A-Z]?\d*|RPLP\d|RPSA(P\d+)?)$'


def chromosome_rank(chrom: pd.Series) -> pd.Series:
    order = {c: i for i, c in enumerate(AUTOSOMES + ['chrX'])}
    return chrom.map(order)


def make_bins(genes, positions: pd.DataFrame, gene_mean, exclude=frozenset(),
              bin_genes=10, min_mean=0.02, chromosomes=tuple(AUTOSOMES), exclude_pattern=RIBOSOMAL,
              *, genome_build='hg38', centromeres=None, max_gap_bp=None):
    """Returns (bins table, gene_bin array with -1 for unused genes, used-gene table).

    Genes named in `exclude` or matching `exclude_pattern` (None: none) are never features."""
    if not isinstance(bin_genes, (int, np.integer)) or bin_genes <= 0:
        raise ValueError('bin_genes must be a positive integer')
    if not np.isfinite(min_mean) or min_mean < 0:
        raise ValueError('min_mean must be finite and nonnegative')
    if max_gap_bp is not None and (not isinstance(max_gap_bp, (int, np.integer)) or max_gap_bp < 0):
        raise ValueError('max_gap_bp must be a nonnegative integer or None')
    namespace = positions.attrs.get('gene_namespace', 'symbol')
    genes = pd.Index(validate_feature_ids(genes, namespace))
    pos = positions.copy()
    pos['gene'] = validate_feature_ids(pos.gene, namespace)
    pos = pos.set_index('gene')
    table = pd.DataFrame({'gene': genes, 'column': np.arange(len(genes)),
                          'mean': np.asarray(gene_mean, dtype=float)})
    key = table.gene
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
    cen = read_centromeres(genome_build, centromeres)
    boundaries = table.chromosome.map(cen)
    if boundaries.isna().any():
        raise ValueError('Centromere table lacks a retained chromosome')
    # A gene spanning the boundary cannot support either arm unambiguously.
    table = table[~((table.start < boundaries) & (table.end > boundaries))].copy()
    table['arm'] = arm_labels(table, centromeres, genome_build)
    previous_end = table.groupby('chromosome', sort=False).end.transform(lambda x: x.cummax().shift(1))
    gap = table.start.to_numpy() - previous_end.to_numpy()
    starts = table.arm.ne(table.arm.shift()).to_numpy()
    if max_gap_bp is not None:
        starts |= gap > max_gap_bp
    table['block'] = np.cumsum(starts) - 1
    gene_bin = np.full(len(genes), -1, dtype=np.int32)
    rows, used = [], []
    for _, frame in table.groupby('block', sort=True):
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
                             last_gene=block.gene.iloc[-1], genomic_block=int(block.block.iloc[0])))
            used.append(block.assign(bin=k))
    if not rows:
        raise ValueError('No genomic bins: check gene names, positions and expression filter')
    bins = pd.DataFrame(rows)
    used = pd.concat(used)[['gene', 'column', 'chromosome', 'start', 'end', 'mean', 'bin']]
    return bins, gene_bin, used


def chromosome_starts(bins: pd.DataFrame) -> np.ndarray:
    """HMM restarts at chromosomes and declared arm/gap blocks."""
    chrom = bins.chromosome.to_numpy()
    reset = np.ones(len(chrom), dtype=bool)
    reset[1:] = chrom[1:] != chrom[:-1]
    if 'genomic_block' in bins:
        block = bins.genomic_block.to_numpy()
        reset[1:] |= block[1:] != block[:-1]
    return reset


def read_centromeres(genome_build='hg38', centromeres=None):
    """Centromere boundary in internal zero-based coordinates for the declared build."""
    if centromeres is None:
        if genome_build != 'hg38':
            raise ValueError('Only hg38 has bundled centromeres; supply a centromeres table for this build')
        centromeres = Path(__file__).with_name('data') / 'centromeres_hg38.tsv'
    table = pd.read_csv(centromeres, sep='\t', comment='#', keep_default_na=False)
    if not {'chromosome', 'p_arm_end'} <= set(table):
        raise ValueError('Centromere table requires chromosome and p_arm_end')
    if table.chromosome.duplicated().any() or table.chromosome.astype(str).str.strip().eq('').any():
        raise ValueError('Centromere chromosomes must be nonempty and unique')
    values = pd.to_numeric(table.p_arm_end, errors='coerce').to_numpy(dtype=float)
    if not np.isfinite(values).all() or np.any(values <= 0) or np.any(values != np.floor(values)):
        raise ValueError('Centromere boundaries must be finite positive integers')
    return pd.Series(values, index=table.chromosome)


def arm_labels(bins: pd.DataFrame, centromeres=None, genome_build='hg38') -> np.ndarray:
    """'chrNp' / 'chrNq' from explicit zero-based centromere boundaries."""
    cen = read_centromeres(genome_build, centromeres)
    boundaries = bins.chromosome.map(cen)
    if boundaries.isna().any():
        raise ValueError('Centromere table lacks a retained chromosome')
    q = bins.start.to_numpy() >= boundaries.to_numpy()
    return (bins.chromosome + np.where(q, 'q', 'p')).to_numpy()


def arm_starts(bins: pd.DataFrame, centromeres=None, genome_build='hg38') -> np.ndarray:
    """chromosome_starts plus the first q-arm bin of each chromosome (bin start past the p-arm end).

    centromeres: table with chromosome, p_arm_end; default is the bundled hg38 table.
    Chromosomes missing from the table are rejected."""
    reset = chromosome_starts(bins)
    labels = arm_labels(bins, centromeres, genome_build)
    return reset | np.r_[True, labels[1:] != labels[:-1]]
