"""Bulk mode on synthetic bulk samples: arm gains and losses at 40% tumour fraction."""
import numpy as np
import pandas as pd

from haplodose.bulk import run_bulk


def test_bulk_recovers_arm_events(tmp_path):
    rng = np.random.default_rng(7)
    n_chrom, per = 22, 300
    chrom = np.repeat([f'chr{c}' for c in range(1, n_chrom + 1)], per)
    within = np.tile(np.arange(per), n_chrom)
    genes = [f'G{i:05d}' for i in range(len(chrom))]
    # Gene starts past every hg38 p-arm end (> 125 Mb) so each chromosome is one q arm.
    pos = pd.DataFrame({'gene': genes, 'chromosome': chrom, 'start': 130_000_000 + within * 100_000,
                        'end': 130_000_000 + within * 100_000 + 50_000})
    base = rng.lognormal(3, 1.5, len(genes))
    cols, roles = {}, []
    for i in range(20):
        cols[f'N{i}'] = rng.poisson(base * rng.lognormal(0, 0.15, len(genes)) * 20)
        roles.append('normal')
    dose = np.ones(len(genes))
    dose[chrom == 'chr7'] = 1.5
    dose[chrom == 'chr18'] = 0.5
    p, b = 0.4, 0.8
    for i in range(5):
        rate = base * rng.lognormal(0, 0.25, len(genes)) * ((1 - p) + p * dose ** b)
        cols[f'T{i}'] = rng.poisson(rate * 20)
        roles.append('tumour')
    counts = pd.DataFrame(cols, index=genes)
    counts.index.name = 'gene'
    counts.to_csv(tmp_path / 'counts.tsv', sep='\t')
    pd.DataFrame({'sample': list(cols), 'role': roles}).to_csv(tmp_path / 'samples.tsv', sep='\t', index=False)
    pos.to_csv(tmp_path / 'positions.tsv', sep='\t', index=False)
    out = run_bulk(tmp_path / 'counts.tsv', tmp_path / 'samples.tsv', tmp_path / 'positions.tsv', tmp_path / 'out')
    arms = pd.read_csv(out / 'arms.tsv', sep='\t')
    t = arms[arms['sample'].str.startswith('T')].set_index(['sample', 'arm'])
    assert t.xs('chr7q', level='arm').call.eq('gain').all()
    assert t.xs('chr18q', level='arm').call.eq('loss').all()
    assert (t.xs('chr7q', level='arm').gain_fraction > 0.8).all()     # HMM segments agree
    other = t[~t.index.get_level_values('arm').isin(['chr7q', 'chr18q'])]
    assert other.call.ne('neutral').mean() < 0.05
    n = arms[arms['sample'].str.startswith('N')]
    assert n.call.ne('neutral').mean() < 0.05
