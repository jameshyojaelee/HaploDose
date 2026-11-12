"""Bulk mode on synthetic bulk samples: arm gains and losses at 40% tumour fraction."""
import json

import numpy as np
import pandas as pd
import pytest

from haplodose.bulk import run_bulk, _read_bulk_counts

CONTRACT = dict(counts_kind='raw', genome_build='hg38', coordinate_system='zero-based-half-open',
                gene_namespace='custom')


def _write(tmp_path, n_normals=20, n_chrom=22, per=300):
    rng = np.random.default_rng(7)
    chrom = np.repeat([f'chr{c}' for c in range(1, n_chrom + 1)], per)
    within = np.tile(np.arange(per), n_chrom)
    genes = [f'G{i:05d}' for i in range(len(chrom))]
    # Gene starts past every hg38 p-arm end (> 125 Mb) so each chromosome is one q arm.
    pos = pd.DataFrame({'gene': genes, 'gene_symbol': genes, 'chromosome': chrom, 'start': 130_000_000 + within * 100_000,
                        'end': 130_000_000 + within * 100_000 + 50_000})
    base = rng.lognormal(3, 1.5, len(genes))
    cols, roles = {}, []
    for i in range(n_normals):
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
    return [tmp_path / 'counts.tsv', tmp_path / 'samples.tsv', tmp_path / 'positions.tsv']


@pytest.mark.parametrize('bin_stat', ['sum', 'median'])
def test_bulk_recovers_arm_events(tmp_path, bin_stat):
    """Median mode rebuilds its pseudo counts without each held-out normal, so normals keep a clean null."""
    out = run_bulk(*_write(tmp_path), tmp_path / 'out', bin_stat=bin_stat, **CONTRACT)
    arms = pd.read_csv(out / 'arms.tsv', sep='\t')
    t = arms[arms['sample'].str.startswith('T')].set_index(['sample', 'arm'])
    assert t.xs('chr7q', level='arm').call.eq('gain').all()
    assert t.xs('chr18q', level='arm').call.eq('loss').all()
    assert (t.xs('chr7q', level='arm').gain_fraction > 0.8).all()     # HMM segments agree
    other = t[~t.index.get_level_values('arm').isin(['chr7q', 'chr18q'])]
    assert other.call.ne('neutral').mean() < 0.05
    n = arms[arms['sample'].str.startswith('N')]
    assert n.call.ne('neutral').mean() < 0.05
    run_json = json.loads((out / 'run.json').read_text())
    assert run_json['normal_panel_size'] == 20
    assert 'normal_panel_size_adequate' not in run_json
    assert 'arm_calls_calibrated' not in run_json
    assert run_json['calibration_status'] == 'unvalidated'
    assert run_json['calibration_method'] == 'none'
    assert run_json['output_schema_version'] == 2
    assert 0.0 <= run_json['normal_arm_call_rate_in_sample'] < 0.05
    units = 'median_derived_pseudocounts' if bin_stat == 'median' else 'input_count_units'
    assert run_json['binned_count_units'] == units
    segments = pd.read_csv(out / 'segments.tsv', sep='\t')
    assert segments['count_units'].eq(units).all()
    samples = pd.read_csv(out / 'samples.tsv', sep='\t').set_index('sample')
    assert samples['binned_count_units'].eq(units).all()
    assert samples['input_binned_count_units'].eq('input_count_units').all()
    used = pd.read_csv(out / 'genes.tsv.gz', sep='\t')['gene']
    counts = pd.read_csv(tmp_path / 'counts.tsv', sep='\t', index_col=0)
    expected = counts.loc[used, samples.index].sum(axis=0)
    np.testing.assert_array_equal(samples.input_binned_counts, expected)
    if bin_stat == 'sum':
        np.testing.assert_array_equal(samples.binned_counts, expected)
    else:
        assert not np.array_equal(samples.binned_counts, expected)


def test_small_normal_panel_is_flagged(tmp_path):
    with pytest.warns(UserWarning, match='No panel size establishes calibration'):
        out = run_bulk(*_write(tmp_path, n_normals=3, n_chrom=4, per=200), tmp_path / 'out', **CONTRACT)
    metadata = json.loads((out / 'run.json').read_text())
    assert metadata['normal_panel_size'] == 3
    assert 'normal_panel_size_adequate' not in metadata
    assert metadata['calibration_status'] == 'unvalidated'


@pytest.mark.parametrize('header', ['N0\tN0\tN1', 'N0\t\tN1', 'N0\t \tN1'])
def test_bulk_count_headers_reject_ambiguous_identity(tmp_path, header):
    path = tmp_path / 'counts.tsv'
    path.write_text('gene\t' + header + '\nG0\t3\t7\t5\n')
    with pytest.raises(ValueError, match='nonempty, unique sample IDs'):
        _read_bulk_counts(path)


@pytest.mark.parametrize('compressed', [False, True])
def test_bulk_count_headers_preserve_distinct_literal_ids(tmp_path, compressed):
    expected = pd.DataFrame([[3, 7, 5]], index=['G0'], columns=['N0', 'N0.1', '001'])
    expected.index.name = 'gene'
    path = tmp_path / ('counts.tsv.gz' if compressed else 'counts.tsv')
    expected.to_csv(path, sep='\t')
    pd.testing.assert_frame_equal(_read_bulk_counts(path), expected, check_exact=True)
