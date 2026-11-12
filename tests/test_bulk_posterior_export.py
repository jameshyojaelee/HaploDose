"""The optional boundary-averaged export preserves native count inference."""
import json

import numpy as np
import pandas as pd
import pytest

from haplodose import bulk, cli
from haplodose.bins import chromosome_starts, arm_starts
from haplodose.posterior import posterior_rna_profile


CONTRACT = dict(counts_kind='raw', genome_build='hg38',
                coordinate_system='zero-based-half-open', gene_namespace='custom')


def test_bulk_posterior_export_preserves_native_outputs_and_identified_gene_inputs(tmp_path):
    rng = np.random.default_rng(20261110)
    genes = ['G'+str(i) for i in rng.permutation(120)]
    samples = ['001', 'N2', 'N1', 'N3']
    counts = pd.DataFrame(rng.poisson(30, (120, 4)), index=genes, columns=samples)
    counts.to_csv(tmp_path/'counts.tsv', sep='\t')
    meta = pd.DataFrame(dict(sample=samples, role=['tumour', 'normal', 'normal', 'normal']))
    meta.to_csv(tmp_path/'samples.tsv', sep='\t', index=False)
    starts = np.r_[np.arange(60)*1000, 10_000_000+np.arange(60)*1000]+130_000_000
    pd.DataFrame(dict(gene=['G'+str(i) for i in range(120)], gene_symbol=['G'+str(i) for i in range(120)], chromosome='chr1',
        start=starts, end=starts+500)).to_csv(tmp_path/'positions.tsv', sep='\t', index=False)
    args = [tmp_path/'counts.tsv', tmp_path/'samples.tsv', tmp_path/'positions.tsv']
    options = dict(bin_genes=4, min_cpm=0, n_programmes=0, max_gap_bp=100_000, **CONTRACT)
    with pytest.warns(UserWarning, match='No panel size establishes calibration'):
        before = bulk.run_bulk(*args, tmp_path/'before', **options)
    with pytest.warns(UserWarning, match='No panel size establishes calibration'):
        after = bulk.run_bulk(*args, tmp_path/'after', rna_profile_smoothing='posterior-regions', **options)
    for path in before.iterdir():
        if path.suffix == '.npy':
            np.testing.assert_array_equal(np.load(path), np.load(after/path.name))
        elif path.name.endswith(('.tsv', '.tsv.gz')):
            pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'),
                pd.read_csv(after/path.name, sep='\t'), check_exact=True)
    mapping = pd.read_csv(after/'genes.tsv.gz', sep='\t')
    bins = pd.read_csv(after/'bins.tsv', sep='\t')
    blocks = np.cumsum(chromosome_starts(bins) | arm_starts(bins, genome_build='hg38'))-1
    assert len(np.unique(blocks)) == 2
    effect = np.load(after/'gene_rna_log_ratio_programme.npy')
    expected = np.load(after/'gene_fitted_expected.npy')
    observed = np.load(after/'rna_log_ratio_posterior_regions.npy')
    assert observed.shape == (4, len(bins))
    assert np.isfinite(observed[0]).all() and np.isnan(observed[1:]).all()
    # Independently reconstruct the declared gene effects from source counts
    # and fitted exposures; ordered gene identity is consequential here.
    mapped = counts.loc[mapping.gene].to_numpy(float).T
    raw = np.log((mapped+.5)/(expected+.5))
    for i in range(4):
        centre = np.median([np.median(raw[i, mapping.bin == b]) for b in range(len(bins))])
        np.testing.assert_array_equal(effect[i], raw[i]-centre)
    reference = json.loads((after/'run.json').read_text())['rna_profile_smoothing']['reference']
    independent, _ = posterior_rna_profile(effect[0], expected[0], mapping.bin.to_numpy(),
        effect[1:], expected[1:], blocks, reference_id=reference['identity'], counts_kind='raw')
    np.testing.assert_array_equal(observed[0], independent)
    assert reference['samples'] == samples[1:]
    assert pd.read_csv(after/'samples.tsv', sep='\t', dtype=str)['sample'].tolist() == samples
    decisions = json.loads((after/'rna_posterior_regions.json').read_text())
    assert decisions[0]['sample'] == '001' and decisions[0]['calibration_status'] == 'unvalidated'
    assert [r['status'] for r in decisions[1:]] == ['reference_only']*3
    assert json.loads((after/'run_manifest.json').read_text())['status'] == 'complete'


def test_bulk_posterior_rejects_median_pseudocounts_before_creation(tmp_path):
    out = tmp_path/'out'
    with pytest.raises(ValueError, match="requires bin_stat='sum'"):
        bulk.run_bulk('unused', 'unused', 'unused', out, bin_stat='median',
            rna_profile_smoothing='posterior-regions', **CONTRACT)
    assert not out.exists()


def test_cli_forwards_explicit_posterior_export_and_count_units(monkeypatch):
    calls = []
    def capture(*args, **kwargs):
        calls.append((args, kwargs))
        return 'fresh-output'
    monkeypatch.setattr(bulk, 'run_bulk', capture)
    cli.main(['bulk', '--counts', 'counts.tsv', '--samples', 'samples.tsv',
        '--positions', 'positions.tsv', '--out', 'fresh-output',
        '--counts-kind', 'raw', '--genome-build', 'hg38',
        '--coordinate-system', 'zero-based-half-open', '--gene-namespace', 'custom',
        '--rna-profile-smoothing', 'posterior-regions'])
    assert len(calls) == 1
    assert calls[0][1]['rna_profile_smoothing'] == 'posterior-regions'
    assert calls[0][1]['bin_stat'] == 'sum' and calls[0][1]['counts_kind'] == 'raw'
