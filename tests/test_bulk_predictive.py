"""Opt-in descriptive profiles leave every existing inferential output intact."""
import json

import numpy as np
import pandas as pd
import pytest

from haplodose.bulk import run_bulk
from haplodose.smoothing import predictive_rna_profile


@pytest.mark.parametrize('smoothing', ['predictive', 'native-regions'])
def test_predictive_export_preserves_native_inference_and_literal_axes(tmp_path, smoothing):
    rng = np.random.default_rng(20261109)
    genes = ['G'+str(i) for i in rng.permutation(120)]
    samples = ['T', '001', 'N2', 'N1']
    counts = pd.DataFrame(rng.poisson(30, (120, 4)), index=genes, columns=samples)
    counts.to_csv(tmp_path/'counts.tsv', sep='\t')
    pd.DataFrame(dict(sample=samples, role=['tumour', 'normal', 'normal', 'normal'])).to_csv(
        tmp_path/'samples.tsv', sep='\t', index=False)
    # A long gap divides the same arm into two smoothing blocks.
    starts = np.r_[np.arange(60)*1000, 10_000_000+np.arange(60)*1000]+130_000_000
    pd.DataFrame(dict(gene=['G'+str(i) for i in range(120)], gene_symbol=['G'+str(i) for i in range(120)], chromosome='chr1',
                      start=starts, end=starts+500)).to_csv(tmp_path/'positions.tsv', sep='\t', index=False)
    args = [tmp_path/'counts.tsv', tmp_path/'samples.tsv', tmp_path/'positions.tsv']
    options = dict(bin_genes=4, min_cpm=0, n_programmes=0, max_gap_bp=100_000,
                   counts_kind='raw', genome_build='hg38', coordinate_system='zero-based-half-open',
                   gene_namespace='custom')
    with pytest.warns(UserWarning, match='No panel size establishes calibration'):
        before = run_bulk(*args, tmp_path/'before', **options)
    with pytest.warns(UserWarning, match='No panel size establishes calibration'):
        after = run_bulk(*args, tmp_path/'after', rna_profile_smoothing=smoothing, **options)
    for path in before.iterdir():
        if path.suffix == '.npy':
            np.testing.assert_array_equal(np.load(path), np.load(after/path.name))
        if path.name.endswith(('.tsv', '.tsv.gz')):
            pd.testing.assert_frame_equal(pd.read_csv(path, sep='\t'),
                                          pd.read_csv(after/path.name, sep='\t'), check_exact=True)
    mapping = pd.read_csv(after/'genes.tsv.gz', sep='\t')
    bins = pd.read_csv(after/'bins.tsv', sep='\t')
    from haplodose.bins import chromosome_starts, arm_starts
    blocks = np.cumsum(chromosome_starts(bins) | arm_starts(bins, genome_build='hg38'))-1
    assert len(np.unique(blocks)) == 2
    if smoothing == 'predictive':
        observed = np.load(after/'rna_log_ratio_predictive.npy')
        for row, values in enumerate(np.load(after/'gene_rna_log_ratio.npy')):
            expected, _ = predictive_rna_profile(values, mapping.bin.to_numpy(), blocks)
            np.testing.assert_array_equal(observed[row], expected)
    else:
        import statistics
        observed = np.load(after/'rna_log_ratio_native_regions.npy')
        gene_effect = np.load(after/'gene_rna_log_ratio_programme.npy')
        gene_expectation = np.load(after/'gene_fitted_expected.npy')
        # Compute directly from exported input identities, expectations and
        # scalar medians; the production regional helper is not the oracle.
        mapped_counts = counts.loc[mapping.gene].to_numpy(float).T
        uncentred = np.log((mapped_counts+0.5)/(gene_expectation+0.5))
        for row, path in enumerate(np.load(after/'rna_state_index.npy')):
            centre = statistics.median([statistics.median(uncentred[row, mapping.bin == b])
                                        for b in range(len(bins))])
            np.testing.assert_array_equal(gene_effect[row], uncentred[row]-centre)
            boundaries = np.r_[0, np.flatnonzero((path[1:] != path[:-1]) | (blocks[1:] != blocks[:-1]))+1,
                               len(path)]
            for a, b in zip(boundaries[:-1], boundaries[1:]):
                expected = statistics.median(gene_effect[row, mapping.bin.between(a, b-1)])
                np.testing.assert_array_equal(observed[row, a:b], np.full(b-a, expected))
    metadata = json.loads((after/'run.json').read_text())['rna_profile_smoothing']
    assert metadata['axes'] == ['samples.tsv', 'bins.tsv']
    assert metadata['units'] == 'natural_log_relative_rna_rate'
    assert metadata['calibration_status'] == 'unvalidated'
    assert pd.read_csv(after/'samples.tsv', sep='\t', dtype=str)['sample'].tolist() == samples
