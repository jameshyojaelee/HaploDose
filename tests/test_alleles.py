"""Allelic imbalance with phase switch errors; cnLOH clones found only through alleles."""
import numpy as np
import pandas as pd
import pytest

from haplodose.alleles import AlleleParams, allele_transitions, clone_imbalance, imbalance_runs, switch_probs
from haplodose.pipeline import run
from haplodose.simulate import write_dataset


def _snps(n, chrom='chr1', step_cm=0.01):
    return pd.DataFrame({'chromosome': chrom, 'pos': np.arange(n) * 10_000 + 1, 'cM': np.arange(n) * step_cm,
                         'gene': ''})


def test_transitions_are_stochastic():
    p, _ = switch_probs(_snps(50), 1.0)
    T = np.exp(np.asarray(allele_transitions(p, 1e-3, 4)))
    assert np.allclose(T[1:].sum(axis=2), 1.0, atol=1e-5)


def test_imbalance_survives_phase_switches():
    rng = np.random.default_rng(0)
    n = 400
    snps = _snps(n)
    flip = np.cumsum(rng.random(n) < 0.1) % 2 == 1          # input phase flips about every 10 SNPs
    D = rng.poisson(20, n).astype(float)
    b_true = np.where(np.arange(n) < 200, 0.5, 0.8)          # second half: (2,0)-like imbalance
    a_true = rng.binomial(D.astype(int), b_true)
    A = np.where(flip, D - a_true, a_true).astype(float)
    second = slice(200, None)
    assert abs(A[second].sum() / D[second].sum() - 0.5) < 0.1   # naive pooling cancels the imbalance
    post, B, nu, _ = clone_imbalance(A[None], D[None], snps, 200.0, AlleleParams())
    runs = [r for r in imbalance_runs(post[0], snps, D, AlleleParams()) if r['resolved']]
    assert len(runs) == 1 and runs[0]['first_snp'] >= 190 and runs[0]['last_snp'] == n - 1
    assert (1 - post[0, :190, 0]).mean() < 0.1


@pytest.fixture(scope='module')
def cnloh(tmp_path_factory):
    root = tmp_path_factory.mktemp('al')
    spec = [dict(name='A', n_cells=600, normal_frac=0.15, clones=[
        dict(frac=0.5, events=[('chr3', 0.0, 1.0, 3), ('chr8', 0.0, 1.0, 1)]),
        dict(frac=0.5, events=[('chr3', 0.0, 1.0, 3), ('chr8', 0.0, 1.0, 1), ('chr12', 0.0, 1.0, 2, 'loh')])])]
    data = write_dataset(root / 'data', specimens=spec, alleles=True, ref_patient='A', ref_cells_per_type=400)
    run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', root / 'with', alleles=data / 'alleles.tsv.gz')
    run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', root / 'without')
    return root


def test_cnloh_clone_separated_only_with_alleles(cnloh):
    def tumour(out):
        cells = pd.read_csv(cnloh / out / 'cells.tsv.gz', sep='\t')
        return cells[cells.truth_clone.isin(['clone1', 'clone2']) & cells.malignant_call]
    t = tumour('with')
    assert t.clone.nunique() == 2
    table = pd.crosstab(t.truth_clone, t.clone)
    acc = max(table.values.diagonal().sum(), table.values[::-1].diagonal().sum()) / table.values.sum()
    # chr12 carries ~3 SNP reads per cell; a balanced cell whose reads all match the LOH
    # orientation looks LOH (~20% of cells), so the Bayes-optimal accuracy is ~0.89.
    assert acc > 0.75
    assert tumour('without').clone.nunique() == 1


def test_joint_states(cnloh):
    seg = pd.read_csv(cnloh / 'with' / 'segments.tsv', sep='\t')
    tum = seg[~seg.clone.str.endswith('_c0') & seg.resolved]
    assert set(tum.loc[tum.chromosome.eq('chr8'), 'joint_state']) == {'loss'}
    allele = pd.read_csv(cnloh / 'with' / 'allele_segments.tsv', sep='\t')
    chr12 = allele[allele.resolved & allele.chromosome.eq('chr12')]
    assert chr12.clone.nunique() == 1                         # only the cnLOH clone
    no = pd.read_csv(cnloh / 'without' / 'segments.tsv', sep='\t')
    assert 'joint_state' not in no.columns


def test_normals_not_called_by_alleles(cnloh):
    ref = pd.read_csv(cnloh / 'with' / 'reference_check.tsv', sep='\t')
    assert (ref.fraction_called_aneuploid <= ref.expected_fraction + 0.05).all()
