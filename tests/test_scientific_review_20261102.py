"""Regression tests for docs/scientific_code_review_2026-11-02.md (verdicts in the _response.md file):
#1 untestable depth-matched pairs, #2 allele veto exemptions and mirrored imbalance, #6 zero usable SNPs,
#7 evidence-level joint-state label."""
import numpy as np
import pandas as pd
import pytest

from haplodose.alleles import (AlleleParams, SpecimenAlleles, balanced_cells, clone_imbalance, joint_state,
                               switch_probs)
from haplodose.clones import CloneParams, Rows, SpecimenModel
from haplodose.profile import ProfileParams, difference_segments


INPUT_CONTRACT = dict(counts_kind='raw', genome_build='hg38',
                      coordinate_system='one-based-closed', gene_namespace='symbol')


def _reset(J, starts=(0,)):
    r = np.zeros(J, bool)
    r[list(starts)] = True
    return r


# 1. Depth-matched pairs with no depth overlap ------------------------------------------------------------
def test_zero_umi_pair_is_untestable_not_nan():
    """Before: a pair with L = 0 gave a segment with mean_diff, z and posterior all NaN (NaN < threshold is
    False, so it was kept). Now it is reported untestable with no segments; other pairs are unaffected."""
    rng = np.random.default_rng(1)
    J = 200
    Ya = rng.poisson(50.0, (2, J)).astype(float)
    Yb = rng.poisson(50.0, (2, J)).astype(float)
    Ya[1, 40:100] *= 2
    La = np.array([0.0, Ya[1].sum()])
    Ya[0] = 0.0
    Lb = Yb.sum(axis=1)
    res = difference_segments(Ya, La, Yb, Lb, _reset(J), 0.8, ProfileParams())
    assert res[0] == dict(segments=[], distance=float('inf'), tau=res[0]['tau'], testable=False)
    assert res[1]['testable'] and np.isfinite(res[1]['distance'])
    assert any(s['direction'] == 'higher' for s in res[1]['segments'])
    for r in res:
        for s in r['segments']:
            assert np.isfinite([s['mean_diff'], s['z'], s['posterior']]).all()


def test_matched_differences_disjoint_depths_untestable():
    rng = np.random.default_rng(0)
    J = 240
    depth = np.r_[np.full(150, 800.0), np.full(150, 20000.0), np.full(150, 820.0)]
    Y = rng.poisson(depth[:, None] / J * np.ones((1, J))).astype(np.float32)
    m = SpecimenModel(Rows(Y, np.arange(len(Y))), Y.sum(axis=1), np.full((J, 1), 1.0 / J), _reset(J, (0, 120)),
                      np.full(J, 0.05), CloneParams(n_programmes=0), gamma=np.zeros(J))
    labels = np.repeat([1, 2, 3], 150)
    disjoint, overlapping = m.matched_differences(labels, [(1, 2), (1, 3)], report_depth_balance=True)
    assert not disjoint['testable'] and disjoint['segments'] == [] and disjoint['distance'] == float('inf')
    assert disjoint['untestable_reason'] == 'no_shared_depth_stratum_mass'
    assert disjoint['weighted_log_depth_wasserstein'] is None
    assert overlapping['testable'] and np.isfinite(overlapping['distance'])
    assert np.isfinite(overlapping['weighted_log_depth_wasserstein'])


# 2. Allele veto and mirrored imbalance --------------------------------------------------------------------
def _alleles(S=100):
    snps = pd.DataFrame({'chromosome': 'chr1', 'pos': np.arange(S) * 1000 + 1, 'cM': np.arange(S) * 0.001,
                         'gene': ''})
    return SpecimenAlleles(np.zeros((2, S)), np.ones((2, S)), snps, 200.0, AlleleParams()), np.full(S, 60.0)


def test_veto_keeps_one_copy_logic_but_exempts_ratio_preserving_changes():
    al, D = _alleles()
    idx = np.arange(len(D))
    bal, two_one, loh = 0.5 * D, 2 / 3 * D, 0.98 * D
    # Unchanged, well-powered imbalance still vetoes a one-copy-sized difference (unchanged behaviour).
    assert al.veto(idx, bal, D, bal, D)
    assert al.veto(idx, two_one, D, two_one, D, log_dose=np.log(1.5))
    # LOH in both groups: 1:0 vs 2:0 leaves b unchanged, so alleles cannot rule a dosage change out.
    assert not al.veto(idx, loh, D, loh, D)
    # A doubling-sized difference with unchanged balance (1:1 vs 2:2) is not vetoed, in either direction.
    assert not al.veto(idx, bal, D, bal, D, log_dose=np.log(2.0))
    assert not al.veto(idx, bal, D, bal, D, log_dose=-np.log(2.0))


def test_differ_sees_mirrored_imbalance():
    al, D = _alleles()
    run = [dict(first_snp=0, last_snp=len(D) - 1, resolved=True)]
    assert al.differ(run, 0.9 * D, D, run, 0.1 * D, D) > 0          # was 0: same |b - 0.5|, opposite haplotypes
    assert al.differ(run, 0.9 * D, D, run, 0.5 * D, D) > 0
    assert al.differ(run, 0.9 * D, D, run, 0.9 * D, D) == 0


def test_differ_null_false_positive_rate():
    """Same beta-binomial allele fraction in both groups (with phase flips shared by both): the signed contrast
    must not add differences (z >= 4) beyond the magnitude test."""
    al, D = _alleles()
    run = [dict(first_snp=0, last_snp=len(D) - 1, resolved=True)]
    rng = np.random.default_rng(3)
    kappa = al.kappa
    hits = 0
    for _ in range(300):
        b = np.where(rng.random(len(D)) < 0.5, 0.8, 0.2)              # input phase switches, same in both groups
        draw = lambda: rng.binomial(D.astype(int), rng.beta(kappa * b, kappa * (1 - b))).astype(float)
        hits += al.differ(run, draw(), D, run, draw(), D) > 0
    assert hits <= 3


# 6. Zero usable SNPs --------------------------------------------------------------------------------------
def test_empty_snp_set_gives_empty_allele_channel():
    """Before: switch_probs indexed cm[0] and raised IndexError when filtering left no SNPs."""
    empty = pd.DataFrame({'chromosome': pd.Series(dtype=str), 'pos': pd.Series(dtype=int),
                          'cM': pd.Series(dtype=float), 'gene': pd.Series(dtype=str)})
    p, reset = switch_probs(empty, 1.0)
    assert p.shape == (0,) and reset.shape == (0,)
    post, B, nu, ll = clone_imbalance(np.zeros((3, 0)), np.zeros((3, 0)), empty, 50.0, AlleleParams())
    assert post.shape[:2] == (3, 0) and B.shape == (3, 0) and np.isnan(nu) and np.all(ll == 0)
    al = SpecimenAlleles(np.zeros((5, 0)), np.zeros((5, 0)), empty, 50.0, AlleleParams())
    Ak, Dk = al.pool(np.ones((5, 2)))
    post, B = al.fit(Ak, Dk)
    assert al.runs(post[0], Dk[0]) == [] and not al.has_events(post[0], Dk[0])
    assert al.differ([], Ak[0], Dk[0], [], Ak[1], Dk[1]) == 0
    assert not al.veto(np.array([], dtype=int), Ak[0], Dk[0], Ak[1], Dk[1])
    assert al.cell_loglik(slice(0, 5), B).shape == (5, 2)
    assert len(balanced_cells(al)[0]) == 0


def test_pipeline_runs_when_specimen_has_no_snps(tmp_path):
    """Allele reads only in reference cells: every SNP is uncovered in the tumour specimen. The run must finish
    expression-only for it and record n_snps = 0."""
    from haplodose.pipeline import run
    from haplodose.simulate import write_dataset
    spec = [dict(name='A', n_cells=300, normal_frac=0.15,
                 clones=[dict(frac=1.0, events=[('chr3', 0.0, 1.0, 3)])])]
    data = write_dataset(tmp_path / 'data', specimens=spec, alleles=True, ref_patient='A', ref_cells_per_type=200)
    cells = pd.read_csv(data / 'cells.tsv', sep='\t')
    al = pd.read_csv(data / 'alleles.tsv.gz', sep='\t')
    ref = set(cells.loc[cells.role.eq('reference'), 'cell_id'])
    al[al.cell.isin(ref)].to_csv(tmp_path / 'ref_only_alleles.tsv.gz', sep='\t', index=False)
    out = run(data / 'counts.h5', data / 'cells.tsv', data / 'positions.tsv', tmp_path / 'out',
              alleles=tmp_path / 'ref_only_alleles.tsv.gz', allele_genome_build='hg38', **INPUT_CONTRACT)
    specs = pd.read_csv(tmp_path / 'out' / 'specimens.tsv', sep='\t')
    assert (specs.n_snps == 0).all()
    assert out is not None


# 7. Evidence-level labels ---------------------------------------------------------------------------------
def test_neutral_imbalanced_is_not_called_cnloh():
    p = AlleleParams()
    assert joint_state(2, 0.99, 50, 500.0, p) == 'neutral_allele_imbalanced'
    assert joint_state(2, 0.01, 50, 500.0, p) == 'neutral'
