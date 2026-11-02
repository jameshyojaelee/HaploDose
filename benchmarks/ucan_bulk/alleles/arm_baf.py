"""Arm-level allelic imbalance from bulk RNA allele counts at common SNPs, without germline genotypes.

Per sample, each site is a mixture over the unknown germline genotype with Hardy-Weinberg priors from the
1000G EUR allele frequency f: hom-ref (1-f)^2 and hom-alt f^2 (alt reads ~ Binomial(n, eps) or (n, 1-eps)),
het 2f(1-f) (alt reads ~ beta-binomial with mean b or 1-b, phase unknown, concentration c). b is the arm's
minor-allele fraction (0.5 = balanced), fitted per arm on a grid. c (allele-specific expression and
other overdispersion) is one value per sample, chosen by profile likelihood over the grid. An arm with
LOH in a pure tumour has b near 0: het sites look homozygous, and the HWE prior is what tells them apart.
Output per sample and arm: b_hat, LLR of b_hat against b = 0.5, sites used, expected het sites.
Usage: arm_baf.py PILEUP_DIR SITES_DIR OUT_TSV [MIN_DEPTH]
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import betabinom, binom

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'ucan_spatial'))
from dna_truth import CEN, arm_of  # noqa: E402

B_GRID = np.round(np.arange(0.01, 0.505, 0.01), 3)
C_GRID = (10.0, 20.0, 40.0, 80.0, 160.0, 1e4)
EPS = 0.005
COLS = ['chrom', 'pos', 'ref', 'alt', 'ref_count', 'alt_count', 'depth']


def site_loglik(k, n, f, c):
    """log p(k | n) per site [S] for every b in B_GRID: [len(B_GRID), S]."""
    hom = np.logaddexp(2 * np.log1p(-f) + binom.logpmf(k, n, EPS), 2 * np.log(f) + binom.logpmf(k, n, 1 - EPS))
    b = B_GRID[:, None]
    het = np.logaddexp(betabinom.logpmf(k, n, b * c, (1 - b) * c), betabinom.logpmf(k, n, (1 - b) * c, b * c)) - np.log(2)
    return np.logaddexp(hom[None], np.log(2 * f * (1 - f)) + het), hom


def fit_sample(t):
    arms = t.arm.to_numpy()
    codes, names = pd.factorize(arms)
    best = None
    for c in C_GRID:
        ll, hom = site_loglik(t.alt_count.to_numpy(), t.n.to_numpy(), t.af.to_numpy(), c)
        per_arm = np.stack([np.bincount(codes, weights=row, minlength=len(names)) for row in ll])  # [B, arms]
        total = per_arm.max(axis=0).sum()
        if best is None or total > best[0]:
            best = (total, c, ll, hom, per_arm)
    _, c, ll, hom, per_arm = best
    i = per_arm.argmax(axis=0)
    b_hat = B_GRID[i]
    # Posterior het probability per site at the arm's b_hat, summed per arm.
    lhet = ll[i[codes], np.arange(len(codes))]
    p_het = 1 - np.exp(hom - lhet)        # hom already carries its genotype priors
    out = pd.DataFrame(dict(arm=names, b_hat=b_hat, llr=per_arm[i, np.arange(len(names))] - per_arm[-1],
                            n_sites=np.bincount(codes), het_sites=np.bincount(codes, weights=np.clip(p_het, 0, 1))))
    return out.assign(concentration=c)


def one_sample(f, af, cen, min_depth):
    sample = f.name.removesuffix('.alleles.tsv.gz')
    t = pd.read_csv(f, sep='\t', header=None, names=COLS, dtype={'chrom': str}).merge(af, on=['chrom', 'pos'])
    t['n'] = t.ref_count + t.alt_count
    t = t[(t.n >= min_depth) & t.af.between(0.01, 0.99)]
    t['chrom'] = 'chr' + t.chrom
    t = t[t.chrom.isin(cen.index)]
    t['arm'] = arm_of(t.chrom, t.pos.to_numpy(), cen)
    print(sample, len(t), flush=True)
    return fit_sample(t).assign(sample=sample)


def main(pileup_dir, sites_dir, out_tsv, min_depth=8):
    min_depth = int(min_depth)
    assert B_GRID[-1] == 0.5
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    af = pd.read_csv(Path(sites_dir) / 'sites.tsv', sep='\t', header=None, usecols=[0, 1, 5], names=['chrom', 'pos', 'af'],
                     dtype={'chrom': str})
    files = sorted(Path(pileup_dir).glob('*.alleles.tsv.gz'))
    with ProcessPoolExecutor(int(os.environ.get('SLURM_CPUS_PER_TASK', 1))) as pool:
        rows = list(pool.map(partial(one_sample, af=af, cen=cen, min_depth=min_depth), files))
    res = pd.concat(rows)[['sample', 'arm', 'b_hat', 'llr', 'n_sites', 'het_sites', 'concentration']]
    res.to_csv(out_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main(*sys.argv[1:])
