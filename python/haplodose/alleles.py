"""Phased heterozygous SNP counts: clone-level allelic imbalance and per-cell evidence.

Input is a Numbat-style allele table (cell, CHROM, POS, cM, AD, DP, GT, gene) with
population-phased genotypes (GT 0|1 or 1|0). Counts are oriented to the input
phase: a = AD if GT is 0|1, else DP - AD. The input phase has switch errors, so
pooling oriented counts over a segment cancels real imbalance; the HMM below
follows the orientation along the genome instead.

Clone level: pooled oriented counts A_ks, D_ks of clone k per SNP s. Hidden
states are balanced (b = 0.5) or imbalanced at level theta_l in either
orientation (b = 0.5 +- theta_l). Orientation flips between neighbouring SNPs with
the phase-switch probability (1 - exp(-2 nu d)) / 2, d in cM (Numbat's form;
nu is chosen by likelihood on the tumour pools). Emissions are beta-binomial
with concentration kappa estimated from same-patient reference cells. This
represents reference overdispersion; context-specific allele expression remains
an alternative to genomic imbalance.

SNPs are dropped when same-patient reference cells show clear imbalance at them
(imprinting, strong eQTL, mapping bias), when their gene or site is excluded, or
when they lie off the autosomes.

Cell level: with b_ks the posterior mean oriented B-allele fraction of clone k,
cell c adds sum_s a_cs log b_ks + (d_cs - a_cs) log(1 - b_ks) to its clone
log-likelihood. Normal-like clone 0 has b = 0.5.

Joint states: expression gives loss/neutral/gain/amplification relative dosage;
alleles add balanced/imbalanced. Relative-neutral dosage with imbalance is labelled
`neutral_allele_imbalanced`: it is consistent with copy-neutral LOH (cnLOH) but does not
establish it (relative dosage cannot fix absolute copy number or ploidy; e.g. 2:1 in a
near-triploid genome reads the same). Expression-only runs never report imbalance.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import jax.numpy as jnp
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.special import betaln

from .io import AUTOSOMES
from .kernels.hmm import forward_backward_batch

LEVELS = np.array([0.08, 0.17, 0.3, 0.45])     # |b - 0.5|: subclonal/ASE-diluted, (2,1) gain, mixed, LOH
NU_GRID = np.array([0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0])
KAPPA_FACTORS = np.array([1.0, 2.5, 5.0, 10.0])


@dataclass
class AlleleParams:
    p_change: float = 5e-4
    min_snps: int = 10
    min_depth: float = 30.0
    posterior: float = 0.95
    bias_alpha: float = 1e-3
    bias_min_depth: int = 10
    cell_shrink: float = 0.9
    default_kappa: float = 50.0
    cell_kappa: float | None = None       # explicit cell-level concentration; None retains legacy binomial


@dataclass
class AlleleData:
    snps: pd.DataFrame          # chromosome, pos (one-based SNP position), cM, gene; sorted along the genome
    A: sparse.csr_matrix        # cells x SNPs, counts on the input-phase haplotype
    D: sparse.csr_matrix        # cells x SNPs, total counts
    excluded: dict = field(default_factory=dict)

    def subset(self, rows):
        return self.A[rows], self.D[rows]


def _chrom(x):
    x = pd.Series(x, dtype=str)
    return np.where(x.str.startswith('chr'), x, 'chr' + x)


def read_alleles(path, cell_ids, exclude_genes=(), exclude_sites=(), chunksize=1_000_000) -> AlleleData:
    """Numbat allele table -> AlleleData over `cell_ids` (matrix order).

    POS and site exclusions 'chrN:pos' retain one-based SNP coordinates. Expression bins and emitted
    allele intervals use zero-based half-open coordinates; a SNP belongs to [start, end) iff start < POS <= end.
    """
    index = pd.Index(np.asarray(cell_ids, dtype=str))
    exclude_genes = {g.upper() for g in exclude_genes}
    exclude_sites = set(exclude_sites)
    parts, snp_parts, excluded = [], [], dict(unphased_or_bad=0, non_autosome=0, gene=0, site=0, unknown_cell=0)
    cols = ['cell', 'CHROM', 'POS', 'cM', 'AD', 'DP', 'GT', 'gene']
    for f in pd.read_csv(path, sep='\t', usecols=cols, chunksize=chunksize, dtype={'gene': str, 'CHROM': str}):
        pos = pd.to_numeric(f.POS, errors='raise').to_numpy(np.float64)
        if (not np.isfinite(pos).all() or (pos < 1).any() or (pos != np.floor(pos)).any()
                or (pos > np.iinfo(np.int32).max).any()):
            raise ValueError('allele table POS must contain finite positive one-based integer positions '
                             'within the supported int32 range')
        f['POS'] = pos.astype(np.int64)
        for c in ('AD', 'DP'):
            v = pd.to_numeric(f[c], errors='raise').to_numpy(np.float64)
            if (not np.isfinite(v).all() or np.any(v < 0) or np.any(v != np.floor(v))
                    or np.any(v > 2 ** 24)):
                raise ValueError(f'allele table {c} must hold finite nonnegative integer read counts '
                                 'within the exact float32 integer range [0, 2**24]')
        if (f.AD > f.DP).any():
            raise ValueError('allele table AD cannot exceed DP')
        ok = f.GT.isin(['0|1', '1|0']) & (f.DP > 0) & (f.AD >= 0) & (f.AD <= f.DP)
        excluded['unphased_or_bad'] += int((~ok).sum())
        f = f[ok]
        chrom = _chrom(f.CHROM)
        auto = np.isin(chrom, AUTOSOMES)
        excluded['non_autosome'] += int((~auto).sum())
        f, chrom = f[auto], chrom[auto]
        bad_gene = f.gene.fillna('').str.upper().isin(exclude_genes).to_numpy()
        bad_site = np.isin(np.char.add(np.char.add(chrom.astype(str), ':'), f.POS.astype(str).to_numpy()),
                           list(exclude_sites)) if exclude_sites else np.zeros(len(f), bool)
        excluded['gene'] += int(bad_gene.sum())
        excluded['site'] += int((bad_site & ~bad_gene).sum())
        keep = ~(bad_gene | bad_site)
        f, chrom = f[keep], chrom[keep]
        row = index.get_indexer(f.cell.astype(str))
        excluded['unknown_cell'] += int((row < 0).sum())
        good = row >= 0
        a = np.where(f.GT.eq('0|1'), f.AD, f.DP - f.AD)
        # Per read-count row only integers are kept (row, chromosome rank, position, counts); SNP metadata
        # (cM, gene) goes to a per-chunk SNP table, so memory scales with SNPs, not with the strings of every row.
        rank_c = pd.Index(AUTOSOMES).get_indexer(chrom[good]).astype(np.int8)
        pos = f.POS.to_numpy(np.int64)[good]
        parts.append(pd.DataFrame(dict(row=row[good].astype(np.int32), rank=rank_c, pos=pos,
                                       a=a[good].astype(np.int32), d=f.DP.to_numpy()[good].astype(np.int32))))
        snp_parts.append(pd.DataFrame(dict(rank=rank_c, pos=pos, cM=f.cM.to_numpy()[good], gene=f.gene.to_numpy()[good]))
                         .groupby(['rank', 'pos'], sort=False).agg(cM=('cM', 'first'), gene=('gene', 'first'))
                         .reset_index())
    t = pd.concat(parts, ignore_index=True)
    snps = (pd.concat(snp_parts, ignore_index=True).groupby(['rank', 'pos'], sort=False)
            .agg(cM=('cM', 'first'), gene=('gene', 'first')).reset_index())
    snps['chromosome'] = np.asarray(AUTOSOMES, dtype=object)[snps['rank'].to_numpy()]
    snps = snps.sort_values(['rank', 'pos']).reset_index(drop=True)
    key_s = snps['rank'].to_numpy(np.int64) * (1 << 32) + snps.pos.to_numpy(np.int64)
    key_t = t['rank'].to_numpy(np.int64) * (1 << 32) + t.pos.to_numpy(np.int64)
    snps = snps[['chromosome', 'pos', 'cM', 'gene']].copy()
    # Missing genetic positions (e.g. before the first map marker): interpolate within the
    # chromosome from physical position, else 1 cM per Mb.
    snps['cM'] = pd.to_numeric(snps.cM, errors='coerce')
    for c, idx in snps.groupby('chromosome').groups.items():
        cm, bp = snps.loc[idx, 'cM'].to_numpy(), snps.loc[idx, 'pos'].to_numpy(np.float64)
        ok = np.isfinite(cm)
        snps.loc[idx, 'cM'] = (np.interp(bp, bp[ok], cm[ok], left=np.nan, right=np.nan) if ok.sum() >= 2
                               else bp / 1e6)
        miss = ~np.isfinite(snps.loc[idx, 'cM'].to_numpy())
        if miss.any():
            ref = bp[ok][0] if ok.any() else 0.0
            base = cm[ok][0] if ok.any() else 0.0
            snps.loc[np.asarray(idx)[miss], 'cM'] = base + (bp[miss] - ref) / 1e6
    col = pd.Index(key_s).get_indexer(key_t)
    shape = (len(index), len(snps))
    def count_matrix(values):
        # Duplicate cell/site rows are summed exactly before the float32 likelihood channel.
        matrix = sparse.csr_matrix((values.to_numpy(np.int64), (t.row.to_numpy(), col)), shape=shape)
        if matrix.data.size and np.any(matrix.data > 2 ** 24):
            raise ValueError('Aggregated allele counts exceed the exact float32 integer range [0, 2**24]')
        return matrix.astype(np.float32)

    A, D = count_matrix(t.a), count_matrix(t.d)
    return AlleleData(snps=snps, A=A, D=D, excluded=excluded)


def normal_snp_filter(A, D, ref_rows, params: AlleleParams):
    """Keep mask over SNPs and beta-binomial concentration kappa from same-patient reference cells.

    A SNP is dropped when the pooled reference cells reject b = 0.5 (two-sided
    binomial, alpha, at least bias_min_depth reads). kappa: method of moments on
    the remaining SNPs with at least 20 reference reads (Var(a/d) = b(1-b)(1 + (d-1)/(kappa+1))/d)."""
    S = A.shape[1]
    if len(ref_rows) == 0:
        return np.ones(S, bool), params.default_kappa, 0
    a = np.asarray(A[ref_rows].sum(axis=0)).ravel()
    d = np.asarray(D[ref_rows].sum(axis=0)).ravel()
    tested = d >= params.bias_min_depth
    lo = stats.binom.cdf(a, d, 0.5)
    hi = stats.binom.sf(a - 1, d, 0.5)
    p = np.minimum(1.0, 2 * np.minimum(lo, hi))
    keep = ~(tested & (p < params.bias_alpha))
    sel = keep & (d >= 20)
    kappa = params.default_kappa
    if sel.sum() >= 20:
        x = a[sel] / d[sel]
        ratio = (x - 0.5) ** 2 * d[sel] / 0.25                         # E = 1 + (d-1) rho, rho = 1/(kappa+1)
        rho = float(np.sum(ratio - 1.0) / np.sum(d[sel] - 1.0))
        kappa = float(np.clip(1.0 / rho - 1.0, 5.0, 1e4)) if rho > 0 else 1e4
    return keep, kappa, int((tested & ~keep).sum())


def switch_probs(snps: pd.DataFrame, nu):
    """Phase-switch probability between each SNP and the previous one (0.5 at chromosome starts)."""
    cm = snps.cM.to_numpy(np.float64)
    chrom = snps.chromosome.to_numpy()
    if len(cm) == 0:                                  # no usable SNPs (e.g. all removed by the bias filter)
        return np.zeros(0), np.zeros(0, bool)
    reset = np.ones(len(snps), bool)
    reset[1:] = chrom[1:] != chrom[:-1]
    d = np.abs(np.diff(cm, prepend=cm[0]))
    p = 0.5 * (1 - np.exp(-2 * nu * d))
    p[reset] = 0.5
    return np.clip(p, 1e-6, 0.5), reset


def allele_transitions(p_switch, p_change, n_levels):
    """[T, S, S] log transitions; S = 1 + 2 n_levels: balanced, then (+, -) per level."""
    L = n_levels
    S = 1 + 2 * L
    macro = np.array([0] + [1 + l for l in range(L) for _ in (0, 1)])
    orient = np.array([0] + [o for _ in range(L) for o in (1, -1)])
    ps = np.asarray(p_switch, np.float64)[:, None, None]
    same_macro = macro[:, None] == macro[None, :]
    imb_to = (orient != 0)[None, :]
    move = np.where(imb_to, p_change / L / 2, p_change / L)
    T = np.where(same_macro, 0.0, move)[None].repeat(len(ps), axis=0)
    stay = (1 - p_change)
    same_orient = (orient[:, None] == orient[None, :])
    within = np.where(same_orient, stay * (1 - ps), stay * ps)
    T = np.where(same_macro[None] & (orient[:, None] != 0)[None], within, T)
    T[:, 0, 0] = stay
    return jnp.asarray(np.log(np.maximum(T, 1e-300)), jnp.float32)


def state_b(levels=LEVELS):
    return np.concatenate([[0.5], np.ravel([[0.5 + t, 0.5 - t] for t in levels])])


def _emissions(Ak, Dk, kappa, b):
    """Beta-binomial log-likelihood [K, S, states] up to a constant shared by states."""
    a = Ak[..., None]
    n = Dk[..., None]
    alpha, beta = kappa * b, kappa * (1 - b)
    return betaln(a + alpha, n - a + beta) - betaln(alpha, beta)


def clone_imbalance(Ak, Dk, snps, kappa, params: AlleleParams, nu=None, levels=LEVELS):
    """Posteriors [K, S, states], oriented B fraction [K, S], chosen nu and per-clone log-likelihoods.

    With no SNPs (S = 0) the allele channel is empty: zero-width posteriors and B, nu NaN (not estimated)
    and log-likelihood 0, so every downstream allele test sees no evidence."""
    b = state_b(levels)
    Ak, Dk = np.atleast_2d(Ak), np.atleast_2d(Dk)
    if Ak.shape[1] == 0:
        K = Ak.shape[0]
        return (np.zeros((K, 0, len(b))), np.zeros((K, 0)), float('nan') if nu is None else nu, np.zeros(K))
    emit = jnp.asarray(_emissions(np.atleast_2d(Ak), np.atleast_2d(Dk), kappa, b), jnp.float32)
    init = np.full(len(b), params.p_change / (len(b) - 1))
    init[0] = 1 - params.p_change
    init = jnp.log(jnp.asarray(init, jnp.float32))
    best = None
    for v in ([nu] if nu is not None else NU_GRID):
        p, reset = switch_probs(snps, v)
        trans = allele_transitions(p, params.p_change, len(levels))
        post, ll = forward_backward_batch(emit, trans, init, jnp.asarray(reset))
        total = float(np.sum(np.asarray(ll)))
        if best is None or total > best[0]:
            best = (total, v, np.asarray(post), np.asarray(ll))
    _, nu, post, ll = best
    return post, post @ b, nu, ll


def imbalance_runs(post, snps, Dk, params: AlleleParams, levels=LEVELS):
    """Contiguous runs with P(imbalanced) > 0.5; emitted intervals are zero-based half-open."""
    p_imb = 1.0 - post[:, 0]
    theta = (post[:, 1:] @ np.repeat(levels, 2)) / np.maximum(p_imb, 1e-12)
    chrom = snps.chromosome.to_numpy()
    new = np.r_[True, chrom[1:] != chrom[:-1]]
    on = p_imb > 0.5
    starts = np.flatnonzero(on & (new | ~np.r_[False, on[:-1]]))
    if not len(starts):
        return []
    last = np.r_[new[1:], True]
    ends = np.flatnonzero(on & (last | ~np.r_[on[1:], False]))
    cs = lambda x: np.r_[0.0, np.cumsum(x)]
    cd, cp, ct = cs(Dk), cs(p_imb), cs(theta * p_imb)
    pos = snps.pos.to_numpy()
    out = []
    for a, e in zip(starts, ends):
        n, depth, sp = e - a + 1, cd[e + 1] - cd[a], cp[e + 1] - cp[a]
        out.append(dict(chromosome=chrom[a], start=int(pos[a]) - 1, end=int(pos[e]), first_snp=int(a), last_snp=int(e),
                        n_snps=int(n), depth=float(depth), p_imbalanced=float(sp / n),
                        theta=float((ct[e + 1] - ct[a]) / max(sp, 1e-12)),
                        resolved=bool(n >= params.min_snps and depth >= params.min_depth
                                      and sp / n >= params.posterior)))
    return out


def _theta2(A, D, kappa, min_depth=20):
    """Weighted estimate of theta^2 = (b - 0.5)^2 over SNPs and its variance; None if too shallow."""
    ok = D >= 2
    if D[ok].sum() < min_depth or ok.sum() < 3:
        return None
    a, d = A[ok], D[ok]
    rho = 1.0 / (kappa + 1.0)
    s2 = 0.25 * (1 + (d - 1) * rho) / d
    y = (a / d - 0.5) ** 2 - s2
    theta2 = max(float(np.average(y, weights=d)), 0.0)
    var = 2 * s2 ** 2 + 4 * theta2 * s2
    w = 1.0 / var
    est = float(np.sum(w * y) / w.sum())
    return est, float(1.0 / w.sum())


def _delta2(Aa, Da, Ab, Db, kappa, min_depth=20):
    """Weighted estimate of delta^2 = (b_a - b_b)^2 over SNPs and its variance; None if too shallow.

    Signed per-SNP contrast between two groups. Both groups' counts are oriented to the same input phase at
    each SNP, so phase switch errors do not affect it, and it sees mirrored imbalance (b_a = 0.9, b_b = 0.1:
    opposite haplotypes) that the magnitude contrast of `_theta2` cannot. The subtracted sampling variance is
    the b = 0.5 one (the largest), so the estimate is conservative away from balance."""
    ok = (Da >= 2) & (Db >= 2)
    if min(Da[ok].sum(), Db[ok].sum()) < min_depth or ok.sum() < 3:
        return None
    rho = 1.0 / (kappa + 1.0)
    s2 = sum(0.25 * (1 + (d - 1) * rho) / d for d in (Da[ok], Db[ok]))
    y = (Aa[ok] / Da[ok] - Ab[ok] / Db[ok]) ** 2 - s2
    delta2 = max(float(np.average(y, weights=1.0 / s2)), 0.0)
    w = 1.0 / (2 * s2 ** 2 + 4 * delta2 * s2)
    return float(np.sum(w * y) / w.sum()), float(1.0 / w.sum())


def cell_loglik(A, D, B, shrink=0.9, kappa=None):
    """[n, K] allele log-likelihood, omitting the clone-independent combinatorial term.

    A finite positive kappa selects a cell-level beta-binomial observation model. It must be supplied
    independently: pooled-clone overdispersion is not a validated cell-level concentration. None uses the
    legacy binomial model. Neither option integrates uncertainty in estimated clone fractions.
    """
    return _cell_loglik(A, D, B, shrink, kappa, np.float32)


def _cell_loglik(A, D, B, shrink, kappa, dtype):
    """Shared observation model; selection contrasts retain float64 precision."""
    B = 0.5 + shrink * (np.asarray(B, np.float64) - 0.5)
    if not np.isfinite(B).all() or ((B <= 0) | (B >= 1)).any():
        raise ValueError('shrunken allele fractions must be finite and strictly between 0 and 1')
    if kappa is not None:
        if not np.isfinite(kappa) or kappa <= 0:
            raise ValueError('cell-level kappa must be finite and positive')
        A, D = sparse.csr_matrix(A, dtype=np.float64), sparse.csr_matrix(D, dtype=np.float64)
        if (A.shape != D.shape or D.shape[1] != B.shape[1] or not np.isfinite(A.data).all()
                or not np.isfinite(D.data).all() or (A.data < 0).any() or (D.data < 0).any()
                or ((D - A).data < 0).any()):
            raise ValueError('allele counts must satisfy finite 0 <= A <= D and match clone fractions')
        D = D.copy()
        D.eliminate_zeros()
        out = np.zeros((D.shape[0], len(B)), np.float64)
        # Iterate cells: allocations scale with covered sites x clones, never cells x sites x clones.
        for c in range(D.shape[0]):
            lo, hi = D.indptr[c:c + 2]
            idx, depth = D.indices[lo:hi], D.data[lo:hi]
            if not len(idx):
                continue
            a = A[c, idx].toarray().ravel()
            if ((a < 0) | (a > depth)).any() or not np.isfinite(depth).all():
                raise ValueError('allele counts must satisfy finite 0 <= A <= D')
            alpha, beta = kappa * B[:, idx], kappa * (1 - B[:, idx])
            out[c] = (betaln(a[None] + alpha, depth[None] - a[None] + beta)
                      - betaln(alpha, beta)).sum(axis=1)
        return out.astype(dtype)
    la, lb = np.log(B).T, np.log1p(-B).T
    return np.asarray(A @ la + (D - A) @ lb, dtype)


def balance_loglik_ratio(A, D, b, kappa=None):
    """Per-cell log likelihood of already-shrunken fractions b versus balance.

    Uses the explicitly supplied cell-level observation model. Pooled-clone
    concentration is a different quantity and must not be substituted here.
    """
    if kappa is None:
        # Preserve the legacy float64 contrast, including its evaluation order.
        return np.asarray(A @ np.log(b / 0.5) + (D - A) @ np.log((1 - b) / 0.5)).ravel()
    likelihood = _cell_loglik(A, D, np.stack((b, np.full_like(b, .5))), 1., kappa, np.float64)
    return likelihood[:, 0] - likelihood[:, 1]


def segment_alleles(post, snps, Dk, bins, first_bin, last_bin, levels=LEVELS):
    """Allele summary of a zero-based half-open expression segment; SNP pos remains one-based."""
    c = bins.chromosome.iloc[first_bin]
    lo, hi = bins.start.iloc[first_bin], bins.end.iloc[last_bin]
    sel = np.flatnonzero((snps.chromosome.to_numpy() == c) & (snps.pos.to_numpy() > lo) & (snps.pos.to_numpy() <= hi))
    if not len(sel):
        return dict(n_snps=0, allele_depth=0.0, p_imbalanced=np.nan, theta=np.nan)
    w = Dk[sel]
    p_imb = 1.0 - post[sel, 0]
    theta = (post[sel, 1:] @ np.repeat(levels, 2)) / np.maximum(p_imb, 1e-12)
    wt = np.maximum(w, 1e-9)
    return dict(n_snps=int(len(sel)), allele_depth=float(w.sum()), p_imbalanced=float(np.average(p_imb, weights=wt)),
                theta=float(np.average(theta, weights=wt * p_imb + 1e-12)))


def joint_state(cn_state, p_imbalanced, n_snps, allele_depth, params: AlleleParams):
    """Combine relative dosage state (1-4) with allele evidence into a joint label.

    Labels describe the evidence (relative RNA dosage x allelic imbalance), not a DNA genotype:
    'neutral_allele_imbalanced' is a cnLOH candidate, not a cnLOH call."""
    if n_snps < params.min_snps or allele_depth < params.min_depth or not np.isfinite(p_imbalanced):
        return {1: 'loss', 2: 'neutral', 3: 'gain', 4: 'amplification'}[cn_state] + '_alleles_uninformative'
    imb = p_imbalanced >= params.posterior
    bal = p_imbalanced <= 1 - params.posterior
    if cn_state == 2:
        return 'neutral_allele_imbalanced' if imb else ('neutral' if bal else 'neutral_alleles_unclear')
    if cn_state == 1:
        return 'loss' if imb else ('loss_allele_balanced' if bal else 'loss_alleles_unclear')
    if cn_state == 3:
        return 'gain' if imb else ('gain_allele_balanced' if bal else 'gain_alleles_unclear')
    return 'amplification_imbalanced' if imb else ('balanced_amplification' if bal else 'amplification_alleles_unclear')


class SpecimenAlleles:
    """One specimen's allele counts (rows aligned with the specimen's cells) and clone-level state."""

    def __init__(self, A, D, snps, kappa, params: AlleleParams | None = None, n_biased=0):
        self.A, self.D = sparse.csr_matrix(A, dtype=np.float32), sparse.csr_matrix(D, dtype=np.float32)
        self.snps = snps.reset_index(drop=True)
        self.kappa, self.params, self.n_biased = float(kappa), params or AlleleParams(), int(n_biased)
        self.nu = None
        self.bins = None                # expression bins (set by the pipeline) for segment-level allele checks

    def pool(self, R, rows=None):
        """Pooled oriented and total counts [K, S] for responsibilities (or one-hot) R over `rows` (default all)."""
        R = np.asarray(R, np.float64)
        A, D = (self.A, self.D) if rows is None else (self.A[rows], self.D[rows])
        return np.asarray((A.T @ R).T), np.asarray((D.T @ R).T)

    def fit(self, Ak, Dk):
        """HMM posteriors and oriented B fractions per pooled group.

        On first use, nu and kappa are chosen by likelihood: nu over NU_GRID, kappa
        over multiples of the reference estimate (reference cells of another type
        can be more overdispersed than the tumour's own normal epithelium)."""
        if self.nu is None:
            best = None
            for f in KAPPA_FACTORS:
                post, B, nu, ll = clone_imbalance(Ak, Dk, self.snps, self.kappa * f, self.params)
                if best is None or ll.sum() > best[0]:
                    best = (ll.sum(), f, nu, post, B)
            _, f, self.nu, post, B = best
            self.kappa *= f
            return post, B
        post, B, _, _ = clone_imbalance(Ak, Dk, self.snps, self.kappa, self.params, nu=self.nu)
        return post, B

    def runs(self, post_k, Dk_k):
        return imbalance_runs(post_k, self.snps, Dk_k, self.params)

    def has_events(self, post_k, Dk_k):
        return any(r['resolved'] for r in self.runs(post_k, Dk_k))

    def differ(self, runs_a, Aa, Da, runs_b, Ab, Db, z_min=4.0, min_theta_diff=0.1):
        """Summed |z| of allele-fraction differences between two pooled groups (0: none resolved).

        Tested over each resolved run of either group, on odd SNPs only (splits are
        proposed on even SNPs). First the imbalance magnitude: per SNP, (x - 0.5)^2
        minus its sampling variance estimates theta^2, with x = a / d of the group.
        If that does not resolve a difference, the signed per-SNP contrast (`_delta2`)
        is tried, which also needs no phase and sees mirrored imbalance (the same
        |b - 0.5| on opposite haplotypes, review 2026-11-02 #2)."""
        total, seen = 0.0, set()
        for r in list(runs_a) + list(runs_b):
            key = (r['first_snp'], r['last_snp'])
            if not r['resolved'] or key in seen:
                continue
            seen.add(key)
            idx = np.arange(r['first_snp'], r['last_snp'] + 1)
            idx = idx[idx % 2 == 1]
            est = [_theta2(A[idx], D[idx], self.kappa) for A, D in ((Aa, Da), (Ab, Db))]
            if any(e is None for e in est):
                continue
            (t1, v1), (t0, v0) = est
            z = (t1 - t0) / np.sqrt(v1 + v0)
            gap = abs(np.sqrt(max(t1, 0.0)) - np.sqrt(max(t0, 0.0)))
            if abs(z) >= z_min and gap >= min_theta_diff:
                total += abs(z)
                continue
            sd = _delta2(Aa[idx], Da[idx], Ab[idx], Db[idx], self.kappa)
            if sd is not None:
                z = sd[0] / np.sqrt(sd[1])
                if z >= z_min and np.sqrt(max(sd[0], 0.0)) >= min_theta_diff:
                    total += z
        return total

    def snps_in(self, bins, first_bin, last_bin):
        """One-based SNP positions inside zero-based half-open expression bins."""
        c = bins.chromosome.iloc[first_bin]
        lo, hi = bins.start.iloc[first_bin], bins.end.iloc[last_bin]
        pos, chrom = self.snps.pos.to_numpy(), self.snps.chromosome.to_numpy()
        return np.flatnonzero((chrom == c) & (pos > lo) & (pos <= hi))

    def veto(self, idx, Aa, Da, Ab, Db, power_z=3.0, null_z=1.0, theta=1 / 6, log_dose=None, loh_theta=0.375):
        """True if SNPs idx could see a one-copy imbalance change (theta 0 -> 1/6) but the groups agree.

        A one-copy change from a state with both haplotypes present (m:M, m >= 1) changes
        |b - 0.5|, so an expression difference with well-powered, unchanged imbalance is not
        taken as copy-number evidence. Changes that keep the allele ratio are invisible to
        alleles and are never vetoed (review 2026-11-02 #2):
        - both groups near LOH (|b - 0.5| >= loh_theta): 1:0 vs 2:0 both read b = 0 or 1;
        - with `log_dose` (the segment's log dosage difference between the groups), a
          difference nearer a doubling than a one-copy gain from two copies
          (|log_dose| >= log 1.75), as 1:1 vs 2:2. The one-copy change of that size
          (1:1 vs 1:0) changes imbalance, so it is not vetoed anyway.
        Expression differences of that size are therefore never vetoed by alleles. This is a
        compatibility check on two summaries, not a joint allele/copy state model."""
        if len(idx) < self.params.min_snps:
            return False
        if log_dose is not None and abs(log_dose) >= np.log(1.75):
            return False
        est = [_theta2(A[idx], D[idx], self.kappa) for A, D in ((Aa, Da), (Ab, Db))]
        if any(e is None for e in est):
            return False
        (t1, v1), (t0, v0) = est
        if min(t1, t0) >= loh_theta ** 2:
            return False
        rho = 1.0 / (self.kappa + 1.0)
        null_var = []
        for D in (Da, Db):
            d = D[idx][D[idx] >= 2]
            s2 = 0.25 * (1 + (d - 1) * rho) / d
            null_var.append(1.0 / np.sum(1.0 / (2 * s2 ** 2 + 4 * theta ** 2 * s2)))
        power = theta ** 2 / np.sqrt(sum(null_var))
        return bool(power >= power_z and abs(t1 - t0) / np.sqrt(v1 + v0) < null_z)

    def cell_loglik(self, rows, B):
        return cell_loglik(self.A[rows], self.D[rows], B, self.params.cell_shrink, self.params.cell_kappa)


def balanced_cells(alleles: SpecimenAlleles, min_theta=0.15, posterior=0.99):
    """Cells whose reads in the specimen's strong imbalance runs are confidently balanced.

    The pooled specimen is fitted once; runs with theta >= min_theta (LOH, single-copy
    gains/losses) orient each SNP. A cell's log-likelihood ratio of that oriented
    imbalance against b = 0.5 gives P(balanced) with a flat prior; cells at or above
    `posterior` are allele-balanced candidates. The explicitly selected cell
    observation model is used; these plug-in scores are not calibrated confidence
    in normal identity. Balanced gains and genome doubling
    remain possible, so balance does not establish copy neutrality or normal identity.
    Returns (row positions, summary dict)."""
    Ak, Dk = alleles.pool(np.ones((alleles.A.shape[0], 1)))
    post, B = alleles.fit(Ak, Dk)
    runs = [r for r in alleles.runs(post[0], Dk[0]) if r['resolved'] and r['theta'] >= min_theta]
    if not runs:
        return np.array([], dtype=int), dict(runs=0, snps=0)
    idx = np.concatenate([np.arange(r['first_snp'], r['last_snp'] + 1) for r in runs])
    b = 0.5 + alleles.params.cell_shrink * (B[0, idx] - 0.5)
    A, D = alleles.A[:, idx], alleles.D[:, idx]
    llr = balance_loglik_ratio(A, D, b, alleles.params.cell_kappa)
    p_bal = 1.0 / (1.0 + np.exp(np.clip(llr, -50, 50)))
    return np.flatnonzero(p_bal >= posterior), dict(runs=len(runs), snps=int(len(idx)))
