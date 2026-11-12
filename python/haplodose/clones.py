"""Per-specimen copy-number clones.

Segmentation happens on pooled cells, where depth supports it; single cells are
only assigned to clone profiles by likelihood. Steps:

1. Baseline: mixture of reference groups fitted to the pooled specimen counts.
2. Pooled specimen profile: centre, dosage slope beta and bin-level noise tau.
3. Proposals: Leiden groups on windowed, depth-regressed log ratios.
4. Merge groups whose pooled profiles differ in no resolved copy-number segment.
5. EM: cell posteriors over clone profiles; re-pool; re-segment.
6. Divisive refinement: re-cluster within clones and look for cells sharing an
   arm-level shift; keep a split only if the pooled difference is resolved.

Sequencing depth changes bin composition (gene-dense chromosomes take a different
share in shallow cells). A per-bin depth slope gamma_j, estimated within clones,
enters the cell likelihood, pooled expectations and difference tests to model
depth-associated variation. Expression programmes learned on normal cells
(`haplodose.programmes`) enter the cell likelihood as per-cell nuisance scores.
Pooled profiles and differences also fit programme amplitudes. These working
corrections do not guarantee that depth or cell state cannot create candidate
groups or that genuine DNA effects are preserved. Each pooled group carries
XL = (sum of L_c x_c, sum of L_c x_c^2), with
x_c the centred log depth: the UMI-weighted mean and variance of depth, from which the pooled expectation is
built (`SpecimenModel.expected`).

Proposals use even bins and difference tests use odd bins. This separates immediate
features, not the complete adaptive search: full-count EM fits memberships and nuisance
quantities, and test bins are reused across cycles. It does not establish independent
confirmation or a specimen-wide error rate. Clone 0 is reference-like, not verified normal.

`CloneParams.split_mode = 'thin'` replaces the even/odd bin split by negative-binomial count splitting
(data thinning; Neufeld et al., for convolution-closed families): each cell x bin count Y is split as
Y_prop | Y ~ BetaBinomial(Y, eps b_j, (1 - eps) b_j), Y_test = Y - Y_prop, with b_j = 1 / phi_j and eps = 0.5.
If Y ~ NB(mean m, size b), the halves are independent NB(eps m, eps b) and NB((1 - eps) m, (1 - eps) b). The
independence is exact only for the true size b; b_j here is the per-bin size from the reference fit (frozen at
model start), so it is approximate. Proposals (Leiden features, event mixtures) use Y_prop on every bin; split
tests, the agglomeration of proposed groups, merge tests and reported differences use Y_test on every bin. EM,
profiles and cell assignment use the full counts. Bin-count settings defined on proposal or test bins
(window_bins, event_windows, the difference test's min_bins and p_change) are converted so they cover the same
genomic span as in 'evenodd' (`_thin_diff_params`). Under the NB model, eps of the counts on every bin carry
about as much information about a region's rate as all counts on half its bins, so the expected gain is finer
localisation (no bins lost at event boundaries), not more counts [Inference]. 'thin' is opt-in: on the Kinker
v3 subclone design (development data, 166 lines) it matched 'evenodd' (mean ARI 0.165 vs 0.167, paired Wilcoxon
p 0.50) with more one-clone single lines (85 vs 70 of 181) (results/vnext_thin_20261102_1600/compare/).

`CloneParams.split_null_draws > 0` adds a model-based filter for divisive splits (`split_null`): a
clone's best split must beat the best split found by the same proposal-and-test search on that many simulations
from the fitted one-clone NB model. It does not gate the initial clustering (`fit_specimen`), which can still
start with groups that only the merge tests can join. Its simulated search conditions on selected parent
cells and fitted nuisance quantities; no specimen-wide false-split control is established. Allele-enabled
runs cannot use this expression-only null until a compatible joint search and null are implemented.

With allele counts (`haplodose.alleles`), each clone also carries pooled phased
SNP counts and an allelic-imbalance HMM. Cells add their SNP evidence to the
clone likelihood; a clone is kept if it has a resolved dosage or allele event
(so copy-neutral LOH clones survive); groups merge only if they differ in
neither; and cells sharing a clone's imbalance are proposed as splits.

Spatial data: a Potts prior on clone labels over a within-specimen nearest-
neighbour graph G (rows sum to 1). The E-step is mean-field:
R_c ~ pi * exp(LL_c + lambda * (G R)_c), iterated. lambda is chosen per E-step
by how well neighbours' data-only posteriors predict each cell's own data
(0 = no smoothing), so smoothing is only as strong as the data's spatial coherence.
Expression is never averaged across neighbours.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.cluster.hierarchy import linkage
from scipy.special import logsumexp

from .kernels.estep import cell_loglik, posteriors, programme_scores
from .profile import (ProfileParams, _subset_reset, difference_segments, fit_pooled,
                      group_centres, group_tau, log_ratio, runs, segment, segments_table,
                      weighted_median, pooled_lograte)
from .alleles import balance_loglik_ratio, state_b
from .programmes import fit_amplitudes


@dataclass
class CloneParams:
    min_cells: int = 30
    min_umi: float = 100.0
    leiden_resolution: float = 1.0
    max_proposal_cells: int = 20000
    window_bins: int = 10
    n_pcs: int = 20
    n_neighbors: int = 15
    max_groups: int = 40
    outer_rounds: int = 8
    cycles: int = 3
    normal_prior: float = 0.05
    event_candidates: int = 5
    event_windows: tuple = ()             # sliding-window region sizes in proposal (even) bins; off: no gain on
                                          # the Kinker subclone test (v2, v3) and 2.7x slower at 100k cells
    in_sample_reference: bool = False     # opt-in: locus-wise allele balance does not establish normal identity
    spatial_neighbors: int = 8
    potts_lambda: float = -1.0           # < 0: chosen from neighbour-predicted cell likelihoods
    n_programmes: int = 10
    programme_window: int = 0            # bins of running mean removed; 0 = arm means (hg38 centromeres)
    depth_source: str = 'query'          # 'reference': depth slopes from fit reference cells, held fixed
    seed: int = 20261031
    split_mode: str = 'evenodd'           # 'evenodd': proposals on even bins, tests on odd bins; 'thin': NB count
                                          # splitting of every bin (module docstring)
    split_null_draws: int = 0             # > 0: accept a clone's best split only if its score beats this many
                                          # fitted-parent NB draws; no specimen-wide guarantee; 0 off
    split_null_model: str = 'nb'          # 'nb': fitted NB only; 'cell_region': plus a per-cell, per-region
                                          # (chromosome half) log-normal rate effect with moment-matched variance;
                                          # 'cov': the same with the full moment-matched region covariance
    split_null_initial: bool = False      # with the null on: pool the initial tumour groups into one clone, so every
                                          # tumour clone has to come from a null-gated split
    profile: ProfileParams = field(default_factory=ProfileParams)


THIN_EPS = 0.5                            # share of each count given to the proposal half in split_mode 'thin'
SPLIT_MODES = ('evenodd', 'thin')


def chromosome_halves(reset):
    """[J] region id per bin: the first and second half of each chromosome's bins (as in event_mixtures)."""
    chrom = np.cumsum(reset) - 1
    out = np.empty(len(reset), dtype=int)
    for c in np.unique(chrom):
        idx = np.flatnonzero(chrom == c)
        out[idx] = 2 * c + (np.arange(len(idx)) >= len(idx) // 2)
    return np.unique(out, return_inverse=True)[1]


def region_sums(X, region):
    """[n, J] -> [n, regions] sums of X over each region's bins."""
    onehot = sparse.csr_matrix((np.ones(len(region)), (np.arange(len(region)), region)))
    return np.asarray(sparse.csr_matrix(X) @ onehot) if sparse.issparse(X) else np.asarray(X @ onehot)


def region_effect_variance(Y, m, phi, region):
    """Moment estimate of s2_r, the variance of a per-cell log-normal rate effect shared by region r's bins.

    Under Y_cj ~ NB(m_cj exp(e_cr), 1 / phi_j), e_cr ~ N(0, s2_r) independent over cells, the region sum O_cr has
    E (O_cr - M_cr)^2 ~ V_cr + M_cr^2 s2_r, with M_cr = sum_j m_cj and V_cr = sum_j (m_cj + phi_j m_cj^2)
    (small s2_r). So s2_r = sum_c [(O_cr - M_cr)^2 - V_cr] / sum_c M_cr^2, floored at 0."""
    M = region_sums(m, region)
    V = region_sums(m + np.asarray(phi)[None] * m * m, region)
    O = region_sums(np.asarray(Y, np.float64), region)
    return np.maximum(((O - M) ** 2 - V).sum(axis=0) / np.maximum((M * M).sum(axis=0), 1e-300), 0.0)


def region_effect_covariance(Y, m, phi, region):
    """Moment estimate of S, the [regions, regions] covariance of per-cell log rate effects e_c shared by each
    region's bins (as `region_effect_variance`, but keeping correlations between regions).

    z_cr = log((O_cr + 0.5) / (M_cr + 0.5)) ~ e_cr + sampling noise of variance ~ V_cr / M_cr^2, independent over
    regions given e_c. So S = Cov_c(z) - diag(mean_c V_cr / M_cr^2), projected onto the positive semidefinite cone
    (negative eigenvalues set to 0). The null in SigClust (Liu et al. 2008, JASA) uses a Gaussian with the data's
    covariance in the same way: any unimodal, correlated cell-by-region structure is in the null, and a split has
    to beat it."""
    M = region_sums(m, region)
    V = region_sums(m + np.asarray(phi)[None] * m * m, region)
    O = region_sums(np.asarray(Y, np.float64), region)
    z = np.log((O + 0.5) / (M + 0.5))
    S = np.cov(z, rowvar=False) - np.diag((V / np.maximum(M + 0.5, 1e-300) ** 2).mean(axis=0))
    w, Q = np.linalg.eigh((S + S.T) / 2)
    return (Q * np.maximum(w, 0.0)) @ Q.T


def thin_counts(Y, keys, size, seed, eps=THIN_EPS):
    """Proposal half of NB count splitting: Y_prop | Y ~ BetaBinomial(Y, eps size_j, (1 - eps) size_j).

    Independence of the two parts requires the NB law with its true fixed size.
    Heterogeneous gene sums or plug-in size estimates may violate that condition.
    Y [n, J] counts; keys [n] global row indices. Each row is drawn from its own Philox stream keyed by
    (seed, key), so a row's draw does not depend on chunking or on which other rows are requested.
    Returns float32 [n, J]; the test half is Y - Y_prop."""
    Y = np.asarray(Y)
    if Y.size and not (np.isfinite(Y).all() and Y.min() >= 0 and np.abs(Y - np.rint(Y)).max() <= 1e-3):
        raise ValueError('thin_counts needs finite, nonnegative, integer-valued counts')
    size = np.asarray(size, np.float64)
    a_all, b_all = eps * size, (1.0 - eps) * size
    base = (int(seed) % 2 ** 64) << 64
    out = np.zeros(Y.shape, np.float32)
    for i, key in enumerate(np.asarray(keys)):
        nz = np.flatnonzero(Y[i] > 0)
        if len(nz) == 0:
            continue
        g = np.random.Generator(np.random.Philox(key=base + int(key)))
        p = g.beta(a_all[nz], b_all[nz])
        out[i, nz] = g.binomial(np.rint(Y[i, nz]).astype(np.int64), p)
    return out


def thin_gene_counts(Xg, gbin, Ybin, Ytest, seed):
    """Test-half gene counts consistent with the bin-level split.

    Given a cell x bin count Y split into Y_prop + Y_test, the Y_test reads of that bin are a uniform random
    subset of its Y reads (beta-binomial thinning is exchangeable over reads), so gene counts given Y_test are
    multivariate hypergeometric. Drawn as sequential conditional hypergeometrics over the bin's genes. Xg [n,
    genes] CSR counts of the used genes, gbin [genes] their bins, Ybin / Ytest [n, J] full and test-half bin
    counts. If a bin's gene counts do not sum to its bin count, the test total is capped at the gene total.
    Returns CSR float64 [n, genes]."""
    X = sparse.csr_matrix(Xg).tocoo()
    keep = (X.data > 0) & (np.asarray(gbin)[X.col] >= 0)
    i, g, v = X.row[keep], X.col[keep], np.rint(X.data[keep]).astype(np.int64)
    j = np.asarray(gbin)[g]
    order = np.lexsort((g, j, i))
    i, g, v, j = i[order], g[order], v[order], j[order]
    if len(v) == 0:
        return sparse.csr_matrix(X.shape, dtype=np.float64)
    start = np.r_[True, (i[1:] != i[:-1]) | (j[1:] != j[:-1])]
    grp = np.cumsum(start) - 1
    first = np.flatnonzero(start)
    pos = np.arange(len(v)) - first[grp]
    gsum = np.bincount(grp, weights=v).astype(np.int64)
    after = gsum[grp] - (np.cumsum(v) - np.r_[0, np.cumsum(v)[first[1:] - 1]][grp])   # later genes' counts
    remaining = np.minimum(np.rint(np.asarray(Ytest)[i[first], j[first]]).astype(np.int64), gsum)
    rng = np.random.default_rng([int(seed) % 2 ** 32, 7])
    out = np.zeros(len(v), np.int64)
    for k in range(int(pos.max()) + 1):
        sel = np.flatnonzero(pos == k)
        h = rng.hypergeometric(v[sel], after[sel], remaining[grp[sel]])
        out[sel] = h
        remaining[grp[sel]] -= h
    return sparse.csr_matrix((out.astype(np.float64), (i, g)), shape=X.shape)


class Rows:
    """Bounded-chunk access to one specimen's rows of the cells x bins matrix.

    fold: None (full counts), 'prop' or 'test' (the halves of NB count splitting, drawn per chunk on the fly
    after `set_thinning`; nothing beyond the current chunk is held in memory)."""

    def __init__(self, Y, rows, chunk=4096):
        self.Y, self.rows, self.chunk = Y, np.asarray(rows), int(chunk)
        self.n = len(self.rows)
        self.thin = None

    def set_thinning(self, seed, size, eps=THIN_EPS):
        self.thin = (int(seed), np.asarray(size, np.float64), float(eps))

    def _fold(self, Yc, keys, fold):
        if fold is None:
            return Yc
        if self.thin is None:
            raise ValueError('fold requested without set_thinning')
        seed, size, eps = self.thin
        prop = thin_counts(Yc, keys, size, seed, eps)
        if fold == 'prop':
            return prop
        if fold == 'test':
            return Yc - prop
        raise ValueError(f'unknown fold {fold!r}')

    def chunks(self, fold=None):
        for a in range(0, self.n, self.chunk):
            b = min(a + self.chunk, self.n)
            yield a, b, self._fold(np.asarray(self.Y[self.rows[a:b]], dtype=np.float32), self.rows[a:b], fold)

    def take(self, idx, fold=None):
        keys = self.rows[np.asarray(idx)]
        return self._fold(np.asarray(self.Y[keys], dtype=np.float32), keys, fold)

    def weighted_sum(self, R, fold=None):
        """R [n, G] dense or scipy sparse -> [G, bins] float64. With a fold, only rows with nonzero weight are
        thinned (the draws are per row, so the result is the same)."""
        out = np.zeros((R.shape[1], self.Y.shape[1]), np.float64)
        if fold is None:
            for a, b, Yc in self.chunks():
                out += np.asarray(R[a:b].T.astype(np.float64) @ Yc)
            return out
        for a in range(0, self.n, self.chunk):
            b = min(a + self.chunk, self.n)
            Rc = R[a:b]
            live = (np.diff(Rc.tocsr().indptr) > 0) if sparse.issparse(Rc) else np.asarray(Rc != 0).any(axis=1)
            live = np.flatnonzero(live)
            if len(live) == 0:
                continue
            Yc = self.take(a + live, fold)
            out += np.asarray(Rc[live].T.astype(np.float64) @ Yc)
        return out


def depth_moments(L, x):
    """Per-cell [n, 2] (L x, L x^2): summed over a group, the UMI-weighted first and second depth moments."""
    L, x = np.asarray(L, np.float64), np.asarray(x, np.float64)
    return np.column_stack([L * x, L * x * x])


def _depth_mean_var(Lk, XLk):
    L = np.maximum(np.atleast_1d(Lk), 1.0)
    XL = np.atleast_2d(XLk)
    xbar = XL[:, 0] / L
    return xbar, np.maximum(XL[:, 1] / L - xbar ** 2, 0.0)


def nb_rate(O, E, r, w, theta=None, iters=100, tol=1e-10):
    """Weighted NB maximum-likelihood rate ratio: theta solving sum_c w_c r_c (O_c - theta E_c) / (theta E_c + r_c) = 0,
    for O_c ~ NB(mean theta E_c, size r_c). Fixed-point iteration theta = sum w r O / (theta E + r) / sum w r E / (theta E + r);
    the pooled ratio sum w O / sum w E (the Poisson answer) is the start and the limit as r grows. O, E, r and w
    broadcast (e.g. cells x bins with per-bin r and per-cell w [n, 1])."""
    O, E = np.asarray(O, np.float64), np.asarray(E, np.float64)
    w = np.broadcast_to(np.asarray(w, np.float64), O.shape)
    r = np.broadcast_to(np.asarray(r, np.float64), O.shape)
    theta = float((w * O).sum() / max((w * E).sum(), 1e-12)) if theta is None else float(theta)
    theta = max(theta, 1e-12)
    wr = w * r
    for _ in range(iters):
        d = theta * E + r
        new = float((wr * O / d).sum() / max((wr * E / d).sum(), 1e-300))
        new = max(new, 1e-12)
        if abs(new - theta) <= tol * theta:
            return new
        theta = new
    return theta


def _event_mixture_objective(x, O, E, r, null_logtheta):
    """Negative NB-mixture gain and gradient, relative to one fitted component.

    The omitted NB constants never enter the calculation. Close component rates
    use log1p/expm1 so subtraction of two large log likelihoods is unnecessary.
    """
    from scipy.special import expit
    with np.errstate(divide='ignore'):
        logE = np.log(E)
    base = logE + null_logtheta
    logr = np.log(r)
    q0 = expit(base - logr)
    values, scores = [], []
    for eta in x[:2]:
        delta = eta - null_logtheta
        if abs(delta) < 0.5:
            denominator_ratio = np.log1p(q0 * np.expm1(delta))
        else:
            denominator_ratio = np.logaddexp(logr, base + delta) - np.logaddexp(logr, base)
        values.append((O * delta - (O + r) * denominator_ratio).sum(axis=1))
        scores.append((O - (O + r) * expit(logE + eta - logr)).sum(axis=1))
    pi = expit(x[2])
    a = values[0] - np.logaddexp(0., x[2])
    b = values[1] - np.logaddexp(0., -x[2])
    total = np.logaddexp(a, b)
    responsibility = np.exp(b - total)
    gradient = np.array([np.dot(1 - responsibility, scores[0]),
                         np.dot(responsibility, scores[1]),
                         (responsibility - pi).sum()])
    return -float(total.sum()), -gradient, responsibility


def _fit_event_mixture(O, E, r, step, *, maxiter=1000):
    """Polish the historical EM proposal with an additive-constant invariant fit.

    This certifies local numerical stationarity, not a global mixture optimum.
    The mixing-weight bounds are the historical event-proposal bounds.
    """
    from scipy.optimize import brentq, minimize
    from scipy.special import expit
    if isinstance(maxiter, (bool, np.bool_)) or not isinstance(maxiter, (int, np.integer)) or maxiter <= 0:
        raise ValueError('Event mixture maxiter must be a positive integer')
    O, E = np.asarray(O, np.float64), np.asarray(E, np.float64)
    r = np.broadcast_to(np.asarray(r, np.float64), O.shape)
    if (O.shape != E.shape or O.ndim != 2 or O.size == 0 or not np.isfinite(O).all()
            or not np.isfinite(E).all() or not np.isfinite(r).all()
            or (O < 0).any() or (E < 0).any() or (r <= 0).any()
            or ((E == 0) & (O > 0)).any()):
        raise ValueError('Event mixtures require finite nonnegative counts and supported means, positive sizes')
    if O.sum() == 0:
        return dict(gain=0., responsibility=np.full(len(O), .5),
                    parameters=None, converged=True, iterations=0, numerical_gain_budget=0.)
    with np.errstate(divide='ignore'):
        logE = np.log(E)
    logr = np.log(r)

    def null_score(eta):
        return float((O - (O + r) * expit(logE + eta - logr)).sum())

    lo, hi = -1., 1.
    while null_score(lo) <= 0:
        lo *= 2
    while null_score(hi) >= 0:
        hi *= 2
    eps = np.finfo(np.float64).eps
    null = brentq(null_score, lo, hi, xtol=eps, rtol=4 * eps)
    theta = np.exp(null + np.array([-.5, .5]) * step)
    pi = .5
    # EM only initializes the numerical optimizer. A bounded warmup does not
    # substitute for the final score and observed-curvature convergence checks.
    for _ in range(20):
        x = np.r_[np.log(theta), np.log(pi) - np.log1p(-pi)]
        _, _, responsibility = _event_mixture_objective(x, O, E, r, null)
        pi = float(np.clip(responsibility.mean(), 1e-3, 1 - 1e-3))
        theta = np.array([nb_rate(O, E, r, (1 - responsibility)[:, None], theta[0]),
                          nb_rate(O, E, r, responsibility[:, None], theta[1])])
    x = np.r_[np.log(theta), np.log(pi) - np.log1p(-pi)]
    bound = np.log(1 - 1e-3) - np.log(1e-3)

    def objective(x):
        value, gradient, _ = _event_mixture_objective(x, O, E, r, null)
        return value / len(O), gradient / len(O)

    fit = minimize(objective, x, jac=True, method='L-BFGS-B',
                   bounds=[(None, None), (None, None), (-bound, bound)],
                   options=dict(maxiter=maxiter, ftol=eps, gtol=eps ** .75, maxls=50))
    solver = 'L-BFGS-B'
    if not fit.success and np.isfinite(fit.x).all():
        # A line search can fail at machine precision without exhausting the
        # information in a proposal. Retry the same bounded objective; the
        # score, curvature and gain-error requirements below remain mandatory.
        fit = minimize(objective, fit.x, jac=True, method='SLSQP',
                       bounds=[(None, None), (None, None), (-bound, bound)],
                       options=dict(maxiter=maxiter, ftol=eps))
        solver = 'SLSQP_retry'
    value, gradient, responsibility = _event_mixture_objective(fit.x, O, E, r, null)
    # Conservative float64 summation/elementary-operation budget for this gain,
    # not a fraction of the arbitrarily large NB log-pmf baseline.
    delta = np.abs(fit.x[:2] - null).max()
    magnitude = float(((2 * O + r) * delta).sum()) + len(O) * (1 + abs(fit.x[2]))
    operations = O.size + 16 * len(O) + 32
    gamma = operations * eps / (1 - operations * eps)
    budget = gamma * max(1., magnitude)
    # Complete-data information bounds the score's natural curvature scale.
    information, scores, curvatures = [], [], []
    for eta, w in zip(fit.x[:2], (1 - responsibility, responsibility)):
        q = expit(logE + eta - logr)
        curvature = ((O + r) * q * (1 - q)).sum(axis=1)
        curvatures.append(curvature)
        scores.append((O - (O + r) * q).sum(axis=1))
        information.append(float(np.dot(w, curvature)))
    pi = expit(fit.x[2])
    information.append(len(O) * pi * (1 - pi))
    projected = gradient.copy()
    if fit.x[2] <= -bound and projected[2] > 0 or fit.x[2] >= bound and projected[2] < 0:
        projected[2] = 0.
    score_budget = np.sqrt(2 * budget * np.maximum(information, 1.))
    # Observed information checks weak curvature and saddle points separately:
    # a tiny score alone is not evidence of an accurately located mixture fit.
    t = responsibility * (1 - responsibility)
    s0, s1 = scores
    hessian = np.array([
        [information[0] - np.dot(t, s0 ** 2), np.dot(t, s0 * s1), np.dot(t, s0)],
        [np.dot(t, s0 * s1), information[1] - np.dot(t, s1 ** 2), -np.dot(t, s1)],
        [np.dot(t, s0), -np.dot(t, s1), information[2] - t.sum()]])
    free = np.array([True, True, projected[2] == gradient[2]])
    curvature, directions = np.linalg.eigh(hessian[np.ix_(free, free)])
    curvature_budget = gamma * max(1., float(np.abs(hessian).sum()), sum(information))
    components = directions.T @ projected[free]
    positive = curvature > curvature_budget
    local_gain_error = .5 * float(np.sum(components[positive] ** 2 / curvature[positive]))
    flat_score_budget = gamma * max(1., float(O.sum() + r.sum()))
    curvature_ok = (np.all(curvature >= -curvature_budget)
                    and np.all(np.abs(components[~positive]) <= flat_score_budget)
                    and local_gain_error <= budget)
    if (not fit.success or not np.isfinite(value) or not np.isfinite(projected).all()
            or np.any(np.abs(projected) > score_budget) or not curvature_ok):
        raise RuntimeError(f'Event mixture did not converge: {fit.message}; '
                           f'projected score={projected}, score budget={score_budget}; '
                           f'curvature={curvature}, local gain error={local_gain_error}, budget={budget}')
    return dict(gain=-value, responsibility=responsibility, parameters=fit.x,
                converged=True, solver=solver, iterations=int(fit.nit), numerical_gain_budget=budget,
                projected_score=projected, score_budget=score_budget, null_logtheta=null,
                observed_information_eigenvalues=curvature, local_quadratic_gain_error=local_gain_error)


def mixture_weights(Ytot, mu, iters=500, *, return_diagnostics=False):
    """Constrained Poisson MLE, certified by a global concave-objective gap.

    Columns of mu are reference compositions. Normalizing counts makes the fit
    invariant to their common scale. For h(w)=sum_j p_j log((mu w)_j), the
    Frank-Wolfe gap max(gradient)-w@gradient bounds h(optimum)-h(w). Active-set
    Newton steps and exact line searches replace a fixed number of EM updates.
    An iteration cap never implies convergence: an uncertified fit raises.
    Weights can be nonidentifiable for redundant references; only the fitted
    mixture and objective, not unique biological proportions, are claimed.
    """
    from scipy.optimize import brentq

    raw_mu = np.asarray(mu)
    precision = (np.finfo(raw_mu.dtype).eps if np.issubdtype(raw_mu.dtype, np.floating)
                 else np.finfo(np.float64).eps)
    Ytot, mu = np.asarray(Ytot, np.float64), np.asarray(mu, np.float64)
    if (Ytot.ndim != 1 or mu.ndim != 2 or mu.shape[0] != len(Ytot) or not mu.shape[1]
            or not len(Ytot) or not np.isfinite(Ytot).all() or (Ytot < 0).any()
            or not np.isfinite(mu).all() or (mu < 0).any()):
        raise ValueError('Mixture counts and reference compositions must be finite, nonnegative and aligned')
    if not np.any(Ytot > 0):
        raise ValueError('Mixture fitting requires positive total counts')
    if int(iters) != iters or iters < 1:
        raise ValueError('Mixture iteration safeguard must be a positive integer')
    columns = mu.sum(axis=0)
    column_error = len(mu) * precision / max(1 - len(mu) * precision, precision)
    if (columns <= 0).any() or not np.isfinite(columns).all() or np.any(np.abs(columns - 1) > column_error):
        raise ValueError('Each reference composition column must sum to one')
    if np.any((Ytot > 0) & (mu.max(axis=1) == 0)):
        raise ValueError('Positive observed counts have no support in any reference composition')
    mu = mu / columns[None]
    counts = Ytot / Ytot.max()
    counts /= counts.sum()
    observed = counts > 0
    p, M = counts[observed], mu[observed]
    K = M.shape[1]
    w = np.full(K, 1.0 / K)
    eps = np.finfo(np.float64).eps
    operations = len(p) + K + 16
    gamma = operations * eps / max(1 - operations * eps, eps)

    for iteration in range(int(iters) + 1):
        mix = M @ w
        if not (mix > 0).all() or not np.isfinite(mix).all():
            raise FloatingPointError('Mixture optimizer lost positive observed-count support')
        gradient = (M / mix[:, None]).T @ p
        multiplier = float(w @ gradient)
        gap = max(0.0, float(gradient.max() - multiplier))
        objective = float(p @ np.log(mix))
        budget = gamma * (max(1., float(p @ np.abs(np.log(mix)))) + max(1., float(gradient.max())))
        if gap <= budget:
            # Undo only the accepted input-column rounding normalization, so
            # callers using their original mu recover the fitted composition.
            answer = w / columns
            answer /= answer.sum()
            diagnostic = dict(converged=True, iterations=iteration, normalized_log_likelihood=objective,
                              normalized_objective_gap_upper=gap, float64_roundoff_budget=budget,
                              kkt_complementarity=float(np.max(np.abs(w * (gradient - multiplier)))),
                              reference_columns=K,
                              interpretation='conditional Poisson composition MLE; weights need not be identifiable')
            return (answer, diagnostic) if return_diagnostics else answer
        if iteration == iters:
            break
        best = int(np.argmax(gradient))
        active = np.flatnonzero((w > 0) | (np.arange(K) == best))
        A = M[:, active]
        # Centre before forming curvature: subtracting nearly equal entries of
        # A.T@diag(p/mix^2)@A would erase information for similar references.
        D = A - A.mean(axis=1, keepdims=True)
        relative = D / mix[:, None]
        H = relative.T @ (p[:, None] * relative)
        reduced_gradient = relative.T @ p
        step = np.linalg.lstsq(H, reduced_gradient, rcond=None)[0]
        step -= step.mean()
        direction = np.zeros(K)
        direction[active] = step
        negative = direction < 0
        upper = min(1., float(np.min(-w[negative] / direction[negative]))) if negative.any() else 0.
        if (not np.isfinite(direction).all() or float(gradient @ direction) <= 0 or upper <= 0):
            # Guaranteed feasible ascent when the Newton active set is singular
            # or a zero-weight boundary blocks its direction.
            direction = -w.copy()
            direction[best] += 1.
            upper = 1.
        change = M @ direction

        def derivative(distance):
            with np.errstate(divide='ignore', invalid='ignore'):
                candidate_mix = M @ np.maximum(w + distance * direction, 0.)
                return float(np.sum(p * change / candidate_mix))

        if derivative(0.) <= 0:
            raise FloatingPointError('Mixture optimum cannot be certified at floating-point resolution')
        if derivative(upper) >= 0:
            distance = upper
        else:
            distance = brentq(derivative, 0., upper, xtol=np.nextafter(0., 1.), rtol=4 * eps, maxiter=int(iters))
        boundary = np.zeros(K, bool)
        if distance == upper:
            decreasing = direction < 0
            boundary[decreasing] = (-w[decreasing] / direction[decreasing]) == distance
        w = np.maximum(w + distance * direction, 0.)
        w[boundary] = 0.
        w /= w.sum()
        new_objective = float(p @ np.log(M @ w))
        if not np.isfinite(new_objective) or new_objective < objective - budget:
            raise FloatingPointError('Mixture line search decreased the objective beyond its roundoff budget')
    raise RuntimeError(f'Mixture optimizer did not converge: objective gap {gap:.6g}, roundoff budget {budget:.6g}')


def depth_design(library, df=4):
    """Cubic truncated-power spline in log10 depth (as CNV_infercnv.depth_design)."""
    x = np.log10(np.maximum(library, 1.0))
    knots = np.quantile(x, np.linspace(0, 1, df + 1)[1:-1])
    return np.column_stack([np.ones_like(x), x] + [np.maximum(x - k, 0) ** 3 for k in knots])


def proposal_features(Y, L, rho, reset, window, mask):
    idx = np.flatnonzero(mask)
    z = np.log((Y[:, idx] + 0.5) / (L[:, None] * rho[idx] + 0.5))
    chrom = np.cumsum(_subset_reset(reset, idx))
    out = np.empty_like(z)
    half = max(1, window // 2)
    for c in np.unique(chrom):
        cols = np.flatnonzero(chrom == c)
        cs = np.concatenate([np.zeros((len(z), 1)), np.cumsum(z[:, cols], axis=1)], axis=1)
        lo = np.clip(np.arange(len(cols)) - half, 0, len(cols))
        hi = np.clip(np.arange(len(cols)) + half + 1, 0, len(cols))
        out[:, cols] = (cs[:, hi] - cs[:, lo]) / (hi - lo)
    D = depth_design(L)
    coef, *_ = np.linalg.lstsq(D, out, rcond=None)
    return (out - D @ coef).astype(np.float32)


def neighbourhood_sum(X, graph, idx):
    """Sum of X over each cell and its graph neighbours within idx (proposals only), and the counts summed."""
    G = graph[idx][:, idx]
    G = (G > 0).astype(np.float32)
    return np.asarray(X + G @ X), 1.0 + np.asarray(G.sum(axis=1)).ravel()


def smooth_features(F, graph, idx):
    """Mean of each cell's proposal features over itself and its neighbours (within idx). Proposals only;
    cell likelihoods always use the cell's own counts."""
    if graph is None:
        return F
    total, n = neighbourhood_sum(F, graph, idx)
    return (total / n[:, None]).astype(np.float32)


def leiden_groups(F, params: CloneParams):
    """PCA, exact kNN graph (symmetrized, unweighted) and igraph Leiden (modularity) on proposal features."""
    import igraph as ig
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    n = F.shape[0]
    if n < 2 * params.min_cells:
        return np.zeros(n, dtype=int)
    X = PCA(n_components=int(min(params.n_pcs, n - 1, F.shape[1] - 1)), svd_solver='randomized',
            random_state=params.seed).fit_transform(np.asarray(F, np.float64))
    k = int(min(params.n_neighbors, n - 1))
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    # Duplicate feature vectors need not return the query itself first. Remove
    # self by row identity so every vertex retains k nonself neighbours.
    idx = np.array([neighbours[neighbours != i][:k] for i, neighbours in enumerate(idx)])
    edges = np.column_stack([np.repeat(np.arange(n), k), idx.ravel()])
    g = ig.Graph(n=n, edges=edges, directed=False).simplify()
    ig.set_random_number_generator(random.Random(params.seed))
    part = g.community_leiden(objective_function='modularity', resolution=params.leiden_resolution, n_iterations=2)
    return consolidate(np.asarray(part.membership), X, params)


def consolidate(g, X, params: CloneParams):
    """Fold groups below min_cells (and beyond max_groups) into the nearest centroid."""
    g = np.unique(g, return_inverse=True)[1]
    while True:
        sizes = np.bincount(g)
        live = np.flatnonzero(sizes)
        small = [k for k in live if sizes[k] < params.min_cells]
        if len(live) > params.max_groups:
            small = [live[np.argmin(sizes[live])]]
        if not small or len(live) == 1:
            break
        k = min(small, key=lambda x: sizes[x])
        cent = np.array([X[g == j].mean(axis=0) for j in live])
        dist = np.linalg.norm(cent - X[g == k].mean(axis=0), axis=1)
        dist[live == k] = np.inf
        g[g == k] = live[np.argmin(dist)]
    return np.unique(g, return_inverse=True)[1]


def pooled_depth_slopes(rows: Rows, L, x, mapc, min_cells, n_q=5):
    """gamma_j from depth quintiles within each group of mapc, pooled over groups (not centred)."""
    q = np.zeros(rows.n, dtype=int)
    for k in np.unique(mapc):
        idx = np.flatnonzero(mapc == k)
        if len(idx) >= 2 * min_cells:
            edges = np.quantile(x[idx], np.linspace(0, 1, n_q + 1)[1:-1])
            q[idx] = np.searchsorted(edges, x[idx])
    group = mapc * n_q + q
    G = int(group.max()) + 1
    onehot = sparse.csr_matrix((np.ones(rows.n, np.float32), (np.arange(rows.n), group)), shape=(rows.n, G))
    Yg = rows.weighted_sum(onehot)
    Lg = np.bincount(group, weights=L, minlength=G)
    xg = np.bincount(group, weights=L * x, minlength=G) / np.maximum(Lg, 1.0)
    num, den = np.zeros(Yg.shape[1]), 0.0
    for k in np.unique(mapc):
        gs = [g for g in range(k * n_q, (k + 1) * n_q) if g < len(Lg) and Lg[g] > 0]
        if len(gs) < 3:
            continue
        share = np.log((Yg[gs] + 0.5) / Lg[gs, None])
        share -= share.mean(axis=0)
        xc = xg[gs] - xg[gs].mean()
        num += xc @ share
        den += xc @ xc
    if den <= 0:
        return np.zeros(Yg.shape[1])
    return np.clip(num / den, -1.5, 1.5)


def reference_depth_slopes(Y, L, groups, min_cells, chunk=4096):
    """Depth slopes from reference cells (list of row arrays, one per reference group). Reference cells
    carry no tumour copy number, so the slope cannot absorb it when depth tracks tumour content, as it
    does for Visium HD bins."""
    rows = np.concatenate(groups)
    mapc = np.repeat(np.arange(len(groups)), [len(g) for g in groups])
    logL = np.log(np.maximum(np.asarray(L, np.float64)[rows], 1.0))
    return pooled_depth_slopes(Rows(Y, rows, chunk), np.asarray(L, np.float64)[rows], logL - logL.mean(), mapc,
                               min_cells)


def _thin_diff_params(pp: ProfileParams):
    """Difference-test settings for split_mode 'thin' (tests on every bin) with the same genomic meaning as on
    odd bins: a run of min_bins odd bins spans 2 min_bins - 1 bins, and a change probability per odd-bin step
    is about half that per bin step."""
    from dataclasses import replace
    return replace(pp, min_bins=2 * pp.min_bins - 1, p_change=pp.p_change / 2)


class SpecimenModel:
    def __init__(self, rows: Rows, L, mu, reset, phi, params: CloneParams, programmes=None, alleles=None,
                 graph=None, gamma=None):
        if params.split_null_draws < 0:
            raise ValueError('split_null_draws must be nonnegative')
        if alleles is not None and params.split_null_draws > 0:
            raise ValueError('Alleles with split_null_draws > 0 are unsupported: the expression-only null '
                             'does not reproduce allele proposals or positive allele evidence. Disable the '
                             'split null or run expression-only; this filter is not a joint calibrated test.')
        self.rows, self.L, self.reset, self.params = rows, np.asarray(L, np.float64), reset, params
        self.pp = params.profile
        self.J = len(reset)
        self.Ytot = rows.weighted_sum(np.ones((rows.n, 1), np.float32))[0]
        self.weights = mixture_weights(self.Ytot, mu)
        base = mu @ self.weights
        self.base = base / base.sum()
        # Per-bin baseline mismatch from reference heterogeneity: large where the fitted mixture
        # blends cell types that disagree (no matched reference), ~0 when one type dominates.
        self.v = (self.weights[None] * (np.log(mu / mu.sum(axis=0, keepdims=True)) - np.log(self.base)[:, None]) ** 2
                  ).sum(axis=1)
        self.centre, self.beta, self.tau = fit_pooled(self.Ytot, self.L.sum() * self.base, reset, self.pp, v=self.v)
        self.phi = np.asarray(phi, np.float64).copy()
        self.offset = np.zeros(self.J)
        logL = np.log(np.maximum(self.L, 1.0))
        self.x = logL - logL.mean()
        self.gamma = np.zeros(self.J)
        # A fixed gamma (from reference cells) is centred once and never refit from query cells.
        self.fixed_gamma = gamma is not None
        self.gamma = (np.asarray(gamma, np.float64) - np.average(gamma, weights=self.base) if self.fixed_gamma
                      else self.depth_slopes(np.ones((rows.n, 1), np.float32)))
        self.programmes = programmes
        self.W = np.zeros((1, self.J)) if programmes is None else programmes.W
        self.U = np.zeros((rows.n, len(self.W)), np.float32)
        if params.split_mode not in SPLIT_MODES:
            raise ValueError(f'split_mode must be one of {SPLIT_MODES}, got {params.split_mode!r}')
        self.thin = params.split_mode == 'thin'
        if self.thin:
            # Every bin is both a proposal and a test bin; the counts are split instead (Rows folds).
            self.propose_mask = np.ones(self.J, dtype=bool)
            self.test_mask = np.ones(self.J, dtype=bool)
            self.prop_fold, self.test_fold = 'prop', 'test'
            self.fold_scale = {'prop': THIN_EPS, 'test': 1.0 - THIN_EPS}
            self.thin_size = 1.0 / np.maximum(self.phi, 1e-3)        # frozen: draws stay fixed while phi is refit
            rows.set_thinning(params.seed, self.thin_size, THIN_EPS)
            self.diff_pp = _thin_diff_params(self.pp)
            self.window_bins = 2 * params.window_bins                 # same genomic span as on even bins
        else:
            self.propose_mask = np.arange(self.J) % 2 == 0
            self.test_mask = ~self.propose_mask
            self.prop_fold = self.test_fold = None
            self.fold_scale = {}
            self.diff_pp = self.pp
            self.window_bins = params.window_bins
        self.log = []
        self.alleles = alleles
        self.Ak = self.Dk = self.B = self.allele_post = None
        self.graph, self.lam = graph, 0.0
        self.genes = None               # (subsample rows, gene counts CSR [sub, genes], gene bin) for gene checks
        self.R_prev = None              # mean-field state carried between E-steps (spatial only)

    # Directed spatial neighbourhood smoothing --------------------------------
    def posterior(self, LL, pi, iters=10):
        logpi = np.log(pi)
        R = posteriors(LL, logpi)
        if self.graph is None or LL.shape[1] < 2:
            return R
        # lambda from the data-only labels (smoothed labels would always look coherent).
        self.lam = self.params.potts_lambda if self.params.potts_lambda >= 0 else self.choose_lambda(LL, R, logpi)
        if self.lam <= 0:
            self.R_prev = R
            return R
        if self.R_prev is not None and self.R_prev.shape == R.shape:
            R = self.R_prev               # Warm start the finite smoothing iterations.
        for _ in range(iters):
            R = posteriors(LL + self.lam * np.asarray(self.graph @ R, np.float32), logpi)
        self.R_prev = R
        return R

    def choose_lambda(self, LL, R0, logpi, grid=(0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)):
        """lambda maximizing how well neighbours predict each cell's own data.

        score = sum_c log sum_k softmax(log pi + lambda (G R0)_c)_k exp(LL_ck), with R0
        the data-only responsibilities. The graph excludes self-edges. The
        selected strength is a predictive heuristic and can be positive by chance."""
        GR = np.asarray(self.graph @ R0)
        ll = LL - LL.max(axis=1, keepdims=True)
        best, lam = -np.inf, 0.0
        for g in grid:
            z = logpi[None] + g * GR
            z = z - logsumexp(z, axis=1, keepdims=True)
            score = float(logsumexp(z + ll, axis=1).sum())
            if score > best + 1e-6:
                best, lam = score, g
        return lam

    # Alleles ----------------------------------------------------------------
    def set_allele_pools(self, Ak, Dk):
        """Clone-level allele state for pooled counts aligned with the current doses (row 0 normal-like)."""
        if self.alleles is None:
            return
        self.Ak, self.Dk = np.asarray(Ak, np.float64), np.asarray(Dk, np.float64)
        post, B = self.alleles.fit(self.Ak, self.Dk)
        B[0] = 0.5
        self.allele_post, self.B = post, B

    def allele_events(self, k):
        return (self.alleles is not None and self.allele_post is not None
                and self.alleles.has_events(self.allele_post[k], self.Dk[k]))

    def _keep_alleles(self, keep):
        if self.B is not None:
            self.Ak, self.Dk, self.B, self.allele_post = (x[keep] for x in (self.Ak, self.Dk, self.B, self.allele_post))

    # Profiles -------------------------------------------------------------
    def depth_slopes(self, R, n_q=5):
        """gamma_j from depth quintiles within each clone (MAP), pooled over clones."""
        gamma = pooled_depth_slopes(self.rows, self.L, self.x, R.argmax(axis=1), self.params.min_cells, n_q)
        return gamma - np.average(gamma, weights=self.base)

    def depth_log_rate(self, Lk, XLk):
        """[G, J] log-rate shift of each pooled group from its cells' depths (up to a per-group constant).

        The pooled expectation is sum_c L_c softmax(log base + gamma x_c). With log Z(x) ~ c x^2 / 2
        (c = sum_j base_j gamma_j^2; gamma is base-centred) and x Gaussian under UMI weighting with mean m and
        variance v, the average is base_j exp((gamma_j m + gamma_j^2 v / 2) / (1 + c v)) up to a constant.
        On simulated depths (sd 0.6, slopes up to 1.5) the worst bin error was 2%, against 24% for the
        mean depth alone (tests/test_numerics.py)."""
        xbar, var = _depth_mean_var(Lk, XLk)
        c = float(np.sum(self.base * self.gamma ** 2))
        k = 1.0 / (1.0 + c * var)
        return k[:, None] * (self.gamma[None] * xbar[:, None] + 0.5 * self.gamma[None] ** 2 * var[:, None])

    def expected(self, Lk, XLk):
        rate = self.base[None] * np.exp(self.depth_log_rate(Lk, XLk))
        return np.atleast_1d(Lk)[:, None] * rate / rate.sum(axis=1, keepdims=True)

    def shared_profiles(self, Yk, Lk, XLk, labels=None):
        """profiles() for every group, with tumour clone doses (groups 1..) replaced by one pooled tumour dose
        plus, on bins where the clone differs from another tumour clone, its measured mean difference from the
        pool.

        profiles() segments each clone on its own, and segmentation shrinks weak evidence to neutral, so a
        small clone's genome-wide dose is noisier and more shrunk than a large clone's. Cells then fit the
        large clone better on hundreds of bins where the clones do not differ, and EM drains the small clone
        (MCF7 with a 15-bin implanted loss in 57/143 cells: the region doses stayed -0.33 vs +0.22 while the
        implanted cells' log-likelihood preference moved from -15 to +2 nats in seven rounds). Here clones
        share every bin outside resolved pairwise differences (`differences`, test bins), and inside them take
        a continuous, run-averaged difference: the step between a clone and a mixed pool is fractional, and
        the segmentation states are whole copy-number steps (anchoring each clone's states to the pool, or
        segmenting against the pooled dose, failed on the same MCF7 case). Adding the deviation on every run,
        not only on resolved differences, rescued a whole-chr12 subclone in test_spatial on one node but cost
        0.02 ARI on Kinker v3 and failed test_recovers_clones_and_normals (dev/hd1_ablate_20261101_2350/).
        Returns profiles() output with 'dose' replaced for groups 1.. and the unanchored doses in 'dose_own'
        (for has_events).

        The deviation is averaged within atomic intervals: the bins are cut at every accepted segment boundary
        of every pair and wherever the pooled path changes, so a third clone's long event cannot join two
        opposite events that separate the first two (review 2026-11-02 #1: a 0.099 contrast fell to 0.005).
        Segments pass the gene-consistency veto when MAP labels are given (alleles are not used here: the
        allele pools are updated after this step), and the deviation has the fitted programme contrast
        removed, as in `differences` (#2)."""
        Yk, Lk, XLk = np.atleast_2d(Yk), np.atleast_1d(Lk), np.atleast_2d(XLk)
        prof = self.profiles(Yk, Lk, XLk)
        t = np.arange(1, len(Yk))
        if len(t) < 2:                                     # one tumour clone: it is the shared profile
            return prof
        Yt, Lt, XLt = Yk[t].sum(axis=0), Lk[t].sum(), XLk[t].sum(axis=0)
        pooled = self.profiles(Yt, Lt, XLt)
        shared = pooled['dose'][0]
        pairs = [(a, b) for i, a in enumerate(t) for b in t[i + 1:]]
        A, B = [a for a, _ in pairs], [b for _, b in pairs]
        res = self.differences(Yk[A], Lk[A], XLk[A], Yk[B], Lk[B], XLk[B])
        mark = np.zeros((len(Yk), self.J), bool)
        cut = self.reset.copy()
        cut[1:] |= pooled['path'][0][1:] != pooled['path'][0][:-1]
        for (a, b), r in zip(pairs, res):
            segs = r['segments'] if labels is None else self.supported(r['segments'], None, None, None, None,
                                                                       labels=labels, a=a, b=b, log=False)
            for seg in segs:
                f, l = seg['first_bin'], seg['last_bin']
                mark[[a, b], f:l + 1] = True
                cut[f] = True
                if l + 1 < self.J:
                    cut[l + 1] = True
        atom = np.cumsum(cut)
        Et = self.expected(Lt, XLt)[0]
        E = self.expected(Lk, XLk)
        dose = prof['dose'].copy()
        for k in t:
            d = np.log((Yk[k] + 0.5) / E[k]) - np.log((Yt + 0.5) / Et)
            w = 1.0 / (1.0 / (Yk[k] + 0.5) + 1.0 / (Yt + 0.5))
            d = d - weighted_median(d, w)
            if self.programmes is not None:
                d = d - fit_amplitudes(d - np.median(d), 1.0 / (0.04 + 1.0 / w), self.W)[0]
            dose[k] = shared
            for a, b, key in runs(2 * atom + mark[k], self.reset):
                if key % 2:
                    dose[k, a:b + 1] = shared[a:b + 1] + np.sum(w[a:b + 1] * d[a:b + 1]) / np.sum(w[a:b + 1])
        prof['dose_own'] = prof['dose']
        prof['dose'] = dose
        return prof

    def profiles(self, Yk, Lk, XLk):
        Yk, Lk = np.atleast_2d(Yk), np.atleast_1d(Lk)
        E = self.expected(Lk, XLk)
        tau = np.maximum(group_tau(Yk, E, np.full(len(Yk), self.centre), self.reset, self.pp.tau_floor),
                         self.tau)
        centre = group_centres(Yk, E, tau)
        post, path, dose, _ = segment(Yk, E, self.reset, centre, self.beta, self.bin_tau(tau), self.pp,
                                      decode=self.programmes is None)
        shift = np.zeros_like(E)
        if self.programmes is not None:
            # Programme amplitudes fitted around the current dose, then re-segmented (twice).
            E0 = E
            for i in range(2):
                w = 1.0 / (tau[:, None] ** 2 + 1.0 / (Yk + 0.5))
                shift = fit_amplitudes(log_ratio(Yk, E0) - centre[:, None] - dose, w, self.W)
                E = E0 * np.exp(shift)
                E *= (E0.sum(axis=1) / E.sum(axis=1))[:, None]
                tau = np.maximum(group_tau(Yk, E, centre, self.reset, self.pp.tau_floor), self.tau)
                centre = group_centres(Yk, E, tau)
                post, path, dose, _ = segment(Yk, E, self.reset, centre, self.beta, self.bin_tau(tau), self.pp,
                                              decode=i == 1)
        return dict(post=post, path=path, dose=dose, tau=tau, centre=centre, E=E,
                    programme_shift=np.sqrt((shift ** 2).mean(axis=1)))

    def bin_tau(self, tau):
        """[G, J] mismatch sd: the group's tau combined with the per-bin reference heterogeneity."""
        return np.sqrt(np.asarray(tau)[:, None] ** 2 + self.v[None])

    def has_events(self, Y, E, post, path, dose, tau):
        rows = segments_table(post, path, dose, Y, E, self.reset, self.beta, tau, self.pp, reference_var=self.v)
        return any(r['resolved'] for r in rows)

    def log_rates(self, doses):
        lr = np.log(self.base)[None] + self.offset[None] + np.asarray(doses)
        return lr - logsumexp(lr, axis=1, keepdims=True)

    # Group merging --------------------------------------------------------
    def differences(self, Ya, La, XLa, Yb, Lb, XLb, *, count_variance=None, weight_scale=None):
        shift = self.depth_log_rate(La, XLa) - self.depth_log_rate(Lb, XLb)
        if self.programmes is not None:
            Ya, Yb = np.atleast_2d(Ya), np.atleast_2d(Yb)
            Qa, Qb = (None, None) if count_variance is None else count_variance
            ca, cb = (None, None) if weight_scale is None else weight_scale
            ra, va = pooled_lograte(Ya, La, Qa, ca)
            rb, vb = pooled_lograte(Yb, Lb, Qb, cb)
            d = ra - rb - shift
            v = va + vb
            shift = shift + fit_amplitudes(d - np.median(d, axis=1, keepdims=True), 1.0 / (0.04 + v), self.W)
        return difference_segments(Ya, La, Yb, Lb, self.reset, self.beta, self.diff_pp, mask=self.test_mask,
                                   shift=shift, count_variance=count_variance, weight_scale=weight_scale)

    def matched_differences(self, labels, pairs, n_strata=8, *, cell_weights=None,
                            min_effective_cells=None, report_depth_balance=False):
        """Differences between groups of cells with depth-matched pooling.

        labels [n]: group per cell (-1: none). For each pair, cells are weighted so
        both groups have the same weight mass in each shared log-depth stratum.
        Depth and RNA composition can still differ within a stratum. A comparison
        without shared stratum mass is untestable.

        Optional nonnegative cell_weights are matched by weight mass, not cell
        count. Squared weights determine the conditional sampling term. These
        calculations condition on the supplied weights; fitted membership weights
        do not acquire an independent-confirmation or calibration claim.
        With report_depth_balance, return the weighted empirical Wasserstein
        distance of the matched log depths for each pair, in natural-log units.
        """
        if cell_weights is not None:
            cell_weights = np.asarray(cell_weights, np.float64)
            if (cell_weights.shape != (self.rows.n,) or not np.isfinite(cell_weights).all()
                    or np.any(cell_weights < 0)):
                raise ValueError('cell_weights must be finite, nonnegative and aligned to specimen rows')
        if min_effective_cells is not None and (not np.isfinite(min_effective_cells) or min_effective_cells < 0):
            raise ValueError('min_effective_cells must be finite and nonnegative')
        # Sparse cells x 2P weights: each cell is in at most the pairs that include its group (a dense
        # matrix was cells x 2P float32 plus a float64 copy, ~17 GiB at 1M cells and 780 pairs).
        r_idx, c_idx, val = [], [], []
        for i, (a, b) in enumerate(pairs):
            ia, ib = np.flatnonzero(labels == a), np.flatnonzero(labels == b)
            edges = np.quantile(np.r_[self.x[ia], self.x[ib]], np.linspace(0, 1, n_strata + 1)[1:-1])
            sa, sb = np.searchsorted(edges, self.x[ia]), np.searchsorted(edges, self.x[ib])
            wa = np.ones(len(ia)) if cell_weights is None else cell_weights[ia]
            wb = np.ones(len(ib)) if cell_weights is None else cell_weights[ib]
            na = np.bincount(sa, weights=wa, minlength=n_strata)
            nb = np.bincount(sb, weights=wb, minlength=n_strata)
            m = np.minimum(na, nb)
            r_idx += [ia, ib]
            c_idx += [np.full(len(ia), 2 * i), np.full(len(ib), 2 * i + 1)]
            val += [wa * np.divide(m, na, out=np.zeros_like(m), where=na > 0)[sa],
                    wb * np.divide(m, nb, out=np.zeros_like(m), where=nb > 0)[sb]]
        W = sparse.csr_matrix((np.concatenate(val), (np.concatenate(r_idx), np.concatenate(c_idx))),
                              shape=(self.rows.n, 2 * len(pairs)), dtype=np.float64)
        Y, L, XL = self.pooled(W, self.test_fold)
        W2 = W.multiply(W)
        Q = self.rows.weighted_sum(W2, self.test_fold)
        total_weight = np.asarray(W.sum(axis=0)).ravel()
        c = np.divide(np.asarray(W2.sum(axis=0)).ravel(), total_weight,
                      out=np.zeros_like(total_weight), where=total_weight > 0)
        result = self.differences(Y[0::2], L[0::2], XL[0::2], Y[1::2], L[1::2], XL[1::2],
                                  count_variance=(Q[0::2], Q[1::2]), weight_scale=(c[0::2], c[1::2]))
        if report_depth_balance:
            from scipy.stats import wasserstein_distance

            for i, r in enumerate(result):
                ia, ib = r_idx[2 * i:2 * i + 2]
                wa, wb = val[2 * i:2 * i + 2]
                shared = bool(np.all(total_weight[2 * i:2 * i + 2] > 0))
                r['weighted_log_depth_wasserstein'] = (
                    float(wasserstein_distance(self.x[ia], self.x[ib], wa, wb)) if shared else None)
                if not r['testable']:
                    r['untestable_reason'] = 'unusable_test_contrast' if shared else 'no_shared_depth_stratum_mass'
        if min_effective_cells is not None:
            Wc = W.tocsc()
            effective = []
            for j in range(W.shape[1]):
                w = Wc.data[Wc.indptr[j]:Wc.indptr[j + 1]]
                w = w[w > 0].astype(np.longdouble)
                # Constant weights have exactly the number of positive cells;
                # avoid a rounded ratio just below an integer support boundary.
                effective.append(0. if not len(w) else float(len(w)) if np.all(w == w[0])
                                 else float(w.sum() ** 2 / (w * w).sum()))
            for i, r in enumerate(result):
                r['effective_cells'] = effective[2 * i:2 * i + 2]
                if min(r['effective_cells']) < min_effective_cells:
                    r.update(segments=[], distance=float('inf'), testable=False,
                             untestable_reason='insufficient_weighted_cell_support')
        return result

    def supported(self, segments, Aa, Da, Ab, Db, labels=None, a=None, b=None, log=True):
        """Difference segments that pass the gene-consistency test and are not vetoed by alleles."""
        keep = []
        for s in segments:
            where = f"bins{s['first_bin']}-{s['last_bin']}"
            if labels is not None:
                gene_support = self.gene_consistent(s, labels, a, b)
                if gene_support is False:
                    if log:
                        self.log.append(dict(action='gene_veto', context=where))
                    continue
                if gene_support is None and self.genes is not None and log:
                    self.log.append(dict(action='gene_veto_uninformative', context=where))
            if (self.alleles is not None and getattr(self.alleles, 'bins', None) is not None and Aa is not None
                    and self.alleles.veto(self.alleles.snps_in(self.alleles.bins, s['first_bin'], s['last_bin']),
                                          Aa, Da, Ab, Db, log_dose=s['mean_diff'] / self.beta)):
                if log:
                    self.log.append(dict(action='allele_veto', context=where))
                continue
            keep.append(s)
        return keep

    def gene_consistent(self, seg, labels, a, b, min_genes=8, alpha=0.01):
        """Return concordance, veto, or None when the gene contrast is unavailable.

        This heuristic uses unweighted hard assignments on a proposal subsample,
        which can differ from the depth-matched population. The selected interval
        and correlated gene signs limit the nominal sign test. Coordinated cell
        state and variable dosage response both affect concordance. In thinning
        mode the gene counts come from the test half (`thin_gene_counts`).
        """
        if self.genes is None:
            return None
        sub, Xg, gbin = self.genes
        la = labels[sub]
        ia, ib = la == a, la == b
        if ia.sum() < 10 or ib.sum() < 10:
            return None
        ya, yb = np.asarray(Xg[ia].sum(axis=0)).ravel(), np.asarray(Xg[ib].sum(axis=0)).ravel()
        total_a, total_b = ya.sum(), yb.sum()
        if total_a == 0 or total_b == 0:
            return None
        r = np.log((ya + 0.5) / total_a) - np.log((yb + 0.5) / total_b)
        informative = (ya + yb) >= 10
        if informative.sum() < min_genes:
            return None
        r = r - np.median(r[informative])
        sel = (informative & (gbin >= seg['first_bin']) & (gbin <= seg['last_bin'])
               & self.test_mask[np.clip(gbin, 0, self.J - 1)])
        n = int(sel.sum())
        if n < min_genes:
            return None
        k = int((np.sign(r[sel]) == (1 if seg['direction'] == 'higher' else -1)).sum())
        from scipy.stats import binom
        return bool(binom.sf(k - 1, n, 0.5) < alpha)

    def merge_unsupported(self, R, doses, pi):
        """Merge tumour clones with no depth-matched dosage difference and no allele difference."""
        merged = False
        statistical_merge = False
        # Exact duplicate RNA profiles have identical expression likelihoods. Hard MAP
        # assignments can hide all but one of them from the pairwise tests below.
        # Preserve the normal component and any allele or spatial distinctions.
        if self.alleles is None and self.graph is None:
            keep = [0]
            duplicate = []
            for k in range(1, len(doses)):
                target = next((i for i in range(1, len(keep))
                               if np.all(np.isfinite(doses[k]))
                               and np.array_equal(doses[k], doses[keep[i]])), None)
                if target is None:
                    keep.append(k)
                else:
                    duplicate.append((k, target))
            if duplicate:
                original_R, original_pi = R, pi
                R, pi = R[:, keep], pi[keep]
                for k, target in duplicate:
                    R[:, target] += original_R[:, k]
                    pi[target] += original_pi[k]
                    self.log.append(dict(action='merge_identical_profiles', a=keep[target], b=k,
                                         criterion_scope='exact_expression_profile',
                                         interpretation='RNA_indistinguishable_not_genomic_equivalence'))
                doses = [doses[k] for k in keep]
                merged = True
        while len(doses) > 2:
            assignment_R = R
            labels = R.argmax(axis=1)
            if self.thin:
                # Groups judged on the test half are formed from the proposal half only, so the counts that
                # test a pair did not assign its cells (review 2026-11-02 #3). Data-only posteriors: no Potts
                # state is touched. EM itself keeps the full counts.
                assignment_R = posteriors(self.loglik(doses, self.prop_fold), np.log(pi))
                labels = assignment_R.argmax(axis=1)
            tum = [k for k in range(1, len(doses)) if (labels == k).sum() >= self.params.min_cells]
            pairs = [(a, b) for i, a in enumerate(tum) for b in tum[i + 1:]]
            if not pairs:
                break
            res = self.matched_differences(labels, pairs)
            cand = []
            for (a, b), r in zip(pairs, res):
                if not r.get('testable', True):
                    self.log.append(dict(action='untestable_pair', a=a, b=b, context='no depth overlap'))
                ad = (None,) * 4 if self.Ak is None else (self.Ak[a], self.Dk[a], self.Ak[b], self.Dk[b])
                if self.supported(r['segments'], *ad, labels=labels, a=a, b=b):
                    continue
                if self.allele_post is not None:
                    ra, rb = (self.alleles.runs(self.allele_post[k], self.Dk[k]) for k in (a, b))
                    if self.alleles.differ(ra, self.Ak[a], self.Dk[a], rb, self.Ak[b], self.Dk[b]):
                        continue
                cand.append((r['distance'], a, b))
            if not cand:
                break
            _, a, b = min(cand)
            pair_cells = (labels == a) | (labels == b)
            pair_mass = assignment_R[pair_cells, a] + assignment_R[pair_cells, b]
            pair_probability = np.divide(assignment_R[pair_cells, a], pair_mass,
                                         out=np.full(pair_mass.shape, .5), where=pair_mass > 0)
            from scipy.special import xlogy
            entropy = -xlogy(pair_probability, pair_probability) - xlogy(1 - pair_probability, 1 - pair_probability)
            self.log.append(dict(action='merge_matched', a=a, b=b, n_a=int((labels == a).sum()),
                                 n_b=int((labels == b).sum()),
                                 mean_pair_entropy_nats=float(entropy.mean()),
                                 criterion_scope='assigned_cell_pools',
                                 interpretation='no_resolved_difference_not_genomic_equivalence'))
            keep = [k for k in range(len(doses)) if k != b]
            R = R.copy()
            R[:, a] += R[:, b]
            R = R[:, keep]
            pi = pi.copy()
            pi[a] += pi[b]
            pi = pi[keep]
            if self.Ak is not None:
                self.Ak[a] += self.Ak[b]
                self.Dk[a] += self.Dk[b]
                self.set_allele_pools(self.Ak[keep], self.Dk[keep])
            doses = [doses[k] for k in keep]
            merged = True
            statistical_merge = True
        if statistical_merge:
            Yk, Lk, XLk = self.pooled(R)
            prof = self.profiles(Yk, Lk, XLk)
            doses = [np.zeros(self.J)] + [prof['dose'][k] for k in range(1, len(doses))]
            if self.graph is not None:
                self.R_prev = R
        return doses, pi / pi.sum(), merged

    def agglomerate(self, Yg, Lg, XLg, labels=None, context='initial', AD=None):
        """Merge pooled groups that differ in no resolved segment (dosage or alleles); most similar first.

        AD: optional pooled allele counts (Ag, Dg) [G, S]. Returns members, pooled arrays and AD."""
        members = [[i] for i in range(len(Yg))]
        Yg, Lg, XLg = [y for y in Yg], list(Lg), list(XLg)
        use_alleles = AD is not None and self.alleles is not None
        if use_alleles:
            Ag, Dg = [x for x in AD[0]], [x for x in AD[1]]
            post = list(self.alleles.fit(np.array(Ag), np.array(Dg))[0])
            runs = [self.alleles.runs(q, d) for q, d in zip(post, Dg)]
        # Pair results depend only on the two groups, so they are cached by group id and only pairs involving
        # a newly merged group are tested again (40 -> 1 groups: 1,521 pair tests instead of 10,660).
        ids, next_id, cache = list(range(len(members))), len(members), {}
        while len(members) > 1:
            pairs = [(a, b) for a in range(len(members)) for b in range(a + 1, len(members))]
            todo = [(a, b) for a, b in pairs if (ids[a], ids[b]) not in cache]
            if todo:
                A, B = [a for a, _ in todo], [b for _, b in todo]
                res = self.differences(np.array([Yg[a] for a in A]), np.array([Lg[a] for a in A]),
                                       np.array([XLg[a] for a in A]), np.array([Yg[b] for b in B]),
                                       np.array([Lg[b] for b in B]), np.array([XLg[b] for b in B]))
                for r, (a, b) in zip(res, todo):
                    differ = bool(r['segments']) or (use_alleles and self.alleles.differ(
                        runs[a], Ag[a], Dg[a], runs[b], Ag[b], Dg[b]))
                    cache[(ids[a], ids[b])] = (r['distance'], differ)
            same = [(cache[(ids[a], ids[b])][0], (a, b)) for a, b in pairs if not cache[(ids[a], ids[b])][1]]
            if not same:
                break
            dist, (a, b) = min(same)
            ids[a] = next_id
            ids.pop(b)
            next_id += 1
            self.log.append(dict(action='merge', context=context, a=_name(members[a], labels), b=_name(members[b], labels),
                                 umi_a=float(Lg[a]), umi_b=float(Lg[b]), distance=dist))
            members[a] += members.pop(b)
            Yg[a] = Yg[a] + Yg.pop(b)
            Lg[a] = Lg[a] + Lg.pop(b)
            XLg[a] = XLg[a] + XLg.pop(b)
            if use_alleles:
                Ag[a] = Ag[a] + Ag.pop(b)
                Dg[a] = Dg[a] + Dg.pop(b)
                post.pop(b), runs.pop(b)
                post[a] = self.alleles.fit(Ag[a][None], Dg[a][None])[0][0]
                runs[a] = self.alleles.runs(post[a], Dg[a])
        AD = (np.array(Ag), np.array(Dg)) if use_alleles else None
        return members, np.array(Yg), np.array(Lg), np.array(XLg), AD

    # EM ---------------------------------------------------------------------
    def loglik(self, doses, fold=None):
        """[n, K] cell log-likelihoods. fold: None (full counts, with alleles) or a count-splitting half
        (expression only): a half of an NB(m, size) count is NB(s m, s size), s its share, so L is scaled by s
        and phi (1 / size) by 1 / s."""
        lr = self.log_rates(doses)
        out = np.empty((self.rows.n, len(doses)), np.float32)
        allele = self.B is not None and len(self.B) == len(doses) and fold is None
        s = 1.0 if fold is None else self.fold_scale[fold]
        for a, b, Yc in self.rows.chunks(fold):
            out[a:b] = cell_loglik(Yc, s * self.L[a:b], lr, self.phi / s, self.x[a:b], self.gamma, self.U[a:b],
                                   self.W)
            if allele:
                out[a:b] += self.alleles.cell_loglik(slice(a, b), self.B)
        return out

    def update_scores(self, R, doses):
        """Per-cell programme scores around each cell's responsibility-weighted clone rate."""
        if self.programmes is None:
            return
        lr = self.log_rates(doses)
        for a, b, Yc in self.rows.chunks():
            self.U[a:b] = programme_scores(Yc, self.L[a:b], R[a:b].astype(np.float32), lr, self.phi, self.W,
                                           self.programmes.prior_prec, self.x[a:b], self.gamma)

    def pooled(self, R, fold=None):
        """Pooled counts, UMIs and depth moments per column of R. With a count-splitting fold, counts are that
        half and UMIs and moments are scaled by its expected share."""
        Rt = R.T.astype(np.float64)
        Y = self.rows.weighted_sum(R, fold)
        L, XL = np.asarray(Rt @ self.L).ravel(), np.asarray(Rt @ depth_moments(self.L, self.x))
        if fold is not None:
            L, XL = L * self.fold_scale[fold], XL * self.fold_scale[fold]
        return Y, L, XL

    def em(self, doses, pi, *, stop_on_map=True):
        """Budgeted profile/nuisance refinement; clone 0 keeps dose 0.

        NB assignments, pooled PLN segmentation and residual nuisance updates
        are separate working steps with no established common objective.
        Stable MAP labels do not certify convergence of responsibilities or
        profiles. The historical MAP-stability shortcut is a heuristic, not a
        convergence certificate. ``stop_on_map=False`` runs the configured budget
        for controlled diagnostics. Returns doses, pi, responsibilities, loglik,
        profiles; neither stopping rule certifies convergence. Final scores
        use the last refinement-time allele fit. Its parameters and fit pools
        remain frozen during final scoring, rather than refitting after LL/R
        have been calculated. Final responsibilities need not be a fixed point.
        """
        p = self.params
        prev = None
        termination = 'budgeted'
        completed_rounds = 0
        for iteration in range(p.outer_rounds):
            completed_rounds = iteration + 1
            LL = self.loglik(doses)
            for _ in range(3):
                R = self.posterior(LL, pi)
                pi = np.maximum(R.mean(axis=0), 1e-6)
                pi[0] = max(pi[0], p.normal_prior / 5)   # keep the normal-like option open
                pi /= pi.sum()
            Yk, Lk, XLk = self.pooled(R)
            prof = self.shared_profiles(Yk, Lk, XLk, labels=R.argmax(axis=1))
            if self.alleles is not None:
                self.set_allele_pools(*self.alleles.pool(R))
            new = [np.zeros(self.J)]
            keep = [0]
            for k in range(1, len(doses)):
                if R[:, k].sum() < p.min_cells:
                    self.log.append(dict(action='drop_small', clone=k, n_eff=float(R[:, k].sum())))
                    continue
                if not (self.has_events(Yk[k], prof['E'][k], prof['post'][k], prof['path'][k],
                                        prof.get('dose_own', prof['dose'])[k], prof['tau'][k])
                        or self.allele_events(k)):
                    self.log.append(dict(action='to_normal', clone=k, n_eff=float(R[:, k].sum())))
                    continue
                new.append(prof['dose'][k])
                keep.append(k)
            self._keep_alleles(keep)
            if self.R_prev is not None:
                Rk = self.R_prev[:, keep]
                self.R_prev = Rk / np.maximum(Rk.sum(axis=1, keepdims=True), 1e-12)
            self._update_offset(R, new, keep)
            self._update_phi(R[:, keep], new)
            if not self.fixed_gamma:
                self.gamma = self.depth_slopes(R[:, keep])
            self.update_scores(self._renormalise(R, LL, pi, keep), new)
            changed = len(keep) != len(doses)
            doses, pi = new, pi[keep] / pi[keep].sum()
            mapc = R[:, keep].argmax(axis=1)
            if stop_on_map and not changed and prev is not None and np.mean(mapc != prev) < 0.005:
                termination = 'map_stability'
                break
            prev = None if changed else mapc
        LL = self.loglik(doses)
        R = self.posterior(LL, pi)
        # The last profile/nuisance update can move a component below the same
        # support minimum checked inside refinement. Enforce it on the actual
        # returned state too, retaining clone0 and the fitted parameter snapshot.
        # Removing columns and renormalising is not another refinement round or
        # a convergence claim. Spatial re-scoring can change support again;
        # strictly decreasing component count bounds these checks.
        for _ in range(len(doses)):
            mass = R.sum(axis=0)
            keep = [0] + [k for k in range(1, len(doses)) if mass[k] >= p.min_cells]
            if len(keep) == len(doses):
                break
            for k in range(1, len(doses)):
                if mass[k] < p.min_cells:
                    self.log.append(dict(action='drop_small_final', clone=k, n_eff=float(mass[k]),
                                         minimum=p.min_cells))
            if self.R_prev is not None:
                self.R_prev = self._renormalise(self.R_prev, LL, pi, keep)
            self._keep_alleles(keep)
            doses, pi = [doses[k] for k in keep], pi[keep] / pi[keep].sum()
            LL = LL[:, keep]
            R = self.posterior(LL, pi)
        else:
            raise RuntimeError('Final component support checks did not terminate')
        Yk, Lk, XLk = self.pooled(R)
        self.log.append(dict(action='em', n_eff=' '.join(f'{x:.0f}' for x in R.sum(axis=0)), potts_lambda=self.lam,
                             refinement_rounds=completed_rounds, convergence_status='not_certified_' + termination))
        return doses, pi, R, LL, Yk, Lk, XLk

    @staticmethod
    def _renormalise(R, LL, pi, keep):
        """Responsibilities over the kept clones. A cell whose kept responsibilities all underflowed to 0 (it
        was confidently in a dropped clone) is recomputed in log space from its kept log-likelihoods and log
        mixture weights; a cell with no finite kept log-likelihood falls back to the kept mixture weights."""
        Rk = R[:, keep]
        tot = Rk.sum(axis=1, keepdims=True)
        bad = ~(np.isfinite(tot[:, 0]) & (tot[:, 0] > 0))
        out = Rk / np.where(bad[:, None], np.ones_like(tot), tot)     # unchanged arithmetic for other cells
        if bad.any():
            logpi = np.log(np.asarray(pi, np.float64)[keep])
            z = np.asarray(LL[bad][:, keep], np.float64) + logpi[None]
            ok = np.isfinite(z).any(axis=1)
            fallback = np.exp(logpi - logsumexp(logpi))
            zb = np.where(np.isfinite(z), z, -np.inf)
            out[bad] = np.where(ok[:, None], np.exp(zb - logsumexp(zb, axis=1, keepdims=True)), fallback[None])
        if not np.isfinite(out).all():
            raise FloatingPointError('nonfinite responsibilities after removing clones')
        return out

    def _update_offset(self, R, doses, keep):
        """Legacy clone-only pooled RNA residual for the kept latent population.

        Observed and expected pools must represent the same retained responsibilities. This heuristic
        excludes depth/programme terms; it is not a conditional-cell-mean residual or an exact NB
        likelihood maximization. Changing that mean alone fails the current recovery contract.
        """
        Rk = np.asarray(R[:, keep], np.float64)
        N = Rk.T @ self.L
        lr = np.log(self.base)[None] + np.asarray(doses)
        rates = np.exp(lr - logsumexp(lr, axis=1, keepdims=True))
        expected = (N[:, None] * rates).sum(axis=0)
        observed = self.rows.weighted_sum(Rk.sum(axis=1)[:, None])[0]
        o = np.log((observed + 0.5) / (expected + 0.5))
        self.offset = o - np.median(o)

    def _update_phi(self, R, doses):
        lr = self.log_rates(doses)
        rates = np.exp(lr)
        num = np.zeros(self.J)
        den = np.zeros(self.J)
        for a, b, Yc in self.rows.chunks():
            k = R[a:b].argmax(axis=1)
            rc = rates[k] * np.exp(self.x[a:b, None] * self.gamma[None] + self.U[a:b] @ self.W)
            m = self.L[a:b, None] * rc / rc.sum(axis=1, keepdims=True)
            num += ((Yc - m) ** 2 - m).sum(axis=0)
            den += (m ** 2).sum(axis=0)
        phi = np.clip(np.divide(num, den, out=np.full(self.J, 0.05), where=den > 0), 1e-3, 10.0)
        self.phi = 0.5 * phi + 0.5 * np.median(phi)

    # Divisive refinement -----------------------------------------------------
    def split_candidates(self, R, doses):
        """Within each tumour clone: Leiden subgroups and cells sharing an arm-level shift.

        Large clones are subsampled to max_proposal_cells (bounded memory); EM assigns the rest."""
        p = self.params
        mapc = R.argmax(axis=1)
        lr = self.log_rates(doses)
        out = []
        rho = self.Ytot / self.Ytot.sum()
        rng = np.random.default_rng(p.seed + len(doses))
        for k in range(1, len(doses)):
            cells = self._proposal_cells(np.flatnonzero(mapc == k), rng)
            if len(cells) < 2 * p.min_cells:
                continue
            out.extend(self._clone_candidates(k, cells, lr[k], rho))
        if self.allele_post is not None:
            for k in range(len(doses)):                    # clone 0 too: an allele-only (cnLOH) tumour
                cells = self._proposal_cells(np.flatnonzero(mapc == k), rng)
                if len(cells) >= 2 * p.min_cells:
                    out.extend((k, cells, g, src) for g, src in self.allele_splits(k, cells))
        return out

    def _clone_candidates(self, k, cells, lr_k, rho):
        """Leiden and event-mixture split proposals for clone k's cells (proposal bins or proposal half)."""
        p = self.params
        out = []
        Yc = self.rows.take(cells, self.prop_fold)
        Lc = self.L[cells]
        F = proposal_features(Yc, Lc if self.prop_fold is None else Lc * self.fold_scale['prop'], rho,
                              self.reset, self.window_bins, self.propose_mask)
        g = leiden_groups(smooth_features(F, self.graph, cells), p)
        if g.max() > 0:
            out.append((k, cells, g, 'leiden'))
        out.extend((k, cells, g, src) for g, src in self.event_mixtures(Yc, Lc, self.x[cells], lr_k,
                                                                     self.U[cells] @ self.W, cells))
        return out

    def _proposal_cells(self, cells, rng):
        n = self.params.max_proposal_cells
        return cells if len(cells) <= n else np.sort(rng.choice(cells, n, replace=False))

    def allele_splits(self, k, cells):
        """Per resolved imbalance run of clone k: cells whose reads there favour the imbalance over balance."""
        out = []
        A_all, D_all = self.alleles.A[cells], self.alleles.D[cells]
        for r in self.alleles.runs(self.allele_post[k], self.Dk[k]):
            if not r['resolved']:
                continue
            idx = np.arange(r['first_snp'], r['last_snp'] + 1)
            idx = idx[idx % 2 == 0]                          # propose on even SNPs; differ() tests odd SNPs
            b = self.allele_post[k][idx] @ state_b()
            b = 0.5 + self.alleles.params.cell_shrink * (b - 0.5)
            A, D = A_all[:, idx], D_all[:, idx]
            llr = balance_loglik_ratio(A, D, b, self.alleles.params.cell_kappa)
            g = (llr > 0).astype(int)
            if min(g.sum(), len(g) - g.sum()) >= self.params.min_cells:
                out.append((g, f"alleles_{r['chromosome']}:{r['start']}-{r['end']}"))
        return out

    def event_mixtures(self, Yc, Lc, xc, lr, uw, cells=None):
        """Cells sharing a shift of a chromosome, a chromosome half, or a sliding window of event_windows
        proposal bins (step half a window).

        Per region, a two-component mixture over cells, with each proposal bin's count
        O_cj ~ NB(E_cj * theta_k, size 1/phi_j) and free rate ratios theta_k shared by the region's bins,
        is fitted by EM initialization and numerical likelihood polishing on proposal bins only.
        Per-bin dispersions keep one overdispersed gene
        (FTH1, MALAT1, MT2A) from dominating a region, which a model of region totals with one
        median dispersion did not. The log-likelihood gain over one component ranks regions; the best
        few with distinct cell splits become split candidates, which odd bins then test.
        In split_mode 'thin', Yc is the proposal half on every bin, the NB size is eps / phi_j, and E is built
        from each cell's own proposal-half total, so it already carries the eps scaling.
        """
        p = self.params
        cols = np.flatnonzero(self.propose_mask)
        Yp = np.asarray(Yc[:, cols], np.float64)
        Lp = Yp.sum(axis=1)
        A = lr[cols][None] + xc[:, None] * self.gamma[cols][None] + uw[:, cols]
        Ep = Lp[:, None] * np.exp(A - logsumexp(A, axis=1, keepdims=True))   # per-cell, depth and programmes
        r_bin = 1.0 / np.maximum(self.phi[cols], 1e-3)
        thin = p.split_mode == 'thin'
        if thin:
            r_bin = THIN_EPS * r_bin                          # a thinned NB count has size eps / phi
        windows = tuple(2 * w for w in p.event_windows) if thin else p.event_windows   # same genomic span
        # Per-cell NB size factor over all proposal bins (the row-wise nb_rate fixed point): a cell's total
        # carries the noise of its overdispersed bins, which would otherwise shift every region of that cell.
        s = np.ones(len(Yp))
        for _ in range(20):
            d = s[:, None] * Ep + r_bin[None]
            s = (r_bin * Yp / d).sum(axis=1) / np.maximum((r_bin * Ep / d).sum(axis=1), 1e-300)
        Ep = Ep * np.maximum(s, 1e-12)[:, None]
        Rnb = np.broadcast_to(r_bin[None], Yp.shape)
        if self.graph is not None and cells is not None:
            # Spatial: counts summed over each cell's neighbourhood (proposal only). A sum of independent NB
            # counts with unequal means m_i and common dispersion phi has variance sum m + phi sum m^2, so its
            # moment-matched size is (sum m)^2 / (phi sum m^2); n times one cell's size overstated it whenever
            # depths differ (review 2026-11-02 #7: means 10 and 100, phi 0.2, variance 2,130 vs 1,320).
            J0 = len(cols)
            OE, _ = neighbourhood_sum(np.hstack([Yp, Ep, Ep ** 2]), self.graph, cells)
            Yp, Ep, E2 = OE[:, :J0], OE[:, J0:2 * J0], OE[:, 2 * J0:]
            Rnb = r_bin[None] * Ep ** 2 / np.maximum(E2, 1e-300)
        chrom = np.cumsum(self.reset)[cols]
        step = self.beta * np.log(1.5)
        found = []
        for c in np.unique(chrom):
            idx = np.flatnonzero(chrom == c)
            parts = [idx[: len(idx) // 2], idx[len(idx) // 2:], idx]
            for w in windows:
                if len(idx) > w:
                    parts += [idx[s: s + w] for s in range(0, len(idx) - w // 2, max(1, w // 2))]
            for part in parts:
                if len(part) < 2:
                    continue
                O, E = Yp[:, part], Ep[:, part]
                r_nb = Rnb[:, part]

                try:
                    fit = _fit_event_mixture(O, E, r_nb, step)
                except RuntimeError as error:
                    if getattr(self, 'log', None) is not None:
                        self.log.append(dict(action='event_mixture_unresolved', chromosome=int(c),
                                             first_bin=int(cols[part[0]]), last_bin=int(cols[part[-1]]),
                                             reason=str(error)))
                    continue
                gain = fit['gain']
                if 0 < gain <= fit['numerical_gain_budget']:
                    if getattr(self, 'log', None) is not None:
                        self.log.append(dict(action='event_mixture_numerically_ambiguous', chromosome=int(c),
                                             first_bin=int(cols[part[0]]), last_bin=int(cols[part[-1]]),
                                             gain=float(gain), numerical_gain_budget=fit['numerical_gain_budget']))
                    continue
                sel = fit['responsibility'] > 0.5
                if gain > 0 and p.min_cells <= sel.sum() <= len(sel) - p.min_cells:
                    found.append((gain, sel.astype(int),
                                  f'event_chrom{c}_bins{cols[part[0]]}-{cols[part[-1]]}'))
        found.sort(key=lambda x: -x[0])
        out, kept = [], []
        for _, g, src in found:                             # overlapping windows of one event give one candidate
            # Compare minority cell sets by Jaccard (labels may be flipped). Raw agreement >= 0.9 called two
            # disjoint 5% groups duplicates; 0.7 Jaccard is about 0.9 agreement for a 30% group.
            m = g.astype(bool)
            m = ~m if m.mean() > 0.5 else m
            if all((m & k).sum() / max((m | k).sum(), 1) < 0.7 for k in kept):
                out.append((g, src))
                kept.append(m)
            if len(out) == p.event_candidates:
                break
        return out

    def try_splits(self, R, doses, pi):
        """Accept at most one split per clone: the candidate with the strongest resolved difference."""
        best = {}
        for k, cells, g, source in self.split_candidates(R, doses):
            scored = self._score_candidate(k, cells, g, source)
            if scored is None:
                continue
            score, members, Ym, Lm, XLm, ADm = scored
            if score > best.get(k, (0,))[0]:
                labels = np.zeros(len(cells), dtype=int)
                for j, m in enumerate(members):
                    labels[np.isin(g, m)] = j
                best[k] = (score, cells, labels, Ym, Lm, XLm, ADm, source)
        if self.params.split_null_draws > 0:
            for k in list(best):
                passed, draws, null_max = self.split_null(k, best[k][1], doses, best[k][0])
                self.log.append(dict(action='split_null', clone=k, source=best[k][7], score=best[k][0],
                                     draws=draws, null_max=null_max, passed=passed,
                                     null_scope='selected_parent_expression_search',
                                     nuisance_refitted=False, initial_groups_pooled=self.params.split_null_initial,
                                     specimen_error_control=False, confirmation_status='not_independent',
                                     monte_carlo_resolution=1.0 / (self.params.split_null_draws + 1)))
                if not passed:
                    del best[k]
        if not best:
            return doses, pi, False
        return self._apply_splits(best, R, doses, pi)

    def _score_candidate(self, k, cells, g, source):
        """Agglomerate a proposed partition on test data and score its weakest depth-matched part pair.
        Returns None (fewer than two parts, or no difference) or (score, members, Ym, Lm, XLm, ADm), the pooled
        parts on full counts."""
        Yc = self.rows.take(cells, self.test_fold)
        onehot = np.eye(g.max() + 1, dtype=np.float32)[g]
        Yg = onehot.T.astype(np.float64) @ Yc
        Lg = onehot.T.astype(np.float64) @ self.L[cells]
        XLg = onehot.T.astype(np.float64) @ depth_moments(self.L, self.x)[cells]
        AD = None if self.Ak is None else self.alleles.pool(onehot, rows=cells)
        if self.test_fold is None:
            members, Ym, Lm, XLm, ADm = self.agglomerate(Yg, Lg, XLg, context=f'split_test_clone{k}_{source}',
                                                         AD=AD)
        else:
            # Agglomerate on the test half; the parts' profiles (new clone doses) use the full counts.
            s = self.fold_scale[self.test_fold]
            members, _, _, _, ADm = self.agglomerate(Yg, s * Lg, s * XLg, context=f'split_test_clone{k}_{source}',
                                                     AD=AD)
            Yg = onehot.T.astype(np.float64) @ self.rows.take(cells)
            Ym, Lm, XLm = (np.array([X[m].sum(axis=0) for m in members]) for X in (Yg, Lg, XLg))
        if len(members) < 2:
            return None
        # Every pair of parts must differ after depth matching (all parts are added as clones below);
        # the split's score is its weakest pair.
        part = np.full(self.rows.n, -1)
        for j, m in enumerate(members):
            part[cells[np.isin(g, m)]] = j
        pairs = [(a, b) for a in range(len(members)) for b in range(a + 1, len(members))]
        post = None if ADm is None else self.alleles.fit(ADm[0], ADm[1])[0]
        scores = []
        for (a, b), res in zip(pairs, self.matched_differences(part, pairs)):
            if not res.get('testable', True):
                self.log.append(dict(action='untestable_pair', clone=k, source=source, a=a, b=b,
                                     context='no depth overlap'))
            ad = (None,) * 4 if ADm is None else (ADm[0][a], ADm[1][a], ADm[0][b], ADm[1][b])
            s = sum(abs(x['z']) for x in self.supported(res['segments'], *ad, labels=part, a=a, b=b))
            if ADm is not None:
                ra, rb = self.alleles.runs(post[a], ADm[1][a]), self.alleles.runs(post[b], ADm[1][b])
                s += self.alleles.differ(ra, ADm[0][a], ADm[1][a], rb, ADm[0][b], ADm[1][b])
            scores.append(s)
        score = min(scores)
        if score <= 0:
            self.log.append(dict(action='reject_split', clone=k, source=source,
                                 context=f'no depth-matched difference in {sum(s <= 0 for s in scores)}/{len(pairs)} part pairs'))
            return None
        return score, members, Ym, Lm, XLm, ADm

    def split_null(self, k, cells, doses, observed):
        """Model-based filter of a selected parent's expression proposal-and-test search.

        Each null draw simulates counts for the same proposal cells from the fitted one-clone model, Y*_cj ~
        NB(mean L_c softmax(log rate_k + gamma x_c + U_c W)_j, size 1 / phi_j), reruns the clone's proposals
        (`_clone_candidates`) and scoring (`_score_candidate`), and keeps the best score (0 if none passes).
        The split passes if its score beats all `split_null_draws` null scores; the first
        null score >= observed stops the filter (same decision). The rank resolution is 1 / (draws + 1),
        but fitted parameters, parent selection and repeated searches preclude an exact p-value claim.
        Draws are seeded by (seed,
        clone, number of clones, draw). The null knows no cell-state structure beyond the fitted programme
        scores and has no gene veto (an approximation). Allele runs are rejected at construction.
        split_null_model 'cell_region' adds a log-normal rate effect per cell and chromosome half, e_cr ~ N(-s2_r / 2,
        s2_r), with s2_r estimated from the clone's own cells (`region_effect_variance`): a continuous,
        unimodal cell-by-region variation of the observed size, which a split must beat. A real subclone in the
        clone inflates s2_r, so this null costs power for true splits.
        split_null_model 'cov' draws e_c ~ N(-diag(S) / 2, S) with the full region covariance S
        (`region_effect_covariance`), so correlated shifts across several regions are in the null too (SigClust
        logic: a split must beat a Gaussian with the clone's own covariance). Same power cost, larger.
        Returns (passed, draws run, largest null score)."""
        p = self.params
        m = self.null_means(k, cells, doses)
        size = 1.0 / np.maximum(self.phi, 1e-3)
        lr = self.log_rates(doses)[k]
        if p.split_null_model in ('cell_region', 'cov'):
            region = chromosome_halves(self.reset)
            Yk = self.rows.take(cells)
            if p.split_null_model == 'cell_region':
                s2 = region_effect_variance(Yk, m, self.phi, region)
                S_half = np.diag(np.sqrt(s2))
            else:
                S = region_effect_covariance(Yk, m, self.phi, region)
                s2 = np.diag(S).copy()
                w, Q = np.linalg.eigh(S)
                S_half = Q * np.sqrt(np.maximum(w, 0.0))
            del Yk
        elif p.split_null_model != 'nb':
            raise ValueError(f"split_null_model must be 'nb', 'cell_region' or 'cov', got {p.split_null_model!r}")
        rho = self.Ytot / self.Ytot.sum()
        logL_mean = np.log(np.maximum(self.L, 1.0)).mean()
        local = np.arange(len(cells))
        null_max = 0.0
        for b in range(p.split_null_draws):
            rng = np.random.default_rng([int(p.seed) % 2 ** 32, 31, int(k), len(doses), b])
            mb = m
            if p.split_null_model in ('cell_region', 'cov'):
                e = rng.normal(0.0, 1.0, (len(cells), len(s2))) @ S_half.T - 0.5 * s2[None]
                mb = m * np.exp(e[:, region])
            Ysim = rng.negative_binomial(size[None], size[None] / (size[None] + mb)).astype(np.float32)
            nm = self._null_model(cells, Ysim, logL_mean, b)
            best = 0.0
            for _, c, g, source in nm._clone_candidates(k, local, lr, rho):
                scored = nm._score_candidate(k, c, g, source)
                if scored is not None:
                    best = max(best, scored[0])
            null_max = max(null_max, best)
            if best >= observed:
                return False, b + 1, null_max
        return True, p.split_null_draws, null_max

    def null_means(self, k, cells, doses):
        """[cells, J] expected counts of `cells` under clone k alone (depth slope and programme scores included)."""
        lr = self.log_rates(doses)[k]
        A = lr[None] + self.x[cells, None] * self.gamma[None] + self.U[cells] @ self.W
        return self.L[cells, None] * np.exp(A - logsumexp(A, axis=1, keepdims=True))

    def _null_model(self, cells, Ysim, logL_mean, b):
        """Shallow copy of this model over simulated counts for `cells` (local rows 0..n-1); fitted profiles,
        slopes, dispersions and programme loadings are shared, data-dependent state is replaced."""
        import copy
        nm = copy.copy(self)
        nm.rows = Rows(Ysim, np.arange(len(cells)), self.rows.chunk)
        if self.thin:
            nm.rows.set_thinning((int(self.params.seed) * 7919 + 37 + b) % 2 ** 63, self.thin_size, THIN_EPS)
        nm.L = Ysim.sum(axis=1, dtype=np.float64)
        nm.x = np.log(np.maximum(nm.L, 1.0)) - logL_mean
        nm.U = self.U[cells]
        nm.graph = None if self.graph is None else sparse.csr_matrix(self.graph[cells][:, cells])
        nm.genes = None
        nm.alleles = None
        nm.Ak = nm.Dk = nm.B = nm.allele_post = None
        nm.log, nm.R_prev = [], None
        return nm

    def _apply_splits(self, best, R, doses, pi):
        doses, pi = list(doses), list(pi)
        pools = None if self.Ak is None else ([a for a in self.Ak], [d for d in self.Dk])
        R_new = [R[:, j] for j in range(R.shape[1])]
        for k, (score, cells, labels, Ym, Lm, XLm, ADm, source) in best.items():
            prof = self.profiles(Ym, Lm, XLm)
            share = pi[k]
            parent = R_new[k].copy()
            for j in range(1, len(Ym)):
                doses.append(prof['dose'][j])
                pi.append(share * (labels == j).mean())
                col = np.zeros_like(parent)
                col[cells] = parent[cells] * (labels == j)
                R_new.append(col)
                if pools is not None:
                    pools[0].append(ADm[0][j])
                    pools[1].append(ADm[1][j])
            R_new[k] = parent.copy()
            R_new[k][cells] = parent[cells] * (labels == 0)
            if k > 0:                                       # clone 0 stays copy-neutral
                doses[k] = prof['dose'][0]
            pi[k] = share * (labels == 0).mean()
            if pools is not None:
                pools[0][k], pools[1][k] = ADm[0][0], ADm[1][0]
            self.log.append(dict(action='split', clone=k, source=source, parts=len(Ym), score=score,
                                 n_eff='/'.join(str(int((labels == j).sum())) for j in range(len(Ym)))))
        if pools is not None:
            self.set_allele_pools(np.array(pools[0]), np.array(pools[1]))
        if self.graph is not None:
            Rn = np.column_stack(R_new)
            self.R_prev = Rn / np.maximum(Rn.sum(axis=1, keepdims=True), 1e-12)
        pi = np.asarray(pi)
        return doses, pi / pi.sum(), True


def _name(members, labels):
    return '+'.join(str(m if labels is None else labels[m]) for m in members)


def _finalize_search(model, doses, pi):
    """Validate a budget-exhausted search without accepting another split.

    A last-cycle split needs the same merge check as an earlier split. Refitting
    after a merge can change the other groups' support, so check again. Strictly
    decreasing group counts bound this finalization even when search oscillates.
    """
    for iteration in range(len(doses)):
        fitted = model.em(doses, pi)
        if len(fitted[0]) > len(doses):
            raise RuntimeError('Terminal EM must not add groups')
        doses, pi, merged = model.merge_unsupported(fitted[2], fitted[0], fitted[1])
        if not merged:
            model.log.append(dict(action='search_stop', reason='cycle_budget_exhausted',
                                  terminal_merge_checked=True, terminal_merge_rounds=iteration))
            return fitted
        if len(doses) >= len(fitted[0]):
            raise RuntimeError('Terminal merge must reduce the number of groups')
    raise RuntimeError('Terminal merge validation did not terminate')


def fit_specimen(rows: Rows, L, mu, reset, phi, params: CloneParams, programmes=None, alleles=None, graph=None,
                 gene_loader=None, gamma=None):
    """Fit one specimen. Returns a dict of arrays and tables (bin indices, not coordinates).

    gene_loader(local_rows) -> (gene counts CSR, gene bin): enables gene-consistency checks."""
    if gene_loader is not None and params.split_null_draws > 0:
        raise ValueError('Gene-veto-enabled split null is unsupported: simulated bin counts cannot reproduce '
                         'the observed gene-consistency veto. Use the exploratory workflow without this filter; '
                         'bin-only null diagnostics do not establish specimen error control.')
    model = SpecimenModel(rows, L, mu, reset, phi, params, programmes, alleles, graph, gamma)
    rng = np.random.default_rng(params.seed)
    n = rows.n
    # Proposals on a bounded subsample, cross-fitted bins (or count halves in split_mode 'thin').
    sub = np.sort(rng.choice(n, min(n, params.max_proposal_cells), replace=False))
    Ys, Ls, XLs = rows.take(sub), model.L[sub], depth_moments(model.L, model.x)[sub]
    if model.thin:
        Yp = rows.take(sub, model.prop_fold)
        Yt = Ys - Yp
    if gene_loader is not None:
        Xg, gbin = gene_loader(sub)
        if model.thin:
            # Gene-consistency checks use the test half of each gene, nested in the bin-level split.
            Xg = thin_gene_counts(Xg, gbin, Ys, Yt, params.seed)
        model.genes = (sub, Xg, gbin)
    if not model.thin:
        F = proposal_features(Ys, Ls, model.Ytot / model.Ytot.sum(), reset, params.window_bins,
                              model.propose_mask)
    else:
        F = proposal_features(Yp, Ls * model.fold_scale['prop'], model.Ytot / model.Ytot.sum(), reset,
                              model.window_bins, model.propose_mask)
    g = leiden_groups(smooth_features(F, graph, sub), params)
    onehot = np.eye(g.max() + 1, dtype=np.float32)[g]
    AD = None if alleles is None else alleles.pool(onehot, rows=sub)
    if not model.thin:
        members, Yk, Lk, XLk, AD = model.agglomerate(onehot.T.astype(np.float64) @ Ys,
                                                     onehot.T.astype(np.float64) @ Ls,
                                                     onehot.T.astype(np.float64) @ XLs, AD=AD)
    else:
        s = model.fold_scale['test']
        members, _, _, _, AD = model.agglomerate(onehot.T.astype(np.float64) @ Yt, s * (onehot.T.astype(np.float64) @ Ls),
                                                 s * (onehot.T.astype(np.float64) @ XLs), AD=AD)
        Yk, Lk, XLk = (np.array([X[m].sum(axis=0) for m in members])
                       for X in (onehot.T.astype(np.float64) @ Ys, onehot.T.astype(np.float64) @ Ls,
                                 onehot.T.astype(np.float64) @ XLs))
    prof = model.profiles(Yk, Lk, XLk)
    gpost = None if AD is None else alleles.fit(*AD)[0]
    doses, shares, normal_share = [np.zeros(model.J)], [], 0.0
    pools = None if AD is None else ([np.zeros(AD[0].shape[1])], [np.zeros(AD[0].shape[1])])
    tum = []
    for k in range(len(Yk)):
        frac = float(np.isin(g, members[k]).mean())
        allele_event = gpost is not None and alleles.has_events(gpost[k], AD[1][k])
        if allele_event or model.has_events(Yk[k], prof['E'][k], prof['post'][k], prof['path'][k], prof['dose'][k],
                                            prof['tau'][k]):
            doses.append(prof['dose'][k])
            shares.append(frac)
            tum.append(k)
            if pools is not None:
                pools[0].append(AD[0][k])
                pools[1].append(AD[1][k])
        else:
            normal_share += frac
            if pools is not None:
                pools[0][0] = pools[0][0] + AD[0][k]
                pools[1][0] = pools[1][0] + AD[1][k]
    if params.split_null_draws > 0 and params.split_null_initial and len(doses) > 2:
        # Initial tumour groups are not null-gated: pool them into one tumour clone and let null-gated splits
        # rebuild any structure.
        Yt, Lt, XLt = Yk[tum].sum(axis=0), Lk[tum].sum(), XLk[tum].sum(axis=0)
        doses = [doses[0], model.profiles(Yt, Lt, XLt)['dose'][0]]
        shares = [float(sum(shares))]
        if pools is not None:
            pools = ([pools[0][0], np.sum(pools[0][1:], axis=0)], [pools[1][0], np.sum(pools[1][1:], axis=0)])
        model.log.append(dict(action='pool_initial', groups=len(tum)))
    if pools is not None:
        model.set_allele_pools(np.array(pools[0]), np.array(pools[1]))
    pi = np.array([max(normal_share, params.normal_prior)] + shares)
    pi /= pi.sum()
    budget_exhausted = False
    for cycle in range(params.cycles):
        doses, pi, R, LL, Yk, Lk, XLk = model.em(doses, pi)
        # Merge tumour clones without depth-matched dosage or allele differences.
        doses, pi, merged = model.merge_unsupported(R, doses, pi)
        if merged:
            continue
        doses, pi, split = model.try_splits(R, doses, pi)
        if not split:
            break
    else:
        budget_exhausted = True
        doses, pi, R, LL, Yk, Lk, XLk = _finalize_search(model, doses, pi)
    prof = model.profiles(Yk, Lk, XLk)
    evidence = dict(group_interpretation='RNA_inferred_candidate',
                    confirmation_status='not_independent', specimen_error_control=False,
                    split_filter='selected_parent_expression_search' if params.split_null_draws else 'disabled',
                    initial_groups_pooled=bool(params.split_null_draws and params.split_null_initial),
                    search_termination='cycle_budget_exhausted' if budget_exhausted else 'no_accepted_split',
                    terminal_merge_checked=True,
                    query_derived_reference_opt_in=params.in_sample_reference)
    return dict(model=model, doses=np.asarray(doses), pi=pi, R=R, LL=LL, Yk=Yk, Lk=Lk, XLk=XLk, prof=prof,
                evidence=evidence)


def neighbour_graph(xy, k, *, cell_ids=None):
    """Row-normalized k-nearest-neighbour graph over one specimen's coordinates, self always excluded.

    With explicit cell_ids, search uses canonical ID order and restores the input
    order afterward, so arbitrary distance ties do not depend on row permutation.
    Ties have no biological preference. Calls without IDs retain legacy row-order
    tie handling. Within returned candidates neighbours are ordered by (distance, index), and the cell itself is removed by index rather than by
    position: at duplicated coordinates the first neighbour returned need not be the cell (review 2026-11-02 #7:
    60 cells on 20 triplicated points gave 20 self-edges and lambda 1.0 for random labels). The graph is
    directed and row-normalized, so G @ R is a neighbourhood-average smoothing heuristic, not the mean-field
    gradient of a symmetric Potts energy."""
    from scipy import sparse
    from sklearn.neighbors import NearestNeighbors
    n = len(xy)
    order = None
    if cell_ids is not None:
        ids = np.asarray(cell_ids, dtype=object)
        if (ids.ndim != 1 or len(ids) != n
                or any(not isinstance(value, str) or not value.strip() for value in ids)
                or len(set(ids)) != n):
            raise ValueError('spatial cell IDs must be aligned, unique, nonempty strings')
        order = np.argsort(ids, kind='stable')
        xy = np.asarray(xy)[order]
    k = int(min(k, n - 1))
    if k < 1:
        return None
    dist, idx = NearestNeighbors(n_neighbors=k + 1).fit(xy).kneighbors(xy)
    dist = np.where(idx == np.arange(n)[:, None], np.inf, dist)
    nearest = np.lexsort((idx, dist), axis=-1)[:, :k]
    cols = np.take_along_axis(idx, nearest, axis=1).ravel()
    rows = np.repeat(np.arange(n), k)
    graph = sparse.csr_matrix((np.full(len(cols), 1.0 / k, np.float32), (rows, cols)), shape=(n, n))
    if order is not None:
        inverse = np.argsort(order)
        graph = graph[inverse][:, inverse]
    return graph


def clone_order(R):
    """Clone 0 stays first; tumour clones ordered by decreasing MAP size."""
    mapc = R.argmax(axis=1)
    sizes = np.bincount(mapc, minlength=R.shape[1])
    return [0] + sorted(range(1, R.shape[1]), key=lambda k: -sizes[k])


def tree_table(dist_genes, names):
    """Average-linkage tree over clones; distance = genes in resolved differing segments."""
    k = len(names)
    if k < 2:
        return []
    cond = dist_genes[np.triu_indices(k, 1)]
    Z = linkage(cond, method='average')
    nodes = {i: names[i] for i in range(k)}
    rows = []
    for i, (a, b, h, size) in enumerate(Z):
        node = f'node{i + 1}'
        rows.append(dict(node=node, left=nodes[int(a)], right=nodes[int(b)],
                         genes_differing=float(h), n_clones=int(size)))
        nodes[k + i] = node
    return rows
