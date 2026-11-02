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
enters the cell likelihood, pooled expectations and difference tests, so clones
are not split by depth. Expression programmes learned on normal cells
(`haplodose.programmes`) enter the cell likelihood as per-cell nuisance scores,
so a cell with, say, an immune-like gene cluster switched on is not called a
clone by that cluster. Pooled profiles and pooled differences get a free
amplitude per programme, so only changes no programme combination explains can
become segments or separate clones. Each pooled group carries XL = (sum of L_c x_c, sum of L_c x_c^2), with
x_c the centred log depth: the UMI-weighted mean and variance of depth, from which the pooled expectation is
built (`SpecimenModel.expected`).

Proposals use even bins and difference tests use odd bins (cross-fitting), so a
split cannot be supported only by the noise that created it. Clone 0 is copy
neutral relative to the reference baseline (normal-like).

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
                      weighted_median)
from .alleles import state_b
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
    in_sample_reference: bool = True      # with alleles: allele-balanced cells as an extra baseline group
    spatial_neighbors: int = 8
    potts_lambda: float = -1.0           # < 0: chosen from neighbour-predicted cell likelihoods
    n_programmes: int = 10
    programme_window: int = 0            # bins of running mean removed; 0 = arm means (hg38 centromeres)
    depth_source: str = 'query'          # 'reference': depth slopes from fit reference cells, held fixed
    seed: int = 20261031
    profile: ProfileParams = field(default_factory=ProfileParams)


class Rows:
    """Bounded-chunk access to one specimen's rows of the cells x bins matrix."""

    def __init__(self, Y, rows, chunk=4096):
        self.Y, self.rows, self.chunk = Y, np.asarray(rows), int(chunk)
        self.n = len(self.rows)

    def chunks(self):
        for a in range(0, self.n, self.chunk):
            b = min(a + self.chunk, self.n)
            yield a, b, np.asarray(self.Y[self.rows[a:b]], dtype=np.float32)

    def take(self, idx):
        return np.asarray(self.Y[self.rows[np.asarray(idx)]], dtype=np.float32)

    def weighted_sum(self, R):
        """R [n, G] dense or scipy sparse -> [G, bins] float64."""
        out = np.zeros((R.shape[1], self.Y.shape[1]), np.float64)
        for a, b, Yc in self.chunks():
            out += np.asarray(R[a:b].T.astype(np.float64) @ Yc)
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


def mixture_weights(Ytot, mu, iters=500):
    """Poisson MLE of reference-group weights (each mu column sums to 1)."""
    w = np.full(mu.shape[1], 1.0 / mu.shape[1])
    for _ in range(iters):
        mix = mu @ w
        w = (Ytot[:, None] * mu * w / mix[:, None]).sum(axis=0)
        w /= w.sum()
    return w


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
            random_state=params.seed).fit_transform(F)
    k = int(min(params.n_neighbors, n - 1))
    _, idx = NearestNeighbors(n_neighbors=k + 1).fit(X).kneighbors(X)
    edges = np.column_stack([np.repeat(np.arange(n), k), idx[:, 1:].ravel()])
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
    onehot = np.zeros((rows.n, group.max() + 1), np.float32)
    onehot[np.arange(rows.n), group] = 1
    Yg = rows.weighted_sum(onehot)
    Lg = onehot.T.astype(np.float64) @ L
    xg = (onehot.T.astype(np.float64) @ (L * x)) / np.maximum(Lg, 1.0)
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


class SpecimenModel:
    def __init__(self, rows: Rows, L, mu, reset, phi, params: CloneParams, programmes=None, alleles=None,
                 graph=None, gamma=None):
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
        self.propose_mask = np.arange(self.J) % 2 == 0
        self.test_mask = ~self.propose_mask
        self.log = []
        self.alleles = alleles
        self.Ak = self.Dk = self.B = self.allele_post = None
        self.graph, self.lam = graph, 0.0
        self.genes = None               # (subsample rows, gene counts CSR [sub, genes], gene bin) for gene checks
        self.R_prev = None              # mean-field state carried between E-steps (spatial only)

    # Spatial Potts prior (mean-field) ------------------------------------------
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
            R = self.R_prev               # variational EM: continue from the previous labels
        for _ in range(iters):
            R = posteriors(LL + self.lam * np.asarray(self.graph @ R, np.float32), logpi)
        self.R_prev = R
        return R

    def choose_lambda(self, LL, R0, logpi, grid=(0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)):
        """lambda maximizing how well neighbours predict each cell's own data.

        score = sum_c log sum_k softmax(log pi + lambda (G R0)_c)_k exp(LL_ck), with R0
        the data-only posteriors: the cell's own labels never enter its prior, so
        spatially random clones give lambda = 0."""
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
        rows = segments_table(post, path, dose, Y, E, self.reset, self.beta, tau, self.pp)
        return any(r['resolved'] for r in rows)

    def log_rates(self, doses):
        lr = np.log(self.base)[None] + self.offset[None] + np.asarray(doses)
        return lr - logsumexp(lr, axis=1, keepdims=True)

    # Group merging --------------------------------------------------------
    def differences(self, Ya, La, XLa, Yb, Lb, XLb):
        shift = self.depth_log_rate(La, XLa) - self.depth_log_rate(Lb, XLb)
        if self.programmes is not None:
            Ya, Yb = np.atleast_2d(Ya), np.atleast_2d(Yb)
            d = (np.log((Ya + 0.5) / np.atleast_1d(La)[:, None]) - np.log((Yb + 0.5) / np.atleast_1d(Lb)[:, None])
                 - shift)
            v = 1.0 / (Ya + 0.5) + 1.0 / (Yb + 0.5)
            shift = shift + fit_amplitudes(d - np.median(d, axis=1, keepdims=True), 1.0 / (0.04 + v), self.W)
        return difference_segments(Ya, La, Yb, Lb, self.reset, self.beta, self.pp, mask=self.test_mask,
                                   shift=shift)

    def matched_differences(self, labels, pairs, n_strata=8):
        """Differences between groups of cells with depth-matched pooling.

        labels [n]: group per cell (-1: none). For each pair, cells are weighted so
        both groups have the same distribution over log-depth strata (coarsened
        exact matching); depth-driven composition then cancels in the difference,
        whatever its functional form. Groups without depth overlap have no power."""
        # Sparse cells x 2P weights: each cell is in at most the pairs that include its group (a dense
        # matrix was cells x 2P float32 plus a float64 copy, ~17 GiB at 1M cells and 780 pairs).
        r_idx, c_idx, val = [], [], []
        for i, (a, b) in enumerate(pairs):
            ia, ib = np.flatnonzero(labels == a), np.flatnonzero(labels == b)
            edges = np.quantile(np.r_[self.x[ia], self.x[ib]], np.linspace(0, 1, n_strata + 1)[1:-1])
            sa, sb = np.searchsorted(edges, self.x[ia]), np.searchsorted(edges, self.x[ib])
            na, nb = np.bincount(sa, minlength=n_strata), np.bincount(sb, minlength=n_strata)
            m = np.minimum(na, nb)
            r_idx += [ia, ib]
            c_idx += [np.full(len(ia), 2 * i), np.full(len(ib), 2 * i + 1)]
            val += [(m / np.maximum(na, 1))[sa], (m / np.maximum(nb, 1))[sb]]
        W = sparse.csr_matrix((np.concatenate(val), (np.concatenate(r_idx), np.concatenate(c_idx))),
                              shape=(self.rows.n, 2 * len(pairs)), dtype=np.float64)
        Y, L, XL = self.pooled(W)
        return self.differences(Y[0::2], L[0::2], XL[0::2], Y[1::2], L[1::2], XL[1::2])

    def supported(self, segments, Aa, Da, Ab, Db, labels=None, a=None, b=None, log=True):
        """Difference segments that pass the gene-consistency test and are not vetoed by alleles."""
        keep = []
        for s in segments:
            where = f"bins{s['first_bin']}-{s['last_bin']}"
            if labels is not None and not self.gene_consistent(s, labels, a, b):
                if log:
                    self.log.append(dict(action='gene_veto', context=where))
                continue
            if (self.alleles is not None and getattr(self.alleles, 'bins', None) is not None and Aa is not None
                    and self.alleles.veto(self.alleles.snps_in(self.alleles.bins, s['first_bin'], s['last_bin']),
                                          Aa, Da, Ab, Db)):
                if log:
                    self.log.append(dict(action='allele_veto', context=where))
                continue
            keep.append(s)
        return keep

    def gene_consistent(self, seg, labels, a, b, min_genes=8, alpha=0.01):
        """Do most genes of a difference segment move with it? (one-sided sign test, test bins only)

        Copy number shifts every gene of a segment; cell-state programmes or a few highly
        expressed genes (e.g. ribosomal clusters) shift a minority. Genes are pooled over
        the proposal subsample; log ratios are centred by the genome-wide median."""
        if self.genes is None:
            return True
        sub, Xg, gbin = self.genes
        la = labels[sub]
        ia, ib = la == a, la == b
        if ia.sum() < 10 or ib.sum() < 10:
            return True
        ya, yb = np.asarray(Xg[ia].sum(axis=0)).ravel(), np.asarray(Xg[ib].sum(axis=0)).ravel()
        r = np.log((ya + 0.5) / ya.sum()) - np.log((yb + 0.5) / yb.sum())
        informative = (ya + yb) >= 10
        if informative.sum() < min_genes:
            return True
        r = r - np.median(r[informative])
        sel = (informative & (gbin >= seg['first_bin']) & (gbin <= seg['last_bin'])
               & self.test_mask[np.clip(gbin, 0, self.J - 1)])
        n = int(sel.sum())
        if n < min_genes:
            return True
        k = int((np.sign(r[sel]) == (1 if seg['direction'] == 'higher' else -1)).sum())
        from scipy.stats import binom
        return bool(binom.sf(k - 1, n, 0.5) < alpha)

    def merge_unsupported(self, R, doses, pi):
        """Merge tumour clones with no depth-matched dosage difference and no allele difference."""
        merged = False
        while len(doses) > 2:
            labels = R.argmax(axis=1)
            tum = [k for k in range(1, len(doses)) if (labels == k).sum() >= self.params.min_cells]
            pairs = [(a, b) for i, a in enumerate(tum) for b in tum[i + 1:]]
            if not pairs:
                break
            res = self.matched_differences(labels, pairs)
            cand = []
            for (a, b), r in zip(pairs, res):
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
            self.log.append(dict(action='merge_matched', a=a, b=b, n_a=int((labels == a).sum()),
                                 n_b=int((labels == b).sum())))
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
        if merged:
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
    def loglik(self, doses):
        lr = self.log_rates(doses)
        out = np.empty((self.rows.n, len(doses)), np.float32)
        allele = self.B is not None and len(self.B) == len(doses)
        for a, b, Yc in self.rows.chunks():
            out[a:b] = cell_loglik(Yc, self.L[a:b], lr, self.phi, self.x[a:b], self.gamma, self.U[a:b], self.W)
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

    def pooled(self, R):
        Rt = R.T.astype(np.float64)
        return self.rows.weighted_sum(R), np.asarray(Rt @ self.L).ravel(), np.asarray(Rt @ depth_moments(self.L, self.x))

    def em(self, doses, pi):
        """Clone 0 keeps dose 0. Returns doses, pi, responsibilities, loglik, profiles."""
        p = self.params
        prev = None
        for _ in range(p.outer_rounds):
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
            self.update_scores(R[:, keep] / R[:, keep].sum(axis=1, keepdims=True), new)
            changed = len(keep) != len(doses)
            doses, pi = new, pi[keep] / pi[keep].sum()
            mapc = R[:, keep].argmax(axis=1)
            if not changed and prev is not None and np.mean(mapc != prev) < 0.005:
                break
            prev = None if changed else mapc
        LL = self.loglik(doses)
        R = self.posterior(LL, pi)
        Yk, Lk, XLk = self.pooled(R)
        if self.alleles is not None:
            self.set_allele_pools(*self.alleles.pool(R))
        self.log.append(dict(action='em', n_eff=' '.join(f'{x:.0f}' for x in R.sum(axis=0)), potts_lambda=self.lam))
        return doses, pi, R, LL, Yk, Lk, XLk

    def _update_offset(self, R, doses, keep):
        Rk = R[:, keep]
        N = Rk.T.astype(np.float64) @ self.L
        lr = np.log(self.base)[None] + np.asarray(doses)
        rates = np.exp(lr - logsumexp(lr, axis=1, keepdims=True))
        expected = (N[:, None] * rates).sum(axis=0)
        o = np.log((self.Ytot + 0.5) / (expected + 0.5))
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
            Yc = self.rows.take(cells)
            Lc = self.L[cells]
            F = proposal_features(Yc, Lc, rho, self.reset, p.window_bins, self.propose_mask)
            g = leiden_groups(smooth_features(F, self.graph, cells), p)
            if g.max() > 0:
                out.append((k, cells, g, 'leiden'))
            out.extend((k, cells, g, src) for g, src in self.event_mixtures(Yc, Lc, self.x[cells], lr[k],
                                                                         self.U[cells] @ self.W, cells))
        if self.allele_post is not None:
            for k in range(len(doses)):                    # clone 0 too: an allele-only (cnLOH) tumour
                cells = self._proposal_cells(np.flatnonzero(mapc == k), rng)
                if len(cells) >= 2 * p.min_cells:
                    out.extend((k, cells, g, src) for g, src in self.allele_splits(k, cells))
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
            llr = np.asarray(A @ np.log(b / 0.5) + (D - A) @ np.log((1 - b) / 0.5)).ravel()
            g = (llr > 0).astype(int)
            if min(g.sum(), len(g) - g.sum()) >= self.params.min_cells:
                out.append((g, f"alleles_{r['chromosome']}:{r['start']}-{r['end']}"))
        return out

    def event_mixtures(self, Yc, Lc, xc, lr, uw, cells=None):
        """Cells sharing a shift of a chromosome, a chromosome half, or a sliding window of event_windows
        proposal bins (step half a window).

        Per region, a two-component mixture over cells, with each proposal bin's count
        O_cj ~ NB(E_cj * theta_k, size 1/phi_j) and free rate ratios theta_k shared by the region's bins,
        is fitted by EM on proposal bins only. Per-bin dispersions keep one overdispersed gene
        (FTH1, MALAT1, MT2A) from dominating a region, which a model of region totals with one
        median dispersion did not. The log-likelihood gain over one component ranks regions; the best
        few with distinct cell splits become split candidates, which odd bins then test.
        """
        from scipy.special import xlogy
        p = self.params
        cols = np.flatnonzero(self.propose_mask)
        Yp = np.asarray(Yc[:, cols], np.float64)
        Lp = Yp.sum(axis=1)
        A = lr[cols][None] + xc[:, None] * self.gamma[cols][None] + uw[:, cols]
        Ep = Lp[:, None] * np.exp(A - logsumexp(A, axis=1, keepdims=True))   # per-cell, depth and programmes
        r_bin = 1.0 / np.maximum(self.phi[cols], 1e-3)
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
            for w in p.event_windows:
                if len(idx) > w:
                    parts += [idx[s: s + w] for s in range(0, len(idx) - w // 2, max(1, w // 2))]
            for part in parts:
                if len(part) < 2:
                    continue
                O, E = Yp[:, part], Ep[:, part]
                r_nb = Rnb[:, part]

                def logp(theta):
                    # Per-cell NB log pmf summed over the region's bins, without the gammaln terms: constant
                    # per cell, they cancel in responsibilities and in the gain. xlogy: m = 0 and O = 0 give 0.
                    m = E * max(theta, 1e-12)
                    return (xlogy(O, m / (m + r_nb)) + r_nb * np.log(r_nb / (m + r_nb))).sum(axis=1)
                theta1 = nb_rate(O, E, r_nb, 1.0)
                ll1 = logp(theta1).sum()
                theta = theta1 * np.exp(np.array([-0.5, 0.5]) * step)
                pi, prev = 0.5, -np.inf
                for _ in range(200):
                    la, lb = np.log1p(-pi) + logp(theta[0]), np.log(pi) + logp(theta[1])
                    total = np.logaddexp(la, lb)
                    ll = float(total.sum())
                    r = np.exp(lb - total)
                    if ll - prev < 1e-6 * max(abs(ll), 1.0):
                        break
                    prev = ll
                    pi = float(np.clip(r.mean(), 1e-3, 1 - 1e-3))
                    theta = np.array([nb_rate(O, E, r_nb, (1 - r)[:, None], theta[0]),
                                      nb_rate(O, E, r_nb, r[:, None], theta[1])])
                gain = float(total.sum() - ll1)
                sel = r > 0.5
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
            Yc = self.rows.take(cells)
            onehot = np.eye(g.max() + 1, dtype=np.float32)[g]
            Yg = onehot.T.astype(np.float64) @ Yc
            Lg = onehot.T.astype(np.float64) @ self.L[cells]
            XLg = onehot.T.astype(np.float64) @ depth_moments(self.L, self.x)[cells]
            AD = None if self.Ak is None else self.alleles.pool(onehot, rows=cells)
            members, Ym, Lm, XLm, ADm = self.agglomerate(Yg, Lg, XLg, context=f'split_test_clone{k}_{source}', AD=AD)
            if len(members) < 2:
                continue
            # Every pair of parts must differ after depth matching (all parts are added as clones below);
            # the split's score is its weakest pair.
            part = np.full(self.rows.n, -1)
            for j, m in enumerate(members):
                part[cells[np.isin(g, m)]] = j
            pairs = [(a, b) for a in range(len(members)) for b in range(a + 1, len(members))]
            post = None if ADm is None else self.alleles.fit(ADm[0], ADm[1])[0]
            scores = []
            for (a, b), res in zip(pairs, self.matched_differences(part, pairs)):
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
                continue
            if score > best.get(k, (0,))[0]:
                labels = np.zeros(len(cells), dtype=int)
                for j, m in enumerate(members):
                    labels[np.isin(g, m)] = j
                best[k] = (score, cells, labels, Ym, Lm, XLm, ADm, source)
        if not best:
            return doses, pi, False
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


def fit_specimen(rows: Rows, L, mu, reset, phi, params: CloneParams, programmes=None, alleles=None, graph=None,
                 gene_loader=None, gamma=None):
    """Fit one specimen. Returns a dict of arrays and tables (bin indices, not coordinates).

    gene_loader(local_rows) -> (gene counts CSR, gene bin): enables gene-consistency checks."""
    model = SpecimenModel(rows, L, mu, reset, phi, params, programmes, alleles, graph, gamma)
    rng = np.random.default_rng(params.seed)
    n = rows.n
    # Proposals on a bounded subsample, cross-fitted bins.
    sub = np.sort(rng.choice(n, min(n, params.max_proposal_cells), replace=False))
    if gene_loader is not None:
        model.genes = (sub, *gene_loader(sub))
    Ys, Ls, XLs = rows.take(sub), model.L[sub], depth_moments(model.L, model.x)[sub]
    F = proposal_features(Ys, Ls, model.Ytot / model.Ytot.sum(), reset, params.window_bins,
                          model.propose_mask)
    g = leiden_groups(smooth_features(F, graph, sub), params)
    onehot = np.eye(g.max() + 1, dtype=np.float32)[g]
    AD = None if alleles is None else alleles.pool(onehot, rows=sub)
    members, Yk, Lk, XLk, AD = model.agglomerate(onehot.T.astype(np.float64) @ Ys, onehot.T.astype(np.float64) @ Ls,
                                                 onehot.T.astype(np.float64) @ XLs, AD=AD)
    prof = model.profiles(Yk, Lk, XLk)
    gpost = None if AD is None else alleles.fit(*AD)[0]
    doses, shares, normal_share = [np.zeros(model.J)], [], 0.0
    pools = None if AD is None else ([np.zeros(AD[0].shape[1])], [np.zeros(AD[0].shape[1])])
    for k in range(len(Yk)):
        frac = float(np.isin(g, members[k]).mean())
        allele_event = gpost is not None and alleles.has_events(gpost[k], AD[1][k])
        if allele_event or model.has_events(Yk[k], prof['E'][k], prof['post'][k], prof['path'][k], prof['dose'][k],
                                            prof['tau'][k]):
            doses.append(prof['dose'][k])
            shares.append(frac)
            if pools is not None:
                pools[0].append(AD[0][k])
                pools[1].append(AD[1][k])
        else:
            normal_share += frac
            if pools is not None:
                pools[0][0] = pools[0][0] + AD[0][k]
                pools[1][0] = pools[1][0] + AD[1][k]
    if pools is not None:
        model.set_allele_pools(np.array(pools[0]), np.array(pools[1]))
    pi = np.array([max(normal_share, params.normal_prior)] + shares)
    pi /= pi.sum()
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
        doses, pi, R, LL, Yk, Lk, XLk = model.em(doses, pi)
    prof = model.profiles(Yk, Lk, XLk)
    return dict(model=model, doses=np.asarray(doses), pi=pi, R=R, LL=LL, Yk=Yk, Lk=Lk, XLk=XLk, prof=prof)


def neighbour_graph(xy, k):
    """Row-normalized k-nearest-neighbour graph over one specimen's coordinates, self always excluded.

    Neighbours are ordered by (distance, index), and the cell itself is removed by index rather than by
    position: at duplicated coordinates the first neighbour returned need not be the cell (review 2026-11-02 #7:
    60 cells on 20 triplicated points gave 20 self-edges and lambda 1.0 for random labels). The graph is
    directed and row-normalized, so G @ R is a neighbourhood-average smoothing heuristic, not the mean-field
    gradient of a symmetric Potts energy."""
    from scipy import sparse
    from sklearn.neighbors import NearestNeighbors
    n = len(xy)
    k = int(min(k, n - 1))
    if k < 1:
        return None
    dist, idx = NearestNeighbors(n_neighbors=k + 1).fit(xy).kneighbors(xy)
    dist = np.where(idx == np.arange(n)[:, None], np.inf, dist)
    order = np.lexsort((idx, dist), axis=-1)[:, :k]
    cols = np.take_along_axis(idx, order, axis=1).ravel()
    rows = np.repeat(np.arange(n), k)
    return sparse.csr_matrix((np.full(len(cols), 1.0 / k, np.float32), (rows, cols)), shape=(n, n))


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
