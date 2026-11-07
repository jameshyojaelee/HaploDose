# HaploDose: separating relative RNA dosage, DNA copy-number change and clone identity in single-cell, spatial and bulk RNA-seq

Draft 2026-11-03. Working draft for internal review, not a submission. Every number carries an HTML comment that
names its source file and line or JSON key. "Development" marks data used to choose model changes; "held-out"
marks data first analysed on frozen code under rules written before the run. Claims 4 and 5 were run once
on frozen code and filled in on 2026-11-04 and 2026-11-06. Reframed on 2026-11-06 after all prespecified claims
were run or withdrawn: 2 of 7 claims held, claim 8 was withdrawn. The paper is now a
failure-inclusive methods report, with bulk focal-loss detection as its supported result. It makes no
state-of-the-art claim.

## Abstract

Methods that infer copy number from RNA usually report one profile per cell or per group of cells. That profile
mixes three quantities: relative RNA dosage, DNA copy-number change and tumour clone membership. Each needs its
own DNA truth at its own resolution, and agreement at one level does not validate another. HaploDose keeps the
three apart. It segments pooled profiles with a Poisson-lognormal hidden Markov model against a fitted mixture
of reference types, assigns single cells to those profiles with a negative binomial likelihood that carries
depth, expression-programme and optional allele terms, and tests clone differences on genomic bins not used to
propose them. We fixed eight benchmark claims and their decision rules before seeing the data used to test them,
and report all of them.
On a held-out MIX-seq set of 127 cell lines <!-- docs/heldout_protocol.md:165 -->, HaploDose arm profiles were
more accurate than infercnvpy, CopyKAT and SCEVAN, tied with R inferCNV (Holm p 0.15
<!-- results/heldout_20261102_1424/analysis/summary.json e1.tests[b=infercnv_r].p_holm = 0.1529 -->) and were
below a pseudobulk baseline, so the prespecified profile claim failed. Peak memory was the lowest of all tools
from 50,000 cells up but not at 20,000 cells, so the prespecified memory claim also failed. In two independent
bulk RNA sets, focal-loss recall at a 5% false-positive rate met the prespecified rule against a naive gene
ratio and infercnvpy (CCLE, 822 lines; CPTAC-3, 618 tumours) <!-- results/ccle_bulk_20261102_1551/score/primary.json all.comparisons[0].n_lines; results/cptac_bulk_20261103_0618/score/primary.json all.comparisons[0].n_lines -->,
but the margin over the naive ratio was small and its bootstrap interval for the mean crossed zero in both.
Clone-difference recovery in ovarian tumours with DLP+ single-cell DNA (SPECTRUM, claim 4, 21 patients) failed:
mean per-patient AUROC was 0.647 for HaploDose and 0.601 to 0.646 for four comparators, and no Holm-corrected
paired test reached significance <!-- results/spectrum_run_20261103_2028/primary/primary.json -->. In the external benchmark of Schmid et al. 2025 (claim 5, 4 datasets), HaploDose had the worst mean rank of 10
outputs (8.42; best inferCNV expression output 3.33) <!-- results/schmid_claim5_rank_20261106_0003/summary_combined.tsv -->.
On development data this traced to bin-to-bin noise in the exported profile: arm-level accuracy was unaffected. A
smoothed export closed the gap to CopyKAT but stayed below a smoothed pseudobulk <!-- docs/STATUS.md:1324 -->. A third
bulk cohort (claim 7, papillary thyroid) had too few focal losses to evaluate, and a single-cell profile claim
(claim 8) was withdrawn before any data were read because development data predicted it would fail. Across
development data, the clone step over-split DNA-clonal lines and changed clone counts between random seeds.
HaploDose is therefore a competitive method for focal-loss detection in bulk RNA at a strict false-positive rate.
It is not the best available method for single-cell dosage profiles or clone identity. We report these failures
with the method.

## 1. Introduction

Copy-number change shifts the expression of the genes it covers, roughly in proportion to dosage. Tools such as
inferCNV, infercnvpy, CopyKAT, SCEVAN and Numbat use this to infer copy number and tumour subclones from RNA
[CITATION NEEDED: inferCNV; infercnvpy; CopyKAT; SCEVAN; Numbat]. Recent benchmarks compare these tools on
profile accuracy against DNA [CITATION NEEDED: Chang 2026; Song 2025; Schmid et al. 2025, Nat Commun 16:8777].

A single RNA profile answers three different questions, and each needs different evidence.

- Relative RNA dosage: how much a genomic region's expression is shifted relative to a reference. Bulk or
  pseudobulk DNA copy number relative to the same reference is the matching truth.
- DNA copy-number change: whether a region is gained, lost or allelically imbalanced. Discrete DNA calls
  (ABSOLUTE, ASCAT, FACETS) are the matching truth. Expression alone cannot see a uniform doubling of the genome,
  because every bin moves together (`python/haplodose/bulk.py:20`).
- Clone identity: which cells share a genotype. Only per-cell DNA, same-cell DNA and RNA, or known synthetic
  cell sets can test this. Aggregate DNA agreement does not.

HaploDose is built so that each output maps to one of these estimands. Following the structure of Wang et al.
(Cell Systems 2026, doi 10.1016/j.cels.2026.101762), we state the model's assumptions, write all components
in one notation, run targeted analyses for each estimand and list where the method fails. We also fixed the
benchmark claims and their decision rules in advance, so that a reader can see which results were tested once
on data never used to tune the method.

## 2. Methods

### 2.1 Estimands and outputs

| Estimand | HaploDose output | DNA truth | Resolution scored |
|---|---|---|---|
| Relative RNA dosage | `profiles.tsv` `log_ratio_gene`, `log_ratio`; `arms.tsv` `median_log_ratio_gene`; bulk log ratio and segmented dose | Line or tumour DNA copy number relative to the same reference | Arm (Spearman rho per line or tumour) |
| DNA copy-number change | `segments.tsv` states (loss, neutral, gain, amplification on a relative diploid scale); bulk arm calls; allele columns and `joint_state` | ABSOLUTE (DepMap 24Q4, TCGA), AscatNGS (CPTAC-3), FACETS (U-CAN) | ~10 Mb windows (focal loss); arms; per-tumour WGD |
| Clone identity | `cells.tsv.gz`, `posteriors.npy`, `clones.tsv`, `differences.tsv` | Synthetic implanted cell sets; per-arm clone differences from scDNA (Andor, SPECTRUM DLP+) | Line (ARI), patient (clone-difference AUROC) |

Segment states are relative copy number 1 to 4 on a diploid scale, not absolute copy number
(`python/haplodose/profile.py:1-7`). The label `neutral_allele_imbalanced` is consistent with copy-neutral LOH
but does not establish it (`python/haplodose/alleles.py` module docstring). `tree.tsv` is an average-linkage
similarity tree, not a phylogeny (`docs/scientific_code_review_2026-11-03.md:158`).

### 2.2 Notation

Cell $c$, gene $g$, bin $j = 1, \dots, J$, reference type $t$, clone $k$, programme $f = 1, \dots, F$, SNP $s$.
A bin holds 10 consecutive expressed genes by default <!-- python/haplodose/cli.py:25 --bin-genes default 10 -->;
one hidden Markov chain runs per chromosome. Clone $k = 0$ is copy-neutral relative to the reference baseline.
Defaults below are from `ProfileParams` (`python/haplodose/profile.py:34-41`) and `CloneParams`
(`python/haplodose/clones.py:82-113`). The frozen 14:24 code used for claims 1, 2, 3 and 6
(`results/heldout_20261102_1424/code/`) has the same `profile.py` constants and likelihood kernels as the
current repository; the differences are a guard for untestable pairs and device selection (checked by `diff`
against `python/haplodose/`).

### 2.3 Reference baseline

Each reference type has bin rates $\mu_{jt}$ with $\sum_j \mu_{jt} = 1$. The specimen baseline is
$\mu_j = \sum_t w_t \mu_{jt}$, with $w$ the Poisson maximum-likelihood mixture weights fitted to the pooled
specimen counts (`mixture_weights`, `python/haplodose/clones.py:311-318`). Per-bin heterogeneity between
reference types gives a mismatch variance $v_j$. Genes whose expression differs from every reference type by
more than a factor of 5 <!-- python/haplodose/pipeline.py:222 max_log_ratio=np.log(5.0) --> are removed before
binning, as are ribosomal protein genes and a user mask. In bulk mode, $\mu_j$ is the pooled normal panel and
$v_j$ the robust variance of the normals' log ratios (`python/haplodose/bulk.py:46-60`).

Assumption: the reference types span the non-malignant expression states in the specimen. If a tumour dosage
change lines up with a contrast between reference types, the mixture can absorb part of it
(`docs/scientific_code_review_2026-11-03.md:155`).

### 2.4 Pooled profile and segmentation

For a pooled group (a clone, a specimen or a bulk sample) with bin counts $Y_j$ and expected counts $E_j$:

$$Y_j \sim \mathrm{Poisson}(\lambda_j), \qquad \log \lambda_j = \log E_j + \text{centre} + \beta \log(n_j/2) + o_j, \qquad o_j \sim N(0, \tau^2 + v_j).$$

- $n_j \in \{1, 2, 3, 4\}$ is a hidden state (loss, neutral, gain, amplification) (`profile.py:23-27`).
- $\beta$ is the dosage slope. It combines dosage response and purity, and is not a purity estimate. Start 0.6,
  bounds 0.3 to 1.0 for clones <!-- python/haplodose/profile.py:36-37 -->; bounds 0.05 to 1.0 in bulk, where
  the start is the value on a grid (0.05 to 1.0) with the highest HMM likelihood
  <!-- python/haplodose/bulk.py:38-39 BULK_PROFILE, BETA_GRID -->.
- $\tau$ is the bin-level extra-Poisson standard deviation, estimated from the MAD of first differences within
  chromosomes after removing sampling variance, floor 0.02 <!-- python/haplodose/profile.py:41, 54-67 -->.
- The integral over $o_j$ uses adaptive Gauss-Hermite quadrature with 12 nodes centred at the posterior mode
  <!-- python/haplodose/kernels/likelihood.py:17 hermgauss(12) -->, written in a deviance form that avoids
  cancellation at high depth (`likelihood.py:37-69`).
- Transitions: $T_{ii} = 1 - p$, $T_{ij} = p\,q_j/(1 - q_i)$, with $p = 2 \times 10^{-3}$ and
  $q = (0.2, 0.5, 0.2, 0.1)$ <!-- python/haplodose/profile.py:27, 38; kernels/hmm.py:148-153 -->. The long-run
  neutral occupancy is therefore 0.379, not 0.5 <!-- docs/STATUS.md:710 -->.
- A segment is reported as resolved if it is non-neutral, spans at least 3 bins and has mean posterior support
  of at least 0.95 <!-- python/haplodose/profile.py:39-40, 154-155 -->.

The centre, $\beta$ and $\tau$ are refitted for 4 rounds from the pooled specimen (`fit_pooled`,
`profile.py:97-122`).

### 2.5 Depth

Each cell has centred log depth $x_c$, and each bin a depth slope $\gamma_j$, estimated from depth quintiles
within clones (or held fixed from reference cells with `--depth-source reference`;
`python/haplodose/clones.py:555-558`). For a pooled group with UMI-weighted depth mean $m$ and variance $v$,

$$E_j \propto \mu_j \exp\{(\gamma_j m + \gamma_j^2 v/2)/(1 + c v)\}, \qquad c = \textstyle\sum_j \mu_j \gamma_j^2$$

(`depth_log_rate`, `clones.py:560-571`). On simulated depths the worst bin error of this approximation was 2%,
against 24% for the mean depth alone <!-- python/haplodose/clones.py:566-567 docstring -->.

### 2.6 Expression programmes

Loadings $W \in \mathbb{R}^{F \times J}$, $F = 10$ <!-- python/haplodose/clones.py:100 n_programmes -->, are the
top singular vectors of residuals of fit reference cells. Each arm's mean is removed from every loading, so no
programme can explain a whole-arm change (`python/haplodose/programmes.py:1-17`). Cells receive scores
$u_{cf}$ shrunk towards 0 by 5 damped Fisher steps (`python/haplodose/kernels/estep.py:36-61`); pooled profiles
receive one free amplitude per programme. The converse limit is stated in the code: the arm-mean part of a
programme is left to dosage, so on short arms a strong programme can look like dosage
(`programmes.py:14-17`).

### 2.7 Cell likelihood

Clone $k$ has log rate $\ell_{kj}$ (baseline plus clone dose). For cell $c$ with binned total $L_c$:

$$\rho_{ckj} = \mathrm{softmax}_j\Big(\ell_{kj} + x_c\gamma_j + \textstyle\sum_f u_{cf} W_{fj}\Big), \qquad Y_{cj} \mid z_c = k \sim \mathrm{NB}(L_c \rho_{ckj}, \phi_j),$$

with variance $\mu + \phi_j \mu^2$ (`python/haplodose/kernels/likelihood.py:72-77`). Because the rates are a
softmax, only compositional differences between clones enter. Cell log-likelihoods are computed relative to
clone 0 from float64 log-rate differences, so small contrasts are not lost in float32
(`python/haplodose/kernels/estep.py:18-33`).

### 2.8 Alleles (optional)

Oriented counts $a_{cs}$, $d_{cs}$ come from population-phased SNPs. At clone level an HMM runs over a balanced
state ($b = 0.5$) and imbalanced states $b = 0.5 \pm \theta$, $\theta \in \{0.08, 0.17, 0.3, 0.45\}$
<!-- python/haplodose/alleles.py:44 LEVELS -->, with phase-switch probability $(1 - e^{-2\nu d})/2$ for $d$ in cM
(`alleles.py:166-175`) and beta-binomial emissions whose concentration $\kappa$ is a method-of-moments estimate
from same-patient reference cells, clipped to 5 to $10^4$ <!-- python/haplodose/alleles.py:141-163 -->. At cell
level the clone fractions are first shrunk towards 0.5,
$\tilde b_{ks} = 0.5 + 0.9\,(b_{ks} - 0.5)$ <!-- python/haplodose/alleles.py:57 cell_shrink, 302-304 -->, and the
cell adds $\sum_s a_{cs}\log\tilde b_{ks} + (d_{cs} - a_{cs})\log(1 - \tilde b_{ks})$. This binomial term does not
carry the beta-binomial overdispersion used at clone level (Limitations, F11).

When alleles are supplied, cells that are confidently balanced in the specimen's strong imbalance runs are added
as an in-sample normal reference by default (`alleles.py:453-471`; `CloneParams.in_sample_reference`,
`clones.py:97`). Limitation 10 explains why this is an assumption, not a validated step (F4).

### 2.9 Clone mixture and spatial prior

With mixture weights $\pi_k$ and summed cell log-likelihood $LL_{ck}$, the non-spatial E-step is
$R_{ck} \propto \pi_k \exp(LL_{ck})$. For spatial data, a mean-field update over a row-normalised
8-nearest-neighbour graph $G$ <!-- python/haplodose/clones.py:98 spatial_neighbors -->:

$$R_{ck} \propto \pi_k \exp\big(LL_{ck} + \lambda (G R)_{ck}\big).$$

$\lambda$ is chosen from $\{0, 0.5, 1, 2, 4, 8, 16\}$ <!-- python/haplodose/clones.py:519 --> to maximise how well
neighbours' data-only posteriors predict each cell's own likelihood, so spatially random clones give
$\lambda = 0$ (`clones.py:519-535`). $G$ is directed, so this is a smoothing heuristic, not the gradient of a
symmetric Potts energy (`neighbour_graph` docstring, `docs/STATUS.md:619`). Expression is never averaged across
neighbours for the cell likelihood.

### 2.10 Clone search and difference tests

Proposals come from Leiden clustering of windowed, depth-regressed log ratios, from two-component mixtures over
chromosome halves, and (with alleles) from imbalance runs (`clones.py:1-24`). By default proposals use even bins
and even SNPs; tests use odd bins and odd SNPs (`split_mode = 'evenodd'`, `clones.py:104`).

A difference between two pooled groups is tested after depth-matched pooling (coarsened exact matching on 8
log-depth strata, `clones.py:691-713`). The log-rate difference $d_j$ on test bins, centred by its weighted
median, enters a 3-state Gaussian HMM with levels $\{-1, 0, +1\} \times \beta \log 1.5$ and prior
$(0.15, 0.7, 0.15)$ <!-- python/haplodose/profile.py:30-31, 190 -->; the long-run neutral occupancy of this chain is
0.452 <!-- docs/STATUS.md:710 -->. A resolved difference segment must then pass a one-sided sign test that most of
its genes move with it (at least 8 genes, alpha 0.01; `clones.py:733-760`) and, with alleles, an allele veto
(`clones.py:715-731`). Tumour clones share one pooled tumour dose except on resolved pairwise differences
(`shared_profiles`, `clones.py:577`). Clones below 30 cells are dropped
<!-- python/haplodose/clones.py:83 min_cells -->.

**v-next split null (claim 4 code only).** The code frozen for claim 4 adds a Monte Carlo acceptance step for
divisive splits: a clone's best split must beat the best split found by the same proposal-and-test search on 19
simulations from the fitted one-clone model <!-- results/vnext_freeze_20261103_0744/README.md:3 -->. In the
`cov` model each simulated cell gets a log-normal rate effect per chromosome half drawn from the clone's own
moment-matched region covariance (`clones.py:1199-1255`). The null does not gate the initial groups in this
configuration (`split_null_initial` false; `clones.py:111`). It is a model-based split filter, not a
calibrated specimen-level false-split rate (F2).

### 2.11 Gene-level output

`log_ratio_gene` is the median over a bin's genes of $\log((y_g + 0.5)/(E_g + 0.5))$, where $E_g$ splits the
bin expectation by gene shares from the fitted reference mixture; genes with $E_g < 1$ are left out
<!-- python/haplodose/pipeline.py:259-272 min_expected=1.0 -->. This column was chosen after a decomposition on
Kinker development data <!-- docs/STATUS.md:500-506 -->, and was fixed as the claim 1 endpoint before the
held-out run <!-- docs/heldout_protocol.md:49-50 -->.

### 2.12 Bulk mode

Each sample is one pooled profile against a normal panel (`python/haplodose/bulk.py:1-20`). Programmes are
singular vectors of the normals' log ratios, high-passed within arms. Arm calls are z scores against
leave-one-out normals, scaled by the sample's $\tau$, called at $|z| > 3$
<!-- python/haplodose/bulk.py:40 ARM_Z -->. With fewer than 10 normals the run records that the null scale is not
calibrated <!-- python/haplodose/bulk.py:45 SMALL_PANEL -->; that flag reflects panel size, not observed
calibration (F11). Reported log ratios are centred by the plain median over bins.

### 2.13 Comparators

| Method | Reference model | Error model | Segmentation | Clone step as benchmarked | Needs reads |
|---|---|---|---|---|---|
| HaploDose | Fitted mixture of reference types | PLN pooled bins; NB per cell; beta-binomial alleles | 4-state HMM on pooled profiles | Proposals, EM, cross-fit split and merge tests | Optional (alleles) |
| infercnvpy | Reference cell mean [Unverified] | None documented [Unverified] | 100-gene running window | Tool output plus our malignant and clone rule | No |
| R inferCNV | Reference cell mean [Unverified] | i6 HMM, Bayesian mixture filter | HMM per subcluster | Leiden subclusters plus our rule | No |
| CopyKAT | Diploid cells inferred or supplied [Unverified] | Poisson-gamma MCMC per window | MCMC windows | No subclone output; our Leiden post-processing | No |
| SCEVAN | Normal cells inferred or supplied [Unverified] | [Unverified] | Mumford-Shah | Native clones | No |
| Numbat | Reference mixture by constrained least squares | PLN expression; binomial alleles | Joint expression-allele HMM | Iterative phylogeny | Yes |
| Pseudobulk / naive gene ratio | Reference pool / pooled normals | None | None; arm median of gene ratios | None | No |

Rows other than HaploDose are taken from `docs/manuscript_outline.md:127-135`, which cites
`docs/related_work.md`; cells marked [Unverified] were not confirmed against the tools' papers or source. Tool
versions for the field comparators: CopyKAT 1.2.5 and SCEVAN 1.0.3 in env `cnv_field`, R 4.4
<!-- docs/STATUS.md:542-543 -->. In every Kinker and MIX-seq run, all tools received other cancer lines as the
reference <!-- docs/STATUS.md:526-527 -->. None of these tools was built for that.

Pseudobulk and the naive gene ratio are baselines, not tools. For a single-line specimen, which is one
population, pseudobulk against the same reference is close to the best possible profile
<!-- docs/STATUS.md:408-409 -->.

### 2.14 Benchmarks, prespecification and freezing

**Development data.** Kinker 2020, Andor 2020, U-CAN (bulk and Visium HD), TCGA, Moorman/HTAN, SMC and YCC were
used to choose model changes. Their results are development evidence and are expected to be optimistic
(`docs/heldout_protocol.md:5-8`).

**Prespecified claims.** On 2026-11-02 the user fixed SOTA target claims and decision rules before any held-out
output existed <!-- docs/STATUS.md:54-58 -->. A claim holds only if its rule is met on data not used in
development, and each is reported whatever the outcome. Claim 6 was added on 2026-11-03 before any CPTAC data
was downloaded <!-- docs/STATUS.md:863-867 -->; claim 7 on 2026-11-06 before any REBC-THYR data was downloaded <!-- docs/STATUS.md:1296 -->;
claim 8 on 2026-11-06 before any Tahoe-100M expression was read, and withdrawn the same day <!-- docs/STATUS.md:1395, 1446 -->.
Claims 7 and 8 were added after the earlier results were known, to seek a third win. Readers should weigh the
count of wins with that in mind.

| Claim | Data (role) | Primary endpoint and rule | Source |
|---|---|---|---|
| 1 Profile accuracy | MIX-seq, 127 lines (held-out) | Per-line arm rho; paired Wilcoxon, Holm over 5 (4 tools plus pseudobulk); superior if Holm p < 0.05 and 95% CI > 0 against all 4 tools; plus common-gene-universe CI > 0 (amendment 5) | `docs/STATUS.md:60-65`; `docs/heldout_protocol.md:47-54, 126-143` |
| 2 Peak memory | Kinker cells resampled, 5k to 1M (controlled) | Lowest median peak memory among finishing tools at every size from 20k cells | `docs/STATUS.md:66-69` |
| 3 Focal loss, cell-line bulk | CCLE vs DepMap 24Q4 ABSOLUTE (held-out) | Focal-loss recall at 5% realised FPR; paired Wilcoxon over lines; Holm over 2; both Holm p < 0.05 with HaploDose higher | `docs/STATUS.md:70-84` |
| 4 Clone differences | SPECTRUM ovarian, DLP+ scDNA (held-out, blind) | Per-patient clone-difference AUROC; paired Wilcoxon over 21 patients; Holm over 4; every Holm p < 0.05 and bootstrap CI of mean > 0; each comparator scored for >= 50% of patients | `docs/STATUS.md:94-105, 952-975, 997-1013` |
| 5 External benchmark | Schmid et al. 2025 pipeline, 4 primary datasets (blind) | Best mean rank over 3 metrics x 4 datasets against 9 published outputs; descriptive | `docs/STATUS.md:121-153` |
| 6 Focal loss, tumour bulk | CPTAC-3 vs AscatNGS (held-out) | As claim 3, unit = tumour; both Holm p < 0.05 and mean diff > 0 | `docs/STATUS.md:868-897` |
| 7 Focal loss, thyroid bulk | REBC-THYR vs AscatNGS (held-out) | As claim 6; not evaluable if < 20 tumours have a focal loss | <!-- docs/STATUS.md:1296 --> |
| 8 1 Mb profile, pooled lines | Tahoe-100M plate1 DMSO (held-out; withdrawn before data) | Per-line 1 Mb Pearson; smoothed HaploDose vs CopyKAT, infercnvpy, smoothed pseudobulk; Holm over 3 | <!-- docs/STATUS.md:1395 --> |

**Freezing.** The method for claims 1, 2, 3 and 6 is the package as of 2026-11-02 14:24, copied to
`results/heldout_20261102_1424/code/` with SHA-256 of every file; every held-out job checks the hashes
(`docs/heldout_protocol.md:149-160`). Evaluation code was hashed before data were opened for claim 3
(`results/ccle_bulk_20261102_1551/README.md:1-5`) and claim 6 <!-- docs/STATUS.md:943 -->. Claim 4 uses a later
freeze, `results/vnext_freeze_20261103_0744/`, adopted even though it failed two of five development gate
criteria <!-- docs/STATUS.md:930-937 -->. Its scorer was validated on Andor development data and hashed before any
SPECTRUM RNA output was read (`results/spectrum_scorer_20261103_0613/README.md`), then amended twice before any
CNV tool ran on SPECTRUM RNA <!-- docs/STATUS.md:997-1022 -->. Every amendment is dated and listed with its
reason in `docs/heldout_protocol.md:70-143` and `docs/STATUS.md`.

**MIX-seq design (claim 1).** Control pools from McFarland et al. 2020 (Nat Commun 11:4296), one control sample
per pool composition (`docs/heldout_protocol.md:72-78`). Within each pool, lines are split in half at random;
one half is the reference and each line of the other half is a single-line specimen, then the halves swap
(`docs/heldout_protocol.md:38-45`). Truth per gene is line DepMap copy number minus the UMI-weighted mean of the
pool's fit reference lines; arm truth is the median over the run's used genes. 127 lines met the rules; 44 were
also in Kinker, so the prespecified sensitivity set has 83 lines <!-- docs/heldout_protocol.md:147 -->. E2
implants one event (50, 150 or 500 genes; gain or loss; 15% or 40% of cells) by binomial thinning in each line
with at least 100 cells, 47 specimens <!-- docs/heldout_protocol.md:119 -->.

**Bulk focal-loss design (claims 3 and 6).** Gene-anchored windows of 10 Mb with at least 10 genes; thresholds
set at a 5% realised FPR on DNA-neutral windows of the same group; random tie-break, seed 20261101
<!-- results/ccle_bulk_20261102_1551/score/params.json window_bp, min_genes, fpr, seed -->. A focal window is a
segment shorter than 25% of its arm with $|TCN - \text{ploidy}| \ge 1$; DNA-neutral is
$|TCN - \text{ploidy}| < 0.5$ <!-- docs/STATUS.md:875-876 -->. Claim 3 normal panel: GTEx v10 tissue-matched
normals; truth DepMap 24Q4 ABSOLUTE segments <!-- docs/STATUS.md:74-80 -->. Claim 6 normal panel: every Solid
Tissue Normal file of the CPTAC-3 cohort, at least 10 per cohort; truth GDC AscatNGS segments
<!-- docs/STATUS.md:874-886 -->. The same panel is used for all three methods in each claim.

**SPECTRUM design (claim 4).** Truth per arm is 1 if any two DLP+ tumour clones differ on at least 30% of the
arm's bins by at least one copy; a tool's score per arm is the largest absolute difference between any two of
its RNA tumour clones <!-- docs/STATUS.md:96-100 -->. 28 of 39 patients with clone tables were DNA-eligible
<!-- docs/STATUS.md:902-903 -->; 21 also have unsorted or CD45-negative tumour scRNA <!-- docs/STATUS.md:956-958 -->.
Reference cells are the authors' published non-cancer labels, up to 300 per type, the same for every tool
<!-- docs/STATUS.md:960-963 -->. A tool that fails for a patient scores AUROC 0.5 <!-- docs/STATUS.md:970-971 -->.

### 2.15 Statistics

The unit of replication is the cell line (Kinker, MIX-seq, CCLE), the tumour or patient (U-CAN, TCGA, CPTAC-3,
SPECTRUM) or the run (scaling). Cells are never the unit of a test. Paired comparisons use the two-sided
Wilcoxon signed-rank test. Each family for Holm adjustment is named next to its result. Raw p values are reported
from the JSON key `p` and adjusted values from `p_holm`; "Holm p" in the text always means `p_holm`. Effect
intervals are 95% percentile bootstrap intervals over the biological unit (10,000 resamples; seed 20261103 for
MIX-seq and CPTAC-3) <!-- docs/heldout_protocol.md:52-53; docs/STATUS.md:889 -->. For MIX-seq E1 the protocol
interval is for the median paired difference; for claims 3 and 6 it is for the mean. A signed-rank test and an
interval for the mean answer different questions, and the two can disagree (F10). No model in these tests has
covariates; tests are paired within the biological unit. MIX-seq lines share 6 reference sets
(3 pools x 2 rotations), so line-level tests assume an independence that does not hold; we report a cluster
sign-flip test and a mixed model as descriptive checks <!-- docs/heldout_protocol.md:138-140 -->. Seeds are fixed
in each run directory. Environments are recorded as `pip_freeze*.txt` and `sessionInfo_cnv_field.txt` in each
run directory.

## 3. Results

### 3.1 Evidence matrix

| Dataset (unit, n) | Role | Code | Relative RNA dosage | DNA copy-number change | Clone identity |
|---|---|---|---|---|---|
| MIX-seq (line, 127) | Held-out | Frozen 14:24 | Claim 1: fails (3.2) | Not scored | E2 implant: no superiority (3.2) |
| CCLE bulk (line, 822 with focal loss) | Held-out | Frozen 14:24 | n/a | Claim 3: holds, small margin over naive (3.4) | n/a |
| CPTAC-3 bulk (tumour, 618 with focal loss) | Held-out | Frozen 14:24 | n/a | Claim 6: holds, small margin over naive (3.5) | n/a |
| REBC-THYR bulk (tumour, 6 with focal loss) | Held-out | Frozen 14:24 | n/a | Claim 7: not evaluable, below the prespecified floor of 20 (3.5) | n/a |
| Tahoe-100M (line) | Held-out; withdrawn before data | Frozen 14:24 + 5-bin export | Claim 8: withdrawn (development predicted a loss to smoothed pseudobulk) | n/a | n/a |
| SPECTRUM (patient, 21) | Held-out, blind | v-next freeze | n/a | n/a | Claim 4: fails, parity with CopyKAT and SCEVAN (3.6) |
| Schmid 2025 (dataset, 4 primary) | External, blind | v-next freeze (14:24 as sensitivity) <!-- docs/STATUS.md:1034 --> | Claim 5: fails, last of 10 on mean rank (3.7) | Claim 5 gain/loss AUC: included in the mean rank (3.7) | n/a |
| Scaling (run, 5k to 1M cells) | Controlled compute | Frozen 14:24 | n/a | n/a | Claim 2: fails as prespecified (3.3) |
| Kinker 2020 (line, 181 single-line; 166 implant) | Development | Repo 2026-11-02 05:29 | Median rho 0.861 vs pseudobulk 0.863 | Not scored | Implant ARI 0.167, tie with infercnvpy |
| Andor 2020 (line; SNU-601 n = 1 full-depth scDNA) | Development | Frozen 14:24 | Not re-scored here | Not scored | Clone-difference AUROC 0.78 to 0.84 by seed; over-splitting of DNA-clonal lines |
| U-CAN bulk (patient, 1,034) | Development | Frozen 14:24 | Dose 0.649 vs naive 0.630 | LOH AUC 0.93; WGD AUC 0.81 to 0.89 (older code) | n/a |
| TCGA bulk (tumour, 6,243) | Development | 2026-11-01 code | Negative: naive higher | Focal-loss recall at 5% FPR 0.531 vs naive 0.472 | n/a |
| U-CAN Visium HD (patient, 20) | Development, secondary | Several | Median rho 0.62 to 0.75 by reference | Arm calls (3.7) | Negative: subclone overlap not above chance |
| Moorman / SMC normals (patient) | Development | Frozen 14:24 | n/a | n/a | 0 false splits in held-out normal cells (easy null) |

Sources for each cell are given in the subsections below.

### 3.2 Claim 1: arm-profile accuracy on held-out MIX-seq (fails)

All 10 tool jobs completed and every tool returned output for all 127 lines
<!-- docs/heldout_protocol.md:164; results/heldout_20261102_1424/analysis/summary.json e3.lines_without_output -->.
Median per-line arm Spearman rho was 0.832 for HaploDose `log_ratio_gene`, 0.840 for pseudobulk, 0.823 for R
inferCNV, 0.799 for infercnvpy, 0.754 for SCEVAN and 0.699 for CopyKAT
<!-- results/heldout_20261102_1424/analysis/summary.json e1.median_rho -->.

| HaploDose vs | n lines | Median paired diff | 95% CI (median) | HaploDose higher / lower | p | Holm p (family of 5) |
|---|---|---|---|---|---|---|
| infercnvpy | 127 | +0.027 | 0.017 to 0.045 | 91 / 36 | 5.5e-6 | 1.6e-5 |
| R inferCNV | 127 | +0.011 | -0.002 to 0.022 | 73 / 54 | 0.153 | 0.153 |
| CopyKAT | 127 | +0.115 | 0.091 to 0.136 | 117 / 10 | 2.6e-20 | 1.3e-19 |
| SCEVAN | 127 | +0.076 | 0.062 to 0.089 | 114 / 13 | 7.3e-18 | 2.9e-17 |
| pseudobulk | 127 | -0.010 | -0.015 to -0.0002 | 52 / 75 | 0.0028 | 0.0056 |

<!-- all rows: results/heldout_20261102_1424/analysis/summary.json e1.tests[*].{n, median_diff, ci95, a_higher, b_higher, p, p_holm} -->

HaploDose met the superiority rule against infercnvpy, CopyKAT and SCEVAN, conditional on independence of lines.
It did not against R inferCNV, and it was below pseudobulk. The common-gene-universe sensitivity required by
amendment 5 gave the same pattern: against R inferCNV the mean difference was +0.006 with a 95% CI of -0.004 to
0.016 <!-- results/heldout_20261102_1424/analysis/sensitivity/summary.json a_e1_common_universe.tests[b=infercnv_r].mean_diff, mean_ci95 -->,
and pseudobulk was further ahead (median difference -0.019, Holm p 3.7e-8)
<!-- results/heldout_20261102_1424/analysis/sensitivity/summary.json a_e1_common_universe.tests[b=pseudobulk].median_diff, p_holm -->.
The decision recorded by the analysis is `claim_holds: false`
<!-- results/heldout_20261102_1424/analysis/sensitivity/summary.json claim1_decision.claim_holds -->.

Sensitivity and dependence checks, all descriptive:

- 83 lines never in Kinker: same pattern; R inferCNV Holm p 0.18 and pseudobulk Holm p 0.18
  <!-- results/heldout_20261102_1424/analysis/summary.json e1_sensitivity_not_in_kinker.tests[b=infercnv_r].p_holm = 0.179, [b=pseudobulk].p_holm = 0.179 -->.
- Cluster sign-flip over the 6 reference sets: p 0.031 (the smallest attainable) against infercnvpy, and 0.97
  against R inferCNV <!-- results/heldout_20261102_1424/analysis/sensitivity/summary.json e1_primary.dependence[*].cluster_signflip.p -->.
  A mixed model with a random intercept per reference set gave an intercept of +0.022 (95% CI 0.004 to 0.040)
  against infercnvpy and +0.006 (-0.007 to 0.019) against R inferCNV
  <!-- results/heldout_20261102_1424/analysis/sensitivity/summary.json e1_primary.dependence[*].mixed_model.intercept, ci95 -->.
  The mixed model reported an estimate on the boundary of the parameter space for both.

The secondary implant endpoint (E2, n = 47 lines) did not favour HaploDose. Mean `ari_all` was 0.136 for
HaploDose, 0.223 for R inferCNV, 0.185 for infercnvpy, 0.137 for CopyKAT and 0.001 for SCEVAN
<!-- results/heldout_20261102_1424/analysis/summary.json e2.mean_ari_all -->. With Holm over 4, no comparison met the
superiority rule; against SCEVAN Holm p was 0.027 but the CI of the median difference had a lower bound of 0
<!-- results/heldout_20261102_1424/analysis/summary.json e2.tests[b=scevan].p_holm, ci95 -->. Below 30 implanted
cells HaploDose's mean ARI was 0.0006, against 0.156 for infercnvpy and 0.185 for R inferCNV
<!-- results/heldout_20261102_1424/analysis/summary.json e2.by_floor -->, which reflects the 30-cell minimum clone size.
HaploDose called one tumour clone in 93 of 127 single-line specimens, two in 30 and three in 4
<!-- results/heldout_20261102_1424/analysis/summary.json e3.tumour_clones.haplodose.distribution -->.

Interpretation. On held-out lines, HaploDose's arm dosage was as accurate as R inferCNV and more accurate than
three other tools, but it was not better than a plain pseudobulk ratio against the same reference. The
single-line design favours pseudobulk, since each specimen is one population. The comparison that would show
whether HaploDose's model adds anything over pseudobulk is a mixed specimen with a diploid normal reference,
which no held-out benchmark here has.

### 3.3 Claim 2: peak memory at scale (fails as prespecified)

Peak memory is the Slurm MaxRSS, median of 3 repeats (1 at 1M cells), on v6 CPU nodes with 16 CPUs; unit = run
<!-- docs/STATUS.md:833-835 -->. HaploDose runs used the frozen 14:24 code; comparator runs came from an earlier
batch with unchanged comparator code (`results/scaling_combined_20261103_0102/README.md`).

| Cells | HaploDose GiB (wall min) | infercnvpy | R inferCNV | CopyKAT | SCEVAN |
|---|---|---|---|---|---|
| 20k | 2.47 (1.73) | 2.44 (0.42) | 11.67 (13.95) | 64.74 (114.48) | 128.01 (53.57) |
| 50k | 3.52 (3.90) | 4.67 (0.54) | 28.86 (65.72) | 171.0 (646.6) | out of memory at 200 GB |
| 100k | 4.50 (7.65) | 9.28 (1.27) | 63.72 (245.07) | did not finish | out of memory |
| 1M | 11.06 (59.06) | 103.51 (4.47) | not run | not run | not run |

<!-- all cells: results/scaling_combined_20261103_0102/summary.tsv peak_rss_gib, wall_min; failures: not_finished.tsv; "not run" and SCEVAN 50k/100k OOM: docs/STATUS.md:840-842 -->

At 20k cells infercnvpy's median peak was 0.025 GiB lower than HaploDose's
<!-- results/scaling_combined_20261103_0102/runs.tsv, medians of peak_rss_gib: haplodose 2.4653, infercnvpy 2.4403 -->.
The per-run HaploDose values were 2.27, 2.47 and 2.47 GiB, and infercnvpy 2.44, 2.44 and 2.42 GiB
<!-- results/scaling_combined_20261103_0102/runs.tsv cells=20000 -->. The gap is within HaploDose's run-to-run spread,
but the rule compared medians and we do not relax it after the fact. HaploDose had the lowest median peak from
50k cells up, about 9 times lower than infercnvpy at 1M cells. infercnvpy was faster at every size, for example
4.47 vs 59.06 min at 1M cells. The GPU build of HaploDose used more host memory than the CPU build (12.86 GiB at
1M) <!-- results/scaling_combined_20261103_0102/summary.tsv haplodose_gpu 1000000 -->. Larger inputs contain
resampled duplicate cells, so this tests computation, not biological diversity
(`docs/scientific_code_review_2026-11-03.md:180`).

### 3.4 Claim 3: focal-loss detection in CCLE bulk RNA (holds; small margin over naive)

1,049 lines in 18 lineage groups were scored <!-- results/ccle_bulk_20261102_1551/score/type_info.tsv, sum of n_tumours over 18 rows -->;
822 had at least one focal-loss window and enter the test (unit = line)
<!-- results/ccle_bulk_20261102_1551/score/primary.json all.comparisons[*].n_lines -->.

| HaploDose dose vs | Mean recall HaploDose / comparator | Mean diff | Bootstrap 95% CI (mean) | Higher / lower / tie | p | Holm p (family of 2) |
|---|---|---|---|---|---|---|
| naive gene ratio | 0.509 / 0.484 | +0.024 | -0.001 to 0.049 | 354 / 268 / 200 | 0.0032 | 0.0032 |
| infercnvpy | 0.509 / 0.423 | +0.085 | 0.062 to 0.108 | 361 / 209 / 252 | 8.2e-15 | 1.6e-14 |

<!-- results/ccle_bulk_20261102_1551/score/primary.json all.comparisons[*].{mean_hd_dose, mean_comparator, mean_diff, ci_low, ci_high, hd_higher, comparator_higher, ties, p, p_holm}; claim_holds = true -->

The claim holds by its rule. The interval for the mean difference against the naive ratio crosses 0, so the
advantage there is small and the signed-rank result does not establish a positive mean. Descriptive subsets:

- Excluding lines whose DepMap copy number was used in development (n = 646): same direction; vs naive Holm p
  0.011, CI of the mean -0.005 to 0.053
  <!-- results/ccle_bulk_20261102_1551/score/primary.json excluding_dev_lines.comparisons[0] -->.
- Unstranded RNA protocol (n = 554 lines with focal loss): HaploDose ahead of naive (+0.047, Holm p 4.8e-5)
  <!-- results/ccle_bulk_20261102_1551/score/primary.json by_stranded.False.comparisons[0] -->. Stranded protocol
  (n = 268): no advantage over naive (-0.023, Holm p 0.69) or infercnvpy (+0.025, Holm p 0.69)
  <!-- results/ccle_bulk_20261102_1551/score/primary.json by_stranded.True.comparisons[*] -->. The GTEx normals are
  unstranded <!-- results/ccle_bulk_20261102_1551/README.md:8 -->, so the margin depends on protocol match.

### 3.5 Claim 6: focal-loss detection in CPTAC-3 tumour bulk RNA (holds; small margin over naive)

1,553 tumours in 9 cohorts were scored <!-- docs/STATUS.md:981-983; results/cptac_bulk_20261103_0618/score/type_info.tsv sum of n_tumours -->;
618 had at least one focal-loss window (unit = CPTAC-3 case)
<!-- results/cptac_bulk_20261103_0618/score/primary.json all.comparisons[*].n_lines -->.

| HaploDose dose vs | Mean recall HaploDose / comparator | Mean diff | Bootstrap 95% CI (mean) | Higher / lower / tie | p | Holm p (family of 2) |
|---|---|---|---|---|---|---|
| naive gene ratio | 0.542 / 0.514 | +0.028 | -0.003 to 0.060 | 229 / 184 / 205 | 0.025 | 0.025 |
| infercnvpy | 0.542 / 0.407 | +0.136 | 0.109 to 0.163 | 279 / 107 / 232 | 2.1e-21 | 4.3e-21 |

<!-- results/cptac_bulk_20261103_0618/score/primary.json all.comparisons[*]; claim_holds = true -->

The claim holds by its rule (both Holm p < 0.05, both mean differences > 0). As for CCLE, the bootstrap interval
for the mean against the naive ratio crosses 0. Excluding pooled-aliquot RNA (n = 508), the direction was the
same (vs naive Holm p 0.033; vs infercnvpy Holm p 1.5e-18)
<!-- results/cptac_bulk_20261103_0618/score/primary.json excluding_pooled_rna.comparisons[*].p_holm -->.

Per-cohort results are descriptive and heterogeneous. Against the naive ratio HaploDose was higher in CM (+0.174,
n = 95) and lower in PDAC (-0.263, n = 20) and HNSCC (-0.101, n = 51)
<!-- results/cptac_bulk_20261103_0618/score/primary.json by_cohort.{CM,PDAC,HNSCC}.comparisons[0].mean_diff, n_lines -->.
In CCRCC the difference was -0.011 (n = 71) <!-- results/cptac_bulk_20261103_0618/score/primary.json by_cohort.CCRCC.comparisons[0] -->.
Claims 3 and 6 test the same endpoint in two cohorts. Together they replicate one bulk task; they do not show a
second capability (`docs/scientific_code_review_2026-11-03.md:196`).

Claim 7 repeated the claim 6 rule on REBC-THYR, the only other GDC cohort with open AscatNGS and at least 10
solid-normal RNA samples <!-- docs/STATUS.md:1296 -->. It was prespecified with a power floor: at least 20 tumours with a
focal-loss window. One amendment was made before any score existed: RNA was matched to DNA by case, because REBC RNA
sample IDs never match the DNA sample IDs <!-- docs/STATUS.md:1354 -->. Of 372 tumours, only 6 had a focal loss
(median ploidy 2.00). The claim is therefore not evaluable and counts as not won <!-- docs/STATUS.md:1368 -->. Papillary
thyroid cancers are close to copy-number neutral, so this cohort cannot test focal-loss detection.

### 3.6 Claim 4: clone-difference recovery in SPECTRUM ovarian tumours (blind)

Design as prespecified (Section 2.14). Code: v-next freeze with flags `--split-mode evenodd --split-null-draws 19
--split-null-model cov`, seed 20261031, no alleles <!-- docs/STATUS.md:966-967 -->. Comparators: infercnvpy,
CopyKAT, SCEVAN, Numbat <!-- docs/STATUS.md:1002 -->. Primary: paired Wilcoxon over 21 patients, Holm over 4;
claim holds only if every Holm p < 0.05, every bootstrap 95% CI of the mean difference is above 0, and each
comparator scored at least 50% of patients <!-- docs/STATUS.md:972-973, 1004-1005 -->. Descriptive: the same tests
over patients where both tools scored, a common-arm sensitivity (amendment 2) and the count of non-finite
posteriors in each HaploDose `cells.tsv.gz` (review F1) <!-- docs/STATUS.md:1006, 1017-1022 -->.

- Primary result: claim 4 fails <!-- results/spectrum_run_20261103_2028/primary/primary.json -->. Mean AUROC: HaploDose 0.647, CopyKAT 0.646, SCEVAN 0.645, Numbat
  0.626, infercnvpy 0.601. HaploDose minus comparator (mean difference, bootstrap 95% CI, raw p, Holm p over 4):
  infercnvpy +0.047 (-0.003 to 0.097), p 0.100, Holm 0.402; CopyKAT +0.001 (-0.043 to 0.048), p 0.973, Holm 1.0;
  SCEVAN +0.003 (-0.115 to 0.100), p 0.374, Holm 1.0; Numbat +0.021 (-0.025 to 0.070), p 0.573, Holm 1.0.
- Completion per tool: every comparator passed the 50% rule. Patients scored 0.5 for missing output: infercnvpy 1
  (OV-031), SCEVAN 1 (OV-009, out of memory), Numbat 1 (OV-031), CopyKAT 0, HaploDose 0.
- Both-scored and common-arm sensitivities: same conclusion; the smallest Holm p was 0.685 (infercnvpy,
  both-scored, n = 20). HaploDose AUROC was below 0.5 in OV-002 (0.305) and OV-050 (0.047); OV-050 has one
  DNA-subclonal arm of 44, so its AUROC rests on a single positive arm.
- Non-finite HaploDose posteriors per patient: 0 in all 21 patients, for both the primary and the `--alleles` runs.
  `assignment_margin` is empty in 5 single-clone specimens, where it is undefined.
- HaploDose with `--alleles` (descriptive, allele path not gated, no test): mean AUROC 0.661 vs 0.647 without;
  higher in 11 of 21 patients, lower in 9.

Development context, not evidence for claim 4: on Andor SNU-601 the per-patient clone-difference AUROC was 0.779
for frozen HaploDose, 0.845 for R inferCNV and 0.766 for infercnvpy
<!-- results/spectrum_scorer_20261103_0613/README.md:10-12 -->. The v-next code passed the SNU-601 gate
(median AUROC 0.818 over 5 seeds) <!-- results/vnext_gate2_20261103_0612/analysis/summary.json cov.c_median -->
but failed seed stability and the MKN-45 false-split gate (Section 3.9). This endpoint scores which arms differ
among clones within a patient. It does not test which cells belong to which clone, so even a positive result
would not validate per-cell clone labels (`docs/scientific_code_review_2026-11-03.md:159`).

### 3.7 Claim 5: external benchmark of Schmid et al. 2025 (blind)

Design as prespecified: HaploDose run through the authors' pipeline (commit 319ea34) next to the 9 outputs the
paper ranks (inferCNV, SCEVAN and Numbat each as CNV and expression output, plus copyKat, CONICSmat and CaSpER),
all comparators rerun by us <!-- docs/STATUS.md:121-133 -->. Primary datasets never used in our development:
COLO320-HSR, BCC06, BCC06post, MM199 (n = 4) <!-- docs/STATUS.md:124-125 -->. Rule: best mean rank over Pearson,
AUC gains and AUC losses across the 4 datasets; with n = 4 no test can pass Holm, so the rule is descriptive
<!-- docs/STATUS.md:127-133 -->. Gate before HaploDose is scored: each rerun comparator output within 0.05 of the
published value on all three metrics; outputs that miss are still ranked and flagged
<!-- docs/STATUS.md:146-149 -->. A375 and HCT116 are plate-based full-length data (DNTR-seq), an untested modality
for HaploDose, and are reported separately with MCF7 as a line-overlap set <!-- docs/STATUS.md:125-126, 138-139 -->.

- Comparator reproduction gate: BCC06post 27/27 rows within 0.05; BCC06 26/27; COLO320 26/27 (CaSpER AUC losses,
  gap 0.052); MM199 control {8} 21/27 and {5,6,8} 25/27; MCF7 21/21 non-Numbat rows, Numbat failed (3 normal cells)
  <!-- docs/STATUS.md claim 5 gate entries -->. The gate needed one fully passing primary dataset (BCC06post); misses
  were ranked and flagged.
- Primary mean rank, 3 metrics x 4 datasets: claim 5 fails <!-- results/schmid_claim5_rank_20261106_0003/summary_combined.tsv -->. inferCNV (Expr) 3.33, Numbat (Expr) 4.04,
  copyKat 4.12, Numbat (CNV) 4.25, inferCNV (CNV) 5.25, CaSpER 5.42, CONICSmat 6.25, SCEVAN (CNV) 6.83, SCEVAN (Expr)
  7.08, HaploDose 8.42. HaploDose Pearson with DNA truth was 0.14 to 0.28, about half of copyKat's on three of the
  four datasets; it ranked last on all three metrics in BCC06post.
- Normalised overall score of the paper (descriptive): HaploDose 0.246, lowest; comparators 0.278 to 0.364.
- Sensitivities: per-method scorer 8.00 (9th of 10); MM199 control {5,6,8} 8.50; 14:24 code 8.42. MCF7 line overlap:
  7.67, 8th of 10 (Numbat outputs missing and ranked last). HaploDose gene coordinates were checked against the
  scorer's hg38 annotation; Pearson and AUC do not depend on scale.

Recorded deviations: CellRanger 8.0.1 instead of 7.0.0; a locally built chr1-22 Numbat panel; external normal
references chosen to match the authors' unnamed files <!-- docs/STATUS.md:140-153 -->. During transcription of
published comparator values, an agent printed the first 4 rows of the paper's published truth-track Source Data;
no HaploDose run or tuning followed <!-- docs/STATUS.md:150-152 -->.

### 3.8 Development evidence: relative RNA dosage

**Kinker single-line specimens (181 lines; repo code of 2026-11-02 05:29; Holm over 5 pairs with HaploDose).**
Median arm rho: HaploDose `log_ratio_gene` 0.861, pseudobulk 0.863, R inferCNV 0.835, infercnvpy 0.813, CopyKAT
0.765, SCEVAN 0.710 <!-- results/kinker_field_20261102_0529/single_line_with_infercnv_r/summary.json median_rho -->.
HaploDose was higher than each published tool (Holm p 1.4e-6 vs R inferCNV and smaller for the others) and lower
than pseudobulk (74 vs 106 lines, Holm p 0.032)
<!-- results/kinker_field_20261102_0529/single_line_ci/summary.json paired_tests[*].p_holm, a_higher, b_higher -->.
The `log_ratio_gene` column was chosen on these data, so they are not independent evidence for it
<!-- docs/STATUS.md:505-506 -->. Arm dosage was stable across 5 seeds, two CPU generations and a GPU run: median
per-line rho 0.853 to 0.856 across 8 runs <!-- docs/STATUS.md:170 -->.

**U-CAN bulk (1,034 tumours; frozen 14:24 code; Holm over 6 method pairs).** Median per-tumour arm rho against
FACETS: HaploDose segmented dose 0.649, HaploDose median-centred log ratio 0.627, naive gene ratio 0.630,
infercnvpy 0.566 <!-- results/ucan_bulk_frozen1424_20261102_1612/score_20261102_162020/summary.json median_rho -->.
Dose was higher than naive (587 of 1,034 tumours, Holm p 1.1e-9); the log ratio did not differ from naive (Holm p
0.81) <!-- results/ucan_bulk_frozen1424_20261102_1612/score_20261102_162020/paired_tests.tsv rows haplodose_dose-naive, haplodose-naive -->.

**TCGA bulk, arm level (6,243 tumours in 16 types; 2026-11-01 code; Holm over 6 method pairs within pooled).**
This is a negative result. The naive gene ratio was more accurate than HaploDose: median per-tumour rho 0.670 vs
0.637, with HaploDose higher in only 2,054 of 6,243 tumours (Holm p 4.8e-191)
<!-- results/tcga_bulk_20261101_1405/score_centred/paired_tests.tsv row pooled haplodose-naive -->. Naive was also
higher within 12 of 16 types <!-- docs/STATUS.md:331-332 -->. At arm level, HaploDose bulk does not beat a
gene-level median log ratio.

**TCGA bulk, sub-arm (3,013 tumours with focal loss; Holm over 10 method pairs within pooled).** Focal-loss
recall at a 5% realised FPR was 0.531 for HaploDose dose and 0.472 for naive (Holm p 1.1e-24; dose higher in
1,349, naive in 871) <!-- results/tcga_bulk_20261101_1405/subarm_20261101_1431/score_fprmatched_random_20261101_1442/paired_tests.tsv row pooled recall_loss hd_dose naive; Holm family: benchmarks/tcga_bulk/subarm.py:30 -->.
The ordering reversed at looser thresholds (0.603 vs 0.630 at 10%, 0.689 vs 0.781 at 20%; descriptive, not
tested) <!-- .../score_fprmatched_random_20261101_1442/curve_pooled.tsv patient_mean_recall_loss -->. Focal-gain
recall did not differ (0.580 vs 0.589; Holm p 1.0) <!-- .../paired_tests.tsv row pooled recall_gain hd_dose naive -->.
This development result motivated claims 3 and 6.

**U-CAN Visium HD (20 patients; secondary).** With in-section immune and stromal references, median arm rho was
0.62 <!-- docs/STATUS.md:240 -->. In-section references leaked tumour RNA: the immune and stromal bins' own profiles
tracked the patient's tumour DNA (median rho 0.43 and 0.54) <!-- docs/STATUS.md:286-288 -->. With an
epithelium-only reference, median rho was 0.75 vs 0.63 for infercnvpy (15 of 20 patients higher, nominal p 0.033,
Holm p 0.13 over 4 comparisons) <!-- docs/STATUS.md:293-297 -->. Reference choice mattered more than the method.

**Kinker at 1 Mb resolution (181 lines; profile-gap diagnosis, after claim 5; Holm over 3 per metric).** Arm rho
does not show bin-level noise, because arm medians average it out. Against 1 Mb DNA bins, the exported
`log_ratio_gene` had median Pearson 0.619, below CopyKAT 0.672. Its roughness along the genome was 0.37, against
0.29 for truth and 0.05 for CopyKAT. A running median over 5 HaploDose bins raised it to 0.696, above CopyKAT
(+0.016, Holm p 1.8e-6) and infercnvpy (+0.028, Holm p 1e-6), but below a pseudobulk smoothed the same way (0.708;
HaploDose higher in 43 of 181 lines, Holm p 5e-18) <!-- docs/STATUS.md:1324, 1436 -->. The window was chosen
after seeing these data. At 1 Mb, simple smoothed pseudobulk is the baseline to beat, and HaploDose does not beat it.

### 3.9 Development evidence: DNA copy-number change and alleles

U-CAN bulk allele features against FACETS (1,034 tumours; 2026-11-01 code): median per-arm AUC 0.93 for LOH and
0.88 for any imbalance; WGD out-of-fold AUC 0.81 for dosage features, 0.88 for allele features and 0.89 for both
(label ploidy >= 2.7) <!-- docs/STATUS.md:373-379 -->. Adding allele features to dosage raised WGD AUC by 0.06 to
0.10 (95% paired bootstrap over patients) <!-- docs/STATUS.md:379 -->. Expression alone cannot see a uniform genome
doubling; the allele features carry most of the WGD signal here. These numbers predate the leave-one-out normals
fix and were not rerun on frozen code. No single-cell allele benchmark with DNA truth exists yet.

### 3.10 Development evidence: clone identity and false splits

**Kinker implants (166 lines; Holm over 10 tool pairs).** Mean `ari_all`: HaploDose 0.167, infercnvpy 0.165, R
inferCNV 0.163, CopyKAT 0.089, SCEVAN 0.027; HaploDose vs infercnvpy p 0.91 and vs R inferCNV p 0.51, both ties
<!-- docs/STATUS.md:528 -->. An oracle mixture told the implanted region reached ARI 0.58 to 0.84 at 150 genes and
0.80 to 0.95 at 500 genes, against 0.00 to 0.34 for HaploDose at 50 to 150 genes <!-- docs/STATUS.md:470-472 -->.
Locating the region, not the acceptance test, limits small events.

**Seed and node stability of the frozen 14:24 code (Kinker 181 lines; Andor 9 lines).** Clone counts were not
stable: only 47.5% of lines had the same clone count in all 5 v6 seeds, and 21% switched between one and more
than one tumour clone by seed alone <!-- docs/STATUS.md:172-173 -->. At a fixed seed, clone counts agreed for
99.4% of lines between v4 and v6 CPUs and 93.4% between GPU and CPU <!-- docs/STATUS.md:173 -->. On Andor, MKN-45
(1 DNA clone) had 2 RNA clones in every run, NUGC-4 (1 DNA clone) had 4 to 6, and SNU-16 (2 DNA clones) had 7 to
10 <!-- docs/STATUS.md:176-177 -->. SNU-601 (15 DNA clones) gave clone-difference AUROC 0.78 to 0.84 across seeds
<!-- docs/STATUS.md:175 -->.

**v-next gates (development).** Two split-null designs were tested against five criteria fixed in advance.
Neither passed all. The frozen `cov` arm had seed agreement 0.519 (needed 0.80), MKN-45 with 2 clones in 5 of 5
seeds and NUGC-4 with 1, 4, 5, 1, 5 clones, and implant ARI 0.150
<!-- results/vnext_gate2_20261103_0612/analysis/summary.json cov.a_same_clones_frac, b_MKN-45_rna_clones, b_NUGC-4_rna_clones, d_mean -->.
Also gating the initial groups raised seed agreement to 0.917 but cut implant ARI to 0.074
<!-- docs/STATUS.md:919-922 -->. Stability came from calling fewer clones, not from more accurate ones.

**Held-out normal cells (frozen 14:24 code).** Normal epithelial cells from SMC (536 cells, 10 patients) and
Moorman (750 cells, 20 patients) gave 0 malignant calls and 0 false splits <!-- docs/STATUS.md:204-207 -->. This is
an easy null: the query cells are the same cell type as the reference.

**Visium HD subclones (negative).** 109 of 166 subclonal DNA segments overlapped a tumour clone split, against
95.5 expected by chance, and 3 patients had 21 to 36 clones <!-- docs/STATUS.md:245-247 -->. Spatial subclone
detection was not better than chance.

## 4. Discussion

HaploDose estimates relative RNA dosage about as well as the best published tool in a held-out test, and better
than three others. It does not beat pseudobulk against the same reference, and in TCGA bulk at arm level a naive
gene ratio beat it. Where segmentation helped was detection of focal losses in bulk RNA at a strict
false-positive rate: development data from TCGA suggested it, and two independent cohorts, CCLE lines and
CPTAC-3 tumours, met the prespecified rule. Even there the gain over a naive ratio was small, the interval for the
mean crossed zero, and it disappeared in stranded CCLE libraries and in some CPTAC cohorts. A fair reading is
that the HMM gives cleaner calls for strong losses at a strict threshold, and that the advantage depends on how
well the normal panel matches the samples.

The clone estimand is the weakest. On held-out implants HaploDose was not better than any comparator, and on
development lines the frozen clone step over-split DNA-clonal lines and changed clone counts between seeds.
Adding a split null made counts stable only by calling fewer clones. These results agree with the 2026-11-03 code
review, which concluded that the evidence supports a focused RNA-dosage methods paper more than a claim of
calibrated cell-level clone identification (`docs/scientific_code_review_2026-11-03.md:5`) [Inference]. The
SPECTRUM claim will test arm-level clone differences, not cell membership, so even a positive result there leaves
per-cell clone accuracy unvalidated.

Three design choices are features of the model, not demonstrated benefits. Expression programmes, depth slopes
and cross-fitting each have a stated rationale and passed targeted development checks (`docs/STATUS.md:699-710`),
but none has yet been ablated on frozen code against all three estimands with pseudobulk and the naive ratio as
floors. That ablation is the main missing analysis for a methods paper.

Of the eight prespecified claims, two hold (3 and 6, which are the same endpoint in two cohorts) and four fail
(1, 2, 4 and 5). Claim 7 is not evaluable because thyroid tumours are nearly copy-number neutral, and claim 8 was
withdrawn before any data were read <!-- docs/STATUS.md:1368, 1446 -->. No result here supports a claim of
general superiority. The goal set at the start, wins on at least three prespecified benchmarks, was not met.

What the evidence supports is narrower. In bulk RNA with a matched normal panel, HaploDose's segmented dose
detects focal losses at a strict false-positive rate slightly better than a naive ratio and clearly better than
infercnvpy, in two independent cohorts. For single-cell dosage profiles, a smoothed pseudobulk against the same
reference is at least as accurate, at arm level and at 1 Mb. For clone identity, HaploDose is at parity with
existing tools on clone differences and unstable on clone counts. The paper's contribution is the separation of
estimands, the unified model, and a benchmark design that reports failures. A gain in accuracy across the
field is not part of it.

## 5. Limitations

1. **Cancer-line references.** In Kinker and MIX-seq every tool's reference was other cancer lines, and the
   truth is relative to those lines <!-- docs/STATUS.md:526-527 -->. None of the tools was built for this, and
   results may not transfer to tumours with diploid normal cells. No scRNA benchmark here uses a diploid normal
   reference.
2. **Node and seed dependence.** Clone counts change with the random seed (47.5% of lines stable across 5 seeds)
   and, less often, with CPU generation or GPU <!-- docs/STATUS.md:172-173 -->. A spatial unit test passes on v6
   and fails on v4 CPUs <!-- docs/STATUS.md:481-485 -->. Every JAX entry point casts to float32, so
   `JAX_ENABLE_X64=1` cannot test whether float precision causes this <!-- docs/STATUS.md:194-198 -->. Dosage
   outputs are stable (Section 3.8); clone outputs are not.
3. **Pseudobulk is the ceiling for single-population profiles**, and HaploDose was below it on both Kinker and
   MIX-seq (Section 3.2).
4. **Development vs held-out.** Kinker, Andor, U-CAN, TCGA, Moorman, SMC and YCC were used to choose changes.
   Only MIX-seq, CCLE, CPTAC-3, SPECTRUM and the Schmid primary datasets support claims. Several development
   results come from code versions before the 14:24 freeze (TCGA, U-CAN alleles, Visium HD, Kinker field).
5. **Shared references.** MIX-seq lines share 6 reference sets; CCLE lines and CPTAC-3 cases share normal panels
   within groups, and the 5% FPR threshold is set on the same group's DNA-neutral windows. Bootstrap intervals
   condition on the fitted reference and threshold (F10, `docs/scientific_code_review_2026-11-03.md:108-129`).
6. **Test and interval disagree.** The bulk rules combine a signed-rank test with a positive observed mean. Both
   bulk claims meet that rule while the interval for the mean against the naive ratio spans 0 (F10).
7. **NaN cell assignments after dropping a small clone (F1).** Removing a clone below `min_cells` can leave cells
   with zero probability in every surviving clone, which produced NaN responsibilities in a synthetic probe
   (`docs/scientific_code_review_2026-11-03.md:13-21`). How often this happens in real runs is unknown
   [Unverified]. The frozen code cannot change. In the 42 SPECTRUM HaploDose runs no assigned cell had a non-finite posterior
   (Section 3.6).
8. **The split null is not a calibrated false-split rate (F2).** It replays only the proposal-and-score search
   for one clone, does not gate initial groups in the claim 4 configuration, has no specimen-level correction,
   and reuses parameters fitted on the selected cells (`docs/scientific_code_review_2026-11-03.md:23-40`).
9. **Allele scores vs an expression-only null (F3).** With alleles, observed split scores include allele
   evidence while null scores do not. This does not affect the expression-only SPECTRUM primary
   (`docs/scientific_code_review_2026-11-03.md:42-50`).
10. **Allele-balanced cells as normal reference (F4).** A tumour clone that lacks another clone's private LOH can
    pass the balance filter while carrying shared truncal changes, and allele balance cannot tell 1:1 from 2:2
    (`docs/scientific_code_review_2026-11-03.md:52-60`). The CLI has no switch to turn this off.
11. **Missing predictions can raise the clone-difference AUROC (F5).** A missing arm contrast is scored as 0
    difference; a synthetic example moved AUROC from 0.0 to 1.0. The SPECTRUM primary keeps this rule; a
    common-arm sensitivity was added before any SPECTRUM output (`docs/STATUS.md:1011-1021`).
12. **Comparator family (F6).** The original SPECTRUM test could declare success against a subset of comparators.
    Amended before any run: the family is fixed at 4, and each comparator must score >= 50% of patients
    <!-- docs/STATUS.md:999-1005 -->. A failing comparator still scores 0.5, which mixes reliability with accuracy.
13. **Input contracts (F7, F8).** The reader accepts fractional counts, and thinning can then produce negative
    counts; a specimen name reused across patients or datasets is silently pooled
    (`docs/scientific_code_review_2026-11-03.md:82-96`). No benchmark cohort here was found to contain either
    case.
14. **Bulk window truth (F9).** A window's DNA truth needs only one gene with a DNA value, so windows with mostly
    missing DNA are scored. The 10 Mb endpoint measures window discrimination, not breakpoints or very small
    events such as tandem duplications (`docs/scientific_code_review_2026-11-03.md:98-106`).
15. **Uncertainty labels (F11).** `arm_calls_calibrated` reflects panel size, not observed calibration;
    single-cell segment standard errors use a scalar $\tau$ and omit per-bin reference variance; the cell-level
    allele term omits overdispersion (`docs/scientific_code_review_2026-11-03.md:131-139`).
16. **Assembly and modality (F12).** Centromeres default to hg38 without an assembly declaration, bins can span
    large gaps, and HMM transitions depend on bin count, not distance. UMI scRNA, full-length reads, mixed
    spatial spots and bulk RNA have different noise and need separate validation
    (`docs/scientific_code_review_2026-11-03.md:141-147`).
17. **Small and rare events.** Events under about 60 genes cannot pass the 3-test-bin minimum; balanced 2:2
    gains escape the allele veto; the default 30-cell floor misses rare clones <!-- docs/STATUS.md:778-782 -->.
    In MIX-seq E2, HaploDose found almost nothing below 30 implanted cells (Section 3.2).
18. **Implanted subclones exist only in counts.** Thinning changes counts, not reads, so Numbat and allele mode
    cannot be scored on implants <!-- docs/STATUS.md:552-553 -->; Numbat is absent from MIX-seq
    <!-- docs/heldout_protocol.md:79-82 -->.
19. **Cell-level clone truth is thin.** Andor gives one line with full-depth scDNA from separate cells; DNA clones
    in MKN-45 cannot rule out a small RNA subgroup (771 cells, arm sd 0.30) <!-- docs/STATUS.md:926-927 -->.
20. **SPECTRUM query cells.** How the authors labelled cancer cells was not checked [Unverified]; if CNV signals
    were used, the query set favours CNV-rich cells for every tool <!-- docs/STATUS.md:974-975 -->.
21. **1 Mb profiles.** The exported per-bin profile is noisy (Section 3.8). Smoothing helps, but a smoothed
    pseudobulk remains more accurate on development data. The 5-bin window was never tested on held-out data.
22. **Claim multiplicity.** Claims 6, 7 and 8 were added after earlier claims had been run, and claim 8 was
    withdrawn on development evidence. Counting wins across a growing claim set overstates evidence, so each
    claim should be read on its own.
23. **Power in near-diploid cohorts.** Claim 7 shows that the bulk focal-loss endpoint needs cohorts with frequent
    focal losses. A win there would not have shown that HaploDose works in near-diploid tumours.

## 6. Data and code availability

- Package: `python/haplodose/`; Rust aggregation `crates/core/src/`.
- Frozen code, claims 1, 2, 3, 6: `results/heldout_20261102_1424/code/` (`code_sha256.txt`).
- Frozen code, claim 4: `results/vnext_freeze_20261103_0744/` (`code_sha256.txt`).
- MIX-seq held-out: input `results/heldout_mixseq_input_20261102_1159/`; analysis
  `results/heldout_20261102_1424/analysis/summary.json` and `analysis/sensitivity/summary.json`; protocol
  `docs/heldout_protocol.md`.
- Scaling: `results/scaling_combined_20261103_0102/` (`summary.tsv`, `runs.tsv`, `not_finished.tsv`).
- CCLE: `results/ccle_bulk_20261102_1551/score/primary.json`.
- CPTAC-3: `results/cptac_bulk_20261103_0618/score/primary.json`; survey `dev/cptac_survey_20261103_0612/survey.md`.
- REBC-THYR (claim 7): `results/rebc_bulk_20261106_0711/score/primary.json`; wrappers `benchmarks/rebc_bulk/`;
  first run with the sample-ID match failure `results/rebc_bulk_20261106_0649/` (unscored).
- Tahoe-100M (claim 8, withdrawn): frozen scorer `benchmarks/tahoe/score_profile.py`, Kinker check
  `dev/claim8_scorer_fixture_20261106_1058/`; metadata only in `data/external/tahoe100m_20261106/metadata/`.
- Profile-gap diagnosis: `dev/profile_gap_20261106_0007/out/`.
- SPECTRUM scorer: `results/spectrum_scorer_20261103_0613/`, amendments `results/spectrum_scorer_20261103_amend_1328/`
  and `results/spectrum_scorer_20261103_amend2_1338/`; run plan `dev/spectrum_run_plan_20261103_0750/plan.md`.
- Schmid: `benchmarks/schmid/PLAN.md`; pipeline `external/benchmark_scrnaseq_cnv_callers` at commit 319ea34.
- Development: `results/kinker_field_20261102_0529/`, `results/kinker_seeds_frozen1424_20261102_1545/`,
  `results/andor_frozen1424_20261102_1542/`, `results/ucan_bulk_frozen1424_20261102_1612/`,
  `results/tcga_bulk_20261101_1405/`, `results/vnext_gate2_20261103_0612/`,
  `results/normals_null_frozen1424_20261102_1554/`.
- Public data: MIX-seq Figshare 10298696 v3 <!-- docs/heldout_protocol.md:25 -->; Kinker GSE157220
  <!-- docs/STATUS.md:394 -->; SPECTRUM GSE180661 / SRP329496 and Synapse syn66366960 <!-- docs/STATUS.md:90-91 -->;
  Schmid truth Zenodo 20260649 <!-- docs/STATUS.md:123 -->; DepMap 24Q4; GTEx v10; GDC release 46.0
  <!-- dev/cptac_survey_20261103_0612/survey.md:3 -->. U-CAN, SMC, YCC and Moorman access terms: [CITATION NEEDED].
- Environments: `pip_freeze*.txt` and `sessionInfo_cnv_field.txt` in each run directory.
- Benchmark scripts still contain absolute workspace paths (`docs/scientific_code_review_2026-11-03.md:198`);
  a portable release is needed before external reproduction.

## 7. Open items before submission

1. Done 2026-11-04 (claim 4 fails, Section 3.6). Was: run claim 4 (SPECTRUM) once and fill every "RESULT PENDING (blind until run)" in 3.6, the abstract and the
   evidence matrix, including tool completion, both sensitivities and the non-finite posterior count.
2. Done 2026-11-06 (claim 5 fails, Section 3.7). Was: run claim 5 (Schmid): comparator reproduction gate, primary mean rank, normalised score, sensitivities.
3. Done 2026-11-03 (v-next primary, 14:24 sensitivity) <!-- docs/STATUS.md:1034 -->. Was: fix which frozen code claim 5 uses before it runs. `benchmarks/schmid/PLAN.md:14-16, 362` says the v-next
   freeze; `results/vnext_freeze_20261103_0744/README.md:6` says that freeze is for claim 4 only.
4. Ablations on frozen code, one component off at a time (programmes, depth slopes, depth matching, gene veto,
   cross-fitting, alleles, Potts, split null), scored on all three estimands with pseudobulk and naive floors.
5. A false-split analysis of the whole frozen clone workflow on independently justified one-clone specimens,
   with binomial intervals (F2).
6. A diploid-normal-reference scRNA benchmark for relative dosage.
7. Rerun U-CAN allele and WGD scoring on frozen code, or label it as pre-freeze in every figure.
8. Fix or disclose F1 (NaN after clone drop), F7 and F8 (input contracts) in a post-freeze release; keep frozen
   results as run.
9. Citations: inferCNV, infercnvpy, CopyKAT, SCEVAN, Numbat, Chang 2026, Song 2025, Kinker 2020, Andor 2020,
   SPECTRUM, DepMap, GTEx, CPTAC-3, TCGA, ABSOLUTE, ASCAT, FACETS, Neufeld et al. count splitting. Confirm the
   [Unverified] comparator-table cells against each tool's paper or source.
10. Figures (outline section 8) are not yet made for the held-out results.
11. Reconcile the source inconsistencies listed in the report that accompanies this draft (heldout_protocol
    common-universe Holm p; sign-flip vs pool Wilcoxon; STATUS U-CAN label; outline stale U-CAN and scaling
    rows) before numbers are copied into figures.
12. Done 2026-11-06 (user: "go"). The framing is a failure-inclusive methods report: estimand separation, a unified
    model, and prespecified benchmarks. The supported result is bulk focal-loss detection. There is no SOTA claim.
