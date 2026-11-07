# HaploDose manuscript outline

Draft 2026-11-02. Structure follows the assumptions-first model of Wang et al., Cell Systems 2026
(doi 10.1016/j.cels.2026.101762): state assumptions, put all methods in one notation, run targeted
analyses, list limitations. Every number below is copied from `docs/STATUS.md` (STATUS) or
`docs/HANDOFF_2026-11-02.md` (HANDOFF) with its section named. Where the two differ, the STATUS section
"Takeover review, 2026-11-02" wins. All current numbers are development evidence unless the row says
held-out. Nothing here is a claim of readiness or a journal choice.

## 1. Working title and abstract skeleton

Working title: *HaploDose: separating relative RNA dosage, DNA copy-number change and clone identity
in single-cell, spatial and bulk RNA-seq*

Abstract skeleton (4 sentences):

1. Tools that infer copy number from RNA report one profile per cell or clone, but that profile mixes
   three different quantities (relative RNA dosage, DNA copy-number change and clone membership), and
   each needs its own truth at its own resolution.
2. HaploDose segments pooled clone profiles with a Poisson-lognormal hidden Markov model against a
   fitted reference mixture, assigns single cells to those profiles with a negative binomial likelihood
   that carries depth, expression-programme and allele terms, and tests clone differences with
   cross-fitting, so each output maps to one estimand.
3. On a held-out MIX-seq set of 127 cell lines, frozen before evaluation, HaploDose arm profiles
   [E1 result pending] relative to infercnvpy, R inferCNV, CopyKAT and SCEVAN; peak memory at scale was
   [claim 2 result pending]; in held-out CCLE bulk RNA, focal-loss recall was [claim 3 result pending];
   and subclone recovery was [claim 4 result pending].
4. Simple baselines match or beat HaploDose for some tasks (a naive gene ratio at arm level in TCGA
   bulk; pseudobulk for single-population profiles), and the false-split rate on data with no true
   subclone was [false-split result pending]; we report both as limits of what RNA can establish.

## 2. Three estimands

The paper keeps these separate in every table and figure. Agreement at one level is not evidence for
another: aggregate DNA agreement does not validate cell-level clone labels.

| Estimand | HaploDose output | Truth that validates it | Resolution scored | Datasets with that truth |
|---|---|---|---|---|
| Relative RNA dosage | `profiles.tsv` `log_ratio_gene` (median per-gene log ratio in each bin) and `log_ratio`; `arms.tsv` `median_log_ratio_gene`; bulk log ratio | DNA copy number relative to the reference (line CN minus reference-line CN; FACETS or ABSOLUTE log2(TCN/ploidy)), scored by rank correlation because RNA scales DNA by the dosage slope beta | Arm (Spearman per line or tumour); sub-arm windows of about 10 Mb (bulk) | Kinker, MIX-seq, U-CAN bulk, TCGA, Moorman, U-CAN Visium HD, CCLE |
| DNA copy-number change | `segments.tsv` states loss / neutral / gain / amplification (relative, diploid scale); bulk arm calls (z vs leave-one-out normals, abs(z) > 3); allele columns and `joint_state` (allelic imbalance; `neutral_allele_imbalanced` = cnLOH candidate); WGD from allele features (benchmark script `alleles/score_alleles.py`, not a package output) | Discrete DNA calls: FACETS TCN, LCN, ploidy (U-CAN); ABSOLUTE (TCGA, DepMap 24Q4 for CCLE); MSK-IMPACT segments (Moorman) | Arm calls (AUC, altered-right, neutral-called); focal loss at a matched FPR; per-arm LOH and imbalance AUC; per-tumour WGD | U-CAN bulk (incl. alleles, WGD), TCGA (arm, sub-arm), CCLE (held-out), U-CAN Visium HD, Moorman |
| Tumour clone identity | `cells.tsv.gz` labels, `posteriors.npy`, `clones.tsv`, `differences.tsv`, `tree.tsv` | Per-cell DNA clone (scDNA), or known implanted cell sets (synthetic), or one population (single-line specimens) | Cell (ARI, difference-arm precision/recall, AUROC); clone count per DNA-clonal specimen | Andor SNU-601 (full-depth scDNA, n = 1); Kinker and MIX-seq implant designs (synthetic counts); second held-out set (planned) |

Notes for the text:
- `bulk.py` docstring: bulk output is relative RNA dosage, not DNA copy number; a uniform genome
  doubling is invisible to expression. WGD is therefore claimed only through allele features.
- Segment states are relative copy number 1-4 on a diploid scale (`profile.py`), not absolute CN.
- RNA clone groups are not validated DNA clones (STATUS, "Rules this work follows").

## 3. Unified notation

Taken from the docstrings of `profile.py`, `clones.py`, `programmes.py`, `alleles.py`, `bulk.py`,
`pipeline.py` and `kernels/{likelihood,hmm,estep}.py`. Defaults from `ProfileParams` and `CloneParams`.

Indices: cell $c$, gene $g$, bin $j = 1..J$ (10 expressed genes per bin, one chain per chromosome),
reference type $t$, clone $k$ ($k = 0$ is copy-neutral relative to the baseline, normal-like),
programme $f = 1..F$ ($F = 10$), SNP $s$.

**Reference baseline.** Each reference type has bin rates $\mu_{jt}$ with $\sum_j \mu_{jt} = 1$. The
specimen baseline is a mixture $\mu_j = \sum_t w_t \mu_{jt}$, with $w$ the Poisson MLE on the pooled
specimen counts (`mixture_weights`). Per-bin reference heterogeneity gives a mismatch variance $v_j$.
In bulk mode $\mu_j$ is the pooled normal panel and $v_j$ the robust variance of normals' log ratios.

**Pooled profile (segmentation).** For a pooled group (clone, specimen or bulk sample) with bin counts
$Y_j$ and expected counts $E_j$:

$$Y_j \sim \text{Poisson}(\lambda_j),\quad \log \lambda_j = \log E_j + \text{centre} + \beta \log(n_j/2) + o_j,\quad o_j \sim N(0, \tau^2 + v_j)$$

- $n_j \in \{1, 2, 3, 4\}$ is a hidden Markov state (loss, neutral, gain, amplification).
- $\beta$ is the dosage slope (dosage compensation times purity); default start 0.6, bounds [0.3, 1.0]
  for clones and [0.05, 1.0] for bulk, where the start is chosen by HMM likelihood on a grid.
- $\tau$ is the bin-level extra-Poisson sd, estimated from the MAD of first differences within
  chromosomes, floor 0.02.
- The integral over $o_j$ uses adaptive Gauss-Hermite quadrature (`pln_logpmf`).
- Transitions: $T_{ii} = 1 - p$, $T_{ij} = p\, q_j / (1 - q_i)$, $p = 2 \times 10^{-3}$,
  $q = (0.2, 0.5, 0.2, 0.1)$. Long-run neutral occupancy is 0.379, not 0.5 (STATUS, "Decisions").
- A segment is reported if its posterior support is at least 0.95 over at least 3 bins.

**Depth.** Centred log depth $x_c$ and a per-bin depth slope $\gamma_j$ estimated within clones
(or from reference cells with `--depth-source reference`). The pooled expectation uses the UMI-weighted
depth mean $m$ and variance $v$ of the group:
$E_j \propto \mu_j \exp\{(\gamma_j m + \gamma_j^2 v/2)/(1 + c v)\}$, $c = \sum_j \mu_j \gamma_j^2$.

**Expression programmes (nuisance).** Loadings $W \in \mathbb{R}^{F \times J}$ are the top SVD
components of fit-reference-cell residuals, high-passed by removing each arm's mean, so no programme can
explain a whole-arm change. Cells get scores $u_{cf}$ shrunk towards 0; pooled profiles get one free
amplitude per programme.

**Cell likelihood.** Clone $k$ has log rate $\ell_{kj}$ (baseline plus clone dose). For cell $c$ with
binned total $L_c$:

$$\rho_{ckj} = \operatorname{softmax}_j\Big(\ell_{kj} + x_c \gamma_j + \sum_f u_{cf} W_{fj}\Big),\quad Y_{cj} \mid z_c = k \sim \text{NB}(L_c \rho_{ckj},\ \phi_j)$$

with variance $\mu + \phi_j \mu^2$. Only compositional differences between clones matter.

**Alleles (optional).** Oriented counts $a_{cs}, d_{cs}$ from population-phased SNPs. Clone level:
HMM over balanced ($b = 0.5$) or imbalanced ($b = 0.5 \pm \theta_l$, $\theta \in \{0.08, 0.17, 0.3,
0.45\}$) states; orientation flips with probability $(1 - e^{-2\nu d})/2$, $d$ in cM; beta-binomial
emissions with concentration $\kappa$ from same-patient reference cells. Cell level adds
$\sum_s a_{cs} \log b_{ks} + (d_{cs} - a_{cs}) \log(1 - b_{ks})$.

**Clone mixture and Potts prior.** $\pi_k$ are clone weights; $LL_{ck}$ the summed cell
log-likelihood. Non-spatial E-step: $R_{ck} \propto \pi_k \exp(LL_{ck})$. Spatial E-step (mean field
over a row-normalised 8-nearest-neighbour graph $G$):

$$R_{ck} \propto \pi_k \exp\big(LL_{ck} + \lambda (G R)_{ck}\big)$$

$\lambda$ maximises how well neighbours' data-only posteriors predict each cell's own data, so spatially
random clones give $\lambda = 0$. $G$ is directed, so this is a smoothing heuristic, not the gradient of
a symmetric Potts energy (`neighbour_graph` docstring).

**Clone search and tests.** Proposals (Leiden groups, region mixtures, allele runs) use even bins and
even SNPs; difference tests use odd bins and odd SNPs. A difference between pooled groups is a 3-state
HMM on the log-rate difference (levels $-1, 0, +1$ copy step, prior 0.15/0.7/0.15) after depth-matched
pooling (coarsened exact matching on log depth), then a gene-consistency sign test and an allele veto.
Tumour clones share one pooled tumour dose except on resolved pairwise differences
(`shared_profiles`).

**Gene-level output.** `log_ratio_gene` = median over the bin's genes of
$\log((y_g + 0.5)/(E_g + 0.5))$, where $E_g$ splits the bin expectation by gene shares from the fitted
reference mixture; genes with $E_g < $ `min_expected` are left out.

### Comparator table

Sources: STATUS sections named, and `docs/related_work.md` (checked 2026-10-31 against papers or
source). Cells marked [Unverified] were not confirmed in the repo.

| Method | Null / reference model | Error model | Normalisation | Segmentation | Clone step (as benchmarked) | Needs reads / alleles |
|---|---|---|---|---|---|---|
| HaploDose | Fitted mixture of reference types; clone 0 copy-neutral | PLN for pooled bins; NB per cell; beta-binomial alleles | Library composition (softmax rates), depth slope, programme scores | 4-state HMM on pooled clone profiles | Proposals, EM, cross-fit split and merge tests; Potts prior in spatial mode | Optional alleles (phased SNP table) |
| infercnvpy | Reference cell mean [Unverified] | No probabilistic error model documented [Unverified] | Log fold-change, clipped (related_work) | None; 100-gene running window (related_work) | Tool output plus our malignant/clone rule (STATUS, head-to-head) | No |
| R inferCNV | Reference cell mean [Unverified] | Moving-average residuals; i6 HMM, Bayesian mixture filter (related_work) | [Unverified] | i6 HMM per subcluster (related_work) | Leiden subclusters (related_work) plus our rule (STATUS) | No |
| CopyKAT | Diploid cells inferred or supplied [Unverified] | Poisson-gamma MCMC per window (related_work) | Freeman-Tukey and DLM smoothing (related_work) | MCMC windows merged [Unverified] | No subclone output; our Leiden post-processing `cnv_post` (STATUS, Kinker field) | No |
| SCEVAN | Normal cells inferred or supplied [Unverified] | [Unverified] | Total-variation smoothing (related_work) | Mumford-Shah (related_work) | Native clones `--clones native`, per pool (STATUS, Kinker field) | No |
| Numbat | Reference mixture by constrained least squares (related_work) | PLN expression; binomial alleles (related_work) | [Unverified] | Joint expression-allele HMM on pseudobulks (related_work) | Iterative phylogeny (related_work) | Yes: per-line pileup and phasing (STATUS, Kinker field) |
| Pseudobulk / naive gene ratio | Reference pool (scRNA); pooled normals (bulk) | None | log2 CPM ratio (bulk) | None; arm median of gene ratios (bulk) | None (single profile) | No |

In every Kinker and MIX-seq run all tools received other cancer lines as the normal reference; none
was built for that (STATUS, Kinker field 2026-11-02).

## 4. Evidence matrix

Deviation from the requested layout: each dataset has one role (development or held-out), so the role
is a column and the three estimands are columns, rather than six estimand x role columns with half
empty. Status: done (number, code version, STATUS section), running, queued or planned. "Queued" and
"running" are from `squeue` on 2026-11-02 afternoon.

| Dataset (unit, n) | Role | Relative RNA dosage | DNA copy-number change | Clone identity |
|---|---|---|---|---|
| Kinker 2020 (line, 181 single-line; 166 implant) | Dev | Done: median arm rho, `log_ratio_gene` 0.861 vs pseudobulk 0.863, R inferCNV 0.835, infercnvpy 0.813, CopyKAT 0.765, SCEVAN 0.710; Holm over 5; run 2026-11-02 05:29, repo code clones variant C (Takeover; Kinker field). Numbat: pool 6 only (19 lines, 5 failed), all pools running | Not scored (truth is relative to reference lines) | Done: implant ari_all 0.167 vs infercnvpy 0.165 (p 0.91), R inferCNV 0.163 (p 0.51), CopyKAT 0.089, SCEVAN 0.027; synthetic counts (Kinker field). Seed/node stability running |
| Andor 2020 (line; SNU-601 full depth n = 1; 8 lines low depth) | Dev | Done: SNU-601 profile rho 0.765 (bins), 0.815 (gene column) (HANDOFF section 4) | Not scored | Done: SNU-601 precision 0.77, recall 0.89 (17/19), difference AUROC 0.76-0.87 by node and column; 00:42 code, not the frozen version (Takeover; Andor rerun). Over-splitting of DNA-clonal lines (MKN-45, NUGC-4). Seed/node refits running |
| U-CAN bulk (patient, 1,034 tumours, 120 normals) | Dev | Done: median per-tumour rho, segmented dose 0.649, log ratio 0.627, naive 0.630, infercnvpy 0.566; dose > naive 587/1,034, Holm p 1.1e-9 (Holm over 6 method pairs); frozen 14:24 code (`results/ucan_bulk_frozen1424_20261102_1612/score_20261102_162020/paired_tests.tsv`) | Done: arm AUC gain/loss, HD log ratio 0.887/0.904, naive 0.887/0.907; alleles LOH AUC 0.93, imbalance 0.88; WGD out-of-fold AUC dosage 0.81, alleles 0.88, both 0.89 (Bulk RNA mode, allele work) | n/a |
| TCGA bulk (tumour, 6,243 in 16 types) | Dev | Done, negative: median per-tumour rho HD log ratio 0.637 vs naive 0.670; naive higher pooled (Holm p < 1e-4) and in 12/16 types (Takeover; Bulk RNA mode) | Done: arm AUC gain/loss HD 0.909/0.890 vs naive 0.926/0.915; focal-loss recall at 5% FPR dose 0.531 vs naive 0.472 (Holm < 1e-4, n 3,013); gains 0.580 vs 0.589, no difference (Bulk RNA mode, sub-arm) | n/a |
| Moorman / HTAN (specimen, 9 pairs, 8 patients) | Dev | Done: median rho 0.76 vs MSK-IMPACT; all tools 0.76-0.79 (Validation against DNA; head-to-head) | Done: 75/99 altered arms right; 77/252 neutral arms called altered (Validation against DNA) | No DNA clone truth. Held-out normals called 0.007 vs 0.002 expected (debug code, jobs 22593864-74). Normals null on frozen code queued |
| SMC (specimen, 22) | Dev | No DNA truth | No DNA truth | Multi-clone 8/22; held-out normals 0.027 vs 0.020 expected (debug code). Normals null on frozen code queued |
| YCC (sample, 5) | Dev, comparator only | No DNA truth | No DNA truth | Agreement with Numbat only: malignant 0.513-0.985; clone ARI vs Numbat 0-0.04 (Latest debug results) |
| U-CAN Visium HD (patient, 20) | Dev, secondary | Done: median rho 0.62 (11-01 09:11 code); gene column 0.649 (11-02 08:44 repo, Holm 0.15); epithelium-only reference 0.75 vs infercnvpy 0.63 (p 0.033 nominal, Holm 0.13) (Validation; head-to-head; Review response) | Done: altered arms right 147/371, neutral called 34/423 (09:11); 158/371, 31/423 (08:44) | Done, negative: 109/166 subclonal DNA segments overlap a split vs 95.5 by chance; 21-36 clones in 3 patients (Validation against DNA) |
| GSE149224 (line, 3) | Dev, robustness only | HCT116 rho 0.25, truth-limited; not used for ranking | n/a | Line separation median ARI 0.89; batch/dose ARI <= 0.003 |
| Scaling (Kinker cells, 5k-1M) | Dev (08:19 code) | n/a | n/a | n/a. Compute, done on 08:19 code: HD 1M 46 min / 10.3 GiB vs infercnvpy 4.5 min / 104 GiB; 100k HD 6.2 / 4.1, infercnvpy 1.3 / 11.1, R inferCNV 245 / 64 (Compute). Frozen 14:24 rerun queued; CopyKAT 50k/100k running |
| MIX-seq (line, 127; 83 never in Kinker) | Held-out | E1 pending (Haplodose jobs queued; comparator resubmission after the amendment 4 fix check) | n/a (truth relative to reference lines) | E2 implant pending (47 implanted specimens, per protocol amendment 4f) |
| CCLE bulk (line, about 1,000) | Held-out | Reported alongside claim 3 [design detail pending] | Claim 3 focal-loss recall pending; input build running | n/a |
| Second held-out set with reads | Held-out (v-next) | Planned | Planned | Claim 4 planned. Candidate Gambardella 2022 (GSE173634, "to verify" in HANDOFF); whether it has reads for Numbat is [Unverified] |

## 5. Prespecified SOTA claims

Copied verbatim from STATUS "Prespecified SOTA target claims (user decision, 2026-11-02, recorded before
any held-out output existed)". A claim holds only if its decision rule is met on data not used to
develop the method; each is reported whatever the outcome.

1. **Profile accuracy vs published tools.** "MIX-seq held-out E1 (`docs/heldout_protocol.md`):
   HaploDose `log_ratio_gene` vs infercnvpy, R inferCNV, CopyKAT, SCEVAN; superior to a tool if Holm
   p < 0.05 (family of 5, including pseudobulk) and the 95% CI is above 0. SOTA = superior to all four
   published tools. Pseudobulk is a baseline, reported alongside, not a tool."
   Status: pending (held-out jobs queued).
2. **Peak memory at scale.** "Controlled scaling benchmark rerun on the frozen version (14:24 freeze
   code; the 08:19 run used earlier code): same cells, one node generation, 16 CPUs, 3 repeats where
   feasible. SOTA = lowest median peak memory among all tools that finish, at every size from 20k cells
   up. Wall time reported alongside (infercnvpy is currently about 10x faster; no runtime claim is
   made)."
   Status: pending (frozen rerun queued). The 08:19 run is development evidence only.
3. **Focal (sub-arm) loss detection in bulk RNA.** Decision rule as fixed 2026-11-02: "Primary:
   focal-loss recall at 5% realised FPR on DNA-neutral windows, HaploDose dose vs naive and vs
   infercnvpy, paired Wilcoxon over lines with >= 1 focal-loss window, Holm over these 2. SOTA = both
   Holm < 0.05 with HaploDose higher." Normal panel: GTEx v10 tissue-matched normals; truth DepMap 24Q4
   ABSOLUTE segments; unit = line.
   Status: pending (input build running; results blind).
4. **Subclone recovery.** "After v-next (count-splitting proposals, calibrated split null), frozen, then
   evaluated once on a second unopened held-out set that has reads, so Numbat (and, if feasible, XClone
   and CalicoST) can be included. Endpoint ari_all as E2; SOTA = superior to every comparator, Holm over
   the comparisons, CI above 0."
   Status: planned. Development evidence is a tie (Kinker 0.167 vs infercnvpy 0.165).

Known risks stated in STATUS: Numbat ranks first in recent benchmarks and is absent from MIX-seq; TCGA
arm level is a loss to the naive gene ratio (not a target claim).

## 6. Targeted analyses

Each analysis answers one question about one estimand. Status as of 2026-11-02.

**6.1 False-split calibration (clone identity).** Question: how often does HaploDose split a population
that has no true subclone?
- Held-out normal cells, SMC and Moorman, frozen 14:24 code: queued
  (`results/normals_null_frozen1424_20261102_1554/`). Debug-code values: 0.027 vs 0.020 (SMC), 0.007 vs
  0.002 (Moorman), observed vs expected.
- DNA-clonal Andor lines: development evidence of over-splitting (MKN-45 1 to 3 RNA clones, NUGC-4 1 to
  4; Andor rerun). Planned: report per-line false splits on the frozen code.
- Seed and node stability (Kinker and Andor refits over several seeds on v4 and v6, plus one GPU Kinker
  run): running and queued.
  Needed because the Potts test passes on v6 and fails on v4, and `JAX_ENABLE_X64=1` cannot test this
  since every JAX entry point casts to float32 (Takeover).
- Calibrated split null for small (< 30 cell) clones: planned, v-next.

**6.2 Ablations (which component helps which estimand).** Done (development): clones variants A-H on
Kinker (shared doses, run-wise deviation, Jaccard dedup, targeted region test, windows, permutation
scan); arm-mean programme high-pass (synthetic grid); ribosomal exclusion; depth matching; cell-type
mismatch filter; low-dosage gene exclusion (not adopted); `log_ratio_gene` vs summed bins. Planned on
the frozen version: one switch off at a time (programmes, depth slopes, depth matching, gene veto,
cross-fitting, alleles, Potts), scored on all three estimands, with pseudobulk and naive gene ratio as
floor baselines.

**6.3 Diploid-normal reference arm (relative RNA dosage).** No scRNA benchmark yet uses a diploid normal
reference; all used other cancer lines (Takeover; HANDOFF section 7). Planned. The dataset is not yet
chosen. U-CAN Visium HD gives a related result: an external normal colon epithelium reference raised
median rho from 0.62 to 0.71, and in-section immune and stromal references leak tumour RNA (U0438).

**6.4 CPU vs GPU agreement.** Question: do JAX scan (CPU) and associative scan (GPU) give the same
outputs? GPU smoke test and a GPU Kinker seed run and a 1M-cell GPU scaling run are running
(`dev/gpu_smoke_20261102_1539/`). Planned endpoint: profile and label agreement against the CPU run on
the same input.

**6.5 Scaling (claim 2).** Frozen-version rerun queued (`results/scaling_frozen1424_20261102_1541/`);
CopyKAT 50k and 100k from the 08:19 run still running. Report median wall time and peak memory over 3
repeats, peak = max(sacct MaxRSS, /usr/bin/time).

**6.6 Other targeted analyses, done (development).**
- Oracle headroom: a mixture told the implanted region reaches ARI 0.58-0.84 at 150 genes and 0.80-0.95
  at 500, against about 0 to 0.34 for every tool at 50-150 genes (HANDOFF section 4). Proposal
  localisation is the bottleneck.
- Proposal floor: below 30 implanted cells, ARI 0.051 vs infercnvpy 0.093; at or above, 0.192 vs 0.180.
- Rare-clone sensitivity by implant size, direction and share (Kinker v3 design).
- Reference leak and depth slope tracking tumour DNA on Visium HD.

## 7. Limitations

1. **Cancer-line references.** In Kinker and MIX-seq, the reference is other cancer lines. None of the
   tools was built for that, and the truth is relative to those lines. Results may not transfer to
   tumours with diploid normal cells.
2. **TCGA arm level is a negative result.** A naive gene-level median log ratio beats HaploDose bulk at
   arm level (median per-tumour rho 0.670 vs 0.637; 12/16 types). HaploDose's bulk advantage is limited
   to focal losses at strict FPR (0.531 vs 0.472 at 5%), and reverses at 10-20% FPR (descriptive).
3. **Pseudobulk is the ceiling for single-population profiles.** On single-line specimens pseudobulk is
   slightly ahead (0.863 vs 0.861, Holm p 0.03).
4. **Node dependence.** Borderline clones and the Potts test change between v4 and v6 CPUs (Andor
   SNU-601 AUROC 0.76 vs 0.82; Kinker about 0.004 ARI). Stability claims wait for the seed and node
   refits.
5. **Development vs held-out.** Kinker, Andor, U-CAN (bulk and spatial), TCGA, Moorman, SMC and YCC were
   used to choose changes. Only MIX-seq, CCLE and the second held-out set support claims.
6. **Kinker lines may carry real subclones.** Splits of single-line specimens are not all errors, so
   one-clone counts are not a clean false-split rate.
7. **Implanted subclones exist only in counts.** Thinning changes counts, not reads, so Numbat and allele
   mode cannot be scored on the implant design.
8. **Cell-level clone validation is thin.** One line (Andor SNU-601) with full-depth scDNA, scored on
   00:42 code, not the frozen version. The 8 low-depth Andor lines cannot rank tools.
9. **Allele mode validated in bulk only.** LOH, imbalance and WGD results are from U-CAN bulk RNA; no
   scRNA allele benchmark with DNA truth.
10. **Small events.** Events under about 60 genes cannot pass the 3-test-bin minimum; balanced 2:2 gains
    escape the allele veto; the default 30-cell floor misses rare clones.
11. **Visium HD subclones are not better than chance** (109/166 vs 95.5), and in-section references can
    leak tumour RNA.
12. **U-CAN bulk allele, LOH and WGD numbers predate the leave-one-out normals fix** (2026-11-01 code); not
    rerun. The dosage rho was rerun on frozen 14:24 code (STATUS, "U-CAN bulk rerun on frozen 14:24 code").
13. **Bootstrap CIs ignore pool clustering.** MIX-seq has 6 reference sets, too few for a cluster
    bootstrap (protocol amendment 3a).

## 8. Figure plan

Development panels stay labelled "development"; held-out panels are new and separate.

- **Figure 1. Estimands and model.** (a) Schematic: relative RNA dosage, DNA copy-number change, clone
  identity, with the truth for each. (b) Generative model in the section 3 notation. (c) Cross-fitting
  (even bins propose, odd bins test).
- **Figure 2. Held-out profile accuracy (claim 1).** (a) MIX-seq E1 per-line rho, all tools and
  pseudobulk. (b) Paired differences with 95% bootstrap CIs and Holm p. (c) Sensitivity: 83 non-Kinker
  lines; pool-level medians.
- **Figure 3. Development profile accuracy.** (a) Kinker single-line rho by tool (existing A1/A2).
  (b) U-CAN bulk and TCGA per-tumour rho, including the TCGA arm-level loss. (c) TCGA focal-loss recall
  vs FPR.
- **Figure 4. Bulk DNA copy number and alleles.** (a) CCLE focal-loss recall at 5% FPR (claim 3).
  (b) U-CAN LOH and imbalance AUC per arm. (c) WGD AUC, dosage vs alleles vs both.
- **Figure 5. Clone identity.** (a) Implant ARI by size, direction, share, Kinker and MIX-seq E2.
  (b) Floor strata. (c) Oracle headroom (existing C). (d) Andor SNU-601 precision and recall (existing D).
- **Figure 6. False splits and stability.** (a) Held-out normals, observed vs expected. (b) Clone count
  on DNA-clonal lines. (c) Seed and node agreement (v4, v6, GPU).
- **Figure 7. Ablations.** Component off vs full model, one row per estimand.
- **Figure 8. Compute (claim 2).** Wall time and peak memory vs cells, frozen version (replaces E1/E2).
- **Supplement.** Visium HD (reference choice, leak diagnosis), Moorman, YCC agreement with Numbat,
  Andor low-depth table, GSE149224 truth-quality check.

## 9. What would falsify the main claims

The profile claim fails if, on MIX-seq E1, HaploDose does not reach Holm p < 0.05 with a CI above 0
against any one of infercnvpy, R inferCNV, CopyKAT or SCEVAN; it also weakens if the gap to pseudobulk
grows beyond the development value or the 83 non-Kinker lines disagree with the full set. The memory
claim fails if any finishing tool has a lower median peak at any size from 20k cells on the frozen
rerun. The bulk claim fails if CCLE focal-loss recall at 5% FPR is not higher than both naive and
infercnvpy at Holm < 0.05, which TCGA arm-level results make a real possibility. The claim that HaploDose
separates the three estimands fails if the held-out normals or DNA-clonal lines show false splits well
above the expected rate, or if seed and node refits change clone counts on the same data, because then
clone labels would carry run noise rather than biology. Claim 4 fails if Numbat, or any comparator, is
not beaten on the second held-out set; the development tie makes this the least likely claim to hold
[Inference].
