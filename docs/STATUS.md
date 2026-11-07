# Haplodose: status, decisions and open work

## Membership uncertainty and merging (2026-11-07)

The [membership investigation](MEMBERSHIP_2026-11-07.md) traces the three remaining
failures and tests two proposed corrections. The supported chr12 profile contrast
survives refinement in seed 20261108, but reassignment and hard-pool merging remove
the distinction. Stable MAP labels do not certify convergence; using all eight
refinement rounds nevertheless reduces recovery from 2/5 to 1/5 and worsens the
observed error/power tradeoff. That default change is rejected and preserved as a
negative result. The retained implementation records actual refinement rounds,
explicitly uncertified stopping status, and conditional pair-membership entropy
for merge decisions. It does not interpret a merge as genomic equality.

Integrating uncertain memberships improves conditional prediction on the selected
cases, but the complete binned-search event gate fails in both split modes. In the
mode-matched thinning experiment it has 10/200 state-null flags, a simultaneous
30-endpoint upper bound of 11.26%, and zero detections in 200 alternatives. This is
an RNA event statistic, not a validated clone caller; initial group contrasts and
public gene/reference preprocessing are not fully covered. No replacement rule
or biological calibration claim is promoted.

Compute job 22850847 passes 315 Python tests (one GPU-only skip) and four Rust
tests; three targeted statistic checks also pass (two overlap the full suite).
The unchanged five-case panel remains 2/5, with all three failures preserved.
Independent audits verify input identity, counts, confidence bounds and observer
invariance. A different CPU generation changes experimental proposal eligibility
in two specimens; the original-host replay and first failed hash audit are both
retained. Cross-hardware numerical certification remains unresolved.

## Recovery causes and error/power replay (2026-11-07)

The [cause investigation](RECOVERY_CAUSES_2026-11-07.md) repairs a failed optimizer
status that discarded a numerically recoverable proposal. The identical objective
is retried with SLSQP and must pass all previous numerical certificates. Compute
job 22846715 passed 311 Python tests (one GPU-only skip) and four Rust tests.
The original scientific assertions, generators and frozen results are preserved.

Complete recovery remains 2/5: all three previously failing cases still fail.
One restored proposal is later merged; two others fail confirmation because
assigned pools dilute the contrast. Post-fit truth-group probes and a gene-level
known-nuisance oracle distinguish this from complete absence of measurement
information, but do not provide deployable labels or nuisance estimates.
Neither reference-derived depth slopes nor removing programmes improves recovery.
Normal-cell errors in one case worsen from 9/90 to 10/90; this is retained.

The paired 998-fit regression replay and independent count/bound audit leave the
false-split/power result unchanged: the state-calibrated gate has 8/200 errors
(simultaneous upper bound 8.96%) and 65/200 detections (lower bound 24.34%).
The 5% target remains unmet. RNA-equivalent transcriptional-state and CNA
explanations also prevent a biological identity claim without independent evidence.
No additional scientific default or calibration claim is promoted.

## Recovery and uncertainty follow-up (2026-11-07)

The original thinning recovery regression is repaired without changing its test,
seed or scientific criteria. Final candidate 5 passed 310 Python tests (one skipped)
and four Rust tests in compute job 22844951. On five frozen simulation inputs,
the preserved baseline passed 0/5 original recovery assertions and the correction
passed 2/5, without a new failing criterion. The other failures remain evidence
of limited recovery. This is not a biological validation panel.

The [follow-up report](RECOVERY_2026-11-07.md) records the numerical corrections,
complete binned-search false-distinction experiment and the separate optional
systematic-RNA interval contract. The latter requires independently fixed
expectations, known covariance and fixed genomic eligibility. Existing approximate
SEs remain scale diagnostics; automatic inference does not acquire a calibration
claim from the new helper. Historical evidence and decisions below are preserved.
State-calibrated global testing still fails the clone criterion: 8/200 false
distinctions, simultaneous upper bound 8.96%, and 65/200 detections. The optional
known-covariance intervals covered all 1,000 simulations at each of six endpoints,
but were wide and often insensitive. Neither result establishes biological validity.

## Schema-2 correctness and scientific validation (2026-11-06)

Current development follows the [scientific contract](SCIENTIFIC_CONTRACT.md) and
[0.2 migration](MIGRATION_0.2.md). The historical entries below are preserved as
records of decisions and outcomes, not renewed claims for the revised method.
Former evaluation cohorts are regression-only. No suitable unopened biological
panel has been verified for a 5% specimen false-split guarantee.

The initial dirty source was captured before edits under
`dev/hardening_20261106_2119/baseline`, with its Git state, file hashes and selected
frozen-result hashes alongside it. The isolated baseline wheel passed Rust 3/3
and Python 99 passed/1 skipped on compute job 22811740. Independent baseline
counterexamples (job 22812061) nevertheless reproduced invalid constant-reference
rank, an incorrect responsibility-weighted programme objective, accepted
fractional/duplicate inputs and lost integer aggregation precision. Candidate
validation and scientific falsifiers are recorded in new task directories;
passing these does not establish biological clone identity or error control.

Written 2026-11-01. Covers work from 2026-10-31 to 2026-11-01. Plan file:
`~/.claude-teams/plans/t3-infercnv-like-caller-handoff-dapper-cat.md` (milestones M0-M4).
Git: one commit `d6ed692` (2026-11-01 08:35), made while 2 tests were failing.
Current plan: `~/.claude-teams/plans/inspect-the-codebase-clearly-vast-lampson.md`.

## Manuscript direction (user decision, 2026-11-02)

HaploDose is intended to become a full scientific manuscript: a methods paper
with a focused benchmark and a clear account of what RNA can reliably establish
about genomic dosage and tumour clones. The reference for presentation is
Wang et al., *Categorization and analysis of 16 computational methods for
estimating cell potency from single-cell RNA-seq data* (Cell Systems, 2026;
[DOI](https://doi.org/10.1016/j.cels.2026.101762)). That paper is a Synthesis with
targeted reanalyses; borrow its treatment of assumptions, mathematical notation,
implementation and validation limits. This direction does not establish current
readiness, acceptance prospects or a target journal.

The manuscript should distinguish relative RNA dosage, DNA copy-number
alterations and clone identity as separate estimands. Bulk DNA agreement supports
aggregate profile fidelity, not cell-level clone accuracy. Explain which
HaploDose components improve which task through ablations, including cases where
simple gene-ratio or pseudobulk baselines perform as well or better. Existing CNV
benchmarks are documented in [related work](related_work.md); the contribution
must go beyond another broad tool comparison.

Evidence requirements before manuscript claims are finalized:

- Freeze the method, parameters, scoring and primary endpoints before independent
  evaluation. Record datasets already used to select outputs or tune the method
  as development evidence, even when no model parameters were fitted on them.
- Validate at the resolution claimed, using matched DNA where available and
  documenting truth quality, sample correspondence and uncertainty. Synthetic
  implants complement biological validation; they do not replace it.
- Use patients or independently derived cell lines as biological replication
  units, account for repeated samples and shared pools, and specify covariates,
  effect sizes, uncertainty and the multiple-testing correction and family.
- Include simple baselines, fair comparator inputs and post-processing, component
  ablations, reference mismatch, expression-program confounding, rare-clone
  sensitivity and false splits. Report failures and negative results.
- Produce reproducible figures and tables from versioned code, inputs, seeds and
  environment inventories. Resolve material numerical and node-dependent
  behavior before making stability or reproducibility claims; compare resource
  use under controlled conditions.

Current results below are development evidence from multiple code versions, not
a frozen manuscript result set. Broad superiority and reliable subclone recovery
are not established. Bulk and scRNA-seq remain the primary benchmarks; spatial
RNA-seq provides a secondary evaluation. The package and manuscript scope remain
general-purpose; the CRC/IGF2 motivation is a research question, not an established
mechanism.

## Prespecified SOTA target claims (user decision, 2026-11-02, recorded before any held-out output existed)

User goal: HaploDose is to be state of the art in at least three key benchmarks. To keep that claim testable,
the claims are fixed here before results are seen; each is reported whatever the outcome, losses included.
A claim holds only if its decision rule is met on data not used to develop the method.

1. **Profile accuracy vs published tools.** MIX-seq held-out E1 (`docs/heldout_protocol.md`): HaploDose
   `log_ratio_gene` vs infercnvpy, R inferCNV, CopyKAT, SCEVAN; superior to a tool if Holm p < 0.05 (family of 5,
   including pseudobulk) and the 95% CI is above 0. SOTA = superior to all four published tools. Pseudobulk is a
   baseline, reported alongside, not a tool. Tightened by protocol amendment 5 (before any held-out job started):
   also requires the common-gene-universe E1 sensitivity to give a CI above 0 against each published tool;
   `superior` is stated as conditional on independence of lines (lines share 6 reference sets).
2. **Peak memory at scale.** Controlled scaling benchmark rerun on the frozen version (14:24 freeze code; the
   08:19 run used earlier code): same cells, one node generation, 16 CPUs, 3 repeats where feasible. SOTA = lowest
   median peak memory among all tools that finish, at every size from 20k cells up. Wall time reported alongside
   (infercnvpy is currently about 10x faster; no runtime claim is made).
3. **Focal (sub-arm) loss detection in bulk RNA.** Independent set: CCLE bulk RNA vs DepMap CN (about 1,000
   lines; not used in development). Endpoint as TCGA `subarm.py` at a matched 5% false-positive rate; HaploDose
   bulk vs naive gene ratio and infercnvpy, paired over lines, Holm over the comparisons. Open design point, to be
   fixed before any run: the normal panel for bulk mode (CCLE has no matched normals).
   **Fixed 2026-11-02 (user decision, before any CCLE data was downloaded):** normal panel = GTEx v10 tissue-
   matched normals (RNA-SeQC 2.4.2 gene reads, the same quantifier as DepMap 24Q4
   `OmicsExpressionRNASeQCGeneCountProfile.csv`); each CCLE lineage gets the GTEx tissue(s) of its organ; lineages
   with no GTEx match are excluded; the same panel for all three methods. Truth: DepMap 24Q4
   `OmicsAbsoluteCNSegmentsProfile.csv` (ABSOLUTE TCN and ploidy, the TCGA `subarm.py` truth type). Unit = line
   (default model profile). Primary: focal-loss recall at 5% realised FPR on DNA-neutral windows, HaploDose dose
   vs naive and vs infercnvpy, paired Wilcoxon over lines with >= 1 focal-loss window, Holm over these 2.
   SOTA = both Holm < 0.05 with HaploDose higher. Lines whose DepMap CN was used as truth in development (Kinker,
   MIX-seq) are kept (their bulk RNA was never used) and reported as a sensitivity subset. The lineage-to-tissue
   map and the frozen code are written to the run directory before scoring. Caveat known in advance: tissue
   normals contain stroma and immune cells that cultured lines lack.
4. **Subclone recovery.** After v-next (count-splitting proposals, calibrated split null), frozen, then evaluated
   once on a second unopened held-out set that has reads, so Numbat (and, if feasible, XClone and CalicoST) can be
   included. Endpoint ari_all as E2; SOTA = superior to every comparator, Holm over the comparisons, CI above 0.
   Development evidence today is a tie (Kinker 0.167 vs infercnvpy 0.165).

   **Data fixed 2026-11-02 (user decision, before any download):** primary = MSK SPECTRUM ovarian (tumour scRNA
   reads GSE180661/SRP329496; DLP+ scDNA clones Synapse syn66366960; clone truth is profile-matched, as Andor);
   secondary, cell-level = Bian 2018 CRC scTrio (GSE97693, same-cell DNA+RNA; DNA CN re-derived). Both stay
   unopened (download and read QC only) until v-next is frozen.
   **Endpoint amended 2026-11-02, before any download** (dataset survey, `dev/dataset_survey_20261102/survey.md`):
   `ari_all` needs a truth clone per RNA cell, which SPECTRUM cannot give (DNA and RNA from different cells;
   mapping RNA cells to DNA clones would be circular). SPECTRUM primary endpoint = per-patient clone-difference
   AUROC as in `benchmarks/andor/score.py` item 2: truth per arm = 1 if any two DLP+ tumour clones differ on
   >= 30% of the arm's bins by >= 1 copy; tool score per arm = largest |difference| between any two of its RNA
   tumour clones (clones of the cells the tool calls malignant, >= max(20, 2%) cells). Patients: tumour scRNA
   from an unsorted or CD45-negative sample and >= 2 DLP+ tumour clones with >= 1 differing arm (clone counts
   and arm flags only; fixed from the Synapse clone tables before any RNA tool runs). Test: paired Wilcoxon over
   patients, HaploDose v-next vs each comparator that runs (Numbat, infercnvpy, CopyKAT, SCEVAN; XClone and
   CalicoST if feasible), Holm over those; SOTA = every Holm p < 0.05 with bootstrap 95% CI of the mean
   difference above 0. Secondary, descriptive: precision/recall of called differences, RNA vs DNA clone counts.
   The SPECTRUM scorer is written and hashed on Andor development data before any SPECTRUM RNA output is read.
   Bian scTrio: `ari_all` against same-cell DNA clone labels (DNA clusters from the re-derived CN), descriptive.
   Bian provenance note (2026-11-02): the CRC project (~/CRC src/46-49, results/bian_*) has processed Bian DNA
   (raw DNA alignment, DNA clones, alleles, chr11p) for the IGF2 question; no RNA CNV caller or HaploDose was run on
   Bian there (grep of src/46-49 and results/bian_*). HaploDose development has not used Bian. Bian stays secondary
   and descriptive.
   Claim 3 amendment (same day, before the score stage ran): descriptive per-protocol (Stranded) sensitivity
   added (`results/ccle_bulk_20261102_1551/README.md`).
   **Claim 3 result (read once, 2026-11-02, frozen 14:24 code, `results/ccle_bulk_20261102_1551/score/primary.json`):
   holds by the prespecified rule.** n = 822 lines with >= 1 focal-loss window (unit = line; 18 lineage groups).
   Focal-loss recall at 5% FPR, mean: HaploDose 0.509, naive 0.484, infercnvpy 0.423. Paired Wilcoxon, Holm over
   2: vs naive diff +0.024, Holm p = 0.0032, bootstrap 95% CI -0.001 to 0.049 (crosses 0; descriptive); vs
   infercnvpy diff +0.085, Holm p < 1e-10, CI 0.062 to 0.108. Excluding development lines (n = 646): same
   direction, vs naive Holm p = 0.011, CI -0.005 to 0.053. Stranded sensitivity (descriptive): unstranded lines
   HaploDose ahead of both (Holm p < 0.001); stranded lines no advantage over naive (diff -0.023, Holm p = 0.69).
   Report the margin over naive as small and protocol-dependent.
5. **Outside benchmark (Schmid et al. 2025, Nat Commun 16:8777; added 2026-11-02 by user decision, before any
   run).** Frozen HaploDose run through the authors' pipeline (github.com/colomemaria/benchmark_scrnaseq_cnv_callers,
   truth Zenodo 20260649) next to their six tools (inferCNV, copyKat, SCEVAN, CONICSmat, CaSpER, Numbat), all
   rerun by us in their pipeline. Primary set = datasets never used in our development: COLO320-HSR, BCC06,
   BCC06post, MM199 (n = 4). Line-overlap set, reported separately: A375 and MCF7 (Kinker lines), HCT116
   (GSE149224). The 9 Andor gastric lines and the mouse set are excluded (development; mouse unsupported).
   Metrics = the paper's Pearson correlation and gain/loss AUC per dataset. With n = 4, no test can pass Holm,
   so the rule is descriptive: SOTA = best mean rank over the 4 datasets on both metrics, against all six tools.
   Reported whatever the outcome.
   Claim 5 amendment (2026-11-02 ~21:00, before any Schmid run or Schmid content was read; from reading the
   authors' code at commit 319ea34, `benchmarks/schmid/PLAN.md`): (a) "both metrics" is three numbers: Pearson,
   AUC gains, AUC losses; SOTA = best mean rank over 3 metrics x 4 primary datasets; the paper's normalised score
   reported descriptively. (b) Ranked against all 9 outputs the paper ranks (inferCNV, SCEVAN and Numbat each as
   CNV and Expr, plus copyKat, CONICSmat, CaSpER). (c) Primary scoring = the authors' combined scorer with
   HaploDose added (it intersects bins across methods; frozen HaploDose is autosomal, so chrX drops for all);
   sensitivity = the per-method scorer and the comparators without HaploDose. (d) Gate before HaploDose is scored:
   our comparator reruns reproduce the published comparator values on one primary dataset; if not, report the gap.
   (e) A375 and HCT116 are DNTR-seq (plate-based full-length): HaploDose there is an untested modality, reported
   as such. (f) External normal references fixed (best match to the authors' unnamed files; deviation recorded):
   COLO320 = gut cell atlas Full_obj_raw_counts_nosoupx_v2.h5ad, Diagnosis "Healthy adult", Region "LargeInt",
   category "Epithelial", atlas row order, then R set.seed(1); sample(1:n, 5000) (n = 38,612); REC and APD excluded,
   no rectum sensitivity. MCF7 = Tabula Sapiens TS_Mammary raw_counts, method "10X", compartment "epithelial"
   (5,630 cells). The selected cell IDs are saved and used by every tool. CellRanger 8.0.1 here vs 7.0.0 in the
   paper (same refdata-gex-GRCh38-2020-A), recorded as a deviation. Comparators run in a separate env pinned to
   the paper's versions (`schmid_r41`, user-approved 2026-11-02).
   (g) Gate pass rule for (d): our rerun of each comparator output is within 0.05 of the published value on Pearson,
   AUC gains and AUC losses (values transcribed from the paper's Source Data before the rerun). An output that misses
   is still ranked, flagged, and the gap reported. SCEVAN pinned to commit 92761f6 (last 1.0.1); if SCEVAN misses the
   gate, it is rerun once at the v1.0.1 tag 90b45b9 and the closer of the two is used, both reported.
   Exposure note (2026-11-02 ~22:20): while transcribing published comparator values, an agent's glob printed the
   first 4 rows (chr1, ~1 Mb) of the paper's published truth-track Source Data (supp fig 19A-C: MM, BCC06,
   BCC06post). Published data; no HaploDose run or tuning followed. Numbat 1000G hg38 panel: official pklab file
   unreachable (bot check); CRC-built chr1-22 panel used (their phasing uses chr1-22 only), recorded as a deviation.

Known risks: Numbat is ranked first in recent benchmarks (Chang 2026, Song 2025) and is absent from MIX-seq;
TCGA arm level is currently a loss to the naive gene ratio (not a target claim).

## Held-out MIX-seq result (read once 2026-11-02 ~22:30; full numbers in docs/heldout_protocol.md "Results")

Claim 1 does NOT hold. n = 127 lines. Haplodose median arm rho 0.832 beats infercnvpy, CopyKAT and SCEVAN (Holm
over 5, all p_holm < 1e-4) but ties R inferCNV (0.823, p_holm 0.15) and is below pseudobulk (0.840, p_holm 0.006).
E2 implant clone recovery (n = 47): no superiority to any tool; R inferCNV highest (0.223 vs 0.136).
SOTA count after this read: claim 3 holds; claim 1 fails; claims 2, 4, 5 open. Three wins now need all of
claims 2, 4 and 5, or a new prespecified independent benchmark.

## Seed, node and GPU stability of the frozen 14:24 code (2026-11-02, development data)

Kinker single-line design, 181 line specimens, `results/kinker_seeds_frozen1424_20261102_1545` (v6 seeds 20261031, 1-4;
v4 seeds 20261031, 1; GPU b6k seed 20261031):
- Arm dosage is stable: median per-line rho 0.853-0.856 across all 8 runs (pseudobulk 0.863); per-line rho range
  across the 5 v6 seeds median 0.014, max 0.053. v4 vs v6 at the same seed: |rho diff| max 0.004; GPU vs CPU max 0.015.
- Clone counts are NOT stable: only 47.5% of lines get the same clone count in all 5 v6 seeds; 21% switch between
  1 and > 1 tumour clones by seed alone. Same seed, other hardware: clone count equal for 99.4% (v4) and 93.4% (GPU).
Andor 9 gastric lines, `results/andor_frozen1424_20261102_1542` (same 7 CPU runs):
- SNU-601 (15 DNA clones): difference AUROC 0.78-0.84 across seeds, precision 0.70-0.76, recall 0.84.
- Over-splitting where DNA is clonal or near clonal: MKN-45 (1 DNA clone) 2 RNA clones in every run; NUGC-4 (1 DNA
  clone) 4-6; SNU-16 (2 DNA clones) 7-10 RNA clones with 21-29 arms called different (precision 0.03-0.04).
  HGC-27 difference AUROC ranges 0.29-0.71 by seed.
Consequence: the frozen clone step over-splits and is seed-sensitive; the clone-identity estimand (claim 4) depends
on v-next (count splitting plus a calibrated split null). Dosage claims are unaffected. Node dependence is small
at fixed seed; seed dependence is the larger effect and is reported as a limitation of the frozen version.

## Takeover review, 2026-11-02 ~15:30 (supersedes "Current state" below where they differ)

One session now owns the whole repo (HaploDose1/2 lanes in `docs/HANDOFF_2026-11-02.md` no longer apply). Plan:
`~/.claude/plans/review-docs-handoff-2026-11-02-md-and-ta-robust-pascal.md`.
- **Held-out run repaired before it started** (protocol amendment 4). The 14:24 job scripts would have failed: the
  comparator run glob matched their own stderr files, so `evaluate.py` and the `afterok` analysis would never have
  run, and any tool failure crashed the pipeline instead of scoring NaN/0. Comparator jobs 22611772-75, 22611777-80 and analysis
  22611794 were cancelled while pending (no held-out output exists; Haplodose jobs 22611771/76 kept, unaffected) and
  are resubmitted with the fixed `benchmarks/heldout/run_tool.sbatch`,
  `analyse.py` and `results/heldout_20261102_1424/analyse_v2.sbatch`. Pre-run check on Kinker pools 16/18:
  `dev/heldout_fixcheck_20261102_1528/` (job 22612248). Frozen Haplodose code is unchanged (all 41 hashes match).
- **The float64 v4/v6 experiment cannot answer its question.** Every JAX entry point casts to float32
  (`kernels/likelihood.py:41-43`, `kernels/hmm.py:139-150`, `profile.py:84-87`, `kernels/estep.py`), so
  `JAX_ENABLE_X64=1` changes nothing there. v4 with x64 still fails the Potts test (22611904/05). Reading
  [Inference]: the fixture sits on hard thresholds (posterior 0.95, split score > 0, Jaccard 0.7) and BLAS/XLA
  rounding tips it. Node dependence is to be measured (seed and subsample refits) and stated as a limitation.
- **U-CAN bulk rerun on frozen 14:24 code** (`results/ucan_bulk_frozen1424_20261102_1612`, development data;
  replaces the older table): n = 1,034 tumours (28 skipped, NA purity/ploidy), per-tumour arm Spearman rho vs
  FACETS. Median rho: HaploDose arm calls 0.627, HaploDose dose 0.649, naive 0.630, infercnvpy 0.566. Paired
  Wilcoxon, Holm over all 6 pairs: HaploDose calls vs naive p = 0.81 (no difference); HaploDose dose > naive
  p_holm = 1e-9; both HaploDose outputs > infercnvpy p_holm < 1e-15.
- **Normal-cell false-split null** (`results/normals_null_frozen1424_20261102_1554`, frozen code, development
  data): held-out normal-epithelium cells as query, pooled across patients. SMC 536 cells / 10 patients, Moorman
  750 cells / 20 patients: 0 malignant calls, 1 clone (normal), 0 false splits in both. Easy null (same cell type
  as the reference); harder nulls (DNA-clonal lines, seed refits) still pending. Summary recomputed with the
  corrected parser: `null_summary_v2.tsv` (the original `null_summary.tsv` miscounted tumour-clone cells).
- Gambardella 2022 (GSE173634) is no longer the second held-out set: rejected (lines from different
  donors, bulk truth only); claim 4 uses SPECTRUM and Bian 2018. Manuscript outline: `docs/manuscript_outline.md`.
- **Test count**: 52 test functions, 61 cases with parametrisation, plus 3 Rust tests. The "70/70" in the handoff
  and "61" here count differently; quote 61 Python cases. README is not edited (standing rule) and is stale.
- **Evidence state**: every number in this file is development evidence. Gaps that bear most on the manuscript:
  no calibrated false-split rate on data with no true subclone; cell-level clone validation only on Andor SNU-601
  (n = 1, scored on 00:42 code, not the frozen 14:24 version); no diploid-normal reference in any scRNA benchmark;
  allele mode validated against DNA only in bulk. TCGA arm level is a negative result (naive gene ratio 0.670 vs
  Haplodose 0.637 median per-tumour rho). Canonical Kinker single-line figure: gene column 0.861 (field run
  2026-11-02 05:29, Holm over 5); the 0.863 elsewhere is the 2026-11-01 decomposition run.

## Current state (read this first)

- **2026-11-02: 61 Python tests pass on v6 (jobs 22601385, 22601594; 3 Rust in 22601012); Potts fails on v4** (node-dependent,
  see the Kinker subclone section). Repo `clones.py` = ablation variant C (shared tumour doses, Jaccard event
  dedup, all-pairs split test).
- Kinker vs the field: subclone design v3 Haplodose 0.167-0.169 vs infercnvpy 0.165 ARI (tie; region found
  0.36-0.40 vs 0.67); single-line rho 0.854 (log_ratio_gene 0.863) vs pseudobulk 0.863. CopyKAT, SCEVAN, R inferCNV and Numbat runs
  in progress.
- Published CRC results (`~/CRC/results/haplodose_benchmark/{SMC_full,Moorman_full}/`)
  are from jobs 22589677 and 22589678 (2026-10-31). They predate ribosomal-gene
  exclusion, alleles, spatial mode, depth-matched tests, vetoes and the
  cell-type mismatch filter. All later real-data runs are debug runs under
  `logs/runs/` and are not published.

## Validation against DNA (2026-11-01, code after the event_mixtures NaN fix)

Scorers: `benchmarks/ucan_spatial/`, `benchmarks/moorman/`, `benchmarks/cell_lines/`, `benchmarks/ycc/`.
Comparator runners (HaploDose2): `benchmarks/comparators/`. Datasets: `docs/benchmark_datasets.md`.

**U-CAN Visium HD vs WGS FACETS** (20 patients, one section each; `results/ucan_spatial_20261101_0911/`):
- Spearman RNA vs DNA log2(TCN/ploidy) across arms: median 0.62 (0.19 to 0.93). U0396 and U0991 have no
  tumour clone (all cells neutral). U0438 (IGF2 TD) 0.20 with 0/27 altered arms called.
- Discrete arm calls: 147/371 DNA-altered arms in the right direction; 34/423 neutral arms called altered.
- Continuous profile: right sign on 287/336 altered arms. Arm-level RNA SD on DNA-neutral arms 0.085 vs
  median +0.10 on gained arms. The noise floor is per-gene baseline mismatch, not sampling.
- Subclones: 109/166 subclonal DNA segments (CF < 0.8 x purity, >= 10 Mb) overlap a tumour clone split,
  vs 95.5 by chance (random placement). Subclone detection is not better than chance. 21 to 36 clones in
  U0007, U0046, U2955 (over-splitting).
**Moorman vs MSK-IMPACT** (9 same-sample pairs, 8 patients; `results/moorman_20261101_0911/`): median rho 0.76;
  75/99 altered arms right; 77/252 neutral arms called altered (DNA threshold |log2| 0.2, not purity adjusted).
**Cell lines GSE149224**: GEO SW480 labels do not match DepMap SW480 (expression Spearman 0.07; HCT116 and RKO
  0.42, 0.43), so only HCT116 (rho 0.25, same as DepMap's own bulk RNA ceiling 0.29), batch/dose ARI (<= 0.003)
  and line separation (median ARI 0.89) are usable.
**Fix**: `event_mixtures` gave 0 x log 0 = NaN for cells with no expected UMIs in a region, silently dropping
  region split candidates (Moorman). Now `xlogy`; tests green (job 22594015). Effect on all benchmarks: none
  material (scale subclone ARI still 0.19).

## Head-to-head and reference experiments (2026-11-01)

U-CAN Visium HD vs FACETS, 20 patients (comparators by HaploDose2, `results/comparators_20261101_0915/`,
table `collected_20261101_102632.tsv`; malignant/clone calls for both inferCNVs are tool + our rule):

| Tool | median rho | no result | altered right | neutral called | subclonal splits (chance) | median time / RSS |
|---|---|---|---|---|---|---|
| Haplodose | 0.62 | 2 | 147/371 | 34/423 | 109 (95.5)/166 | 84 s / 1.9 GB |
| infercnvpy | 0.62 | 1 | 168/355 | 32/416 | 121 (111.3) | 51 s / 1.0 GB |
| T3 prototype | 0.61 | 0 | 21/331 | 3/386 | 57 (32.7), cn_state runs | 40 s |
| R inferCNV | 0.64 | 13 (16/27 jobs failed) | 66/129 | 62/187 | 16 (15.4) | 678 s / 4.5 GB |

Paired per-patient rho, Haplodose vs infercnvpy / prototype / R inferCNV: higher in 12/18, 14/18, 6/6;
Wilcoxon p 0.23, 0.021, 0.031; Holm over these 3: none below 0.05. Moorman same-sample median rho 0.76-0.79
for all tools. Numbat on YCC: 10-34 min per sample vs Haplodose with alleles 1-6 min.

Experiments on the 20 U-CAN patients (none used for fitting anything):
- **External normal colon epithelium reference** (bin2cell cells from JGCRC4NT, JGCRC14NT mucosa, P3_NAT,
  P5_NAT; 40 um tiles; `benchmarks/ucan_spatial/external_ref/`, `results/ucan_spatial_extref_20261101_1156/`):
  tumour clone found in 20/20 (was 18); median rho 0.71 vs 0.62 on 18 common patients; altered arms right
  200/364 vs 147/371; neutral called 46/416 vs 34/423; tau 0.227 vs 0.232; paired Wilcoxon p 0.23 (10/18
  higher). Caveats: different segmentation (bin2cell), platform/sequencer batch, 5-14x deeper tiles.
- **infercnvpy with the same external reference** (`results/comparators_extref_20261101_1228/`): median rho
  0.64 (0.62 with its baseline reference; higher in 14/19, Wilcoxon p 0.029, nominal). Haplodose vs infercnvpy,
  both with the external reference: median 0.71 vs 0.64, Haplodose higher in 14/20, p 0.18; altered arms right
  200/364 vs 187/357; neutral called 46/416 vs 44/413. Haplodose still worse on U0247 (0.09 vs 0.41) and U0438
  (IGF2 TD; 0.35 vs 0.75).
- **U0438 diagnosis: in-section immune/stromal references leak tumour RNA.** A plain pool of all query
  bins vs the external epithelium scores median rho 0.72 (U0438 0.88 vs Haplodose 0.35;
  `results/ucan_pool_diagnostic_20261101_1329/`), so cell selection was not the loss. The in-section
  immune and stromal bins' own profiles track the patient's tumour DNA (median rho 0.43 and 0.54; U0438
  immune 0.86), and U0438's fitted baseline was 98% immune, which cancels the tumour signal
  (`results/ucan_reference_leak_20261101_1341/`). Also: the depth slope learned on query bins tracks DNA
  (median rho 0.39 vs -0.16 on reference bins; `results/ucan_depth_slope_20261101_1331/`), since deeper
  Visium HD bins hold more tumour.
- **Epithelium-only reference** (in-section immune/stromal set to role other;
  `results/ucan_spatial_epionly_20261101_1342/`, infercnvpy on the same inputs
  `results/comparators_epionly_20261101_1342/`): Haplodose median rho 0.75 (U0438 0.89, U2964 0.91),
  altered right 237/364, neutral called 60/416; infercnvpy 0.63, 213/357, 77/413. Haplodose higher in 15/20,
  p 0.033 nominal (Holm over the 4 comparisons run on these sets: 0.13). Versus Haplodose with in-section
  references: 10/20 higher, p 0.87. Adding `--depth-source reference` (new flag, default query): 0.76,
  11/20 vs without, p 0.67; not made default. Recommendation for Visium HD: no in-section reference groups.
- **Excluding 1,682 low-dosage genes** (U-CAN/TCGA bulk slope < 0.25; `results/dosage_response_20261101_0931/`,
  slopes replicate U-CAN vs TCGA r 0.83): worse, median rho 0.60, altered right 128/371. Not adopted.

## Bulk RNA mode and U-CAN bulk benchmark (2026-11-01)

User direction (2026-11-01): bulk RNA and scRNA-seq are the main benchmarks; Visium HD is secondary; datasets
outside the in-house CRC sets are allowed.

`haplodose bulk` (`python/haplodose/bulk.py`, test `tests/test_bulk.py`): one profile per sample against a
normal panel. Pooled-normal rate, per-bin normal variance, normal-sample programmes high-passed within arms,
beta start chosen by HMM likelihood on a grid (fit_pooled alone never moved beta off 0.6 in low-purity
samples), arm calls as z against leave-one-out normals scaled by sample tau (|z| > 3). Reported log ratios
are centred by the plain median over bins.

U-CAN (`benchmarks/ucan_bulk/`): 1,062 tumours + 120 normals (exon counts, 18,522 protein-coding genes),
truth WGS FACETS per arm log2(TCN/ploidy); 1,034 tumours scored (28 FACETS purity/ploidy NA). Unit =
patient; Holm over method pairs. `results/ucan_bulk_20261101_1412/`:

| Method | median per-tumour rho | median arm AUC gain | median arm AUC loss |
|---|---|---|---|
| Haplodose segmented dose | 0.649 | 0.851 | 0.846 |
| Haplodose log ratio (median-centred) | 0.627 | 0.887 | 0.904 |
| naive gene log2 CPM ratio, arm median | 0.630 | 0.887 | 0.907 |
| infercnvpy (samples as cells) | 0.565 | 0.810 | 0.773 |

Dose beats naive per tumour (589/1,034 higher, Holm p < 1e-4); both beat infercnvpy. Diagnostics: per-gene
median bins, no programmes, mean-of-logs reference and 1-CPM pseudocount changed nothing material
(`results/ucan_bulk_variants_*`, `results/ucan_bulk_diag_ref_*`); the AUC gap (0.854 loss before) was the
HMM-fit centre (`results/ucan_bulk_diag_centre_*`), now replaced by the plain median for reported values.
TCGA pan-cancer (HaploDose2, `benchmarks/tcga_bulk/`, `results/tcga_bulk_20261101_1405/score_centred/`):
16 types with >= 10 solid normals, ABSOLUTE truth, 6,243 paired tumours. Median per-tumour rho / arm AUC
gain / loss: Haplodose log ratio 0.637 / 0.909 / 0.890; Haplodose dose 0.638 / 0.889 / 0.845; naive 0.670 /
0.926 / 0.915; infercnvpy 0.532 / 0.823 / 0.753. Naive > Haplodose pooled (Holm p < 1e-4) and in 12/16 types
(Holm < 0.05 within type); not significant in KICH, KIRP, READ, THCA. Raw STAR counts vs Xena pseudo-counts
(COAD+READ, n 527) change nothing (0.732 vs 0.737). Conclusion: at arm level, Haplodose bulk does not beat
a gene-level median log ratio.
Sub-arm (HaploDose2, `benchmarks/tcga_bulk/subarm.py`, `results/tcga_bulk_20261101_1405/subarm_20261101_1431/`):
gene-anchored ~10 Mb windows, DNA at hg19 gene positions (Ensembl GRCh37 r75). Pooled medians, dose / HD
log ratio / naive / naive smoothed (101 genes) / infercnvpy: within-arm rho (arm medians removed, n 5,055)
0.335 / 0.257 / 0.253 / 0.243 / 0.181; focal-loss recall (mean, n 3,013) 0.53 / 0.46 / 0.47 / 0.48 / 0.27;
focal-gain recall 0.58 / 0.59 / 0.59 / 0.51 / 0.50. Dose vs each method Holm < 1e-4 on within-arm rho and
loss recall; within-arm rho dose > naive in 12/16 types, Holm < 1e-4 in 9 (BRCA, COAD, HNSC, KIRC, KIRP,
LIHC, LUAD, READ, UCEC); naive higher in STAD (Holm < 1e-4). Thresholds at the 5th/95th percentile of DNA-neutral
windows. FPR-matched rerun (seeded random tie-break; `score_fprmatched_random_20261101_1442/`, pooled realised
FPR equal to target for all methods): focal-loss recall dose vs naive 0.254 vs 0.174 at 1% FPR, 0.379 vs
0.276 at 2%, 0.531 vs 0.472 at 5% (Holm < 1e-4, n 3,013, dose higher in 1,349, naive in 871), 0.603 vs 0.630
at 10%, 0.689 vs 0.781 at 20% (10-20% descriptive, not tested). Gains: no difference at any FPR (0.580 vs
0.589 at 5%). Segmentation gives cleaner strict-threshold loss calls; at loose thresholds the HMM's zeroed
weak events cap its recall.

Andor 2020 SNU-601, full-depth scDNA truth (HaploDose2, `results/andor_20261101_1224/scdna_fulldepth/SNU-601`):
1,527 cells, 15 DNA clones, 19/39 arms differ between some clone pair (chance precision 0.49). RNA clones /
difference AUROC / arms called / precision / recall: Haplodose 5 / 0.66 / 13 / 0.77 / 0.53; infercnvpy 6 /
0.77 / 24 / 0.50 / 0.63; R inferCNV 5 / 0.84 / 16 / 0.56 / 0.47; prototype 4 / 0.65 / 6 / 0.67 / 0.21.
One line, 39 arms: no ranking is supported. Other 8 lines (10 specimens) have only low-depth cell-tagged scDNA (median 4.1k fragments per cell). Simulated depth check on SNU-601 (`benchmarks/andor/simulate_lowdepth.py`, seed 20261101; `depth_check_SNU-601/compare/`): low-depth truth gives 2 clones vs 15, ARI 0.12, difference-arm recall 2/19. The 8-line table (`scores_combined/lowdepth_pooled.tsv`; 11 positives; precision 0.06-0.10 for all tools, a lower bound) cannot rank tools; no tests run. Andor result = SNU-601 only; 8-line table is a supplement sanity check.

Andor rerun on repo code 2026-11-02 (clones.py variant C + HaploDose2 fix b; frozen `dev/hd2_20261101/repo_0042`;
`results/andor_rerun_20261102_0042/`), vs the 2026-11-01 12:24 code (all changes since then, not one fix).
SNU-601 full depth (1 specimen, no test): RNA clones 5 -> 12 (DNA 15), difference AUROC 0.66 -> 0.82, arms called
13 -> 22, true 10 -> 17, precision 0.77 -> 0.77, recall 10/19 = 0.53 -> 17/19 = 0.89, profile rho 0.77 -> 0.78.
Comparators (AUROC / precision / recall): infercnvpy 0.77 / 0.50 / 0.63, R inferCNV 0.84 / 0.56 / 0.47; Haplodose
now has the best precision and recall, R inferCNV a slightly higher AUROC. Low depth, 8 lines (supplement):
pooled arms called 54 -> 97, true 4 -> 10 of 11 positives, precision 0.07 -> 0.10, median AUROC 0.89 (n 7) ->
0.97 (n 8). Cost: over-splitting of DNA-clonal lines (MKN-45 1 -> 3 RNA clones, 2 false arms; NUGC-4 1 -> 4, 6
false arms; HGC-27 4-6 clones vs DNA 2; SNU-16 6-7 vs DNA 2). Open: split calibration (see Kinker subclone).

Allele work (WGD, LOH from bulk RNA): 246,275 common exonic SNPs (1000G GRCh38, EUR AF 0.05-0.95,
biallelic; `results/ucan_allele_sites_20261101_1403/`), pileup `benchmarks/ucan_bulk/alleles/pileup.sh`.
Streaming mode (`-T`) 3.6 min / 0.5 GB per sample vs > 29 min indexed (`-R`); full cohort 1,182 samples,
`results/ucan_allele_pileup_20261101_1518/` (jobs 22596652-57). Arm minor-allele fraction without germline
genotypes: `alleles/arm_baf.py` (HWE genotype mixture from EUR AF, beta-binomial het sites, per-arm b on a
grid, per-sample overdispersion). Two test tumours: per-tumour rho of b vs FACETS purity-mixed expected minor
fraction 0.91 (U0001) and 0.86 (UM052); every FACETS LOH arm imbalanced. Scoring (`alleles/score_alleles.py`:
LOH and imbalance AUC per arm, WGD classifier allele vs dosage, IGF2 TD vs non-TD), job 22596673).
Results (`results/ucan_alleles_20261101_1523/score_20261101_185027/`; 1,034 tumours, 28 skipped; unit = patient):
- Per arm vs FACETS: median per-tumour rho 0.73; median per-arm AUC 0.93 for LOH, 0.88 for any imbalance. At a
  call threshold set to the 99th percentile of normal-sample arms, 289/23,740 DNA-balanced arms are called
  (1.2%) and 3,600/16,555 DNA-imbalanced arms (22%): strict calls catch only clear imbalance.
- WGD (out-of-fold AUC, logistic regression, 5-fold x 10 repeats, seed 0): label FACETS ploidy >= 2.7 (295
  WGD): dosage features 0.81, allele features 0.88, both 0.89; label major copy >= 2 over > 50% of genome (328):
  0.83 / 0.90 / 0.92. Adding allele features to dosage: +0.06 to +0.10 (95% paired bootstrap over patients)
  and +0.07 to +0.11. Alleles carry the WGD signal that dosage alone misses.
- IGF2 TD (78) vs non-TD (956), Holm within each source: RNA broad gain chr2 9.0% vs 1.8% (Holm p 0.008);
  RNA arm log ratio chr2 and chr20 higher in TD (Holm p < 1e-4 each); RNA WGD score 0.21 vs 0.13 (nominal p
  0.016, Holm p 0.078, not significant). DNA for comparison: gain chr2 27% vs 10% and chr20 64% vs 38% (Holm
  p < 1e-3 each), FACETS ploidy Holm p 0.13. RNA reproduces the slide-14 chr2 and chr20 gain enrichment; the
  WGD difference is in the same direction but not significant in RNA or in DNA ploidy at n = 78 TD.

GSE149224 HCT116 arm rho (Haplodose 0.25, lowest of four) is truth-limited, not a method ranking. Truth is
DepMap line minus RKO: only 7/39 arms differ for HCT116 (|diff| >= 0.2), so rho is mostly rank noise on
neutral arms; Haplodose gets 5 of the 7 (9q, 8q, 7q, 10p down, 17q up) and misses chr20. For SW480 (22
informative arms), all four tools show no 7p, 3q, 13q gain (DepMap +0.6 to +0.74 vs RKO) and none sees 18q
loss; the tools agree with each other and not with DepMap, so the GSE149224 stocks likely differ from DepMap's
[Inference]. Not used for ranking from here on.

Kinker 2020 (GSE157220; `benchmarks/kinker/`): 53,513 cells, 198 lines in 8 pools + an 8-line custom pool.
GEO has linear CPM and no cell labels: counts recovered exactly (library = 1e6 / min nonzero CPM, max
deviation 2e-11); lines assigned per cell by correlation with DepMap expression of the pool's lines
(`results/kinker_input_20261101_1535/`). 181/198 lines kept (`prepare.eligible`): 7 lines lack DepMap
expression and their cells went to neighbours (custom pool dropped; pool 10 NCIH2077 -> NCIH1048, RD), so
lines with assigned cells off the published count by > 10% are dropped. Design (`results/kinker_input_v2_20261101_1551/`):
two rotations per pool, half the lines as reference, each other line a single-line specimen; truth = DepMap CN
minus the UMI-weighted reference-line CN. Jobs 22596804 (Haplodose) and 22596806 (infercnvpy per pool);
`results/kinker_20261101_1536/` and `kinker_infercnvpy_20261101_1537/` are from cancelled runs on the
uncleaned design and should not be used.
Results (`results/kinker_compare_20261101_*/`, `benchmarks/kinker/compare.py`; n = 181 lines, one query each,
paired Wilcoxon on arm-level Spearman rho, Holm over 3 tool pairs): median rho Haplodose 0.854, pseudobulk
0.863, infercnvpy 0.813; sign on informative arms (|truth| >= 0.2, median 15 per line) 0.964, 0.967, 0.931.
Haplodose > infercnvpy on 130/181 lines (Holm p 2e-11); pseudobulk > Haplodose on 119/181 (Holm p 2e-6, median
gap 0.01, same in single- and multi-clone specimens). A single-line specimen is all one population, so its
pseudobulk is the best possible estimate; Haplodose cannot beat it here, only match it. Haplodose gives 1 tumour
clone in 115/181 specimens, 2-8 in 66 (Kinker report CNA subpopulations in some lines [Unverified], so not all
splits are errors).

Kinker subclone design (`benchmarks/kinker/implant.py`, `score_implant.py`, `run_implant.sbatch`): in each
single-line query (166 lines, >= 100 cells) one event is implanted in a random share of cells (binomial thinning;
size 50/150/500 expressed genes x loss/gain x 15%/40%, 13-14 lines per setting). Same expression state, so only
the event separates the two groups. infercnvpy runs with `--rule all` (every query cell clustered; the
held-out-reference malignancy null is meaningless when the reference is other cancer lines). Mean ARI vs the
implant split, unit = line, paired Wilcoxon, Holm over 1 pair: before the review fixes Haplodose 0.126 vs
infercnvpy 0.157 (p 0.086; `results/kinker_implant_20261101_1619/`), after 0.135 vs 0.157 (p 0.19;
`results/kinker_implant_v2code_20261101_1717/`). Haplodose leads only on 500-gene losses (0.83 vs 0.67 at 40%,
0.47 vs 0.30 at 15%); it finds almost nothing at 50-150 genes or for gains, where infercnvpy reaches 0.07-0.33.
[Inference, later found incomplete] Event proposals are whole chromosomes or halves, so a 15-bin event is diluted in its region.

Diagnosis (2026-11-01 evening). (1) Region ranking: the same chr11 (FTH1, MALAT1) and chr16 (MT2A) regions
topped the event proposals in every specimen; region totals with one median dispersion let one overdispersed
gene dominate. `event_mixtures` now uses per-bin NB likelihoods with a shared rate ratio, a per-cell NB size
factor, sliding windows (`CloneParams.event_windows` = (8, 24) proposal bins; CLI `--event-windows`) and drops
candidates whose cell split agrees >= 90% with a better one; `tests/test_numerics.py` has a test that fails on
the old code. Result (`results/kinker_implant_windows_20261101_1854/`): proposals now hit the implanted region
(TOV21G chr12, SKMEL30 chr3), but ARI fell to 0.116 vs infercnvpy 0.157 (Holm p 0.003); 500-gene losses lost
0.115 ARI on average from more accepted false splits that then starve a small true subclone below min_cells
(CAOV3: chr12 and chr10 splits at score ~10 before the true chr11 split at 48). Single-line: rho 0.854 and
115/181 one-clone lines (116 before). Scale 100k: 404 s (was 148 s), ARI unchanged 0.19. (2) Gain design flaw:
v2 gains thinned implanted cells outside the region to 2/3, so they had a third fewer UMIs; Haplodose models
depth as a covariate and reads that as depth. Design v3 (`implant.py` GAIN_MODE depth_neutral, same regions and
cells as v2) thins the whole specimen to 2/3 except the region in implanted cells. Reruns with and without
windows: jobs 22598244 / 22598245, `results/kinker_implant_v3[_nowin]_20261101_2021/`. Open: the split
acceptance test is not calibrated (review: even/odd reuse), which more candidates make worse.

v3 results (n = 166 lines, paired Wilcoxon, Holm over 1 pair): windows on 0.147 vs infercnvpy 0.165 (p 0.087),
off 0.153 (p 0.18); windows are now off by default (`event_windows = ()`). Cause of the remaining loss (MCF7,
15-bin loss in 57/143 cells): EM drained the small clone because each clone was segmented on its own and a
small clone's profile is more shrunk. Fix (`SpecimenModel.shared_profiles`, used in `em`): tumour clones take
the pooled tumour dose plus their weighted mean deviation from the pool on each run where the pooled path, own
path and resolved pairwise differences are constant. Kept version (deviation only on resolved differences):
0.169 vs 0.165 (p 0.93; `results/kinker_implant_v3_shared_20261101_2109/`), region found 0.40 vs 0.67. Also kept:
event dedup is Jaccard >= 0.7 on minority cell sets (raw agreement >= 0.9 dropped disjoint rare groups), and
`try_splits` requires every pair of parts to differ, not only the first two. Ablation on frozen copies, all on
v4 (`dev/hd1_ablate_20261101_2350/`; 166 lines implant, 181 single-line):

| variant | implant ARI | region found | single-line rho | one-clone lines | tests |
|---|---|---|---|---|---|
| A shared doses | 0.169 | 0.37 | 0.854 | 70 | 58/59 |
| B A + deviation on every run | 0.147 | 0.39 | 0.854 | 75 | 57/59 |
| C A + Jaccard + all-pairs (repo) | 0.167 | 0.36 | 0.854 | 70 | 58/59 (v4), 59/59 (v6) |
| D B + C | 0.149 | 0.38 | 0.854 | 76 | 57/59 |

infercnvpy: 0.165, region found 0.67. Where implanted events are lost (variant C, `trace_region.py`): an event
proposal overlaps the implanted region in 49% of specimens, and most of those are merged back by the odd-bin
agglomeration; 22% accepted. Tried, not adopted: variant E, a targeted held-out test of each event proposal's
region (robust z vs same-size windows, |z| >= 4, gene-consistency veto), also guarding merge_unsupported:
ARI 0.173 (noise-level vs C), one-clone lines 70 -> 61, rejected most proposals that hit real implants
(150-gene 40% gains: 0.64), and failed `test_expression_programme_is_not_a_clone` (no programme correction; a
co-regulated cluster shifts even and odd bins alike). Reading: the proposed cell groups are too mixed for any
held-out test to have power; proposal quality, not the acceptance test, limits small events.
Proposal diagnostic (`dev/hd1_proposals_20261102_0600/`, r0, 150- and 500-gene settings, 51 specimens, up to 50
event candidates kept): Leiden proposals have ARI ~0 with the implanted cells; an event proposal overlapping a
150-gene region exists in nearly every specimen, but the best one has ARI 0.15-0.29 (windows on: 0.26-0.51) and
ranks 3-9 (5 are kept). Fate follows quality: ARI > 0.6 accepted 11/13, ARI < 0.1 mostly merged. Oracle bound
(`oracle.py`; region known, 2-component binomial mixture on each cell's region UMIs, fitted without labels, real
counts): ARI 0.32-0.63 at 50 genes, 0.58-0.84 at 150, 0.80-0.95 at 500, vs Haplodose 0.00-0.34 and infercnvpy
0.00-0.32 at 50-150. So there is large headroom, and it is in localising the region. Tried, not adopted:
variant F, windows (8, 24) on and candidates ranked by gain as a robust z within their kind (half, whole,
window size): ARI 0.171 (noise vs C 0.167), better rank for 150-gene 40% gains (1.3), worse for losses (10-17),
500-gene 15% gains 0.094 -> 0.016. Same oracle on Haplodose's proposal genes only (used genes in even bins,
`oracle_even.tsv`): 150 genes 0.40-0.62 (median region UMIs ~60 vs ~190 with all genes), 500 genes 0.52-0.89.
So the even/odd split plus the gene filter costs ~0.2 ARI at 150 genes, and localisation another ~0.1-0.2.
The even/odd design is what makes the odd-bin test valid (labels chosen on even bins); using all bins needs a
different calibration, e.g. a null for the best window's mixture gain from a parametric bootstrap.
One-clone lines were 115/181 before shared doses (truth unknown: Kinker
lines can carry real subclones). `test_spatial::test_potts_improves_coherent_clones` fails on v4 for the old
own-profile code too (accuracy 0.75 < 0.85, job 22600052) and passes on v6 (22597830, 22601012): the fixture's
per-cell signal is weak by design (plain run: the copy-neutral clone takes 334 tumour cells; the chr12 split
proposal separates the subclones barely above chance, `dev/hd1_potts_20261101_2250/dbg_chr12_22600993.out`), so
the outcome flips with CPU float differences. Open. A/B noise: identical code on v4 vs v5 differed by 0.004 ARI
(HaploDose2); `--constraint` on the sbatch line does not override `--prefer` in the script, so pin in the script.

Tried, not adopted: variant G (`dev/hd1_ablate_20261101_2350/G`), event proposals scanned on every bin and accepted when the region's mixture gain beat the largest gain of its kind in 19 bin-order permutations (per-bin dispersion and cross-gene programme correlation kept, adjacency destroyed); accepted regions treated as known differences in shared doses and merges. Kinker v3, v4: ARI 0.158 vs C 0.167 (56/63, p 0.27); at/above floor 0.181 vs 0.192, below 0.050 vs 0.051; single-line gene rho 0.861, one-clone lines 71; passes the programme test; 17x slower (6 h vs 21 min). Together with E and F, three acceptance or ranking schemes on top of the current proposals did not improve small-event recovery.

Proposal floor (#5). Field rescore by floor stratum (HaploDose2): at/above 30 implanted cells (137) Haplodose
0.192 vs infercnvpy 0.180 ARI; below (29) 0.051 vs 0.093, so the gap to infercnvpy is in the below-floor
specimens. `--min-cells` experiment on repo code (clones H, v4, `dev/hd1_ablate_20261101_2350/out_R1_*_mc{15,20}`):
min-cells 30 / 20 / 15: below-floor ARI 0.031 / 0.081 / 0.150, at/above 0.195 / 0.183 / 0.184, all 166 0.167 /
0.165 / 0.178 (mc15 vs 30: 58/64, p 0.80); single-line one-clone lines 70 / 46 / 38 of 181, gene rho 0.860 in
all. Split scores of small (< 30 cell) splits do not separate helpful ones (implant below-floor, ARI +0.1:
scores 5-25) from single-line small splits (median 11.6, 90th percentile 21.4), so a size-dependent score
threshold is not justified without a calibrated null. Default stays 30; `--min-cells 15` is the option when rare
subclones matter more than extra splits.

Single-line profile gap to pseudobulk (HaploDose2, `results/singleline_decomp_20261101_2354/`): summed-count
bins let a few highly expressed genes dominate a bin (about 0.01 rho); the median of per-gene ratios within the
same bins recovers it (0.862 vs pseudobulk 0.863). Added (HaploDose2, add-only): `profiles.tsv` log_ratio_gene
and `arms.tsv` median_log_ratio_gene (per-gene log ratio, genes with E_g >= 1 in the clone, plain median in bin
or arm; bulk.py not touched); existing columns unchanged. Kinker single-line rho 0.854 -> 0.863 (181 lines,
109/72, Holm p 0.001 over 2 gates), but the column was chosen after the Kinker decomposition, so that is not
independent evidence. Independent check, Andor SNU-601 full depth (1 specimen, no test;
`results/andor_gene_20261102_0113/`): profile rho 0.765 -> 0.815, difference AUROC 0.76 -> 0.87; low depth
profile rho median 0.71 -> 0.77 (9/10 specimens higher). Costs: U-CAN spatial rho 0.64 -> 0.63 (20 patients,
p 0.52), small clones' region contrast 0.419 -> 0.411 (65 clones < 100 cells, p 0.009). Tests 61/61 on v6.
Andor numbers move with node type (SNU-601 v4: 11 clones, AUROC 0.76; v6: 12, 0.82, same code): quote ranges.
Remaining gap to pseudobulk (field run, 181 lines; d = rho gene - rho pseudobulk): no clone-count effect
(pseudobulk ahead in 39/71 one-clone vs 67/110 multi-clone lines, Fisher p 0.44; median d -0.004 in both); a
small shift over all lines, about 0.004 median, of which about 0.0015 is bin coarsening (pb_gene 0.8634 vs
pb_binmed 0.8619). Not pursued.

Kinker field comparison 2026-11-02 (`results/kinker_field_20261102_0529/`; Haplodose = repo code, clones.py
variant C; unit = cell line; paired Wilcoxon, Holm over the pairs against Haplodose (single-line) or over all
tool pairs (subclone)). Single-line, 181 lines, median rho vs DepMap arms / mean sign on informative arms:
Haplodose log_ratio_gene 0.861 / 0.966, Haplodose bins 0.854 / 0.964, pseudobulk 0.863 / 0.967, infercnvpy
0.813 / 0.931, CopyKAT 0.765 / 0.915, SCEVAN 0.710 / 0.896. Haplodose beats infercnvpy (138/43 lines, Holm p
9e-15), CopyKAT (169/12, 7e-29) and SCEVAN (174/7, 5e-30); pseudobulk is slightly ahead (106/74, Holm p 0.03).
Subclone v3, 166 lines, mean ARI / region found: Haplodose 0.167 / 0.36 (v4 run `out_C_implant`), infercnvpy
0.165 / 0.67, CopyKAT 0.089 / 0.48, SCEVAN 0.027 / 0.17. Haplodose vs infercnvpy p 0.91 (tie); Haplodose beats
CopyKAT (Holm p 0.0015) and SCEVAN (7.5e-6). Caveats: CopyKAT has no subclone output, so its clones are our
Leiden post-processing (`cnv_post`); SCEVAN's clones are its own (`--clones native`), fitted per pool, and it
called a median of 1 clone per line. All tools got other cancer lines as the normal reference (none was built
for that). infercnvpy single-line scores are from the 2026-11-01 16:09 evaluate run on the same input.
R inferCNV (job 22598581, `results/kinker_infercnv_r_20261101_2110/`; comparison with all tools in `results/kinker_field_20261102_0529/single_line_with_infercnv_r/`): single-line median rho 0.835, sign 0.954; Haplodose higher in 118/181 lines (Holm p 1.4e-6 over 5 pairs). Subclone design v3, all five tools with the 2026-11-02 scorer (`results/kinker_field_20261102_0529/implant_all5/`, ari_all, 166 lines, paired Wilcoxon, Holm over 10 tool pairs): Haplodose 0.167, infercnvpy 0.165, R inferCNV 0.163, CopyKAT 0.089, SCEVAN 0.027. Haplodose vs infercnvpy p 0.91, vs R inferCNV 0.51 (ties); Haplodose beats CopyKAT (Holm 0.0025) and SCEVAN (1.2e-5). By floor stratum (>= 30 / < 30 implanted cells; 137 / 29): Haplodose 0.192 / 0.051, infercnvpy 0.180 / 0.093, R inferCNV 0.175 / 0.107. Effect sizes for single-line (median paired rho difference, 95% bootstrap CI): vs R inferCNV +0.018 [0.009, 0.026], infercnvpy +0.033 [0.025, 0.048], CopyKAT +0.095 [0.083, 0.110], SCEVAN +0.144 [0.132, 0.161], pseudobulk -0.004 [-0.009, -0.001] (`single_line_ci/`). These are development-benchmark numbers; confirmatory claims follow `docs/heldout_protocol.md`.

Compute (2026-11-02). Rough timings from the Kinker runs, not a controlled comparison (different nodes, CPU
counts and batching): Haplodose 46 min for both rotations (2 runs of ~53k cells, 8 CPUs, v3 node), 3.1 GiB;
infercnvpy 13 min (16 per-pool runs), 4.6 GiB; SCEVAN 6.3 h, 16 GiB; CopyKAT 9.2 h, 22 GiB (both 16 CPUs, v4).
Controlled scaling benchmark: `results/scaling_20261102_0819/` (`benchmarks/scaling/`: prepare.py, run_one.sbatch,
summarise.py). Real Kinker r0 query and reference cells (45,977 available): 5k / 20k without replacement, 50k
(30,431 unique) / 100k / 1M with replacement (duplicates renamed); 3 repeats up to 50k, 1 at 100k, 1M for
Haplodose and infercnvpy only; every run on cpu,
`--constraint=v6`, 16 CPUs, 72 h cap, 200G memory cap (SCEVAN 5k and 20k first ran with 32G / 96G and ran out of
memory; rerun at 200G as repeats 11-13). Wall time from /usr/bin/time -v; peak memory = max(sacct MaxRSS of the
batch step, /usr/bin/time), since /usr/bin/time misses forked R workers and sacct samples every 30 s. Haplodose
code frozen in `code/` (08:19). Jobs in `jobs.txt`. Results so far (`summary.tsv`, median of repeats, wall min / peak GiB): Haplodose 5k 0.4 / 1.9, 20k 1.8 / 2.4, 50k 3.7 / 3.3, 100k 6.2 / 4.1, 1M 46 / 10.3; infercnvpy 0.4 / 1.3, 0.4 / 2.4, 0.5 / 4.7, 1.3 / 11.1, 1M 4.5 / 104; R inferCNV 2.7 / 4.1, 14 / 11.7, 66 / 28.9, 100k 245 / 64; CopyKAT 5k 12 / 14.7, 20k 114 / 65 (50k, 100k running); SCEVAN 5k 9.6 / 31, 20k 54 / 128, 50k and 100k out of memory at 200 GB. infercnvpy is ~10x faster than Haplodose at every size; Haplodose has the lowest peak memory from 20k up (1M: 10 vs 104 GiB). Panels: `results/figures_headline_20261102_1438/E{1,2}_scaling_*.pdf`.

Field comparators (env `cnv_field`, `slurm/build_cnv_field.sbatch`: CopyKAT 1.2.5 @01924c0, SCEVAN 1.0.3
@5a49b88, R 4.4; runner `benchmarks/comparators/run_r_field.py`, `benchmarks/kinker/run_field.sbatch`). Both use
the fit reference lines as normal cells, as for inferCNV. SCEVAN's `plotConsensusCNA` fails with >= 10
subclones (grep matches subclone1 and subclone10); `scevan_run.R` replaces it with a no-op. Runs:
`results/kinker_{copykat,scevan}_20261101_2246/`, `results/kinker_implant_v3_{copykat,scevan}_20261101_2246/`.
R inferCNV: `benchmarks/kinker/run_infercnv_r.sbatch` (pinned `--constraint=v6`: SIGILL on a v4 node).
Numbat on Kinker needs the ENA BAMs (`benchmarks/kinker/alleles/fetch_bams.sbatch`; HTTPS gave 403 pages and
truncations, now FTP with size-checked resume, job 22601406) and an hg19 1000G panel (`build_hg19_panel.sbatch`,
done: `data/external/numbat_hg19_panel/`). Cell IDs are `<barcode>-<channel>-<pool>`; ENA titles are
`Pool<channel>_<pool>`, so each line's barcodes per channel BAM are known. Pileup and phasing must be per line
(each line is a different person). Implanted subclones exist only in thinned counts, not in reads, so Numbat and
allele mode can be scored on the single-line design (vs DepMap) only, not on the subclone test.

## Efficiency review response (2026-11-01, `docs/efficiency_review_2026-11-01.md`)

All five numerical findings reproduced and were fixed; each has a test in `tests/test_numerics.py` that fails on
the old code. PLN log pmf rewritten in deviance form (was +16 at y = mu = 1e7); HMM messages renormalised per step
in scan, associative scan and Viterbi (posteriors moved up to 0.12 under a constant emission offset); programme
score gradient includes the softmax term and runs 5 Fisher steps (3 did not converge; off by 0.29 with unequal
dispersions); event mixtures solve the weighted NB score equation (`clones.nb_rate`) with a convergence stop;
pooled expectation uses UMI-weighted depth mean and variance (`SpecimenModel.depth_log_rate`; worst bin error 2%
vs 24% for the mean alone). Efficiency: Viterbi skipped where paths are unused, sparse depth-matched weights,
agglomeration pair cache, randomized SVD for programmes, power-of-two JAX chunk shapes, sparse bulk binning,
block-aligned gene reads, compact allele reading (new test). Bulk normals are now fully leave-one-out (rate,
variance, programmes); U-CAN bulk not rerun yet. Not done: even/odd independence (needs null simulation),
regional dispersion calibration, neighbour search, bulk beta-grid batching, Rust kernels. The review lists
scikit-learn 1.9.1; the env has 1.7.2. Code snapshot `results/code_snapshot_20261101_1717.tar.gz`.
Effect (new code, jobs 22597228-30): Kinker single-line median rho 0.854 (unchanged); subclone ARI 0.126 -> 0.135;
scale 100k: 148 s, 3.45 GiB, ARI 0.19, normals called malignant 0.20 (unchanged). The fixes are correct but do
not explain the subclone weakness.

## Review response 2026-11-02 (`docs/efficiency_review_2026-11-02.md`)

Split: HaploDose1 #1, #2, #7 (clones.py) and #5 (reporting strata); HaploDose2 #3, #6, #8, #9; #4 a decision (see
Decisions). HaploDose2 part, in the repo; full suite 65/66 on v4 (Potts, v4 only; job 22604409):
- #9 `kernels/estep.py`: cell log-likelihoods are computed relative to clone 0 per bin from float64 log-rate
  differences (`cell_loglik` now returns values up to a per-cell constant; every consumer uses within-cell
  differences). The review's case (+0.01 in 20 of 2,000 bins) returned 0.0; now within 2% of float64 (new tests).
- #3 `pipeline.py`: log_ratio_gene splits each bin's expectation by gene rates mixed with the specimen's fitted
  reference weights (and the in-sample normal type), not raw pooled reference counts. The review's pure-type-A
  query: max |log_ratio_gene| 0.02, was -1.0 (new test). No effect where a pool has one reference group (Kinker).
- #6 `benchmarks/kinker/score_implant.py`: ari is NaN below 2 assigned cells (was 1.0 for an empty specimen);
  primary endpoint ari_all (unassigned cells one extra label, missing specimen 0), coverage, n_no_output. The field
  implant rescore reproduces every tool's mean ARI exactly (coverage 1.0 for all; `results/review2_gates_20261102_0844/implant_rescore`).
- #8 `bulk.py`: median mode rebuilds pseudo counts without each held-out normal; fewer than 10 normals warns and
  sets run.json `arm_calls_calibrated` false (threshold unchanged). `benchmarks/ucan_spatial/dna_truth.py` also
  reports rho_log_ratio_gene.
Gate, U-CAN spatial, 20 patients (18 with tumour clones), v4, base = code of 01:13 (clones C, old E-step, pooled
gene shares) vs new = 08:44 repo (`results/review2_gates_20261102_0844/`; paired Wilcoxon, Holm over 3): bins
rho 0.635 vs 0.635 (1/4 patients higher/lower, Holm 0.16; #9 alone does not move U-CAN); gene rho 0.636 -> 0.649
(10/7, p 0.05, Holm 0.15; U0354 0.47 -> 0.64, U2761 0.62 -> 0.70); new gene vs new bins 0.649 vs 0.635 (10/8,
p 0.50). Neutral arms called 31/423 and altered right 158 vs 156/371 in both. The earlier U-CAN cost of the gene
column (7/11 against) is gone after #3. Kinker effect of #9 (+#3, a no-op with one reference group per pool),
v4, `dev/hd1_ablate_20261101_2350/out_{Hbase,C}_*` (Hbase = 08:40 repo with clones.py C): v3 implant ARI
0.1655 vs 0.1672 (166 lines, 24/17 higher/lower, p 0.92); single-line bins rho 0.8540 vs 0.8540 (13/18, p 0.38),
one-clone lines 72 vs 70. No effect on Kinker.
Also (HaploDose2, on request): `profiles.tsv` assignment_dose (the dose that set the final cell likelihoods,
`fit_specimen` doses; normal-like clone 0) next to log_dose (the clone's own HMM fit, as segments.tsv) for #2;
test in test_pipeline.py (15/15). #5 strata in score_implant.py (`below_floor` = n_implanted < 30, the default
min_cells; all 166 stay in the primary endpoint). Field implant runs (`results/review2_gates_20261102_0844/
implant_rescore_floor`), mean ari_all / region found, descriptive: 137 at or above the floor: Haplodose 0.192 /
0.42, infercnvpy 0.180 / 0.72, CopyKAT 0.102 / 0.51, SCEVAN 0.032 / 0.20; 29 below: Haplodose 0.051 / 0.07,
infercnvpy 0.093 / 0.48, CopyKAT 0.025 / 0.34, SCEVAN 0.004 / 0.03. Haplodose's deficit against infercnvpy sits
in the below-floor specimens; there it almost never reports the region (2/29).

HaploDose1 part, in the repo (clones.py = variant H; new `tests/test_review_20261102.py`):
- #1 `shared_profiles`: deviations are averaged within atomic intervals, cut at every accepted segment boundary of
  every tumour pair and at every pooled-path change. The review's probe (tumour 1 vs 2 opposite events, tumour 3
  spanning both): contrast 0.005 before, now at least 0.9 x the independent fit's 0.099 with opposite signs, and
  invariant to clone order. Test fails on the old code (job 22604390).
- #2 the deviation has the fitted programme contrast removed (as in `differences`), and segments pass the
  gene-consistency veto with MAP labels (`supported(..., log=False)`); alleles are not used there because the
  allele pools are updated after this step in `em`. Reporting: `fit_specimen` returns both `doses` (assignment)
  and `prof` (independent fit); profiles.tsv still writes the independent fit. Open: write the assignment dose
  as its own column.
- #7 `neighbour_graph` removes the cell by index and orders ties by (distance, index): 0 self-edges at
  triplicated coordinates (was 20 in the review's probe; test fails on old code). The directed row-normalized
  graph is documented as a smoothing heuristic, not a Potts mean-field gradient. Spatial event proposals use the
  moment-matched NB size (sum m)^2 / (phi sum m^2) for neighbourhood sums (was n x one cell's size).
Gate, v4, H vs Hbase: implant ARI 0.1668 vs 0.1655 (50/48, p 0.68), region found 0.39 vs 0.37; single-line bins
rho +0.0003 (47/46, p 0.43); gene rho -0.0013 (36/64 lines, p 0.019 nominal, Holm 0.056 over the 3 H tests);
one-clone lines 70 vs 72. Adopted for correctness; the benchmark effect is at most -0.001 rho.

## Package layout

| Path | Content |
|---|---|
| `crates/core/src/` | Rust (PyO3): streaming per-cell bin sums, grouped column sums (`aggregate.rs`) |
| `python/haplodose/io.py` | count readers (10x H5, h5ad, CSR H5), output guard (`HAPLODOSE_PROTECTED_PATHS`), gene and site lists |
| `python/haplodose/bins.py` | fine bins (10 genes), ribosomal exclusion, chromosome and arm boundaries, arm labels |
| `python/haplodose/data/centromeres_hg38.tsv` | p-arm ends from UCSC cytoBand (CRC copy, md5 in file header) |
| `python/haplodose/profile.py` | pooled Poisson-lognormal HMM segmentation, specimen fit (centre, beta, tau), clone difference test |
| `python/haplodose/programmes.py` | normal-cell expression programmes as per-cell nuisance terms |
| `python/haplodose/alleles.py` | phased SNP reading, normal-cell SNP filter, allelic-imbalance HMM, cell allele likelihood, vetoes |
| `python/haplodose/clones.py` | per-specimen clone model: proposals, merges, EM, splits, Potts prior |
| `python/haplodose/pipeline.py` | end-to-end run, reference check, outputs |
| `python/haplodose/kernels/` | JAX: HMM (scan on CPU, associative scan on GPU, per-position transitions), NB/PLN likelihoods, cell E-step, programme scores |
| `python/haplodose/simulate.py` | simulator: clones, depth effects, programmes, alleles with phase switches, spatial bands |
| `python/haplodose/cli.py` | `haplodose run`, `haplodose simulate` |
| `benchmarks/window_grid.py` | programme high-pass grid (done) |
| `benchmarks/scale.py` | streaming large simulation + run, records wall time and peak memory |
| `slurm/build_test.sbatch` | cargo test, editable install, pytest (compute node) |
| `slurm/window_grid.sbatch`, `slurm/scale.sbatch` | benchmarks |
| `slurm/gpu_smoke.sbatch` | prepared, **never submitted** (each GPU run needs explicit approval) |
| `logs/dbg/` | debug scripts and their logs (see last section) |
| `logs/runs/` | debug real-data runs (not published) |
| `results/` | window grid and scale benchmark outputs |

CRC side (outside this repo): `~/CRC/src/44_haplodose/` (adapter: `masks.txt`,
`allele_exclude_genes.txt`, `allele_exclude_sites.txt` with chr11:2130588),
`~/CRC/slurm/haplodose_cpu.sbatch` (staged publish with archive). Journal entries
H.2 and H.3 appended to `~/CRC/docs/notes/2026-10-30_meeting_followups.md`.

## What the model does now

1. **Binning.** Genes are summed into bins of 10 expressed genes in one streaming pass.
   Never features: user mask (CRC: IGF2 and aliases), ribosomal protein genes
   (`bins.RIBOSOMAL`), and genes no reference type matches within log 5
   (`pipeline.type_mismatch_genes`, listed in `type_mismatch_genes.tsv`).
2. **Baseline.** Mixture of reference types fitted to the pooled specimen. Per-bin
   mismatch variance from reference heterogeneity (`SpecimenModel.v`) is added to tau.
   With alleles, cells that are confidently balanced in the specimen's strong
   imbalance runs become an extra in-sample normal reference (`alleles.balanced_cells`).
3. **Programmes.** Top 10 SVD components of normal-cell residuals. Arm means are removed
   from the loadings (default `programme_window=0`), so arm-level copy number cannot be absorbed.
4. **Clone search.**
   - Leiden proposals (direct igraph Leiden; scanpy was removed from this path because
     of a scanpy/scipy `getnnz` bug), merged unless their pooled profiles differ.
   - EM over clones with NB cell likelihood plus depth slopes, programme scores and
     allele evidence.
   - Divisive splits proposed by Leiden, per-region two-component mixtures and per-run
     allele mixtures.
   - Even bins propose and odd bins test; even SNPs propose and odd SNPs test.
5. **Split and merge tests.**
   - Depth-matched pooling (coarsened exact matching on log depth) for split tests,
     merges after EM and reported differences.
   - Gene-consistency veto: one-sided sign test that most genes in a segment move.
   - Allele veto: expression differences with well-powered, unchanged allelic imbalance.
6. **Alleles.**
   - Numbat-format phased counts.
   - HMM over balanced and imbalanced (4 levels × 2 orientations) states.
   - Phase-switch probability (1 − e^(−2νd))/2 with d in cM; ν and the beta-binomial
     κ are chosen by likelihood.
   - SNPs biased in same-patient reference cells are dropped.
   - Joint states per segment, including `neutral_allele_imbalanced` (a cnLOH candidate, not a cnLOH call).
7. **Spatial.**
   - Potts prior over an 8-nearest-neighbour graph, mean-field.
   - λ chosen by how well neighbours' data-only posteriors predict each cell's own data.
   - The mean-field state is carried across E-steps and splits.
   - Proposals use neighbourhood-summed features. Cell likelihoods always use the
     cell's own counts.
8. **Outputs.** `cells.tsv.gz`, `posteriors.npy`, `clones.tsv`, `segments.tsv` (with
   allele columns and `joint_state` when alleles are given), `allele_segments.tsv`,
   `differences.tsv`, `profiles.tsv.gz`, `arms.tsv` (new, untested), `tree.tsv`,
   `log.tsv`, `reference_check.tsv`, `specimens.tsv`, `bins.tsv`, `genes.tsv.gz`,
   `type_mismatch_genes.tsv`, `bin_counts.npy`, `qc.json`.

## Decisions and the evidence for each

| Decision | Evidence |
|---|---|
| Expression programmes as nuisance terms | SMC held-out normals called aneuploid 0.121 observed vs 0.063 expected before; 0.050 vs 0.043 after (job 22577293) |
| Arm-mean removal instead of a 25-bin window | Synthetic grid (job 22589519, `results/window_grid_20261031_203023/`, 20 datasets per setting). Mean ARI / tumour clones (truth 2) / false programme differences / normals obs vs exp: off 0.14 / 3.85 / 126 / 0.013 vs 0.004; w25 0.17 / 2.45 / 25 / 0.010 vs 0.004; arm 0.27 / 1.35 / 0 / 0.004 vs 0.005. On SMC, every high-pass option could absorb at most about 5% of any real clone difference (job 22589503) |
| Exclude ribosomal protein genes | Whole-chr19 splits (SMC14/23/25, KG146M2, KG74P, KG136P) were driven by RP genes: SMC14 chr19 UMI-weighted shift +0.28, without RP genes +0.08; SMC18 chr8 control unchanged (job 22593201). Moorman clone ARI vs infercnvpy rose from 0.111 to 0.195 |
| Depth-matched difference tests | YCC155 clones were depth strata (median UMIs 7,700 vs 2,400 vs 1,770) separated by whole chromosomes |
| Exclude cell-type mismatched genes | Malignant agreement with Numbat: YCC155 0.56 → 0.985, YCC230 0.81 → 0.977; YCC230 clones 20 → 11; SMC normals 0.047/0.027 → 0.027/0.020 (jobs 22593864-74) |
| Potts λ from neighbour-predicted likelihood | λ from pseudo-likelihood on smoothed labels ran away (λ = 16 with shuffled coordinates). The current rule passes both spatial tests |
| Per-run allele splits with phase-free test | Summing all runs gave 457/54 splits driven by truncal events; requiring "clearly balanced" in the other group rejected true cnLOH splits |
| HMM transitions keep equal exit rates (no change; review 2026-11-02 #4) | `sticky_transition` uses T_ii = 1-p, T_ij = p q_j/(1-q_i), so the long-run state occupancy is proportional to q_i(1-q_i), not q_i: neutral 0.379 for the copy-number HMM (q 0.2/0.5/0.2/0.1) and 0.452 for the difference HMM (q 0.15/0.7/0.15), against the 0.5 and 0.7 written in STATE_PRIOR/DIFF_PRIOR (recomputed by HaploDose2). q sets where a state change goes, p how often; every state is left at the same rate. Valid as is; switching to an occupancy-preserving redraw would change state exit rates and need recalibration on null false calls, so it is not done to gain sensitivity |

## Tried and reverted

- `cell_doses='resolved'` (only resolved segments in the cell likelihood): lost the
  expression subclone in S1.
- Seven-state (fractional step) difference HMM: more false splits on real data (Moorman
  12 of 19 multi-clone, U0438 24 clones, depth test failed) and did not recover the
  scale subclone.
- Proposal refinement by even-bin reassignment: broke the S1 subclone tests and did not
  change the scale result.

## Latest debug results (unpublished; code as of jobs 22593864-74, last green tests 22593863)

**YCC with alleles vs native Numbat** (`logs/runs/ycc6_1101_0815/`):

| Sample | Malignant agreement | Numbat CNV segments matched | Tumour clones | Clone ARI vs Numbat |
|---|---|---|---|---|
| YCC230 | 0.977 | 15 of 19 | 11 | 0.04 |
| YCC061 | 0.513 | 20 of 26 | 3 | 0 |
| YCC097 | 0.549 | 0 of 2 | 2 | NA |
| YCC155 | 0.985 | 11 of 14 | 5 | 0 |
| YCC162 | 0.957 | 10 of 21 | 1 | 0 |

**SMC and Moorman, expression only** (`logs/runs/*_v6_1101_0815/`):

| Dataset | Multi-clone specimens | Held-out normals, observed vs expected |
|---|---|---|
| SMC | 8 of 22 | 0.027 vs 0.020 |
| Moorman | 9 of 19 | 0.007 vs 0.002 |

**Visium HD, spatial mode** (`logs/runs/spatial6_1101_0815/`). Spearman of arm-level
RNA against patient DNA log ratio, Haplodose arm median vs prototype:

| Sample | Haplodose | Prototype |
|---|---|---|
| U0354 | 0.565 | 0.564 |
| U0440 | 0.610 | 0.638 |
| U0247 | 0.230 | 0.227 |
| U0438 | 0.064 | 0.638 |

**Simulations:**
- **cnLOH clone found only with alleles.** Accuracy 0.83 with about 3 SNP reads per cell
  on the event [Estimate: the best achievable is about 0.89].
- **Spatial.** At 800 UMIs Potts reaches about 0.99 accuracy vs 0.66 without
  (debug job 22593332, earlier code).
- **Scale, 100k cells.** 11,000 genes, 4 specimens, depth 1,500 (job 22593859):
  145 s, 3.35 GiB peak. The subclone (a 250-gene half-chromosome gain) was not found
  (median ARI 0.19). 20% of normals are called malignant [Inference: information limit
  at d′ ≈ 3 with a 90:10 prior; the reference check agrees].

## Open problems (not done)

1. **Failing tests** (see the top section).
2. **Scale subclone missed.** Region mixtures rank large regions above the real event
   (top-5 candidates were chr1, chr11, chr12, chr16; the event was on chr3).
   The per-region moment fix broke S1. Needs a calibrated ranking or more candidates.
3. **Over-splitting where references do not match the query** (YCC230 11 clones vs
   Numbat 4; clone ARI vs Numbat about 0 everywhere).
4. **YCC061 and YCC097 malignant calls disagree with Numbat** (0.51, 0.55).
   No in-sample normals were found there.
5. **U0438 spatial profile does not track DNA** (0.06 vs prototype 0.64).
   The baseline is 0.997 immune and β sits at its bound.
6. **Borderline small clones are unstable across CPU types** (SMC19 on v4 Intel vs
   v6 AMD; floating-point differences of about 1e-16). A clone-stability report is
   planned, not built.
7. **`programme_overlap` is computed per clone vs normal**, so truncal events dominate it.
   It should be per clone difference.
8. **Known model limits.**
   - Events under about 60 genes cannot pass the 3-test-bin minimum.
   - Balanced 2:2 gains escape the allele veto.
   - A programme's arm-mean part is left to dosage on short arms.
   - Potts absorbs randomly scattered normal cells.

## Planned next steps

1. Restore green tests (revert or fix the region NB size). Run the closest tests, then
   the full suite.
2. Fix candidate ranking for region mixtures (scale subclone), then rerun the 100k and
   1M scale benchmarks (`slurm/scale.sbatch`, target ≤ 8 GiB peak). Record wall time
   and memory per 100k cells.
3. Rerun all real data (YCC ×5, SMC, Moorman, Visium HD ×4) on green code. Publish SMC
   and Moorman through `~/CRC/slurm/haplodose_cpu.sbatch`. Extend the CRC adapter
   for allele runs (YCC) and spatial runs (Visium HD, U0247 positions for samples
   without `positions.tsv`).
4. Clone-stability output (refit under perturbed seeds and subsamples; report how
   often each clone recurs).
5. Per-difference `programme_overlap`; move `arms.tsv` from untested to tested.
6. Benchmarks still owed from the plan:
   - comparison with the frozen prototype and infercnvpy on the synthetic grid;
   - patient-mixture benchmark;
   - Moorman arms vs MSK-IMPACT.
7. GPU (M4): submit `slurm/gpu_smoke.sbatch` only after explicit approval
   (`--prefer=b6k` after a smoke test).
8. Rewrite `README.md` (still describes M1 only) and add journal entries for this work.

## Rules this work follows

- Show the sbatch header before submitting.
- Account `innovation`, partition `cpu`, `--prefer=v6`, 90 h wall time.
- No arrays.
- Tests on compute nodes.
- CRC source data read only.
- Outputs to new directories.
- Never edit the frozen prototype (`~/CRC/src/43_cnv_caller/`).
- Each GPU run needs approval.
- RNA clone groups are not validated DNA clones; patient is the biological unit.

## Debug scripts (`logs/dbg/`)

| Script | Purpose |
|---|---|
| `dsrun.sbatch` + `summarise.py` | run a CRC prepared dataset, summary table |
| `yccrun.sbatch` + `numbat_compare.py` | YCC with alleles vs native Numbat |
| `spatrun.sbatch` + `dna_compare.py`, `dna_only.sbatch` | Visium HD run and DNA arm comparison |
| `agreement.py/.sbatch` | malignant and clone agreement with CRC infercnvpy calls |
| `chr19_genes.py` | gene-level view of whole-chromosome splits |
| `smc_window_check.py`, `window_shapes.py` | programme high-pass checks |
| `cnloh_dbg.py`, `spatial_dbg.py`, `ycc_alleles.py`, `prog_test_seeds.py`, `smc_subset.py` | simulation and allele debugging |

## Claim 2 (peak memory) result, 2026-11-03

Rule (fixed above): lowest median peak memory among all tools that finish, at every size from 20k cells up.
Summary: `results/scaling_combined_20261103_0102/summary.tsv` (benchmarks/scaling/summarise.py over symlinks: frozen 14:24 HaploDose runs plus
08:19 comparator runs; pre-freeze HaploDose runs excluded). Peak = sacct MaxRSS, median of 3 repeats (1 at 1M), cpu
v6, 16 CPUs. Unit = run.

| cells | HaploDose GiB | infercnvpy | R inferCNV | CopyKAT | SCEVAN |
|---|---|---|---|---|---|
| 20k | 2.47 | **2.44** | 11.7 | 64.7 | 128 |
| 50k | **3.52** | 4.67 | 28.9 | 171 | OOM at 200G |
| 100k | **4.50** | 9.28 | 63.7 | did not finish | OOM |
| 1M | **11.1** | 103.5 | not run | not run | not run |

**Claim 2 does NOT hold as prespecified**: at 20k, infercnvpy is 0.025 GiB (1%) lower (corrected 2026-11-03 from 0.026; medians 2.4653 vs 2.4403) (per-run HaploDose 2.27,
2.47, 2.47; infercnvpy 2.44, 2.44, 2.42). HaploDose is lowest from 50k up, and about 9x lower at 1M. This is
reportable as a descriptive result, not as a prespecified win. The 20k gap is within run-to-run spread for
HaploDose, but the rule was a point comparison of medians and is not relaxed after the fact. The GPU arm
(`haplodose_gpu`) uses more host memory than the CPU arm (20k 3.29, 100k 5.53, 1M 12.9 GiB).

## v-next gate result, 2026-11-03

`results/vnext_gate_20261102_2154/analysis/summary.json`. Both null arms fail criteria (a) and (b) and pass
(c), (d) and (e). Seed stability (a): evenodd_null 0.475, thin_null 0.464, against 0.80 needed. That matches the
frozen code (0.475). DNA-clonal lines (b): MKN-45 has 2 RNA clones in every evenodd_null seed, and NUGC-4 has
4-6. In evenodd_null, 3,645 of 3,659 null tests passed, and the null maximum was 0 in 3,466 of them. The likely
cause [Inference] is a one-clone NB null that is too clean: it lacks cell-level and regional variation beyond
the fitted programmes. v-next is not adopted. The diagnosis and a revised null are with the vnext-splits agent,
on development data only, with the same gate criteria.

SOTA count: claim 3 holds; claims 1 and 2 fail; claims 4 (needs a v-next that passes the gate) and 5 are open.
Three wins now require both 4 and 5.

## Claim 6 prespecified 2026-11-03 (user decision: CPTAC bulk backup), before any CPTAC data was downloaded

Survey from metadata only: `dev/cptac_survey_20261103_0612/survey.md` (GDC release 46.0). No CPTAC expression or
copy-number file had been downloaded when this was written.

6. **Focal (sub-arm) loss detection in tumour bulk RNA, CPTAC-3.** This replicates claim 3 on primary tumours that
   were not used in development (0 of 2,164 CPTAC case ids are in TCGA).
   - Code: HaploDose = frozen 14:24 (`results/heldout_20261102_1424/code`, hash-checked), with the same bulk flags as
     CCLE. infercnvpy = `benchmarks/tcga_bulk/subarm_infercnvpy.py`. Naive gene ratio = as in `subarm.py`.
   - Scorer: `benchmarks/tcga_bulk/subarm.py` unchanged, driven like `benchmarks/ccle_bulk/score.py`: hg38 gene
     positions (Ensembl 112), hg38 UCSC arms, random tie-break.
   - Truth: GDC CPTAC-3 WGS AscatNGS "Allele-specific Copy Number Segment" (GRCh38). `Copy_Number` gives total CN.
     Ploidy = length-weighted mean total CN over chr1-22 (the claim 3 rule). Definitions are as in claim 3: focal =
     segment shorter than 25% of its arm with |TCN - ploidy| >= 1; DNA-neutral = |TCN - ploidy| < 0.5.
   - RNA: STAR - Counts unstranded column, primary tumour only. Per case, use the RNA file whose sample list
     contains the ASCAT tumour sample. If several qualify, take the smallest GDC file_id. If none qualify, exclude
     the case. Pooled-aliquot files (2-5 tumour samples) are kept if they contain the ASCAT sample, and a
     sensitivity analysis excluding them is reported.
   - Groups: one per cohort: CCRCC, UCEC, LUAD, GBM, CM, STAD, PDAC, HNSCC, LSCC. AML is excluded (no normals);
     CPTAC-2 is excluded (no CN).
   - Normal panel: every Solid Tissue Normal STAR - Counts file of that cohort, one per case (smallest file_id),
     including matched and normal-only cases. Minimum 10 normals per cohort, otherwise the cohort is excluded
     (GBM has 17). The same panel is used for all three methods.
   - Unit: tumour (case), counting tumours with >= 1 focal-loss window. Primary: focal-loss recall at 5% realised
     FPR, HaploDose dose vs naive and vs infercnvpy. Paired Wilcoxon signed-rank; Holm over these 2.
   - SOTA criterion: both Holm p < 0.05 and mean paired difference > 0. Also reported: bootstrap 95% CI of the mean
     difference (10,000 resamples of tumours, seed 20261103). Descriptive only: per cohort, and excluding
     pooled-aliquot RNA.
   - Caveats stated in advance:
     - Tumour RNA contains stroma, and ASCAT TCN is purity-corrected while RNA is not, which lowers recall for
       every method.
     - GDC marks AscatNGS as deprecated, so the downloaded files are checksummed and kept.
     - This is a second bulk claim (same endpoint as claim 3), not a new modality.
   The scorer, its input builder and their hashes are written before any CPTAC data is opened. Data are fetched
   only after this section exists.

## Claim 4 SPECTRUM amendment, 2026-11-03, before any SPECTRUM scRNA was opened or any tool was run on it

Scorer frozen at `results/spectrum_scorer_20261103_0613/`: `eval_code_sha256.txt` passes; the Andor match is 20/20,
with patient AUROC equal to stored andor item 2 for 5 SNU-601 runs. Eligible patients = 28 of 39
(`eligible_patients.tsv`), from the DLP+ clone tables only.
- Excluded: 9 patients with fewer than 2 tumour clones; OV-007 with no differing arm.
- Added exclusion: OV-081. Its WGD clone differs from the non-WGD clones on all 44 arms, so no arm is truth 0 and
  the AUROC is undefined. This is a mechanical consequence of the endpoint, decided from DNA only.
- Compound leaf labels such as '0/2' are single published clones (cluster_info and branch_info) and are kept.
  Cells with an empty leaf_id are dropped.
Development context (Andor SNU-601, development data, not evidence for claim 4): patient AUROC is 0.779 for
HaploDose frozen, 0.845 for R inferCNV and 0.766 for infercnvpy.

## v-next gate 2 result, 2026-11-03

`results/vnext_gate2_20261103_0612/analysis/summary.json`: covariance (SigClust-style) split null, `--split-null-model
cov`. Criteria (a)-(e) are the same as gate 1. Neither arm passes.

| criterion | cov (null on splits) | cov_initial (also gates the initial groups) |
|---|---|---|
| (a) same clone count across 5 seeds, >= 0.80 of 181 lines | 0.519 fail | **0.917 pass** |
| (b) MKN-45 / NUGC-4 one clone in >= 4 of 5 seeds | MKN-45 2,2,2,2,2; NUGC-4 1,4,5,1,5: fail | MKN-45 2,2,2,2,2 (fail); NUGC-4 1,1,1,1,1 (pass) |
| (c) SNU-601 diff AUROC median >= 0.80 | 0.818 pass | 0.839 pass |
| (d) implant mean ARI >= 0.147 (n = 166) | 0.150 pass | **0.074 fail** (paired Wilcoxon vs evenodd p 3e-4, nominal) |
| (e) dosage median rho >= 0.85 every seed | pass | pass |

Gating the initial groups makes clone counts stable, but it halves subclone recovery: the stability comes from
calling fewer clones. MKN-45 gets a split-born second clone (about 8% of cells) under every null tried. That
clone cannot be checked against the scDNA, which has 771 cells and an arm sd of 0.30. This was predicted in the
gate 2 README before the run. Runtime is at most 46 min per run, with a peak of 2.8 GiB.

## v-next frozen for claim 4, 2026-11-03 (user decision, before any SPECTRUM scRNA was opened)

`results/vnext_freeze_20261103_0744/` (read-only code + SHA-256). This is the gate 2 'cov' arm:
`--split-mode evenodd --split-null-draws 19 --split-null-model cov`. It was frozen even though it failed gate
criteria (a) seed stability (0.519) and (b) MKN-45. Those failures are reported as limitations of the clone
estimand, together with the gate 2 table above. Claim 4 runs this code once on the 28 eligible SPECTRUM patients,
with the frozen scorer `results/spectrum_scorer_20261103_0613/`. The reference-cell rule and tool settings are
written below before any tool runs.
   Amendment to claim 6, 2026-11-03, made before any input, group or score stage ran (no CPTAC value read; only GDC
   metadata and md5 checks):
   (i) a case with several ASCAT segment files uses the smallest file_id;
   (ii) cohort = GDC primary site + disease type mapped to the 9 cohorts; cases fitting none are excluded (1 melanoma
   case with a urinary primary site).
   The eval-code hashes were written before the download (`results/cptac_bulk_20261103_0618/eval_code_sha256.txt`).
   Of 4,528 files, 4,508 passed md5 on the first fetch; the remaining 20 are being re-fetched (job 22626342).
   Cases whose files never verify are excluded and counted.
Caveat on the frozen 'cov' null, from the vnext-splits diagnostics (`dev/vnext_null3_20261103_0123/out/`): it keeps
a simulated true chr12 subclone only narrowly (null max 7.40 vs observed 7.86). It only partly stops false splits:
it rejects the first splits in NUGC-4 and SNU-16, but not in MKN-45. The changed code is
`python/haplodose/clones.py` (`region_effect_covariance`), `cli.py` and `tests/test_split_null.py`. pytest on v6
(job 22624614): 95 passed, 1 skipped. Not committed.

## Claim 4 SPECTRUM run rules, 2026-11-03 (fixed before CellRanger or any CNV tool ran on SPECTRUM scRNA)

Plan: `dev/spectrum_run_plan_20261103_0750/plan.md` (code SHA-256 in `code_sha256.txt` there, checked at every job
start). Only metadata has been read: the GEO cell table, SRA run lists and DLP+ clone tables.
- **Patients:** n = **21**, not the 28 above. Seven DNA-eligible patients (OV-004, 006, 044, 046, 129, 133, 139)
  have no CD45-negative or UNSORTED tumour scRNA, so they fail the "tumour scRNA" part of the rule.
  `paired_test.py` gets `patients_with_scrna.tsv`. Samples: 71 of 74. Three samples with no published cell labels
  are excluded; each of those patients keeps other samples.
- **Reference cells:** the same for every tool. They are the authors' published non-cancer labels (T, B, plasma,
  myeloid, dendritic, mast, endothelial, fibroblast), pooled over the patient's samples, up to 300 per type
  (seeded), split 50/50 fit/holdout.
- **Query cells:** cells labelled Ovarian.cancer.cell, up to 3,000 per sample (seeded). Each tool's own
  malignant call then applies as in the frozen scorer: comparators use `--rule cluster`, HaploDose uses its own call.
- **Counts:** CellRanger 8.0.1, GRCh38-2020-A, at the authors' barcodes.
- **Tools:** HaploDose v-next frozen (`vnext_freeze_20261103_0744`, cov flags, seed 20261031, no alleles) = primary.
  Comparators: infercnvpy, CopyKAT, SCEVAN, Numbat (GRCh38 1000G phase 3 panel; this is a deviation from "same
  panel as Kinker", which is hg19). XClone and CalicoST are not installed, so they are not run.
  HaploDose with `--alleles` is descriptive only (the allele path was not gated).
- **Failures:** a tool that fails or makes no call for a patient scores AUROC 0.5. Failures are counted per tool
  (frozen `paired_test.py`). Expected in advance: Numbat on OV-031, which has 43 reference cells.
- **Test:** unchanged. Paired Wilcoxon over the 21 patients, HaploDose vs each of the 4 comparators, Holm over 4.
  SOTA = every Holm p < 0.05 with the bootstrap 95% CI of the mean difference above 0.
- **Caveat [Unverified]:** how the authors labelled cancer cells (markers vs CNV inference) is not checked. If CNV
  signals were used, the query set favours CNV-rich cells for every tool alike.

## Claim 6 (CPTAC-3 bulk) result, 2026-11-03

Run: results/cptac_bulk_20261103_0618 (score job 22627305, frozen 14:24 code, eval code hashed before any CPTAC
file was opened). Primary as prespecified: focal-loss recall at 5% FPR, paired Wilcoxon over tumours with >= 1
focal-loss window, Holm over 2 (hd_dose vs naive, hd_dose vs infercnvpy). n = 618 tumours (unit: CPTAC-3 case) of
1,553 scored across 9 cohorts; all 1,553 input tumours were scored (the survey listed 1,579; 26 dropped at input
building by the aliquot/ASCAT rules).

| comparator | mean HD | mean comp | diff | 95% CI (bootstrap) | HD>/</= | Holm p |
|---|---|---|---|---|---|---|
| naive | 0.542 | 0.514 | +0.028 | -0.003 to 0.060 | 229/184/205 | 0.025 |
| infercnvpy | 0.542 | 0.407 | +0.136 | 0.109 to 0.163 | 279/107/232 | 4.3e-21 |

Claim 6 HOLDS under the prespecified rule (both Holm p < 0.05, both mean diff > 0). Caveat to report: the margin
over the naive gene ratio is small and its bootstrap CI of the mean crosses 0 (the Wilcoxon tests ranks, the CI
tests the mean). Per-cohort results are descriptive (e.g. CCRCC vs naive -0.011). Same sign excluding pooled-RNA
cases (n = 508, Holm p 0.033 and < 1e-15). Test recomputed independently from per_patient.tsv.gz: p match.
Correction 2026-11-03 (from docs/scientific_code_review_2026-11-03.md F10): the claim 6 table first gave 2e-21 for
infercnvpy under "Holm p"; that was the raw p. Holm-adjusted p = 4.3e-21. Decision unchanged.

## Claim 4 amendment 2026-11-03 13:30 (before any CNV tool ran on SPECTRUM scRNA)

Trigger: docs/scientific_code_review_2026-11-03.md F6. The run rules above say "Holm over 4", but the frozen
paired_test.py put a comparator in the Holm family only if it scored at least one patient, so a comparator that
failed everywhere dropped out and SOTA could be declared against a subset. Reconciled to the stricter written rule:
- Family fixed at infercnvpy, CopyKAT, SCEVAN, Numbat, whatever they output. XClone and CalicoST removed from the list
  (not installed, never part of the claim).
- New: SOTA also requires each of the 4 to have scored output for >= 50% of the 21 patients. Otherwise claim 4 is
  "not established", since a comparator that fails everywhere sits at 0.5 and is beaten trivially.
- New, descriptive: the same tests over patients where both tools scored (no 0.5 fill). Completion per tool is
  reported next to the test.
Not changed: endpoint, 0.5 fill, Wilcoxon, CI rule, patients. Old script kept as
dev/spectrum_run_plan_20261103_0750/paired_test_20261103_0613.py. Validation and new hashes:
results/spectrum_scorer_20261103_amend_1328.
Review F5 (a missing arm counts as 0 difference, which can raise AUROC) is not changed in the primary. A
common-arm sensitivity is to be added before any SPECTRUM tool output exists, or else disclosed as post hoc.
Pointer: benchmarks/spectrum/run_tool.sbatch SCORER now results/spectrum_scorer_20261103_amend_1328; plan hashes redone (old: code_sha256_20261103_1340_pre_amend.txt).

## Claim 4 amendment 2, 2026-11-03 13:38 (before any CNV tool ran on SPECTRUM scRNA)

Review F5 addressed as a descriptive sensitivity, primary unchanged. score.py writes a 'measured' flag per arm
(>= 2 kept tumour clones with a value); paired_test.py adds `sensitivity_common_arms`: per patient, AUROC over the
arms both tools measured, Wilcoxon over patients where both are defined, Holm over the 4. Validation 22/22:
results/spectrum_scorer_20261103_amend2_1338 (job 22632075). run_tool.sbatch SCORER now points there; plan hashes
redone (old: code_sha256_20261103_1345_pre_amend2.txt). Old scripts kept in dev/spectrum_run_plan_20261103_0750/.
Review F1 (NaN assignments after a small clone is dropped): the frozen code cannot change, so after the SPECTRUM HaploDose runs, count non-finite posteriors in each cells.tsv.gz and report them next to the result. Fixed before any output.

## Kinker Numbat comparator, development data, 2026-11-03 (job 22612086)

results/kinker_numbat_20261102/score_all_20261103_1424 (log: score_all_22612086.out). Numbat on 181 single-line
specimens: 146 ok, 26 no clones after the size filter, 9 failed. Line-level dosage rho (unit = cell line), HaploDose
(frozen field run) vs Numbat over the 146 scored lines (35 missing): median 0.861 vs 0.818 (the same HaploDose median 0.861 appears in every pair, so it may be over all 181 lines; not checked), median paired diff +0.031 (95% CI 0.027 to 0.050), HaploDose higher in 116, Numbat in
30, Holm p 8.5e-15 (Holm over the 6 pairs with haplodose). Clone false splits: Numbat reports a median of 4 tumour
clones (IQR 3-6) per single-line specimen where 1 is expected; one clone in 7 of 146. Development evidence only
(Kinker was used in development); not a SOTA claim.
Claim 4 run started 2026-11-03 20:28 in waves (operational only, analysis unchanged): results/spectrum_run_20261103_2028/README.md. Wave 1 = 7 patients; watcher job 22642637 submits the rest and the paired test.

## Claim 5 frozen code, 2026-11-03 ~20:40 (user: "use whatever is best"; before any HaploDose run on Schmid data)

Primary: v-next freeze results/vnext_freeze_20261103_0744/code with the claim 4 flags
(--split-mode evenodd --split-null-draws 19 --split-null-model cov), seed 20261103, other settings per
benchmarks/schmid/PLAN.md item 4 (value log_ratio_gene, reference holdout 0.3, CPU). Reason: same method as claim 4,
so the paper reports one current method for the single-cell claims. Descriptive sensitivity, not part of the claim:
the 14:24 freeze (results/heldout_20261102_1424/code) with CLI defaults. The v-next README line "claim 4 only" is
superseded by this entry. Gate (amendment d/g) still precedes any HaploDose scoring.
Clarification (before any run started): in the v-next primary arm, seed 20261103 is passed both to the HaploDose CLI (--seed) and to the wrapper's reference split. This differs from claim 4, which used the CLI default 20261031; seeds are arbitrary and both were fixed before any Schmid run. The 14:24 sensitivity arm keeps the CLI default seed. Jobs: v-next 22642692-95, 14:24 22642696-99.

### Claim 5 gate deviations, 2026-11-04 (written before the gate rerun; no comparator or HaploDose value read)

First gate runs 22642665 (BCC06), 22642666 (BCC06post), 22642667 (MCF7) ended COMPLETED in Slurm, but snakemake
failed in 4 rules on all three, so no evaluation ran and there is no gate table. copyKat, SCEVAN, BAFExtract and the
Numbat pileup succeeded. Fixes (team-lead approved 1, 3, 4; 2 waits for the user):
1. inferCNV: reticulate loaded `~/.local/share/r-miniconda` python (no igraph). `run_gate.sbatch` now sets
   `RETICULATE_PYTHON` to the schmid_r41 python (igraph 0.9.9, leidenalg 0.8.9). No user-level config changed.
2. Numbat: `run_numbat` fails in plotting (`Ops.data.frame`), most likely ggplot2 3.5.0 vs their 3.4.1 [Inference].
   Not fixed; the reinstall waits for the user. Until then the gate reruns with `NO_NUMBAT=1`: their evaluation rule is
   copied with only the numbat_obj input removed (`benchmarks/schmid/gate/make_nonumbat_rule.py`; their
   `evaluate_results.R` skips Numbat when that input is NULL). Numbat bins then do not enter the combined bin intersection,
   so those rows are provisional; the gate (amendment g) is final only after a run with Numbat.
3./4. CONICSmat and CaSpER annotations: biomaRt 2.50.3 cannot query Ensembl now (HTTP 404; BiocFileCache/dbplyr error).
   Release: Ensembl 113 (`oct2024.archive.ensembl.org`). Basis: their commit 41dfa38 (2025-03-25) set the CONICSmat host
   to www.ensembl.org "#currently: GRCh38.p14", and 113 was current on that date; the CaSpER script uses the www default
   with no date, so 113 is used for both [Inference]. Raw BioMart tables (query XML, URL, date, md5) are downloaded once and
   subset to the genes as `getBM(filters = "hgnc_symbol")` would; the code after getBM is copied verbatim
   (`benchmarks/schmid/annot/`). Assumption: the unfiltered table keeps the server's row order that getBM would return,
   which matters only for duplicated symbols [Inference].
Also: `run_gate.sbatch` now exits non-zero when snakemake fails and writes `gate_comparison.tsv` itself.

Update 2026-11-04 (after the no-Numbat reruns 22645303/4 failed; before the final gate runs; no gate value read):
- Fix 1 replaced. reticulate 1.27 ignores RETICULATE_PYTHON's conda handling here: for a conda python it takes the conda
  binary from the env's `conda-meta/history` (`~/.local/bin/conda`, absent). The gate now points RETICULATE_PYTHON at a
  job-local symlink to the env python (`OUT/pyshim/bin/python`), which reticulate treats as a plain python. Tested in
  job 22647316 (igraph 0.9.9 and leidenalg import). No user-level config changed.
- Fix 2 (user-approved): the Numbat error is patchwork 1.1.2 (their yml) with ggplot2 3.5.0 (reproduced in 22647316
  with a guides-collect plot). ggplot2 3.4.1 was approved first and then withdrawn, because ggraph 2.2.0 (needed by
  numbat) requires ggplot2 >= 3.5.0. Installed instead: patchwork 1.2.0 from the Posit CRAN snapshot 2024-02-29, in
  schmid_r41 only; ggplot2 stays 3.5.0. Package list and sessionInfo before and after are in the fix run dir. Deviation
  from their yml (patchwork 1.1.2).
- Fix 5 (approved): CaSpER's `casper.R` reads `@assays$RNA@data`, which Seurat 5.0.1 (Assay5) lacks. The gate's working
  copy of `casper.R` gets one line after `library(Seurat)`: `options(Seurat.object.assay.version = "v3")` (tested in
  22647306). The original is kept as `casper.R.orig` and the diff is written to each gate dir.
- The final gate runs include Numbat and their unchanged evaluation rule (not NO_NUMBAT). The no-Numbat route is dropped.

### Claim 5 gene-position file fix, 2026-11-04 ~09:30 (team-lead decision; no evaluation output had been written)
- Gate 22647673 (BCC06post) failed only in `evaluate_cnv_prediction_wes`; every method rule had finished. Their
  `lib.R` reads gene positions with `data.table::fread` (1.14.8, as in their yml). fread stops at line 82 of
  `external/annotations/hg38_gencode_v27.txt` (md5 0c08b487706b2b8362470570d5b9eeb8). That line (TNFRSF18) has one
  trailing tab: 5 fields, while the other 56,608 lines have 4. fread keeps 81 genes, then IRanges fails on NA.
- New copy `external/annotations/hg38_gencode_v27_fixtab_20261104.txt` (md5 53909a1213b591c12894222cc47b6b03) strips
  that tab. The diff is that one line only, with 56,609 lines in both files. md5s are in
  `hg38_gencode_v27_fixtab_20261104.md5`. The original file is unchanged.
- Only the evaluation rules use the copy: the 14 `infercnv_gene_pos` inputs in the gate's working-copy `workflow.sm`.
  The method rules (inferCNV `gene_pos`) keep the original. R `read.table` (header=FALSE, row.names=1, sep="\t"), the
  call inferCNV makes [Inference: from the inferCNV gene-order reader, not re-read here], gives identical data frames
  for both files (checked 2026-11-04). So the finished method outputs are unaffected. HaploDose read the original
  file and is left as run.
- Repair: `benchmarks/schmid/gate/rerun_eval.sbatch` reruns only the evaluation rule in an existing gate dir. It
  stops unless a snakemake dry run plans exactly that one job. BCC06post: job 22650125. BCC06, MCF7 and COLO320 get the
  same repair if their evaluation fails the same way. `run_gate.sbatch` now applies the copy from the start (diff
  written to OUT/workflow_fixtab.diff). Gates submitted earlier keep the script they were submitted with. So MM199
  gate 22649741 is replaced by 22650126, which is held like its predecessor.

### Claim 5 MM199 control-cluster rule, 2026-11-04 ~09:00 (team-lead decision; written before prep 22649737 ran; no marker table, comparator output, HaploDose output or truth file read)
- Their `prepare_MM_sciCNV_scran.Rmd` (lines 249-250) labels Seurat clusters `c(4,6,8)` at resolution 0.5 as
  "control" and all other clusters as "MM". It states no numeric threshold. It lists markers only: "Specific MM marker
  genes: B2M, CCND1, Immunoglobulins (IGHM, IGKC)", "Additional marker used in the Numbat paper for MM: MZB1",
  "Additional marker used in the Casper paper for MM: SCD1, CD38" (lines 227-233). So the thresholds below are ours.
- Our clustering (CellRanger 8.0.1, Seurat 5.0.1) may number clusters differently. Rule, fixed now: the paper's control
  is the non-plasma immune clusters. A cluster is control if its fraction of PTPRC-positive cells (count > 0) is
  >= 0.5 and its SDC1-positive and MZB1-positive fractions are each < 0.2. These fractions come from
  `MM_cluster_marker_check.tsv` written by `benchmarks/schmid/prep/prepare_MM.R`, computed on the cells kept after
  their QC.
- If this marker set equals `{4,6,8}` exactly, their labels stand. Otherwise the marker-defined set is used for all
  comparator tools and both HaploDose arms. This is recorded as a deviation, and the HaploDose jobs (22649739, 22649740)
  and the gate (22649741) are rerun on the relabelled input in new dirs. If the marker set is empty, nothing is run on
  MM199 and the case is reported.
- No truth file is read for this decision.
- Outcome, 2026-11-04 10:22 (prep 22649737; `results/schmid_inputs_MM_20261104_0849/MM_cluster_marker_check.tsv`):
  the rule gives control = {8} (32 cells; PTPRC 0.78, SDC1 0.06, MZB1 0.13). Their {4,6,8} does not match. Cluster 4
  is plasma-like (MZB1 1.00, SDC1 0.70, PTPRC 0.01), so the numbering has shifted. Clusters 5 and 6 are PTPRC-high and
  SDC1-low but fail MZB1 < 0.2 (0.38, 0.25), so the rule labels them MM. Applied as written: control = {8} for all
  tools and both HaploDose arms. Deviation. The rule is not changed after seeing the table.
  Rerun: prep 22651203 (`results/schmid_inputs_MM_ctrl8_20261104_1023`), then determinism check 22651205 (count
  matrix byte-identical and the same marker table as 22649737). Then HaploDose 22651206 (vnext) and 22651207
  (sens1424), and gate 22651208. The annotation from 22650110 is reused: it depends on genes only, and the check
  requires an identical matrix. Superseded held jobs 22649739, 22649740, 22650126 (with 22649738, 22649741) need a
  user `scancel`.

### Claim 4 run incident, 2026-11-03 ~21:00 (no outputs read)

- Scoring failed for the 14 wave-1 haplodose and infercnvpy jobs. `results/spectrum_scorer_20261103_amend2_1338/` had no
  `dna_clones/` (the amendment validation jobs skipped `dna_clones.py`). The tools had finished. Fix: `dna_clones/` in both
  amend dirs is now a symlink to `results/spectrum_scorer_20261103_0613/dna_clones` (`dna_clones.py` unchanged, DNA-only).
  Job 22642879 (`results/spectrum_run_20261103_2028/rescore_wave1.sbatch`) runs `score.py` on the recorded run dirs; no
  tool is rerun. Later waves use the symlink directly.
- infercnvpy on SPECTRUM-OV-031 failed by the frozen runner rule "at least 20 fit and 20 held-out reference cells" (the
  patient has 43 reference cells; HaploDose ran with 24/19). By the protocol it scores as no_output. This rule is in our
  wrapper, not in infercnvpy itself, so the result is reported with `sensitivity_both_scored` alongside the primary. The rule
  is not changed after the run.

### Package fixes from review 2026-11-03 (F1, F7, F8), development code only

Not in any frozen copy: claims 4 and 5 run `results/vnext_freeze_20261103_0744/code` (SHA-256 checked), which these edits
do not touch.
- F1: in `clones.py` EM, cells whose responsibilities sum to 0 or are non-finite after a clone is dropped are
  renormalised in log space. A cell with no finite likelihood falls back to the mixture weights. All other cells are
  unchanged. These cells skip the spatial smoothing term.
- F7/F8: fractional counts are rejected in `io.py`, `clones.thin_counts` and allele AD/DP. `bulk.py` checks sample IDs,
  negative, non-finite and zero-total samples. `pipeline.check_identities` rejects empty query IDs and specimens that
  span more than one patient or dataset (reference cells are not checked).
- Tests: `tests/test_review_20261103.py` (4 new). Job 22642671, ne1dc6-010: 99 passed, 1 skipped; the new tests
  fail 4/4 on the pre-fix copy. Logs: `dev/review_fixes_20261103_2032/`.
- 2026-11-04 ~08:35: patchwork 1.2.0 installed in schmid_r41 by job 22649569 (`external/schmid_env/fix_patchwork_20261104_0828/`;
  `r_packages.diff` shows patchwork 1.1.2 -> 1.2.0 as the only change; ggplot2 3.5.0; patchwork smoke test passed). The job
  exited 1 only because its Numbat smoke test called snakemake with a wildcard target (a test bug, not an env failure). The
  four final gates 22647672-75 had their dependency removed by hand; the Numbat rule inside each gate is the real check.
- **Claim 5 gate, BCC06post: PASS** (2026-11-04, `results/schmid_gate3_BCC06post_20261104_0641/gate_comparison.tsv`,
  evaluation rerun 22650125 with the fixtab gene file). Rule (g): every comparator output within 0.05 of the published value
  on Pearson, AUC gains and AUC losses. 27/27 within 0.05 (9 outputs x 3 metrics); largest gap CONICSmat Pearson 0.041
  (0.218 vs 0.177); all other gaps <= 0.009. HaploDose for BCC06post may now be scored; that requires the other primary
  lines for the claim itself (best mean rank over 3 metrics x 4 primary datasets).
- **Claim 5 MM199 control sensitivity, declared 2026-11-04 ~10:35, after the marker table and before any MM199 gate or
  HaploDose number exists.** The prespecified rule (PTPRC >= 0.5, SDC1 < 0.2, MZB1 < 0.2) gave control = {8}, 32 cells,
  not the authors' three clusters. Clusters 5 (178 cells) and 6 (152) are PTPRC-high and SDC1-low but MZB1 0.38 and 0.25;
  MZB1 is also expressed in some non-plasma B cells, so the MZB1 cutoff I set was probably too strict [Inference]. The
  primary stays as prespecified: control = {8}. One sensitivity is added, fixed now: control = {5,6,8} (PTPRC >= 0.5 and
  SDC1 < 0.2, MZB1 dropped). Both run through the gate and HaploDose with everything else identical, and both are reported.
  Neither is chosen by its result: the primary is {8} whatever the gate shows.
- Exposure note 2026-11-04 ~10:40: while checking the file format, an agent printed the header and first 2 gene rows of the
  BCC06post HaploDose vnext gene_by_clone output (two clone log ratios). No truth file was read, no metric computed, and
  no decision depends on it. Scoring job prepared, not submitted: `benchmarks/schmid/score/run_score_all.sbatch`
  (combined scorer primary, per-method sensitivity; refuses to run unless every gate row is within 0.05 or a written override
  exists). Fixture test 22651292 replaces HaploDose with copyKat's pseudobulk and must match copyKat's rows to 1e-9.
- **Claim 5 gate, BCC06: 26/27 within 0.05** (2026-11-04, `results/schmid_gate3_BCC06_20261104_0641/gate_comparison.tsv`,
  evaluation rerun 22651591 with the fixtab gene file). One miss: Numbat CNV Pearson 0.184 vs published 0.235 (gap 0.051,
  just past the 0.05 rule). Per rule (g), numbat_cnv on BCC06 is still ranked, flagged gate_miss, and the gap reported.
  No rerun of Numbat is made to close the gap (the rule gives a second run only to SCEVAN). Cause not investigated before
  scoring; candidates are Numbat run-to-run variation and our env deviations (patchwork 1.2.0, CellRanger 8.0.1) [Speculation].
- **Claim 5 MCF7 Numbat failure and missing-output ranking rule, 2026-11-04 ~16:45 (team-lead decision; written before any
  HaploDose number on any Schmid line exists; MCF7 is line-overlap, not primary).** In gate 22647674 `run_numbat` failed
  on MCF7: "Found 3 normal cells.. Running HMMs on 7 cell groups.. Error in smooth_segs(., min_genes = min_genes) : No
  segments containing more than 10 genes for CHROM 13,18,21,22." Numbat is not rerun. Deviation: MCF7 is evaluated with
  a copy of the authors' rule `evaluate_cnv_prediction_wgs` with only the `numbat_obj` input removed
  (`benchmarks/schmid/gate/make_nonumbat_rule.py`, output to `evaluation_nonumbat/`); their scorer handles a NULL Numbat
  object. The evaluation-only guard is the same as `rerun_eval.sbatch`: a dry run must plan exactly one job of the
  no-Numbat rule, so no method reruns. The HaploDose MCF7 scoring also drops the Numbat input. In `gate_check`, Numbat
  (CNV) and Numbat (Expr) on MCF7 are gate_miss with no output (ours = NA; published values exist: Pearson 0.300 and
  0.615). Ranking rule, new (PLAN.md and the claim 5 text above define none; the claim 4 "failure scores AUROC 0.5" rule
  is a different claim and is not used here): an output that is missing for a dataset and metric ranks last on that
  dataset and metric; ties among several missing outputs share the last ranks (average rank). Correction to the 10:40
  note above: the scorer no longer requires every gate row within 0.05; per (d)/(g) it stops only if no primary line
  passes fully, keeps misses ranked and flagged, and skips a line with no gate table.
- **Claim 5 gate, MM199 sensitivity control {5,6,8}: 25/27 within 0.05** (2026-11-04, `results/schmid_gate3_MM_ctrl568_20261104_1027/gate_comparison.tsv`).
  Misses: CONICSmat Pearson 0.391 vs 0.539 (gap 0.148); CaSpER AUC gains 0.629 vs 0.693 (gap 0.064). Both flagged and ranked
  under rule (g). This is the sensitivity arm; the primary MM199 arm (control {8}, gate 22651208) is still running. The
  primary stays {8} whichever arm reproduces better.
- Claim 4 run, 2026-11-04 ~16:50: SCEVAN on SPECTRUM-OV-009 (job 22651351) was OOM-killed at its frozen 128G after 1.6 min
  (R child killed, 5 oom_kill events). score.py wrote status no_output, which paired_test scores 0.5 by the protocol. Not rerun
  with more memory: the protocol runs each tool once at fixed resources, and a rerun chosen after the failure would be a
  post-hoc change. It counts against SCEVAN's completion (MIN_COMPLETION 0.5) and is excluded in `sensitivity_both_scored`.
- **Claim 5 gate, MM199 primary control {8}: 21/27 within 0.05** (2026-11-04, `results/schmid_gate3_MM_20261104_1023/gate_comparison.tsv`).
  All misses are CONICSmat and CaSpER, on every metric: CONICSmat Pearson 0.109 vs 0.539, AUC gains 0.794 vs 0.868, losses
  0.665 vs 0.792; CaSpER Pearson 0.039 vs 0.391, gains 0.640 vs 0.693, losses 0.665 vs 0.773. inferCNV, copyKat, SCEVAN
  and Numbat (both outputs) pass. Under rule (g) both are ranked and flagged gate_miss. The sensitivity arm {5,6,8} misses
  less (CONICSmat Pearson only 0.148, CaSpER gains 0.064), consistent with a 32-cell reference hurting the two tools that
  use the reference most directly [Inference]. The primary stays {8}, as fixed before either gate result existed.
- **Claim 5 gate, MCF7 (line-overlap): 21/21 non-Numbat rows within 0.05** (largest gap SCEVAN CNV AUC losses 0.012;
  `results/schmid_gate3_MCF7_20261104_0641/gate_comparison.tsv`, no-Numbat evaluation 22653989). Numbat failed (3 normal
  cells; smooth_segs error), so its 6 rows are NA and count as gate misses, ranked last per the 16:45 rule.

### Claim 4 result (SPECTRUM, read 2026-11-04 22:10, paired test job 22652272)

Prespecified outcome: **claim 4 fails.** HaploDose is not superior to any of the 4 required comparators.
Source: `results/spectrum_run_20261103_2028/primary/{primary.json,per_patient.tsv}`.

- Endpoint: clone-difference AUROC per patient (score.py level=patient). n = 21 patients; unit = patient.
- Test: paired two-sided Wilcoxon over patients; Holm over the 4 required comparators (infercnvpy, CopyKAT, SCEVAN, Numbat). Missing output scored 0.5. Bootstrap 95% CI of the mean difference, 10,000 patient resamples, seed 20261103.
- Completion: all 4 comparators passed MIN_COMPLETION 0.5. Failures scored 0.5: infercnvpy 1 (OV-031), SCEVAN 1 (OV-009, OOM), Numbat 1 (OV-031), CopyKAT 0, HaploDose 0.
- Mean AUROC: HaploDose 0.647, CopyKAT 0.646, SCEVAN 0.645, Numbat 0.626, infercnvpy 0.601.

| vs | mean diff | 95% CI | HD higher / lower | raw p | Holm p |
|---|---|---|---|---|---|
| infercnvpy | +0.047 | -0.003 to 0.097 | 14 / 6 | 0.100 | 0.402 |
| CopyKAT | +0.001 | -0.043 to 0.048 | 10 / 11 | 0.973 | 1.0 |
| SCEVAN | +0.003 | -0.115 to 0.100 | 13 / 8 | 0.374 | 1.0 |
| Numbat | +0.021 | -0.025 to 0.070 | 11 / 8 | 0.573 | 1.0 |

- Descriptive sensitivities (both-scored only; common arms) give the same conclusion: no Holm p below 0.68.
- Per-patient AUROC below 0.5 for HaploDose: OV-002 (0.305) and OV-050 (0.047; SCEVAN 0.977 on the same patient). Reported, not investigated for tuning; SPECTRUM stays an evaluation set.
- Review F1 check: non-finite clone_posterior or p_normal among assigned cells = 0 in all 42 HaploDose runs (21 primary, 21 alleles). assignment_margin is empty for 8,870 cells in 5 single-clone specimens (OV-009 1, OV-065 3, OV-075 1); the margin is undefined with one clone, so these are not F1 NaNs.
- Interpretation: on this endpoint HaploDose is at parity with CopyKAT and SCEVAN and numerically above infercnvpy and Numbat, without statistical support. Claims status now: holds 3 (CCLE bulk), 6 (CPTAC-3 bulk); failed 1, 2, 4; claim 5 pending.

### Claim 5 scorer fixture, 2026-11-05 00:30 (operational, no HaploDose output read)

- Fixture v5 (22653951) hit its 8 h limit in case A; no verdict. Each scored run takes about 2.3 h per HaploDose arm (2 arms), so case A alone needs about 18 h.
- Fixture v6 (22690386) submitted through `submit_frozen.sh` with the current scorer (adds the Numbat-NA case and case C, MCF7 line_overlap), `--time=72:00:00` on the command line. Code copy and sha256: `dev/schmid_claim5_20261103/score_fixture_v6_20261105_0030/code/`.
- The all-lines scoring job (6 runs x 2 arms, about 28 h [Estimate]) will be submitted with `--time=90:00:00` on the command line; the 24 h in the `run_score_all.sbatch` header is too short. Timing only; no rule changes.

### Claim 5 scoring fan-out, 2026-11-05 (user: "full fan out jobs as much as you need"), written before submission

- Scoring is submitted now, in parallel with fixture v6, one `run_score_all.sbatch` job per run_id, each from its own
  `submit_frozen.sh` copy. The scorer tree is byte-identical to the fixture v6 copy (`diff -r` clean at submission).
- Blindness: no shard output is opened until fixture v6 reports FIXTURE_PASS on that same code. If v6 fails, every
  shard is discarded unread and rerun after the fix; that is recorded here.
- Gate rule (d)/(g) is unchanged. The stop rule needs one primary run passing every gate row. For shards whose own run
  does not pass every row, the map adds an anchor row `BCC06post_gateanchor` (the real BCC06post gate table, which
  passes 27/27) with HaploDose paths `none`, so the anchor is listed in `skipped_runs.tsv` as haplodose_output_missing
  and never scored. The run's own gate rows are flagged exactly as in a single job.
- COLO320 gets its shard when its gate table exists.
- 2026-11-05 11:26: COLO320 gate 22647675 finished every method rule, then its evaluation stopped on the GENCODE line 82
  trailing tab (same as BCC06/BCC06post; the job predates the fixtab fix). Evaluation-only rerun 22724842
  (`gate/rerun_eval.sbatch`, frozen copy in the gate dir); scoring shard 22724843 afterok on it. Shard 22721980 had
  started afterany on the gate, found no gate table and exited (NO_RUN_READY); nothing scored.
- 2026-11-05 ~13:50: COLO320 scoring shard 22724843 stopped in its first step ("object 'loss_wgs' not found"): score_one.R
  passed no params for COLO320, but their rule evaluate_cnv_prediction_cnvkit_COLO320 sets genomic_format = "CNVkit"
  (workflow.sm at 319ea34). No COLO320 score was produced. Fixed in score_one.R (backup dev/score_one_pre_cnvkitfix_20261105.R);
  the change touches only the COLO320 branch, so the five running shards (BCC06, BCC06post, MM x2, MCF7) get identical
  params. Fixture v6 does not exercise COLO320, so a COLO320 fixture (22725812, score/test/run_fixture_colo320.sbatch,
  copyKat stand-in, HaploDose == copyKat, CaSpER flagged) runs on the same frozen code as the new shard 22725813.
  COLO320 scores are read only after both fixture v6 and the COLO320 fixture pass.

### Claim 5 ranking script frozen, 2026-11-05 ~20:20 (before any HaploDose score file was opened)

- `benchmarks/schmid/score/rank_claim5.R` sha256 87a3f5c4f3d81f6c5f006aa150824157c24ec93398550001b2956cffebd9e5e2. Per dataset and metric (Pearson vs truth from the corr table,
  AUC gains, AUC losses), 10 outputs ranked (1 = best), ties average, missing output last (average among missing).
  SOTA = HaploDose strictly lowest mean rank over 3 metrics x 4 primary datasets (BCC06, BCC06post, COLO320, MM_ctrl8);
  a tie for first is not SOTA. Descriptive: normalised overall score (Pearson floored 0, (AUC-0.5)x2 floored 0, max F1).
  Sets reported: primary (combined scorer, vnext arm); sensitivities: MM_ctrl568 in place of MM_ctrl8, per-method
  scorer, sens1424 arm; MCF7 line-overlap separately.
- Tested on fixture tables (copyKat stand-in; no HaploDose output): `dev/schmid_claim5_20261103/rank_test_20261105_2010/`.
  HaploDose ties copyKat exactly; a removed Numbat output ranks 10th.
- 2026-11-05 ~22:00: COLO320 fixture 22725812 FIXTURE_PASS (58 checks; HaploDose stand-in == copyKat in every flagged table, CaSpER flagged, CNVkit truth used). dev/schmid_claim5_20261103/score_fixture_colo320_20261105_1303/check_fixture.log
- 2026-11-06: fixture v6 22690386 FIXTURE_PASS (236 checks; case exits A 0, B 3 as expected, C 0). Scores may now be read.

### Claim 5 result (Schmid et al. 2025 benchmark, read 2026-11-06 ~00:10, ranking job 22737628)

Prespecified outcome: **claim 5 fails.** HaploDose has the worst mean rank of the 10 outputs on the primary set.
Source: `results/schmid_claim5_rank_20261106_0003/{summary,verdict,per_dataset}_{combined,permethod}.tsv`
(frozen `rank_claim5.R` sha256 87a3f5c4…; scores from `results/schmid_score_claim5_20261105_0836/`, after fixture v6 and
the COLO320 fixture passed).

- Unit = dataset; n = 4 primary (BCC06, BCC06post, COLO320, MM199 control {8}); 3 metrics (Pearson vs DNA truth, AUC
  gains, AUC losses); 10 outputs. Descriptive rule, no test (n = 4), as prespecified.
- Primary (combined scorer, v-next arm) mean rank: inferCNV (Expr) 3.33, Numbat (Expr) 4.04, copyKat 4.12, Numbat (CNV)
  4.25, inferCNV (CNV) 5.25, CaSpER 5.42, CONICSmat 6.25, SCEVAN (CNV) 6.83, SCEVAN (Expr) 7.08, **HaploDose 8.42**.
  Normalised overall score (descriptive): HaploDose 0.246, lowest; comparators 0.278 to 0.364.
- Per dataset HaploDose Pearson / AUC gains / AUC losses: BCC06 0.164 / 0.582 / 0.583; BCC06post 0.138 / 0.598 / 0.540
  (last on all three); COLO320 0.187 / 0.652 / 0.490; MM199 0.282 / 0.645 / 0.652. Pearson is about half of copyKat's
  on BCC06post, COLO320 and MM199.
- Sensitivities, same conclusion: MM199 control {5,6,8} 8.50 (last); 14:24 arm 8.42; per-method scorer 8.00
  (9th of 10; best inferCNV (Expr) 3.62). MCF7 line-overlap (separate): HaploDose 7.67, 8th of 10 (Numbat outputs
  missing, ranked last).
- Checks before accepting: HaploDose gene coordinates match the hg38 annotation used by the scorer (EGFR, MYC, GAPDH,
  TP53 identical); Pearson and AUC are scale-free, so the +2 offset in read_haplodose cannot change ranks; the two
  HaploDose arms differ in values but not in ranks. Not investigated for tuning: Schmid stays an evaluation set.
- Interpretation: on external single-cell data with DNA truth, HaploDose's clone-level dosage profile tracks DNA copy
  number less well than every comparator output. This agrees with claim 1 (MIX-seq profile below a pseudobulk baseline).
- Claims status: hold 3 (CCLE bulk), 6 (CPTAC-3 bulk); fail 1, 2, 4, 5. The prespecified SOTA goal (>= 3 wins) is not met.

## Claim 7 prespecified 2026-11-06 (user "do all"), before any REBC-THYR data was downloaded

Survey from GDC metadata only: `dev/claim7_survey_20261106/survey.txt`. Every TCGA cohort with >= 10 solid normals
was used in development (`benchmarks/tcga_bulk`), and CPTAC-3 is claim 6. The only other GDC project with open
AscatNGS, primary-tumour STAR - Counts and >= 10 Solid Tissue Normal RNA files is REBC-THYR (radiation-related
papillary thyroid cancer; 449 AscatNGS WGS segment files, open; 437 primary tumour and 400 solid normal RNA files).

7. **Focal (sub-arm) loss detection in tumour bulk RNA, REBC-THYR.** Third replication of the claim 3/6 endpoint on
   tumours not used in development.
   - Everything is as claim 6 (code, scorer, truth, RNA file rule, normal-panel rule, pooled sensitivity, test,
     Holm family, bootstrap seed), except: project REBC-THYR; one group THYR (primary_site thyroid gland, disease
     type adenomas and adenocarcinomas).
   - Unit: tumour (case) with >= 1 focal-loss window. Primary: focal-loss recall at 5% realised FPR, HaploDose dose vs
     naive and vs infercnvpy, paired Wilcoxon signed-rank, Holm over these 2. SOTA criterion: both Holm p < 0.05 and
     mean paired difference > 0.
   - Power floor, stated in advance: papillary thyroid cancers carry few copy-number changes, so many tumours may
     have no focal-loss window. If fewer than 20 tumours qualify, claim 7 is reported as not evaluable, which counts
     as not won. The number of qualifying tumours is not looked at before the run.
   - Code: `benchmarks/rebc_bulk/` wraps `benchmarks/cptac_bulk/` (project and cohort only). Eval code hashes are
     written before the fetch.
   - Caveats in advance: as claim 6 (stroma, purity-corrected truth vs uncorrected RNA, AscatNGS deprecated). Most
     normals are matched to a tumour, as in CPTAC. A thyroid win would replicate claims 3 and 6 on a third bulk
     cohort; it is not a new modality.
   - Run 2026-11-06: eval code hashed in `results/rebc_bulk_20261106_0649/eval_code_sha256.txt` (claim 6 hashes re-checked OK) before any
     download. Manifest `data/external/rebc_thyr_20261106/files.tsv`: 442 AscatNGS and 837 RNA files after the THYR
     map (survey count 449 included cases outside that map), 3.6 GB, all open. Jobs: fetch 22759847 (io), then
     input 22759934, group 22759935, score 22759936 (afterok chain, cpu v6 preferred, 8 CPU 48G each).

## Profile-gap diagnosis, 2026-11-06 (development data only: Kinker; Schmid not reopened)

`dev/profile_gap_20261106_0007/` (job 22737862; `out/summary.json`, `out/paired_tests.tsv`). Question: why does
HaploDose rank last on Schmid's bin-level Pearson when its arm-level rho on Kinker beats CopyKAT? n = 181 Kinker
single-line specimens (unit = cell line); truth = DepMap relative CN as in `benchmarks/kinker/evaluate.py`.
Paired Wilcoxon vs `hd_lrg`, Holm within each metric over 6 pairs.
- Arm rho (median): hd_lrg 0.861, pseudobulk 0.863, CopyKAT 0.812. No gap at arm level.
- 1 Mb Pearson (Schmid-like): hd_lrg 0.619 (the output exported to Schmid, log_ratio_gene), hd_dose 0.659,
  CopyKAT 0.672, hd_lrg after a 5-bin running median 0.696 (higher than hd_lrg in 152/181 lines, Holm p < 1e-20).
- Roughness along the genome: hd_lrg 0.369, truth 0.290, CopyKAT 0.052, hd_lrg_sm5 0.091.
- Reading: the bin-level gap comes mostly from bin-to-bin noise in the exported output, not from a wrong
  arm-level signal. A smoothed export would likely close it on Kinker [Inference for Schmid; not tested there,
  and Schmid stays closed for tuning]. Any smoothed export is a new version that must be frozen and evaluated only
  on unopened data.
- Agent report details (same run): smoothed HaploDose is above CopyKAT in 114/181 lines (median +0.017); arm rho
  stays 0.863. CopyKAT here is Schmid's input (mean of CNA_raw_results_gene_by_cell), hence 0.812, not 0.765. Winsorising
  does not change the gap (0.614), so outliers are not the cause. The holdout on/off comparison was not run. The 5-bin window was picked
  after seeing these results, so a smoothed export is a new version. It would be frozen on Kinker, checked on
  Andor, and evaluated only on unopened data.

## Claim 8 (Bian scTrio) survey, 2026-11-06: not a confirmatory claim

`dev/bian_claim8_survey_20261106_0654/survey.md` (no RNA tool run).
- There are 609 same-cell trio cells in 6 patients. DNA CN exists for only 140 of them (CRC01 and CRC11), and those cells were selected with RNA markers. The existing clone labels use RNA-informed features. At the existing cut, each patient has one dominant DNA clone.
- At most 5 evaluable patients. The smallest exact Wilcoxon p at n = 5 is 0.0625, so no per-patient test can pass Holm. A cell-level test would be pseudoreplication.
- GEO RNA is FPKM/TPM, so it would need STAR + featureCounts. Only 20 normal cells are available as a reference.
- Blinding: the survey agent read one TPM line and two cells' marker values by accident. The CRC project had already used Bian RNA.

Decision: Bian is not a SOTA claim. At most it is a descriptive cell-level check, after RNA-blind DNA CN is re-derived for all trio cells.

## Claim 7 amendment 1, 2026-11-06, before any claim 7 score was computed

The first run (`results/rebc_bulk_20261106_0649/`) matched only 5 of 383 ASCAT cases to tumour RNA. REBC RNA files
carry sequencing-centre sample IDs (`SC*_merged`), never the ASCAT tumour sample ID, so the claim 6 rule (RNA file
whose sample list contains the ASCAT tumour sample) fails on metadata alone. Its score stage crashed (subarm.py
groupby on an empty curve table, 5 tumours). No score, recall or per-tumour output of that run was read. The 5-tumour
HaploDose and infercnvpy outputs there are left unread.
Amendment: cases with no sample-ID match fall back to a case-level match: the case's Primary Tumor STAR - Counts file;
several -> smallest file_id (8 cases have 2). `selection.tsv` gains column `match` (sample | case).
Dry run on metadata: 372 tumours (5 sample-matched, 367 case-matched), 400 normals, 0 pooled. Caveat: for the case-matched tumours,
the RNA and DNA may come from different aliquots of the same tumour (one primary tumour sample per ASCAT case).
Everything else is unchanged (test, Holm over 2, power floor 20). Rerun: `results/rebc_bulk_20261106_0711/` (eval code re-hashed; builder
before amendment kept at `dev/rebc_make_input_pre_amend1_20261106.py`). Jobs: input 22760991, group 22760992, score 22760993.

## Claim 7 result (read once, 2026-11-06, `results/rebc_bulk_20261106_0711/score/primary.json`): not evaluable, so not won

Run directory is `results/rebc_bulk_20261106_0711/` (the amendment text above names it via the same variable). 372 tumours, 400 normals.
Truth: median ploidy 2.00 (IQR 1.99-2.00). Only 6 tumours have >= 1 focal-loss window (10 have a focal gain), below
the prespecified floor of 20. Descriptive only (n = 6, unit = tumour; Wilcoxon p not computable, < 5 non-zero
differences; Holm over 2 = 1): focal-loss recall at 5% FPR, mean HaploDose 0.50, naive 0.58, infercnvpy 0.17; 4 of 6
tumours tie on both comparisons. The power risk stated in advance came true: papillary thyroid cancers are nearly
copy-number neutral.
Claims status: hold 3 (CCLE bulk), 6 (CPTAC-3 bulk); fail 1, 2, 4, 5; 7 not evaluable. The prespecified SOTA goal
(>= 3 wins) is not met. No further unopened GDC cohort meets the claim 6 rules (survey
`dev/claim7_survey_20261106/survey.txt`), and Bian cannot support a confirmatory claim.

## Claim 8 data survey, 2026-11-06 (metadata only; no claim prespecified)

`dev/claim8_data_survey_20261106_0727/survey.md`.
- Candidates for a profile claim, with n = cell lines that have DepMap truth and no development overlap:
  - Gambardella GSE173634 (Drop-seq): 17 lines.
  - Tahoe-100M DMSO wells (Parse chemistry, new to HaploDose): up to 51 lines [Unverified: which of these were sequenced].
  - CPTAC-3 scRNA: 30 cases, but 28 have DNA truth already exposed in claim 6.
- Power from development: on Kinker, smoothed HaploDose beat CopyKAT in 114/181 lines, a median of +0.017 [Inference]. At that rate, 17 lines are unlikely to pass Holm. About 50 lines might.
- DepMap DNA truth for these lines was likely already used in claim 3 [Inference]. The RNA data are unopened.
- Nothing is prespecified. It awaits the user's decision.
- Tahoe-100M metadata (`dev/tahoe_meta_20261106/`):
  - The cell-line table lists 102 names (99 with DepMap IDs); the obs table, not yet read, says which were sequenced.
  - There are 14 plates and 31 DMSO samples, at about 2,000 transcripts per cell.
  - Lines are pooled in each well, which is the MIX-seq design where claim 1 failed against pseudobulk. A fair claim needs a smoothed pseudobulk comparator too [Inference]. Odds of a win look low.

## Claim 8 prespecified 2026-11-06 (user: "Tahoe claim 8"), before any Tahoe expression value was read

Read so far: Tahoe-100M `cell_line_metadata.parquet` and `sample_metadata.parquet`, plus the obs_metadata schema
(`dev/tahoe_meta_20261106/`). No expression row has been read.

Multiplicity, stated up front: this is the eighth prespecified claim. It was chosen after claims 1-7 were known,
and after a development diagnosis (profile gap) that motivated its endpoint. All eight are reported.

8. **Relative RNA dosage profile at 1 Mb, pooled cell lines, Tahoe-100M (Parse chemistry, new to HaploDose).**
   - Cells: DMSO_TF control samples of plate1 (the 2 samples with drug == 'DMSO_TF'; trametinib-solvate samples
     excluded). Only cells with pass_filter 'full' are used [if that value does not exist, the authors' filter flag is used as is,
     recorded before any expression is read]. Line = obs cell_name (Tahoe SNP demultiplexing).
   - Lines: >= 50 cells. At most 300 cells per line, sampled with seed 20261106. DepMap 24Q4 OmicsCNGene available. The
     design is benchmarks/heldout/prepare_mixseq.py (benchmarks/kinker/prepare.py write_design): one pool, two
     rotations, half the eligible lines as reference, 70% fit / 30% holdout, seed 20261106.
   - Primary population: query lines not used in development or earlier claims as scRNA. Excluded: the Kinker,
     MIX-seq, Andor and Schmid (MCF7, HCT116, A375, COLO320) lines. They can still serve as reference lines. DepMap DNA
     for most lines was used as truth in claim 3. The RNA is new.
   - Truth per gene: DepMap CN of the line minus the UMI-weighted mean of the pool's fit reference lines
     (benchmarks/kinker/evaluate.py). It is averaged into 1 Mb bins. Only bins shared by all methods for that line are kept.
   - Methods (all from one run of each on the same input):
     - HaploDose frozen 14:24 code (`results/heldout_20261102_1424/code`, claim 1 flags) with a fixed export:
       log_ratio_gene, a running median over 5 HaploDose bins within each chromosome, clones weighted by assigned
       cells. The window was picked on Kinker after seeing the profile-gap result.
     - CopyKAT: mean over the line's query cells of the gene-by-cell CNA_raw_results.
     - infercnvpy: mean over the line's query cells.
     - Smoothed pseudobulk: the evaluate.py pseudobulk, as medians over HaploDose's 10-gene bins, then the same
       5-bin running median.
   - Endpoint per query line: Pearson between the method's 1 Mb values and 1 Mb truth. Primary test: HaploDose vs
     each of the 3 comparators, paired Wilcoxon over lines (unit = cell line), Holm over 3. SOTA criterion: all 3 Holm
     p < 0.05 and mean paired difference > 0, with bootstrap 95% CI of the mean difference (10,000 resamples of lines,
     seed 20261106). A comparator that fails on a line scores NaN; that line is dropped from that pair only.
   - Power floor: fewer than 20 primary lines means not evaluable, which counts as not won.
   - Scorer: written from `dev/profile_gap_20261106_0007/profile_gap.py`, hashed, and checked on Kinker first. On
     Kinker it must reproduce that run's hd_lrg_sm5, copykat and pb_bin10 medians to 3 decimals. It is frozen
     before any Tahoe method output exists.
   - Caveats stated in advance:
     - Parse chemistry is new to every method here.
     - All 50 lines share one well, so the reference is other cancer lines, as in MIX-seq, where claim 1 failed
       against pseudobulk.
     - Development evidence predicts small margins (Kinker median +0.017 vs CopyKAT).
- Scorer fixture on Kinker passed, 2026-11-06 (`dev/claim8_scorer_fixture_20261106_1058/`, job 22774343).
  - The 1 Mb medians match the profile-gap run to 4 decimals: hd_sm5 0.6964, copykat 0.6725, pb_bin10 0.6353, hd_lrg 0.619.
  - The infercnvpy reader works (median 0.654). Scorer frozen (read-only): `eaa5c8eaadded527db125149a2dd10525c60f0222ce7529de388dd1213708f6e  benchmarks/tahoe/score_profile.py`.
- Development prediction, from the same fixture (n = 181 Kinker lines; Holm over 3):
  - The smoothed pseudobulk baseline beats smoothed HaploDose: median 0.708 vs 0.696, mean diff -0.014, HaploDose
    higher in only 43/181 lines, Holm p 5e-18.
  - HaploDose beats CopyKAT (+0.016, Holm p 2e-6) and infercnvpy (+0.028, Holm p 1e-6).
  - Under the prespecified rule (all 3 must pass), claim 8 would fail on development data. The comparator set is
    not changed: dropping a baseline because it wins would bias the claim.

## Claim 8 withdrawn 2026-11-06 (user decision), before any Tahoe expression value was downloaded or read

Reason: on development data (Kinker, n = 181 lines), the frozen scorer predicts a fail. Smoothed pseudobulk beats
smoothed HaploDose (median 0.708 vs 0.696, Holm p 5e-18), and the rule needs HaploDose to beat all 3 comparators.
Only metadata was fetched (`data/external/tahoe100m_20261106/metadata`, checksummed). The fetch job was cancelled before the
shard download. Claim 8 is reported as prespecified and then withdrawn, not as a test.
Final claims status: hold 3 (CCLE bulk) and 6 (CPTAC-3 bulk); fail 1, 2, 4 and 5; 7 not evaluable; 8 withdrawn. The
prespecified SOTA goal (>= 3 wins) is not met.
Development finding for the paper: at 1 Mb resolution, a smoothed pseudobulk is a strong baseline that HaploDose's
profile output does not beat (Kinker). Arm level: HaploDose 0.861 vs pseudobulk 0.863.
- Tahoe fetch stopped after the footer scan: 76 of 3388 shards would have been needed (~5.4 GB). A partial 8 MB shard download (never read) was deleted.
