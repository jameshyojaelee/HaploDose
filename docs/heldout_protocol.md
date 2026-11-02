# Held-out validation protocol (written 2026-11-02, before any method was run on these data)

## Why

Kinker 2020, Andor 2020 and U-CAN were used to choose model changes, so their results are development
results and are optimistic. Publication claims about Haplodose versus other tools rest on this held-out
evaluation, run once on frozen code, with the endpoints and tests below fixed in advance. Any deviation from
this protocol is reported as a deviation, with its reason.

## Frozen method

- Haplodose: the repository package at the freeze point, copied to `results/heldout_<date>/code/` with
  SHA-256 of every file in `code_sha256.txt`, default parameters (CLI defaults at the freeze point; no
  dataset-specific options). The freeze point is recorded in this file (section "Freeze record") before the
  first held-out run. Development may continue afterwards, but held-out numbers are only ever reported for the
  frozen version.
- Comparators, as run on Kinker (`benchmarks/comparators/`): infercnvpy (`run_infercnvpy.py`), R inferCNV
  (`run_infercnv_r.py`), CopyKAT and SCEVAN (`run_r_field.py`, env `cnv_field`), all with their wrapper
  defaults; `--rule cluster` for the single-line design, `--rule all` for the implant design. Pseudobulk is
  computed by `benchmarks/kinker/evaluate.py`. ~~Numbat is added if its Kinker pipeline is complete at the
  freeze point (single-line design only).~~ (struck: amendment 2)

## Data

Primary: MIX-seq (McFarland et al. 2020, Nat Commun 11:4296; Figshare 10298696 v3, CC-BY 4.0), control pools
only (DMSO and untreated; files and md5 in `benchmarks/heldout/mixseq_files.tsv`). Selection rules, fixed now:
- One control experiment per distinct pool composition, so a line enters once per pool: the DMSO 24 h sample
  where a pool has one, otherwise its other control; ties broken by more QC-passing cells.
- Cells: `cell_quality == normal` (the authors' QC-passing singlets).
- Lines: DepMap 24Q4 copy number available and at least 50 QC-passing cells in the chosen experiment. A line in
  several pools is a query once: in the pool with the most cells of that line.
- Primary population: all lines meeting these rules. Prespecified sensitivity population: lines not in the
  Kinker 2020 pools (never seen in development in any form).

Secondary: Gambardella et al. 2022 (Nat Commun; GEO GSE173634, to verify), breast lines, used only if its counts
and line labels can be read without changing any method; same rules and endpoints, reported separately.

## Design

As `benchmarks/kinker/prepare.py`: within each pool, eligible lines are split in half at random (seed
20261103); in rotation r, half r is the reference (one reference group per pool, 70% fit / 30% holdout cells)
and each line of the other half is its own single-line specimen; every line is a query exactly once. Truth per
gene = line DepMap CN minus the UMI-weighted mean of the pool's fit reference lines; arm truth = median over the
run's used genes. Implant design: `benchmarks/kinker/implant.py` design v3 (depth-neutral gains), its settings
unchanged (50/150/500 genes, gain/loss, 15%/40% of cells), seed as in that script.

## Endpoints and tests

E1 (primary): per-line Spearman rho between the tool's arm profile and DepMap arm truth (evaluate.py; arms with
>= 5 genes). Haplodose's endpoint column is `log_ratio_gene` (fixed now; `log_ratio` reported as secondary).
Comparisons: Haplodose vs each of infercnvpy, R inferCNV, CopyKAT, SCEVAN and pseudobulk; paired two-sided
Wilcoxon signed-rank over lines; Holm over these 5 comparisons. Effect: median paired difference with a 95%
percentile bootstrap CI (10,000 line resamples, seed 20261103). Claim of superiority over a tool only if its
Holm-adjusted p < 0.05 and the CI excludes 0 in Haplodose's favour.

E2 (secondary): implant design, per line ARI between tool clones and implanted cells (`score_implant.py`
`ari_all`); Haplodose vs each of infercnvpy, R inferCNV, CopyKAT, SCEVAN; paired Wilcoxon, Holm over 4; same
bootstrap CI. Reported overall and by floor stratum (implanted cells < 30 vs >= 30) and by setting
(descriptive).

E3 (secondary, descriptive): sign agreement on DNA-informative arms (|truth| >= 0.2); number of tumour clones
per single-line specimen (one expected); coverage (cells with output / cells sent).

Unit of replication: the cell line. Lines share a pool and its reference, so pool is reported as a cluster;
a sensitivity analysis repeats E1 with pool-level medians (Wilcoxon over pools) and is descriptive.

Failures: a tool that fails or returns no output for a line scores rho = NaN for E1, and that line is dropped
from that pair's test only (count reported); for E2, `ari_all` scores a missing specimen 0.

## Amendments (each dated, made before any method was run on the held-out data)

1. 2026-11-02, before any method run, after reading only the authors' cell metadata (classifications.csv:
   line and QC per cell). The control pools are three compositions: experiment 1 (DMSO 24 h, DMSO 6 h,
   untreated 6 h; 24 lines), experiment 3 (DMSO 24 h, DMSO 6 h; 97 lines) and experiment 10 (DMSO; 93 lines).
   The rule "DMSO 24 h first" would use experiment 3's 24 h sample with 53 lines at >= 50 cells, while its 6 h
   sample has 83. Amended rule: per composition, the control sample with the most lines at >= 50 QC-passing
   cells, ties broken by more QC-passing cells. Selected: DMSO_24hr_expt1 (23 lines; ties untreated 6 h at
   23, 2,623 vs 2,417 cells), DMSO_6hr_expt3 (83), DMSO_expt10 (34). Reason: sample size only; no outcome seen.
2. 2026-11-02, before any method run. Numbat is not included in the held-out evaluation: it needs per-line
   allele counts from the held-out reads (pileup and phasing on the matching genome build), and only count
   matrices are prepared for MIX-seq. Its Kinker single-line result is reported as development only. (Replaces
   the conditional sentence under "Frozen method".)
3. 2026-11-02, before any method run. (a) The E1 and E2 bootstrap CIs resample lines and ignore that lines share
   a pool and its reference (3 compositions x 2 rotations = 6 reference sets, too few for a cluster bootstrap);
   this is stated with every CI, and the pool-level sensitivity analysis stays descriptive. (b) Next to each E1
   comparison, the number of lines where either tool returned no profile is reported (`compare_tools.py`
   n_missing), with the tool that failed. (c) The freeze record also lists SHA-256 of the evaluation code:
   `benchmarks/kinker/{evaluate,compare_tools,score_implant,implant}.py`, `benchmarks/comparators/*`, and
   `benchmarks/heldout/prepare_mixseq.py`.
4. 2026-11-02 (~15:30), before any method run: the 8 comparator jobs submitted at the 14:24 freeze
   (22611772-75, 22611777-80) and the analysis job 22611794 were cancelled while pending, with no held-out output
   written, because a review of the job scripts found faults that would have stopped or corrupted the run. The two
   Haplodose jobs (22611771, 22611776) were kept: their path through the 14:24 script is not affected by (a)-(c). Technical fixes only; endpoints, tests,
   families, data and the frozen Haplodose code are unchanged.
   (a) `run_tool.sbatch` wrote each comparator's stderr to `<tool>_time.txt` inside the pool directory, and the
   run glob `<tool>_*` matched that file, so `evaluate.py` would have failed for all four comparators and the
   `afterok` analysis would never have run. stderr now goes to `time_<tool>.txt`; runs are collected as
   `<tool>_2*/` directories.
   (b) Failures now follow the "Failures" rule above instead of stopping the pipeline. Only runs with `qc.json`
   (the last file every wrapper and Haplodose write) are scored; incomplete runs, failed pools and failed
   evaluations are listed in `<tool>_failures.txt`; tool jobs exit 0 after recording a failure and the analysis is
   submitted with `afterany`. In `analyse.py`, a tool with no output keeps its place in each Holm family (5 for E1,
   4 for E2): E1 lines it did not return are NaN and dropped from that pair (counted); a tool with no E2 output
   scores `ari_all` 0 on every line; a comparison with no testable pairs (or all differences 0) enters Holm with
   p = 1 and is flagged `testable: false`.
   (c) The comparator jobs of one design shared per-pool cell files written with an unlocked check-then-write;
   they are now written to a job-private temporary file and renamed.
   (d) Each resubmitted job checks the frozen Haplodose package against `code_sha256.txt` at start (the two kept
   Haplodose jobs do not; all 41 hashes were verified at ~15:30 and `code/` was then made read-only); the skip guard for an
   existing Haplodose run tests `qc.json` (the pipeline writes no `run.json`). `code/` and `eval_code/` are made
   read-only.
   (e) `analyse.py` now also computes what this protocol already lists but the 14:24 version left out: the
   pool-level Wilcoxon (descriptive, nominal p; with 6 reference sets the smallest two-sided p is 0.03125), E3
   coverage (query cells with a clone label / query cells sent; a line with no output scores 0) and tumour clones
   per specimen, and E2 by setting.
   (f) Clarifications, no change to the code that defines them: E2 covers the query lines with at least 100 cells
   (`implant.py` MIN_CELLS, as on Kinker); the implant log of 2026-11-02 lists 47 implanted specimens. The
   per-pool line splits use seed 20261103 + k in pool order (`prepare_mixseq.py`), not 20261103 itself. The
   conditional Numbat sentence under "Frozen method" is struck (amendment 2).
   The 14:24 `analyse.py` and `analyse.sbatch` stay in `eval_code/` and `analysis_sha256.txt`; the amended versions
   are `eval_code/analyse_v2.py`, `analyse_v2.sbatch` and `analysis_sha256_v2.txt`. Before resubmission the fixed
   scripts were run end to end on Kinker pools 16 and 18 (development data) with an incomplete run directory and
   three tools without output (`dev/heldout_fixcheck_20261102_1528/`).
5. 2026-11-02 ~16:30, before any held-out job had started (sacct: 22611771, 22611776 and all comparator jobs
   pending; no held-out output exists). Prompted by an external code review
   (`docs/scientific_code_review_2026-11-02.md`, findings 4 and 5; verdicts in `..._response.md`). Primary
   endpoints, tests and frozen scripts are unchanged. Added, all fixed now and computed by a separate script
   (`benchmarks/heldout/sensitivity.py`, hashed into `eval_code/` with `analyse_v2.py` before resubmission):
   (a) Targets differ between tools: arm truth is the median over each run's used genes, and HaploDose uses an
   expression-filtered gene set while comparators use all positioned genes. Sensitivity: E1 rescored for every
   tool on a common gene universe per pool and rotation (genes used by every tool that produced output there),
   arms eligible at >= 5 common genes, same truth formula, same tests and Holm family of 5.
   (b) E3 coverage restricted to query cells present in each tool's own matrix; n_in_matrix / n_sent per tool.
   (c) E2 repeated with cells absent from a tool's matrix left unlabelled instead of clustered from zeros.
   (d) Lines share 6 reference sets (3 compositions x 2 rotations), so line-level tests assume an independence
   that does not hold. `superior` is reported as "conditional on independence of lines". Descriptive: a cluster
   sign-flip permutation over the 6 reference sets of per-set median differences (64 patterns, p >= 1/32) and a
   linear mixed model of the per-line difference with a random intercept per reference set.
   Decision rule tightened (can only make the claim harder): prespecified SOTA claim 1 (`docs/STATUS.md`) holds
   only if the primary E1 rule is met AND the common-universe sensitivity (a) gives a mean difference > 0 with
   bootstrap 95% CI above 0 against each of the four published tools.

## Held-out input (built 2026-11-02, no method run)

`results/heldout_mixseq_input_20261102_1159/` (`benchmarks/heldout/prepare_mixseq.py`, job 22607080): DMSO_24hr_expt1 2,623, DMSO_6hr_expt3 8,398, DMSO_expt10 4,452 QC-passing cells; 32,643 gene symbols (18,097 with a GRCh38 position); median 13,549 UMIs per cell. 127 lines meet the rules (each a query once over the two rotations); 44 were also in Kinker, so the sensitivity population has 83 lines. Per-line table: `lines.tsv`.

## Freeze record

Frozen 2026-11-02 14:24, before any held-out run. Root: `results/heldout_20261102_1424/`.
- Haplodose package: `code/haplodose/` (repo python/haplodose at 14:24: clones.py with review fixes #1/#2/#7,
  estep #9, gene shares #3, assignment_dose column); clones.py SHA-256 974f5573..., full list in `code_sha256.txt`
  (41 files, package and evaluation code). CLI options: `cli_run_help.txt`; defaults from the frozen source: min_cells 30 (cli.py:30), event windows none (cli.py, CloneParams.event_windows = ()).
- Evaluation and comparator code: copies in `eval_code/` (evaluate, compare_tools, score_implant, implant,
  prepare_mixseq, comparators/*); every held-out job checks the repo files it runs against these hashes and stops
  on a mismatch (`benchmarks/heldout/run_tool.sbatch`).
- Environments: `pip_freeze_haplodose.txt`, `pip_freeze_scrna_cnv.txt` (infercnvpy and wrappers),
  `sessionInfo_cnv_field.txt` (CopyKAT, SCEVAN); R inferCNV from module R/4.4.1 as on Kinker.
- Hardware: every held-out run on cpu `--constraint=v6`, 16 CPUs.
