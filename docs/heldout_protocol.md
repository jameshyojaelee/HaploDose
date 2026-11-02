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
  computed by `benchmarks/kinker/evaluate.py`. Numbat is added if its Kinker pipeline is complete at the
  freeze point (single-line design only).

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

## Held-out input (built 2026-11-02, no method run)

`results/heldout_mixseq_input_20261102_1159/` (`benchmarks/heldout/prepare_mixseq.py`, job 22607080): DMSO_24hr_expt1 2,623, DMSO_6hr_expt3 8,398, DMSO_expt10 4,452 QC-passing cells; 32,643 gene symbols (18,097 with a GRCh38 position); median 13,549 UMIs per cell. 127 lines meet the rules (each a query once over the two rotations); 44 were also in Kinker, so the sensitivity population has 83 lines. Per-line table: `lines.tsv`.

## Freeze record

(to be filled before the first held-out run: date, package SHA-256 list path, CLI defaults, comparator env
versions)
