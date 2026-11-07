# Response to the scientific and code review of 2026-11-02

Response to `docs/scientific_code_review_2026-11-02.md` (not edited). Written 2026-11-02 by a review subagent.
Each finding was checked against the working tree (line numbers below are current, not the review's) and,
where cheap, reproduced on a compute node. Changes to `python/haplodose` described here are **v-next**. They are
not the frozen held-out version in `results/heldout_20261102_1424/code`, which is what gets reported on MIX-seq.
No held-out output, MIX-seq input contents or `results/ccle_bulk_*` outputs were read.

Reproduction: `dev/review_response_20261102/repro.py`. It was run on a snapshot of the package before these
fixes (`dev/review_response_20261102/before/`; job 22612795 failed on an output path, rerun as job 22612820)
and again after the fixes (first block of job 22612872's log).

| # | Verdict | Action |
|---|---|---|
| 1 | CONFIRMED | Fixed in `profile.py`. A clones.py logging diff is proposed. |
| 2 | CONFIRMED (both parts) | Fixed in `alleles.py`. A clones.py diff is proposed to pass `log_dose` to the veto. |
| 3 | CONFIRMED | No code change. A thin-mode diff is proposed; null calibration is needed in all cases. |
| 4 | CONFIRMED (both parts), size unknown | Frozen code, so not edited. Proposed as (b), a sensitivity analysis declared now. |
| 5 | PARTLY | Frozen code, so not edited. Proposed as (b), a sensitivity analysis plus wording. |
| 6 | CONFIRMED | Fixed in `alleles.py`. |
| 7 | CONFIRMED for cnLOH; PARTLY for bulk | Label changed in `alleles.py`. A note was added in `bulk.py` (no column or key renamed). |
| 8 | Already fixed | Nothing remains beyond rerunning the tests (see below). |

## Findings

### 1. No depth overlap gives NaN difference segments: CONFIRMED

Evidence. `matched_differences` (`clones.py:631`) gives both groups zero weight when no depth stratum holds
cells from both groups (`m = min(na, nb) = 0`). `difference_segments` then computes `log(0.5 / 0)`. The HMM
runs on NaN emissions, and `support < params.posterior` is False for NaN, so the run is kept. Before the fix,
two groups at 800 vs 20,000 UMIs gave 2 segments with `mean_diff`, `z` and `posterior` all NaN (job 22612820).

Consequences. In `merge_unsupported` the NaN segment keeps the pair unmerged. In `try_splits`, `min(scores)`
with a NaN score depends on its position in the list. NaN first gives NaN and rejects the split. NaN later is
skipped, so the split is accepted on the other pairs.

Fix. `profile.py:182` and `profile.py:200`: a pair with a zero UMI total or non-finite log ratios is
`testable=False`, with no segments and distance inf. Callers already treat this as "no difference": merge
last, reject the split. Every result now also carries `testable`. After the fix the same input gives 0
segments and distance inf, and the overlapping pair is unchanged (job 22612872).

Tests:
- `tests/test_scientific_review_20261102.py::test_zero_umi_pair_is_untestable_not_nan`
- `::test_matched_differences_disjoint_depths_untestable`

Not done. There is no minimum-overlap threshold beyond zero. Small overlap still gives a valid but
low-powered test.

### 2. Allele logic discards real differences: CONFIRMED

Veto (`alleles.py:416`). A one-copy change from a state with both haplotypes present (m:M, m >= 1) always
changes |b - 0.5|. So the veto's premise holds there, and the 1:1 -> 2:2 exception in its docstring was never
implemented. Two cases preserve the allele ratio and are invisible to alleles:
- Any change inside LOH: 1:0 vs 2:0 both read b = 1. This is the realistic case (loss vs cnLOH subclones).
- Proportional changes, such as 1:1 vs 2:2.

Before the fix, both identical-balanced and identical-LOH counts vetoed (job 22612820).

Fix:
- No veto when both groups' estimated |b - 0.5| >= 0.375 (LOH-like).
- New optional `log_dose`. No veto when |log_dose| >= log 1.75, which is nearer a doubling than a one-copy gain.
  A one-copy change of that size (1:1 vs 1:0) changes imbalance, so the veto would not fire on it anyway.

A one-copy-sized difference (log 1.5) with unchanged imbalance still vetoes. The 0.375 and log 1.75 cut
points are design choices, not calibrated [Inference]. This is a compatibility check on two summaries. A
joint allele and copy-state model, as the review suggests, is a larger v-next design item.

`differ` (`alleles.py:377`). It compared (b - 0.5)^2, so 0.9 vs 0.1 scored 0.0 (job 22612820). A signed
per-SNP contrast `_delta2` (`alleles.py:279`) needs no phase either, because both groups share the
orientation at each SNP. It is now tried when the magnitude test does not resolve a run. After the fix, 0.9 vs
0.1 scores 26.9 (job 22612872). In a null check (300 replicates, beta-binomial, shared phase flips) it may
add at most 3 false differences. That bound is the test's threshold, not a measured rate.

Tests:
- `::test_veto_keeps_one_copy_logic_but_exempts_ratio_preserving_changes`
- `::test_differ_sees_mirrored_imbalance`
- `::test_differ_null_false_positive_rate`

`log_dose` takes effect only once `clones.py` passes it. See the proposed diff.

### 3. Split-validation independence does not survive EM: CONFIRMED

Code trace:
- `em` (`clones.py:826`) assigns cells with `loglik` (`clones.py:797`) on all bins and all SNPs, and in
  thin mode on the full counts (module docstring: "EM, profiles and cell assignment use the full counts").
- `merge_unsupported` (`clones.py:702`) then takes `labels = R.argmax(axis=1)` and tests pairs on odd bins or
  the test half.
- The test data therefore helped form the groups it judges. That biases towards keeping clones apart.
- `has_events` in `em` (keep a clone vs fold it into normal) uses the same full-count pools.
- Odd bins or the test half are reused across every split and merge round, with no correction.

Count splitting (`split_mode='thin'`) makes the proposal and the test independent at proposal time only. It
does not fix any of the reuse above.

Proposed (thin mode only; diff below). Merge-test groups come from a data-only E-step on the proposal half,
and EM keeps the full counts. Even-odd mode would need an even-bin-only likelihood; that is not proposed.

Required regardless. A 0.95 posterior is not false-split control. Measure the whole adaptive procedure's
false-split rate on specimens with no subclone, crossing depth, cell state, reference mismatch and spatial
autocorrelation. Report that rate, not the posterior.

### 4. Methods not scored against identical targets: CONFIRMED (both parts), size unknown

Evidence:
- `benchmarks/kinker/evaluate.py:41-59`: per-gene truth and the arm median use the run's own
  `genes.tsv.gz`. The protocol states this ("arm truth = median over the run's used genes").
- Comparator wrappers write `genes_df`, all positioned autosomal genes (`cmp_common.write_layout`).
  HaploDose writes its expression-filtered set. So arm truth and the arms eligible at `MIN_GENES` differ
  between tools. Within one run, HaploDose and pseudobulk share genes, so that pair is fair.
- `run_r_field.py:115-116`: unmapped genes and NaN are set to 0.
- `run_r_field.py:126-127`: under `--rule all`, used for the implant design (E2), cells absent from the
  tool's matrix are clustered from zero profiles. They then count as covered in E3
  (`analyse.py::coverage`). Under `--rule cluster` (E1) absent cells are not analysed.

The effect on rankings has not been estimated [Unverified].

Proposed: **(b) a post-hoc sensitivity analysis, declared now before any comparator output exists.** Leave
E1/E2/E3 and the frozen scripts unchanged. Add a dated protocol note specifying:
- (i) E1 rescored for every tool on a common gene universe per pool and rotation. This is the intersection of
  the genes used by all six methods, with arms eligible at >= 5 common genes, the same truth formula and the
  same tests.
- (ii) E3 coverage restricted to query cells present in each tool's native matrix (`native/cells.txt`). The
  n_in_matrix / n_sent fraction is reported per tool.
- (iii) E2 repeated with absent cells as unlabelled. It is descriptive.

Why not (a), a pre-run amendment. The HaploDose held-out jobs were kept at the freeze and may already have
written output. Changing the primary endpoint then looks like tuning, even though the change is symmetric. The
issue is a fairness and comparability weakness with an unknown direction, not an invalid primary.

Option (c), the v-next eval. Fix the scoring universe inside evaluate.py and add a `native_present` column to
wrapper outputs.

### 5. Superiority tests ignore shared reference pools: PARTLY

Real: Wilcoxon and the line bootstrap (`analyse.py:43-58`) assume independent lines. Lines share 6 reference
sets (3 compositions x 2 rotations).

Already handled:
- Amendment 3a discloses this for the CIs.
- The pool-level Wilcoxon is computed, with a minimum two-sided p of 0.03125.
- Holm handles multiplicity only, as the review says.

Mitigating [Inference]:
- E1 tests paired differences between tools that use the same reference. A reference-composition error moves
  both tools' rho, and the paired difference partly cancels it.
- The residual within-set correlation of the differences is unknown.

Proposed: **(b) a post-hoc sensitivity analysis plus wording, no change to the primary.**
- Report `superior` as "conditional on independence of lines".
- Add a cluster sign-flip permutation over the 6 reference sets of per-set median differences
  (2^6 = 64 patterns, so p >= 1/32) and a line-level linear mixed model with a random intercept per
  reference set. Both are descriptive.

Stronger confirmation needs more independent datasets, as the review says. (a) is not warranted: with 6
clusters there is no better confirmatory test to swap in.

### 6. Zero usable SNPs crash allele inference: CONFIRMED

Evidence. Before the fix, `switch_probs` and `SpecimenAlleles.fit` on an empty SNP set raised `IndexError`
(job 22612820). The pipeline builds `SpecimenAlleles` even when `snp_cols` is empty (`pipeline.py:139-141`).

Fix:
- `switch_probs` (`alleles.py:170`) returns empty arrays.
- `clone_imbalance` (`alleles.py:218`) returns zero-width posteriors, B of shape [K, 0], nu NaN and
  log-likelihood 0.
- Every downstream test then sees no allele evidence.

The reason is recorded in `specimens.tsv`: `n_snps` = 0 and `phase_switch_nu` = NaN. pipeline.py was not
edited.

Tests:
- `::test_empty_snp_set_gives_empty_allele_channel`
- `::test_pipeline_runs_when_specimen_has_no_snps`, end to end with allele reads only in reference cells.

### 7. Output labels assert too much: CONFIRMED (cnLOH), PARTLY (bulk)

cnLOH. Relative-neutral dosage plus imbalance is also 2:1 in a near-triploid genome, or a subclonal mixture.
This was already raised in `docs/efficiency_review_2026-11-01_followup.md`.

Change. `joint_state` (`alleles.py:329`) now returns `neutral_allele_imbalanced` (a cnLOH candidate), and the
docstrings were updated. No benchmark script or test reads the value `cnLOH` (grep, excluding `dev/` and
`results/`). `docs/STATUS.md:611` and `docs/manuscript_outline.md:40` still say "cnLOH". The docs are left to
the lead.

Bulk. The flag has a basis: a Gaussian simulation of the median/MAD scale, cited in the `bulk.py:122`
warning. It is still a panel-size rule, not calibration on the data. The key was kept because
`tests/test_bulk.py` reads it. Added in `bulk.py:192-197`:
- `arm_calls_calibration_basis`, stating panel size only.
- `normal_arm_call_rate_in_sample`, labelled optimistic because the normals define the null.

Test: an assertion was added to `tests/test_bulk.py::test_bulk_recovers_arm_events`.

`normal_like`. The column name was kept. It means "assigned to the reference-like clone 0", and the
manuscript should word it that way.

### 8. GPU init errors: already fixed

`device.select` (`device.py:39-46`) now catches any Exception only when device is 'gpu'. It re-raises as a
RuntimeError chained with `from err`, so the cause is kept. Errors on the 'cpu' and 'auto' paths propagate
unchanged. Nothing else remains.

`test_gpu_without_visible_gpu_fails_clearly` and `test_auto_falls_back_to_cpu_without_gpu` are in the full-suite
run below.

## Computational efficiency

| Item | Verdict | Action |
|---|---|---|
| Depth pooling one-hot (`clones.py:339`, `pooled_depth_slopes`) | CONFIRMED. The arithmetic is right for 40 *clones* (200 strata): 1e6 x 200 x (4 + 8) B = 2.24 GiB. With 40 strata it would be 0.45 GiB. | Sparse one-hot plus `np.bincount`; in the proposed diff |
| `shared_profiles` pairs (`clones.py:517`) | PARTLY. 780 pairs for 40 clones run as one batched HMM; about 47 MB at 3,000 test bins [Estimate]. Pools change every EM round, so a cache would rarely hit. | No change; low priority |
| Gene-level reporting (`pipeline.py`) | Not assessed | None |
| Bulk SVD (`bulk.py:59`) | NOT MATERIAL. The SVD is of a normals x bins matrix with `full_matrices=False`. | None |

## Manuscript direction

This section agrees with AGENTS.md and STATUS: calibration over expansion, and Numbat as precedent. The
CRC/IGF2 point is right: relative RNA composition cannot detect uniform doubling. No action.

## Proposed changes for clones.py, cli.py and pipeline.py

The full text is `dev/review_response_20261102/clones_proposed.diff`. It is against `clones.py` as of
2026-11-02 about 16:30, which already contains the split_mode work.

The diff has four parts:
- (1) Finding 1: log `untestable_pair` in `merge_unsupported` and `try_splits`.
- (2) Finding 2: pass `log_dose=s['mean_diff'] / self.beta` to `alleles.veto` in `supported`. `mean_diff` is
  in beta-attenuated log-ratio units (`difference_segments` levels are multiples of beta log 1.5).
- (3) Finding 3, thin mode only:
  - `loglik(doses, fold=None)` scores a count-splitting half. L is scaled by s and phi by 1/s, because a half of
    NB(m, size) is NB(s m, s size).
  - The half's likelihood is expression only.
  - `merge_unsupported` takes its test groups from `posteriors(loglik(doses, prop_fold), log pi)`. These are
    data-only, so no Potts state is touched.
  - Caveat: the thinning size is frozen at model start while phi is refit, so independence is approximate,
    as the module docstring already says.
- (4) Efficiency: sparse one-hot and `bincount` in `pooled_depth_slopes`.

No change is needed in `cli.py` or `pipeline.py`.

How the diff interacts with count splitting:
- (1) and (2) are mode-independent.
- (3) is a no-op in even-odd mode.
- In thin mode, (3) costs merge-test power, because labels come from half the counts. Pair this with the
  null calibration in finding 3.

The diff was run on a package copy (`dev/review_response_20261102/proposed/`). See the test results.

## Tests run

All jobs ran on the cpu partition with `--prefer=v6`.

| Job | What | Node | Result |
|---|---|---|---|
| 22612820 | `repro.py` on the pre-fix snapshot | v6 | Findings 1, 2a, 2b and 6 reproduced |
| 22612872 | `repro.py` after the fixes, then `test_scientific_review_20261102.py`, `test_alleles.py`, `test_bulk.py` and `test_profile.py` | v6 | 19 passed |
| 22612912 | `slurm/build_test.sbatch`: Rust tests, then the full Python suite on the working tree | ne1dc4-050 (v4) | Rust 3 passed. Python 88 passed, 1 skipped, 1 failed. |
| 22612913 | Full Python suite with the proposed `clones.py` diff (`dev/review_response_20261102/proposed/`) | v6 | 89 passed, 1 skipped |

The working-tree failure in job 22612912 is `tests/test_spatial.py::test_potts_improves_coherent_clones`. It is the
known failure on v4 nodes. The same test passed on v6 in job 22612913.

Both runs include `tests/test_device.py`, and it passed. The review's reported failure of
`test_gpu_without_visible_gpu_fails_clearly` no longer occurs.

The skipped test is the GPU comparison, skipped because no GPU was visible.

Not run:
- A GPU run of the device tests.
- Real-cohort reruns. The changes to the veto and `differ` may change clone counts on development data.

## Amendment 5 implementation (lead decision: option b)

`benchmarks/heldout/sensitivity.py ROOT INPUT_DIR` runs after `analyse.py` and writes
`ROOT/analysis/sensitivity/`. It imports the following and does not copy them:
- From evaluate.py: `score_run` and `load_counts`.
- From score_implant.py: `score_run`.
- From analyse.py: `e1`, `paired`, `family` and `query_cells`.

It uses seed 20261103 and 10,000 bootstrap resamples.

How the gene universe is applied. `evaluate.score_run` keeps only genes that are in the DepMap table. Passing
it a DepMap table cut down to the common universe therefore applies the universe to every part of the score:
truth, pseudobulk, the tool profiles and the arm eligibility rule (5 or more genes).

End-to-end check on Kinker pools 16 and 18 (development data) used
`dev/review_response_20261102/sens_check/sens_check.sbatch`:
- Job 22613001: first run.
- Job 22617092: rerun after the mixed-model fix.

The fake root:
- HaploDose: the frozen-14:24 seed run.
- infercnvpy and CopyKAT in every pool.
- SCEVAN only in r0_pool16.
- R inferCNV absent.

Results:
- (a) The common universe has 10,096 to 10,116 genes per pool and rotation. Median rho is 0.860 for
  HaploDose and 0.808 for infercnvpy, close to the primary.
- (a) R inferCNV is untestable (p = 1 in Holm). SCEVAN has n = 9.
- (b) n_in_matrix / n_sent: CopyKAT 8674 of 8676; SCEVAN 4007 of 8676, because it ran in one pool only.
- (b) Under `--rule cluster` no cell is labelled without being in the tool's matrix, as expected.
- (c) 0 to 2 query cells per CopyKAT or SCEVAN run were relabelled to NA.
- (d) The sign-flip test has 4 sets here, so p >= 0.125.

Mixed model. With few sets, the set variance sits at its 0 boundary. statsmodels' default optimizer then returned
one fit with intercept 0 and SE NaN, which caused the sqrt warning. The lbfgs optimizer reported "converged" at
intercept 0 with a lower likelihood.

The script now fits with lbfgs, powell, nm and bfgs. It keeps the converged fit with a finite SE and the highest
REML likelihood, and reports the other optimizers' results. Powell and nm agree. At variance about 0 the
intercept equals the mean paired difference (CopyKAT, (a): 0.0819 both ways).

Interpretation chosen for the amendment-5 decision rule. "Mean paired difference > 0 with bootstrap 95% CI
above 0" is read as a 95% CI of the **mean** (`mean_ci95`). The median-difference CI from `analyse.paired` is
also reported. The lead should confirm this reading before hashing.
