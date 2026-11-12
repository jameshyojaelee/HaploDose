# A stronger cell-specific profile candidate, with a remaining REFCON gap

[Observed pattern, not a guarantee] On the exposed scONE specimen, the fixed
gene-trunk / expected-log-abundance candidate raises mean whole-profile Pearson
from the existing HaploDose stack's 0.419 to 0.574, versus pseudobulk 0.561.
Mean cell-specific Pearson rises from 0.227 to 0.298. Reference-refined REFCON
still leads at 0.616 and 0.376. Whole-profile improvement over pseudobulk occurs
in only 161/390 cells, despite the higher mean; cell-specific improvement occurs
in 370/390. This is a measured development improvement, not population superiority,
a release promotion, or completion of the SOTA objective.

## Measurements, estimand and identifying assumptions

The specimen is the previously exposed scONE release fixture: 390 annotated
tumour cells, 450 reference cells and processed corresponding same-cell DNA.
RNA consists of rounded deposited expression, not verified raw UMI/read counts.
Keep the original identifiers, gene mapping, axes, reference roles and 572
5 Mb/arm/gap regions. Candidate fitting is conditional on these exposed features
and roles; target-excluded prior fitting is not an independent complete-workflow
evaluation. No additional biological specimen is supplied.

The initial regional estimand is relative RNA log odds under a working binomial
regional margin with variable observed library size. The anchor follow-up uses
relative log abundance, matching a multiplicative RNA-response model. Neither
quantity establishes absolute DNA copy number. The source DNA remains a processed
relative-depth reference measurement, not error-free copy-state or clone truth.
Corresponding DNA supports this regional compatibility comparison; it does not
validate membership probabilities, 5% false-split control or segment intervals.

The anchored model assumes RNA rate proportional to `reference_g * A_g * D_ib`,
where A_g is shared within the query population and D_ib is a cell-specific
regional response. Within-query count contrasts can cancel A_g; gene-median
pseudobulk supplies a robust shared RNA profile relative to the external reference.
Gene-by-cell state interactions, reference composition, mapping bias, purity and
ploidy can violate that interpretation. Independent RNA replication, gene-level
concordance and suitable DNA measurements would distinguish these explanations.
The present test measures compatibility and does not resolve all of them.

## Inference changes and matched controls

1. Aggregate RNA counts exactly into the frozen regions. Gene/library totals
   must be conserved. For each target-excluded fold, estimate latent probability
   variance from the mean of `Y(Y-1)/(N(N-1))` minus cross-cell products of `Y/N`.
   This uses variable libraries, here 4,889–61,245 counts. It assumes independent
   exchangeable latent probabilities independent of the conditioned library size.
   The finite ratio is not unbiased and depth-associated state can violate it.
2. Apply the existing half-count-regularized beta mean-probability or expected
   log-odds action. Moment variance is constrained to its probability range;
   negative estimates give complete pooling. All unmodified moments are saved.
   There is no DNA-selected variance floor or threshold. Three folds use seed
   20261110; each prior excludes the entire target fold.
3. The first transfer reaches 0.519 whole and 0.285 cell-specific Pearson.
   Regional count pseudobulk itself reaches only 0.451 whole after smoothing,
   against gene-median pseudobulk 0.561. This motivates separating the shared
   gene profile from within-specimen count deviations.
4. Construct a gene-median trunk using the other two folds and original reference.
   Add `log(p_cell) - log(pooled_training_region)` with coefficient exactly one.
   Compare unshrunk counts, the same-prior mean-probability action and `E log p`.
   The trunk-only control is mandatory. For Beta(a,b), the last action is
   `digamma(a) - digamma(a+b)`. Check Beta(2,3) against the exact value −13/12.
5. Smooth each complete combined log profile with the existing five-region
   within-block median. Components are not smoothed separately. No width or
   coefficient is fitted to DNA. The score grid is 5 Mb; smoothing spans up to
   five regions and does not demonstrate 5 Mb event detection resolution.
6. Add three direct log-abundance controls with the **same saved priors**, counts,
   units and smoothing. These refit nothing and isolate the fixed trunk/reference
   correction from the earlier change from log odds to log abundance.

Every candidate is finite on the original support. No method disappears to
improve its verdict. The original five methods, eight direct regional exports,
four anchored/trunk exports and three matched-unit controls are all retained.
RNA predictions and priors are hashed before each executed DNA scoring step.
The specimen and historical scores were already exposed; this ordering does
not turn the development exercise into prospective evaluation.

## Every method on identical DNA support

All values are means over the same 390 cells and 11,260,470 common gene pairs,
summarized to the original regions. Biological replication is **one specimen**.
The cell-specific score normalizes each cell to mean one on callable regions
and subtracts the other cells' mean, using the frozen scorer. All historical
unscaled residual scores are also retained in the output tables.

| Method | Whole Pearson | Spearman | Cell-specific Pearson |
|---|---:|---:|---:|
| Released REFCON | 0.432251 | 0.417711 | 0.371134 |
| REFCON with reference | 0.615822 | 0.557644 | 0.376007 |
| Existing HaploDose stack, dual | 0.419387 | 0.412188 | 0.226796 |
| Existing naive individual | 0.387827 | 0.393342 | 0.198905 |
| Existing pseudobulk | 0.561258 | 0.399008 | −0.048242 |
| Regional individual, raw | 0.098032 | 0.430623 | 0.161464 |
| Regional individual, sm5 | 0.489112 | 0.490268 | 0.280798 |
| Regional pseudobulk, raw | 0.351111 | 0.341724 | −0.212092 |
| Regional pseudobulk, sm5 | 0.450583 | 0.371971 | −0.201401 |
| Regional shrinkage, raw | 0.381558 | 0.434208 | 0.233104 |
| Regional shrinkage, sm5 | 0.511370 | 0.490435 | 0.283596 |
| Regional expected log odds, raw | 0.393159 | 0.437616 | 0.238888 |
| Regional expected log odds, sm5 | 0.518899 | 0.493106 | 0.284654 |
| Cross-fitted gene trunk, sm5 | 0.561442 | 0.403155 | −0.010205 |
| Anchored individual, sm5 | 0.532637 | 0.517566 | 0.293051 |
| Anchored shrinkage, sm5 | 0.568466 | 0.527571 | 0.296602 |
| Anchored expected log abundance, sm5 | 0.573712 | 0.528837 | 0.297868 |
| Direct individual log abundance, sm5 | 0.489067 | 0.490252 | 0.280762 |
| Direct shrinkage log abundance, sm5 | 0.511320 | 0.490430 | 0.283561 |
| Direct expected log abundance, sm5 | 0.518868 | 0.493094 | 0.284618 |

The matched-unit expected-action comparison attributes +0.054844 whole and
+0.013250 cell-specific mean Pearson to the fixed anchor/reference correction
in this fixture. Within the anchored family, expected log abundance improves
whole Pearson by +0.041075 over unshrunk counts, with 381/390 cell wins;
its cell-specific gain is +0.004817 with 229/390 wins. Shrinkage contributes
beyond anchoring and the common smoother here. The gain over the trunk-only
mean is +0.012271 whole and +0.308073 cell-specific.

These are descriptive paired effects with **no patient-level p-values or
superiority confidence interval**. The mean improvement over pseudobulk is
concentrated: 161 cell wins and 229 losses for whole Pearson. The corresponding
REFCON comparison has 84/390 whole and 42/390 cell-specific wins. All negatives
remain. Do not replace this evidence with a claim of universal superiority.

## Independent checks and preservation

The separately running **default-primal** fit 23188436 subsequently completes
`0:0` in 2h37m17s; its nslab scorer 23188477 completes `0:0` in 58s. Its mean
whole/cell-specific Pearson is 0.419387457983/0.226795785970, agreeing at the
displayed precision with the dual baseline. The full check 23190558 verifies
all input/profile identities and geometry, both simple-baseline profile matrices
exactly, and all 1,560 non-stack comparator cell records exactly. The stack
profiles are **not byte-identical**: maximum difference is 1.63×10⁻⁸, with
maximum cell-score discrepancy 1.24×10⁻⁹. Both artifacts remain preserved.

Against the actual default primal, the anchored candidate has mean paired
whole/cell-specific gains +0.154325/+0.071072, with 385/390 and 302/390
cell wins. Thus the reported development gain does not depend on substituting
the optional dual solver for the default. This is still one exposed specimen.

- The first transfer exactly reproduces all 1,950 original cell records.
  The anchor reproduces all 5,070 preceding records; the matched-unit controls
  reproduce all 6,630 preceding records, five fields each. These include the
  unchanged callable-window counts and all four numerical score fields.
- Exact rational enumeration verifies the unequal-depth moment expectations;
  scalar fixtures check every moment component, shrinkage limits and target-fold
  isolation. Independent 70-digit calculations of 18 real-data moments use
  positive cross-cell products. Maximum discrepancies for mean, factorial second
  and cross moment are 3.47×10⁻¹⁸, 2.54×10⁻²¹ and 1.69×10⁻²⁰, within
  operation/magnitude-derived allowances.
- Twelve independent integral calculations evaluate the same expected log
  probability using the convergent digamma-difference integral, with separate
  adaptive pieces and a sixteenfold refinement. Maximum discrepancy is
  1.78×10⁻¹⁵. Acceptance uses integration errors, consecutive-refinement
  discrepancy and arithmetic allowance, not a biological confidence bound.
- Two rational derivatives independently check the leave-one-read weighted
  solver. They agree within the short-operation bounds. That separate pilot's
  mixed clean-case result remains negative; the real-data transfer does not
  retroactively make its 64-comparison gate pass.

No production package/native file, dependency, original test, scientific threshold
or frozen result is changed. The new code is a runnable research candidate.
Full package/GPU suites were **not run** for these research-only additions.
The three original recovery failures, specimen false-split criteria, power and
uncertainty failures remain open. No membership-refinement transfer is presumed;
the preceding failed profile-to-assignment experiment remains applicable.

## Every compute attempt

All jobs use nslab account/QOS on CPU nodes, one CPU, no array, existing Conda
environments and new per-job output directories. Running jobs are untouched.

| Job | Purpose | RAM/time requested | State / scheduler elapsed | Timed peak RSS |
|---:|---|---|---|---:|
| 23189745 | Locked leave-one-read pilot | 8 GB / 1 h | 0:0 / 1m02s | 87,552 KiB |
| 23189931 | Variable-depth regional transfer | 8 GB / 1 h | 0:0 / 53s | 919,040 KiB |
| 23190102 | Fixed gene-trunk anchor | 8 GB / 1 h | 0:0 / 43s | 920,040 KiB |
| 23190118 | Independent numerical reader | 4 GB / 15 min | 1:0 / 15s | 306,404 KiB |
| 23190127 | Reader serialization correction | 4 GB / 15 min | 0:0 / 19s | 306,272 KiB |
| 23190151 | Matched log-abundance controls | 4 GB / 1 h | 0:0 / 39s | 1,079,532 KiB |
| 23190558 | Completed default-primal comparison | 4 GB / 15 min | 0:0 / 3s | 58,692 KiB |

The first reader completes its assertions and paired table, then fails when a
NumPy integer region ID reaches JSON serialization. The original code, manifest,
log and partial output are preserved. Version 2 changes only that serialized
index to a standard integer, reruns every unchanged check and completes. It
does not relax a numerical criterion. All executed source manifests remain valid.

[First design](../dev/scone_regional_shrink_20261110/DESIGN.md),
[anchor design](../dev/scone_regional_shrink_20261110/ANCHOR_DESIGN.md),
[matched-unit design](../dev/scone_regional_shrink_20261110/CONTROL_DESIGN.md),
[variable-depth estimator](../dev/scone_regional_shrink_20261110/estimator.py),
[first transfer](../dev/scone_regional_shrink_20261110/study.py),
[anchor implementation](../dev/scone_regional_shrink_20261110/anchored.py),
[matched controls](../dev/scone_regional_shrink_20261110/controls.py),
[complete 20-method report](../dev/scone_regional_shrink_20261110/run_23190151/REPORT.json),
[all cell scores](../dev/scone_regional_shrink_20261110/run_23190151/scored_cells.tsv),
[all old scores preserved](../dev/scone_regional_shrink_20261110/run_23190151/PREVIOUS_EXACT.json),
[locked RNA profiles and provenance](../dev/scone_regional_shrink_20261110/run_23190151/RNA_PREDICTIONS_LOCKED.json),
[saved priors and RNA lock](../dev/scone_regional_shrink_20261110/run_23190102/RNA_PREDICTIONS_LOCKED.json),
[independent reader v2](../dev/scone_regional_shrink_20261110/summarize_v2.py),
[all independent checks](../dev/scone_regional_shrink_20261110/summary_23190127/results/SUMMARY.json),
[completed default-primal check](../dev/scone_regional_shrink_20261110/run_23190558/REPORT.json),
[default-primal score report](../dev/refcon_regression_20261110/score_23188477/SCONE_REPORT.json),
[default-primal prediction and resource record](../dev/refcon_regression_20261110/hd_23188436/resources.txt),
[paired cell descriptions](../dev/scone_regional_shrink_20261110/summary_23190127/results/paired_descriptive.tsv),
[preserved reader failure](../dev/scone_regional_shrink_20261110/summary_23190118/study.log),
[control resources](../dev/scone_regional_shrink_20261110/run_23190151/resources.txt),
[post-job environment inventory](../dev/scone_regional_shrink_20261110/environment_inventory_after_jobs.json),
[all job accounting](../dev/scone_regional_shrink_20261110/accounting.tsv),
[preceding same-cell comparator study](REFCON_REGRESSION_2026-11-10.md),
[leave-one-read pilot](LEAVE_ONE_READ_STACK_2026-11-10.md).
