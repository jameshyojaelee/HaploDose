# Resolve numerical stacking failures before judging the estimator

The correction resolves all 1,113 original optimizer failures across the
full panel. Three other fits remain numerically unresolved and retain pooled
abstention. The scientific gate remains failed at 47/64. Clean low-depth
wins against pseudobulk increase from 33/48 to 34/48, below the unchanged
required 36. The research correction is retained without a production promotion.

## Selected numerical mechanisms

Selection was locked before new runs: the lexical first failed and successful
regional fit in each of six condition/depth strata, ordered by seed, phase,
fold and region. The exact same generator, cell IDs, molecule halves, priors,
expert forecasts and validation counts reconstruct the saved fit. Both the
old weights and diagnostics must match exactly. All 12 do.

The estimand at this stage is the best convex combination of six supplied
regional RNA probabilities under Bernoulli cross-entropy on reserved
molecules. It is a finite numerical objective conditional on these matrices.
Its identification condition is convex prediction support; individual mixture
weights can remain weakly identified when experts coincide. A biological
interpretation of the resulting RNA response remains separate: transcriptional
state can predict held-out RNA without being a DNA alteration.

Pairwise coordinate updates preserve the two-weight mass and minimize each
line by bracketed root finding or its boundary solution. The simplex convex
gradient gap bounds remaining objective suboptimality in exact arithmetic.
Stopping uses a matrix-specific floating arithmetic estimate based on dot
products, signed gradient sums and local derivative sensitivity. This is
neither a confidence interval nor a universal epsilon. The 1,000-cycle ceiling
is a computational cap; a fit reaching it without the numerical criterion
remains unresolved. The independent rational two-expert optima are exactly
one-half and three-quarters. Both pass.

All six failed SLSQP returns in the probe are finite and feasible after the same explicit
simplex projection. Independent refinement reaches the numerical criterion
for all 12 matrices in one or two cycles. The old failure messages include
iteration limits and positive directional derivatives. Successful old fits
also change by small amounts, up to 1.81×10⁻⁹ in validation probability.
Four of the six old failed fits change validation predictions; the largest
absolute change is 0.006272. Two already have the same pooled optimum.
Equal objective values do not imply equal target predictions.

Falsifiers are failed exact reconstruction, an incorrect known optimum,
infeasible weights, increased loss beyond the arithmetic allowance or an
unresolved gradient gap. The new full-panel path keeps pooled abstention
for unresolved fits. It does not accept a solver status as numerical evidence.

## Full-panel correction

The separate full-panel refinement preserves all 288 previously exposed
simulated specimens, all eight baselines, all 64 simultaneous scientific
comparisons, target-fold exclusion and the 96/192-iteration sensitivity
analysis. All 2,304 baseline scores must replay exactly. It is a matched
regression analysis, not prospective evaluation. The original 47/64 result
remains unchanged in its historical directory.

Main refinement job 23189593 completes with all 2,304 baseline scores exact,
25,920 regional fits and 2,592 score rows. All 1,113 old optimizer failures
resolve. Three previously successful SLSQP fits reach the 1,000-cycle limit
without meeting the new numerical criterion and therefore abstain:

| Condition | Depth | Seed | Phase/fold/region | Gap | Arithmetic estimate |
|---|---:|---:|---|---:|---:|
| Homogeneous | 3,000 | 9641032 | 1/2/7 | 6.13×10⁻¹⁵ | 1.54×10⁻¹⁵ |
| Two populations | 300 | 9642015 | 1/1/17 | 2.30×10⁻¹³ | 2.32×10⁻¹⁵ |
| Two populations plus state | 3,000 | 9643024 | 0/1/6 | 6.73×10⁻¹⁵ | 1.68×10⁻¹⁵ |

The original six-condition simulation has 48 independent generated specimens
per condition, 120 cells and 24 regions. These are zero biological specimens.
The exact same 64 Bonferroni-adjusted one-sided specimen-win comparisons
still pass 47/64. Seven of eight pseudobulk endpoints pass. Clean low-depth
whole-profile comparison has 34/48 wins and simultaneous lower bound
0.471874, failing the required bound above one-half. The exact threshold
remains 36/48 wins. No failed comparator disappears.

[Observed pattern, not a guarantee] Mean RNA errors change as follows.
These are matched exposed regression results, not independent validation.

| Condition | Reads | Prior whole | Refined whole | Pseudobulk whole | Prior minority | Refined minority |
|---|---:|---:|---:|---:|---:|---:|
| Two populations | 300 | 0.003884 | 0.003882 | 0.004175 | 0.058452 | 0.058093 |
| Two populations | 3,000 | 0.000501 | 0.000494 | 0.003482 | 0.005858 | 0.005660 |
| Two populations plus state | 300 | 0.008888 | 0.008894 | 0.016113 | 0.066538 | 0.066521 |
| Two populations plus state | 3,000 | 0.002503 | 0.002498 | 0.015395 | 0.013343 | 0.013191 |

Homogeneous whole-profile errors are 0.001456/0.000149 at 300/3,000 reads,
versus pseudobulk 0.000792/0.000079. Doubling expert EM iterations still
changes probabilities by up to 0.001983. Refinement is not a convergence
proof for those component fits. Better RNA prediction also retains the
genomic interpretation problem: with state at high depth, state-free error
is 0.007740 versus pseudobulk 0.003714. The three original recovery
failures, biological false-split control and selected uncertainty remain
outside this correction.

## Independent assessment and execution

The completed reader independently evaluates saved probe gradients using
70-digit decimal arithmetic on the exact binary64 input values. It requires
the discrepancy to fit the recorded arithmetic estimate. It also independently
recounts every scientific decision by integer binomial tails and checks the
32 posterior-action ablation decisions. All 12 decimal checks, all 64 full-panel
decisions and all 32 ablation decisions pass. Maximum decimal gradient-gap
discrepancy is 2.46×10⁻¹⁸. This validates the stated finite computations;
it does not certify the three unresolved fits or biological claims.

Probe job 23189476 completes `0:0` in 21 seconds on nslab CPU, one CPU,
8 GB and one hour requested. Its timed payload uses 5.42 seconds and
104,180 KiB peak RSS. Main refinement completes `0:0` in 13m10s with
192,000 KiB timed peak RSS, requesting one CPU, 8 GB and one hour.
Reader 23189603 completes `0:0` in 11s with 124,816 KiB timed peak RSS,
requesting one CPU, 4 GB and 15 minutes. All run through
the validated isolated environment. The source and actual native/backend
identities are checked. Existing running jobs are untouched.

The first file-only preparation command used unavailable `python` on the
login node. Repeating that preparation with `python3` succeeded. Scientific
execution remains confined to compute nodes. No production files, frozen
outputs, comparator families or scientific thresholds change.

[Locked probe design](../dev/stack_weight_numerics_20261110/DESIGN.md),
[full-panel design](../dev/stack_weight_numerics_20261110/REFINEMENT_DESIGN.md),
[independent line solver](../dev/stack_weight_numerics_20261110/solve.py),
[exact reconstruction](../dev/stack_weight_numerics_20261110/probe.py),
[completed probe report](../dev/stack_weight_numerics_20261110/run_23189476/results/report.json),
[full-panel implementation](../dev/stack_weight_numerics_20261110/refined_study.py),
[refinement submission](../dev/stack_weight_numerics_20261110/refinement_submission.json),
[full-panel report](../dev/stack_weight_numerics_20261110/run_23189593/results/report.json),
[every scientific comparison](../dev/stack_weight_numerics_20261110/run_23189593/results/comparisons.tsv),
[all fits, including unresolved cases](../dev/stack_weight_numerics_20261110/run_23189593/results/decisions.json),
[independent reader](../dev/stack_weight_numerics_20261110/summarize_refinement.py),
[reader submission](../dev/stack_weight_numerics_20261110/summary_submission.json),
[independent results and all means](../dev/stack_weight_numerics_20261110/summary_23189603/results/SUMMARY.json),
[job accounting](../dev/stack_weight_numerics_20261110/accounting.tsv).
