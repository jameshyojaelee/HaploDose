# Match the regional posterior action to its declared loss

[Observed pattern, not a guarantee] The expected-log-odds action improves seven of eight simultaneous
comparisons against the original same-prior regional estimator. It still
loses clean low-depth whole-profile accuracy to pseudobulk. This is a
mechanism ablation, not a successful release gate or a production change.

## Scientific contract

The estimand is cell-specific regional RNA log odds relative to the same
reference. The observation model is the binomial margin of the existing
fixed-total multinomial generator. The empirical beta prior, reference,
training folds, counts and truth are identical between posterior actions.
The new action is `digamma(a) - digamma(b)`, the posterior mean of log odds;
the earlier action is the log odds of posterior mean probability. Expected
log odds minimizes posterior squared log-odds loss under this beta model.
That algebra does not establish frequentist accuracy under an estimated
or misspecified prior. Continuous RNA state remains an alternative to a
genomic cause, and this experiment has zero biological specimens.

Independent beta integration distinguishes an incorrect implementation
from a loss mismatch. Matched generated errors distinguish a useful action
correction from a prior or information limitation. Independently corresponding
DNA and biological replication remain necessary for genomic claims.

## Numerical result

For Beta(2,3), the exact expected log odds is −0.5. The implementation
returns −0.4999999999999998; independent weighted adaptive integration
returns −0.5. Their allowed discrepancy, 7.92×10⁻¹⁴, combines the recorded
quadrature error and an operation-count floating arithmetic allowance.
Independent squared-risk integrals give an improvement of 0.008936845785,
agreeing with the exact identity `(log(2/3) + 1/2)**2` within 9.20×10⁻¹⁴.
Zero and complete shrinkage limits and exact target-fold isolation pass.
All 1,152 preceding regional-shrinkage scores replay exactly.

## Matched performance

There are 48 new independent simulated specimens per condition, with
paired 300/3,000 read depths, 120 cells, 24 regions, 20% minority and
randomized event positions. Seeds are `9721000 + 1000*condition + replicate`.
Cells and regions are not treated as independent specimens. Values below
are mean squared relative RNA log-odds errors.

| Condition | Reads | Expected log odds, whole | Original action, whole | Pseudobulk, whole | Expected log odds, minority | Original action, minority |
|---|---:|---:|---:|---:|---:|---:|
| Two populations | 300 | 0.004701 | 0.004707 | 0.004205 | 0.066174 | 0.067861 |
| Two populations | 3,000 | 0.001271 | 0.001276 | 0.003483 | 0.015840 | 0.016155 |
| Two populations plus state | 300 | 0.010832 | 0.010970 | 0.015881 | 0.070819 | 0.071689 |
| Two populations plus state | 3,000 | 0.002362 | 0.002374 | 0.015129 | 0.013605 | 0.014032 |

All 32 prespecified ablation contrasts are reported with exact one-sided
Clopper–Pearson bounds and Bonferroni correction over this 32-test family.
Twenty-seven bounds exceed 0.5: seven of eight against pseudobulk, seven
against global shrinkage, six against individual ratios and seven against
the original regional action. The independent reader verifies all 32
decisions by exact integer binomial tails. Clean low-depth whole-profile superiority
over the original action has 28/48 wins and fails. Against pseudobulk
that endpoint has 1/48 wins and negative mean improvement. Homogeneous
whole-profile error is also higher than pseudobulk at both depths.

The smaller ablation family does not replace the existing 56/64-comparison
promotion gates or remove any release comparator. The full release gate
is **not assessed** here. The expected-action correction alone does not
resolve the profile failure. No genomic membership, false-split, calibrated
interval, DNA state or biological validation claim follows.

## Execution and preservation

Job 23189477 completes `0:0` in 20 seconds on nslab CPU with one CPU,
8 GB and one hour requested. The timed study uses 4.40 seconds and
82,176 KiB peak RSS. The validated isolated Conda environment and actual
CPU backend are recorded. No frozen artifact or production source changes.

[Prespecified design](../dev/posterior_logodds_20261110/DESIGN.md),
[executed code](../dev/posterior_logodds_20261110/study.py),
[source manifest](../dev/posterior_logodds_20261110/source.sha256),
[submission](../dev/posterior_logodds_20261110/submission.json),
[numerical and replay report](../dev/posterior_logodds_20261110/run_23189477/results/report.json),
[all ablation contrasts](../dev/posterior_logodds_20261110/run_23189477/results/paired_ablation.tsv),
[all scores](../dev/posterior_logodds_20261110/run_23189477/results/scores.tsv),
[environment identity](../dev/posterior_logodds_20261110/run_23189477/identity.json),
[resource record](../dev/posterior_logodds_20261110/run_23189477/resources.txt),
[independent decision recount](../dev/stack_weight_numerics_20261110/summary_23189603/results/SUMMARY.json).
