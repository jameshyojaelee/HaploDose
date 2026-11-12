# Local regional mixtures improve the state challenge and lose clean borrowing

The local count mixture passes 38 of 56 simultaneous comparisons. On
high-depth simulations with overlapping continuous RNA state, minority-region
RNA error falls from 0.06914 for the whole-profile mixture to 0.01797.
Without added state, the same change raises minority error from 0.00122
to 0.01038. The full gate fails. The fit remains experimental.

## Estimand, mechanism and competing explanation

The estimand is each cell's regional relative RNA log odds. The observation
law is the binomial margin of a fixed-total multinomial count vector. Each
region learns its own two- or three-probability mixture, selected by
held-out cell likelihood against a pool. The fit uses the original two
adjacent partitions, chromosome boundaries, half-count regularization and
96 iterations. The probabilities are normalized across grouped regions,
allocated within groups by training count shares, and averaged over phases.

The tested mechanism is that a whole-profile mixture can allocate its
limited components to continuous expression state. Local components can
represent marginal variation separately. Their cost is losing membership
information shared across distant CNA-shaped regions. This is compatible
with the measured tradeoff. Finite EM fitting is another contributor,
because doubling the budget changes some selected models. A further
predictive combination is being tested using reserved molecule counts.

State and genomic dosage can produce the same RNA probabilities. These
calculations identify predictive RNA behavior under the generator. Suitable
corresponding DNA would be needed for a genomic interpretation.

## Complete conditional comparison

Each condition contains 48 independently generated specimens, with paired
depths, 120 cells, 24 regions and 20% minority prevalence. Event starts are
randomized within chromosomes. The main endpoints are whole-profile and
minority-region RNA mean squared errors. All seven comparators remain in
the 56-comparison family. One-sided exact binomial confidence bounds use
Bonferroni correction. The simultaneous lower win-probability bound must
exceed 0.5 for every comparison. Ties count against the candidate.

| Condition | Depth | Local whole error | Global-mixture whole error | Local minority error | Global-mixture minority error | Individual minority error |
|---|---:|---:|---:|---:|---:|---:|
| Two populations | 300 | 0.003983 | 0.004166 | 0.063181 | 0.078914 | 0.113505 |
| Two populations | 3,000 | 0.000645 | 0.000277 | 0.010378 | 0.001217 | 0.010390 |
| Two populations and state | 300 | 0.009270 | 0.010970 | 0.065295 | 0.082989 | 0.107294 |
| Two populations and state | 3,000 | 0.003700 | 0.006405 | 0.017967 | 0.069138 | 0.010593 |

At low depth without extra state, the local model improves whole-profile
error over pseudobulk in 34/48 specimens. The family requires 36/48.
The local model retains 72.1% of the generated high-depth clean contrast
and 73.7% under continuous state. The global model retains 97.5% and 23.4%
on those same draws. The local high-depth state minority error still exceeds
the individual ratio's error. Regional beta shrinkage also has lower state
whole-profile error, 0.002268. The new method has no uniform superiority.

In the homogeneous controls, local whole-profile error is 0.001333 at depth
300 and 0.000136 at depth 3,000. Pseudobulk gives 0.000765 and 0.000081.
Local fits often select additional prediction components there. These
components are outside the genomic false-split event and carry no clone
interpretation. There are zero biological specimens.

## Numerical and execution evidence

The rational two-probability posterior and direct binomial enumeration pass.
Target-fold perturbation preserves fitted choices and an unchanged peer's
prediction exactly. Canonical marginal count ordering makes row-permutation
predictions exact. All outputs are positive compositions within an
operation-derived binary64 summation allowance. All 2,016 previous score
rows reproduce exactly.

An independent reader recomputes all wins and decisions. At 48 specimens,
the threshold follows from the exact integer inequality
`sum(comb(48,k), k=wins..48) * 1120 < 2**48`. This verifies the simultaneous
decision without relying on inverse-beta floating-point output.

Doubling the iteration budget changes 42 of 108 checked fold/phase selection
records. The largest predicted-probability change is 0.006043. The budget
does not establish a converged regional mixture. Both the failed joint gate
and this qualification prevent promotion.

Initial job 23189258 stops before scoring because its fold-isolation
comparison included a target-dependent mass-correction diagnostic. The
corrected comparison isolates fitted choices and retains exact peer-prediction
checking. It changes no estimator calculation, seed or scientific criterion.
The original failed source is retained. Retry 23189266 completes `0:0` in
1m40s on nslab CPU with 1 CPU, 8 GB and a one-hour limit. Independent reader
23189274 completes with 1 CPU, 4 GB and a 15-minute limit. The exact runtime,
memory and environment inventories are retained with the runs.

[Frozen design](../dev/local_count_mixture_20261110/DESIGN.md),
[executed method](../dev/local_count_mixture_20261110/study_v2.py),
[failed initial check](../dev/local_count_mixture_20261110/run_23189258/study.log),
[initial submission](../dev/local_count_mixture_20261110/submission.json),
[retry submission](../dev/local_count_mixture_20261110/retry_submission.json),
[main report](../dev/local_count_mixture_20261110/run_23189266/results/report.json),
[all comparisons](../dev/local_count_mixture_20261110/run_23189266/results/comparisons.tsv),
[all scores](../dev/local_count_mixture_20261110/run_23189266/results/scores.tsv),
[iteration study](../dev/local_count_mixture_20261110/run_23189266/results/iterations.tsv),
[independent reader](../dev/local_count_mixture_20261110/summarize.py),
[independent means and decisions](../dev/local_count_mixture_20261110/summary_23189274/results/SUMMARY.json),
[accounting](../dev/local_count_mixture_20261110/accounting.tsv).
