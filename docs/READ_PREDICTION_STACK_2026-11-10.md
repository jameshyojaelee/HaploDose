# Select regional pooling by reserved-molecule prediction

The stack passes 47 of 64 simultaneous comparisons and fails the full
gate. It beats pseudobulk on seven of eight tested endpoints under the
declared simultaneous criterion. Clean low-depth whole-profile recovery
has 33/48 wins, below the required 36. It also loses several comparisons
with the clean mixture, adjacent ratio and regional shrinkage. All 2,304
prior scores replay exactly. The method remains experimental.

## Results

Values are mean squared regional relative RNA log-odds errors over 48
independently generated specimens in each condition. These are conditional
simulation results with zero biological specimens.

| Condition | Depth | Stack whole | Pseudobulk whole | Stack minority | Pseudobulk minority | Whole-profile mixture minority | Individual minority |
|---|---:|---:|---:|---:|---:|---:|---:|
| Two populations | 300 | 0.003884 | 0.004175 | 0.058452 | 0.084563 | 0.075994 | 0.111857 |
| Two populations | 3,000 | 0.000501 | 0.003482 | 0.005858 | 0.084650 | 0.002369 | 0.010277 |
| Two populations and state | 300 | 0.008888 | 0.016113 | 0.066538 | 0.114646 | 0.084399 | 0.106568 |
| Two populations and state | 3,000 | 0.002503 | 0.015395 | 0.013343 | 0.113629 | 0.067475 | 0.010667 |

Every minority comparison with pseudobulk has 48/48 wins. At high depth,
the stack also wins the whole-profile comparisons with pseudobulk in 48/48
specimens in both heterogeneous conditions. The low-depth state whole-profile
comparison has 48/48 wins. In clean high-depth specimens it loses to the
ungrouped mixture in 48/48 whole-profile and minority comparisons. Under
high-depth state it loses minority error to the individual ratio in 37/48
specimens. Stronger simple or specialized baselines remain consequential.

The stack's homogeneous whole errors are 0.001449 and 0.000148 at the two
depths, compared with pseudobulk's 0.000792 and 0.000079. At high depth it
retains 79.9% of the clean generated contrast and 71.1% with state. Increased
RNA accuracy can also increase state-free error. For example, at high depth
with state, the stack's state-free error is 0.007736 and pseudobulk's is
0.003714. The better RNA score does not resolve the cause of the response.

## Mechanism and statistical contract

The local mixture loses clean shared-population information while reducing
continuous-state error. This experiment chooses convex combinations of six
RNA probability forecasts separately for each grouped region. The experts
are a pooled composition, two whole-profile mixtures, two regional mixtures
and a regional beta-shrinkage model. The mixture sizes remain two and three.
All experts use the original count observations, half-count regularization,
chromosome boundaries and 96-iteration budget.

The estimand is each cell's relative regional RNA log odds. Under the fixed-total
multinomial sampling experiment, a random fixed-size allocation gives two
molecule halves that are independent conditional on the latent composition.
Each expert uses one half to forecast the other, with prior parameters
estimated from another cell fold. Swapping both molecule halves and cell-fold
roles supplies the weight-fitting observations. SLSQP minimizes the regional
binomial cross-entropy subject to nonnegative weights summing to one.
The output records numerical gradient-gap diagnostics and solver fallbacks.

For a target cell fold, the experts are refitted on the two other folds.
Its complete counts determine posterior RNA forecasts. Grouped probabilities
are normalized, assigned to member regions by training count shares, and
averaged over the two original boundary phases. The method conditions on
the supplied folds, cell IDs, observation law, references and genomic order.

The competing explanation is that improved RNA prediction follows cell
state. That possibility is explicit in the stress generator and state-free
error endpoint. The molecule-splitting identity does not justify thinning
biological NB counts with plug-in dispersion. Suitable corresponding DNA
and biological replication remain necessary for genomic claims.

## Falsifiers and independent checks

The original 48 independent specimens per condition and paired depths are
preserved. New seeds randomize event starts, with 20% minority prevalence
and overlapping continuous state. Whole-profile and minority RNA errors
define 64 simultaneous comparisons against all eight preceding methods.
Every exact one-sided Bonferroni lower bound must exceed 0.5. Ties count
against the stack. No comparator is removed. Homogeneous errors, state-free
errors, contrast retention, selected weights, numerical failures and
iteration doubling on the first three seeds per condition are retained.

Rational arithmetic independently enumerates fixed-size splitting for a
four-trial example. It must equal the product of two Binomial(2,p) laws.
The independent reader inspects the actual optimization callback and checks
its analytic derivative against exact rational derivatives at two weight
vectors. Its 32-epsilon arithmetic bound covers the short products, sums and
divisions on that fixture. An asymmetric known optimum tests that the
optimizer moves from its starting weights. Strong convexity and the recorded
gradient gap determine the parameter allowance. These numerical checks make
no statistical or biological calibration claim.

The reader independently recounts all specimen wins and verifies the
64 decisions by the exact binomial-tail inequality
`sum(comb(48,k), k=wins..48) * 1280 < 2**48`. The required count is 36/48.
Iteration sensitivity matters because the preceding local-mixture experiment
changed 42 of 108 checked selection records when the budget doubled. A
performance score cannot establish convergence of its component fits.

The two independently derived gradient discrepancies are zero and
5.55×10⁻¹⁷, within the 7.11×10⁻¹⁵ operation allowance. On the asymmetric
fixture the first weight is 0.750000000048 versus exact 0.75. The objective
and strong-convexity calculation give a 0.0000132 numerical parameter bound.
These are numerical checks on specified finite examples.

The maximum gradient-gap diagnostic over the new study is 1.32×10⁻⁸.
There are 1,113 optimizer failures that trigger the predefined pooled
fallback. Doubling the component-fit iteration budget changes predicted
probability by up to 0.002675. These are additional unresolved numerical
limitations. Neither numerical success nor the 47 passing contrasts
establishes a release-wide accuracy or uncertainty guarantee.

## Remaining decisive questions

The optimizer failures need a saved-prediction-matrix comparison with an
independent convex solver before another complete performance run. Retain
the current feasible weights, solver messages and fallback behavior. A
solver success flag alone cannot establish that its forecasts minimize
the chosen loss or are stable across equally good weights.

The exported loss also deserves a matched analysis. The study fits
probability forecasts and exports their log odds. Conditional mean
probability and conditional mean log odds are different estimands. For a
Beta(2,3) posterior, the former has log odds `log(2/3)` while the latter is
exactly `-1/2`. Expected log odds minimizes posterior squared log-odds
loss. An independent beta-integral calculation followed by a matched
posterior-action ablation can distinguish a loss mismatch from lack of
measurement information. Both actions must retain the same data, prior,
reference and scoring populations. The current results supply no evidence
that this proposed correction passes the full performance gate.

## Execution

Main job 23189322 completes `0:0` in 13m02s on nslab CPU with 1 CPU, 8 GB
and one hour. Reader 23189332 completes `0:0` in 10s with 1 CPU, 4 GB and
15 minutes. The scheduler samples 80,332 KiB peak RSS for the main job.
Both use the validated read-only Conda environment. Source hashes
cover the executed stack and its exact local, adjacent and shrinkage
dependencies. Frozen datasets and earlier study outputs supply regression
checks. There are zero biological specimens.

[Frozen design](../dev/read_prediction_stack_20261110/DESIGN.md),
[executed stack](../dev/read_prediction_stack_20261110/study.py),
[submission](../dev/read_prediction_stack_20261110/submission.json),
[independent reader](../dev/read_prediction_stack_20261110/summarize.py),
[reader submission](../dev/read_prediction_stack_20261110/summary_submission.json),
[main report](../dev/read_prediction_stack_20261110/run_23189322/results/report.json),
[all comparisons](../dev/read_prediction_stack_20261110/run_23189322/results/comparisons.tsv),
[all specimen scores](../dev/read_prediction_stack_20261110/run_23189322/results/scores.tsv),
[iteration study](../dev/read_prediction_stack_20261110/run_23189322/results/iterations.tsv),
[independent derivatives and summaries](../dev/read_prediction_stack_20261110/summary_23189332/results/SUMMARY.json),
[job accounting](../dev/read_prediction_stack_20261110/accounting.tsv),
[earlier local tradeoff](LOCAL_COUNT_MIXTURE_2026-11-10.md).
