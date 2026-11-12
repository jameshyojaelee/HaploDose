# Region-excluded residual precision repairs preservation, with two comparisons unresolved

The new normal-only conditional residual weights pass all12 original signal-
preservation safeguards and94/96 simultaneous profile comparisons. The complete
gate still fails. Both remaining comparisons have positive mean improvement but
23/32 strict seed wins against Poisson precision at lower depth, below the fixed
simultaneous requirement. No real-data expansion, default promotion, changed
threshold or biological claim follows.

This makes a concrete advance over the
[previous weighting study](WITHIN_REGION_WEIGHTING_2026-11-11.md): the method now
retains Poisson efficiency while handling extra gene noise without failing the
clean, dosage-only or noise-only safeguards. The estimand remains planted regional
natural-log RNA response under an identity dosage-response simulation. Independent
DNA and biological replication are still required for genomic interpretation.

## Correction and identifying assumptions

Total normal variance mixes state and sampling effects. For each target region,
the new learner predicts its normal gene log contrasts using gene contrasts from
all other regions. The whole target region is excluded. With Gamma rates shared
across molecule halves, retaining the target gene would let it predict its own
fluctuation and underestimate its variability. Region exclusion blocks that direct
route; it does not remove arbitrary cross-region gene dependence.

Fit normal identities0–39 and use40–79 to select the unchanged14-member ridge
path. Average the opposite-half mean squared residual by gene, invert and normalize
within regions. Use these same weights in the full-count regional RNA summary
and its state-prediction target. Query values, DNA and labels choose no weights,
penalties or state coefficients. Residuals remain working estimates that include
state-prediction error and calibration-sample selection. No calibrated variance
or uncertainty interpretation is supplied.

State must have informative within-region markers and a normal/query relationship
that transports. Gene-specific dosage response, internal breakpoints, unmarked
coherent state and transport reversal remain competing explanations or explicit
limits. Noiseless dosage protection is not a finite-count genomic guarantee.

## Frozen evaluation and complete outcome

The unchanged generator uses32 new seeds64–95 across7 laws and2 depths:448 new
condition datasets,32 independent base draws, zero biological specimens. There
are80 normals/80 queries,24 regions and20 genes/region. Retained-gene depths are
1e6/1e7; the balanced arm contrast is log(1.5). Fourteen old cases cover one old
seed's complete law/depth grid and remain regression-only. Every one of their140
original and84 preceding weighted/raw profile arrays, with the associated14+14
parameter files, reproduces exactly. The earlier complete448-case replay remains
unchanged; these14 cases are a bounded new integration preflight.

All462 predictions lock before a separate scorer opens any new truth. Twelve
comparators cross2 state laws,2 depths and2 errors:96 primary comparisons. Each
requires positive mean improvement and a one-sided exact95% lower probability of
strict improvement above1/2, Bonferroni alpha=.05/96. All12 original point-mean
no-increase safeguards against raw gene means remain. No best-per-law choice is
permitted. The candidate must pass all108 decisions.

| Endpoint family | Passed | Total | Decision |
|---|---:|---:|---|
| State-challenge comparisons |94|96|Fail complete family|
| Clean/pure-dosage/noise-only safeguards |12|12|Pass|

Both failures compare with Poisson precision in shared Poisson state at depth1e6.
Whole/cell-specific mean improvements are4.79425724e-7/4.42686886e-7. Both have
23/32 wins and simultaneous lower bound0.41623676, below1/2. Positive point means
cannot override that requirement. Every other comparison passes. The strongest
simple, previous within-region and all three weighted competitors remain present.

Selected depth1e6 whole/cell-specific MSE:

| Condition | Previous within-region | New conditional residual |
|---|---:|---:|
| Clean | .000063248 / .000061070 | .000026012 / .000025080 |
| Shared state | .000079811 / .000076885 | .000034560 / .000033216 |
| Pure dosage | .000066918 / .000064326 | .000028940 / .000027728 |
| Noise only | .00263217 / .00253922 | .00201315 / .00194965 |
| Noisy state | .00291372 / .00281248 | .00230540 / .00222485 |

These are simulated within-law improvements, not observed tumour DNA results.
The planted arm contrast remains around log(1.5). Both depths and all19 methods
remain in the full score table, including raw-summary ablations and the explicitly
truth-assisted oracle. Unmarked/reversed-state cases are reported as limits and
are not silently counted as supported biological regimes.

## Independent checks and compute

The independent reader reconstructs all17,024 error scores in extended precision;
maximum relative discrepancy is4.08e-16. All108 decisions agree with rational sums
and integer binomial tails. All33,264 normalized region weights,331 selected exact
no-action profiles, recorded target-region exclusions and both query perturbations
pass. Perturbing queries leaves all weights, residual variance fits, ridge choices
and selection losses identical. Independent noiseless dosage/state and primal-dual
fixtures pass; excluded-region predictors are exactly unchanged by perturbing that
region's counts. Numerical allowances never relax a scientific decision.

| Job | Work | Progress | State | Elapsed |
|---|---|---|---|---|
|23206195|Numerics, preflight and new fit/score|462/462 ██████████|Complete0:0|3m42s|
|23206219|Independent reader and isolation|108/108 ██████████|Complete0:0|12s|

Both use nslab account/QOS, CPU partition, one CPU,8GB,1h and original Intel v4
node003 with existing Conda Python/isolated imports. Running0;pending0;completed2;
job errors0. No arrays, new dependencies, shared environment changes, public code
or original-test edits occur. All91 validated nonbinary public/data/test identities
remain unchanged. Frozen outputs, historical decisions and scientific criteria
remain unchanged. Original recovery2/5, state false splits14/200 and correlated
selected-interval coverage54/1,000 are not improved by this study.

The remaining two comparisons warrant one separately frozen, fixed-horizon
replication of the *unchanged* method using entirely new simulation draws. The
32 cases above stay failed development evidence and cannot be pooled into that
replication or used to tune weights. Increasing simulation precision will not
supply missing biological replication, suitable DNA truth or transport evidence.
A passing replication would only justify separately designing an appropriate
real-data challenge; it would not establish SOTA or promote a default.

[Frozen design](../dev/within_region_conditional_20261111/DESIGN.md),
[estimator](../dev/within_region_conditional_20261111/conditional.py),
[executed workflow](../dev/within_region_conditional_20261111/study.py),
[source identities](../dev/within_region_conditional_20261111/source.sha256),
[submission](../dev/within_region_conditional_20261111/submission.json),
[full96+12 report](../dev/within_region_conditional_20261111/run_23206195/results/REPORT.json),
[every individual score](../dev/within_region_conditional_20261111/run_23206195/results/scores.tsv),
[all means](../dev/within_region_conditional_20261111/run_23206195/results/means.tsv),
[independent reader](../dev/within_region_conditional_20261111/reader_23206219/RECOUNT.json),
[reader submission](../dev/within_region_conditional_20261111/reader_submission.json),
[public-source preservation](../dev/within_region_conditional_20261111/current_source_check.json).
