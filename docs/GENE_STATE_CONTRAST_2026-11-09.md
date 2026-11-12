# Repair a missed continuous-state null, then reject the insufficient safeguard

An actual optimizer defect produced false positive support for a discrete RNA
mixture: fitting shared-state standard deviation could stop at0 even though
positive variance improved the null likelihood. The repaired research helper
optimizes variance>=0 with an analytic boundary derivative. Independent exact
and SciPy checks pass. Corrected comparisons reject the original107 and109
recovery contrasts, so the proposed continuous-state safeguard still cannot
satisfy recovery plus error control. Do not promote it, retune the penalty or
run its larger workflow expansion. Preserve all provisional/failed fits.

## Decisive comparison and competing explanations

The [five-case candidate](CONDITIONED_GENE_SCORE_2026-11-09.md) recovers5/5
but falsely splits4/10 fresh state nulls. The earlier
[NB state safeguard](NB_STATE_RECOVERY_2026-11-09.md) lost107 and recovered3/5.
Test whether replacing approximate NB bin sums with individual gene likelihoods
recovers information that safeguard discarded. Unequal-success-probability
NB genes do not sum exactly to another NB; SciPy convolution provides an
independent counterexample with maximum PMF discrepancy.0155329 against a
derived6.39e-12 rounding allowance.

[Initial design](../dev/gene_state_contrast_20261109/DESIGN.md),
[bounded runner](../benchmarks/hardening/gene_state_contrast.py).
Lock originally lost107 and difficult109 before new outcomes. Reconstruct raw
proposal genes from the exact saved gene-thinning CSR, identify rows/genes by
the existing confirmation record, and use frozen proposal-only8-stratum offsets,
expected rates, common-support cells and reference-only gene NB sizes. The
confirmation count array does not supply fitted observations. No truth is used.
Compare moment-bin and gene-product likelihood on exactly the same observations.

Estimand: discrete regional RNA mixture conditional on Gaussian continuous state
and depth nuisance. Null has1 location/state SD; alternative2 locations, mixture
weight and shared SD. Same BIC2-log(n) penalty, starts and24/48/96/192 quadrature
convergence rule. Conditional independent NB genes/shared Gaussian state are
assumptions, not established biology. Non-Gaussian state, dependent genes and
reference mismatch remain alternatives. Independently corresponding DNA is
needed for CNA/clone identity. Biological n=0 and no new validation cohort.

Both initial arms appear to pass both cases: bin/gene BIC107=21.8503/22.8276,
109=.393272/2.12181. But **all four null SD fits collapse to0**. Both bin and
gene arms now pass107, so gene resolution alone cannot explain improvement
over the prior negative result; conditioning and population also changed.
Do not treat the positive BICs as a recovered biological mechanism.

## Independent optimum falsification

[Boundary design](../dev/gene_state_contrast_20261109/BOUNDARY_DESIGN.md),
[independent probe](../benchmarks/hardening/gene_state_null_boundary.py).
With variance v, at0 the log marginal derivative is mean(S_i^2+H_i)/2, where
S/H are a cell's NB log-rate score/Hessian. Negative loss derivatives at0
falsify a reported boundary minimum. Freeze the fitted mean and evaluate a
prespecified SD grid scaled by observed log-ratio spread at24/48/96 nodes.
This grid is a numerical counterexample, not a biological threshold.

All4 fits miss a better null by much more than quadrature/rounding discrepancy.
For native genes,107 improves.0305826 negative-log-likelihood per cell, enough
to bound old-alternative BIC above by-8.48895;109 improves.00866842, bound-6.58128.
The variance-boundary derivatives are-6.13151/-2.69962. Identical quadrature
signs do not establish optimizer correctness. This rejects the provisional
passing comparison before any further recovery/null/power expansion.
[All grid values and derived budgets](../dev/gene_state_contrast_20261109/boundary_v2/report.json).

## Actual correction and unchanged tests

[Correction design](../dev/gene_state_contrast_20261109/CORRECTION_DESIGN.md),
[repaired research helper](../benchmarks/hardening/nb_state_guard.py).
Optimize v>=0 instead of SD>=0. At positive v use the chain rule. At0 use
negative mean(S_i^2+H_i)/2, with mixture responsibilities for the alternative.
The likelihood at SD0 has no integral: compute its analytic value and exactly
zero SD derivative. No positive variance floor, new starts, new penalty or
favourable threshold. Return SD in the existing parameter arrays and explicitly
record optimized_state_parameter=variance for provenance.

The first new regression23154489 fails its **unchanged exact-zero** SD derivative
assertion because old quadrature leaves rounding residue. Keep that assertion;
fix the analytic boundary rather than add a tolerance. The failed output and
source_v4 remain; its two unstarted blocked descendants are canceled. Source_v5
runs the SAME regression and the original independent NB likelihood/gradient
suite before corrected case fits. No original package tests are changed.

On heterogeneous counts[0,0,8,8], mean4,size1000, exact variance-loss derivative
is -374500/63001=-5.94435008968. Error8.88e-16 is within independently derived
1.69e-12 rounding budget. The original SD trap is reproduced; starting the
corrected optimizer at variance0 escapes to positive variance. Null/mixture
interior finite differences and one-sided boundary convergence pass. The
existing4 independent SciPy PMF/adaptive-quadrature fixtures, gradients and
identical-component identity also pass.
[Exact boundary regression](../dev/gene_state_contrast_20261109/boundary_regression_v5/report.json),
[original independent suite](../dev/gene_state_contrast_20261109/original_numerics_v5/report.json).

## Corrected performance decision

| Original seed | Same-cell moment-bin BIC | Same-cell native-gene BIC | Gene null SD | Result |
|---|---:|---:|---:|---|
|20261107|-9.04472|-9.43792|.116154|reject safeguard as sufficient|
|20261109|-10.15029|-9.90241|.101994|reject safeguard as sufficient|

Corrected null objectives are explicitly checked against the already proven
per-order fixed-mean grid; all pass. Quadrature24/48 resolves the negative
signs. The alternative optimizations retain their previous best likelihoods;
the false advantage disappears because the null is fitted better. This is an
optimization correction, not a newly weakened clone test. The gene-resolution
change does not recover these distinctions under the shared-state model.
[Corrected107](../dev/gene_state_contrast_20261109/corrected_20261107_v5/report.json),
[corrected109](../dev/gene_state_contrast_20261109/corrected_20261109_v5/report.json).

The other3 cases, full guarded workflow and fresh error/power expansion are NOT
RUN after the declared gate fails. The previous full5/5 regression and4/10
state-null failure remain intact, as do all older scientific negative results.
No biological error rate, power improvement, absolute copy number or membership
validation is established by the corrected numerical result.

## Every submitted job

| Job | Role | Account/QOS; CPU /GB /time | Outcome |
|---|---|---|---|
|23154450|initial107 bin/gene comparison|innovation;1 /8 /2h|complete0:0,27s; optimum later falsified|
|23154451|initial109 bin/gene comparison|nslab;1 /8 /2h|complete0:0,26s; optimum later falsified|
|23154481|independent variance boundary/grid|nslab;1 /4 /1h|complete0:0; all4 optima falsified|
|23154489|first exact boundary regression|nslab;1 /4 /1h|FAILED1:0; exact-zero assertion retained|
|23154490|dependent corrected107|innovation;1 /8 /2h|canceled unstarted, failed dependency|
|23154491|dependent corrected109|nslab;1 /8 /2h|canceled unstarted, failed dependency|
|23154495|same exact regression after analytic fix|nslab;1 /4 /1h|complete0:0|
|23154496|existing independent NB numerical suite|nslab;1 /4 /1h|complete0:0|
|23154497|corrected107 with grid assertion|innovation;1 /8 /2h|complete0:0; scientific gate fails|
|23154498|corrected109 with grid assertion|nslab;1 /8 /2h|complete0:0; scientific gate fails|

All CPU compute nodes with captured runtime003, neutral working directory and
read-only validated noneditable wheel. Source_v1 initial comparisons,v2 grid,
v4 failed exact regression,v5 corrected computation. Source_v3 is unexecuted.
[All immutable records/commands/resources](../dev/gene_state_contrast_20261109/).
Production/native/original tests, dependencies, shared installations and frozen
outputs remain unchanged. Full package/GPU suites are not rerun for this
research-only helper correction.
