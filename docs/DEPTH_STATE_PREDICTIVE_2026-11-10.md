# Depth-conditioned state prediction still fails one required recovery

[Observed pattern, not a guarantee] The new depth-conditioned predictive
safeguard retains four original recoveries and rejects ALL four known false
splits, but still rejects original20261109. It fails the frozen screen, so
reject expansion and threshold tuning. All11 CPU jobs complete0:0; no
production, original test,5% error/power criterion or frozen result changes.

## Implemented correction and scientific contract

The preceding shared-state comparison used one global mixture prevalence; the
passing five-case recovery candidate conditions mixture composition on depth.
This correction uses the SAME locked discovery-only depth strata, with
stratum-specific means and soft cell prevalence pi_s=mean(discovery q in s).
This is a cell-prior prevalence, rather than a gene/exposure share. No new
stratum grid or depth threshold is fitted to these outcomes.

The null is a Gaussian continuous regional log-rate state with mean alpha_s
and common variance v. The alternative retains that continuous variance and
adds a shared contrast delta, with component means alpha_s-c and alpha_s+delta-c,
c=log(1-pi_s+pi_s*exp(delta)). Both have arithmetic RNA mean exp(alpha_s+v/2);
delta0 nests the null. Each family estimates its own stratum means/variance
from discovery only; the alternative also estimates delta. pi_s is learned
from the original discovery q, not independently known. This is a predictive
working model, not an independent confidence or identification claim.

Fit the corrected variance parameterization/analytic boundary with no variance
floor, the exact zero-variance NB likelihood and prespecified multistart
L-BFGS-B. Successful finite multistart fits are recorded; global optimality is
not certified. Fit at quadrature24/48/96, save models before opening any
confirmation-count entry, then evaluate joint(D,C) minus discovery(D) likelihood
integrals using the SAME cell latent state. Full counts=D+C exactly.

Estimand: conditional confirmation RNA log-predictive advantage of a selected
regional mixture over a depth-conditioned continuous-state alternative.
Assumptions: adequate stratum structure/latent distributions, transferable
reference-only NB sizes and the resulting count-fold law, with upstream
selection fixed. Regional transcription/reference mismatch is the strongest
alternative to genomic clones and can be RNA-equivalent. Same-cell DNA tests
membership; independently corresponding DNA at matched resolution tests
profile compatibility. Biological n=0 here. This screen establishes neither.

[Design frozen before outcomes](../dev/depth_state_predictive_20261110/DESIGN.md),
[all exact cases](../dev/depth_state_predictive_20261110/CASES.json),
[implemented model and executable screen](../benchmarks/hardening/depth_state_predictive.py).

## Independent numerical checks

Exact Fraction NB mixture probabilities agree with the conditional prediction
to4.44e-16, within the operation-derived8.73e-11 budget. Adaptive SciPy
integration of FULL NB PMFs agrees with the likelihood objective to3.33e-16.
An independent adaptive integral of Gaussian heat-equation derivatives agrees
with all mean/contrast/variance gradients to2.78e-16. Integration-error plus
operation/magnitude budget is1.16733e-9. The12-SD integration tail is explicitly
bounded using PMF<=1, score<=sum(y+r) and the corresponding score-squared/
Hessian bound. No arbitrary universal epsilon is used.

Nested-null, component swap, cell permutation and stratum permutation pass,
along with the preceding independent predictive-integral checks. These
establish numerical agreement for tested fixtures, not sampling calibration.
[Numerical report](../dev/depth_state_predictive_20261110/numerics_v1/report.json).

## Every selected success and failure

ALL5 original cases plus ALL4 known phase6 false splits are locked before
outcomes. Same populations, genes, reference-only sizes, selected chr12 extent
and original common-support depth strata; no favorable region/cell selection.
These are exposed mechanism probes, not population-performance estimates.

| Case | Cells | Genes | Mixture predictive advantage | Numerical diagnostic | Proposed safeguard |
|---|---:|---:|---:|---:|---|
|original_20261031|508|100|4.348085251|5.41296385e-05|retain|
|original_20261107|512|93|4.713200383|6.71703197e-05|retain|
|original_20261108|514|94|4.527262944|6.82762039e-05|retain|
|original_20261109|502|93|-3.396704732|3.98525397e-05|reject|
|original_20261110|510|94|2.863254628|4.93640217e-05|retain|
|state_202611096001|599|92|-2.277543057|9.51967163e-05|reject|
|state_202611096002|600|94|-1.836995519|6.19891822e-05|reject|
|state_202611096004|599|93|-4.167722939|9.54887426e-05|reject|
|state_202611096006|600|94|-10.889295687|6.69733722e-05|reject|

All nine signs resolve beyond the existing two-times consecutive-order
discrepancy plus rounding diagnostic. That is a convergence heuristic, not
a proven quadrature enclosure or biological interval. Positive/negative signs
are a bounded-screen requirement, NOT a5% hypothesis-test threshold.

The strongest required recovery still fails:109 moves from preceding score
-3.923530592 to-3.396704732, not to a passing sign. Previously unresolved
state202611096001 now has a resolved negative score. Removing global prevalence
misspecification therefore does not suffice for the original recovery/error
objective. Do not tune a threshold, delete109, promote the four recovered cases
or equate rejection with proof that the generating clone is biologically absent.
This falsifies the implemented safeguard as the requested complete correction;
it does not prove universal impossibility of recovering109.

[All order-specific models, predictions and decisions](../dev/depth_state_predictive_20261110/summary_v1/report.json).
No full guarded workflow or fresh calibration is run after the failed screen.
Production remains2/5, and the earlier5/5 research candidate still has its
unchanged4/10 state-null error failure. Whole-specimen one-sided simultaneous
95% error upper<=5% with power remains required and unestablished.

## Compute and preservation

All jobs CPU/v6&amd, no arrays/interactive QOS: numerical1CPU/4GB/1h nslab;
case workers1CPU/8GB/2h alternating innovation/nslab; summary1CPU/4GB/1h nslab.
The prior isolated noneditable wheel passes456 Python/4 Rust tests,1 GPU-only
skip; it is reused READ-ONLY, not rebuilt or modified. Source_v1 is immutable;
imports are isolated, wheel/native/source hashes asserted, outputs new and
frozen paths protected. No new dependency.

| Job | Role | Account/QOS | Runtime | Outcome |
|---|---|---|---|---|
|23156747|numerics|nslab|00:00:19|complete0:0|
|23156758|case_00|innovation|00:01:30|complete0:0|
|23156759|case_01|nslab|00:01:34|complete0:0|
|23156760|case_02|innovation|00:00:40|complete0:0|
|23156761|case_03|nslab|00:01:08|complete0:0|
|23156762|case_04|innovation|00:01:48|complete0:0|
|23156763|case_05|nslab|00:01:45|complete0:0|
|23156764|case_06|innovation|00:02:01|complete0:0|
|23156765|case_07|nslab|00:01:05|complete0:0|
|23156766|case_08|innovation|00:01:12|complete0:0|
|23156767|summary|nslab|00:00:05|complete0:0|

Peak scheduler RSS 121060K (23156761.batch). Exact sbatch commands, actual backend,
versions, all model locks, identities and resource records remain in the
[task directory](../dev/depth_state_predictive_20261110/),
[job journal](../dev/depth_state_predictive_20261110/jobs.json),
[accounting](../dev/depth_state_predictive_20261110/accounting_final.tsv).
Only one new research module is added. Original production/native/test code,
thresholds and frozen artifacts stay unchanged. The next performance experiment
tests the residual reference-noise assumption using disjoint adjacent-depth
reference pairs, without choosing a scalar dispersion from the simulator.
