# Reference depth must identify the intercept, not only the slope

A deterministic counterexample shows that the current query-centred depth
expectation can leave substantial regional effects when the whole query has a
different depth from its reference. With known slopes, a fourfold depth shift
leaves centred residuals of minus/plus log(2), even with one unchanged genotype.
An added true RNA gain contrast0.243 becomes negative,−0.450, in that expectation
comparison. This is not a measured package DNA call or a biological false-split rate.

A conditional reference-intercept calculation removes the known depth effect
while retaining the gain. This is research code with known nuisance inputs;
**a production replacement is not implemented or validated**. No inference default,
historical artifact or scientific criterion changes.

## What is and is not identified

The estimand is normalized RNA rate at the query's exposure distribution after
the specified depth component, relative to a recorded reference population. The
working mean is L*softmax(b+gamma*log(L/anchor)). The deterministic model declares
gamma to be a non-genomic effect. Neither that declaration nor a correct optimizer
establishes the origin of depth association in biological samples.

Identifying assumptions are common depth response between reference and query,
stable reference genotype/composition, known exposure and an identified slope.
The strongest biological alternative is depth-linked cell state or CNA content,
whose removal could erase real genomic signal. Corresponding DNA and independently
controlled technical/reference measurements would distinguish those explanations.
No such biological validation is supplied by these calculations. Reference labels
alone do not certify a homogeneous diploid control.

The pipeline forms base rates by pooling reference counts but centres log depth
within each query. `SpecimenModel` receives no reference-depth distribution.
Its existing expectation-variance test checks averaging around a supplied point
intercept; it does not establish that a pooled reference rate is that intercept.
Knowing a slope alone is therefore insufficient to transport a reference pool.
Fitted offsets/programmes and actual clone search may change downstream behavior;
they were not executed or validated by this expectation-level counterexample.

## Locked forward-model cases

Forty equal-baseline bins have slopes−0.5 on the first10, zero on the middle20
and+0.5 on the last10. Eight equal-depth query cells are represented by their
known mean rates. No observed counts are fabricated or admitted through the raw
count interface. These are four deterministic cases, biological n=0.

| Case | Reference depths | Query depth | Slope |
|---|---|---:|---|
|Matched|9000|9000|specified gamma|
|Shifted|9000|36000|specified gamma|
|Reference spread|4500,9000,18000|36000|specified gamma|
|Zero slope|9000|36000|zero|

The null has a single unchanged genotype. The positive control adds RNA contrast
delta=0.6*log(1.5)=0.243279 to the first10 bins. This fixed response model is a
numerical control, not independent evidence of DNA-to-RNA response. In particular,
the constant-depth references cannot themselves estimate gamma; it is known here.

## Conditional correction

Use reference exposures L_i, their weights w_i=L_i/sum(L), and recorded pooled
reference rate mu. Estimate the anchored intercept b by minimizing
sum_i w_i*logsumexp(b+gamma*x_i)−mu*b. The gradient is the predicted pooled
reference rate minus mu. Anchoring one logit removes the additive constant.
This convex calibration reproduces the reference distribution, then predicts
query expectations at their actual depths. It requires more information than
the current expectation API receives.

Compare the installed current expectation, a reference-origin shift alone, and
the calibrated expectation. For null inputs, maximum centred log-ratio error is:

| Case | Current | Origin shift only | Calibrated |
|---|---:|---:|---:|
|Matched|0|0|0|
|Shifted|0.693147|3.3e-16|0|
|Reference spread|0.579682|0.176620|6.8e-15|
|Zero slope|0|0|0|

Thus changing the origin alone does not solve the pooled-reference problem when
reference depths vary. In the shifted gain control, the first region is−0.449868
under the current expectation and+0.243279 after calibration. Across all calibrated
null/gain cases, maximum error is below7e-15. These are mean-function errors,
not uncertainty intervals, detection power or estimates of population performance.

## Independent computational evidence

The design was frozen before execution in
`dev/depth_reference_anchor_20261109/DESIGN.md`. Existing installed `expected()`
and `depth_moments()` are executed directly in the current-model comparison.
The candidate is `benchmarks/hardening/depth_reference_anchor.py`; it uses SciPy
convex optimization and an analytic-Hessian root polish. Two numerical refinements,
known true intercepts, scalar forward rates and an origin change by log(2) verify
its parameterization. Maximum origin-induced prediction change is3.8e-15;
maximum refinement change is9.4e-15.

Independent audit recomputes320 query-bin means at80 decimal digits, checks all24
contrast rows, independently reproduces reference means and calibration equations,
and verifies positive anchored Hessian eigenvalues in all four cases. Smallest
eigenvalues are approximately0.000625–0.000725. Reference-rate residual is at most
1.91e-16. The shifted-depth contrast also has a separate analytic log(2) check.
Numerical allowances depend on operation counts, scale and known-intercept error;
no scientific cutoff is introduced.

| Job | Stage | Account/QOS | CPUs / RAM | Elapsed | Process peak RSS | Exit |
|---|---|---|---|---:|---:|---:|
|23079034|forward-model/candidate|innovation|2 /4GB|16s|189652KiB|0|
|23079491|initial audit|nslab|1 /4GB|6s|79648KiB|1|
|23079815|label-preserving audit|nslab|1 /4GB|5s|79684KiB|0|

The first audit stopped because pandas treated the literal event identifier
`null` as missing. The corrected reader preserves literal labels and asserts
both expected events before recounting. Failed output/source remain frozen;
no case or numerical/scientific criterion was removed or relaxed.

All execution used CPU compute nodes and the validated noneditable environment
`dev/profile_pool_export_20261109/environment`, wheel SHA256
`2ba2a00107a89bc0aa899c8afaf274134544c47f287336521cce53008280cd42`.
Evidence: `panel_v1/{forward_models.npz,contrast_errors.tsv,report.json}` and
`audit_v2/report.json`. Execution inventories/resource records and failed artifacts
are retained under the task; `jobs.json` records the chain. No dependencies,
shared installations, production/native/test/README files or frozen results changed.
The unchanged full package suite was not rerun.

## Decision and next decisive test

Retain this open model-correctness issue and conditional calibration prototype.
Do not describe reference-derived slopes as a validated correction for a
specimen-wide reference/query depth difference. The next test must estimate slopes
from fit reference counts, reserve heldout reference evidence, record depth support
and challenge slope/model misspecification. A zero slope from inadequate depth
variation is not evidence that no depth effect exists. Query depths outside
supported reference exposure remain an extrapolation.

Only after that test should the reference context be carried consistently through
mixture fitting, pooled initialization, profile expectation, cell likelihoods and
merge/search calibration. An expected() patch alone is not a complete corrected
workflow. Shared changes must pass both profile and clone tracks, including the
original recovery,5% false-split/power and uncertainty criteria. Those failures,
the original thinning2/5 and the absent biological null panel remain unresolved.
No SOTA or biological-validation claim follows.
