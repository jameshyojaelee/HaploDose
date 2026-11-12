# Better reference prediction is not necessarily the right dosage reference

The constrained NB rate candidate improves held-out composite NB score in both
reference types of the existing conditional fixture. It nevertheless increases
error for the generating pooled epithelial reference. A post hoc expected-score
check identifies an estimand mismatch under depth-dependent composition, rather
than an optimizer defect. Do not replace the reference estimator or claim better
genomic dosage from this result. No production code, default or criterion changes.

## Estimands, observation model and replication

Use the same already exposed raw-count parity fixture. Each reference type has
30 fit and30 heldout cells, from one conditional simulation, biological n=0.
Reference types/cells/bins are not independent donors. The48-bin universe was
already selected using fit references. This is not an unopened evaluation.

The prediction target is summed heldout log working-NB probabilities using observed
library offsets and fit-only dispersion. Bins and observed totals are dependent,
so the product is a composite score, not a correct joint sampling density or
calibrated inference. The reference-profile target is instead the centred log of
the pooled generating reference rate. Fit-population and holdout-population true
rates are evaluated separately. These targets need not have the same optimizer.

The identifying assumption behind a common reference rate is adequate control
of cell-specific composition and exposure. The strongest alternative is a
depth-dependent composition change: variance weighting then changes which
reference population is represented. Independent, representative reference RNA
and covariates would address that assumption biologically. Corresponding DNA and
biological replication would still be needed to establish CNA accuracy.

## Fixed candidates and information separation

Three methods were fixed before execution: the existing pooled count rate with
its existing0.1 count floor; mean cell-normalized rates; and constrained NB
composite-likelihood rates summing to1. Each type uses one shared fit-only moment
dispersion across all methods, with existing production bounds0.001 and10. This
type-specific diagnostic dispersion is not the full pipeline's across-type fit.
No dispersion, feature or method setting was chosen against holdout.

NB fitting uses anchored softmax log rates, pooled/equal starts, a2000-iteration
budget and two numerical gradient tolerances. The objective is normalized by fit
count total. All eight recorded fits report solver success. Function stagnation
can terminate a refined fit; this is not a global-optimality certificate. Recorded
gradients, convergence messages, iterations and both refinements remain available.

Predictions, dispersions, heldout scores and selection are frozen before generating
means are loaded. The independent audit replays the complete analysis with the
original means and then a deliberately distorted mean matrix in a new file.
All rates, dispersions, optimization records, numerical checks, scores and selected
methods are identical; only truth-dependent errors change. This verifies code-level
isolation, not prospective scientific blinding of an already exposed fixture.

## Results

Larger composite score is better. These are totals over30 heldout cells per type,
not significance tests or biological-unit effect estimates.

| Reference | Pooled | Equal cell | NB |
|---|---:|---:|---:|
|type0|-4755.141718|-4760.547757|-4754.256207|
|type1|-4681.216971|-4679.125211|-4678.877271|

NB wins both coarse/fine heldout selections. The complete rankings and all scores
are retained; no candidate is dropped. Five type0 and two type1 dispersion
estimates reach the fixed lower bound; none reaches the upper bound.

After freezing predictions, generating-rate error gives a different result:

| Reference target | Pooled MSE | Equal-cell MSE | NB MSE |
|---|---:|---:|---:|
|type0 fit population|0.0065120|0.0087893|0.0069811|
|type0 holdout population|0.0063917|0.0085828|0.0068325|
|type1 fit population|0.0041799|0.0041027|0.0042401|
|type1 holdout population|0.0058143|0.0049053|0.0056086|

MSE is centred log reference-rate error, not DNA error. NB does not consistently
improve this declared reference estimand. A prediction-score gain cannot be
silently substituted for improved dosage accuracy.

## Why the score can target a different population

This **post hoc oracle diagnostic** uses generating mean counts M and their known
exposure totals D, after selection is frozen. At the pooled rate q=sum(M)/sum(D),
compute the expected NB objective gradient with the already fixed dispersion.
The objective is linear in counts, so replacing counts by M calculates its
conditional expectation for known exposures. This is not conditioning on random
observed library totals or establishing a biological population law.

The projected gradient is resolved nonzero:1.5184e-4 for type0 and1.4783e-4 for
type1, compared with an operation-count numerical allowance2.05e-11. Thus that
expected working objective is not stationary at the declared pooled rate here.
This is not a general claim that correctly specified NB fitting is biased.

Controls preserve or explicitly redefine their reference target:

- With common cell rates and the SAME depths/pooled target, gradients are near1e-17.
- With equal depths and original cell compositions, gradients at the resulting
  pooled target are near2e-17.
- The analytic Poisson score limit at each pooled target is near2e-17; no small
  dispersion constant is invented for this check.

Independent scalar sums of r*(D*q-M)/(r+D*q), analytic softmax gradients and
complex-step scalar objectives agree to at most2.45e-17 across all six conditions.
The interaction between heterogeneous composition and exposure can therefore
change the target of NB weighting in this conditional example. Known means are
unavailable to deployment and never enter the candidate fit.

## Numerical gates and compute evidence

Equal-library analytic-optimum and exact cell-order invariance fixtures pass.
Analytic gradients agree with independently written complex-step scalar arithmetic;
full SciPy NB probability differences independently check optimizer objectives.

The audit recomputes eight heldout score rows at60 and80 decimal digits using
an independent integer-count rising-factorial identity for NB coefficients and
standard-library Decimal arithmetic. All score choices are numerically resolved;
maximum recorded score discrepancy is below1e-10, far below the selected score
gaps. Precision refinement is an arithmetic check, not a biological confidence
interval or proof of formal interval bounds. Twelve generating-rate MSEs are
independently recounted.

| Job | Stage | Account/QOS | CPUs / RAM | Elapsed | Process peak RSS | Exit |
|---|---|---|---|---:|---:|---:|
|23076470|fit/predict/reveal|innovation|2 /4GB|9s|119040KiB|0|
|23077071|precision/leakage audit|nslab|2 /4GB|21s|119088KiB|0|
|23077489|expected-score initial|innovation|1 /4GB|9s|118392KiB|1|
|23077638|expected-score report fix|innovation|1 /4GB|4s|118340KiB|0|

The initial expected-score run passed numeric checks but failed report serialization
on a NumPy boolean. Only conversion to a native boolean changed before rerun;
the failed source/output and traceback remain frozen. No criterion was relaxed.

All executable work used CPU compute nodes and the validated noneditable wheel
from `dev/profile_pool_export_20261109/environment`, SHA256
`2ba2a00107a89bc0aa899c8afaf274134544c47f287336521cce53008280cd42`.
No new dependencies or shared installations. Full package tests were not rerun
because package/native/test source was unchanged.

Task `dev/reference_rate_prediction_20261109` contains frozen designs/source,
`panel_v1/{predictions.npz,prediction_freeze.json,heldout_scores.tsv,
selection_before_truth.json,post_freeze_truth_errors.tsv,optimizations.json}`,
`audit_v1/{high_precision_scores.tsv,report.json,replay_original,replay_perturbed}`,
`target_v2/report.json` and the preserved failed `target_v1`. Complete jobs/resources
are in `jobs.json`. No dependent work remains pending.

## Decision

Retain this as a bounded reference-model diagnostic. No blanket NB weighting
substitution, default or accuracy claim is supported. Any further correction must
specify which reference population it represents and test exposure/composition
adequacy, as well as prediction. More optimization of the wrong working target
does not repair identification. Original thinning2/5, three recovery failures,
failed false-split/power and uncertainty criteria, frozen results and historical
decisions remain unchanged. No DNA-membership, biological-calibration or SOTA
claim follows.
