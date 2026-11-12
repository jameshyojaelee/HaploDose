# A simpler prevalence model preserves recovery and rejects known false splits

[Observed pattern, not a guarantee] One shared mixture prevalence, with
depth-specific nuisance means, passes the fixed eight-case screen:all5 original
recoveries resolved-positive and ALL3 known state-null false splits resolved-
negative. Every24 old per-cell predictive arrays reproduces EXACTLY. All10 CPU
jobs complete0:0. This is a development mechanism result, not a whole-workflow
recovery, false-split rate, biological guarantee or production fix. It permits
the separately frozen full-workflow follow-up, including merging.

## Why this simplification was tested

The [preceding unrestricted prevalence correction](LIKELIHOOD_PREVALENCE_2026-11-10.md)
fits discovery better but loses all5 true recoveries. Its eight separately fitted
stratum prevalences can approach0/1. That justifies testing a smaller model,
without choosing a regularization strength or adding a prior. The observation
does not uniquely prove that overfitting explains the failure; state/count-model
misspecification remains an alternative.

Use one pi across all depth strata, retaining their individual location
parameters, common contrast delta and within-component state variance. This
reduces7 free parameters relative to8 stratum prevalences. It is distinct from
the earlier pooled model that omitted depth-specific means. Fixed proposal
q-derived stratum weights, unconstrained stratum weights, and this shared-weight
model remain separate reported experiments. [Design frozen before outcomes](../dev/pooled_prevalence_20261110/DESIGN.md).

Estimand:conditional confirmation RNA predictive advantage of two regional
populations over one continuous regional state. Observation model:NB genes with
a shared latent state per cell, depth-specific nuisance locations and common
mixture prevalence. Additional identifying assumption:clone prevalence is
independent of the outside-event depth stratum in this comparison. Selection,
estimated exposure and biological genotype-depth association may violate that
restriction. It is not established for tumours by these simulations.

Regional transcription/reference mismatch and RNA-equivalent state/genotype
laws remain biological competitors. Independently measured depth/state and
same-cell DNA can test the identifying assumption and clone memberships.
Corresponding DNA profiles at matching resolution test genomic compatibility.
These8 exposed mechanism cases have biological n=0; they do not replace
biological replication or a suitable independent truth panel.

## Exact intervention, numerical checks and unchanged decision

Keep the same counts, cells, genes, proposal q, reference noise, depth strata,
offsets, genomic interval and old discovery-only continuous null. For fixed
locations/contrast/variance, combine all cells' component-likelihood scores and
maximize the one concave prevalence likelihood to exact boundary or adjacent
floats. Same7 starts, variance0 treatment, L-BFGS-B options and24/48/96 quadrature.
No prior, tuning grid, new cutoff or reference selected by confirmation.

This restricted family does not contain every old fixed pi_s model. Its
discovery likelihood may therefore be lower; require no loss of its OWN feasible
shared-prevalence anchor. Record that distinction rather than disabling an
unfavorable inclusion check for the same model. Write all models before C
decode, conserve full=D+C EXACT, predict C conditional on D using one shared
latent state across both folds. No all-cell descriptive refit overwrites the
independent predictive statistic.

The SAME scientific screen requires every true score resolved-positive and
every known false score resolved-negative. The diagnostic remains2x successive-
order difference plus operation rounding, explicitly not a formal quadrature
enclosure or biological interval. Failure rejects expansion/tuning; success only
permits full original workflow and separately frozen fresh error/power tests.
All original .95/bin/min5bins/min30support/gene.01 assertions and simultaneous
one-sided95% specimen false-split upper<=5% with power stay unchanged.
[Parameter and inclusion contract](../dev/pooled_prevalence_20261110/PROFILE_DESIGN.md).

Numerical gate preserves the exact Fraction1/6 optimum, boundary/flat/label
checks, independent full-PMF parameterization/heat-equation checks and convergent
five-point gradients. The shared nuisance score is summed over strata, because
only one prevalence is fitted. These are computational checks, not calibration.
[Numerical report](../dev/pooled_prevalence_20261110/numerics_v1/report.json).

## Every fixed case and result

| Fixed case | Cells | Genes | Predictive mixture advantage | Numerical diagnostic | Screen |
|---|---:|---:|---:|---:|---|
|original20261031|508|100|4.996070058|0.0000518216|retain|
|original20261107|512|93|6.279982730|0.0000716558|retain|
|original20261108|514|94|6.886305412|0.0000709500|retain|
|original20261109|502|93|0.435499301|0.0000413460|retain|
|original20261110|510|94|1.888685950|0.0000527720|retain|
|state202611097002|600|93|-1.960328864|0.000130082|reject|
|state202611097005|600|93|-0.160509711|0.0000945250|reject|
|state202611097007|598|94|-1.217669050|0.0000967301|reject|

No case or sign is omitted. All8 signs resolve under the fixed diagnostic.
Original109, previously rejected by fixed and unconstrained prevalence
safeguards, is now positive without a cutoff change. New discovery prevalence
for the five originals is approximately0.504,0.641,0.403,0.547,0.580, subject to
component-label orientation. Those are learned quantities, not supplied true
clone fractions. Model prevalence is not a DNA cellular-purity estimate.
[All order-specific fits and scores](../dev/pooled_prevalence_20261110/summary_v1/report.json).

The3 negative cases were selected because the previous candidate falsely split
them. Their rejection is NOT an independent false-split estimate or3 independent
biological controls. All5 positives are proposal-level results; they do not
prove final memberships survive the adaptive search and merge. The additional
depth-independence assumption and finite-sample model-selection uncertainty
remain explicit. Production recovery remains2/5 until justified full validation.

## Execution, scope and next gate

Research sources:
[profiled_prevalence.py](../benchmarks/hardening/profiled_prevalence.py),
[likelihood_prevalence.py](../benchmarks/hardening/likelihood_prevalence.py).
The shared mode leaves unrestricted-mode mathematics available; immutable old
snapshots preserve the negative experiment. No production/native/original tests,
default, dependency, shared installation, frozen result or prior claim changes.
The unchanged installed wheel identity is checked; no new full package/GPU suite
is claimed for these research edits.

| Job | Role | Account/QOS | CPU / RAM / requested time | Actual time |
|---|---|---|---|---|
|23157746|shared-prevalence numerical gate|nslab|1 /4GB /1h|9s|
|23157747|original1031|innovation|1 /8GB /3h|147s|
|23157748|original107|nslab|1 /8GB /3h|128s|
|23157749|original108|innovation|1 /8GB /3h|118s|
|23157750|original109|nslab|1 /8GB /3h|152s|
|23157751|original110|innovation|1 /8GB /3h|215s|
|23157752|known false7002|nslab|1 /8GB /3h|247s|
|23157753|known false7005|innovation|1 /8GB /3h|323s|
|23157754|known false7007|nslab|1 /8GB /3h|244s|
|23157755|dependent all8 summary|nslab|1 /4GB /1h|11s|

All CPU/v6&amd, original005,10complete0:0, no arrays/interactiveQOS. Peak process
RSS181132KB, peak scheduler RSS124544K.
[Exact commands](../dev/pooled_prevalence_20261110/jobs.json),
[accounting](../dev/pooled_prevalence_20261110/accounting_final.tsv),
[resources](../dev/pooled_prevalence_20261110/resource_summary_final.json),
[preservation](../dev/pooled_prevalence_20261110/preservation_final.json).

The [full-workflow follow-up](../dev/pooled_prevalence_workflow_20261110/DESIGN.md)
fits the SAME model during discovery preparation, before confirmation, then
retains all existing confirmation/gene/merge checks. Captured-runtime exact
numeric and full unguarded controls must pass before candidate inference.
Its job journal is separate; a submitted/running workflow is not a passed
recovery or fresh power test. No fresh phase8 panel is run unless original
full recovery first passes. The old phase7 error/final100 futility and all
historical negative results remain unchanged. Full goal remains active.
