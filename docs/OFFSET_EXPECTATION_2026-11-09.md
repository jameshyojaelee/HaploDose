# Retained populations and the shared RNA offset

The shared-offset numerator included all specimen counts while its denominator
excluded dropped components. The production correction matches the retained
responsibility mass in both pools. This is a population-accounting correction
to the legacy heuristic, not a new conditional-likelihood optimizer. Broader
conditional-mean replacement passes arithmetic checks but fails recovery and
remains an explicitly unpromoted research helper.

## Scientific contract

The estimand is a median-centred regional residual in RNA log-abundance units,
conditional on retained responsibilities, a reference and fitted clone doses.
The legacy expected pool is the retained exposure times normalized
reference-plus-dose rates. It omits cell depth and programme terms. Its output
is neither DNA dosage nor an independently calibrated uncertainty estimate.

For retained components K, observed counts must be
sum_i(sum_k_in_K R_ik)Y_ij. The old numerator was sum_i Y_ij. This difference is
an arithmetic population mismatch even with no nuisance effects. The repair
does not renormalize retained weights or change the half-count or centring.

The strongest alternative explanation for a regional residual remains a depth,
reference or transcriptional effect. A conditional cell mean includes these
effects; a coherent full inference procedure and independent corresponding DNA
would be needed to interpret its residual as a candidate genomic effect.
Correcting one residual expression alone does not supply that evidence.

## Locked experiments and numerical evidence

Task: `dev/offset_expectation_20261109`. Baseline source, candidate wheels,
environments and every intermediate output have separate identities. Existing
frozen artifacts are unchanged.

The first24-case internal covariate probe has equal exposures but supplied
opposite log-depth covariates. This is not constructor-compatible depth data.
It remains a mathematical probe only. A separate24-case follow-up uses exposures
12240/3060 and computes centred log-depth from those exposures, matching the
constructor. Both probes cross nuisance(None,depth,programme,both), retained
membership(hard,dropped,soft), and exact null or a fixed+17/-17count perturbation.
Counts and mean pools have exact rational representations. Tolerances derive
from an operation-count floating-point budget, with independent Decimal80 logs.

Without nuisance effects, dropping a component creates a false centred offset
of approximately±0.168222 under the old method. Four new production tests check
the retained null, partial soft mass, component duplication and chunking. Six
research checks retain the conditional candidate's constructor-compatible
nulls, perturbations, component symmetries, partial mass and logit gauges.
Those six are scoped to the rejected research method; they do not pretend that
the legacy production mean became conditional. All preexisting tests, including
the original full-workflow recovery assertions, remain byte-identical.

## Recovery falsifies the broader change

All arms use the same five frozen raw-input development cases, original
parameters and six original criteria. Truth joins occur after inference.

| Offset arm | Passed cases | Passing seeds |
|---|---:|---|
| Original actual wheel |2/5|20261031,20261107|
| Full conditional actual wheel |1/5|20261108|
| Retained-population ablation, legacy mean |2/5|20261031,20261107|
| Conditional-mean ablation, whole-specimen numerator |1/5|20261108|
| Rebuilt retained-population actual wheel |2/5|20261031,20261107|

The full conditional candidate's full suite records390 passed,1 failed and1
GPU-only skip; the original primary thinning recovery assertion fails. Its four
Rust tests pass. This is a failed validation, preserved at
`validation_23089312`. The two dependent analyses were explicitly reclassified
as diagnostics with afterany dependence on that terminal failure. They were
not treated as successful validation or promoted.

Physical follow-up23091096 passes48 method/probe checks and six strengthened
research tests, while explicitly recording release_validation_passed=false.
Numerical correctness does not override the failed recovery assertion.

The ablation pattern implicates coupling of nuisance means with the rest of
profile refinement [Inference]. The profile path uses an approximate pooled
depth expectation, group programme removal and shared tumour-profile anchoring;
the cell likelihood uses normalized cell-specific rates including offset.
Changing the offset expectation alone does not make these updates a common
likelihood optimizer. This is an open modelling issue, not a demonstrated
biological mechanism or justification to relax recovery criteria.

## Remaining validation and decisions

Rebuilt population-only wheel validation23092399 passes395 Python and4 Rust
tests, with1 GPU-only skip. Its Python test process peak RSS is2528340KiB;
full validation takes2m34s. A second isolated environment preserves the failed
conditional wheel and its test report. Independent Fraction/Decimal80 audits
23092730/23092731 each pass24 old/new comparisons,48 exact count pools and96
scalar offsets, including the constructor-compatible follow-up. These validate
the conditional arithmetic, not the rejected method's full scientific behavior.

Actual-wheel recovery23092729 completes all five cases and retains2/5 passes;
20261108/109/110 remain failures. Independent audit23092732 reconstructs all150
criteria across25 arm/case combinations from cell identities and difference
tables. This supports preserving the narrow correction, not a recovery claim.
Six paired-regression jobs23092897–23092902 complete600 fits; independent
audit23093073 verifies all paired input hashes,1200 old/new exhaustive membership
matchings and168 simultaneous confidence bounds at Decimal60/80 precision.

| Condition/endpoint | Events | One-sided simultaneous95% bound |
|---|---:|---:|
| Matched NB null: extra groups |0/200|upper0.0397854|
| Correlated-state single CNA: extra groups |190/200|upper0.987177|
| Two CNA populations: multigroup detection |200/200|lower0.960215|
| Two CNA populations: exactly two groups |0/200|lower0|
| Two CNA populations: exactly two with accuracy>.65 |0/200|lower0|

All five paired endpoint event counts are unchanged from the old MAP arm; none
of these binary endpoints is gained or lost. This is not proof that fitted
responsibilities or all numerical outputs are identical. Correlated expression
state still defeats the error target, and multigroup detection is not correct
clone recovery. The numerical audits pass while these scientific gates fail.

The regression reuses exposed simulation inputs and known references; it does
not add independent biological specimens. Retain the5% specimen error target
and report exactly-two recovery and all-cell membership accuracy alongside
detection. Five new simultaneous bounds retain all163 earlier endpoints in a
new168-bound record. The narrower accounting repair cannot be advertised as
solving the three recovery failures, correlated-state splitting or uncertainty.
Certified biological null specimens remain n=0.

## Reproduction and disposition

All22 CPU jobs are terminal:21 complete and the original full-conditional
validation23089312 remains a preserved failure. The13-job population-repair
chain completes without execution errors. Resource headers, accounts, runtimes,
step memory and dependencies are in `dev/offset_expectation_20261109/jobs.json`
and `accounting_final.tsv`. Python full-suite process peak is2528340KiB;
the final independent tradeoff audit takes13.26s with124780KiB process peak.
No GPU was requested for these CPU studies; the inference backend is asserted
inside the runs. No production dependency or shared installation was changed.

Build/test command: `slurm/hardening_validate.sbatch TASK
TASK/candidate_population_v1 TASK/population_environment`, submitted with8CPUs,
32GB,2h,nslab/cpu, after successful environment preparation. Analysis command:
`slurm/hardening_analysis.sbatch TASK TASK/population_environment SNAPSHOT SCRIPT
OUTPUT [arguments]`, with `HAPLODOSE_VALIDATED_ARTIFACT` set to
`TASK/validation_23092399/installed_artifact.json`. Exact argv/source/wheel/native
identities and allocation headers are recorded in each `analysis_JOB/execution.json`.
The final wheel SHA256 is
`711325c2020268d339822afc896efbd36c1f1370751317e32368e7eabb55b93f`.

Keep the minimal accounting repair. Retain `offset_conditional.py` and its six
arithmetic checks as rejected-candidate research evidence. Do not change clone
thresholds, confidence names, recovery criteria or historical outputs to disguise
the failures. The next justified modelling investigation must make reference,
profile, cell likelihood and nuisance updates mutually consistent; swapping
one expected mean has already been falsified as a sufficient repair.
