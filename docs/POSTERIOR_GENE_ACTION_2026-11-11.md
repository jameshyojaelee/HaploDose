# Posterior gene-response actions: numerical qualification and matched transfer

The scientific goal remains accurate RNA-derived regional genomic profiles and
supported CNA populations at an independently validated resolution. This
experiment changes the point estimator under an existing working RNA model.
It does not introduce a genomic posterior, change inference criteria or replace
the original clone recovery and false-split/power endpoints.

**Final outcome,2026-11-11:** the primary makes small exposed mean improvements
but fails the fixed majority-cell gate,186/390 against required196/390. The
secondary also fails. All independent numerical, DNA-reconstruction, exact
regression and numerical-boundary checks pass. Reject promotion. All six jobs
are terminal:four complete,two failed reference attempts preserved. No public
default,scientific criterion,frozen result or biological validation claim changes.

Historical execution checkpoint,2026-11-11: numerical qualification passes. The matched
exposed scONE fit23207574 is running; reader23207584 awaits its successful
completion. **No new DNA-performance result is available at this checkpoint.**
All prior production recovery2/5, continuous-state false splits14/200 and
correlated selected-interval coverage54/1000 remain unchanged.

## Question and fixed scientific contract

The preceding gene-NB response estimator finds a mode in eta=log(t). With a
Gamma(a,rate=a) prior on response t and a product of gene NB terms, its score is

`sum_g (y_g-E_g*t)/(1+phi_g*E_g*t) + a*(1-t)`.

The density in eta includes the change-of-variable Jacobian. The root is not
generally E[log(t)|y] or log E[t|y]. In the Poisson submodel only, its exponential
equals the posterior mean of t. This distinction is an action/estimand issue;
it does not make the saved mode calculation erroneous.

The primary new action E[log(t)|y] minimizes squared log-response error under
the working model. The prespecified secondary action log E[t|y] targets the
mean amplitude. Keep the same likelihood, gene noise, regional prior, exposure,
reference, trunk, folds and five-region block smoother. The mode control must
reproduce the saved estimator exactly. No DNA-selected action or combination.

The identifying assumptions are fixed exposures, adequate conditional gene NB
noise, prior transport and a meaningful relative RNA reference. Observed-library
conditioning, gene dependence, rounded deposited expression, regional state and
dosage compensation are competing explanations. These remain working/composite
RNA quantities. Corresponding callable DNA can test profile compatibility;
independent biological specimens are required for transport and biological
validation. Neither numerical accuracy nor a favourable RNA objective resolves
the genomic identification problem.

[Frozen design](../dev/posterior_gene_action_20261111/DESIGN.md),
[posterior kernel](../dev/posterior_gene_action_20261111/posterior.py),
[byte-identical original root helper](../dev/posterior_gene_action_20261111/frozen_response.py),
[earlier estimator and negative majority-cell gate](SCONE_GENE_RESPONSE_2026-11-10.md).

## Numerical method, qualification and failed attempts

The log-response density is strictly concave for positive a. Tangent lines bound
omitted mass, absolute log-response moment and amplitude moment. Boundaries expand
until these tails are below arithmetic diagnostics derived from operation counts
and binary64 precision. Two-sided Gauss-Legendre integration records every order
32/64/128/256/512, consecutive changes, cancellation and tail terms. Diagnostics
are numerical checks, not rigorous libm enclosures or confidence intervals.

Independent qualification uses twelve one-gene Tricomi-U moment/derivative
references at80/120digits, two heterogeneous full-PMF references at60/90digits,
conjugate Poisson Gamma/digamma identities and row/gene permutations. Prior
shapes0.02/0.5/1/35 approximate the retained real-data range, with zero counts, unequal
means and zero/nonzero gene noise. All14 cases and conjugate/permutation checks
pass in job23207565. High-precision references must round identically at both
precisions. The final512-vs256 change and independent discrepancy must fit the
unchanged arithmetic/tail diagnostics. No new dependency is installed;
mpmath1.3.0 was already present in the existing Conda environment.

The exact RNA-derived shape extrema are0.01880632681780952 and35.02594314353668;
the scalar fixtures approximate rather than strictly enclose that range. Every
real-input integration separately requires512-vs256 convergence within its
recorded arithmetic/tail diagnostic. The independent real-input reader adds
three complete heterogeneous gene-PMF probes. These checks do not establish a
universal numerical error guarantee for untested inputs.

Job23207529 stops after three passing cases when a special-function intermediate
fails to converge. An exactly zero Kummer-series term needs an explicit finite
zero bound. The repair bounds that intermediate below twice the requested bit
precision while preserving the final80/120-digit agreement. Job23207558 exposes
my incorrect precision attribute `mp.prec` before any case completes; the
correct existing API is `mp.mp.prec`. Both failures and complete executed source
snapshots remain. No scientific gate was loosened.

The passing run records Poisson identity maximum discrepancy6.128431095930864e-14
and exact row/gene permutation agreement in its fixture. Exposure rescaling with
a fixed prior changes the model, so no unsupported exposure-invariance claim
is made. Unequal-exposure references and conjugate identities check that scale
dependence. The large-positive-z evaluation branch uses stable log-addition to
avoid overflow; it is a numerical representation choice, not an effect cutoff.

[Initial failed source](../dev/posterior_gene_action_20261111/attempt1_manifest.json),
[first failure](../dev/posterior_gene_action_20261111/numerics_23207529/output.txt),
[second failed source](../dev/posterior_gene_action_20261111/attempt2_manifest.json),
[second failure](../dev/posterior_gene_action_20261111/numerics_23207558/output.txt),
[qualified source identities](../dev/posterior_gene_action_20261111/source_manifest.json),
[all independent cases and convergence records](../dev/posterior_gene_action_20261111/numerics_23207565/results/REPORT.json).

## Matched performance experiment and falsification rule

Keep all24 previous methods,9360 existing cell records,390 tumour cells,572
regions,three frozen folds,reference identities,gene support and smoothing.
Use existing RNA-only nuisance/prior fits without refitting them. Lock all new
RNA predictions before the unchanged scorer opens corresponding processed DNA.
The complete new family includes primary, secondary and exact mode control.

The primary must improve both whole and cell-specific mean Pearson against the
preceding gene-NB-moment estimator and beat matched pseudobulk in a strict
majority of whole-profile cells. Retain the old anchored-estimator continuation
gate and both REFCON endpoints separately. A secondary-only win cannot promote
the primary. Preserve all failures and abstentions. No prior/noise/seed/reference/
smoothing/threshold sweep follows a failure.

The independent reader reconstructs full gene PMFs, verifies selected real-input
moments at50/80digits, computes DNA correlations in long double, recounts means
and paired decisions using exact rational sums, and requires exact preservation
of every old cell record and the mode-control profile. Real-input probes use the
first target cell and first eligible region per fold, selected without DNA.

Biological n=1 exposed specimen; cells and regions are nested observations.
Reconstructed RNA and uncertain processed DNA limit interpretation. No
independent evaluation, raw-count uncertainty, calibrated genomic memberships
or biological false-split guarantee is claimed. Numerical qualification has
biological n=0. A scientific improvement is not established until the complete
paired results and reader are available.

[Executed transfer](../dev/posterior_gene_action_20261111/study.py),
[frozen transfer identities](../dev/posterior_gene_action_20261111/performance_manifest.json),
[independent reader design](../dev/posterior_gene_action_20261111/READER_DESIGN.md),
[reader](../dev/posterior_gene_action_20261111/reader.py).

## Compute and preservation

The first five jobs use account=nslab,qos=nslab,partition=cpu,one CPU,8GB. The
subsequent decision recount uses the same account/QOS/partition and2GB. No arrays,
shared installations, production/default/test changes or frozen-result writes.
The initial dirty tree,286 preceding ledger values and91 previously validated
nonbinary source/data/test identities are recorded before implementation.

| Job | Work | Wall-time request | Status at checkpoint |
|---|---|---|---|
|23207529|Initial numerical qualification|1h|Failed1:0,11s; special-function reference|
|23207558|First reference repair|1h|Failed1:0,4s; wrong precision attribute|
|23207565|Corrected reference and complete qualification|1h|Complete0:0,17s;65956KiB process RSS|
|23207574|Matched profile action transfer|3h|Running; no new performance result yet|
|23207584|Independent moment/scoring/decision reader|1h|Dependency afterok:23207574|

Existing read-only Conda environment:Python3.11.16,NumPy2.4.6,mpmath1.3.0 in
the passing numerical job. Each job uses isolated `-I`, neutral fresh output
directory, one BLAS thread, hash checks and process resource measurements.
Old/shared production environments and frozen outputs remain untouched.

[Initial working tree](../dev/posterior_gene_action_20261111/initial_git_status.txt),
[initial ledger](../dev/posterior_gene_action_20261111/claims_before.json),
[prior source identities](../dev/posterior_gene_action_20261111/prior_source_identities.json),
[first submission](../dev/posterior_gene_action_20261111/numerical_submission.json),
[second submission](../dev/posterior_gene_action_20261111/numerical_retry_submission.json),
[qualified retry](../dev/posterior_gene_action_20261111/numerical_retry_v3_submission.json),
[transfer submission](../dev/posterior_gene_action_20261111/performance_submission.json),
[dependent reader submission](../dev/posterior_gene_action_20261111/reader_submission.json).

## Completed matched fit and fixed-table rejection

Transfer23207574 completes0:0 in31m56s. The exact mode-control profile and all
9360 previous cell records reproduce. Process peak RSS is1312312KiB;
scheduler sampled peak1021844KiB. All1716 fold-region fits finish, with four
complete-pooling cases and1712 supported-density cases. The latter require
256-vs512 quadrature agreement for their target cells.

| Method | Whole-profile mean Pearson | Cell-specific mean Pearson | Whole wins over pseudobulk |
|---|---:|---:|---:|
|Primary expected log response|0.603150638|0.308352828|186/390|
|Secondary log mean amplitude|0.601251439|0.307195375|182/390|
|Unchanged mode control|0.602558136|0.308224663|181/390|
|Pseudobulk|0.561258387|-0.048241612|Comparator|
|Reference-refined REFCON|0.615822434|0.376007356|Retained comparator|

The primary improves both means by about0.000593 and0.000128 and adds five
cell-level wins, but fails the strict majority requirement of at least196/390.
The secondary declines on both means and also fails majority. These are tiny
conditional mean changes, not a validated end-to-end upgrade or a new specimen
replication. REFCON remains ahead on both means. Reject unchanged primary and
secondary for promotion without prior/noise/reference/fold/smoother sweeps.

Decision job23207853 independently recounts all27 methods and390 cells with
exact rational sums; it completes0:0 in1s using1CPU/2GB/1h on nslab CPU.
It runs separately while the original full-PMF reader continues, to avoid
serializing fixed-table decisions behind slow high-precision integration.
It does not replace any original reader check. At this checkpoint two of three
real-input PMF probes pass; final DNA/correlation and propagated numerical
decision checks remain pending in23207584. No additional scientific promotion
or change to frozen criteria follows from this partial reader execution.

The reader was prepared before execution with propagation of log-action
numerical diagnostics through centering, medians, exponentiation, normalization
and residual/correlation calculations. Mean comparisons and majority decisions
must resolve outside those diagnostics. The initial and intermediate source
versions are preserved. Environment inventory is written by the reader itself;
the submitted batch launcher remains byte-identical. One temporary inventory
batch edit was not executed and is retained explicitly as such.

[Completed predictions and scores](../dev/posterior_gene_action_20261111/run_23207574/REPORT.json),
[exact preservation](../dev/posterior_gene_action_20261111/run_23207574/PREVIOUS_EXACT.json),
[all quadrature diagnostics](../dev/posterior_gene_action_20261111/run_23207574/POSTERIOR_NUMERICS.json),
[prediction lock before DNA](../dev/posterior_gene_action_20261111/run_23207574/RNA_PREDICTIONS_LOCKED.json),
[independent exact decisions](../dev/posterior_gene_action_20261111/decisions_23207853/results/REPORT.json),
[decision submission](../dev/posterior_gene_action_20261111/decisions_submission.json),
[initial reader source](../dev/posterior_gene_action_20261111/reader_initial_manifest.json),
[numerical-boundary reader source](../dev/posterior_gene_action_20261111/reader_numerical_boundary_manifest.json),
[executed reader source](../dev/posterior_gene_action_20261111/reader_manifest.json),
[unexecuted inventory batch edit](../dev/posterior_gene_action_20261111/unexecuted_inventory_batch_manifest.json).

## Final independent verification and disposition

Reader23207584 completes0:0 in9m27s. All three full121-gene PMF probes pass
50/80-digit identical-rounding requirements, with maximum candidate discrepancy
7.595313267216852e-14. Independently reconstructed processed DNA uses11260470
common gene-cell pairs at572 regions. All2340 primary whole/cell-specific
correlations across three new/control methods agree,maximum discrepancy
8.049116928532385e-16. All108 method/metric means and32 paired comparisons are
recounted with rational score sums. Every9360 preceding record remains exact.

Propagating numerical diagnostics leaves exactly186–186 primary and182–182
secondary whole-profile wins. All mean-comparison margins resolve, including
the primary0.000128165 cell-specific mean improvement against a propagated
diagnostic of0.000001268. These are numerical diagnostics, not statistical
confidence bounds. Both original anchored-estimator continuation gates and
both REFCON superiority decisions also fail. The independent exact decision
job agrees with the full reader. No floating-point near-tie changes a verdict.

The model-conditional expected-log action gives a tiny exposed mean gain; it
does not meet the prespecified performance criterion. Retain it as experimental
evidence and reject unchanged promotion. No hyperparameter,noise,reference,
fold or smoother sweep follows. Biological n remains one exposed specimen;
processed DNA and reconstructed RNA remain uncertain reference measurements.
Computational correctness does not establish genomic identification,clonal
membership validity,biological false-split control or interval calibration.

| Job | Final result | Elapsed | Process peak RSS,KiB | Node |
|---|---|---|---:|---|
|23207529|Failed1:0; special-function reference|11s|52868|ne1dc3-026|
|23207558|Failed1:0; wrong precision attribute|4s|54016|ne1dc3-009|
|23207565|Complete0:0; numerical qualification|17s|65956|ne1dc4-038|
|23207574|Complete0:0; matched transfer|31m56s|1312312|ne1dc4-068|
|23207584|Complete0:0; complete independent reader|9m27s|948848|ne1dc4-020|
|23207853|Complete0:0; exact decision recount|1s|Recorded in resource file|ne1dc4-063|

Scheduler sampling differs from process peaks:1021844KiB for transfer and
287900KiB for the reader. Zero scheduler RSS for the short jobs is an absent
sample,not proof of zero use. Every job uses nslab; no running job was disturbed.
The full environment inventory records Python3.11.16,NumPy2.4.6,SciPy1.17.1,
pandas2.3.3,anndata0.12.19 andmpmath1.3.0. No dependency installation or upgrade.

Final preservation verifies all91 prior nonbinary package/native/data/test
identities,all287 preceding ledger values,current execution manifests and both
failed source snapshots. Add one final result entry; do not revise earlier
checkpoint values. Original production recovery2/5,false splits14/200 and
correlated selected-interval coverage54/1000 remain. The full genomic/clone/
spatial/SOTA goal remains unachieved and active.

[Full independent result and numerical decision margins](../dev/posterior_gene_action_20261111/reader_23207584/results/SUMMARY.json),
[all paired effects](../dev/posterior_gene_action_20261111/reader_23207584/results/paired_cells.tsv),
[environment inventory](../dev/posterior_gene_action_20261111/reader_23207584/results/environment_inventory.json),
[complete scheduler accounting](../dev/posterior_gene_action_20261111/accounting_complete.tsv),
[resource inventory](../dev/posterior_gene_action_20261111/resources_complete.json),
[preservation checks](../dev/posterior_gene_action_20261111/preservation_complete.json).
