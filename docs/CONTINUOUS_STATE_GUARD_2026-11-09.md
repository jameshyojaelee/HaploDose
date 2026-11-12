# A direct state safeguard suppresses genuine recovery and is rejected

The new continuous-state safeguard was implemented and tested on all five
original recovery cases. It returns one S1 tumour group in every case and
recovers0/5, versus4/5 for its unguarded research parent and2/5 production.
Reject it as a sufficient correction. No large error/power campaign, package
promotion or biological claim follows from suppressing true distinctions.

## Why this candidate was tested

The [depth-aware gene candidate](DEPTH_STRATIFIED_GENE_NOISE_2026-11-09.md)
fails the frozen fresh state-null experiment13/100, versus1/100 production,
without improving two-profile power(95/100 versus96/100). Its matched null
passes0/100. This pattern implicates the continuous regional expression state
alternative rather than success under a correctly homogeneous working null.
Splitting counts does not remove a biological state shared by both halves.

The [prespecified design](../dev/continuous_state_guard_20261109/DESIGN.md)
tests one direct model comparison after an event has passed all previous
weighted HMM, depth-match and gene requirements. It adds a veto, never lowers
the original.95 support, segment length, gene or cell criteria. Original
tests and production defaults remain unchanged.

## Observation model and interpretation

The estimand is evidence for two discrete locations in a regional RNA-effect
distribution after allowing continuous state variation. It is not genomic clone
identity. Use the actual proposal's cell-by-bin O/E/NB-shape arrays, with source
and cell ordering retained through the executed event function. Define each
cell's approximate regional log rate aslog((sumO+.5)/(sumE+.5)), with delta-method
measurement variance sum(E+E²/r)/(sumE+.5)² from the native NB moments.

Compare one Gaussian latent state(location and nonnegative state variance,
2 parameters) against two locations with a common continuous-state variance
and mixing weight(4 parameters). Both convolve the same cell-specific working
measurement errors. Shared variance prevents a component's variance collapsing
onto one cell. Two fixed initializations use proposal-group means and20th/80th
percentiles. Historical mixing-weight bounds.001/.999 remain. Local SciPy
optimization is not certified global optimization. Optimization failure or an
unresolved numerical model-comparison margin means abstention.

The new veto requires2*(LL_alternative-LL_null)>2log(n), the BIC penalty for the
two extra parameters. It is a working model-selection rule, not a calibrated
5% whole-specimen test; mixtures are nonregular and all prior selection/nuisance
issues remain. Gaussian state, the aggregate log-rate approximation and delta
variance can be wrong. A non-Gaussian transcriptional state may imitate a
mixture, while a real genomic mixture may have nearly Gaussian RNA. Same-cell
DNA would distinguish genomic identity. Biological n=0.

[Guard implementation](../benchmarks/hardening/continuous_state_veto.py)
independently recomputes both fitted likelihoods using SciPy Normal densities.
The operation/scale accumulation budget is about2.1e-8 on these five fits;
observed BIC disadvantages are6.2–10.1, so a numerical near-boundary decision is
not the explanation. These checks do not validate the observation approximation,
optimizer globality, adaptive error control or genomic labels.

## Full recovery outcome

| Seed |S1 groups, required2|Original assignment, required>.65|Pass all|
|---|---:|---:|---|
|20261031|1|.602000|no|
|20261107|1|.598020|no|
|20261108|1|.596040|no|
|20261109|1|.602000|no|
|20261110|1|.599602|no|

All cases lose the required chr12 distinction.20261108 also has normal FPR.1,
failing the unchanged strict<.1 criterion. The existing109 case never reaches
this veto because it already fails initial weighted confirmation. The guard
executes five times across the other workflow paths; all five BIC decisions
prefer the continuous-state model. Log-likelihood gains for the more complex
model are1.5629,1.2129,2.9399,3.1361,2.0772, insufficient for the fixed penalty.
The larger model fitting slightly better does not imply supported discreteness.

[Full original criteria](../dev/continuous_state_guard_20261109/recovery_v1/report.json),
[every guard decision](../dev/continuous_state_guard_20261109/recovery_v1/guard_records.json)
and the per-seed stage tables retain the complete evidence. The guard decision
list identifies source windows/cell rows; source matrices and case provenance
are bound by the enclosing per-seed execution and stage records rather than a
standalone case field in each guard row. No simulated false-split or power rate
for this guard has been established. Zero errors by abstaining everywhere would
not meet the simultaneous recovery objective.

## Execution, changes and next decision

Job23133178 completes0:0 on nslab CPU4/16GB/2h, actual1m13s process time and
1515920KiB peak. The validated noneditable wheel and immutable source_v1 execute
from a neutral directory; no arrays, shared installs or frozen-output changes.
[Command](../dev/continuous_state_guard_20261109/jobs.json),
[accounting](../dev/continuous_state_guard_20261109/accounting.tsv) and
[process resources](../dev/continuous_state_guard_20261109/resource_summary.json)
are recorded. No new full package suite/GPU parity was run for research-only code.

Research helpers add an optional candidate-filter callback to weighted stages
and a guard flag through the reference-noise/recovery runner; ordinary defaults
retain the previous path. Production/native/test sources are unchanged. Keep
the failed guard and its full recovery results; do not relax its BIC or original
scientific criteria to manufacture recovery. This negative result rejects this
particular approximation and rule as a useful correction. It is not proof that
every better RNA model is impossible, nor evidence for withdrawing independent
DNA-supported clones.
