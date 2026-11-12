# Reference-composition proposals regress recovery and are rejected

Both reference-pooled and depth-interpolated gene proposal candidates recover
3/5, versus4/5 for the retained research control. They lose110 and fail109.
Reject both at the unchanged5/5 continuation gate; their larger null/power
experiment is not run. Production remains2/5 and the5% specimen error/power
criterion remains unmet. A node-matched replay passes the complete exact control;
the earlier failed replay and cancelled candidate jobs remain preserved below.

## Intended correction and identifying assumptions

Previous gene proposals estimate one within-bin composition from pooled query
proposal counts. The new model would allocate each native expected bin count to
genes using fit-reference composition, either pooled per type or interpolated
over the existing five reference depth strata. Native fitted type contributions
mix those gene fractions. This targets heterogeneous gene composition within bins,
not a new posterior threshold or a different merge criterion.

The estimand is distinct regional RNA populations under conditional gene NB
working models. Reference composition must transport to query cells. Gene-specific
response, reference programmes, mapping artifacts and continuous expression state
are alternatives. Corresponding independent DNA and suitable biological replication
would distinguish a genomic interpretation. Biological specimens added=0.
Native bin expectations and reference mixture still depend on existing full-count
fitting; reference-only gene fractions do not make the complete method independent
confirmation or biologically calibrated.

A half-count per gene/type is a declared positive mean regularizer, not observed
data or a calibrated NB/Dirichlet uncertainty model. It is distributed over strata
by reference exposure, so the total prior mass matches the pooled arm. Within-bin
fractions interpolate linearly between exposure-weighted log-depth knots. Outside
those knots the endpoint composition is held fixed and tails are recorded;
transport beyond reference support is unestablished. No query counts, holdout
reference counts, true noise, class labels or DNA fit these fractions.

New code is `reference_gene_composition.py`, `reference_composition_recovery.py`
and `reference_composition_numerics.py`. The existing research
`gene_event_proposal.py` gains a default-off region-provider argument. No package,
native, original test or dependency changes. Immutable source_v1 is retained.

## Initial numerical gate and failed control

Numerical23146198 passes independent Fraction mixed-type proportions, scalar
depth-stratum proportions/knots, interpolation/endpoints, reference row/gene
permutation and explicit failure for unsupported strata. The operation-based
rounding allowance is3.55e-14. This verifies those calculations, not transport or
recovery. Candidate/reference tables and all numerical fixtures remain saved.

Control23146199 completes every original inference run, with the same4/5 result
as the preceding gene-proposal research procedure, then fails exact comparison
of `log.tsv`. The failure is retained verbatim. Gains around1e-14 and their
approximately1e-8 numerical budgets differ. The native workflow already labels
these entries `event_mixture_numerically_ambiguous`; some such records appear or
disappear when the computed near-zero gain changes sign. This is not evidence
that an accepted scientific distinction changed.

Static byte inspection establishes50/50 exact matches: ten arrays/tables per
seed, including raw binned counts, posteriors, cell assignments, clone profiles,
segments, differences, specimen summaries, bins and genes. After removing only
the explicitly numerically-ambiguous action, all other log records are identical
in every seed. That characterization does **not** replace or pass the frozen
exact-log test. The source of the floating-point diagnostic variation has not
been isolated. No particular hardware or library cause is asserted.

Historical/new log row counts are245/245,224/222,260/260,167/165 and248/249 in
seed order31,107,108,109,110. All30 scientific criterion values remain unchanged;
109 still fails. These are regression cases, not false-split or power estimates.

## Initial jobs and decision, before node-matched continuation

| Job | Role | Resources | Outcome |
|---|---|---|---|
|23146198|independent numerical gate|nslab CPU1/4GB/1h|complete0:0,16s|
|23146199|full exact control|innovation CPU4/16GB/3h|failed1:0,3m51s|
|23146200|reference-pooled candidate|nslab CPU4/16GB/3h|cancelled before execution after failed dependency|
|23146201|reference-depth candidate|innovation CPU4/16GB/3h|cancelled before execution after failed dependency|

All attempted execution uses the isolated noneditable validated wheel and neutral
compute directory. Its prior412 Python/4 Rust passes and1 GPU-only skip remain
separate evidence; full suites were not rerun for research-only changes. The two
candidate recovery rates and their null/power behavior were **not estimated**.
No numerical tolerance, scientific criterion, frozen result or failed control is
altered to obtain a passing verdict. The full goal remains unresolved.

[Frozen design](../dev/reference_gene_composition_20261109/DESIGN.md),
[failed exact-control execution](../dev/reference_gene_composition_20261109/analysis_23146199/output.txt),
[bounded static failure characterization](../dev/reference_gene_composition_20261109/control_failure_characterization.json).

## Single node-matched continuation and decisive recovery failure

The original gene control ran on ne1dc6-006, while the failed replay ran on
ne1dc6-010. The [amendment](../dev/reference_gene_composition_20261109/NODE_MATCHED_REPLAY.md)
prespecifies ONE replay on the original node, with identical frozen source_v1,
environment, exact log/array/table/fold checks and all30 scientific criteria.
If it fails, dependent candidates stop again. It changes no tolerance or test.
This control23148129 passes in full. The earlier failure remains a failure;
passing on one node does not isolate the floating-point cause or prove portability.

The two already frozen interventions then execute on the same node. Their
complete original criterion records, saved source counts, references, fractions,
expectations, event likelihood inputs and stage traces remain under pooled_v2
and depth_v2. Both fail the5/5 prerequisite and neither is promoted.

| Seed | Exact control accuracy | Pooled-reference accuracy | Depth-reference accuracy | Pooled / depth original recovery |
|---|---:|---:|---:|---|
|20261031|.836032|.836032|.838057|pass / pass|
|20261107|.738095|.738095|.738095|pass / pass|
|20261108|.756917|.756917|.756917|pass / pass|
|20261109|.602000|.602000|.602000|fail / fail|
|20261110|.811623|.599602|.599602|fail / fail|

Accuracy uses the original tumour-call-conditional denominator and does not
replace the independent two-group, normal-FPR, tumour-TPR, one-group-S2 and
chr12-only difference criteria. All original rules are retained, including
failed exactly-two-group and chr12-support criteria for109/110.

The actual first-half chr12 proposal209–217 still reaches confirmation. For109,
all three paths return zero initial weighted-HMM segments. For110, control has
one initial and one matched supported segment, accepted score4.648711. Pooled
reference retains one initial segment but zero matched segments. Depth reference
has zero initial segments. The loss therefore precedes the final merge; a
merge-only intervention cannot recover these particular candidate failures.
Mean membership margins increase from .579744 control to .596473 pooled and
.607645 depth for110 while recovery worsens. Increased fitted certainty is not
improved or calibrated membership. These records are extracted without fitting
or rescoring from executed stage_trace.json files in
[failed_contrast_stage_records.json](../dev/reference_gene_composition_20261109/failed_contrast_stage_records.json).

| Job | Role | Account/QOS, CPU resources | Actual outcome | Batch peak RSS |
|---|---|---|---|---:|
|23148129|exact original-node control|innovation,4CPU/16GB/3h|complete0:0,4m12s|1,586,836KiB|
|23148130|pooled-reference intervention|nslab,4CPU/16GB/3h|complete0:0,3m22s;3/5|1,737,976KiB|
|23148131|depth-reference intervention|innovation,4CPU/16GB/3h|complete0:0,3m14s;3/5|1,797,148KiB|

Across all seven attempts: four completed jobs, one failed control, two cancelled
before execution. Larger fresh null/power runs are **not run** because recovery
fails, rather than because error criteria changed. No new biological specimens,
production clone change, frozen-output edit or scientific criterion relaxation
follows. Preserve both negative interventions for future mechanism decisions.
