# Fixed-platform contrast for the failed clone replay

The unchanged full control reproduced all30 scientific criteria but failed
historical table equality. A bounded paired-node experiment now tests the
execution platform as a cause. This does not change the control or any recovery,
gene-veto, merge or false-split criterion. The depth-conditioned full candidate
remains unrun. Biological n=0; no biological or performance gain follows here.

## Specific causal question

The first differing log includes a near-zero event mixture gain2.44e-15, below
its existing3.88e-10 rejection budget. Some other rejected-mixture gains/budgets
also differ slightly. These observations do not prove that all scientific
outputs are unchanged; the exact full control still fails. Previous same-host
2/4-thread repeats matched529 numerical calls but did not isolate the historical
cause. Do not repeat that experiment unchanged.

Scheduler metadata shows the original node003 running Linux5.15 and recent
node010 running Linux6.8. Installed package RECORD manifests previously matched;
that does not establish which transitive system libraries were actually loaded.
The hypothesis is platform-dependent numerical behavior. This is unverified.
The relevant estimand is an exact numerical input/output difference at fixed
code/data/parameters, not DNA dosage or clone accuracy. Process history,
transitive libraries and a still-unrecorded execution difference are alternatives.

[Design fixed before results](../dev/recovery_platform_20261109/DESIGN.md).
Run the original seed108 observer on exactly the historical003 and recent010
nodes, with the same validated noneditable wheel/native identity, source,
input hashes, chunk97,2 requested CPUs, one BLAS thread, seed and hash seed0.
Capture actual mapped shared-library hashes, ELF interpreter, OS/hardware and
an archive of the loaded system libraries. No shared installation is modified.
The library archive is evidence for a possible later causal intervention, not
an installed replacement or a demonstrated fix.

[Implementation](../benchmarks/hardening/recovery_platform.py) reuses the existing
observer without inserting new fitter caller frames. The dependent comparison
will compare numerical traces and every retained array/table against BOTH the
historical and recent outputs. Different outputs on identical numerical inputs
would localize the discrepancy; a library difference alone would not prove its
causality. Agreement across hosts would reject this tested explanation. No node
search or selective passing control is allowed. Original failed gates remain.

## Execution checkpoint

| Job | Role | Account/QOS; CPUs / RAM / requested time | Status at checkpoint |
|---|---|---|---|
|23152722|fixed historical node003|innovation after amendment;2 /16GB /1h|pending, Resources then Priority|
|23152723|fixed recent node010|innovation;2 /16GB /1h|complete0:0,48s|
|23152724|dependent exact comparison|nslab;1 /4GB /1h|pending afterok of BOTH worker IDs|

The node003 job was initially submitted on nslab. Its projected wait exceeded a
day; amend only its pending account/QOS to already authorized innovation. This
does not free a node whose48 CPUs are now all allocated. The revised projected
wait remains long; preserve the actual pending job rather than duplicate it.
[Allocation amendment](../dev/recovery_platform_20261109/ALLOCATION_AMENDMENT.md).

A read-only SSH attempt to inspect node003 system-library metadata is rejected
by pam_slurm_adopt because there is no active allocation on that node. No
off-scheduler analysis or file copying occurs; do not bypass the scheduler.
The scheduled worker remains the route to that evidence.

The completed recent worker retains the actual library inventory/archive and
unchanged single-case trace in
[node010_v1](../dev/recovery_platform_20261109/node010_v1/).
The paired result is not yet available; no platform cause or correction is
claimed. Exact job commands and immutable source are in the
[task directory](../dev/recovery_platform_20261109/).
No production/native/original-test/dependency/frozen change. The SOTA objective,
three production recovery failures and5% specimen false-split/power objective
remain unachieved.

Later allocation checkpoint: innovation did not improve the queue estimate;
restore the original nslab account/QOS on23152722, with the same fixed node and
job ID. The table above and ledger checkpoint1 preserve the intervening
innovation state. Do not duplicate the job or repeatedly change accounts.

## Completed contrast and causal correction

Both original workers and their comparison complete. There are529 calls in
each trace,104 differences and no unmatched calls. The first differing event
mixture has EXACTLY equal O/E/r/step inputs but different numerical output;
one reported score-budget difference is1.01e-14. All scientific criteria agree.
Actual OS releases are Ubuntu22.04.4 on003 and Ubuntu24.04.3 on010. All seven
mapped system-library hashes differ, while no additional library basename
appears on either side. This is stronger than package-version matching, but an
association alone does not establish library causality.

The [causal design](../dev/recovery_platform_20261109/RUNTIME_DESIGN.md) fixes
recent node010, Python environment, native wheel, source, input hashes, seeds
and threads, and changes ONLY the dynamic loader/system-library bundle to the
seven captured003 files. They are copied into a new hashed task-local directory;
no shared installation or production dependency is changed. The Python
processes use that explicit loader/library path. Actual mapped library hashes
are checked through the captured process inventory.

Snapshot_v2 corrects only diagnostic comparison/argument handling in the
platform wrapper; inference and the numerical observer remain byte-identical.
The whole source-tree snapshot hashes differ and remain recorded. The causal
claim concerns the unchanged numerical computation, not identical documentation
or unused research files across snapshots.

This intervention reproduces ALL529 original-node numerical calls EXACTLY,
with zero input/output differences and zero mapped-library hash differences.
It supports the captured system runtime as a cause of the platform discrepancy,
not merely a different successful node. It does not identify which one of the
seven libraries is sufficient, nor imply a biological or power gain.
[Causal comparison](../dev/recovery_platform_20261109/compare_runtime_v2/report.json).

The initial comparison helper used np.array_equal, which incorrectly marks
matching NaNs in the posterior array unequal relative to the original test.
Correct it to the SAME np.testing.assert_array_equal used by the strict
historical control, with no tolerance. Preserve compare_v1 and write new
compare_corrected_v2/compare_runtime_v2 outputs. Under corrected exact semantics,
native010 matches all16 recent files. Two-thread003 and the runtime intervention
match each other but differ from the historical full run in three numeric
tables; the single-case traced2-thread run is not the original full4-thread
execution. Do not call that partial comparison a passing historical gate.

Finally run the ORIGINAL unmodified full control, all five cases and original4
threads, under the captured runtime on fixed010, without numerical observers.
Job23153243 completes0:0 in79s and passes EVERY original scientific criterion
record, exact array/table comparison and all15 fold comparisons. It reproduces
the historical4/5 result; this is computational reproducibility, not new clone
recovery. All earlier failed controls remain unchanged. The runtime PLUS original
full execution now supplies a passing prerequisite without weakening it.
[Full control report](../dev/recovery_platform_20261109/control_runtime003_v1/report.json),
[full-control design](../dev/recovery_platform_20261109/FULL_CONTROL_DESIGN.md).

| Job | Role | Account/QOS; CPUs / RAM / requested time | Final outcome |
|---|---|---|---|
|23152722|original003 worker|nslab after restored allocation;2 /16GB /1h|complete0:0,53s|
|23152723|native010 worker|innovation;2 /16GB /1h|complete0:0,48s|
|23152724|initial comparison|nslab;1 /4GB /1h|complete0:0,6s; NaN helper semantics later corrected|
|23153184|captured003 runtime on010|innovation;2 /16GB /1h|complete0:0,36s|
|23153185|corrected native-pair comparison|nslab;1 /4GB /1h|complete0:0,10s|
|23153186|causal runtime comparison|nslab;1 /4GB /1h|complete0:0,7s;529 calls exact|
|23153243|original strict full control|innovation;4 /16GB /3h|complete0:0,79s; all original checks pass|

This permits a separately snapshotted depth-conditioned full candidate ONLY
after its new mode's control path also passes. No production, original-test,
scientific criterion, frozen artifact or confidence claim is changed. Runtime
capture fixes execution provenance/replay; the original3 production recovery
failures and5% false-split/power requirement still need scientific improvements.
