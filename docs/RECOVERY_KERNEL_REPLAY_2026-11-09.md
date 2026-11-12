# Fixed-condition replay does not isolate the historical numerical difference

Four fresh seed108 research fits pass the unchanged recovery assertion and have
identical recorded inputs/outputs across529 calls per worker. This includes370
event mixtures,105 pooled segment calls,31 specimen profile calls,21 cell
likelihood calls, reference estimation and programme learning. The2/4-thread
comparison does not reproduce a numerical divergence. Final scientific files
match the recent failed control exactly. The historical discrepancy is still
unisolated; do not invent a numerical correction or pass its failed exact gate.
Production recovery remains2/5, and the state-null false-split/power failure remains.

## Question, observations and identification limits

The adjacent-depth experiment stopped at a failed full exact control: all30
scientific criteria agree with the old4/5 research result, but logs and three
seed108 profile/segment tables differ slightly. The aim here is to find the first
numerical operation with identical inputs and differing outputs, not to accept
more clones or weaken that comparison. Biological n=0. The observed quantities
are exact numerical function inputs/outputs; genomic identifiability and
biological validation are separate questions.

The [design](../dev/recovery_kernel_replay_20261109/DESIGN.md) fixes seed20261108,
the old quintile reference estimator, gene count splitting, weighted confirmation
and merging, effective-gene veto, chunk97 and all original recovery assertions.
Only observers are added. Wrappers copy exact inputs before computation and
outputs after computation. Content-addressed arrays preserve dtype, shape and
values without rounding. Event results are captured through the existing
profiling callback, preserving the fitter's direct caller. A source or a declared
package version alone cannot establish identical execution.

The four fresh workers use the requested thread sequence2,4,2,4 and fixed
PYTHONHASHSEED=0. BLAS remains1. Comparisons include two same-thread repetitions
and two cross-thread pairs. Each has529 matched calls, no missing calls and zero
different inputs or outputs. Each worker retains3,851 unique arrays totalling
44,789,299 uncompressed bytes. Exact-input replay is not an uncertainty interval,
power experiment or independent biological replication.

Every worker passes the original seed108 assertion. All16 compared final arrays/
tables, including log.tsv, equal the recent failed-control seed108 files byte for
byte. Relative to the historical depth-quintile fit, the same three numeric
tables and log still differ. The observer has not forced a new scientific result
or made the old exact gate pass. This bounded negative result argues against the
tested requested-thread/fresh-process mechanism on this host. It does not establish
portability across CPUs, environments, process histories or unspecified settings.
Stop this probe without an unsupported production change.

## Correct historical node identification

The original depth-quintile both_v1 output is job23131740 on ne1dc6-003 with4
CPUs. The later experiment initially requested006 because a different preceding
control had succeeded there, and actually executed010 with2 CPUs. The new design's
shorthand "006 to010" is superseded by the
[execution-record clarification](../dev/recovery_kernel_replay_20261109/EXECUTION_NODE_CLARIFICATION.md).
The original designs and records are preserved. Neither node identity nor the
present negative thread result isolates the original cause.

All21 historical inference Python modules are byte-identical to the current
source; only bulk.py and cli.py have changed. The native binary SHA is identical
in the old and current execution records. NumPy, SciPy, JAX and JAXlib installed
RECORD manifests also agree in every hashed entry (1043/749/628/132 entries).
That is a comparison of declared installed-file inventories; it is not a new
hash verification of every transitive loaded binary. No environmental cause is
asserted from package versions or inventory agreement.

## Every failed observer attempt and allocation amendment

| Job | Snapshot / role | Requested or final allocation | Outcome |
|---|---|---|---|
|23150305|source_v1, first observer|nslab CPU4/16GB/1h|failed before first specimen fit,41s|
|23150339|source_v2, path correction|nslab CPU2/16GB/1h after pending amendment|failed at first specimen event capture|
|23150379|source_v3, profiler-compatible observer|nslab CPU2/16GB/1h|all four workers and exact comparisons complete0:0|

First failure: the observer wrote worker_N/fit, omitting the seed_<identifier>
output component required by the existing gene-thinning hook. Source_v2 corrects
only that path, retaining the failure and partial output.
[Path correction](../dev/recovery_kernel_replay_20261109/RUNNER_PATH_CORRECTION.md).

Second failure: wrapping _fit_event_mixture inserted a caller frame, so the existing
weighted-stage profiler refused to capture its membership probabilities. This
guard correctly prevented a silent hard-membership fallback. Source_v3 observes
through a composed profiling callback, calling the original callback first and
restoring profile state in finally. It changes no fitter or research stage.
Partial trace saving is added for future exceptions.
[Compatibility correction](../dev/recovery_kernel_replay_20261109/PROFILER_COMPATIBILITY_CORRECTION.md).

Before source_v2 executes, its pending4-CPU allocation has an estimated four-hour
wait. Reduce only the allocation to2CPU/16GB/1h, preserving explicit2,4,2,4 worker
threads. The four-thread workers share the assigned cores; this is not a speedup
benchmark or original four-physical-core replay. Actual affinity is the same four
logical CPUs42,43,90,91 for every successful worker, on AMD EPYC9254. OpenBLAS
reports one thread/Cooperlake kernel, and OpenMP reports the requested2/4 threads.
No node is selected or retried to find a passing scientific result.
[Pending amendment](../dev/recovery_kernel_replay_20261109/TWO_CORE_ALLOCATION.md).

## Evidence and preservation

[Complete comparisons](../dev/recovery_kernel_replay_20261109/probe_v3/report.json),
[final file identity](../dev/recovery_kernel_replay_20261109/final_file_identity.json),
[historical Python identity](../dev/recovery_kernel_replay_20261109/historical_python_identity.json),
[installed RECORD comparison](../dev/recovery_kernel_replay_20261109/environment_record_comparison.json).
Each worker saves arrays.npz, trace.json, report.json, actual hooks, unchanged
criteria, execution backend and per-case outputs. All three jobs, exact commands,
failures, environment/wheel/native identity, process memory/runtime and snapshot
integrity are retained under dev/recovery_kernel_replay_20261109.

The same validated noneditable wheel is used read-only from a neutral directory.
The earlier421 Python/4 Rust passes and one GPU-only skip are separate evidence;
full suites and GPU parity are not rerun for these research-only observers. No
package, native, original test, dependency, shared installation or frozen result
changes. Computational status: bounded exact replay passes, historical control
still fails. Statistical status: no new identifiability or calibration evidence.
Biological status: no new corresponding DNA or independent null specimens. The
adjacent-depth candidates remain cancelled and untested.
