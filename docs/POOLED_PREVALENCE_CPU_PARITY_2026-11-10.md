# Smaller CPU allocation rejected by the exact workflow gate

The later fixed-thread causal replay also fails the SAME exact trace gate.
Job23159255 retains1allocatedCPU but restores OMP/Rayon4 and BLAS1, recording
actual allocation and affinity[16,64]. All5 original assertions/30criterion
values match;108log remains222rows versus231baseline. Complete failed manifest,
exit1:0,2m53s. Holding method thread counts fixed does not restore parity;
the root cause is [Unverified]. No further resource variants are attempted.

This extra contrast was justified by NEW scheduler estimates of multi-day waits,
and separated the two factors changed together by the first experiment. It
uses the same strict gate, not a new tolerance or omitted table. No original
pilot job is canceled, replaced or changed. All six retain4CPU/32GB/90h and
the fixed scientific design. First-attempt findings below remain unchanged.
[Prespecified causal contrast](../dev/pooled_prevalence_cpu_parity_20261110/FIXED_THREAD_DESIGN.md),
[scheduler evidence](../dev/pooled_prevalence_cpu_parity_20261110/scheduler_before_four_threads.txt),
[fixed-thread decision](../dev/pooled_prevalence_cpu_parity_20261110/four_threads_report.json),
[actual environment/affinity](../dev/pooled_prevalence_cpu_parity_20261110/analysis_23159255/execution.json),
[second failed manifest](../dev/pooled_prevalence_cpu_parity_20261110/control_four_threads_v1/run_manifest.json),
[both attempts/accounting](../dev/pooled_prevalence_cpu_parity_20261110/accounting_both_attempts.tsv).

Job23159147 tested whether queued pilot workers could use1CPU instead of4
on the same005 node, captured003 runtime, frozen source and read-only wheel.
It failed the unchanged exact-table comparison after2m53s,exit1:0. No pilot
allocation was changed. The existing worker and queued jobs retain4CPUs,
32GB,90h and the original accounts, source, outputs and scientific design.

This was one bounded execution experiment, not a new biological model or a
fresh error/power look. The original four-CPU guarded recovery remains5/5
with all30 criteria passing. The fresh pilot remains separate and live.

## Decisive comparison and observed failure

The frozen `pooled_prevalence_workflow.py --control` reruns all five original
cases without the new safeguard, then compares with the earlier full unguarded
baseline using literal equality:30 criterion values,10arrays,70tables and15
whole count folds. The gate and source were not edited. Fewer allocated CPUs
also change the launcher's OMP/Rayon thread counts; BLAS remains one-threaded.

All five original recovery assertions pass and criterion values match. The
comparison then fails for seed20261108's `log.tsv`:231x21 baseline versus
222x21 replay. Small-file inspection identifies a changed sequence of
`event_mixture_numerically_ambiguous` records. The exact assertion stops the
comparison; do not claim that all remaining outputs or folds were checked.

The different allocation/thread environment is an execution difference, but
the causal source of the trace difference is [Unverified]. Passing the five
final recovery assertions does not prove equivalence of the adaptive search
under fresh nulls. Do not drop diagnostic rows, add a tolerance, redefine the
gate or assume the difference is harmless to permit the resource change.

The prespecified rule rejects this execution adaptation on any discrepancy.
No more CPU counts or nodes are screened. One worker remains running; five
wait for capacity on the required node. This is not a QOS memory-limit issue.
No interim pilot scientific outcomes are used in the decision.

## Reproducible evidence and preserved state

[Design fixed before submission](../dev/pooled_prevalence_cpu_parity_20261110/DESIGN.md),
[exact command/jobID](../dev/pooled_prevalence_cpu_parity_20261110/jobs.json),
[failed manifest](../dev/pooled_prevalence_cpu_parity_20261110/control_v1/run_manifest.json),
[full failure and original assertions](../dev/pooled_prevalence_cpu_parity_20261110/analysis_23159147/output.txt),
[process resources](../dev/pooled_prevalence_cpu_parity_20261110/analysis_23159147/resources.txt),
[initial source state](../dev/pooled_prevalence_cpu_parity_20261110/initial_git_status.txt).

[Structured decision](../dev/pooled_prevalence_cpu_parity_20261110/report.json),
[scheduler accounting](../dev/pooled_prevalence_cpu_parity_20261110/accounting.tsv),
[unchanged pilot allocations](../dev/pooled_prevalence_cpu_parity_20261110/pilot_allocations_unchanged.json),
[preservation](../dev/pooled_prevalence_cpu_parity_20261110/preservation.json).
Inference process peak RSS1,611,788KB. All734 scientific source files match
the delivered bulk snapshot, and all171 prior ledger values remain exact.

CPU/nslab,1CPU/8GB/1h, original005. Scientific code, installed wheel, original
tests, thresholds, primary family, horizon, frozen results and pilot job handles
remain unchanged. No new package suite is necessary or claimed for this
execution-only experiment. Biological n=0 and no new false-split, power,
calibration or genomic-membership conclusion.
