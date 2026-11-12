# Trace the first fresh membership collapse before proposing another correction

COMPLETE: replay 23160939 exits 0:0 in 30m22s and reader 23161705 exits 0:0
in 25s. Both archived arms reproduce exactly: 4 numerical arrays, 28 tables
INCLUDING logs, and the candidate count-fold archive. All ten captured stage
boundaries are therefore eligible for interpretation. The observer did not
change the measured workflow. These are exposed mechanism results, not new
error/power estimates or biological validation.
Recorded process peaks are1087328KB for the replay and121700KB for the reader;
[resources](../dev/phase8_membership_trace_20261110/resources_completed.json)
and [accounting](../dev/phase8_membership_trace_20261110/accounting_completed.tsv)
retain both measurements and exit states.

| Boundary | Tumour profiles | Target contrast | All-cell correct assignments |
|---|---:|---:|---:|
| First two-profile EM |2|.3246221|445/600|
| Additional split, then EM |3|three pair contrasts .5380745/.2738516/.2642229|375/600|
| Statistical merge, with profile refitting |2|.0077755|377/600|
| Next EM |2|exactly zero; profiles identical|358/600|
| Exact-profile coalescence and final fit |1|no pair|359/600|

The statistical merge joins profiles 1/2, whose target effects are -.2096325
and +.3284420. Its comparison is untestable because one weighted effective
population is 29.6520494, below the unchanged minimum 30. The generic log text
"no depth overlap" is imprecise here. The other comparison involving the small
component has effective support 22.2404222. The two better-supported groups have
effective support 184.6024673/251.6596932 and a resolved segment. The merger's
profile refit also changes the untouched component's profile. All recorded
nuisance arrays remain identical across this merge boundary; the large
attenuation occurs during merging/profile refitting, before the next EM.

This establishes where the contrast disappears, not that the earlier profiles
are genomic clones or that retaining an unsupported group is valid. Removing
exact-profile coalescence, lowering 30 to accommodate 29.65, or retaining three
untestable distinctions would not solve the declared scientific task.

The positive guard's selected 477-cell population has discovery-membership AUC
.6897964 and joint-membership AUC .7454718 against simulator labels. These are
conditional descriptive values on that selected population, not calibrated
genomic probabilities or directly comparable to all600-cell AUCs. Every stage
and component correlation remains in the reader report. Rank AUC agrees exactly
with direct pair enumeration; cell pairs are not biological replication.

Next, test whether rejecting a failed addition after the unchanged full
refinement preserves the earlier supported model. This adds a search constraint
rather than relaxing merge support. The sham must preserve every old output
before the intervention can run. The other under-split/over-split cases, old
recovery gate and fresh false-split/power criteria remain required.
[New causal intervention](REFINEMENT_CHECKED_SPLIT_2026-11-10.md).

[Exact replay report](../dev/phase8_membership_trace_20261110/trace_v1/report.json),
[complete stage analysis](../dev/phase8_membership_trace_20261110/summary_v1/report.json),
[reader manifest](../dev/phase8_membership_trace_20261110/summary_v1/run_manifest.json).
Earlier preparation checkpoints below retain their original context.

Dependent reader23161705 is submitted afterok23160939:1CPU4GB30min,nslab.
It requires the COMPLETE exact replay before computing anything scientific.
All stage profile contrasts,membership masses,all600-cell MAP accuracy,nuisance
changes and pairwise likelihood rank separation are retained. Rank AUC is
independently checked by explicit cell-pair enumeration. Positive guard
memberships align to stages by exact global row IDs. No refit or threshold
selection occurs. The reader can locate a loss boundary,not prove its cause.
[Reader design](../dev/phase8_membership_trace_20261110/SUMMARY_DESIGN.md),
[executed reader](../benchmarks/hardening/phase8_membership_summary.py).
Its909-file snapshot preserves all907 replay source files exactly. Both jobs
are pending/running at this checkpoint; no parity or mechanism result claimed.

Implemented one unchanged full-workflow replay of exposed phase8replicate1,
seed202611098001. Its original candidate selects a long-event cap fallback,
then accepts a short-event guard at6.79947. Later statistical merging and
profile refitting produce identical assignment profiles,which correctly
coalesce. Removing exact-profile coalescence would leave RNA-indistinguishable
components and cannot establish clone recovery. Do not weaken that check.

This follows the COMPLETE rejected pilot:7/10exact candidate recovery versus
production9/10,state-null1/10versus0/10,all-cell assignment.699167versus.708167.
The fixed final100 certification is futile. Neither this replay nor a known
failure panel is new independent evaluation. Biological n=0. Replicate2also
under-splits while4over-splits,so a blanket stricter acceptance rule is not a
demonstrated correction.

The replay retains original production/candidate order,the exact generator,
method,wheel,node005,capturedruntime003 and4CPUallocation. An observer records
input/output profiles,memberships,mixture weights,nuisance state and log
positions around EM and merge calls. It starts AFTER the weighted-stage context
enters and preserves method code identities used by existing caller checks.
Temporary nested event profilers restore it. No extra likelihood evaluations
or parameter updates run during inference.

Before interpreting any captured state,the replay must reproduce ALL archived
numerical arrays,tables INCLUDING logs,and count-fold arrays EXACTLY for both
arms. Any mismatch fails the replay. No tolerance or table exclusion. New
states may locate the first transition that loses the contrast;they do not
prove a biological clone,calibration or a correction's efficacy. Subsequent
replicate2/4 challenges remain required;one-case success cannot qualify a fix.

Job23160939:CPU/nslab,4CPU32GB90h,original005,v6&amd,no array/interactive QOS.
Submitted with original read-only validated environment. Immutable source
snapshot contains907files,all905 original pilot files unchanged,plus the new
observer and design. Initial scheduler status is recorded in the task's job
journal;execution and scientific results remain pending. The current package's
allele/namespace repairs are retained in the worktree;this historical replay
uses its frozen method. No default,scientific criterion or frozen result edit.

[New observer/strict replay](../benchmarks/hardening/phase8_membership_trace.py),
[six saved decision chronologies](../dev/phase8_membership_trace_20261110/saved_chronology_checkpoint.json),
[locked design](../dev/phase8_membership_trace_20261110/DESIGN.md),
[snapshot derivation](../dev/phase8_membership_trace_20261110/snapshot_record.json),
[exact job command](../dev/phase8_membership_trace_20261110/jobs.json),
[original pilot and failure directions](POOLED_PREVALENCE_PILOT_2026-11-10.md).
