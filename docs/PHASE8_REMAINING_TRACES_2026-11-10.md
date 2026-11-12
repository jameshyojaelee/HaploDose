# Trace the remaining under-split and over-split mechanisms

Replicate4 also COMPLETE: replay23166025 exits0:0 in34m03s and reader23166026
exits0:0 in24s. Both arms pass all4array/28table/log/fold exact gates. Eight
EM/merge boundaries show a different defect from replicate2. The initial
two-profile fit has422/600 correct assignments and effects-.083689/+.278354;
the tested pair has effective support366.125/63.860 and one resolved segment.
After another addition, three profiles remain; one has no full-count MAP cells
despite responsibility mass106.377. The terminal refinement drops an earlier
small component at27.3269, but the *new final* responsibilities retain another
component with mass29.91955057367891, below the existing30 minimum and with no
MAP cells. Its target effect+.230947 is close to, but not exactly equal to,
another effect+.223523. The final merge checks only the MAP-eligible1/3pair,
whose effective support349.280/97.499 passes; component2 is never tested there.
Final one-to-one assignments are438/600 with three reported profiles.

This establishes a final-state support gap, not permission to merge close
profiles using an arbitrary tolerance. A correction now enforces the SAME
responsibility-mass minimum after final scoring and recomputes aligned scores/
pools; its old-wheel counterexamples, complete build/tests and actual case4
replay are scheduled separately. No recovery improvement is yet claimed from
that repair. Do not infer that the new refinement-addition constraint fixes4:
its later search can retain extra groups and this support defect has its own
cause. [Final-support repair](FINAL_MEMBERSHIP_SUPPORT_2026-11-10.md),
[complete case4 trace](../dev/phase8_rep4_trace_20261110/trace_v1/report.json),
[complete case4 reader](../dev/phase8_rep4_trace_20261110/summary_v1/report.json).

Replicate2 is COMPLETE: replay23166023 exits0:0 in25m03s and
reader23166024 exits0:0 in25s. Both archived arms match EXACTLY on all4arrays,
28tables including logs and1fold archive. The reader independently checks
rank separation against direct pair enumeration. Replicate4 remains running;
do not interpret its incomplete capture as a validated trace.

The first accepted replicate2 split is the long-event fallback
`event_chrom12_bins210-228`, with338/262 proposal labels. The subsequent
full-count fit has responsibility masses542.4304/56.9210, target-region
effects-.0486907/+.3339804 and404/600 correct one-to-one MAP assignments.
The smaller component's responsibility ranks still separate the imposed
populations (descriptive AUC.807257), but this is not accurate membership.
Proposal-fold assignments used by the merge become565/35. After matching and
classification-margin weights, effective support is366.4130/25.0656, below
the unchanged30-cell minimum. The pair is untestable and merges to one group.
Nuisance arrays do not change during that merge. The second attempted split
repeats the mechanism, with effective support366.4165/25.0666.

Neither discovery/confirmation guard is positive in this case; the accepted
long proposal uses the previously declared exact-state-cap fallback. Thus this
case differs from replicate1's destructive later addition. The proposed
refinement constraint cannot create missing support or turn a negative guard
into evidence. Do not lower30, preserve an untestable distinction, or claim
that returning two profiles before merging fixes recovery. Longer refinement
was already rejected on the original panel in
[the membership study](MEMBERSHIP_2026-11-07.md); no new full-budget sweep is
justified by this trace alone.

Computational status: exact reproduction and stage correspondence pass.
Statistical status: exposed conditional mechanism only, with fitted nuisance
and memberships; no new error/power estimate. Biological status: no specimens
or independent DNA added. The competing explanation is an extreme-tail RNA
partition rather than a recovered CNA population. A correction must preserve
the unchanged final recovery criteria and improve the joint false-split/power
result; descriptive AUC or successful execution alone cannot establish it.

[Complete replicate2 replay](../dev/phase8_rep2_trace_20261110/trace_v1/report.json),
[complete stage reader](../dev/phase8_rep2_trace_20261110/summary_v1/report.json),
[executed merge weights](../dev/phase8_rep2_trace_20261110/trace_v1/replay/seed_202611098002/pooled_prevalence/trace.json).
Earlier preparation checkpoints below retain their original context.

Two further locked phase8 cases now run with the same observation-only capture
and strict output gate used successfully on replicate1. Replicate2 returns one
group; replicate4 returns three. These are exposed mechanism investigations,
not new error/power estimates or biological validation. The method, generator,
wheel, original production/candidate order and scientific criteria are fixed.

| Case | Replay | Account/QOS | Reader | Allocation |
|---|---|---|---|---|
|Replicate2,seed202611098002|23166023|innovation|23166024|replay4CPU32GB90h;reader1CPU4GB30min|
|Replicate4,seed202611098004|23166025|nslab|23166026|replay4CPU32GB90h;reader1CPU4GB30min|

All jobs use CPU; replay jobs retain original005/4CPU/capturedruntime003.
Readers depend afterok on their own replay. No array,interactive QOS,GPU,new
dependency or shared installation. The separate refinement sham/candidate chain
23164505/23164506 continues unchanged. Results remain pending.

The new driver parameters select only the already locked replicate and seed.
Every archived numerical array, table INCLUDING logs, and gene-fold archive
must match EXACTLY before captured states are interpreted. A mismatch fails the
replay; no tolerance or artifact exclusion. The reader preserves all stage
contrasts,memberships,nuisance changes,likelihood rank separation and positively
guarded membership correspondence with exact row IDs. Its terminal group count
must match the archived one/three, respectively.

For replicate2, locate the first loss of a candidate contrast and whether the
unchanged30-effective-cell support gate is consequential. For replicate4,
locate how the additional profile survives support/confirmation and whether it
partitions imposed expression state. These questions do not assume either
case is corrected by the proposed refinement check. No blanket cutoff change
is justified by having both under-splitting and over-splitting.

Static preparation correction: source_v1 initially inherited the replicate1
observer's requirement that at least one outer merge actually occur. The
archived over-split4 log has no such merge. Source_v2 keeps every original
replicate1/2 check and requires EM and merge-check boundaries for all cases,
without requiring a nonexistent merge in4. All exact output gates remain.
The unexecuted source_v1 is retained; no job ran it and no failed scientific
check was discarded. Each source_v2 preserves910/911 source_v1 files exactly;
the earlier base method/package files remain identical throughout.

[Parameterized observer](../benchmarks/hardening/phase8_membership_trace.py),
[reader](../benchmarks/hardening/phase8_membership_summary.py),
[replicate2 design/jobs](../dev/phase8_rep2_trace_20261110/jobs.json),
[replicate4 design/jobs](../dev/phase8_rep4_trace_20261110/jobs.json),
[first completed causal trace](PHASE8_MEMBERSHIP_TRACE_2026-11-10.md),
[separate refinement intervention](REFINEMENT_CHECKED_SPLIT_2026-11-10.md).
