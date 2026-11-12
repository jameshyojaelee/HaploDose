# Estimate proposal exposure outside the candidate genomic block

## Outcome

Outside-block exposure failed the prespecified recovery gate. Job 23185429
completed in 3m35s on node010. All 20 workflows and the input identity checks
completed. The exact-root, independent PMF, perturbation and permutation checks
passed. The sham reproduced both arrays and all 14 tables including logs in
each of the five cases. The scheduler recorded 1,089,972 KiB peak RSS and the
process recorded 1,313,548 KiB.

| Original seed | Native and exact all-bin | Outside block | S1 assignment accuracy, native → outside block | Normal false calls, native → outside block |
|---|---|---|---|---|
| 20261031 | Pass | Fail | 84.44% → 60.20% | 5 → 8 / 90 |
| 20261107 | Pass | Pass | 76.45% → 74.21% | 6 → 6 / 90 |
| 20261108 | Fail | Fail | 59.60% → 59.60% | 10 → 10 / 90 |
| 20261109 | Fail | Fail | 60.20% → 60.20% | 2 → 2 / 90 |
| 20261110 | Fail | Fail | 59.96% → 59.96% | 8 → 8 / 90 |

The exact all-bin fit preserves the original six criteria and S3 group
counts on every case. Outside-block exposure is exercised on all 1,258 window
calls in that arm. The sham and exact all-bin arms each make 1,406 calls.
Each arm records one unresolved window under the existing numerical checks.
The unresolved outside-block window spans bins 114–132 in seed 20261031,
outside the chr12 target. The fitted chr12 target window 210–219 has gain
22.954 under the native exposure, 22.954 under exact all-bin exposure and
29.146 under outside-block exposure. The proposed minority increases from
126 to 131 of 508. Native and exact all-bin workflows accept that split.
The outside-block workflow merges the proposed groups during the split test,
with no accepted S1 split. It loses the two-group, assignment and chr12
criteria. There is no additional S3 group in any case.

Excluding the region can strengthen the proposal likelihood while worsening
the downstream confirmation. The tested change rescues none of the original
three failures. Reject it as a recovery correction. No larger error/power
run, promotion, production edit or criterion change follows. Biological
sample size is zero. The 5% whole-specimen error criterion and uncertainty
calibration remain unresolved.

[Completed report](../dev/outside_block_exposure_20261110/probe_v1/report.json),
[sham controls](../dev/outside_block_exposure_20261110/probe_v1/SHAM_LOCKED.json),
[candidate window trace](../dev/outside_block_exposure_20261110/probe_v1/candidate_trace.json),
[candidate seed 20261031 decisions](../dev/outside_block_exposure_20261110/probe_v1/candidate/seed_20261031/fit/log.tsv),
[completed manifest](../dev/outside_block_exposure_20261110/probe_v1/run_manifest.json),
[job accounting](../dev/outside_block_exposure_20261110/accounting_complete.tsv).

## Qualified alternate node

Job 23185507 completed in 1m07s on node010 with nslab, 1 CPU, 8 GB RAM,
and a 30-minute limit. It reproduced the corrected wheel's frozen fixture
exactly on both arrays and all 14 tables including logs. The scheduler
recorded 674,732 KiB peak RSS and the process recorded 906,872 KiB. While
the main job was still pending, its required node was amended to node010.
Job ID 23185429, innovation account and QOS, 1 CPU, 8 GB RAM, 3-hour limit,
source, wheel, runtime and all gates are preserved.

[Exact fixture report](../dev/outside_block_exposure_20261110/node010_control_v1/report.json),
[control submission](../dev/outside_block_exposure_20261110/node010_control_job.json),
[placement amendment](../dev/outside_block_exposure_20261110/node010_placement.json).

## Original submission and design

Submitted job 23185429 on CPU node003, innovation account and QOS, 1 CPU,
8 GB RAM, 3 hours. The frozen source contains 955 files and preserves all
953 files in the parent manifest. The launcher and driver passed static
syntax checks. The installed wheel was previously validated on compute nodes
and node003 passed the exact fixture check. Inference results are pending.
At 17:18 Eastern the scheduler reported Priority and estimated a 23:13
Eastern start on the same node. The estimate is provisional. No dependent
job was submitted.

The native proposal fit uses the region being tested to estimate each cell's
exposure. A genuine regional RNA shift could contribute to that exposure and
reduce the proposed contrast. Test this on all five original recovery cases
with native, exact sham, exact all-bin and exact outside-block arms. The
all-bin arm isolates changes from solving the NB exposure score to convergence.
The outside-block arm removes the candidate's entire genomic reset block
from the exposure fit. Nuisance terms, proposal windows, confirmation,
merging, and all recovery criteria follow the existing workflow. When exposure
has no outside support in any cell, use the all-bin estimate where available
and the native estimate otherwise.

The exclusion applies to the local exposure fit conditional on existing
nuisance terms. A second alteration or correlated expression state outside
the block can still bias the fit. Exact numerical fixtures and exact sham
arrays and tables including logs must pass before interpreting candidate
outcomes. Advancement requires all five original assertions, every previously
passing individual criterion, and no increase in S3 groups. The exposed
cases are synthetic diagnostics. A positive result would need separate
complete-workflow false-split and power validation.

[Frozen design](../dev/outside_block_exposure_20261110/DESIGN.md),
[execution command](../dev/outside_block_exposure_20261110/jobs.json),
[source snapshot](../dev/outside_block_exposure_20261110/snapshot.json),
[research driver](../benchmarks/hardening/outside_block_exposure.py).
