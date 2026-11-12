# Integrate the qualified partition average without changing old defaults

DELIVERY GATE PASSED:23176769 completes0:0 in19s scheduler time,9.82s measured,
peak165448KiB. ALL181 original ordered cases reproduce exact candidate profiles,
private fit arrays,selection diagnostics,previous gene-Laplace profiles and
diagnostics,and earlier default profiles/diagnostics. Both original source-count
examples also pass exactly. No missing case,alternative population or relaxed
tolerance is used. All four workers complete0:0 in18m28s,19m07s,18m42s,17m25s;
their measured process peaks are236228/228908/227548/223408KiB.
[Exhaustive gate](../dev/native_partition_average_20261110/summary_v1/report.json),
[complete manifest](../dev/native_partition_average_20261110/summary_v1/run_manifest.json),
[all exact commands](../dev/native_partition_average_20261110/jobs.json).

The public option is documented in README and schema2 migration notes ONLY after
this complete gate. Use gene_segments=True,average_segments=True on either
stacked_rna_profile or stacked_count_profile. Previous defaults remain. New
documentation changes no executed scientific code,so full scientific reruns
are not repeated for documentation. The validated wheel/source stays immutable.

Scope:the measured accuracy improvement is for POOLED QUERY profiles from
single-cell RNA datasets against relative bulk DNA. It does NOT establish
per-cell DNA accuracy,clone membership,bulk/spatial transfer or calibrated
uncertainty. Counts may be raw or explicitly reconstructed in this descriptive
API;the latter acquires no raw-count calibration. DNA-reference correspondence
and shared-reference uncertainty remain scientific limits. Field means lead
the existing16-method panel,while all30 formal superiority verdicts remain
untestable under its missing-family rule. No independent/SOTA or biological
false-split guarantee follows. Recovery23170760 and four exposed downstream
tests remain separate unfinished work. Earlier checkpoints retain context.

Full validation23176763 COMPLETE0:0 in3m08s scheduler time. ALL508 Python
tests and4 Rust tests pass,with1GPU-only skip. All10 warnings are the existing
SciPy forwarding of the HiGHS threads option;two arise in the new tests.
Measured pytest peak2131448KiB,scheduler1890452KiB. Candidate wheel SHA256
89d967879423b0fa0d522fa92479d93778e0e72509869ff094c846eaed07e080;
native binary retains e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7.
[Full Python checks](../dev/native_partition_average_20261110/validation_23176763/pytest.txt),
[Rust checks](../dev/native_partition_average_20261110/validation_23176763/cargo_test.txt),
[installed artifact](../dev/native_partition_average_20261110/validation_23176763/installed_artifact.json).

Original count-chain23176768 also COMPLETE0:0 in25s,measured22.57s,process
peak641036KiB. BOTH original identified matrices reproduce exact gene effects,
new profiles and selection diagnostics:stomach103query/2118reference cells/
10133genes;42MGBA517/2467/10116. Cell/gene/reference identity and order are
checked. These are reconstructed measurements and gain NO raw-count calibration.
The old metadata mixed-type warning and HiGHS warning remain visible;no missing
ID,profile or diagnostic is ignored. DNA is not read for reproduction.
[Count report](../dev/native_partition_average_20261110/count_parity_v1/report.json),
[complete count manifest](../dev/native_partition_average_20261110/count_parity_v1/run_manifest.json).

ALLfour installed181-case workers23176764–23176767 are RUNNING. Their first
completed cases match all predictions,fit arrays,diagnostics and old defaults
EXACTLY. This is partial coverage;final delivery gate23176769 still waits for
ALL181 original ordered cases. No full delivered-feature or scientific validation
claim follows until that gate passes. Earlier checkpoints below retain context.

Environment setup23176758 COMPLETE0:0 in3m48s on I/O002. Exact Conda/pip copy
and import checks pass;no project package is installed during setup. The private
environment reports NumPy2.4.6,SciPy1.17.1,JAX0.10.2 from the existing version
set. Build/full-suite23176763 is now RUNNING on CPU005 and has reached pytest.
No full test result or installed reproduction result yet. All environment
warnings remain visible;no safety override. Separate matched field comparison
now has leading point means but ZERO formal superiority verdicts because the
fixed comparator family remains incomplete;see the linked scientific evidence.
Earlier checkpoints below retain their original context.

Resource checkpoint:the SAME environment job23176758 now RUNS on I/O node
ne1di6-002. While it was still PENDING and before the prefix existed,its
partition/features were changed toio/amd for designated dependency-file staging.
Account/QOS,1CPU8GB1h,source and exact dependency versions are unchanged.
No restart or duplicate. Build/tests/scientific reproduction stay on CPU005;
the one-hour request is below the live72-hour I/O limit. Existing package-cache
checksum warnings remain visible;no safety-check override is applied.
[Exact before/after placement](../dev/native_partition_average_20261110/environment_io_resource_amendment.json).
All executable package checks still await this setup;no delivery result yet.
Earlier submission records below retain their original context.

Implementation and submission checkpoint,2026-11-10. The full181-line research
regression passed ALL8 fixed conditional comparisons,mean1Mb DNA Pearson
.713535568 versus current.711350891 andSM5 pseudobulk.667955383. ALL181 old
profiles/1448 old scores remain EXACT. This warrants an opt-in implementation,
not independent superiority,calibrated uncertainty or automatic mode transfer.
[Scientific evidence](PARTITION_AVERAGED_STACK_2026-11-10.md).

The package now implements the same averaging and stacking calculation through
explicit LOCAL interval/fold callbacks. It does not install a process profiler
or modify numerical functions at runtime. Only gene_segments.py and stacking.py
change among the941 preserved pre-edit files;average_segments.py is new. All
old tests and other previous files remain exact. Both existing APIs accept
gene_segments=True,average_segments=True. Nonboolean flags or averaging without
gene_segments=True fail. Existing defaults are unchanged. Public units remain
natural-log relative RNA;Gibbs weights on fitted losses are not integrated
Bayesian evidence,calibrated boundaries or genomic labels. No bulk,spatial,
clone decision or scientific threshold changes.

New independent tests enumerate raw-gene Fraction medians/losses over all8
partitions and compare direct scalar sums,uniform rational weights,singleton,
shift/scale invariance,masked-target exclusion,resets,count-API forwarding,
default preservation and invalid flags. Exact reset invariance uses dyadic
values so the translated observations themselves remain exactly representable.
No convenient numerical epsilon or relaxed scientific test was introduced.
AST and launcher syntax pass;executable checks remain pending at this checkpoint.

The installed-wheel gate runs ALL previous Python and Rust tests plus the new
tests. Four workers then require exact public predictions,private fit arrays,
all common diagnostics and old default/gene-Laplace profiles/diagnostics on
ALL181 cases. A separate job reconstructs BOTH locked examples from the original
identified count matrices. A final gate requires complete ordered coverage and
the successful count report. These reproduction jobs read no DNA. No delivery
claim follows merely from the implementation,build or successful process exit.

| Job | Role | Account/QOS | CPU / RAM / time | Prerequisite |
|---|---|---|---|---|
|23176758|Exact private Conda/pip versions|innovation|1 /8GB /1h|none|
|23176763|Noneditable wheel,full Python/Rust|nslab|2 /32GB /1h|23176758|
|23176764|Installed cases0:46|innovation|1 /8GB /2h|23176763 plus actual successful tests|
|23176765|Installed cases46:92|nslab|1 /8GB /2h|same|
|23176766|Installed cases92:138|innovation|1 /8GB /2h|same|
|23176767|Installed cases138:181|nslab|1 /8GB /2h|same|
|23176768|Original count matrices,two cases|innovation|1 /16GB /2h|same|
|23176769|Complete installed/count coverage|nslab|1 /4GB /1h|all five reproduction jobs plus actual reports|

While environment23176758 was still PENDING and before any environment prefix
or outcome existed,its requested node list was cleared. The SAME job may copy
dependencies and check imports on any compatiblev6/AMD CPU node. All numerical
build/tests and installed reproduction remain pinned005. Source,dependency
versions,account/QOS,CPUs,memory,time and scientific gates are unchanged;no
restart or duplicate. [Before/after scheduler record](../dev/native_partition_average_20261110/environment_resource_amendment.json).

All jobs use CPU/v6&amd,no arrays or interactive QOS. The949-file
source snapshot is immutable. Dependencies are copied at exactly existing
versions into a new private Conda prefix;no new production dependency or shared
installation. The submission logger stopped AFTER safely recording the accepted
environment job;remaining submissions resumed from its ID without duplicating
or restarting it. No scientific execution/source/resource changed.

[Frozen integration design](../dev/native_partition_average_20261110/DESIGN.md),
[package averaging](../python/haplodose/average_segments.py),
[independent tests](../tests/test_average_segments.py),
[installed reproduction](../benchmarks/hardening/partition_average_native_parity.py),
[count reproduction](../benchmarks/hardening/partition_average_count_parity.py),
[exact commands](../dev/native_partition_average_20261110/jobs.json),
[snapshot preservation](../dev/native_partition_average_20261110/snapshot_v1.json),
[static checks](../dev/native_partition_average_20261110/static_patch_check.json),
[submission continuation](../dev/native_partition_average_20261110/SUBMISSION_NOTE.md).

Open scientific limits remain:former evaluation cohorts are regression-only;
biological error calibration/replication and independent truth are absent.
The candidate has no established bulk/spatial transfer,clone recovery benefit,
or whole-specimen5% false-split guarantee. Recovery correction23170760 and its
four exposed downstream cases remain separate,unchanged tests.
