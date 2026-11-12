# Transfer pooled RNA profile accuracy into assignment refinement

**Completed; transfer rejected.** Job23185102exits0:0 in10m10s on003. All15full
fits finish. Independent projection checks pass,all five sham cases exactly
reproduce every baseline array/table including logs,and all actual gene-pool/
tie/order checks pass. Recovery drops2/5→1/5,and normal false calls increase
in every case. No larger null/power study or package integration. The delivered
pooled-profile result remains distinct and unchanged.

| Original seed | Complete recovery baseline / candidate | S1 tumour groups | Normal false calls /90 | Assignment accuracy baseline / candidate |
|---|---|---|---|---|
|20261031|pass / fail|2→1|5→16|.844444→.600394|
|20261107|pass / pass|2→2|6→7|.764471→.760396|
|20261108|fail / fail|1→1|10→13|.596040→.598425|
|20261109|fail / fail|1→1|2→4|.602000→.598020|
|20261110|fail / fail|1→1|8→10|.599602→.598020|

The last column retains the original tumour-call-conditional denominator and
does not substitute for the separate two-group/normal-error/chr12 checks.
Seed31loses normal FPR,exactly-two groups,assignment and chr12-only difference
criteria. Seed110also loses the previously passing normal-FPR criterion. No
additional S3group is produced. All original criteria and failed results remain.

The seed31failure occurs before any accepted split,not as a collapse of two
previously accepted candidate groups. The first three candidate profile updates
each have just one tumour component plus the reference-like component. The
first replacement changes dose by up to.127451natural-log units. Initial EM
effective masses change from about98/502to78/522,with the same three rounds.
Baseline later accepts chr12bins210–219(score6.3437041);candidate accepts no S1
split. Thus the single-tumour profile replacement changes the population before
any multi-tumour isotonic contrast can act. This localizes the observed failure;
it does not attribute it uniquely to a particular gene,programme or depth term.

All67sham and64candidate profile updates exercise the new estimator,without
zero-pool fallback. Maximum gene/bin pooling discrepancy2.73e-11is within
derived arithmetic budgets. The independent exact optima include[-7/20,-7/20,
-7/20,53/20]and[-13/20,-13/20,-13/20,1/10]. SciPy's existing HiGHS option-forwarding
warning remains in logs. Scheduler peak1220520KiB;process peak1478588KiB,measured
10m02.79s. Biological n=0;no new package/GPU suite or genomic calibration claim.

[All original criteria](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/transfer_v1/report.json),
[independent projection checks](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/transfer_v1/numerics.json),
[full exact sham check](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/transfer_v1/SHAM_LOCKED.json),
[all candidate dose/refinement traces](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/transfer_v1/candidate_profile_trace.json),
[complete manifest](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/transfer_v1/run_manifest.json).
Earlier submission checkpoint follows unchanged.

Submitted23185102:CPU003/innovation account and QOS,1CPU16GB6h,one job/no array.
ALLfive original recovery cases run under baseline,exact sham and candidate.
Independent projection checks precede fitting;all baseline/sham arrays,tables
and logs must match exactly BEFORE candidate execution. No outcome or performance
gain is claimed at submission. Public package/defaults/frozen results remain.

The delivered gene-Laplace partition average improves pooled profile regression,
but membership refinement still uses the native HMM/shared assignment effects.
This experiment tests a direct inference transfer,not another confidence rescaling
or a relaxed support threshold. It does not assume the pooled181-line result
transfers automatically;the earlier individual-cell transfer failure remains.

At each shared-profile refinement call,use its actual R to pool full gene counts,
then the SAME native reference-gene allocations and native depth/programme-
adjusted expectations to construct centered RNA effects. Smooth each tumour and
the pooled tumour with the unchanged validated public partition-average stack.
Project the proposed profiles onto the native pairwise order and exact ties in
every bin,with exposure-weighted mean equal to the new pooled tumour profile.
Weighted isotonic regression can pool weak contrasts but cannot separate an
existing tied pair or reverse its native direction. Record all actual updates.

Keep native HMM/post/path/own-profile eligibility,all proposal and confirmation
tests,gene vetoes,merges,minimum support,posterior/core/stopping and all original
assertions. These remain a hybrid RNA working procedure,not exact EM or a genomic
membership model. The new fractional gene pools are descriptive weighted count
summaries;they do not supply independent count calibration. Reference transport,
transcriptional state and initial grouping errors remain competing explanations.
Same-cell DNA and appropriate biological replication establish genomic claims.

Continuation requires5/5complete original recoveries,no lost passing individual
criterion and no additional S3group. Otherwise reject this unchanged transfer
as sufficient and stop expansion. A pass only permits a separately frozen
complete-workflow null/power study,not promotion. The simultaneous5%specimen
false-split criterion remains unchanged and unestablished. Biological n=0.

Independent Fraction enumeration over every order-constrained partition checks
PAVA,along with permutation,exposure scaling,exact ties and ordering. All real
gene pools must aggregate back to native bin pools within operation-derived
rounding budgets. Save R,counts,E,reference fractions,IDs,maps,old/new doses and
every smoother diagnostic. The sham performs all new computations but returns
native profiles,testing that instrumentation itself does not change predictions.

The952-file source preserves all949files of the delivered profile snapshot and
adds this driver/design plus the exact existing output-comparison helper. Native
clones.py/pipeline.py are byte-identical to the support-corrected source. The
read-only isolated508-test native-profile wheel is bound to the frozen source;
all Python/native/wheel hashes and actual CPU backend are asserted from a neutral
compute cwd. No new dependency or shared installation. Full package/GPU suites
are not rerun for research-only code;new projection and inference checks run
on compute nodes. AST and launcher syntax pass.

[Frozen design](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/DESIGN.md),
[research implementation](/gpfs/commons/home/jameslee/Haplodose/benchmarks/hardening/gene_profile_membership.py),
[exact job command](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/jobs.json),
[source preservation](/gpfs/commons/home/jameslee/Haplodose/dev/gene_profile_membership_20261110/snapshot.json).
