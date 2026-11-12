# Test segmentation averaging against the current profile estimator

Matched field23178680 COMPLETE0:0 in1m48s scheduler time,measured1m38.54s,
process peak334564KiB,scheduler241220KiB. It ran on CPUne1dc6-003 after the
pre-outcome placement amendment. ALL5068 preceding scores reproduce EXACT,
including missing outcomes,with all recorded input hashes unchanged. All181
projected RNA outputs lock before DNA. This establishes the executed matched
comparison across these two compatible CPU nodes,not general numerical parity.

| Method | Mean1Mb Pearson | Mean arm Spearman |
|---|---:|---:|
|New averaged-expert stack|.7093163399890323|.8508433045079116|
|Previous gene-Laplace stack|.707030834261393|.850106669252588|
|Previous RNA stack|.702969969234286|.8474459986908777|
|SM5 pseudobulk|.6621648151144185|.8357367749116862|
|MDL pseudobulk|.6649265943612028|.8286443816675824|
|Gene-level pseudobulk|.35155457960290226|.8491078571601732|
|inferCNV R|.6684765572601751|.8132941229730033|
|infercnvpy|.6252941919638846|.7980361782975964|
|CopyKAT|.6094027636176813|.7369578393272872|
|SCEVAN|.5682890161572478|.684472620382208|

Point means use the SAME167 complete-family lines across ALL16 methods. The
new primary has the highest point mean on both endpoints. However CopyKAT13
andSCEVAN14 missing populations remain in the full181 score table,so ALL30
formal comparisons are untestable by the original strict family rule. ZERO
superiority verdicts pass. Do not call the ranking independent SOTA validation,
drop methods,or replace the failed formal gate with these point means. The
new standalone average remains an ablation;all16 outputs/negative results are
retained. Native delivery still requires its separate installed-wheel gates.
[Complete field report](../dev/partition_averaged_stack_20261110/field_v1/report.json),
[completed manifest](../dev/partition_averaged_stack_20261110/field_v1/run_manifest.json),
[all5792 scores](../dev/partition_averaged_stack_20261110/field_v1/scores.tsv),
[resource measurement](../dev/partition_averaged_stack_20261110/analysis_23178680/resources.txt).
Earlier checkpoints below retain their original context.

Matched field comparison frozen and submitted as23178680. Keep ALL14 preceding
methods plus the new primary/ablation,ALL181 populations and both1MbPearson/
armSpearman endpoints. Require ALL5068 previous scores EXACT. Fixed family30
retains original strict missing-family rule:unavailable competitors block ALL
formal superiority verdicts;common-complete point means cannot replace them.
No predictor,penalty,reference,mask or endpoint changes. Native delivery is
undergoing separate validation and may not be inferred from field scoring.

The job was PENDING before any field output when its node pin was cleared;
pure scoring of frozen predictions can run on any compatiblev6/AMD CPU. Same
job,source,runtime,nslab account/QOS,1CPU16GB1h and exact scientific gates.
All installed-package and inference jobs remain pinned005.
[Field design](../dev/partition_averaged_stack_20261110/FIELD_DESIGN.md),
[field driver](../benchmarks/hardening/partition_average_field.py),
[command](../dev/partition_averaged_stack_20261110/field_job.json),
[resource-only amendment](../dev/partition_averaged_stack_20261110/field_resource_amendment.json),
[929-file source preserving all927 predecessor files](../dev/partition_averaged_stack_20261110/snapshot_v3.json).
Earlier checkpoints below retain their original context.

Full181-line regression COMPLETE:all five jobs23176490–23176494 exit0:0.
All181 old profiles and1448 old scores are EXACT;all workers lock RNA before
DNA and keep the full original population. The augmented stack is selected in
169 lines. Mean1Mb DNA Pearson is.7135355679987977 versus current gene-Laplace
.711350890510652,previous RNA stack.7072257052738345,SM5 pseudobulk
.6679553825846039 andMDL pseudobulk.6683618389831502. It wins133/181 versus
current,loses36 and ties12;versus SM5 it wins175 and loses6.

ALL8 prespecified conditional paired-mean contrasts pass BOTH Holm and basic
simultaneous-interval gates. The current-model gain is.002184677488145502,
simultaneous interval[.00144679825675895,.003085100797539899];Holm p
.00039998000099995. Max conditional export-arithmetic Pearson bound is
6.935875433388507e-7,well below even the smallest simultaneous lower bound.
These are approximate eight-reference-group/20,000-draw results on EXPOSED
regression,not calibrated coverage or independent competitor superiority.
No bulk/spatial/clone benefit or default promotion is established.

Workers take9m25s,9m51s,9m29s and9m00s scheduler time;summary26s. First measured
completed worker process peak426084KiB;all resource files remain alongside the
executions. The successful gate permits separately verified opt-in package
integration;research return observation must be replaced by explicit local
data passing and reproduce ALL181 frozen candidate predictions before delivery.
[Full report](../dev/partition_averaged_stack_20261110/summary_v1/report.json),
[completed manifest](../dev/partition_averaged_stack_20261110/summary_v1/run_manifest.json),
[all181 paired scores](../dev/partition_averaged_stack_20261110/summary_v1/paired_profiles.tsv).
Earlier checkpoints below retain their original context.

Two-case panel COMPLETE:23176320 exits0:0 in34s scheduler time,28.82s measured
process time;process peak290704KiB,scheduler146392KiB. BOTH original profiles
and ALL16 old scores are EXACT. All RNA predictions lock before DNA;all recorded
input hashes remain unchanged. The single augmented stack passes the frozen gate:

| Locked case | Current stack | Averaged-expert stack | SM5 pseudobulk | MDL pseudobulk |
|---|---:|---:|---:|---:|
|2313287_STOMACH|.8155510221769611|.8177903729411115|.7487632408673744|.7668579268685296|
|42MGBA_CENTRAL_NERVOUS_SYSTEM|.8372884967941583|.8384877003182553|.7760888357017742|.8103601099901083|

Both certified RNA-loss improvements use positive new-expert weights
.3509835389 and.5771200994. Gains exceed conditional export-arithmetic Pearson
bounds2.79893e-7 and4.52918e-7. Standalone averages score.8088463857 and
.8334682889;they remain ablations. No calibrated biological uncertainty,
independent validation or general superiority follows from two exposed cases.
[Actual panel](../dev/partition_averaged_stack_20261110/panel_v1/report.json),
[completed manifest](../dev/partition_averaged_stack_20261110/panel_v1/run_manifest.json).

The full181-line protocol was separately frozen after this pass and BEFORE
any new full-regression outputs. Four workers23176490–23176493 and dependent
summary23176494 are submitted,initially Priority/Dependency. ALL181 old profiles
and ALL1448 scores must be EXACT. Keep all8 prior comparators in a fixed Holm/
simultaneous-interval family;no comparator or line drops. The same kernel hash,
RNA-before-DNA locks and8-group/20,000-draw conditional inference are enforced.
All8 comparisons must pass before separately verified opt-in integration;
no default promotion or bulk/spatial transfer is implied.
[Full frozen protocol](../dev/partition_averaged_stack_20261110/REGRESSION_DESIGN.md),
[driver](../benchmarks/hardening/partition_average_regression.py),
[commands](../dev/partition_averaged_stack_20261110/regression_jobs.json),
[927-file snapshot preserving all925 earlier files](../dev/partition_averaged_stack_20261110/snapshot_v2.json).
Earlier checkpoints below retain their original context.

Numerical23176319COMPLETE0:0 onne1dc6-009 in20s scheduler time. Independent
uniform-rational error1.78e-15 and exhaustive raw-gene partition error5.33e-15
are within the derived1.43616e-10 effect-arithmetic allowance. One-bin,
translation/scaling and EXACT held-target exclusion checks pass. Time-v13.57s,
process peak116736KiB;scheduler MaxRSS0 is unavailable sampling,not zero usage.
The existing HiGHS threads-option warning remains visible. Real-profile23176320
is now eligible but PENDING(Priority);no DNA performance result yet.
[Numerical report](../dev/partition_averaged_stack_20261110/numerics_v1/report.json),
[execution/resources](../dev/partition_averaged_stack_20261110/analysis_23176319/resources.txt).
Earlier checkpoints below retain their original context.

Pre-execution placement amendment:cleared ONLY pending numerical23176319's
required-node pin,retainingv6&amd,CPU1/4GB/1h,innovation account/QOS and the same
source/runtime. The exact-profile panel remains pinned to005. No new job,
restart,method change or gate relaxation. Original command/design preserved;
[reason and before/after record](../dev/partition_averaged_stack_20261110/RESOURCE_AMENDMENT.md).

Implemented ONE research predictor that averages all contiguous partitions under
the existing gene-Laplace working score. Numerical23176319 and locked two-case
panel23176320 are submitted. No executable result or accuracy gain yet.

| Job | Purpose | Account/QOS | CPU / memory / time | Dependency |
|---|---|---|---|---|
|23176319|Independent numerical references|innovation|1 /4GB /1h|none|
|23176320|Two locked exposed profile cases|nslab|1 /8GB /2h|afterok23176319 plus actual numerical report|

Both use original005/v6&amd,capturedruntime003 and the validated old private
noneditable wheel. Initial numerical wait is Resources;panel waits on Dependency.
No array,interactiveQOS,new dependency,shared installation or frozen/default edit.

The estimand is relative regional RNA abundance in natural-log units against the
specified reference. The mechanism follows the earlier bulk point improvement:
hard segmentation discards uncertainty among plausible boundaries. The strongest
alternative is that averaging merely blurs real dosage-shaped signal,so better
RNA fit need not improve independent DNA prediction. Matched-resolution DNA on
the locked cases is the immediate distinguishing observation;corresponding new
biological replication remains necessary for broader genomic claims.

Read-only observation obtains exact integer interval costs and medians from the
ACTUAL installed old kernel. Forward/backward summation weights every partition
by exp(-existing_score/2),without sampling,truncation or a new temperature.
Old scale,3log(training genes) penalty,genomic resets,folds,reference,targets and
all existing experts remain. These are Gibbs weights on fitted working losses,
not integrated Bayesian evidence or calibrated biological boundary probabilities.

Add one averaged column to the same simplex L1 stack. Masked targets remain
absent from their predictor and scale fit. Retain old output unless improvement
exceeds optimization and arithmetic/averaging budgets. This is conditional fitting,
not unbiased held-out error. Standalone averaging is an ablation only. Weight
estimation stability and biological uncertainty are not supplied by arithmetic
bounds conditional on the fitted float coefficients.

Numerics require independently enumerated Fraction expectations over all8uniform
partitions,independently rebuilt raw-gene median/loss partitions,scalar likelihood
sums,one-bin/translation/scaling checks and exact masked-target perturbation
invariance. Actual block inclusion mass must equalone within a derived bound.
No scientific epsilon,probability clipping or post-hoc renormalization hides a
failure. The panel explicitly reads the passing report before analysis.

Use the SAME exposed cases,in order:2313287_STOMACH and42MGBA_CENTRAL_NERVOUS_SYSTEM.
Require both old profiles and all16old scores EXACT,same callable genes,reference
DNA contrast and1Mb mean projection. Freeze both RNA outputs before parsing DNA.
The new stack must beat current gene-Laplace,SM5 pseudobulk AND MDL pseudobulk
on BOTH cases,beyond the conditional export-rounding Pearson perturbation bound,
with certified RNA-loss gain and positive new-expert weight. An ablation win
cannot rescue a failed stack. A valid negative result stops expansion without
penalty,temperature,fold,case or comparator tuning. A pass permits separately
freezing a full regression;it does not establish independent superiority.

All923 predecessor files remain exact in the925-file immutable snapshot. AST
and launcher syntax pass;executable checks pending. Existing clone jobs continue
unchanged. No new cohort,biological replication,count uncertainty calibration,
DNA identity or false-split guarantee is created. Former evaluation is regression-only.

[Frozen design](../dev/partition_averaged_stack_20261110/DESIGN.md),
[implementation](../benchmarks/hardening/partition_averaged_stack.py),
[source preservation](../dev/partition_averaged_stack_20261110/snapshot_v1.json),
[exact commands/responses](../dev/partition_averaged_stack_20261110/jobs.json),
[bulk mechanism evidence](POSTERIOR_REGION_AVERAGE_2026-11-09.md).
