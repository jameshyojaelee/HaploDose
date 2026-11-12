# Gene-level segmentation improves the single-cell profile regression

The optional gene-level segment expert raises mean1Mb DNA Pearson from
0.7072257053 to0.7113508905 over all181 exposed regression lines, versus
smoothed pseudobulk0.6679553826. All six prespecified conditional comparisons
pass both simultaneous interval and Holm rules. The SAME method fails the bulk
replacement gate:0.6966950787 versus strongest pseudobulk0.7302613671 and retained
posterior-averaged profile0.7648167082. Retain a single-cell descriptive option,
reject bulk transfer, and do not change any default, clone criterion or frozen
artifact. Former evaluation remains regression-only.

## Concrete method and falsifiable mechanism

The existing simplex L1 stack predicts individual gene RNA effects, while its
Gaussian segment expert minimizes squared errors of bin medians. The new
expert instead fits gene medians on contiguous bin segments by absolute error,
retaining gene multiplicity and matching the stack's prediction loss. It is a
working Laplace loss on prepared natural-log RNA contrasts, not a Laplace,
Poisson or NB observation model for counts.

On each training split, estimate b as median(abs(gene−training-bin-median))/log2,
the Laplace MAD identity. Zero/nonfinite b makes the expert unavailable; no
invented floor. Minimize2*sum absolute gene residual/b plus3log(N training genes
in the block) for each ADDED segment. The3log convention is retained as a
working complexity penalty on the actual observations, not a calibrated BIC
claim under dependent genes. Never cross contiguous arm/gap blocks. No DNA
selects genes, reference, scale, penalty, weights or outputs.

All observed float effects, scale and specified float penalty are converted to
exact rational units. Two heaps supply interval medians and integer residual
sums; dynamic programming costs and comparisons are exact integers. Only
fitted medians/diagnostic objectives return to float. Genuine ties choose the
first predecessor and remain recorded. No universal epsilon determines a
scientific choice. The current implementation uses quadratic interval-table
memory within each block; larger-resolution scalability is not established.

The original median-width and Gaussian experts, rank-mod10 folds and singleton
rules remain. Each masked target is excluded from its base predictions AND
scale fit. Add ONE new predictor column to the SAME simplex absolute-error fit,
retaining the old convex hull. Use the augmented output only if it improves
RNA fit over the actual old deployed output beyond both optimization-gap and
operation/magnitude rounding budgets; otherwise return the old output exactly.
Weights themselves use these targets, so the fit is not independent confirmation
or an unbiased CV score. Upstream transformations and references remain fixed.
Standalone gene-Laplace is an ablation, never a selectable second candidate.

Estimand: regional relative RNA abundance in natural-log units against the
specified reference. Identifying assumption: a shared regional gene location
is predictable beyond gene-specific deviations. Coordinated transcription or
reference composition is the strongest alternative and can have identical RNA
observations. Corresponding independently processed DNA at matched resolution
can test candidate genomic dosage; RNA smoothing and Pearson alone cannot
validate cell membership, absolute copy number, purity, WGD or ancestry.

[Initial design](../dev/gene_laplace_profile_20261110/DESIGN.md),
[post-gate regression design](../dev/gene_laplace_profile_20261110/REGRESSION_DESIGN.md),
[research implementation](../benchmarks/hardening/gene_laplace_profile.py),
[native expert](../python/haplodose/gene_segments.py),
[public APIs](../python/haplodose/stacking.py).

## Independent computation, locked cases and all181 regressions

An independently written Fraction calculation enumerates every partition of a
four-bin fixture, scores direct gene residuals without the DP heaps, and agrees
with the exact optimum19965496953594612478567523625077 /
685812462948444219108304617472. Translation by2^20, multiplication by8, within-bin
permutation, masked-target scale exclusion, genomic resets and zero-scale
abstention pass. These are computational invariants, not sampling calibration.
[Numerical results](../dev/gene_laplace_profile_20261110/numerics_v1/report.json).

Before expanding the run, both locked diagnostic lines improve against BOTH
preceding stack and smoothed pseudobulk, with all12 previous scores exact:

| Locked line | Previous stack | Augmented stack |
|---|---:|---:|
|2313287_STOMACH|.8100292815|.8155510222|
|42MGBA_CENTRAL_NERVOUS_SYSTEM|.8296615480|.8372884968|

The estimator remains byte-identical between source_v1 and source_v2. All181
lines replay in original order with identical references, genes, bins, truth
contrasts, masks and1Mb scorer. Every1086 prior scores is EXACT. Predictions,
parameters, decisions and hashes save before DNA is opened. The relative DepMap
truth remains an uncertain reference, potentially from different culture material.
No new biological replication or independent evaluation is created.

| Comparator | Mean1Mb Pearson | Primary gain | Simultaneous basic interval for gain |
|---|---:|---:|---:|
|Raw bin pseudobulk|.5959399404|.1154109501|[.0989428402,.1333497332]|
|Smoothed pseudobulk5|.6679553826|.0433955079|[.0360236451,.0503106698]|
|Gaussian-MDL pseudobulk|.6683618390|.0429890515|[.0358722751,.0502966301]|
|Predictive median|.6802181879|.0311327026|[.0254862420,.0367874589]|
|Predictive partition|.6859045160|.0254463745|[.0183339531,.0324338401]|
|Previous RNA stack|.7072257053|.0041251852|[.0022104254,.0052823583]|

Primary augmented mean.7113508905; standalone ablation.6950356807. The primary
beats the old stack132/181, loses43 and ties6; versus smoothed pseudobulk173 wins
and8 losses. These paired means have the same estimand, tests and intervals.

Fixed family6;20,000 paired whole draws of the8 connected biological/reference
groups, seed2026111003, same draws for every contrast. Each Holm p=.0002999850;
all pass Holm<.05 AND simultaneous Bonferroni interval lower>0. These approximate
conditional intervals hold upstream fitting fixed and do not claim established
coverage with8 groups or erase repeated exposure to this cohort.
[Every score and contrast](../dev/gene_laplace_profile_20261110/summary_v2/report.json).

## Bulk transfer is rejected, including its apparent recall gain

All48 old gene contexts,24 development/24 evaluation patients, normal-reference
inputs, support, loss thresholds and40 preceding methods are retained exactly.
Two new columns make42 methods; standalone remains ablation only. Correlation
is defined on23 evaluation patients, loss recall on13, FPR on23; undefined
outcomes and noncalling patients remain recorded.23 patients are not independent
of shared reference estimation, and these are exposed regression outcomes.

| Method | Mean Pearson | Loss recall | Calling-patient precision | Loss FPR |
|---|---:|---:|---:|---:|
|Augmented stack|.6966950787|.6103021978|.2907557720 /22|.0538906884|
|Standalone ablation|.6728971530|.6345695971|.2839703857 /22|.0590495183|
|Previous RNA stack|.7035111344|see full table|see full table|see full table|
|Strongest pseudobulk|.7302613671|.3996794872|.3585714286 /14|.0412444510|
|Retained posterior average|.7648167082|.6027472527|.4837962963 /18|.0437371741|

The augmented stack wins8/23 against strongest pseudobulk and2/23 against the
retained best profile; it loses15 and21 respectively. RNA fit selects it47/48
contexts yet DNA agreement worsens. This directly falsifies a shared bulk-benefit
claim from better RNA prediction alone. Increased recall includes increased FPR;
it is not a free sensitivity gain. The frozen bulk point gate fails, so no
bulk uncertainty expansion or default promotion follows.
[All42 methods and failures](../dev/gene_laplace_profile_20261110/bulk_v1/report.json),
[complete patient table](../dev/gene_laplace_profile_20261110/bulk_v1/patient_scores.tsv).

## Package delivery and matched field comparison

Both profile APIs accept gene_segments=True. The original default remains exact;
no CLI, bulk, clone, dependency or shared-installation changes. The private
noneditable wheel passes456 Python/4 Rust tests,1 GPU-only skip. Eight warnings
are the existing SciPy forwarding of the HiGHS threads option. Numerical tests
cover exact exhaustive partition scoring, translations/scales, masking, genomic
resets, unavailable scale, invalid mode/identity, unchanged default behavior
and count-API forwarding.

The two original identified count matrices reproduce all new RNA gene effects,
profiles and common decision fields EXACT without reading DNA: stomach103 query/
2118 reference cells/10133 genes,42MGBA517/2467/10116. These source measurements
are reconstructed counts and acquire no raw-count calibration. All181 installed new profiles/common diagnostics and previous default profiles/
diagnostics reproduce EXACT without reading DNA. Four disjoint workers cover
each original case once; the final ordered coverage check passes.
[Full installed parity](../dev/gene_laplace_profile_20261110/native_parity_summary_v2/report.json).
[Count-chain report](../dev/gene_laplace_profile_20261110/count_parity_v1/report.json),
[wheel/test evidence](../dev/gene_laplace_profile_20261110/validation_23156464/).
Wheel SHA256 d0d5183953d4114142047b2c8970269e840502ec076c582d209703e5b710c16a;
native SHA256 e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7,
unchanged from the preceding delivery.

The matched field comparison projects the already frozen new predictions by
gene identity onto the previous common-support RNA tables. All12 previous methods,
all181 populations, masks and4344 scores reproduce EXACT, including NaNs. The
new RNA tables freeze before opening relative DNA truth. Same-cell DNA membership
validation is not involved.

| Method | Mean1Mb Pearson | Mean arm Spearman |
|---|---:|---:|
|Augmented stack|.7070308343|.8501066693|
|Previous RNA stack|.7029699692|.8474459987|
|Smoothed pseudobulk5|.6621648151|.8357367749|
|Gene-level pseudobulk|.3515545796|.8491078572|
|inferCNV R|.6684765573|.8132941230|
|infercnvpy|.6252941920|.7980361783|
|CopyKAT|.6094027636|.7369578393|
|SCEVAN|.5682890162|.6844726204|

These point means use167 complete-family lines under the SAME common support.
The earlier arm disadvantage versus gene-level pseudobulk becomes a small point
gain, not a formal superiority finding. CopyKAT13 and SCEVAN14 unavailable
populations remain missing in the full181-line table. The prespecified26-test
family (13 comparators x2 endpoints),20,000 whole8-group paired draws and
seed2026111004 retain the previous rule: incomplete family makes ALL formal
superiority tests untestable. Zero formal verdicts pass. Do not drop failing
methods, shrink the family or relabel this old cohort independent evaluation.
All14 methods, standalone ablation and conditional intervals remain recorded.
[Field design](../dev/gene_laplace_profile_20261110/FIELD_DESIGN.md),
[all field results](../dev/gene_laplace_profile_20261110/field_v1/report.json).

## Compute, preservation and unresolved clone recovery

[Exact submission commands and all jobs](../dev/gene_laplace_profile_20261110/jobs.json).
All CPU/v6&amd, innovation/nslab, no arrays or interactive QOS. Each analysis
runs from a neutral directory with immutable source manifests, wheel/native
identity checks, explicit actual backend, protected frozen paths and recorded
runtime/memory. Setup copies existing Conda/pip dependencies into a new private
prefix; no shared installation or previous environment is modified.

| Job | Role | Account/QOS | CPUs /GB /time | Outcome |
|---|---|---|---|---|
|23156295|Exact numerical references|nslab|1 /4 /1h|complete0:0,19s|
|23156296|Locked two-line gate|innovation|1 /8 /2h|complete0:0,20s|
|23156297|All48 bulk contexts|nslab|1 /16 /2h|complete0:0,291s; scientific gate fails|
|23156356|SC0:46|innovation|1 /8 /2h|complete0:0,379s|
|23156357|SC46:92|nslab|1 /8 /2h|complete0:0,385s|
|23156358|SC92:138|innovation|1 /8 /2h|complete0:0,380s|
|23156359|SC138:181|nslab|1 /8 /2h|complete0:0,356s|
|23156360|Fixed-family regression summary|nslab|1 /4 /1h|complete0:0,14s|
|23156463|Private environment copy|innovation|1 /8 /1h|complete0:0,235s|
|23156464|Noneditable wheel/full Python and Rust tests|nslab|4 /32 /2h|complete0:0,157s|
|23156465|Unfinished single-worker181-case parity|innovation|1 /8 /2h|cancelled for four-worker replacement; logs kept|
|23156480|All181 matched field comparison|nslab|1 /8 /1h|complete0:0,83s|
|23156492|Two original count-matrix chains|innovation|1 /8 /1h|complete0:0,22s|
|23156594|Installed parity0:46|innovation|1 /8 /2h|complete0:0,00:06:01|
|23156595|Installed parity46:92|nslab|1 /8 /2h|complete0:0,00:06:08|
|23156596|Installed parity92:138|innovation|1 /8 /2h|complete0:0,00:06:01|
|23156597|Installed parity138:181|nslab|1 /8 /2h|complete0:0,00:05:41|
|23156598|Require ordered, exhaustive181-case parity|nslab|1 /4 /1h|complete0:0,00:00:11|

The cancelled parity attempt changed no method or scientific decision. Its frozen
source/logs remain; four disjoint workers execute the SAME exact comparisons,
then a dependent summary requires the original ordered identifiers0..180 with
no omission/duplication and all default/new profile/diagnostic checks passing.
No tests are relaxed and no failed scientific outcome is replaced.

All18 job attempts are terminal:17 complete0:0 and1 deliberate replacement
cancellation, no failed scientific/test job. Peak scheduler RSS 2003368K (23156464.batch);
pytest measured2243152KB. All exact commands, per-step RAM/runtime and statuses
are retained in [accounting](../dev/gene_laplace_profile_20261110/accounting_final.tsv)
and [resource summary](../dev/gene_laplace_profile_20261110/resource_summary_final.json).

No production clone or original recovery test is changed. Production remains2/5;
the earlier population-aligned score recovered5/5 but failed4/10 state-null
false splits. Its later predictive state safeguard lost original109 and is
rejected. This profile improvement does not repair that tradeoff or establish
5% whole-specimen error control. No new calibrated intervals, biological clone
panel or spatial validation have appeared. The full scientific/SOTA objective
remains incomplete; all failed hypotheses and original scientific criteria stay.

## Changed files and reproduction

Production adds [gene_segments.py](../python/haplodose/gene_segments.py) and
changes only the optional branches in [stacking.py](../python/haplodose/stacking.py).
All earlier low-level solvers/experts are byte-identical.
[Test fixtures](../tests/test_gene_segments.py) are new; no earlier test changes.
README documents the opt-in keyword and failed bulk transfer.

Research sources are gene_laplace_profile.py, gene_laplace_panel.py,
gene_laplace_summary.py, gene_laplace_native_parity.py, gene_laplace_field.py,
gene_laplace_count_parity.py and gene_laplace_parity_summary.py, all under
benchmarks/hardening. The task directory preserves DESIGN, REGRESSION_DESIGN,
FIELD_DESIGN, baseline, every immutable executed source, Conda/pip inventory,
wheel, outputs, locks and cancellation reason. jobs.json records executable
argument arrays for every sbatch command. Do not rerun into an existing output
directory; preserve these artifacts and use a new task/output location.

Current STATUS/DEVELOPMENT_LOG and generated indexes add this checkpoint;
the claim ledger appends one field while preserving all157 earlier values.
Of82 previous Python/native/test files, only the intended stacking.py wrappers
change. Original tests, frozen data/results and historical decisions remain.
