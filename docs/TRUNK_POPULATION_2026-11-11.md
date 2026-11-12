# Shared-profile training population: paired performance mechanism

The scientific goal remains accurate RNA-derived genomic profiles and supported
CNA populations,independently validated at the claimed resolution. The preceding
posterior-action candidate fails186/390 whole-profile wins against196 required.
This experiment tests a separate source-observed training-population difference;
it does not retune that action or replace recovery,false-split or uncertainty
criteria. No new biological validation is claimed.

**Final result:** larger shared-profile training improves cell-specific agreement
to0.311803 but reduces whole-profile agreement to0.602302. Primary187/390 wins
fail required196/390. All ablations fail continuation too. Independent counts,
profiles,DNA scoring and numerical decision checks pass;reject promotion. Three
nslab jobs complete,one failed preflight preserved. Frozen criteria and results
remain unchanged;the full scientific goal is still unachieved.

Historical checkpoint:original-context baseline qualification passes exactly. Candidate retry
23207930 and dependent independent reader23207931 are submitted. No new candidate
DNA-performance result is available at this checkpoint.

## Scientific question,estimand and competing explanations

The current shared gene-median trunk uses260 query cells in two training folds.
The fixed pseudobulk comparator uses389 target-excluded query cells. Both use the
same450 external-reference cells. My initial conversational description that
pseudobulk used all390 cells was imprecise:source inspection confirms target
exclusion. Different amounts of shared-profile information could explain part
of the remaining whole-profile gap. RNA/DNA response mismatch,cell state,gene
dependence and reference composition remain competing explanations.

The estimand remains a descriptive relative RNA profile against the specified
external reference,plus a within-specimen cell response. The frozen response
was fitted under a gene-NB/Gamma working model in three folds. Use its existing
point estimates without refitting likelihood,noise or prior. Changing the origin
of a fitted abundance `q=p_fold*exp(delta)` requires
`delta_LOO=delta+log(p_fold)-log(p_LOO)`. This preserves q;it does not produce a
newly fitted posterior under a389-cell prior. Gene medians and regional means
remain different summaries,and counts are rounded deposited expression.

Primary:389-cell LOO trunk plus reanchored frozen expected-log response. Controls:
LOO trunk alone;LOO trunk with unreanchored response;260-cell trunk with reanchored
response;LOO trunk with reanchored mode response. Keep all27 old methods and add
five new outputs. No candidate mode may replace a failing primary.

All390 cells,572 regions,fold identities,reference,genes,common callable DNA
support and five-region within-block median stay fixed. Smooth the combined
profile once. Require both mean Pearson endpoints above the immediately preceding
expected-log candidate0.603150638/0.308352828,not just the older mode. Keep the
196/390 majority gate,original anchored-estimator gates and REFCON comparisons.
Individual-response benefit additionally needs improvement over the matched LOO
trunk;otherwise attribute a gain to shared-profile estimation. A failed primary
rejects the unchanged mechanism without fold/weight/prior/noise/action/smoothing
sweeps. Future independent specimens with matched DNA must establish transport.

[Frozen design](../dev/trunk_population_20261111/DESIGN.md),
[executed candidate](../dev/trunk_population_20261111/study.py),
[preceding complete rejection](POSTERIOR_GENE_ACTION_2026-11-11.md).

## Measurement chain and exact execution contexts

Use the locked scONE RNA matrix,annotated normal/query identities,frozen gene
columns,region/block labels and corresponding processed same-cell DNA. All
target-excluded counts are computed directly from that matrix. No new DNA is
used in fitting,feature choice or action selection. New RNA predictions must
be hash-locked before scoring. Biological n=1 exposed specimen;cells and genes
are not independent biological specimens. Genomic dosage,cellular clone identity
and calibrated uncertainty remain unestablished.

First job23207901 stops before DNA scoring. The reconstructed pseudobulk differs
at9640/223080 values,max2.220446049250313e-16. This is an exact replay failure;
its source and full log remain. No tolerance is relaxed. Inspection traces the
original pseudobulk generation to AMDv6,archived stack_scaling environment and
runtime003;the later scoring and frozen response calculations use Intelv4.

Job23207913 recreates the original row-wise baseline on AMDv6 with that runtime.
All223080 values match exactly. It writes canonical LOO trunks,IDs and hashes,
without opening DNA. The Intelv4 retry imports those trunks,verifies the canonical
baseline/IDs/hash,and independently requires exact old mode,expected-log and
fold-trunk profiles. Cross-context trunk differences are recorded separately;
no universal cross-platform bitwise assertion is made.

[Original failed source](../dev/trunk_population_20261111/attempt1_manifest.json),
[full failed check](../dev/trunk_population_20261111/run_23207901/output.txt),
[source-context qualification](../dev/trunk_population_20261111/baseline_23207913/REPORT.json),
[baseline reconstruction](../dev/trunk_population_20261111/baseline.py),
[execution qualification and reader contract](../dev/trunk_population_20261111/REPLAY_CONTEXT.md).

## Computational and scientific checks

The fit checks source-count conservation,independent exact integer pooling,
own-cell perturbation exclusion,permuted cell identities,all three prior profile
replays and abundance reanchoring. Fitted responses are reused,not recomputed.
The independent reader imports no candidate/scorer. It reconstructs pooled gene
sums and trunks in long double for the first target cell per fold,checks all
regional means as exact rational pseudocount fractions,and reconstructs all five
new profile arrays with independently sorted-window medians. Frozen response
numerical diagnostics are propagated through the new arithmetic;there is no
claim of a new full-PMF run.

It independently reconstructs corresponding DNA,whole/cell-specific correlations,
all32-method means and paired effects. All10530 old cell records must match
exactly. A primary decision overlapping numerical diagnostics is unresolved,
requiring precision resolution. Controls equivalent to pseudobulk may have
unresolved numerical near-ties;they cannot establish superiority. Only the
prespecified primary can pass a continuation gate. No confidence bounds or
patient-level p-values are fabricated from390 cells.

[Independent reader](../dev/trunk_population_20261111/reader.py),
[reader source manifest](../dev/trunk_population_20261111/reader_manifest.json),
[candidate source manifest](../dev/trunk_population_20261111/source_manifest.json).

## Compute and preservation

All jobs use nslab account/QOS,CPU partition,one CPU,8GB,one-hour limit;no arrays.
Baseline qualification uses AMDv6 and the archived read-only runtime;fit/reader
use Intelv4 and the existing read-only allelic_qtl environment. No installation,
public/default/test change,frozen-output overwrite or scientific criterion change.

| Job | Work | Checkpoint state |
|---|---|---|
|23207901|Initial combined-context exact replay|Failed1:0 before DNA|
|23207913|Original AMD-context baseline qualification|Complete0:0,baseline exact|
|23207930|Qualified Intel-context mechanism fit/score|Submitted|
|23207931|Independent reconstruction and decisions|Afterok:23207930|

[Initial working tree](../dev/trunk_population_20261111/initial_git_status.txt),
[288 initial ledger entries](../dev/trunk_population_20261111/claims_before.json),
[91 prior source identities](../dev/trunk_population_20261111/prior_source_identities.json),
[first submission](../dev/trunk_population_20261111/submission.json),
[baseline submission](../dev/trunk_population_20261111/baseline_submission.json),
[qualified retry](../dev/trunk_population_20261111/retry_submission.json),
[dependent reader](../dev/trunk_population_20261111/reader_submission.json).

The frozen production recovery2/5,continuous-state false splits14/200 and
correlated selected-interval coverage54/1000 remain. Complete scientific status
requires the pending paired results and independent reader. The full genomic,
clone,spatial and competitor-superiority goal remains active and unachieved.

## Completed performance results and falsification

| Method | Whole mean Pearson | Cell-specific mean Pearson | Whole wins over pseudobulk |
|---|---:|---:|---:|
|Primary LOO trunk,reanchored expected log|0.602302429|0.311802682|187/390|
|LOO trunk,reanchored mode ablation|0.602235828|0.311709730|187/390|
|LOO trunk,unreanchored expected log|0.602636759|0.312124707|186/390|
|Fold trunk,reanchored expected log|0.602153107|0.307113466|186/390|
|LOO trunk alone|0.561258387|-0.048241612|Numerically equivalent to pseudobulk|
|Previous expected-log primary|0.603150638|0.308352828|186/390|
|Pseudobulk|0.561258387|-0.048241612|Comparator|
|Reference-refined REFCON|0.615822434|0.376007356|Retained comparator|

The primary changes the two means by-0.000848209/+0.003449854 and adds one
whole-profile win. All mean margins and its187–187 numerical win range resolve
outside propagated arithmetic diagnostics. Neither more shared-profile data nor
regional reanchoring repairs the full performance failure. Individual response
still beats the LOO trunk,but this does not satisfy the fixed joint gate.
The higher cell-specific ablation result cannot replace a failing primary.
Reject unchanged training-population adjustment without a fold/weight/action/
prior/noise/smoother sweep. No gain is established on independent biological units.

The trunk-only method is a per-cell rescaling of the original pseudobulk,so
whole and normalized cell-specific correlations are mathematically equivalent.
Its nominal floating-score win count130/390 has no scientific interpretation:
numerical propagation permits0–390. Retain that unresolved control comparison
explicitly;do not count numerical near-ties as meaningful wins or losses.

## Independent verification and preservation

The candidate exactly replays the original mode,expected-log andfold-trunk
profiles and all10530 old cell records. The original-runtime baseline worker
exactly reproduces223080 pseudobulk values. Intel/AMD raw trunk discrepancy is
2.220446049250313e-16 and is recorded without a universal bitwise guarantee.
Count conservation,integer pooling,own-cell perturbation exclusion,identity
permutation and the algebraic change of origin pass. Maximum abundance identity
discrepancy4.439855545237673e-15 is below its64-operation diagnostic.

Reader23207931 independently checks1716 exact rational regional means,three
source-RNA trunks(max2.684e-16),all1115400 new profile values across five arms
(maximum relative discrepancy2.3583220442974745e-15),128 rational method means,
96 paired comparisons and3900 long-double primary DNA correlations(max4.083e-14).
All10530 preceding records remain exact. Primary and all four response-arm
majority decisions resolve;the mathematically equivalent trunk-only comparison
correctly remains numerically unresolved. No confidence or calibrated genomic
membership statement is inferred from these numerical results.

The first failed-source copy used basenames,colliding between two distinct
`study.py` files. The partial copy is retained. The complete original tree is
recovered into `attempt1_reconstructed` with repository-relative paths;each of
all seven original execution identities matches its originally recorded SHA256.
The two subsequently changed files are restored by reversing only their known
execution-context edits and requiring exact original hashes before writing.
This is verified source recovery,not an inferred approximation or a replaced
failed result. Future snapshots must preserve relative paths.

| Job | Result | Elapsed | Node |
|---|---|---|---|
|23207901|Failed1:0,before DNA; exact baseline preflight|45s|ne1dc4-016|
|23207913|Complete0:0; exact original-context baseline|16s|ne1dc6-003|
|23207930|Complete0:0; candidate and exact replays|1m05s|ne1dc4-068|
|23207931|Complete0:0; independent reader|22s|ne1dc4-016|

All use nslab/CPU,1CPU/8GB/1h. Detailed process/scheduler RSS and environments
are retained separately. All91 preceding nonbinary source/data/test identities
and288 prior ledger values remain. Add the new negative result without revising
older values or historical decisions. No public/default/test,shared environment,
frozen-output or scientific criterion change. Original recovery,error/power and
uncertainty failures remain unchanged.

[Full independent result](../dev/trunk_population_20261111/reader_23207931/results/SUMMARY.json),
[paired comparisons](../dev/trunk_population_20261111/reader_23207931/results/paired_cells.tsv),
[source-context recovery record](../dev/trunk_population_20261111/attempt1_reconstruction.json),
[all original source identities](../dev/trunk_population_20261111/attempt1_manifest.json),
[complete accounting](../dev/trunk_population_20261111/accounting_complete.tsv),
[process/scheduler resources](../dev/trunk_population_20261111/resources_complete.json),
[preservation record](../dev/trunk_population_20261111/preservation_complete.json).
