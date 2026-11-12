# Depth-aware reference noise: recovery gain fails the false-split gate

The reference-only depth-quintile candidate improves the exposed original recovery
panel from production2/5 to research4/5, newly recovering20261108 and20261110.
It nevertheless fails the frozen full-workflow error/power experiment:13/100
state-null false splits versus1/100 production, with no power gain. **Reject
promotion.** Production, original tests, frozen results and the5% criterion remain
unchanged. This report retains the promising pilot and its contradictory final result.

## What changed and why

[Estimator](../benchmarks/hardening/reference_gene_noise.py) retains the previous
cross-cell factorial-moment common NB dispersion calculation, now within declared
reference type × depth quintile. Five strata come from the existing reference
depth-model resolution, not a DNA-selected grid. Each stratum estimates rate² by
((sumY)²-sumY²)/((sumL)²-sumL²), then pools excess factorial moments over mean².
Only fit-reference cells enter; no query, confirmation or holdout reference values,
truth labels or simulator noise are used. Nonpositive dispersion fails explicitly.
Weights, effective-count gene veto, proposals, posterior.95/core/gene/cell support
and original six recovery assertions are unchanged. New arguments preserve old
defaults and frozen source snapshots. No production package code changes.

The estimand is recoverable regional RNA-profile population structure, under the
NB and pooled-RNA working models. Locally stable reference gene rates within
depth strata identify the common dispersion; residual depth dependence, cell
state, reference composition and heterogeneous gene noise are alternatives.
Estimated dispersion is not pure technical noise and does not establish independent
halves or biological calibration. Corresponding DNA and independent biological
replication would establish genomic interpretation. Biological n=0.

[Recovery design](../dev/depth_stratified_gene_noise_20261109/DESIGN.md),
[pilot design](../dev/depth_stratified_gene_noise_20261109/PILOT_DESIGN.md) and
[predeclared final continuation](../dev/depth_stratified_gene_noise_20261109/FINAL_DESIGN.md)
record when each decision was frozen. Source_v1 executes the original five-case
arms; source_v2 the pilot; source_v3 the continuation and exact summary.

## All original recovery outcomes

Hard pooling passes3/5; weighted confirmation and merging passes4/5. The fitted
common phi values in seed order are.114465,.119190,.113781,.119622,.111266.
These values are measured reference moments, not tuned toward the simulator.

| Seed | Pass all | Normal FPR (<.1) | Tumour TPR (>.9) | S2 groups (=1) | S1 groups (=2) | Original assignment (>.65) | Difference chromosomes (=chr12) |
|---|---|---:|---:|---:|---:|---:|---|
|20261031|True|0.044444|0.968627|1|2|0.836032|['chr12']|
|20261107|True|0.066667|0.988235|1|2|0.738095|['chr12']|
|20261108|True|0.088889|0.992157|1|2|0.756917|['chr12']|
|20261109|False|0.022222|0.980392|1|1|0.602000|[]|
|20261110|True|0.066667|0.980392|1|2|0.814000|['chr12']|

Assignment is the unchanged tumour-call-conditional original criterion; it does
not replace the separate exactly-two-group and chromosome checks. All30 criterion
values and original assertions, including failures, remain in the
[full report](../dev/depth_stratified_gene_noise_20261109/both_v1/report.json).
The remaining109 case now rejects its chr12 proposal at initial weighted HMM
confirmation, before merging. Its mean membership margin is.658452, but there is
no qualifying segment. Thus a merge-only correction cannot recover this executed
case. No posterior threshold was relaxed to accept it.

## Paired fresh error and power, complete workflow

Replication amendment: conditions share the base RNG seed and reference counts
within a replicate. The old frozen JSON field `independent_simulations=300`
(pilot30) therefore overcounts cross-condition independence. Correct totals are
300 condition datasets from100 independent base seeds, or30/10 in the pilot.
Each condition still has100 independent replicates, so its binomial bounds and
the Bonferroni decisions remain valid. The original JSON stays unchanged; future
summary keys distinguish condition datasets and base seeds. The independence
wording in the following historical description is superseded by this amendment.

Phase3 generates matched, continuous-state and two-population conditions, each
with100 independent simulated specimens. Each is fit by production and candidate:
300 simulations/600 complete workflows, including reference estimation and all
adaptive stages. The first10 are the retained pilot, not extra independent data.
Continuous state has one CNA population; two-profile condition adds a CNA-shaped
signal with the same nuisance. Biological specimens remain0.

The pilot showed0/10 errors for both arms under both nulls, and candidate10/10
versus production9/10 recovery. The exact upper error bound was35.48%, so it
never certified5%. The final100 horizon was fixed before continuation. Two nulls
× two looks use alpha.0125 per one-sided exact bound/excess-error test, preserving
the simultaneous family. Do not extend the horizon after a failed bound.

| Condition | Production | Candidate | Candidate error upper bound | Decision |
|---|---:|---:|---:|---|
| Matched null false splits |0/100|0/100|4.2874%|Conditional matched-null bound passes|
| Continuous-state null false splits |1/100|13/100|22.4272%|Candidate rejected|
| Exact chr12 two-profile recovery |96/100|95/100|Not an error endpoint|No power improvement|
| Mean all-cell two-profile assignment |.719967|.717333|Not an error endpoint|No membership improvement|

The exact candidate state-null excess-error p=.00146435 is below.0125. This
is stronger than merely failing a confidence bound. Production state upper
bound6.2119% also fails the strict bound, although1/100 is not evidence that its
population rate exceeds5%. Matched-null success does not generalize to state.
Two-profile ordinary exact two-sided95% recovery intervals are[.900743,.988996]
production and[.887165,.983568] candidate. They describe this fixed simulated
regime, not universal power; four paired recoveries favor candidate and five
favor production. All-cell denominator600 includes unassigned/incorrect cells.
[All performance endpoints](../dev/depth_stratified_gene_noise_20261109/final_summary_v1/performance.json)
and [paired report](../dev/depth_stratified_gene_noise_20261109/final_summary_v1/report.json)
preserve both arms and every condition. Exact excess tails are independently
calculated as integer binomial sums; bounds use SciPy beta quantiles.

## Execution, files and next action

All28 CPU jobs complete0:0: two recovery arms, six pilot chunks and summary,
eighteen continuation chunks and final summary. Each workflow job requests
4CPU/16GB/2h; summaries1CPU/4GB/2h. Both innovation/nslab, no arrays or interactive
QOS. Validated noneditable wheel/environment unchanged. Actual runtime, memory,
backend, imports, wheel/native/source hashes and dependencies are retained in
[job journal](../dev/depth_stratified_gene_noise_20261109/jobs.json),
[accounting](../dev/depth_stratified_gene_noise_20261109/accounting.tsv) and
[process resource records](../dev/depth_stratified_gene_noise_20261109/resource_summary.json).
Full package/GPU tests: not rerun for research-only changes. Earlier403 Python/
4Rust passes and1GPU-only skip remain separate evidence.

Reject the candidate as a default despite4/5 exposed recovery. A targeted
continuous-state alternative safeguard is the next performance intervention,
with recovery checked before any large error/power expansion. More gene-noise
constants, repetitions of these five cases, or smaller p-values cannot establish
genomic identity. Profile work proceeds separately; neither track can mask the
other’s failed gate. The full scientific goal remains unmet.
