# Fixed-horizon replication of region-excluded residual precision

COMPLETE: the unchanged estimator passes all96 simultaneous comparisons and all12
preservation safeguards on256 new independent base draws. Independent reconstruction
verifies every decision and all136,192 error scores. This is a conditional simulation
milestone, permitting a separately frozen real-count regression challenge. It does
not promote a default or establish biological/SOTA performance.

The earlier execution checkpoint recorded one fixed-horizon simulation replication.
The preceding32-draw study remains failed at94/96 comparisons with12/12 safeguards;
none of those draws enters this replication. No biological or SOTA conclusion is
established by submitting or completing these jobs.

## Method and independent horizon

The scientific target remains accurate genomic gain/loss profiles and supported
CNA-defined populations across cellular, spatial and bulk RNA. This experiment
isolates one high-count profile component: normal-only prediction of transcriptional
state, with regional log weights estimated from gene residuals whose predictors
exclude the gene's entire region. The full estimator source is unchanged from
[the preceding experiment](WITHIN_REGION_CONDITIONAL_2026-11-11.md).

There are256 new base seeds96–351, paired across the unchanged7 laws and2 depths:
3,584 condition datasets, zero biological specimens. The80 normal/80 query
observations,480 genes,24 regions, log(1.5) arm contrast, normal fit/selection
identities, ridge path and all19 methods remain fixed. All eight worker prediction
sets must lock before the pooled scorer opens any new truth. No interim performance
inspection or optional stopping is allowed.

The original96 primary comparisons retain positive mean improvement and exact
one-sided95% simultaneous lower strict-win probability above1/2, Bonferroni
alpha=.05/96. All12 point-mean no-increase safeguards remain. A pass requires all108.
The old14-case preflight is repeated for execution qualification, but it remains
14 unique old datasets, not112 independent specimens or new evaluation draws.

The fixed horizon uses the smallest development win proportion23/32 as a planning
ASSUMPTION, not a confidence bound. Exact binomial calculations require155/256 wins
for a sign criterion. At that assumed win probability, the per-comparison miss
probability is0.0000355466, giving joint sign-criterion power at least0.9965875 by
the union bound across96 comparisons. This covers only the sign criteria under
the assumed law, not mean effects, safeguards, biological power or the5% clone
false-split target. Smaller block-doubling horizons32/64/128 have uninformative
raw union bounds; their negative values are lower bounds, not negative probabilities.

Marker information and reference-to-query transport remain identifying assumptions.
Unmarked coherent state and reversed transport stay in the full report. Gene-specific
dosage response, breakpoints inside a region and biological measurement uncertainty
remain limits. Further simulation precision cannot supply corresponding DNA, suitable
single-CNA null specimens or biological replication. A passing result would only
justify separately freezing a suitable real-data challenge; no automatic default
promotion, genomic confidence claim or clone-recovery improvement follows.

## Compute and execution amendment

All jobs use nslab account/QOS,CPU partition,1CPU/8GB/1h,existing Conda Python and
isolated imports. Eight workers are explicit jobs, not an array. The source and
original design remain byte-frozen. All scripts, full submissions, seeds, environments,
input/prediction/truth hashes and process resource logs are retained.

Workers23206252–23206259 complete on the original node003. Pooled scorer23206283
runs after all eight succeed; independent reader23206286 waits for the scorer.
No new-draw partial performance is inspected.

While seven workers were pending, qualifier23206279 tests node012 on the14 OLD
cases and reproduces all252 arrays and14 complete parameter files exactly. No
replication truth is opened. By qualification completion all pending workers have
started on node003; zero workers are moved and every running job remains untouched.
The unused alternate-node qualification and placement amendment are preserved.

The pooled scorer was submitted before its separate resource-header announcement.
The omission was acknowledged immediately in chat; its recorded nslab1CPU/8GB/1h
allocation and dependency did not change. All other batches received their header
before submission. This execution-record issue does not alter the frozen method,
criteria, fixed horizon or scoring isolation.

Current checkpoint: eight workers and qualification complete; pooled scoring and
its independent reader are not yet reported as complete. The final scientific
decision remains pending; this checkpoint is superseded only by a completed final
report below. The active full goal remains unachieved.

[Frozen protocol](../dev/within_region_replication_20261111/DESIGN.md),
[unchanged estimator identities](../dev/within_region_replication_20261111/frozen_source.sha256),
[worker orchestration](../dev/within_region_replication_20261111/worker.py),
[all eight submissions and seed allocations](../dev/within_region_replication_20261111/WORKERS.json),
[exact planning calculation](../dev/within_region_replication_20261111/worker_0/results/FIXED_SAMPLE_PLANNING.json),
[pooled scorer](../dev/within_region_replication_20261111/pool.py),
[scorer submission](../dev/within_region_replication_20261111/pool_submission.json),
[independent reader](../dev/within_region_replication_20261111/read_results.py),
[reader submission](../dev/within_region_replication_20261111/reader_submission.json),
[placement amendment](../dev/within_region_replication_20261111/RUNTIME_PLACEMENT.md),
[alternate-node qualification](../dev/within_region_replication_20261111/qualification_23206279/results/QUALIFIED.json),
[actual placement outcome](../dev/within_region_replication_20261111/runtime_placement_outcome.json).

## Completed final horizon and independent result

| Requirement | Final result | Decision |
|---|---|---|
| Simultaneous profile comparisons |96/96; weakest190/256 wins, lower bound0.644306923|Pass|
| Original clean/pure-dosage/noise safeguards |12/12|Pass|
| Independent error-score reconstruction |136,192/136,192|Pass|
| Exact rational/integer decision recount |108/108|Pass|
| Query-intervention isolation |Both fixed cases, weights and model choices unchanged|Pass|

The weakest comparisons remain the lower-depth Poisson comparisons with the
Poisson-precision method. Their mean improvements and simultaneous bounds now
pass on the independent fixed horizon. The development94/96 gate remains failed;
those32 seeds are not pooled here. No horizon extension, method change or weaker
criterion was used in this replication.

Selected depth1e6 whole/cell-specific MSE on the new256 draws:

| Condition | Earlier within-region | New residual precision | All-normal count pseudobulk |
|---|---:|---:|---:|
| Clean | .000065874 / .000063512 | .000026458 / .000025529 | .00791661 / .00790177 |
| Shared state | .000082879 / .000080016 | .000034722 / .000033522 | .00908116 / .00790091 |
| Pure dosage | .000068442 / .000065646 | .000029171 / .000028023 | .00845659 / .00790188 |
| Noise only | .00265037 / .00255507 | .00201329 / .00194570 | .00798781 / .00790408 |
| Noisy state | .00287184 / .00276502 | .00228244 / .00220147 | .00913768 / .00790481 |

The planted contrast remains near log(1.5). These baselines pool individual query
observations within the simulated population; they must not be transplanted as
cross-patient bulk pooling. Both depths, all7 laws and all19 methods remain in the
full report. The additional cohort size improves precision of the declared
simulation comparison, not biological replication or measurement validity.

The limits are consequential. Under reversed transport, depth1e6 whole MSE is
0.1966265 versus raw0.0507351; unmarked-state MSE is0.0506444 versus raw0.0505521.
The method does not solve those identifying failures. It needs informative state
markers and a valid reference relationship. Fitting, simulation or comparator
agreement cannot independently verify those biological assumptions.

The separate reader checks266,112 region weight normalizations,2,647 exact
no-correction profiles and all recorded region exclusions. Maximum relative score
reconstruction discrepancy is5.72e-16. Repeated old-case checks total1,120 original
and672 weighted arrays with112 parameter files; these cover the same14 unique old
cases. The separate node012 qualification checks252 arrays and14 files. Every
result is exact at its declared regression gate; no numerical tolerance decides
a scientific pass.

## Completed compute and preservation

| Jobs | Work | Progress | State/errors | Memory | Downstream |
|---|---|---|---|---|---|
|23206252–23206259|Eight workers|8/8 ██████████|Complete0:0;0 errors|8GB each; scheduler peaks76,156–76,748KiB|All prediction locks present|
|23206279|Alternate-node qualification|252/252 ██████████|Complete0:0|8GB; short-run RSS not sampled|Qualified; no placement changes needed|
|23206283|Pooled scoring|3,584/3,584 ██████████|Complete0:0|8GB; scheduler97,284KiB|Reader dependency succeeds|
|23206286|Independent reader|108/108 ██████████|Complete0:0|8GB; scheduler109,772KiB|All checks pass|

Running0;pending0;completed11;job errors0. Worker runtimes are3m26s–3m44s,
qualification14s,pooled scoring3m16s and reader1m47s. Peak process RSS across these
jobs is138,836KiB; process logs and scheduler samples are retained separately.
All use nslab CPU. Across this goal turn, all15 submitted jobs complete without
failure. Existing unrelated jobs are untouched.

All18 execution/source/configuration identities and91 validated nonbinary public,
data and test identities remain unchanged. No public package, original test,
dependency or shared installation changes accompany this milestone. The recorded
521-Python/four-Rust suite remains the applicable package validation; it is not
represented as newly executed. Frozen outputs and historical decisions are intact.

**Next admissible action:** freeze a real-count U-CAN U0001/U0008 regression
challenge using the120 known normal donors and the original count/gene/bin
contracts. Retain the whole fixed native/simple-baseline family and matched
smoothing operators; lock both new RNA profiles before either FACETS DNA read.
Check computational equivalence before caching repeated feature calculations
for that larger gene universe. These two already exposed specimens remain
regression-only and cannot establish independent biological superiority. No
actual real-count fit was run by this replication.

Clone recovery remains2/5; the continuous-state false-split result remains14/200;
correlated selected-interval coverage remains54/1,000. The full genomic profiling,
clone reliability and SOTA objective remains active and unachieved.

[Complete96+12 report](../dev/within_region_replication_20261111/pool_23206283/results/REPORT.json),
[all68,096 score rows](../dev/within_region_replication_20261111/pool_23206283/results/scores.tsv),
[all condition means](../dev/within_region_replication_20261111/pool_23206283/results/means.tsv),
[all predictions locked before truth](../dev/within_region_replication_20261111/pool_23206283/results/RNA_PREDICTIONS_LOCKED.json),
[independent reconstruction](../dev/within_region_replication_20261111/reader_23206286/RECOUNT.json),
[all job outcomes](../dev/within_region_replication_20261111/accounting_complete.tsv),
[process resource inventory](../dev/within_region_replication_20261111/resources_complete.json),
[final source integrity](../dev/within_region_replication_20261111/source_integrity_complete.json).
