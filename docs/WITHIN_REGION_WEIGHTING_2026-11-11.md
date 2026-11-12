# Regional weighting fixes Poisson sampling loss but does not pass the full gate

The matched Poisson-precision arm fixes all clean/pure-dosage safeguards and
improves on the preceding state predictor under Poisson noise. It worsens extra
gene noise. Total normal-log-variance weighting handles noise-only samples but
still loses dosage-only cell-specific accuracy. Count sums introduce additional
nonlinear state error. All three complete candidate gates fail. No arm expands
onto real data or changes production. This is a follow-up to the
[preceding mechanism experiment](WITHIN_REGION_NUISANCE_2026-11-11.md).

The objective remains genomic gain/loss profiles and supported CNA-defined
populations across single-cell, spatial and bulk RNA. This study tests a
specific high-count observation-model component, not biological identification,
clone membership, absolute copy number or uncertainty calibration.

## Mechanism, assumptions and frozen analysis

For independent Poisson genes with known reference means and a common regional
multiplier, the region sum is sufficient for that multiplier. First-order
variance of an equal-gene log average is sum(1/m_g)/G², versus 1/sum(m_g) for
the log sum. Cauchy-Schwarz establishes the ordering. Unequal gene abundance
therefore wastes information in the equal-gene summary. Heterogeneous NB/Gamma
noise breaks the Poisson weighting argument; the experiment explicitly retains it.

Three fixed arms change both the full-count RNA summary and its opposite-half
normal state-prediction target: normal-mean-proportional log weights, inverse
normal log-variance weights, and regional count sums. Weights use the 40 fitting
normals only. Normal identities 40–79 select the unchanged 14-member ridge path.
Queries select no weights or penalties. The variance weights are *total normal
variation*, which includes reference state; they are not known sampling noise.
Count sums have a nonlinear log-sum-exp state response, unlike a linear weighted
log target. These are competing explanations for a failure of simple weighting.

All 448 former datasets remain regression-only. Their 4,480 original usable
profile arrays and 448 parameter files replay exactly before any new scores.
The independent generator uses 32 new seeds, 32–63, under the same seven laws
and two depths: 448 new datasets, zero biological specimens. Each has 80 normals,
80 queries, 24 regions, 20 genes per region and a balanced log(1.5) arm contrast.
Depths 1e6/1e7 apply to the 480 retained genes, not the whole assay. Normal/query
transport and within-region marker information remain identifying assumptions.

The original eight comparators plus the preceding within-region estimator give
72 primary comparisons per arm. The entire three-arm family has 216 one-sided
exact 95% bounds, Bonferroni protected, requiring positive mean improvement and
lower strict-win probability above 1/2. Every arm also retains all 12 original
no-mean-increase safeguards versus raw gene means. No method can qualify by
choosing its best conditions retrospectively. All 896 predictions lock in an
observation-only process before separate new-truth scoring.

## Complete outcome and causal implication

| Candidate | Primary comparisons | Preservation safeguards | Complete gate |
|---|---:|---:|---|
| Poisson precision | 68/72 | 8/12 | Fail |
| Normal log variance | 64/72 | 10/12 | Fail |
| Count sum | 64/72 | 8/12 | Fail |

All arms pass the 64 comparisons against the original baseline family. Their
additional failures are against the preceding within-region predictor. Poisson
precision passes its four Poisson-state comparisons but loses all four noisy-
state comparisons; the other two arms lose all eight. This is not superiority
over the best available candidate.

Selected new-seed mean whole/cell-specific MSE at depth 1e6:

| Condition | Raw mean | Previous within | Poisson precision | Normal log variance | Count sum |
|---|---:|---:|---:|---:|---:|
| Clean | .000063640 / .000061609 | .000063667 / .000061637 | .000025516 / .000024669 | .000026711 / .000025813 | .000025506 / .000024660 |
| Shared state | .050032 / .048316 | .000080865 / .000078123 | .000035045 / .000033788 | .00021763 / .00021042 | .00050056 / .00048308 |
| Pure dosage | .0014029 / .000062308 | .000067441 / .000064720 | .000028573 / .000027460 | .00019705 / .00019036 | .00021437 / .000030646 |
| Noise only | .0027411 / .0026384 | .0027098 / .0026091 | .0054223 / .0052236 | .0020489 / .0019769 | .0051832 / .0049914 |
| Noisy state | .052773 / .050920 | .0028652 / .0027619 | .0056867 / .0054804 | .0032011 / .0030904 | .0061900 / .0059686 |

The Poisson correction's remaining four safeguards fail only in the two
noise-only depths and endpoints, increasing mean error by .00259–.00275.
Normal-log-variance weighting fails pure-dosage cell-specific error at both
depths, by .000128053 and .000036807. Count sums also fail the four noise-only
safeguards. Every failed decision and both depths remain in the full report.

[Inference] Poisson weighting resolves the earlier sampling-efficiency mechanism,
but ignores non-Poisson gene variance. Total-variance weighting includes state
variation, so it cannot isolate that noise either. The count-sum error is
consistent with its nonlinear state response; no nonlinear corrective model was
selected from these outcomes. Reversed transport remains damaging for every arm;
unmarked state remains unresolved. Additional count depth does not establish a
valid reference-to-query bridge.

The next justified question is whether a reference-only conditional residual can
separate state variation from gene noise. A predictor must exclude the target
gene's entire region: Gamma effects shared across molecule halves would otherwise
let a gene predict its own fluctuation. That is a new estimator requiring a new
frozen study and new seeds; it does not rescue these failed gates or reopen their
real-data continuation.

## Executed checks, compute and preservation

Both nslab CPU jobs complete with exit 0:0. One CPU, 8 GB, one hour and the original
Intel v4 node003 are used. No arrays, new dependencies, public package changes,
original-test changes or shared environment modifications occur.

| Job | Work | Progress | State | Elapsed |
|---|---|---|---|---|
| 23206120 | Numerical tests, replay, new generation/fit/score | 896/896 ██████████ | Complete | 3m23s |
| 23206152 | Independent reader and query isolation | 252/252 ██████████ | Complete | 19s |

Running 0; pending 0; complete 2; job errors 0. The dependent reader succeeds.
Independent rational weights/Fisher variance, equal-weight reduction and noiseless
dosage protection pass. All 15,232 new error scores reconstruct in extended
precision, with maximum relative discrepancy 3.80e-16. All 252 scientific
decisions agree with rational sums/integer tails. All 43,008 region-weight sums,
455 exact no-correction profiles and both fixed query-perturbation cases pass.
Query perturbations leave weights, penalties and selection losses unchanged.
Numerical checking allowances never decide or relax scientific gates.

All original frozen outputs, confidence/error criteria, historical decisions and
the failed preceding experiment remain unchanged. Biological n=0. Production
recovery remains 2/5, state false splits 14/200 and correlated selected-interval
coverage 54/1,000. No new biological or uncertainty claim follows. Independent
corresponding DNA and replicated RNA/reference transport evidence would be
required to extend a conditional profile improvement to genomic accuracy.

[Frozen design](../dev/within_region_weighting_20261111/DESIGN.md),
[matched estimators](../dev/within_region_weighting_20261111/weighted.py),
[executed workflow](../dev/within_region_weighting_20261111/study.py),
[source identities](../dev/within_region_weighting_20261111/source.sha256),
[submission](../dev/within_region_weighting_20261111/submission.json),
[complete decisions](../dev/within_region_weighting_20261111/run_23206120/results/REPORT.json),
[all individual scores](../dev/within_region_weighting_20261111/run_23206120/results/scores.tsv),
[all condition means](../dev/within_region_weighting_20261111/run_23206120/results/means.tsv),
[numerical fixtures](../dev/within_region_weighting_20261111/run_23206120/results/NUMERICS.json),
[all independent checks](../dev/within_region_weighting_20261111/reader_23206152/RECOUNT.json),
[reader command](../dev/within_region_weighting_20261111/reader_submission.json),
[all accounting](../dev/within_region_weighting_20261111/accounting_complete.tsv),
[process resource measurements](../dev/within_region_weighting_20261111/resources_complete.json),
[final source identities](../dev/within_region_weighting_20261111/source_integrity_complete.json).

The separately frozen region-excluded residual follow-up subsequently passes all12
preservation safeguards and94/96 comparisons, while retaining the complete-gate
failure. [All outcomes](WITHIN_REGION_CONDITIONAL_2026-11-11.md).
