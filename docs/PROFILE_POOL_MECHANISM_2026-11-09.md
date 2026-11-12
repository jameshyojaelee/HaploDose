# Aggregation preserves an estimand but does not supply a profile advantage

On the existing simulated raw-count fixture, count-pool aggregation slightly
reduces squared error versus the MAP-weighted mean clone log effect, but does
not beat the naive ratio to the known generating reference. All summaries worsen
under the same RNA-selected smoothing operator when scored against generating
RNA means. This is one mechanism fixture, not an accuracy or superiority estimate.
No aggregate, smoothing change, reference correction or default is promoted.

## Primary measurements and target

Use the already executed `profile_pool_export_20261109/parity_v1` fixture, without
refitting. The installed simulator, seed20261109 and all original parameters
reproduce every raw count, gene and cell identity exactly. A generation-boundary
observer retains `L*share` before gamma/Poisson sampling. Independent scalar
recalculation instead starts from each gene's rate, depth and slope, normalizes
with `math.fsum`, then applies the exposure. All67200 gene means agree, with
maximum absolute difference1.14e-13 within an operation-count/scale allowance.
This is a checked numerical reference, not formal interval arithmetic.

The fixture has280 cells:120 reference cells across two types/splits and160 query
cells. Only the30 fit cells of generating type0 define the specified reference.
The48-bin target is the median-centred log ratio of normalized generating pooled
query RNA means to normalized generating reference RNA means. Depth variation,
normal admixture and two generating CNA populations remain present. The reference
identity is supplied by the generator, not independently learned from RNA.

The target is **relative pooled RNA abundance**, not mean DNA copies, calibrated
genomic state or correct clone membership. Reference mismatch, transcriptional
effects and nonlinear RNA response remain competing biological explanations.
All values are conditional on the retained feature universe and this realized
depth structure. Corresponding DNA and biological replication would be needed to
establish the genomic interpretation; neither is added here. Biological n=0,
conditional simulated specimens=1. Genes/bins/cells are not replicate specimens.

## Fixed comparisons and controls

All clone summaries include the reference-like component and use the same
centering convention and48-bin support. Fitted expectations and clone log effects
are those actually exported; no fitted parameter is changed. Mean-log weights
are fixed MAP cell counts, responsibility mass or observed count mass. The pooled
summary sums observed/expected fitted pools before one half-count offset.

Naive uses the specified type0 reference **counts**; oracle-reference-rate uses
the generating type0 reference **means**, retaining query count noise. The latter
is unavailable to a deployed method. Model expectations include fitted reference
mixture/nuisance effects; this is not a matched-reference-estimation ablation or
a fair deployable-comparator benchmark.

The previously fixed RNA-only partition rule estimates first-difference noise
from naive within-chromosome contrasts, with MAD/(sqrt(2)*normal_q75), penalty
3*sigma²*log(n) and the existing minimum-bin rule. Naive RNA selects four chromosome
blocks. The SAME boundaries are applied to every method; truth remains unsmoothed.
No target, DNA or accuracy metric selects boundaries, weights or parameters.

| Method | Unsmoothed RNA MSE | Common-boundary RNA MSE |
|---|---:|---:|
| MAP-cell mean log effect |0.0037666|0.0078331|
| Responsibility-mass mean log effect |0.0036181|0.0077509|
| Count-mass mean log effect |0.0035936|0.0077394|
| Pooled fitted observed/expected counts |0.0036024|0.0077478|
| Naive, known reference identity |0.0034325|0.0075034|
| Oracle reference means |0.0008430|0.0073181|

These are descriptive errors in natural-log RNA units squared. MAP mean-log has
unsmoothed Pearson0.92633, versus0.92322 for naive and0.92325 for pooled fitted
counts: the MSE result is not dominance across metrics. All12 MSE/Pearson/Spearman
rows remain in `metrics.tsv`; no comparator is dropped. Smoothing removes some
true RNA variation in this target; that does **not** establish worse DNA accuracy.

An independent exposure example illustrates the estimand difference: two equally
numerous populations with rates1 and2 but exposures1 and9 have pooled effect
log(1.9), whereas their mean log effect is0.5*log(2). Neither becomes mean DNA
copy number by renaming it. Splitting each exported count/expectation pool into
two equal pieces preserves the pooled effect within binary64 accumulation error;
the half-count correction is applied once after summing, never per new component.

## Independent checks and execution

The audit verifies all original input/output hashes, derives gene means with
scalar rate/depth arithmetic, independently sums13440 bin means, and checks all48
target bins. Independent scalar MSE, centered-product correlation and explicit
tie-average ranks recount all36 metrics; maximum discrepancy3.33e-16. Four common
regions partition the support exactly without crossing chromosomes. Raw pooled
observations agree with directly summing the160 assigned cells' original bins.

All computation ran on compute nodes using the validated noneditable environment
`dev/profile_pool_export_20261109/environment` and its wheel SHA256
`2ba2a00107a89bc0aa899c8afaf274134544c47f287336521cce53008280cd42`.

| Job | Account/QOS | CPUs / RAM | Elapsed | Process peak RSS | Exit |
|---|---|---|---:|---:|---:|
|23075083|innovation|2 /4GB|9s|193036KiB|0|
|23075423|nslab|2 /4GB|3s|84688KiB|0|

Task `dev/profile_pool_mechanism_20261109` contains the prespecified design,
immutable source snapshots, `panel_v1/{primary_means_and_profiles.npz,metrics.tsv,
regions.json,report.json}`, `audit_v1/report.json` and execution/resource manifests.
No package/native/test/README changes or new dependencies occurred. The full385-test
Python/four-test Rust suite was not rerun because production source was unchanged.
Frozen results, original criteria and historical decisions remain intact.

## Decision

Retain exact count-pool exports for their correct interpretation and auditability.
This experiment does not justify presenting pooling or smoothing as the accuracy
repair. The oracle-reference sensitivity points to reference estimation as a
mechanism to investigate, not a demonstrated deployable correction. Reference
uncertainty must be tested without using true rates or DNA outcomes for fitting.

Kinker count-scale and historical correspondence gaps remain unresolved; this
fixture cannot fill them. Original recovery2/5, the three failing seeds, the failed
false-split/power tradeoff and uncertainty failures remain unchanged. No genomic,
biological-calibration or SOTA claim follows from this computational result.
