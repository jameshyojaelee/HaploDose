# Empirical regional prior using every eligible target-excluded cell

The scientific goal is more accurate cell-specific RNA dosage profiles while
preserving minority CNA-shaped contrasts and the strong pseudobulk baseline.
Genomic validation remains a separate requirement. This experiment combines
two diagnosed changes without adding a parameter-selection search.

## Mechanism and fixed test

The [119-cell moment-Beta prior](LEAVE_ONE_CELL_PRIOR_2026-11-11.md) improves
whole-profile MSE over its 60-cell predecessor in all four heterogeneous strata,
but fails 6/32 comparisons and retains high-depth minority bias. The
[finite empirical prior](REGIONAL_EMPIRICAL_PRIOR_2026-11-11.md) reduces minority
bias but uses only 60 eligible cells and fails clean low-depth whole-profile
superiority. The new estimator fits that SAME empirical prior to all 119
eligible cells, excluding each target from its own support and weights.

The observation model is binomial regional count margins at fixed total;
the estimand is individual regional RNA log odds against the unchanged external
reference. Preserve the finite support rule, Dirichlet parameter 1+1/M,
total MAP pseudocount one, convex optimisation and expected-logit action.
No prior-strength, action, support, fold, iteration or region sweep. Regional
transcription remains a competing genomic explanation. RNA MSE and state-free
genomic-response error are separate endpoints.

Use the same 288 exposed generated count cases, all 120 cells and 24 regions,
48 base draws per scenario and paired depths 300/3000. Preserve all baseline
inputs, references, genes and common support. Pseudobulk retains all 119 eligible
cells. Primary family remains all 32 whole/minority comparisons in both
heterogeneous scenarios/depths against pseudobulk, individual, global shrinkage
and original two-fold Beta expected logit. Every simultaneous one-sided 95%
lower strict-win bound must exceed one half, Bonferroni .05/32. Numerical
overlap and ties count against the candidate. A favourable mean is insufficient.
All RNA predictions lock before truth scoring. Biological n=0; this is exposed
mechanism development, not independent population performance.

## Implementation and independent evidence

[Frozen design](../dev/loo_empirical_prior_20261111/DESIGN.md),
[study](../dev/loo_empirical_prior_20261111/study.py),
[unchanged empirical solver](../dev/loo_empirical_prior_20261111/prior.py),
[independent reader](../dev/loo_empirical_prior_20261111/reader.py),
[source identities](../dev/loo_empirical_prior_20261111/SNAPSHOT.json).

Four separate case workers fit exact naive target-deleted counts. No target
value enters its own prior. Saved per-target support, weights, optimisation
residuals, distance bounds and predictions permit independent reconstruction.
Qualifications include exact rational optima, independent SciPy optimisation,
finite differences and target perturbation. The reader reconstructs every
target's support and complete-binomial-PMF action, checks independent fsum
errors, replays all 1440 original baseline/Beta control rows and recounts all
32 decisions with exact integer binomial tails. Numerical tolerances derive
from operation counts and strong-convexity residual bounds, not a scientific
epsilon. Unresolved calculations cannot earn a win.

## Execution

| Job | Purpose | Allocation | Dependency |
|---|---|---|---|
| 23213184 | 288 cases, 829440 target-region fits | CPU AMD-v6, 4 CPUs, 8 GB, 3 h | None |
| 23213185 | Independent reader | CPU AMD-v6, 1 CPU, 4 GB, 2 h | afterok 23213184 |

All account/QOS nslab. Existing isolated Conda environment/runtime003, one BLAS
thread per worker, neutral directories and source hashes before/after. No arrays,
production edit, dependency installation or modification of existing outputs.
Before submission: AST/bash checks, all 91 protected identities and immutable
source snapshot passed.
[Exact submission commands](../dev/loo_empirical_prior_20261111/submissions.json).

Results and independent validation are pending at this checkpoint. Failure
rejects this fixed composed estimator; it does not justify weaker criteria or
a tuning sweep. A pass would permit a separately frozen transfer experiment,
not automatic default promotion or genomic validation.

## Reader I/O correction, with the running original preserved

Fit23213184 completes0:0 in3m54s. Its initial32-way decision is28 passes and
four failures;the complete scientific gate fails. Independent verification
remains required before accepting those computed results.

Source inspection finds unnecessary repeated decompression in the reader's
target/region loop:each `fits[name]` access loads the entire compressed array.
This is a new-reader performance defect;it does not alter the estimator or
scientific test. Running23213185 and its frozen source remain untouched.
New reader23213545 loads each case's arrays once into memory,then performs the
EXACT same complete checks and numerical tolerances. Cache lifetime is one case.
A separate source contract verifies the one-block substitution;compute-node
qualification compares every cached field to direct archive reads. No subset
of cells,regions,cases,comparisons or numerical checks is dropped.

Reader23213545 requests nslab/nslab CPU AMD-v6,1CPU4GB1h,afterok the same
completed23213184. It writes a separate result directory. The original reader's
eventual report remains available for comparison;its observation time is not
treated as a failed check. No running job is cancelled or changed.
[Exact cache correction and contract](../dev/loo_empirical_reader_cache_20261111/DESIGN.md),
[reader source](../dev/loo_empirical_reader_cache_20261111/reader.py),
[compute qualification](../dev/loo_empirical_reader_cache_20261111/qualify_cache.py),
[exact submission](../dev/loo_empirical_reader_cache_20261111/submission.json).

## Completed independent result: 28/32; general gate fails

Cached reader23213545 completes0:0 in1m31s by scheduler accounting;the full
scientific checks take1m27.94s. Its independently reconstructed829440 actions
have maximum discrepancy6.661338147750939e−15,within the unchanged numerical
diagnostics,with zero unresolved actions. Every target-deleted support and
119-cell training count verifies. All1440 historical baseline/Beta control
rows match exactly. Independent fsum endpoints and all32 integer binomial
decisions agree with the fit report. Cache qualification confirms all11 fields
(5736960 bytes) exactly and proves the single source substitution.

| Condition | Depth | Pseudobulk whole RNA MSE | Original Beta whole MSE | New empirical whole MSE |
|---|---:|---:|---:|---:|
| Two CNA populations | 300 | 0.004132589 | 0.005169724 | 0.005449287 |
| Two CNA populations | 3000 | 0.003484575 | 0.001268426 | 0.000688085 |
| Two CNA populations plus state | 300 | 0.016196525 | 0.011178183 | 0.011538307 |
| Two CNA populations plus state | 3000 | 0.015544627 | 0.002313576 | 0.002036112 |

All16 minority-region comparisons pass. The four failures are whole-profile
comparisons at depth300:

| Condition | Comparator | Resolved wins /48 | Simultaneous lower win bound |
|---|---|---:|---:|
| Two CNA populations | Pseudobulk | 1/48 | 0.000032577 |
| Two CNA populations | Global shrinkage | 5/48 | 0.017656 |
| Two CNA populations | Original Beta expected logit | 16/48 | 0.153698 |
| Two CNA populations plus state | Original Beta expected logit | 16/48 | 0.153698 |

No numerical boundary overlaps explain those four failures. More eligible
training cells plus the finite empirical prior do not establish general
superiority;the low-depth whole-profile loss persists. High-depth gains and
minority preservation remain positive conditional results,not a reason to
change the original gate or select a depth cutoff retrospectively. Reject
general promotion without a support/prior/action/fold sweep.

[Inference] Finite training information is not the sole explanation of the
earlier empirical prior failure. Low-depth prior estimation/distribution
regularisation remains a candidate mechanism;this experiment does not isolate
its contribution. Successful complete-PMF reconstruction separates numerical
correctness from this statistical failure. None of these simulated RNA results
establishes individual-cell genomic accuracy or biological clone membership.

Fit process peak92684 KiB;cached reader109700 KiB;qualification has a separate
resource record. Scheduler fit/cached-reader peaks124228/89848 KiB differ due
to sampling and process accounting. Original inefficient23213185 remains
RUNNING and untouched at this checkpoint;it has no downstream dependency.
No job or numerical check has failed. Its eventual output is retained,not
required instead of the completed identical independent check.

[Authoritative independent result](../dev/loo_empirical_reader_cache_20261111/reader_23213545/results/REPORT.json),
[cache/source qualification](../dev/loo_empirical_reader_cache_20261111/reader_23213545/CACHE_QUALIFIED.json),
[all32 original decision records](../dev/loo_empirical_prior_20261111/run_23213184/results/comparisons.tsv),
[all score records](../dev/loo_empirical_prior_20261111/run_23213184/results/scores.tsv),
[fit resources](../dev/loo_empirical_prior_20261111/run_23213184/resources.txt),
[reader resources](../dev/loo_empirical_reader_cache_20261111/reader_23213545/resources.txt).

No production/default or scientific criterion changes. All old results and
decisions remain. This completes the fixed experiment,not the full scientific
goal. Research recovery5/5 remains;state14/200 and selected coverage54/1000
still preclude the corresponding validation claims.
