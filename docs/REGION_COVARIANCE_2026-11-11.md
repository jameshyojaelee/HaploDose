# Cross-region RNA covariance does not improve the exposed cell profiles

[Observed pattern, not a guarantee] The new joint regional covariance estimator
fails its declared continuation and mechanism gates. Whole/cell-specific mean
Pearson is **0.567226/0.256930**, below both the previous gene-NB candidate
**0.602558/0.308225** and the matched diagonal ablation **0.576824/0.301815**.
Whole-profile wins over pseudobulk are **156/390**, failing strict majority.
Reject this unchanged candidate. No public method, default, original test,
scientific criterion or frozen result changes.

This is concrete negative performance evidence from a new implemented estimator,
not a package upgrade or evidence of superiority. The full genomic-profile,
clone-membership and spatial objectives remain active and unachieved. Previous
production recovery 2/5, continuous-state false splits 14/200 and correlated
selected-interval coverage 54/1,000 remain unresolved.

## What was tested and why

The previous goal turn was a status report, without new scientific execution.
The shared-noise count split subsequently recorded in that status already fails
7/9. Earlier whole-cell confirmation also retained state-driven false distinctions.
Neither rejected design is repeated here. The new hypothesis concerns regional
profile accuracy: independent regional shrinkage ignores patterns shared across
regions, which might predict a target cell's profile from its other observations.

The estimand is the cell's retained regional **RNA allocation**, relative to
other-fold query cells and the unchanged gene-median reference contrast. It is
not an absolute DNA state or a genomic membership probability. The working
observation model is multinomial allocation of RNA counts conditional on the
retained library. Training cell covariance combines biological RNA variation
with any unmodelled noise. Transfer requires that this covariance supplies useful
prediction in other cells. The strongest competing biological explanation is
coordinated expression state or lineage rather than DNA dosage. Same-cell DNA
at callable matched resolution can test profile compatibility; independent
specimens and suitable controls are still needed for biological generalisation.

The fixed scONE fixture has 390 tumour cells, 450 annotated reference cells and
572 5-Mb/arm/gap regions. The five-region median can span a larger physical
interval; these results do not establish 5-Mb event resolution. Biological n=1.
RNA is rounded deposited expression, not verified raw molecules. A multinomial
model therefore supplies a descriptive working smoother, with no count-model
calibration claim. Processed corresponding-cell DNA is an uncertain reference
measurement, not error-free copy states or independently certified clone labels.
This exposed cohort remains regression-only.

## Executed estimator and matched controls

For training counts Y_i, observed library N_i and u_i=Y_i/N_i, calculate the
mean m and unbiased sample covariance S. The conditional multinomial noise
covariance has the unbiased estimate

`R_i = (diag(u_i) - u_i u_i') / (N_i - 1)`.

Thus `C_raw = S - mean(R_i)` estimates latent RNA composition covariance.
Its diagonal is algebraically the earlier factorial-minus-cross estimator
before truncation. For nonidentical fixed latent compositions it targets their
sample covariance; the population interpretation additionally requires a
representative independent training sample. Unequal N_i are retained explicitly.

Whiten C_raw by diag(m), remove the sqrt(m) composition direction, symmetrise,
and truncate negative eigenvalues to zero. This is a specified PSD approximation,
not an optimality or coverage guarantee. No rank, ridge or noise multiplier is
selected. The linear prediction is

`m + C (C + diag(m)/N)^(-1) (u - m)`.

An eigensystem implements it with factors `N*lambda/(1+N*lambda)`. On the
composition subspace, with C*1=0, the Poisson diagonal gives the same update as
the corresponding multinomial covariance. The fitted covariance remains a
working RNA object. A Euclidean simplex projection repairs negative linear
predictions; it is part of the fixed method, rather than outcome-specific
clipping. The diagonal ablation removes only off-diagonal entries of the same
PSD covariance before the linear rule. This also removes its compositional
off-diagonal constraint; simplex projection is applied to both arms.

Individual fractions and complete pooling are retained as additional controls.
All four actions become relative log responses through
`log((N*p_hat+1/2)/(N*m+1/2))`, then receive the unchanged other-fold gene-median
trunk, cell centering and five-region within-block median. The half count is a
common boundary regulariser, not a detection threshold or a calibrated prior.
The previous three folds, genes, references, region axes and DNA scorer are
unchanged. Every target fold is excluded from covariance fitting. Shared
upstream preprocessing remains conditioned upon; there is no claim of complete
independent discovery/confirmation separation.

All predictions, covariance fits and counts are locked before new DNA scoring.
All 24 preceding methods are retained with every score record exact. No
population or method disappears from the comparison.

## All new outcomes and the scientific decision

| Method | Whole Pearson | Cell-specific Pearson | Whole wins over pseudobulk |
|---|---:|---:|---:|
| Pseudobulk | 0.561258 | -0.048242 | Reference |
| Previous gene-NB moment candidate | 0.602558 | 0.308225 | 181/390, previous preserved result |
| Reference-refined REFCON | 0.615822 | 0.376007 | See preserved complete score table |
| New joint covariance | 0.567226 | 0.256930 | 156/390 |
| Matched diagonal covariance | 0.576824 | 0.301815 | 163/390 |
| Matched individual fraction | 0.551944 | 0.302780 | See complete paired table |
| Matched complete pooling | 0.561442 | -0.010205 | See complete paired table |

The joint rule loses to its diagonal ablation by mean **-0.009598/-0.044885**;
it wins in 163/390 whole-profile and 74/390 cell-specific comparisons. Relative
to the preceding gene-NB candidate, means decline **-0.035333/-0.051295**, with
75/390 and 63/390 wins. Against reference-refined REFCON it loses both means
and wins only 80/390 and 24/390 cells. There is no competitor-beating result.

The joint method's mean whole-profile gain of +0.005967 over pseudobulk coexists
with only 156/390 cell wins. Do not substitute the mean for the failed majority
criterion. Both the original anchored-candidate continuation gate and the
stronger previous-gene-NB gate fail for both new candidates. The joint-versus-
diagonal mechanism gate also fails. No biological-unit p-values are calculated.

The unprojected joint rule gives at least one negative regional fraction in
360/390 target cells. Maximum simplex displacement by fold is approximately
0.000726, 0.000793 and 0.000623 in fraction L2 units. The diagonal rule has no
negative entries, but its simplex projection still repairs total mass. The
three fitted covariance matrices retain 247, 247 and 246 positive eigenvalues
from 260 training cells. These are observed properties, not proof of a uniquely
identified failure cause. Covariance estimation error, RNA noise-law mismatch,
reference/trunk error and genuine expression state remain alternatives.

The matched ablation establishes that adding these off-diagonal estimates is
harmful here. It does not establish that all possible regional dependence models
are harmful. The fixed design prohibits a post-outcome rank, ridge, noise or
smoother sweep and prohibits expanding this failed candidate's evaluation.
Do not selectively combine its passing cells or endpoints with earlier methods.

## Independent computational evidence and its limits

Numerical prerequisites precede fitting. They include exact rational covariance,
12 enumerated unequal-library binomial configurations, independent fraction-unit
linear solves, zero-covariance pooling, cell/gene permutations and exact simplex
KKT fixtures. The exact unequal-library expected latent covariance is 1/8 on
the two diagonal entries and -1/8 off diagonal.

The independent reader imports no candidate estimator. It verifies:

- All regional counts reconstructed by explicit gene-to-region summation from
  the archived RNA matrix, with exact identity alignment and count conservation.
- All **9,360 preceding cell-score records** exactly preserved; 28 methods and
  390 cells remain, on **11,260,470 common gene/cell measurement pairs**.
- **147 exact Fraction covariance entries**, computed by factorial-minus-cross
  products instead of sample covariance minus noise. Maximum discrepancy
  **1.69e-21**, under the operation/magnitude diagnostic.
- **18 original-unit dense solves** on prespecified first/shallowest/deepest
  target cells, followed by exact Fraction simplex projection. Maximum observed
  prediction discrepancy **3.19e-16**. Covariance matrices in these coordinates
  can be ill-conditioned, up to 5.19e6; their worst conservative propagated
  allowance is 9.66e-4. The small observed disagreement is an agreement check,
  not a uniform high-precision error certificate or an interval calibration.
- Independent callable-DNA reconstruction and long-double whole/cell-specific
  Pearson for all four new methods, **3,120 cell/endpoint checks**. Maximum
  discrepancy **1.78e-14**. Every mean, paired comparison and decision is also
  recounted using exact sums of saved binary64 values. The 76 paired summaries
  retain every declared comparator and endpoint.

The eigen residuals are at most 2.24e-16, orthogonality residuals at most
5.65e-14, and composition residuals at most 3.71e-18. Arithmetic diagnostics
come from operation counts, magnitudes, condition numbers and independent
references. They are not scientific effect cutoffs or confidence intervals.
No reported performance gate is redefined by those diagnostics.

Computational correctness checks pass for the executed calculation. Statistical
identification of genomic variation and biological validity remain unestablished.
No false-split or power experiment is run for this failed profile candidate.
The public package/Python/Rust/GPU suites are not rerun: no public code changed.

## Every attempt, compute resource and preservation record

All three jobs use nslab account/QOS, CPU partition, one CPU, 8 GB, one-hour
limit and v4 Intel placement. They use the existing read-only Conda environment
at `micromamba/envs/allelic_qtl`, with Python 3.11.16 and NumPy 2.4.6 reported
inside execution. Threads are fixed to one; imports use `-I`, no user site,
no PYTHONPATH/PYTHONHOME and no bytecode writes. Transitive project source is
content-hashed. No environment or dependencies are installed or changed.

| Job | Work and progress | State / exit | Elapsed | Process peak RSS |
|---|---|---|---|---:|
| 23207068 | Numeric prerequisite, input gate 1/3 ██ | Failed 1:0 before fitting | 30 s | 365,916 KiB |
| 23207075 | Corrected coordinate comparison, 3/3 folds ████████ | Completed 0:0 | 27 s | 952,052 KiB |
| 23207085 | Independent reader, 3/3 folds ████████ | Completed 0:0 | 21 s | 885,224 KiB |

Two complete, one failed preserved, zero running or pending in this chain.
All dependent analysis completed. Short-job scheduler RSS is often zero;
process peaks above are measured separately and not inferred from that zero.

The first metadata launcher used unavailable `python` on the login host; it
submitted no job. Repeating those stdlib-only metadata operations with installed
`python3` succeeded. No scientific execution occurred on the login node.

The first compute job passed numerical prerequisites but compared raw CSR `.data`
arrays before/after integer conversion. Conversion can canonicalise sparse
storage order. That comparison failed before model fitting. The preserved
eight-file source snapshot identifies exactly what ran. The retry still requires
finite nonnegative integer input and now requires exact **matrix-coordinate**
equality after conversion, then exact regional count conservation. It completes
without changing inputs, model, scoring support or scientific criteria. This
fixes the check's identity estimand; it does not relax count equality.

All 91 prior validated package/native/data/test identities remain exact. All 283
preceding claim-ledger values are preserved; the new candidate is a separate
entry. The preceding arm-scoring correction, frozen outputs and original tests
remain unchanged. Unrelated jobs and installations are untouched. Full genomic
and SOTA claims still require genuinely stronger independent evidence.

## Reproducible source and artifacts

[Fixed scientific design](../dev/region_covariance_20261111/DESIGN.md),
[estimator and exact prerequisites](../dev/region_covariance_20261111/estimator.py),
[executed workflow](../dev/region_covariance_20261111/study.py),
[independent reader design](../dev/region_covariance_20261111/READER_DESIGN.md),
[reader implementation](../dev/region_covariance_20261111/reader.py),
[executed source manifest](../dev/region_covariance_20261111/source_manifest.json),
[failed-attempt identities](../dev/region_covariance_20261111/attempt1_manifest.json),
[reader source manifest](../dev/region_covariance_20261111/reader_manifest.json).

[First command and resource header](../dev/region_covariance_20261111/submission.json),
[preserved failed job](../dev/region_covariance_20261111/run_23207068/output.txt),
[corrected retry command](../dev/region_covariance_20261111/retry_submission.json),
[reader command](../dev/region_covariance_20261111/reader_submission.json),
[numerical prerequisites](../dev/region_covariance_20261111/run_23207075/NUMERICAL_PREREQUISITE.json),
[RNA lock before DNA](../dev/region_covariance_20261111/run_23207075/RNA_PREDICTIONS_LOCKED.json),
[all fitted diagnostics](../dev/region_covariance_20261111/run_23207075/COVARIANCE_DIAGNOSTICS.json),
[complete 28-method score report](../dev/region_covariance_20261111/run_23207075/REPORT.json),
[all cell records](../dev/region_covariance_20261111/run_23207075/scored_cells.tsv),
[prior-record exactness](../dev/region_covariance_20261111/run_23207075/PREVIOUS_EXACT.json),
[independent evidence and failed decisions](../dev/region_covariance_20261111/reader_23207085/results/SUMMARY.json),
[all paired effects](../dev/region_covariance_20261111/reader_23207085/results/paired_cells.tsv),
[accounting](../dev/region_covariance_20261111/accounting_complete.tsv),
[preservation verification](../dev/region_covariance_20261111/preservation_complete.json).
