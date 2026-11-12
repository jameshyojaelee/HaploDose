# RNA-selected segmentation: small profile gain, no justified replacement

The implemented research selector increases mean DNA-profile correlation on the
181-line Kinker regression from .680218 to .685905. Its added benefit over the
retained RNA selector is unresolved across shared-reference groups, and raw-bulk
correlation, loss recall and false-positive rate worsen. Do not promote this
addition. The simpler installed predictive profile remains opt-in; defaults,
frozen results, clone criteria and all prior comparators remain unchanged.

## Model and decisive comparison

The estimand is relative regional RNA abundance on the existing natural-log
scale. The observation criterion is held-out gene absolute prediction error,
conditional on existing feature selection, references, normalization and callable
support. This is neither independent confirmation nor a calibrated DNA model.
Regional transcriptional regulation and reference composition can produce the
same RNA pattern as DNA dosage. Corresponding independent DNA at matching
resolution distinguishes those explanations; no new biological specimens enter.

`rna_predictive_partition.py` retains the installed masked-gene median selector
exactly and adds one piecewise-constant Gaussian mean partition from the existing
dynamic program. Every fold excludes its held-out genes from medians, scale and
partition fitting. Singleton genes remain training-only. First-difference MAD
within genomic blocks estimates scale, divided by .6744897501960817 and sqrt(2).
Unidentified scale makes this candidate ineligible; there is no positive floor.

The fixed penalty is 3 log(block bins) per added segment, inherited from the
existing MDL comparator. Counting a mean and breakpoint gives an approximate
coding rationale on the doubled likelihood scale; this is not an error-control
theorem. The partition is selected only when its masked absolute error, plus
two operation-based rounding allowances, is below every retained median width.
Ties retain the original predictor. Full profiles are then refitted using RNA.
There was no DNA-dependent choice of width, penalty or eligibility.

Before scoring DNA, each panel saves predictions, choices and hashes. Kinker
retains all four previous comparators; bulk retains all twelve originals plus
the three installed gene-profile comparators and this candidate. Common support,
the historical 24/24 tumour calibration/comparison split, truth and thresholds
remain fixed. Missing candidate predictions on informative truth fail explicitly.
Former evaluations remain regression-only.

## Checks and complete results

Numerical checks independently enumerate all eight partitions of a four-value
fixture and verify the optimum, a noiseless median optimum, constant-scale
abstention, and invariance of fold predictions when only held-out values change.
Every retained RNA prediction matches its frozen array exactly for all181 Kinker
lines and all68 bulk samples. The uncertainty job independently verifies all724
old Kinker comparator scores and the first100 ratio-bootstrap estimates using
`math.fsum`. These checks execute on compute nodes.

The two locked diagnostic lines retain the original predictor exactly; their
mean correlation remains .793928. This is a mechanism check, not biological n=2
independent validation.

| Kinker comparison, n=181 lines | Mean paired correlation change | Wins/losses/ties |
|---|---:|---:|
| Retained RNA selector | +.005686 | 33/30/118 |
| Fixed sm5 pseudobulk | +.017949 | 102/47/32 |
| Standalone MDL | +.017543 | 84/34/63 |
| Unsmoothed pseudobulk | +.089965 | 170/11/0 |

The partition is selected for63/181 lines. Mean correlations are .685905 new,
.680218 retained, .667955 fixed sm5, .668362 standalone MDL and .595940 raw
pseudobulk. These point estimates do not establish benefit beyond the retained
selector or independent superiority.

Post-result uncertainty characterization freezes the eight previous shared-query/
reference connected groups, sizes18,21,24,24,24,24,24,22. With seed2026110912,
20,000 group bootstrap draws preserve the line-weighted estimand by taking the
ratio of sampled group sums to sampled group sizes. Two contrasts receive
Bonferroni 97.5% marginal percentile intervals. These are approximate conditional
intervals assuming independent groups and fixed preprocessing, not validated
coverage or a biological calibration guarantee. Four Monte Carlo endpoint
brackets receive two-sided alpha .0125 each.

Against the retained selector the interval is [-.001629,.013435], including zero;
the eight group effects include -.009977 and +.023442. Against fixed sm5 the
interval is [.012502,.026416]. Monte Carlo endpoint brackets preserve these
interpretations. The new component's additional benefit is not established.

| Bulk comparison against retained selector | Paired change | Support | Wins/losses/ties |
|---|---:|---:|---:|
| Correlation | -.000920 | 23 comparison patients | 2/2/19 |
| Loss recall | -.014652 | 13 with loss truth | 1/1/11 |
| Loss precision | +.010683 | 16 with both methods defined | 2/2/12 |
| Loss false-positive rate | +.003012, worse | 23 with callable negatives | 3/2/18 |

FPR wins above mean lower error. Precision is conditional on both procedures
making calls; its small gain is not an unconditional benefit. The partition is
selected in12/68 total tumour/normal samples. Against strong naive segmentation,
correlation changes -.050936 and recall -.102839; against the existing
precision-weighted HaploDose partition mean, changes are -.053743 and -.280266.
All per-patient scores, undefined outputs and comparator denominators are saved.
No threshold or comparator is removed to improve the verdict.

## Decision, execution and preservation

Reject a bulk/shared replacement. Keep the Kinker point gain as a negative or
unresolved component experiment, without adding it to production. Current evidence
supports retaining the simpler installed selector. There is no new clone,
uncertainty-calibration, raw-count calibration or SOTA claim.

| Job | Work | CPUs | RAM | Limit | Outcome |
|---|---|---:|---:|---:|---|
|23141126|independent numerical checks|1|4GB|1h|complete0:0|
|23141127|locked two cases|4|16GB|2h|complete0:0|
|23141128|all181 Kinker lines|4|16GB|2h|complete0:0|
|23141129|raw bulk and all comparators|4|16GB|2h|complete0:0|
|23141957|shared-reference uncertainty|1|4GB|1h|complete0:0|

CPU allocations use innovation/nslab without arrays or interactive QOS. Source_v1
and source_v2 are immutable snapshots, running read-only against the noneditable
validated wheel23138003. Baseline/source manifests, input identifiers, environment,
backend, dependency status, process resources and scheduler sampling are retained.
Zero sampled RSS on short jobs indicates missing sampling, not zero memory use.
The prior412 Python/4 Rust passes and1 GPU-only skip are prior installed-artifact
evidence; full package/GPU suites were not rerun for these research-only files.

[Design and exact outputs](../dev/predictive_partition_20261109/DESIGN.md),
[uncertainty design](../dev/predictive_partition_20261109/UNCERTAINTY_DESIGN.md),
[complete uncertainty report](../dev/predictive_partition_20261109/uncertainty_v1/report.json).
