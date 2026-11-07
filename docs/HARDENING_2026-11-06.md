# Correctness and scientific evidence, 2026-11-06

The scientific goal is to make HaploDose conclusions depend on identifiable measurements,
correct computation and independent evidence at the claimed resolution. Correctness changes
precede equally prioritized profile/bulk and clone investigations. Unsupported biological
interpretations, calibration guarantees and sensitivity claims are withheld.

This is an implementation and evidence record, not publication readiness. The executable claim
contract is `benchmarks/hardening/claims.json`; public interpretation and migration are documented
in `SCIENTIFIC_CONTRACT.md` and `MIGRATION_0.2.md`. Historical result files and decisions are preserved.

## Frozen execution chain

The initial dirty source, including untracked implementation files, was copied before edits to
`dev/hardening_20261106_2119/baseline`. The adjacent baseline manifest, Git status and patches
identify the initial state. Shared installations were not modified. A copied, isolated Conda/Python
environment builds noneditable wheels on compute nodes. Validation executes from a neutral directory
and records the imported Python package, native binary, wheel hash, dependencies, exit status,
runtime and memory. Source snapshots have per-file hashes and a completeness marker.

Selected frozen result hashes are checked before and after the measurement audit. This is a
bounded integrity check of the listed artifacts, not a claim that every byte in all historical
directories has been rehashed. The audit recorded 777 source/measurement/provenance files and eight
diagnostic panels in `dev/hardening_20261106_2119/evidence_v2` (job 22812733). Its four linked records
describe collection, not an inference run. Actual inference/scoring jobs have separate records.

## Correctness evidence

| Check | Preserved baseline | Revised candidate evidence |
|---|---|---|
| Full build/tests | Job 22811740: Rust 3 passed; Python 99 passed, 1 skipped | Candidate v5, job 22814288: Rust 4 passed; Python 280 passed, 1 skipped, 1 unchanged thinning-recovery failure |
| Constant-reference programme rank | Four directions despite zero variation | Zero directions |
| Responsibility-weighted programme objective | Gap to independent SciPy optimum 0.0029103 | Gap 0.000000723; arithmetic comparison budget 0.0007496 |
| Fractional observations | Accepted 1.0001 as a count | Rejected |
| Duplicate feature identities | Accepted | Rejected before indexing |
| Integer aggregation | 2^24 + 1 became 2^24 | Rejected as unrepresentable; exactly representable sums preserved |
| Discovery isolation | No independent confirmation execution contract | Perturbation and unrelated-dataset tests pass; frozen confirmation artifacts unchanged |
| Allele counts | Fractional/nonfinite/overflowing observations could pass or truncate | Exact finite integer validation, bounded exact aggregation, and coordinate-mask tests pass |
| Empty outputs | Headerless tables could not be read | Primary schema columns emitted even when empty |

Before/after counterexamples use the same runner hash and independent expected results:
`baseline_counterexample_results_v2/results.json` and `counterexamples_v2/results.json` under the
task directory. The first baseline probe failed only at JSON serialization; that failed job and
its files are retained. It was corrected and rerun under a new identity.

The latest build evidence is the [v5 Python test report](../dev/hardening_20261106_2119/validation_22814288/pytest.txt),
[Rust test report](../dev/hardening_20261106_2119/validation_22814288/cargo_test.txt) and
[installed-artifact record](../dev/hardening_20261106_2119/validation_22814288/installed_artifact.json).
The formal identification, uncertainty and backend measurements below used the separately
recorded v4 wheel and source snapshot. They are not presented as numerical reruns of v5.

Other repairs cover SVD rank after filtering, responsibility-weighted nuisance optimization,
compositional overlap, segment-specific variance and RNA-effect units, genomic discontinuities,
sample/gene alignment, fixed comparator families, callable support, paired estimands, and complete
or failed execution manifests. New inputs must declare measurement scale, assembly, coordinates
and feature namespace. Input declarations remain assumptions that require source verification.

The public gene-aware pipeline rejects the selected-parent split filter because its bin-only
null omits the observed gene veto. Allele-enabled null filtering is also rejected. Historical
recovery assertions for these withdrawn configurations remain as named benchmark functions and
in the immutable original test snapshot; they were replaced by rejection tests for the new API.
Their earlier failures have not been converted into scientific successes.

The issue-to-check map is separate from biological accuracy:

| Correctness group | Implementation and decisive executable evidence |
|---|---|
| IDs, coordinates, exact counts | `io.py`, `bins.py`, `alleles.py`, Rust/native aggregation; `test_hardening_inputs`, `test_hardening_allele_inputs`, `test_native`, and identity fixtures in `test_hardening_truth` |
| Measurement scale | Required public input declarations and exact integer validation; fractional and undeclared inputs fail before inference; historical source/wheel preserved |
| Dataset/confirmation leakage | `datasets.py`, `confirmation.py`, public pipeline; unrelated-dataset and confirmation-only perturbations in `test_hardening_pipeline`, `test_confirmation`, `test_hardening_confirmation_stability` |
| Numerical models | Rank, weighted NB objective and compositional diagnostics; independent SciPy optimizer, quadrature and exhaustive HMM references in model/backend studies and `test_hardening_model`, `test_hardening_mixture` |
| Clone/null mismatch | Unsupported observed-gene-veto or positive-allele/null combinations rejected; `test_hardening_null_contract`, `test_hardening_clones`; conditional bin-level search remains a limited experiment |
| Truth and scoring | Gap-aware DNA reconstruction, identifier/callable/fixed-family checks, paired estimands; `test_hardening_truth`, `test_hardening_scoring`; corrected scorer siblings preserve read-only historical scorers |
| Execution/provenance | Noneditable isolated wheel, transitive hashes, failed-run manifests, strict backend and exit status; `test_provenance`, `test_device`, compute execution records |

Test names above refer to modules under `tests/`, not assertions that biological claims passed.
Schema 2 preserves historical readers/results, states RNA units and evidence status, and removes
misleading calibration/detection-limit names. No production dependency was added. No branch,
commit, push or publication was performed.

## Scientific results and decisions

**Reference and state identification.** Independent calculations reproduce equivalence of
normalized RNA under constant log-rate shifts, balanced doubling, a purity/copy-number tradeoff,
and regional transcription versus dosage. Equal likelihoods for these examples establish
nonidentifiability in the constructed models, not the cause of any particular specimen's profile.
Fixing a reference from independent metadata can identify a contrast that a free reference mixture
absorbs. Retrospective reference selection using DNA agreement is not permitted.

The first independent mixture comparison also exposed a computational defect: the legacy
500-step reference-weight iteration stopped 75.55 log-likelihood units below the independently
computed optimum. The corrected optimizer reached the independently computed boundary weight
of one in two iterations, with zero recorded objective gap (comparison arithmetic budget
1.61e-7). Its normalized objective-gap and KKT certificates are recorded in the
[v4 identification diagnostic](../dev/hardening_20261106_2119/identifiability_v4/summary.json), job 22813530.
With an independently fixed reference, the constructed regional RNA effect is 0.2433 and every
regional bin is labelled gain; the freely optimized reference makes the effect approximately
zero and every regional bin neutral. Correct optimization therefore exposes reference/dosage
confounding; it does not establish improved biological accuracy. All five tested undeclared,
declared-normalized, fractional, negative or zero-library input cases were rejected. Integer-valued
normalized expression can still mimic raw counts numerically, so provenance declarations
cannot replace source verification.

**Thinning recovery.** Candidates v2 and v5 miss the known simulated chr12 subclone in the original
thinning fixture. The assertion is unchanged and still fails. Four diagnostic runs (job 22812467)
held the data, parameters and assertion fixed: current preprocessing failed; restoring historical
feature selection, historical binning, or both recovered the contrast. This is one sensitivity
counterexample, not an estimate of population power. The correctness fixes remain; the previous
thinning-recovery claim is not promoted. The default even/odd mode passed its existing recovery
checks, which likewise do not establish biological clone validation.

**Conditional uncertainty.** Job 22813531 generated 1,000 independent simulations per condition.
The primary estimand is the unweighted observed log-count-ratio mean against a known reference,
evaluated using the reported approximate SE. It is not DNA copy number or confidence coverage
of the exported HMM posterior mean effect. Fixed intervals were prespecified; selected intervals
were the largest absolute observed mean among actual HMM runs, including neutral runs.

| Generator | Interval | Nominal-95% coverage | Simultaneous 95% interval | Holm undercoverage p |
|---|---|---:|---|---:|
| Independent, homogeneous | Fixed | 945/1000 (94.5%) | 92.3–96.2% | 0.506 |
| Independent, homogeneous | Selected | 807/1000 (80.7%) | 77.2–83.9% | 1.51e-57 |
| Independent, heterogeneous | Fixed | 973/1000 (97.3%) | 95.6–98.5% | 0.9999 |
| Independent, heterogeneous | Selected | 900/1000 (90.0%) | 87.2–92.3% | 2.52e-10 |
| Correlated, heterogeneous | Fixed | 690/1000 (69.0%) | 65.0–72.8% | 2.35e-151 |
| Correlated, heterogeneous | Selected | 48/1000 (4.8%) | 3.2–6.9% | log10(p) = -1156.46 |

The predeclared family contains exactly these six one-sided exact binomial undercoverage tests
against coverage 0.95, with Holm adjustment at family alpha 0.05. The intervals are two-sided
Clopper–Pearson bounds with Bonferroni simultaneous family confidence 0.95. Four tests falsify
nominal coverage in these constructions; nonrejection in the other two does not establish
calibration. Other recorded endpoints have nominal descriptive intervals and are outside this
test family. The saved [primary log-tail table](../dev/hardening_20261106_2119/uncertainty_v4_logtails/primary_logtails.tsv)
and [reporting revision](../dev/hardening_20261106_2119/uncertainty_v4_logtails/corrected_summary.json)
retain the same simulation counts and simultaneous bounds. The reporting-only correction
replaces an underflowed numeric p-value with its log probability; it reruns no inference.

Oracle log-count intervals use the known log-intensity covariance and a Poisson delta
approximation. Their fixed-interval coverage is 95.0%, 94.8% and 94.7% in these conditions;
selected-interval coverage falls to 82.4%, 80.8% and 62.0%. Thus selection
remains consequential even with the generator's covariance. The retained conservative variance
mode gives a mean fixed-interval SE 15.5% above the oracle in independent heterogeneous data,
but 47.0% below it with correlated offsets. Correlated nulls produce at least one resolved HMM
event in 323/1000 simulations (nominal descriptive 95% interval 29.4–35.3%); this is not a clone
false-split rate. These are conditional simulations, not biological replication or clinical
error estimates. The [original simulation summary](../dev/hardening_20261106_2119/uncertainty_v4/summary.json)
and saved counts are preserved. Neither variance subtraction nor wider intervals alone has
been shown to solve both selection and correlated expression.

**Clone error and power.** The [joint summary](../dev/hardening_20261106_2119/clone_summary_v8/summary.json)
combines three procedures (unfiltered, selected-parent filter, conditional complete-search replay)
and three conditions. Its single one-sided 95% Bonferroni family contains six null upper bounds
and three alternative detection lower bounds. Each condition has 100 independent simulation
seeds, 600 cells and 240 bins; there are no biological replicates.

| Condition | Events, each procedure | Joint family bound | Interpretation |
|---|---:|---:|---|
| One CNA population under the NB working generator | 0/100 extra groups | Upper 5.0604% | Even this bound exceeds 5%; the narrower null-only bound was 4.6747% |
| One CNA population with correlated regional state | 91/100 extra groups | Upper 96.7346% | Broad error-control claim falsified |
| Two CNA populations with the same correlated state | 100/100 multigroup detections | Lower 94.9396% | Detection is not correct clone recovery |

In the third condition, median adjusted Rand index against simulated CNA labels is 0.0759
(mean 0.0785; interquartile range 0.0680–0.0865). Only 2/100 runs return exactly two counted
groups; 60 return three and 38 return four. The selected filter changes neither groups nor ARI;
complete-search replay changes no detection decision. Its output is a global gate, not a new
membership solution or calibrated abstention. These results reject the proposed correction's
broader applicability; extra simulation draws do not repair the missing biological assumptions.

The replay repeats bin-level initialization, estimation, proposal/merge/split/stopping operations
conditional on a known reference. It omits the public pipeline's raw feature selection, gene veto,
alleles and spatial graph, so it is not validation of the entire raw-data workflow. The alternative
places CNA and transcriptional state on the same region, deliberately challenging identification.
An earlier clean two-CNA diagnostic detected groups in 67/100 cases; its source/results remain
preserved but it cannot substitute for the same-nuisance alternative. No parameters were changed
to improve these outcomes. A calibrated complete public-workflow procedure remains unestablished.

**Numerical backends.** The formal CPU/L40S comparison, job 22813539, verified the installed v4
wheel against its package files and checked device placement at outputs inside PLN, HMM,
cell-likelihood, programme and profile inference. The
[comparison](../dev/hardening_20261106_2119/backend_gpu_v4/comparison.json) records no changed
profile Viterbi bin, posterior call at 0.95, resolved call or cell MAP assignment. The largest
tested CPU/GPU cell-likelihood difference was 0.000331; HMM posterior differences were below
1.8e-7. Against the input-matched float64 NB reference, maximum disagreement was 2.26e-5 on CPU
and 3.53e-4 on GPU, while reference assignment margins were at least 25.3. These are clear-margin
examples, not a test of arbitrary near-boundary decisions. Independent SciPy quadrature,
optimization, likelihood and exhaustive path results are recorded in the separate
[CPU](../dev/hardening_20261106_2119/backend_gpu_v4/cpu/summary.json) and
[GPU](../dev/hardening_20261106_2119/backend_gpu_v4/gpu/summary.json) summaries.

The GPU summary records JAX allocator peak usage of 134,235,648 bytes (128.02 MiB) and a peak
pool of 130 MiB. These are allocator measurements, not the entire device/process memory peak.
Host peak RSS was 595.27 MiB for the CPU worker and 926.05 MiB for the GPU worker. The combined
two-worker command took 13.97 seconds; this timing is not evidence for GPU speedup.

**Measurement chain.** U-CAN verification matched 2,259,684 values from 120 normal donors and the
two locked tumours against source featureCounts files. Quantifier settings are recorded; original
runtime versions were not saved by the old script. The first downstream reporting attempt failed
on mixed-case versus canonical gene symbols after successful inference; its files are retained.
The matched raw/transformed TCGA comparison is a development diagnostic. Inverted batch-adjusted
RSEM and reconstructed Kinker integers are not evidence for raw-count uncertainty calibration.

The lexical allele selection is `2313287_STOMACH`, with 22 phased chromosomes and two pileup
channels containing 54 and 50 barcodes. The [allele trace](../dev/hardening_20261106_2119/allele_trace_v4/summary.json)
reconciles all 24,576 final AD/DP/GT rows, 104 cells and 3,500 phased sites to the source channel
matrices and donor-specific filtering. No barcode, orientation or row-total discrepancy was found.
This verifies transformation integrity, not donor identity or phasing accuracy.

A real assembly discrepancy justified bounded read inspection: the final VCF header advertised
GRCh38, whereas the executed pipeline used hg19. Both BAM headers have GRCh37 chromosome lengths;
the variant panel declares b37/hs37d5, and both prespecified indexed loci match its REF/ALT.
The [installed Numbat template capture](../dev/hardening_20261106_2119/allele_template_v1/commands.json)
shows `make_vcf_chr` replacing VCF metadata with a fixed GRCh38 template while retaining coordinates.
This explains the false header; it does not independently validate every coordinate or gene/cM
annotation. Historical VCFs remain untouched. New allele-enabled runs require a declared assembly
matching the expression annotation; the declaration must be supported by the execution chain.

Andor DNA reconstruction respects gaps and physical support. In SNU-601, the corrected procedure
changes 15 inferred DNA groups to 16 on the same 1,527 passing cells (partition ARI 0.821), and
the union of differentiated arms from 19 to 39. These are inferred reference measurements.
One simulated low-depth/bias challenge changes 16 groups to six, with ARI 0.226 on 1,473 shared
cells; 10/39 differentiated arms remain. This is not a calibrated clone-event sensitivity curve:
the challenge uses a cancer-derived bias/reference and changes depth plus inference conditions.
MKN-45 and NUGC-4 retain one inferred DNA group, but median fragment counts of 3,220 and 5,665
cannot certify single-CNA biological controls. Full SNU-601 has median 236,973 fragments.

Tracing Andor RNA identities found PRICKLE4 mapped to the wrong raw Ensembl row in five reference
libraries. A unique ID at the exact current-GTF locus identifies the correction without DNA
outcomes. The [new derivative certificate](../dev/hardening_20261106_2119/andor_derivative_v7/input/derivative.json)
records 136 changed reference entries across the three prepared cases (130 changes from 1 to 0,
six from 2 to 0). Every other count and all metadata remain unchanged; all 16,014 selected rows
were checked against raw matrices. These rows overlap across case pools and are not independent
biological replicates. Current-GTF consistency does not recover the historical alignment annotation.
Focused derivative/truth fixtures passed 16/16 in compute job 22814832.

The [corrected three-case RNA analysis](../dev/hardening_20261106_2119/andor_rna_v7/summary.json)
completed with fixed settings (job 22814881). SNU-601 has 12 retained RNA groups versus 16 inferred
DNA groups, MKN-45 one versus one, and NUGC-4 six versus one. Retained groups have at least
`max(20, 2% of query cells)` members; this is a descriptive reporting filter. SNU-601 has DNA
differences on all 39 supported arms, and the other two have none, so each case has only one
truth class and arm AUROC is undefined. RNA/DNA group counts do not establish corresponding
memberships. NUGC-4 subdivisions remain unresolved state-versus-genome contrasts because shallow
DNA cannot certify the null; MKN-45 agreement is descriptive, not a specificity estimate.

**Profile and bulk development.** The [TCGA development analysis](../dev/hardening_20261106_2119/profile_v6/summary.json)
uses 24 calibration patients, 24 reserved development patients and 20 reference donors, with eight
normal-reference resampling refits. Genes, callable windows, smoothing and comparator families
are fixed; 2,130 ambiguous GRCh37 annotation loci are excluded before inference. The unsmoothed
HaploDose ratio does not improve mean Spearman agreement over naive ratios (difference -0.0069).
Applying the same smoothing to both gives a mean difference of 0.0410; the family-adjusted
basic-bootstrap interval is 0.0033 to 0.0755. HMM/export output exceeds the smoothed naive baseline
by 0.0694 (0.0246 to 0.1312). There are 23 evaluable correlation pairs and 13 loss-recall pairs;
missing/undefined endpoints remain reported. None of the six loss-recall contrasts passes Holm
adjustment. These are internal development results, with a finite reference-resampling ensemble;
the interval coverage itself is not validated. They do not establish prospective superiority.

U-CAN raw-source verification and inference complete for both locked tumours. On the same callable
bins, naive five-bin-smoothed / HaploDose five-bin-smoothed / HMM / exported profile Spearman
correlations are 0.629 / 0.715 / 0.817 / 0.809 for U0001 and 0.565 / 0.604 / 0.641 / 0.624 for U0008.
Export smoothing slightly reduces agreement in both. These two contrasts are descriptive and
do not identify the cause of performance in other specimens. Comparison with a rolling filter
alone does not isolate the benefit of an adaptive HMM segmentation operator.

The [locked artifact diagnostic](../dev/hardening_20261106_2119/artifact_diagnostics_v8/summary.json)
applies identical filters and identical frozen-HaploDose-boundary mean projection to both methods.
Boundaries are selected using HaploDose; this is a conditional operator control, not independently
learned naive segmentation. For Kinker `2313287_STOMACH`, pseudobulk Spearman correlation rises
from 0.577 to 0.729 under the fixed five-bin filter, versus 0.506 to 0.667 for the HaploDose gene
ratio. Posterior correlation 0.660 falls to 0.408 after exported segment averaging. For
`42MGBA_CENTRAL_NERVOUS_SYSTEM`, identical boundary projection gives pseudobulk/HaploDose-ratio
correlations 0.806/0.802. These observations reject attributing those gains uniquely to the
HaploDose signal model. They do not identify the cause in other lines.

For transformed-input `TCGA-05-4389-01`, the median signal on 15 callable focal-loss windows is
-0.3683 in the unsegmented HaploDose ratio and -0.00652 in exported dosage. The earlier recall
values 0.50/0.05 remain historical measurements with different callable support; they are not
recomputed as though the denominators were unchanged. Saved per-bin posteriors are absent for
that bulk run, so HMM inference and export averaging cannot be separated. The result establishes
signal compression in their combined stage, not raw-count uncertainty validity.

On the reserved raw TCGA development patients, achieved patient-mean neutral-window loss FPR is
0.0561 for HMM/export, 0.0505 for the unsegmented ratio, 0.0429 for its five-bin smoother, 0.0417
for naive ratios and 0.0400 for smoothed naive ratios. FPR is defined for 23/24 patients and loss
recall for 13/24. Precision averages have method-dependent no-call exclusions and are descriptive.
All methods cover the common classified callable truth. These point estimates do not certify
5% FPR, and neutral-window FPR is distinct from whole-specimen clone false-split probability.
With the same donors and genes, raw-minus-transformed naive loss recall differs by -0.278
(smoothed: -0.172) despite similar FPR. Measurement-scale effects are therefore material.

## Gates that computation cannot supply

The clone target is an unsupported additional CNA distinction within an independently established
single-CNA population at documented DNA resolution. A one-sided exact 95% upper confidence bound
must be at most 5%, with power and abstention reported simultaneously. With no errors and one
prespecified population, the upper bound is `1 - 0.05^(1/n)`, requiring 59 independent biological
specimens. More strata require simultaneous bounds. Simulations, cells, reads, chromosomes and
repeat subsamples are not additional biological specimens.

No suitable unopened biological panel meeting that requirement has been verified. Same-cell DNA/RNA
would support membership validation; separate-cell DNA supports profile compatibility and sampled
clone diversity. Source/donor correspondence, DNA sensitivity and biological replication remain
load-bearing evidence gaps where unverified. Former evaluation datasets are regression-only.
The metadata-only SPECTRUM gate found zero verified same-biopsy correspondences among 148 candidate
RNA/DNA site pairs: 124 unresolved sites, 13 matching site names only and 11 missing sites. These
pairing records are not 148 independent specimens or proof that correspondence cannot be resolved.
Source collection metadata must establish correspondence before a validation claim.

Scientific development results must compare identical support, references and smoothing, retain
simple baselines and failures, and use biological-unit paired effects. Conditional standard errors
remain approximate; a completed workflow, a passing test suite, or comparator agreement does not
establish calibrated uncertainty or biological validity.

## Implementation decision and reproducibility

The deliverable is a corrected exploratory 0.2.0/schema-2 candidate, its isolated wheel and the
reproducible falsification suite. The full test suite retains one scientific recovery failure;
the biological clone guarantee and calibrated intervals have failed their present gates. No new
variance rule, richer mixture or spatial default is promoted. Positive development profile results
are retained alongside negative smoothing, threshold-transfer and clone findings. Addressing the
remaining claims requires a justified observation/selection model and suitable independent truth,
not just more repetitions of the same calculations.

The [job registry](../dev/hardening_20261106_2119/job_registry.json) records actual scheduler states,
exit codes, resource measurements and retained failed attempts. Builds/tests/numerical studies and
analyses ran on CPU, GPU and I/O compute nodes under innovation and nslab; no arrays or interactive
QOS were used. No job needed more than 210 GB, so bigmem was not requested. Per-job execution files
bind exact commands, source snapshots, installed binaries and environment; the wheel build/import
check is distinguished from the failing full test suite.

The [final integrity audit](../dev/hardening_20261106_2119/final_audit/summary.json) and
[changed-file manifest](../dev/hardening_20261106_2119/final_audit/changes.json) compare against the
initial dirty source, rather than attributing all pre-existing Git changes to this task. The
adjacent implementation patch is reviewable independently of the initial uncommitted work.
