# HaploDose scientific contract

Schema 2 separates correct computation, statistical identification and biological
validation. A completed run is not a validated biological result. The executable
claim ledger is `benchmarks/hardening/claims.json`; its diagnostic examples are
development investigations, not a sample for population-level inference.

## Outputs and limits

New profile rows preserve actual observed/expected bin pools and the centering
constant. Count-pool aggregation and mean clone log effects are distinct estimands;
do not interchange them or multiply the half-count offset by clone count. The
pools retain fractional fitted weights and conditional nuisance expectations, not
a newly calibrated count likelihood. [Exact export evidence](PROFILE_POOL_EXPORT_2026-11-09.md).

Gene-sum bins are working observations, not generally sufficient statistics
for heterogeneous gene-level NB counts. In the three original recovery failures,
known-parameter summed-count likelihoods retain substantial assignment information;
moment-matched NB changes few decisions relative to the actual sum law. This
does not validate estimated reference dispersions, nuisance or adaptive grouping,
and does not justify a gene-level replacement by itself. True nominal exposures,
rates and class priors in that diagnostic are unavailable oracle inputs.
[Matched representation and distribution checks](GENE_BIN_INFORMATION_2026-11-09.md).
Marginal class-decision agreement does not validate independence of thinned
proposal/confirmation counts or the complete adaptive null; unequal-gene sums
need not follow the NB law assumed by bin-level thinning.
An exact zero-covariance counterexample demonstrates this joint-law limitation.
Gene-native thinning supports independence only under fixed correctly specified
gene NB laws; full-count nuisance fitting and adaptive grouping do not inherit
it. A known-dispersion research variant improves original recovery2/5→3/5 but
leaves two failures and changes draw realizations. It is not a calibrated or
deployable default. [Law, unchanged gates and negative results](GENE_THINNING_2026-11-09.md).
Under the same fixed research gene splits, weighting both confirmation and merge
recovers4/5 original cases; either stage alone remains3/5.20261109 still loses
accepted distinctions at merging. This interaction does not calibrate fitted
membership weights or remove full-count nuisance leakage. The45 exact matched
fold comparisons exclude changed draws only within this stage ablation.
[Observed transition and unresolved criteria](GENE_THINNING_WEIGHTED_2026-11-09.md).
Aligning gene-veto pools with the actual depth/membership-weighted regional pools
restores5/5 original cases in that research variant. This supports a measurement-
population correction, not calibration of inherited sign-test rules on selected
weighted counts. The public procedure is unchanged; known gene noise, adaptive
selection, nuisance uncertainty and biological interpretation remain unresolved.
[Exact-pool evidence, detection decline and retained gates](GENE_VETO_ALIGNMENT_2026-11-09.md).
The same research veto fails weight-unit invariance: unscaled gene eligibility
and pseudocounts can change its decisions under common weight rescaling, even
when the conditional regional mean/variance are unchanged. The5/5 recovery result
does not justify promotion. [Exact unit falsification and three actual-pool changes](GENE_VETO_UNITS_2026-11-09.md).
Fresh full-workflow gene-count evaluation rejects its state-null applicability:
41/100 candidate false distinctions, simultaneous upper53.4487%, despite98/100
two-profile/chr12 recovery. Matched null0/100 upper4.6747% applies only to that
fixed known-noise regime. Contrast recovery is distinct from cell assignment,
whose all-cell mean is72.8%. [Unchanged5% gate, paired control and independent audit](GENE_VETO_FRESH_2026-11-09.md).

Effective weight units repair one research computational defect and preserve
five exposed recovery cases; they do not establish selected-gene calibration or
5% state-null control. [Unit repair evidence](GENE_VETO_UNIT_FIX_2026-11-09.md).
Hard-MAP eligibility can hide exactly identical expression profiles from merge
checks. Exact redundancy is a computational identification defect; unequal
state-driven profiles remain a scientific failure requiring separate evidence.
[Locked state mechanisms](GENE_STATE_TRACE_2026-11-09.md).

RNA effects describe normalized expression relative to recorded reference cells.
They do not identify absolute DNA copy number, cellular purity or balanced genome
doubling. Regional transcription and reference composition can reproduce a
dosage-shaped expression profile. Candidate groups are RNA-defined; only suitable
independent DNA supports their interpretation as CNA-defined clones. A similarity
dendrogram is not an inferred tumour phylogeny.

A freely fitted shared regional offset and a common tumour RNA effect are
observationally equivalent when the reference component has zero weight:
adding a vector to the offset and subtracting it from each tumour profile
preserves every tumour mean. Fitted reference-like probabilities alone are not
an independently verified normal anchor. Improving the conditional NB offset
objective fails the original recovery gate and does not identify genomic
dosage. [Forward equivalence, numerical checks and workflow falsification](OFFSET_OBJECTIVE_2026-11-09.md).

Conditional-likelihood diagnostics must score the profile update actually
applied, including the fixed normal dose and retained latent population. Raw
segmentation returns can differ from that transition. The corrected production
trace reproduces inference exactly and still finds NB Q decreases; pooled PLN
segmentation is not thereby an NB M-step. A higher NB Q alone does not establish
genomic recovery. [Executed-transition audit](PROFILE_TRANSITION_2026-11-09.md).

Conditioning profile expectations on the current shared RNA offset can match
a known-offset forward example while worsening whole-workflow recovery. The
tested component intervention is not a supported default or an NB M-step;
offsets fitted from the same RNA do not become independent genomic anchors.
[Numerical agreement and failed recovery](PROFILE_OFFSET_2026-11-09.md).

Independent normal anchoring can remove an effect shared with query cells, but
adds reference-sampling noise and requires normal-to-query transport. A conditional
count example reverses a genuine gain with mismatched normals while accurately
estimating RNA relative to those normals. Matched normals still cannot distinguish
observationally equivalent tumour transcription and genomic dosage. Fixing the
shared offset0 also fails all five original recovery cases; neither intervention
is a supported package default. [Normal-anchor and full-workflow evidence](REFERENCE_OFFSET_2026-11-09.md).

Normal-only predictive shrinkage can reduce calibration noise without using
query outcomes, but it does not establish accurate regional amplitudes. In a
fresh conditional experiment it lowers per-bin MSE while worsening regional
contrast MSE; scalar shrinkage chosen on half-sample fits is not automatically
appropriate for a full normal refit. Mismatched-normal and tumour-state ambiguity
also remain. [Selection/leakage audits and loss decomposition](NORMAL_OFFSET_SHRINK_2026-11-09.md).

Carrying a normal-only penalty through the full refit using scalar expected
information partly reduces that regional bias on fresh conditional counts, but
does not match full anchoring's regional accuracy and worsens mismatched-normal
errors. Its exact isotropic quadratic limit does not make it an exact multivariate
NB posterior. It is research evidence, not a supported general correction,
uncertainty estimate or clone recovery guarantee.
[Independent checks and fresh paired failures](NORMAL_OFFSET_INFORMATION_2026-11-09.md).

Depth slopes do not themselves identify a reference-population intercept. Current
query-centred depth correction does not transport the pooled reference depth
distribution. A known-slope forward example retains a whole-specimen depth effect
and can reverse a gain's modeled RNA contrast. Conditional intercept calibration
corrects that example, but estimated-slope and complete-workflow validation remain
open. Reference depth association can also reflect biological state or CNA content.
See the [depth-origin counterexample and limits](DEPTH_REFERENCE_ANCHOR_2026-11-09.md).

Estimated-slope count simulations reduce error with adequate reference depth
support but worsen it with narrow support or true zero slopes. A reference-specific
programme remains nontransportable even with known slopes. Fewer than60 fit cells
per reference group at default min_cells30 supply no quintile slope estimate;
returned zero slopes are not identified absence of a depth effect. No universal
intercept correction or gain-preservation claim is supported. Reference-resampling
bands also fail under absent depth support and reference-specific programmes;
sampling uncertainty does not establish transportability. See the
[estimated-slope falsification panel](DEPTH_REFERENCE_ESTIMATION_2026-11-09.md).

Continuous reference-depth regression with heldout-selected shrinkage improves
some conditional simulations but does not outperform pooling on the exposed
gene-level RNA counterfactual fixture or resolve reference-programme transport.
It is not a supported general correction or release default. See the
[regression and gene-aggregation checks](REFERENCE_DEPTH_SOFTMAX_2026-11-09.md).

Hard membership and latent-weighted dispersion moments have different population
targets. Correcting that ideal-model moment does not establish calibrated RNA
uncertainty or resolve state-driven genomic ambiguity: the tested weighted
prototype fails the original recovery and fresh state-null gates. See the
[dispersion/membership experiment](DISPERSION_MEMBERSHIP_2026-11-09.md).

`clone_posterior` is a plug-in responsibility for a fitted RNA component, not
validated probability of genomic clone membership. Pairwise coassignment computed
from responsibilities additionally conditions on fitted parameters and label
independence. In the exposed joint-mode research experiment, a95% pair-confidence
cutoff retained substantial genomic error under continuous state variation; it
is not a supported correction or release default. This research result is not
a measured biological error rate for the production workflow. See the
[membership-confidence diagnostic](CLONE_MEMBERSHIP_CONFIDENCE_2026-11-09.md).

`n_eff` is responsibility mass. `n_cells_map` counts hard assignments, and the
weighted support calculation uses the squared weight sum divided by the sum of
squared weights. Biological replication is counted at the specimen or donor level.
`cn_state` is a working RNA smoothing state. `malignant_call` is an exploratory
reference-difference flag. `tree.tsv` summarizes resolved RNA profile differences
by average linkage. The specimen's tree is withheld if a required pairwise
contrast is untestable. `qc.json` records each pair and `tree_status` under
`pairwise_depth_diagnostics`.

Depth matching balances coarse log-depth strata. Within-stratum depth and RNA
composition can still differ. `weighted_log_depth_wasserstein` measures the
remaining difference between matched log-depth distributions, in natural-log
units. It is undefined when there is no common stratum mass. It is a descriptive
balance diagnostic and does not enter split, merge, or significance criteria.
The gene sign veto uses unweighted hard assignments
in a proposal subsample. Its population can differ from the depth-matched regional
contrast. `gene_veto_uninformative` records when available gene inputs cannot
produce a usable contrast. Coordinated state can pass this heuristic.

Spatial neighbourhood smoothing is a directed, row-normalized heuristic, not
validated CNA confidence or a symmetric Potts-energy optimization. Adaptive
smoothing may select positive strength under spatially random finite samples;
no false-split guarantee follows. Explicit cell-ID tie handling preserves graph
identity under row permutation, but cannot establish genotype from spatial cell
state. See the [spatial correctness evidence](SPATIAL_IDENTITY_2026-11-09.md).
When metadata includes `section`, edges are restricted to individual sections.
Without section IDs, callers must supply a common coordinate frame per specimen.
Missing or nonfinite coordinates in a requested graph, or missing supplied
section IDs, warn and disable smoothing for that specimen. `qc.json` records
the reason. The assignment model scores competing RNA components. Mixed-spot
RNA contribution and cell fraction require a separate mixture model.

Allele imbalance is conditional on donor identity, phasing, bias and count-noise
assumptions. It does not by itself establish copy-neutral LOH. Allele-balanced
query cells are not established normals and are not added as references by default.
The expression-only selected-parent split filter cannot calibrate positive allele
evidence, so that combined configuration is rejected until its null is implemented.
The same restriction applies to gene-consistency vetoes: the current bin-only simulator
cannot reproduce gene evidence. The public gene-aware pipeline therefore rejects this
split filter. Bin-only selected-parent experiments remain diagnostic-only.

Conditional segment standard errors omit fitted-parameter uncertainty, genomic
correlation and selection unless their metadata explicitly say otherwise. A
three-standard-error diagnostic is not a minimum detectable DNA effect. The
current arm-score scale is descriptive and is not marked calibrated by panel size.
The approximate segment SE is a scale diagnostic for the unweighted observed
log-count-ratio statistic, exported separately from the HMM posterior mean state
effect. An observed statistic is not an unknown estimand; this diagnostic has no
validated parameter-interval target. It is not a confidence interval for the
posterior effect or for DNA copy number.
The preserved earlier simulations falsified nominal 95% coverage in four of six
prespecified fixed/selected interval conditions after Holm correction; the independent
follow-up falsified five of six using the same criteria. Selection and correlated
expression caused substantial undercoverage. These SEs must not be presented as
calibrated confidence intervals; see the dated evidence report for the estimand
and simultaneous coverage bounds.

The optional systematic RNA interval helper requires explicit `counts_kind='raw'`
and count provenance, independently fixed expected intensities, a known Gaussian
latent covariance and a frozen bin/gap universe.
It combines simultaneous Poisson count bounds with all-interval Gaussian bounds.
This is an explicit model contract, not an automatic correction for fitted
references or residual covariance. Provenance declarations cannot establish those
assumptions. The helper does not run automatically in the pipeline; its width,
power and applicability must be assessed separately. See the
[follow-up evidence](RECOVERY_2026-11-07.md).

## Measurement and role contracts

New runs require raw-count, genome-build, coordinate-system and feature-namespace
declarations. Internal genomic intervals are zero-based half-open. Numbat-style
SNP POS and exclusion-site keys retain their one-based convention, converted at
interval comparisons. Count integrality alone cannot prove provenance: normalized
values rounded to integers remain invalid count-model evidence.

Each dataset has independent feature selection and bins. Features are selected
from fit references; held-out references do not determine the feature universe.
Query/reference mismatch flags are diagnostics rather than automatic exclusion of
possible high-amplitude dosage signal. Explicit discovery and confirmation roles
are recorded separately. Disjoint features enforce computational separation but
do not prove independence of biologically correlated counts.
The [shared-state allele experiment](ALLELE_SHARED_STATE_2026-11-08.md) explicitly
falsifies genomic confirmation by reserved SNPs alone: a shared allele-expression
state can reproduce both discovery and reserved contrasts in a single DNA
population. If a genomic alternative has the same RNA/allele observation law,
its detection probability equals the false genomic-split probability for that
state alternative. Independent DNA or independently justified nuisance
constraints are necessary to distinguish those labels.

## Validation decisions

The clone shared-offset residual now matches retained responsibility mass in
observed and expected pools after component removal. Its legacy clone-only
mean remains a heuristic, excluding cell depth and programme terms. A broader
conditional-mean replacement passed exact arithmetic checks but failed primary
recovery and is unpromoted. A passing numerical check cannot override that
failure. The narrower accounting repair passes the full suite and preserves
the existing2/5 development recovery passes; the other three failures remain.
[Derivation, candidate rejection and validation](OFFSET_EXPECTATION_2026-11-09.md).

The clone target is a 5% chance of an unsupported additional distinction in an
independently established single-CNA population, at stated DNA resolution. Passing
requires a one-sided 95% upper confidence bound at or below 5%, together with
reported recovery of genuine contrasts. For one prespecified population, zero
errors still requires 59 independent biological specimens. Simulations, seeds,
cells and repeat subsamples do not add biological replication. Uncertain DNA truth
does not certify a negative control.

Former evaluation cohorts remain regression-only for the revised method. Existing
frozen artifacts and historical decisions are retained. Corrected scores and truth
receive new identities. There is currently no verified unopened clone panel
sufficient for the proposed biological guarantee; candidate outputs remain
exploratory until that evidence exists.
The frozen November 6 bin-level null experiments also fail under correlated expression
state: 91 of 100 single-CNA simulations produce extra groups. Conditional replay
does not correct this failure. This rejects broad error-control claims even
before the missing biological validation panel is considered.
The follow-up state-calibrated global gate reduced accepted distinctions to 8/200,
but its simultaneous upper bound was 8.96%, above the 5% target, with detection in
65/200 alternatives and no membership repair. It is not promoted to a default.

Numerical validation compares against independent float64 optimization,
quadrature and exhaustive HMM calculations. A backend disagreement close enough
to a decision boundary to change its meaning requires refinement or an unresolved
result. More compute cannot resolve observationally equivalent biological states.
