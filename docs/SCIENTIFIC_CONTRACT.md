# HaploDose scientific contract

Schema 2 separates correct computation, statistical identification and biological
validation. A completed run is not a validated biological result. The executable
claim ledger is `benchmarks/hardening/claims.json`; its diagnostic examples are
development investigations, not a sample for population-level inference.

## Outputs and limits

RNA effects describe normalized expression relative to recorded reference cells.
They do not identify absolute DNA copy number, cellular purity or balanced genome
doubling. Regional transcription and reference composition can reproduce a
dosage-shaped expression profile. Candidate groups are RNA-defined; only suitable
independent DNA supports their interpretation as CNA-defined clones. A similarity
dendrogram is not an inferred tumour phylogeny.

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

## Validation decisions

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
