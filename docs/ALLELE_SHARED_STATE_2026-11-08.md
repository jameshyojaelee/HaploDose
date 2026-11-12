# Reserved SNPs reproduce allele-expression state without identifying DNA clones

The physically isolated allele workflow produces positive reserved contrasts in
40/40 shared-state simulations under each existing cell model. Every simulated
specimen has one DNA population by construction. An independent audit reconstructs
all 240 fits and 240000 permutation statistics exactly. This rejects the proposed
interpretation that independently sampled reserved reads alone establish genomic
confirmation. It does not measure the biological prevalence of this countermodel
or the production whole-specimen false-split rate.

## Observation law and competing explanations

The estimand is the probability of a positive candidate held-out allele contrast
after complete discovery, under three explicit RNA observation laws. Each has
120 query cells, 40 RNA bins, 40 phased SNPs, DP20 per SNP, known reference and
fixed RNA dispersion .15. RNA contains the same regional shift in every cell and
no class difference. Allele fractions are either constant .65, independently
chosen .5/.8 per cell/site, or shared .5/.8 per cell across sites. The latter has
60 cells in each expression state, randomly assigned before sampling reads.

All three conditions have expected pooled allele fraction .65. Independent read
draws do not remove the shared latent cell state. The same joint RNA/count law
can be labelled two genomic populations with those allele fractions. Thus the
biological labels have identical likelihoods for every observation; numerical
separation or accurate latent-state assignment cannot select the correct label.
The shared ASE law is an explicit falsification alternative. Its extent and
frequency in real specimens remain [Unverified]. Same-cell DNA copy/allele
measurements would distinguish membership; separate-cell DNA can constrain
profile compatibility at its assay resolution.

## Frozen workflow and working test

Forty independent simulation seeds per condition/model yield 240 full fits.
Metadata-first whole-site hashing selects 13 discovery and 27 reserved SNPs before
counts are inspected. Only discovery observations enter initialization, nuisance
fitting, proposals, splitting, merging and stopping. Existing conditional support
rules remain unchanged, including the 30-cell minimum. This binned experiment
omits raw gene filtering/vetoes, spatial smoothing and reference uncertainty.
The validated package and current unsupported joint-null rejection are unchanged.

After discovery artifacts are frozen, the working statistic takes the maximum
absolute difference in mean per-cell reserved allele fraction across eligible
group pairs and discovery-selected regions. A whole-cell row permutation preserves
within-row SNP dependence; 999 permutations give p=(1+exceedances)/1000, with
inclusive ties. The working .05 cutoff was fixed before execution. Missing
eligible pairs or regions yield abstention and p1. This is a research diagnostic,
not a production clone acceptance rule. Its exchangeability requirement fails
when reserved observations share the state that formed discovery memberships.

## Error, availability and assignment together

| Cell model | Allele condition | Available contrasts | Positive contrasts | Mean latent all-cell matching |
|---|---|---:|---:|---:|
| Binomial | Constant |0/40|0/40|50.00%|
| Binomial | Site-independent state |0/40|0/40|51.17%|
| Binomial | Shared-cell state |40/40|40/40|100.00%|
| Concentration 2 | Constant |0/40|0/40|50.02%|
| Concentration 2 | Site-independent state |0/40|0/40|50.04%|
| Concentration 2 | Shared-cell state |40/40|40/40|85.19%|

All six endpoints share Bonferroni6 simultaneous one-sided 95% exact binomial
bounds. Forty positives have lower bound .8871978258; zero positives have upper
bound .1128021742. All zero-positive conditions abstain, so they do not demonstrate
test calibration. The shared-state group recovery is recovery of a simulated
allele-expression state, not genomic membership. Greater concentration-induced
shrinkage reduces matching without removing the repeated reserved contrast.
Neither setting is selected or promoted.

Under the two observationally identical biological laws, any RNA/allele-only
rule has equal rejection probability. Consequently, if its false genomic-split
probability is bounded by 5% for this shared-state law, its power against the
identical genomic alternative is at most 5%. This conditional identification
limit follows from equal probability laws; it is not a universal bound for other
genomic effects, independent assays or independently justified ASE restrictions.
More count splitting, narrower nominal p-values or extra null draws cannot
separate these two labels. A correction must state what independent information
restricts the nuisance, or abstain from the genomic interpretation.

## Independent evidence and decision

Binomial job 23003288 completes in 32s with 677340 KiB process peak RSS on
innovation CPU. Concentration-2 job 23003289 completes in 40s with 699216 KiB on
nslab CPU. Audit 23003362 completes in 86s. All exit zero. It verifies source and
output hashes, all paired input matrices, independent role hashes, frozen fitted
arrays/metadata/logs, eligible cells/regions, every permutation and exhaustive
assignment matching. Scalar integer count sums and scalar permutation means
reproduce all statistics with zero discrepancy. A floating-operation allowance
is derived from bounded-fraction summation; no available decision intersects
that rounding envelope. Independent 80-digit Decimal roots verify the extreme
exact bounds. Actual resources and immutable outputs remain under
`dev/allele_shared_state_20261108/{binomial_v1,beta_v1,audit_v1}`.

Withhold genomic confirmation from this working rule. Retain the computational
role splitter as a tested research implementation, with raw reserved summaries
separate from any claimed DNA interpretation. No new production dependency,
default, calibration flag, scientific cutoff, original test or frozen output
changes. The unchanged wheel's full suite was not rerun. Original thinning remains
2/5 with three unresolved failures; the broader joint error/power gate remains
failed and certified biological null n remains zero.

Next require an independently observed constraint on shared nuisance or suitable
DNA truth. Any future joint-channel procedure must reproduce positive evidence,
nuisance estimation, selection and stopping under both null and alternative,
and retain power/abstention alongside the same 5% specimen error criterion.
This negative result satisfies a scientific falsifier, not the overall package
accuracy or clone-reliability goal.
