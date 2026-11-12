# Reserved SNP counts change the allele proposals

Changing only odd-index SNP counts changes the current allele proposal's even-SNP
posteriors, discovered runs and cell partition. Independent auditing accounts for
label exchange and confirms different partitions despite unchanged group sizes.
A reference that physically subsets proposal SNPs before fitting is exactly
invariant. This falsifies local confirmation separation; it does not establish a
complete replacement workflow or a biological false-split rate.

## Observation, estimand and isolation

The estimand is computational dependence of a proposal on its designated
confirmation observations, at fixed cell memberships and reference inputs. The
working model is the package's phased beta-binomial HMM and cell likelihood.
There are 120 deterministic fixture cells, not biological specimens. Forty SNPs
share fixed positions and genetic distances. Every SNP depth is 20. Even-index
alternate counts are 10, 12 or 16 in three fixed forty-cell groups. Odd-index
counts are set to 10, 19 or 1 in three prespecified scenarios.

The proposal count arrays have identical hashes across scenarios; only
confirmation alternate counts change. Initial memberships, reference concentration
50, all default scientific thresholds and the cell observation model are fixed.
The current path calls the actual `set_allele_pools` and `allele_splits` methods.
The first fits pooled HMM states and switch/concentration parameters using all
SNPs. The second scores only even SNPs, but inherits that all-SNP fit.

An alternative explanation for a changed output is a mere label exchange or
changed interval name. The independent audit canonicalizes each partition as
unordered cell sets, eliminating both explanations. It verifies literal fixture
counts and array hashes independently of the runner.

## Results and falsification

| Quantity | Odd balanced | Odd same orientation | Odd reverse orientation |
|---|---:|---:|---:|
| Current proposed group sizes |40/80|40/80|40/80|
| Current partition equals baseline, up to labels |yes|no|no|
| Current switch parameter |0.01|0.01|30|
| Current concentration |50|50|50|
| Proposal-only reference partition equals baseline |yes|yes|yes|
| Proposal-only reference switch parameter |0.01|0.01|0.01|
| Proposal-only reference concentration |250|250|250|

Current even-SNP posterior entries differ by up to 0.9999998 and 0.9991521 from
the balanced baseline. The physically subset reference has exactly equal fitted
parameters, posterior bytes, runs and labels in both contrasts. It uses twenty
observed proposal SNPs with the original genomic distances and support thresholds;
it does not count masked zero-depth confirmation SNPs toward minimum support.

The prespecified separation criterion was exact equality under confirmation-only
perturbations. The current path fails it; the subset mechanism reference passes
at fixed memberships. This is a software dependency counterexample, not a null
specimen experiment, uncertainty calibration or estimate of prevalence. No
scientific threshold, confidence family or original recovery criterion changes.

## Correction boundary and biological interpretation

A local proposal-only HMM is an identified remedy for this particular dependency.
It is not yet promoted. Upstream membership estimation, initial clone retention,
run selection and repeated merging can also consume full allele data. A proposal
fit based on those memberships may still depend on confirmation observations.
Similarly, reusing concentration or switch parameters fitted on all SNPs defeats
separation even after masking emissions. These dependencies require a complete
discovery/confirmation contract before any independent-confirmation claim.

The biological alternative remains context-specific allele expression, mapping
bias or phase error rather than DNA alteration. Independent corresponding DNA
and suitable held-out measurement can distinguish those explanations; equal
proposal bytes cannot. Even an independent conditional SNP test would not alone
establish whole-specimen false-split control, supported cell membership or CNA
causation. The mismatched allele-plus-expression-null configuration remains
rejected. Original three thinning failures and false-split/power failures remain.

## Compute and reproducibility

Diagnostic job 23002808 completes in 13s on a CPU node with 529400 KiB process
peak RSS. Audit 23002857 completes in 2s with 79740 KiB peak RSS. Both use the
isolated wheel validated by 348 Python and four Rust tests, with one GPU-only
skip. No production code changed and the full suite was not rerun. No frozen
result, package environment, existing test or dataset role changed.

Implementation is `benchmarks/hardening/allele_confirmation_leakage.py` and its
independent audit. Design, source snapshots, input matrices, fitted posterior
arrays, exact labels, hashes and execution records are under
`dev/allele_confirmation_leakage_20261108/{probe_v1,audit_v1}`. The decision-complete
dependency gate is `dev/allele_confirmation_leakage_20261108/NEXT_GATE.md`.
