# Allele proposals and balance selection honor the chosen cell model

An explicitly supplied cell-level beta-binomial concentration now applies to
allele split proposals and balanced-cell candidate selection, as it already did
to final cell memberships. Previously both selection paths silently used binomial
contrasts. This is a computational correction to an existing option, not a new
scientific default, calibrated clone method or independent biological validation.

## Contract, counterexample and competing explanation

The cell-level observation model is A_cs given D_cs distributed beta-binomial
with mean haplotype fraction b_ks and concentration cell_kappa. With None the legacy
binomial working model remains. The selected estimand is a conditional cell-level
log-likelihood contrast between the proposed fraction and balance, using the same
shrinkage and fixed decision thresholds as before. Pooled-clone concentration is
a different parameter; this correction does not estimate or substitute it.

The old-wheel counterexample fixes pooled fractions/run discovery to isolate the
cell likelihood. With cell_kappa=2 and fractions shrunk from 0.8 to 0.77, three
cells have independent SciPy beta-binomial log contrasts 0.19855, -1.89863 and
-1.04487. Expected proposal labels are [1,0,0]. The old code returned binomial
labels [0,1,0]. It also selected a cell with 500 of 1,000 reads on the haplotype
as balanced at the unchanged 0.99 threshold, even though its selected-model
working balance probability is 0.63085294. The corrected code matches both
independent expectations.

These are fixed numerical fixtures, not biological specimens or a power study.
Conditioning on the proposed fractions and giving both hypotheses equal prior
weight does not account for their estimation, run selection or other adaptive
search. Strong context-specific allele expression, mapping bias and phase error
remain alternatives to genomic imbalance. Independently reserved allele evidence
with an appropriate model, and corresponding DNA at the claimed resolution, are
needed to distinguish them. A finite concentration does not itself supply those
measurements or certify the resulting probabilities.

## Small correction and falsifiers

Both paths now call a shared balance contrast from the same observation-model
implementation as cell membership. Selection contrasts retain float64; the public
membership log-likelihood output retains float32. The None branch keeps the exact
original arithmetic order. No threshold, shrinkage, prior, SNP allocation,
minimum-cell rule, HMM parameter or default changes.

Ten new checks cover independent SciPy PMFs, zero coverage, haplotype reversal,
mixed-depth proposal labels, balance selection, invalid concentrations and legacy
precision. PMF cancellation tolerances scale with the count-dependent log-gamma
terms and floating arithmetic, rather than a universal scientific epsilon. The
decisive counterexamples use the original sign/0.99 criteria with no tolerance.
Different independent labels, changed legacy fixture values or any altered
preexisting test falsify the correction's stated computational contract.

The allele-plus-expression-null configuration remains rejected. This repair does
not add a matched joint-channel null or establish specimen false-split control.
The original three thinning failures and joint false-split/power failures remain;
their frozen results are not rerun or overwritten. Biological null n remains 0.

## Compute validation and preservation

Old-wheel probe 23002432 completes in 8s and reproduces both defects. Environment
setup 23002452 completes in 3m29s on the I/O partition, copying the exact existing
Conda/pip dependencies to an isolated prefix. No new production dependency or
shared installation is introduced. Build/validation 23002524 completes in 2m36s
on a CPU node: 348 Python tests and 4 Rust tests pass; 1 GPU-only test is skipped.
Python-test process peak RSS is 2525528 KiB. GPU execution was not run for this
CPU/SciPy allele calculation.

Corrected probe 23002539 completes in 5s. Independent audit 23002671 completes in
4s, verifies both corrections and exact old/new None fixture equality, hashes all
494 baseline source files, and permits only alleles.py, clones.py and the new test file
among package/native/test changes. Every preexisting test is unchanged. The first
audit 23002576 failed an incorrectly anticipated 346-pass report string; its output
is retained. The revised audit reads the actual 348-pass result without changing
the package, tests, tolerances or any scientific criterion. SeeAUDIT_NOTE.md.

WheelSHA256:c7f4bc04a08d1bef8ed5ebd37969b1f6cb6e247271da1d3cf42c5a16d308d71a.
NativeSHA256:e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7.
Sources, environments, old/corrected outputs, failed/successful audits and actual
backend/resource records are under dev/allele_cell_model_20261108. The authoritative
candidate is candidate_v2; installed artifact is validation_23002524; final audit
is audit_v2/report.json. No claim of competitor superiority follows from this fix.
