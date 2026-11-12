# Binning loses some membership information; the matched NB approximation loses little here

On the three original failed recovery inputs, the known-generator summed-count
oracle correctly assigns439/510,407/510 and403/510 true tumour cells. The matching
gene-count oracle assigns441/510,411/510 and420/510. Moment-matched bin NB changes
only1,1 and0 decisions relative to the actual summed-count distribution. Binning
therefore loses some information, but the retained bin information remains
substantial when parameters and classes are known. These experiments do not
support aggregation alone as a sufficient explanation for the current failures.
They do not establish that estimated pipeline parameters are adequate.

## Contract and controlled comparison

Task `dev/gene_bin_information_20261109` locks `DESIGN.md` before new computation.
The cases are exactly seeds20261108/109/110 from the original frozen recovery
panel. Selected genes, their source-column order and bin memberships come from
the preexisting `recovery_next_20261107/trace_SEED/fit/genes.tsv.gz`. There is no
new feature selection, binning, smoothing, threshold or omitted method.

The estimand is realized classification and true-class log predictive score
under three known-parameter representations of the same raw measurements.
Observation model: independent gene NB counts from the unchanged simulator.
True class-specific gene rates, nominal per-cell exposures, depth slopes,
dispersion and population priors are deliberately supplied only to this
post-fit diagnostic. They are not available to HaploDose. Expected gene counts
retain the original full-gene composition normalization before selecting genes;
exposure is not replaced by the observed retained-gene library total.

All three representations use the same retained genes, cells and parameters:

1. Full independent gene NB likelihood.
2. Actual independent-gene NB sum PMF for each observed bin count.
3. NB likelihood matching each bin's true mean and variance, with
   variance=sum(mu_g+phi*mu_g^2). This supplies the true conditional variance;
   it is not the pipeline's estimated reference dispersion.

Comparisons1-versus2 isolate representation loss under this known model;
2-versus3 isolate the bin-distribution approximation under identical known
parameters. A previously reported67% state-mixture bin oracle uses a different
generator and is not a valid aggregation comparator for these original cases.

Biological n=0. These are three selected development mechanisms, not biological
replication, prospective evaluation or a population superiority test. Known
classes make this a classification diagnostic; it does not discover clone
number. A realized oracle accuracy is neither a per-realization upper bound
nor a deployable performance guarantee.

The strongest competing biological explanation remains regional transcriptional
state rather than DNA. This generator labels its rates with CNA events, but
RNA likelihood agreement cannot independently establish that interpretation.
Independent DNA or justified nuisance restrictions remain necessary in biology.

## Numerical method and falsification

For independent NB genes with sizes r_g and probabilities p_g, let q_g=1-p_g.
The positive PGF coefficient recurrence is

    h_0 = 1
    h_n = sum_{k=1}^n [sum_g r_g*q_g^k] * h_{n-k} / n
    log P(sum_g Y_g=y) = sum_g r_g*log(p_g) + log(h_y)

Truncation at observed y loses no contribution to that coefficient: a
nonnegative gene count above y cannot contribute to total y. Positive long-double
arithmetic avoids FFT cancellation. Operation/magnitude-based allowances retain
both log-zero and log-coefficient magnitudes; they are numerical accumulation
budgets, not confidence intervals. A decision margin overlapping its budget
stops interpretation rather than receiving a relaxed cutoff.

Numerical gate23110287 passes41 independent checks: exact Fraction coefficients,
equal-probability NB closure, direct multiplication of SciPy NB probability
polynomials, and row/gene permutation invariants. A two-gene geometric example
has identical count sums but likelihood ratios160/243 and1000/2187 for the
two allocations. Thus sums are not generally sufficient statistics for
heterogeneous gene NB models. This algebraic example is not an identification
claim about whole-genome doubling under an RNA composition model.

Independent audit23110869 uses raw frozen counts rather than regenerated input,
direct positive probability-polynomial convolution rather than the primary PGF
recurrence, an explicit gene/moment NB formula, and scalar log-score calculations.
It checks25 boundary/class-extreme cells across all retained bins/classes:
31,374 independent bin/class PMFs. Maximum summed log-likelihood discrepancy
is1.82e-12, within the recorded budgets. All1800 cell identities/labels, mappings,
bin sums, class priors, score summaries, paired decision tables and numerical
decision margins are independently checked. Primary gene assignments reproduce
all prior retained-gene oracle decisions exactly. No numerical check fails.

## Matched outcomes

The denominator is all510 true tumour cells, including cells assigned normal;
this differs from conditioning on a fitted malignant call. No new recovery
criterion is substituted for the original full workflow.

| Seed | Gene NB correct | Actual bin sum correct | Moment-NB correct | Gene/bin differing decisions | Actual/moment differing decisions |
|---|---:|---:|---:|---:|---:|
| 20261108 |441/510|439/510|440/510|14/600|1/600|
| 20261109 |411/510|407/510|406/510|24/600|1/600|
| 20261110 |420/510|403/510|403/510|30/600|0/600|

For the first two seeds, every representation has0/90 normal false calls and
detects510/510 tumour cells. For20261110, gene NB has0/90 false calls and510/510
detections; both bin models have3/90 false calls and505/510 detections.

Mean true-class log scores(gene/actual-bin/moment-NB) are
-0.283696/-0.294887/-0.294885,
-0.371763/-0.372712/-0.373146 and
-0.331836/-0.364591/-0.364965 respectively. These are descriptive scores on the
same generated specimens. No superiority p-value or population interval is
attached to them. Finite-realization classification need not obey a universal
ordering; for example, moment matching changes one decision beneficially in
the first case and adversely in the second.

The observations reject the claim that aggregation necessarily destroys the
information needed to distinguish these fixed populations. They do not prove
that a particular estimator can learn the supplied rates/nuisance, recover
unknown groups or control false distinctions. Parameter estimation, proposal
formation, confirmation and merging remain implicated alternatives, not a
uniquely identified cause. A gene-level replacement is not promoted merely
because it retains more information in principle.

Marginal classification agreement does not establish the joint law of proposal
and confirmation counts. An unequal-gene NB sum can depart from NB even when
moment-matched class decisions agree. Thinning that sum as though it were NB
need not produce independent halves; this experiment does not validate the
frozen thinning size, confirmation calibration or complete adaptive null. That
assumption remains a distinct unresolved mechanism, not something excluded by
the small marginal decision differences.

## Implementation, package correction and reproducibility

Research scripts are `gene_bin_reference.py`, `gene_bin_information.py` and
`gene_bin_information_audit.py`. The three case jobs reproduce frozen CSR
values/indices/row pointers plus cell/gene identities exactly, then save selected
identifiers, masks, bins, true parameters, scores, budgets and predictions in
`seed_SEED_v1/evidence.npz`. Source and inputs are hashed in execution/run records.

The sole production-file edit corrects the module docstring in
`python/haplodose/bins.py`: bins are count summaries for working models, not
generally sufficient statistics for heterogeneous gene-level models.
`documentation_ast_parity.json` verifies all functional AST is unchanged.
`measurement_source_integrity.json` establishes that all76 protected files were
unchanged during the measurement runs. A new isolated environment and immutable
snapshot are used to build/test the corrected documentation wheel; prior
environments, wheels, scientific criteria and frozen results are preserved.

All seven CPU jobs complete successfully: numerical gate23110287, measurement/
audit23110866–23110869, isolated environment23111357 and wheel validation23111358.
The new noneditable wheel passes395 Python tests and4 Rust tests;1 GPU-only test
is skipped on CPU. The unchanged original tests remain byte-identical. Wheel
SHA256 is `e30775c72f93251d5105dc0924e9bd0b6797570c98790db418a8d6aae6a5a72d`;
native SHA256 remains
`e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.
Artifact/environment inventories, tests and exact wheel path are recorded in
`validation_23111358/installed_artifact.json` and adjacent reports. The new
environment is `dev/gene_bin_information_20261109/environment`; the old pinned
population environment remains unchanged.

Numerical gate uses CPU1/4GB/1h; cases CPU1/8GB/2h; audit and environment setup
CPU1/8GB/1h; build/tests CPU8/32GB/2h. Accounts are innovation/nslab, without
arrays. Measurement jobs peak below406000 KiB, audit309676 KiB. Python test
execution is117.86s and its measured process peak2529304 KiB. Setup-job actual
peak memory was not recorded; Slurm zero/missing sampled RSS is not interpreted
as a memory measurement. All headers, commands, runtimes, failure status,
available resource records and dependencies are in `jobs.json`,
`accounting_final.tsv`, `analysis_JOB/`, `environment_setup/` and
`validation_23111358/`. AST, JSON and new link/evidence-path checks pass.
`source_integrity.json` verifies the module docstring is the sole changed
protected file; the other75 are unchanged. Full source, parameters, counts,
scores and audit remain under `dev/gene_bin_information_20261109/`.

Computational status: the known-model comparison and independent audits pass.
Identification: these oracle inputs are unavailable to the method; DNA identity
is not inferred from a correct count model. Biological validation: absent.
Original production recovery remains2/5. No new clone recovery, calibrated
uncertainty, specimen false-split/power or package superiority claim is made.
