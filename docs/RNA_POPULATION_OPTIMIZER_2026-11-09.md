# Label-invariant numerical audit preserves the failed clone tradeoff

Correct the comparator's refinement diagnostic:67 cases had changed literal
component labels, but4 were complete reversals of two arbitrary names. There
are63 changed partitions among998 cases, with up to289 changed cells. None of
the32 accepted evaluation distinctions changes partition under refinement.
Independent optimization finds no material supported gain over the saved fine
fits in the corrected eight-case probe. No full comparator replay or HaploDose
method correction is justified by this investigation.

## What changed and what remains fixed

`benchmarks/hardening/rna_population_baseline.py` now compares binary partitions
over both global label permutations, retaining literal label differences as a
separate field. It also reports unordered co-membership disagreements. For n
cells and d literal differences, aligned cell disagreement is min(d,n-d) and
pairwise disagreement is d(n-d). These identities compare repeated fits, not
predictions against true labels; no accuracy is gained by using simulation truth.

Audit23066017 checks all1024 pairs of five-cell binary fixtures by exhaustive
label mappings and direct pairwise comparisons. All998 saved records are
recounted. Eight reruns of the corrected helper require exact counts, loadings,
probabilities, coarse/fine scores and original literal disagreements. All pass.
Only diagnostic metadata changes. Production code, original tests, saved fits,
scientific thresholds, confidence families and error/power events are unchanged.

The older67-case number in the previous report is a literal-label diagnostic,
not the number of changed partitions. Its frozen report remains unchanged;
the new audit supplies the corrected interpretation and explicit raw counts.

| Role | Saved cases | Literal-label changes | Changed partitions | Pure reversals |
|---|---:|---:|---:|---:|
| NB calibration |199|19|17|2|
| State calibration |199|11|11|0|
| NB evaluation |200|17|16|1|
| State evaluation |200|11|10|1|
| Two-population evaluation |200|9|9|0|

## Independent optimization and membership uncertainty

The estimand is a fitted RNA population under the same transformed-RNA Gaussian
working model and known reference as the preceding comparator. It is not a
genomic clone probability. Regional transcription remains a competing biological
explanation; same-cell DNA or independently justified state constraints would
be needed to distinguish it. These simulations add no biological replication.

The first numerical probe selects the four greatest score differences and four
greatest literal MAP differences. It therefore includes four pure label
reversals. Keep this probe; it identifies higher local likelihoods in some
near-collapsed fits but does not certify global optima. The largest recorded
score difference6.59098 is an IMPROVEMENT of the fine solution, not evidence
that refinement discarded a better coarse solution.

The corrected probe replaces literal disagreement with label-invariant
disagreement for selecting the four most unstable partitions. The score-based
selection and lexical tie rules remain unchanged. Both probes independently
reconstruct saved likelihoods and fine MAP assignments from exact frozen counts
and saved PCA loadings, then run L-BFGS-B continuation and independent SciPy
normal-density Nelder-Mead from both saved starts. Each probe has8 cases and32
solver runs. No truth, detection or error outcome chooses cases or solutions.

The corrected probe converges in all32 runs. No gain beyond the operation-scaled
roundoff allowance also meets the existing30-cell MAP support rule. In the two
large supported partition changes (151 and126 cells), fitted scores are only
0.00949 and0.04109, and changed-cell median maximum probabilities are approximately
0.554. Those are weakly separated RNA mixtures, not well-identified genomic
memberships. In the largest change (289 cells), the fine fit's smallest component
has only1 cell; other large score changes similarly favour1–3-cell components
that fail the support rule. A high conditional component probability alone cannot
establish a distinct population or validate DNA membership.

This is a bounded numerical falsification: it provides no supported optimizer
intervention for these cases. It neither proves global optimality nor eliminates
all possible numerical problems. The earlier original-NB proposal investigation
is a separate model and experiment; see
[its independent solver evidence](JOINT_CONFIRMATION_2026-11-08.md#independent-proposal-optimizer-follow-up).
Do not transfer a Gaussian-comparator finding into a HaploDose correction.

## Consequences for the scientific result

Audit23066716 joins all600 fixed evaluation decisions to unique role/replicate
diagnostic records, with no missing or duplicated specimens. All original event
counts agree exactly. No new inference, threshold choice, p-value or fit occurs.

| Evaluation role | Accepted distinctions | Accepted with changed partition |
|---|---:|---:|
| Matched NB null |8/200|0|
| State null |10/200|0|
| Two-population alternative |14/200|0|

The simultaneous upper bounds0.110230 and0.124699 still fail the5% criterion;
power remains14/200. These are conditional simulated results, not biological
guarantees. Numerical stability of specimen decisions does not validate the
model or genomic interpretation. Original thinning recovery2/5, its three
remaining failures, uncertainty failures and missing biological null panel remain.
Do not promote the comparator or claim SOTA. Do not add draws or relax criteria.

## Reproducibility

Immutable designs, source snapshots, selected identifiers, solver parameters,
scores, gradients, posterior-change diagnostics, counts/array identities,
execution environment and resource records are under
`dev/rna_population_optimizer_20261109`. Final consequence evidence is
`summary_v1/report.json`; label proof is `label_v1/report.json`; corrected
numerical probe is `probe_v3/report.json`. Scripts are
`rna_population_optimizer_probe.py`, `rna_population_label_audit.py` and
`rna_population_optimizer_summary.py` under `benchmarks/hardening`.

All jobs run on Intelv4 compute CPUs using innovation/nslab. Failed probe23065331
is retained: its independent reference initially compared an unclipped likelihood
gain to a statistic explicitly defined as max(0,gain). Probe23065619 applies the
same clipping with the unchanged numerical error budget. Successful label and
solver audits preserve all prior fits and decisions. The package suite was not
rerun because production and original tests are unchanged; the corrected helper
and independent audit execute on compute nodes. No dependencies or shared
installations change. Numerical and biological uncertainty remain separate.
