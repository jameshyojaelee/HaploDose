# Confirmation cannot resolve a score tie by candidate order

The frozen confirmation scorer selected the first maximum even when several
candidates had exactly equal scores. Such an assignment can change identity
when the candidate columns are permuted. The correction abstains on exact ties
and carries the reason through chunk aggregation into the cell table. It does
not repair the three thinning recovery failures or establish clone reliability.

## Scientific and computational contract

The immediate estimand is a frozen RNA candidate label, conditional on the
discovery profiles, priors, dispersions and nuisance parameters. Proposal features
provide a composite NB score; test features supply separate nominal contrast
diagnostics. The strongest alternative to an apparent membership is an equally
supported candidate, including a candidate with the same proposal profile but a
different test profile. Candidate array order is not distinguishing evidence.

Independent measurements that distinguish the candidate profiles could resolve
this ambiguity. Reusing confirmation test counts to choose membership would
instead change the role contract. Even a unique proposal maximum does not prove
DNA membership, identify continuous state versus CNA, or calibrate uncertainty.
Near ties and assignments driven only by unequal priors remain unvalidated.

The correction uses exact equality, with no new numerical or scientific cutoff.
It excludes tied cells from either candidate's contrast pool. A candidate-specific
programme score is also marked unavailable for a tied cell. Existing score
arithmetic, priors, discovery fitting, merge rules and scientific criteria remain.

Cell exports add `confirmation_unresolved_reason`. Values distinguish
`tied_candidate_scores`, `no_proposal_depth`, `no_finite_candidate_score` and
`not_fitted`; assigned cells and non-confirmation rows have an empty reason.
The existing `confirmation_status` categories remain available to readers.
Historical result files are unchanged.

## Counterexample and decisive tests

The preceding installed wheel was run on three deterministic grids of 100 cells.
There are 100 tied cells with identical proposal profiles and different test
profiles, ten ties from symmetric observations under two different profiles,
and no ties when identical profiles have unequal priors. The old wheel assigns
the first candidate in all 110 ties. Independent SciPy log-PMF discrepancies are
at most 1.25e-14, reported as numerical diagnostics, not a scientific tolerance.

The final audit uses exact rational NB probabilities to determine all 300
expected assignments, including abstention, without a floating-point epsilon.
It requires byte-value equality of old/new scores and exact preservation of all
190 unique-score assignments and their separately computed contrast diagnostics.
Perturbing test counts must leave assignments and proposal scores identical.
Unit checks additionally cover candidate permutations, a single candidate,
contrast pooling and reason preservation across a 4096-cell chunk boundary.
Existing discovery-leakage and recovery assertions are retained.

Any changed score, changed unique-score assignment, failure of the exact rational
ordering, use of test evidence for membership, or lost reason falsifies the
proposed narrow correction. The audit also verifies the pre-edit snapshot and
unchanged installed modules outside confirmation and pipeline exports.

## Execution

The pre-edit source is preserved under `dev/confirmation_ties_20261108/baseline`.
Environment staging22998195 and old-wheel probe22998221 run on compute nodes.
Candidate source `candidate_v1` is built and tested by22998241. The independent
audit uses `audit_source_v2` in22998307; the earlier pending audit22998245 was
cancelled before execution to add exact rational reference calculations.

Validation22998241 passes 336 Python tests and four Rust tests, with one GPU-only
skip on the CPU allocation. It takes 6m54s including build/install; test process
peak RSS is 2531632 KiB. Audit22998307 passes in 11s, with peak RSS 119104 KiB.
All 300 rational assignment checks pass: 110 tied cells abstain, while all 190
unique assignments and their contrast diagnostics are unchanged. All scores are
exactly equal to the preceding wheel, and 457 baseline files remain unchanged.

The new noneditable wheel has SHA256
`88e26f53b5afe651bba09021cc0f5868248879102708e04c0df979fa73b9f2e9`.
Its isolated environment is `dev/confirmation_ties_20261108/environment`;
the execution records link the wheel, source, native binary, inputs and outputs.
Final audit evidence is `dev/confirmation_ties_20261108/audit_v2/report.json`.
No task jobs remain active. Frozen results and the prior false-split/power
rejection remain unchanged; biological n is zero. Retain this computational
correction without promoting a biological or uncertainty-calibration claim.
