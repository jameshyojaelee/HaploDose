# The three production failures do not share one profile-contraction mechanism

[Observed pattern, not a guarantee] The unchanged five-case production replay
shows different failure locations. Seed20261108 carries a nonzero chr12 contrast
through shared-profile updates before merging removes the distinction.
Seed20261109 never enters an observed S1 update with two tumour components.
Seed20261110 enters its S1 two-tumour update with an already tiny chr12 contrast.
Neither a blanket gene-veto change nor disabling profile sharing follows from
these observations. Original production recovery remains2/5.

All five public tables, posterior arrays and30 criterion values exactly reproduce
the prior corrected production replay. The new trace records67 profile updates
across15 specimen models, with ordinary production thinning and no true-noise
oracle. All54 gene checks in the preceding diagnostic already passed.

## Estimand and identifying assumptions

The diagnostic objective is sum_ck R_ck log NB(Y_c | mu_ck,phi), holding R,
reference, offsets, dispersion, depth effects and programmes fixed. Clone0 remains
zero-dose as in the actual controller. Compare incoming doses, individual pooled
profiles and shared profiles. This is an unpenalized conditional working-NB
objective in natural-log probability units. It is not a genomic likelihood,
independent confirmation or the package's complete profile-prior objective.

Individual and shared profiles also use a PLN working model, segmentation and
regularization. Those are competing explanations for a drop in unpenalized NB
fit; a decrease alone does not establish a computational bug or justify removing
a penalty. Correlated state and inadequate reference modelling remain competing
biological explanations for inferred genomic distinctions. Corresponding DNA,
rather than a favorable objective value, would distinguish those biologically.

The chr12 diagnostic is the maximum absolute pairwise difference in the actual
RNA log-rate vectors over exported bins with chromosome12 and start<10Mb.
It is a descriptive planted-region probe, not a whole-region effect estimate,
DNA-copy-number estimate, breakpoint accuracy measure or inferential threshold.

## Independent numerical evidence

The audit reconstructs count and exposure pools with dense matrix products,
checks depth moments and exact binned-matrix row correspondence, and compares
SciPy NB log PMFs with stable NB likelihood ratios accumulated in long double.
Installed float32 kernel discrepancies are bounded with operation/reduction
lengths and absolute cancellation magnitudes. Values within this numerical
guard remain unresolved; they do not receive a convenient sign or cutoff.
These guards concern arithmetic, not statistical coverage.

| Update comparison | Guarded increase | Guarded decrease | Unresolved |
|---|---:|---:|---:|
| Individual versus incoming |2|20|45|
| Shared versus incoming |3|1|63|
| Shared versus individual |22|0|45|

All67 direct/reference comparisons pass. A universal claim that shared profiles
destroy the working likelihood is unsupported. These descriptive counts span
correlated updates from five simulations, not67 independent biological units.

For S1 seed20261108, the first two-tumour update's chr12 maximum changes
0.201742→0.224162; it reaches0.291774 before a later merge and remains0.265159
near the end of the second multi-component refinement. Thus the contrast is not
simply zeroed by sharing. Its later depth-matched merge and membership population
remain the relevant branch. Seed20261109 has no observed S1 two-tumour update,
so proposal/confirmation precedes any such profile failure. Seed20261110's only
such update changes0.00104960→0 for chr12, while other differences remain; the
incoming proposal already lacks the intended chr12 signal.

These are mechanism localizations, not proof of the remaining biological cause.
No profile-optimizer replacement or relaxed merge threshold is promoted. Any
future correction must retain the original recovery/error/power criteria and
separately address the count-splitting law, selected memberships and nuisance
identification rather than assuming all three failures have one bottleneck.

## Reproduction and remaining limits

Immutable design/source, all input/parameter/output arrays, installed-kernel
scores, independent audit and environment inventories are under
`dev/production_profile_trace_20261109`. Authoritative records are
`trace_v1/report.json` and `audit_v1/report.json`. Jobs23122416 and23122417
complete0:0 on CPU in1m33s and59s. This study changes no production source,
scientific threshold or frozen artifact. The validated wheel remains the
403-Python/4-Rust-test artifact described in the exact-merge report.

Biological n=0. The three production recovery failures, state-null error criterion,
membership calibration and interval coverage remain unresolved. The research
5/5 path still supplies true gene-noise information. A deployable estimate needs
its own observation-model and replication checks. More computation cannot replace
missing corresponding DNA and biological replication.
