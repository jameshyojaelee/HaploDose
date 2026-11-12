# Coherent gene halves do not prevent membership loss

With the fixed known-gene-noise research split,20261109 accepts a chr12
distinction twice and merges it away twice.20261110 loses its proposed chr12
distinction during initial confirmation agglomeration, before matched testing.
These are different failure stages.20261108, now recovered, retains its accepted
distinction through merge. Production recovery remains2/5, research3/5; no default
or scientific criterion changes.

## Locked observational trace

`dev/gene_thinning_stages_20261109/DESIGN.md` fixes the three original failure
cases before execution. A Python call/return profiler copies existing values
from every S1 event fit, candidate score, agglomeration, matched comparison,
support and merge call. It captures all event regions, not a selected successful
window. No operation, seed, proposal, threshold or profile is replaced. The
gene-native prototype still uses unavailable simulator size10. Truth is withheld
until inference finishes and used only for descriptive audit.

Independent audit confirms exact replay of fitted cell/difference tables and all
18 original criterion values, including failures. All568 recorded stages pass
identity/hash checks. Original candidate, support and merge intermediates are
saved, including input profiles before merge mutation. The audit verifies actual
proposal-half argmax labels and compares their postfit truth contingency tables
with the full-count assignments. These are working-model assignments, not
calibrated genomic membership probabilities.

| Seed | Candidate scores | Accepted candidates | Merge calls | Calls reducing group count | Original research recovery |
|---|---:|---:|---:|---:|---|
|20261108|14|2|2|0|pass|
|20261109|12|2|5|2|fail|
|20261110|6|0|1|0|fail|

For20261109, both accepted candidates cover bins209–217. Their matched scores
are5.917695 and5.386993. At the first merge the proposal assignments place
183 clone1 cells in group1 and155 clone2 cells in group2, versus202 and177
using full-count assignments. Penalizing reference-assigned tumour cells gives
338/510 versus379/510 correctly grouped. This is a descriptive common-denominator
comparison, not the final recovery metric conditional on tumour calls.134/600
cells change assignments. The actual matched comparison has no resolved segment;
the second merge repeats this loss with132/600 changed assignments. Neither
comparison establishes genomic equivalence.

For20261110, the implanted-region event fit has gain12.985298 and proposed group
sizes361/149. Its postfit tumour-label table is[[267,34],[87,114]], membership
accuracy381/502 on tumour cells within that proposal, distinct from the510-cell
denominator and the final criterion. Confirmation pools both proposed parts into
one group; no matched comparison is reached for this candidate. Replacing only
terminal merge cannot by itself restore that lost proposal.

For20261108 the first merge changes113/600 assignments, yet retains bins209–217
with minimum retained-bin posterior0.997600 under the existing HMM. Assignment
change alone is not a sufficient failure explanation. Full-count assignments
reuse test observations and are not an independent remedy; earlier joint-data
experiments' state-null failures remain relevant negative evidence.

The estimand is distinguishable regional RNA-profile populations under fixed
research folds. The NB observation model, fitted nuisance and adaptive search
remain working assumptions. A continuous expression state can produce selected
groups. The observation needed to establish a genomic distinction remains
corresponding DNA at suitable resolution; these simulated truth descriptions
cannot supply biological validation or a specimen false-split guarantee.

## Execution and decision

Three CPU traces23119344/345/346 and independent audit23119611 complete0:0.
Trace requests4CPU/8GB/2h, audit1CPU/8GB/1h, innovation/nslab, no arrays.
Elapsed scheduler times41/40/39/4 seconds. Full accounting, process RSS/runtime,
installed wheel and CPU backend identities are preserved in the task records.
All76 protected files match the initial dirty package/native/test/README/build
baseline. Prior frozen results and all scientific criteria remain unchanged.
Full package tests are not rerun for observational research code; the installed
wheel's existing395 Python/4 Rust passes and1GPU-only skip remain separate evidence.

This trace justifies a bounded interaction test of the existing membership-weighted
confirmation/merge procedures under exactly these same gene splits. It does not
justify accepting a split without validation or weakening posterior, cell, gene
or whole-specimen error/power criteria. Biological n=0. No coverage, fresh
false-split/power or independent DNA validation has been established here.
