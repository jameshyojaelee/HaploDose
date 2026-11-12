# Processor-sensitive clone proposals traced to PCA precision

The objective is numerical correctness of the same RNA observation model across
processors, without changing scientific acceptance rules. A bounded investigation
identifies and corrects the recorded neighbour-graph discrepancy in two cases.
It does not establish biological clone validity, universal processor equivalence
or resolution of the remaining recovery failures.

## Primary execution chain and cause

The locked cases are the original two-clone simulation replicates150 and167,
seed keys[20261107,20,2,replicate], with600 cells,240 bins,60-bin genomic resets
and even/odd discovery. The exact same saved counts, references, dispersion and
other arrays are replayed with the validated wheel on Intel v4 and AMD v6.
Both processors reproduce their corresponding historical final responsibilities
exactly. Uninstrumented and instrumented outputs are also exactly equal on each
processor. Thus neither regenerated observations nor tracing explains the change.

The first differing stage is randomized PCA inside initial Leiden proposals.
Input feature matrices are identical. Float32 PCA coordinates differ by up to
3.04e−5 and4.17e−5. Each symmetrized unweighted graph differs by two undirected
edges. Independent double-precision distance calculations agree with the graph
constructed from each processor's own coordinates. Reordering graph edges does
not alter membership. These observations locate the mechanism in PCA numerical
precision rather than graph insertion order or a different distance rule.

One neighbour-cutoff gap in each case lies within a float32 distance-accumulation
scale. This scale is a diagnostic for arithmetic in the distance, not a total
error bound for randomized PCA or a scientific tolerance. Float64 PCA gives
coordinate discrepancies of1.35e−13 and5.54e−14 and identical neighbour graphs
in these probes. All algorithmic dimensions, seeds and clustering parameters stay
fixed. The contrast does not justify rounding away future boundary disagreements.

The old discrepancy propagates into initialization and merging. Replicate150
has maximum final responsibility difference0.07864; replicate167 has five versus
four total components, including the reference-like component. Eligible event
counts are0/1 for Intel/AMD in150 and1/0 in167. This is a decision difference,
not merely different printed decimals. Both historical outcomes remain preserved.

A separate six-row equal-feature fixture exposes another graph defect: dropping
the first nearest-neighbour result retains two self edges when self is not first.
Explicit row-identity exclusion removes them while retaining the requested number
of nonself neighbours. Equal features are permitted; equal cell identifiers are
still rejected by the input contract.

## Candidate and computational evidence

The isolated candidate changes only two operations: use float64 inputs for the
same randomized PCA, and exclude self by row identity before taking k neighbours.
It does not change k, PCA dimensions, Leiden seed/resolution, clone support,
state priors, posterior thresholds, discovery allocation or merge criteria.
Independent distance and unique-distance row-permutation checks pass.

After correction, all nine encountered neighbour graphs across the two complete
replays agree exactly between processors. Final MAP assignments and proposal
eligibility also agree: one eligible event in150 and none in167. Responsibility
values still differ by at most2.52422e−6 and2.37582e−6. These differences are
reported, not accepted through a new tolerance or described as bitwise equivalence.
A future scientific decision affected by such a difference remains unresolved.

Independent audit22969101 verifies the graph identities, self exclusion, MAP
assignments and eligibility from saved arrays. It independently joins frozen truth
and recounts all ten original recovery cases, checks35 recorded input hashes,
410 preserved baseline files and65 unchanged live package/test/native files.
The audit is `dev/hardware_trace_20261108/precision_audit_v1/report.json`.

## Recovery and scientific limits

| Inference path with the precision correction | Original assertions passed | Remaining failed seeds |
|---|---:|---|
| Default thinning path | 2/5 | 20261108,20261109,20261110 |
| Experimental pooled initialization and exact-matched confirmation | 4/5 | 20261109 |

No new recovery is claimed. The remaining experimental failure returns one S1
tumour group, accuracy0.602 and no chr12 distinction. Under the default path,
seed20261108 also has10 normal false calls out of90, failing the original strict
less-than10% criterion. Those assertions are unchanged and remain failed.

The candidate is first implemented as an explicit research hook, then validated
as a noneditable wheel with two added graph-invariant tests and every existing
scientific assertion retained. Job22968997 passes321 Python tests and four Rust
tests; one GPU-only test is skipped on the CPU node. It completes in6m05s; pytest
process peak RSS is2,535,240KiB. The isolated environment copy22968996 completes
in3m27s without new dependencies or changes to shared installations.

Jobs22969581 and22969582 verify exact installed-wheel versus research-hook
responsibility and dose arrays for both locked cases on Intel and AMD separately.
The correction is now applied to live `python/haplodose/clones.py`, with the two
new tests in `tests/test_proposal_graph_precision.py`. All live package, native
source and test hashes match the validated candidate. The wheel identity is
`60fe42a9eabcf9f30a154ea512d5723ee11136b88d230030be136a73f2d7fc15`.
This applies a numerical correction, not a new scientific model or a clone claim.

The running998-fit experiment retains its old wheel and immutable source.
Preservation check22969947 verifies all65 executed package/test/native files and
all five execution-wheel identities unchanged,64 original live files unchanged,
the exact validated replacement of `clones.py`, and the new test hash. The original
pending audit22964195 and dependent tradeoff22968510 were cancelled before execution
and replaced by22969949 and22969969. The new audit records both source versions;
all998-input,1,800-decision,54-bound and paired checks remain unchanged. No running
inference was cancelled, output overwritten or scientific criterion relaxed.

The pending998-fit pooled-initialization study uses its earlier frozen code.
Its eventual result must not be attributed to this numerical correction. Any
scientific promotion still requires error/power evidence for the actual candidate
and independent biological measurements at the claimed resolution. The missing
single-CNA biological validation panel is not supplied by these processor replays.

## Execution records

Stage traces run in jobs22968166/168/192 and22968365/368/370; independent graph
diagnostics in22968575/578; corrected traces in22968754/756/22968806. All complete.
Recovery jobs22968809 and22968810 complete in2m41s and5m36s, each with4 CPUs and8GB
requested. Their scheduler peak RSS is1,113,188KiB and1,131,828KiB. Audit22969101
completes in8s with2 CPUs and4GB requested. Source snapshots, original failures,
processor traces, all arrays and manifests are retained under
`dev/hardware_trace_20261108`. No production dependency or frozen result changes.
