# The improved RNA profile path now accepts identified count matrices

The successful conditional RNA stack can now be called directly on a count
matrix with explicit cell/gene identities and reference membership. This adds
a usable observation-to-profile path, not a new smoother or DNA-state model.
The final version passes448 Python/4 Rust tests and exactly reconstructs both
locked Kinker cases from their source counts, including gene contrasts, profiles
and selection diagnostics. It selects requested observations before float64
conversion. The initial447-test version and its exact replays remain preserved.

## Public interface and scientific meaning

[Implementation](../python/haplodose/stacking.py),
[count/identity tests](../tests/test_stacked_count_profile.py),
[source-count replay](../benchmarks/hardening/count_stacked_parity.py),
[frozen design](../dev/count_stacked_profile_20261110/DESIGN.md),
[memory correction design](../dev/count_stacked_profile_20261110/MEMORY_DESIGN.md).

`stacked_count_profile(counts, cell_ids, gene_ids, features, blocks, query_cells,
reference_cells, counts_kind=...)` returns a bin profile, ordered gene contrasts
and diagnostics. Counts are cells by genes; dense and sparse inputs work.
Features specify unique `gene` and integer `bin` columns in canonical genomic
order. Matrix axes are matched by identity to these features and requested cells.
Query/reference populations must be disjoint, nonempty and present. Missing IDs
are rejected before indexing; ambiguous duplicates are not collapsed.

The estimand is a descriptive regional RNA abundance contrast against the
supplied reference, under the selected-feature library convention. Pool query
counts y and reference counts r, compute E=r/r.sum()*y.sum(), take
log((y+0.5)/(E+0.5)), median-centre the genes, then call the unchanged
`stacked_rna_profile`. The half-count convention is exactly the previously
frozen benchmark, not a new scientific cutoff or count uncertainty model.
Reference membership is caller-supplied and is never selected to match DNA.

Explicit `counts_kind='raw'` or `'reconstructed'` is required. Selected values
must be finite, nonnegative integers; normalized inputs are rejected. Both
labels produce descriptive RNA effects only. Reconstructed CPM does not acquire
raw-count calibration. Inferential `run` and `bulk` retain their raw-only input
contracts. Diagnostics state natural-log RNA units, unvalidated calibration,
input scale, ordered genes, query/reference identities, selected-feature totals,
normalization and all original smoothing/optimization decisions.

Identifying assumptions for a genomic interpretation remain external: suitable
reference composition, correspondence and RNA responsiveness to DNA dosage.
Regional transcription, lineage/reference mixture and spatial cell-state
gradients can mimic these effects. Independently processed corresponding DNA
at the same resolution and prospective locked replication can distinguish the
claimed genomic benefit. This API does not assign individual cells or spots to
DNA clones, estimate absolute copy number or provide calibrated intervals.

## Correctness and memory behavior

All selected counts and their positive pooled totals must be strictly below
2**53. This is the float64 exact-integer boundary, not a convenient tolerance.
Correctly rounded population totals are checked before sparse floating sums;
positive integer subtotals below this boundary remain exactly representable.
Pooling must conserve those totals exactly. Zero-total populations reject.

The feature table defines fold order, so matrix cell/gene permutations preserve
the same gene contrasts and profile. Unselected observations cannot change the
fit. Sparse data remain sparse until two pooled gene vectors are computed.
Version1 unnecessarily converted the full sparse collection to float64. Version2
first selects the actual query/reference populations and features, then validates
and converts those subsets. It explicitly records validation scope as
`selected_cells_and_features`; unrelated observations are unused. Global matrix
shape, numeric dtype and axis uniqueness are still checked.

Tests use an independently calculated scalar count-ratio reference, exact
query/reference totals13/40, dense/sparse parity, identity-aligned permutations,
unrelated-observation changes, missing/duplicate/overlapping IDs, fractional
bin IDs, invalid scales/counts, zero references and inexact totals. Version2
adds one invariance fixture with nonfinite/fractional excluded observations.
No previous test, numerical tolerance or scientific criterion is changed.

## Source-count chain, not just agreement between wrappers

Both locked cases are fixed before this change:
`2313287_STOMACH` and `42MGBA_CENTRAL_NERVOUS_SYSTEM`. Start from the original
count H5, gene/cell header IDs and frozen inference feature tables. Reconstruct
the exact query-plus-reference ordered input rows; assert original matrix
column indices and reference/query membership, and load those counts through
the existing bounded gene loader. The new public API must then reproduce
every gene contrast and frozen stacked bin value **exactly**, together with
all previous selection weights/diagnostics. No DNA is read by this replay.

Both cases pass version1 with the explicit reconstructed-count label.
[Original count-chain report](../dev/count_stacked_profile_20261110/parity_v1/report.json).
This establishes computational reproduction of the input chain on these cases;
it is not a new independent accuracy cohort. The underlying all181-line
regression improvement remains mean1MbPearson0.707226 versus smoothed
pseudobulk0.667955,177 wins/4 losses, with its existing conditional inference.
[All performance and limitations](RNA_STACKED_SMOOTHING_2026-11-10.md).
Bulk transfer remains rejected and gene-level pseudobulk remains competitive
at arm resolution. No default or clone decision is changed.

## All compute attempts

All five jobs use CPU/v6&amd, isolated wheel/private environment, neutral test
directory and exact source/wheel/native import checks. No arrays, new
dependencies, shared installation or frozen result writes.

| Job | Role | Account/QOS | CPUs/RAM/time | Final status |
|---|---|---|---|---|
|23155729|Fresh exact dependency environment|innovation|1/8GB/1h|complete|
|23155758|Version1 full build/tests|nslab|4/32GB/2h|447 Python/4 Rust pass,1 GPU-only skip|
|23155759|Version1 count-to-profile replay|innovation|1/16GB/2h|both exact|
|23156031|Version2 memory correction/full tests|nslab|4/32GB/2h|448 Python/4 Rust pass,1 GPU-only skip|
|23156032|Version2 count-to-profile replay|innovation|1/16GB/2h|both exact|

[Exact commands, immutable source versions and dependencies](../dev/count_stacked_profile_20261110/jobs.json).
The first successful wheel and outputs remain unchanged when the same private
environment is rebuilt after its jobs finish. All five jobs complete0:0.
Final [full tests](../dev/count_stacked_profile_20261110/validation_23156031/pytest.txt)
pass448 Python tests in133 seconds,1 GPU-only skip,7 visible HiGHS option
pass-through warnings; [Rust](../dev/count_stacked_profile_20261110/validation_23156031/cargo_test.txt)
passes4. No GPU inference is claimed for this CPU-only addition.

The final [source-count replay](../dev/count_stacked_profile_20261110/parity_v2/report.json)
passes EXACT equality for103 query/2118 reference cells across10133 genes in
2313287_STOMACH and517 query/2467 reference cells across10116 genes in42MGBA.
Every frozen gene contrast, profile and prior selection diagnostic matches.
There is no DNA read, new biological specimen, altered threshold or unsupported
clone assignment. [Installed wheel/native/environment inventory](../dev/count_stacked_profile_20261110/validation_23156031/installed_artifact.json)
and task accounting record the actual executed artifacts, backend and resources.

Old smoothing functions remain identical; the package change appends the count
API and adds `tests/test_stacked_count_profile.py`. The independent source-count
runner is `benchmarks/hardening/count_stacked_parity.py`. README, the two new
reports, STATUS, DEVELOPMENT_LOG, claims and the generated three development
indexes document the increment. All previous frozen artifacts, claim entries
and tests are retained, including the successful first wheel.

The clone track is separate: a newly implemented predictive continuous-state
safeguard fails one original recovery and is rejected without adjusting its
threshold. [Complete experiment and all eight jobs](STATE_PREDICTIVE_CONFIRMATION_2026-11-10.md).
The5% whole-specimen confidence requirement, simultaneous power assessment and
biological truth gaps remain unresolved. The active SOTA goal is not complete.
