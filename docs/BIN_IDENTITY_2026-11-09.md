# Shared gene starts must not change named bin counts with column order

The pre-correction wheel changes regional RNA counts when two distinct genes
share a genomic start and their input columns exchange order across a bin boundary.
The named observations and annotation are identical. This is computational
identity failure, not evidence of a biological dosage change.

## Counterexample and correction

Probe23054521 executes all720 orders of six named genes on a CPU compute node.
Genes A/B/C/D/E/F start at1/2/3/3/4/5; bins contain three genes. The prior sorter
uses original count-column index after chromosome/start. Half the permutations
give A/B/C versus D/E/F; half give A/B/D versus C/E/F. With counts2/3/5/101/7/11,
the regional observations become10/119 or106/23, although the total stays129.
Two independent raw integer count vectors retain exact total conservation in
every case. All720 unique-start controls have one identical named assignment.
Annotation is deliberately reversed, ruling out its row order as the tie rule.

The correction sorts by chromosome rank, genomic start and canonical unique gene
ID. Count-column index still links genes to their observations; it no longer
decides tied genomic order. No gene is dropped, merged or relabelled. Bin size,
mean filter, arm/gap boundaries, likelihood and scientific thresholds are unchanged.

Seven new test cases cover count-column permutations, unequal ends of genes with
the same start, independently reordered annotation, exact sparse aggregation,
identity-aligned mean filtering and explicit genomic gap resets. Every original
test remains unchanged. This provides reproducibility of a declared binning
rule; it does not claim that lexical gene order reveals finer biological resolution
or that overlapping genes supply independent genomic observations.

## Independent validation gates

Candidate_v1 preserves the pre-correction source and adds only the tie-order
correction and the new test file. Noneditable wheel validation23054611 runs in
a separately copied dependency environment after201-second I/O job23054610.
The copy's unstarted memory request was reduced16GiB to8GiB because both I/O
nodes had16160MB available; its code/dependencies and all test resources remain
unchanged. Build/test allocation is CPU8/32GB/2h. No shared installation changes.

Independent audit23054644 reconstructs all1440 cases by gene-ID lookup and exact
integer summation, compares unique-start controls and verifies frozen source/test
identities. Full public-workflow job23054815 uses one fixed existing development
simulation with a declared G00009/G00010 start tie. It replays original/reversed
column orders on the old and corrected wheels. The original must exhibit changed
bin counts; corrected scientific arrays/tables must agree exactly, including
missingness masks. Input-column metadata and run provenance remain distinct.
The four simulated fits add no biological replication or error-calibration evidence.

Validation23054611 completes in139 seconds:372 Python tests and all4 Rust tests
pass;1 GPU-only test is skipped. Python process peak is2448340KiB; scheduler peak
is2126404KiB. The noneditable wheel SHA256 is
`fc1ec44534629801f08bd4b84bfc326a71ce67e65dae9120065387b329408a60`;
native SHA256 remains
`e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.
Installed package/native paths remain inside the private environment.

Audit23054644 completes in14 seconds, process peak85268KiB, exit0. All1440
cases agree with the independent gene-ID/count reconstruction. All720 corrected
tie cases have one named assignment; all720 unique-start controls exactly retain
their historical membership and counts. Every original test and all593 frozen
baseline source files are preserved. Exactly one existing source file changes,
and only the new test file is added to the candidate inventory.

Workflow23054815 completes in85 seconds, process peak918792KiB, scheduler peak
784304KiB, exit0. All four fits verify isolated source/wheel/native identities and
the actual CPU backend inside inference. The original wheel has different bin
counts under column reversal. Corrected bin-count/posterior arrays have identical
shapes, dtypes, missingness and every observed value. All eight scientific tables
listed in the report agree exactly, as does gene metadata after identity alignment
and removal only of input-column index. This is public-workflow equality, not
matching groups with a favourable relabelling or widening a numerical tolerance.

Exact equality is the criterion for named identities and count observations; no
floating tolerance is introduced. No correctness repair is waived because a
historical metric changes. The separate spatial-wheel GPU job23054254 does not
validate this new wheel. Native kernels are unchanged; no new GPU parity claim
is made for the bin correction.

## Preservation and interpretation

Artifacts and source snapshots are in `dev/bin_identity_20261109`. The old wheel
and earlier frozen evidence remain unchanged. Tied historical bins may change
under the correction; new analyses require new result identities. No historical
output is rewritten. Unique-start behavior must remain exactly identical.

Implementation is `python/haplodose/bins.py`; new tests are
`tests/test_bin_identity.py`. Diagnostic scripts are
`benchmarks/hardening/bin_identity_probe.py`, `bin_identity_audit.py` and
`bin_identity_workflow.py`. Designs, installed wheel/native identities, runtime,
memory, actual CPU backend and dependency status are recorded with their jobs.

The estimand in this check is the regional sum of the same named raw observations.
Its observation model is exact integer aggregation. Incompatible named counts
under a pure column permutation falsify correctness; no stochastic tolerance or
biological assumption repairs that failure. Whether the corrected regional RNA
effects identify genomic CNAs still requires appropriate references, nuisance
control and independent DNA at the claimed resolution. This correction does not
resolve the original three recovery failures, uncertainty failures or the failed
false-split/power experiment. It supports no superiority or SOTA claim.
