# Spatial neighbour ties now retain cell identity

GPU follow-up23054254 completes exit0 in65s on RTX PRO6000 Blackwell. The existing
CPU/GPU simulated-workflow test passes with0 skips and actual backend assertions
inside both inference processes. Process peak1655296KiB, scheduler peak1779264KiB,
sampled owned GPU memory870MiB;5s sampling can miss between-sample peaks. Wheel
SHA256 fd507a0bc21643547084e47ebb8462fc0eddd2f6f11849a6dc5840ceade66546
and native identities are verified. Reports are under
`dev/spatial_identity_20261109/gpu_parity_v1`. The simulation is nonspatial and
supplies no spatial biological validation. It does not cover the later bin-identity
wheel, whose CPU validation is reported separately.

The public pipeline constructs specimen spatial graphs in canonical cell-ID order
and restores the input order afterward. Reordering identical cells/coordinates no
longer changes the graph solely through distance-tie handling. Tied neighbours
have no known genomic preference; this corrects computation without validating
biological CNA membership or improving the statistical observation model.

## Counterexample and correction

Probe23030369 finds70/80/22/8 changed directed entries on the same64-cell8x8 grid
under row permutation at k=1/2/4/8. A60-cell triplicated-coordinate fixture changes
58/0/156/0 entries. All eight cases and zero-change controls remain. Canonical-ID
search yields exact invariance, zero self edges and exactly k valid nonself
neighbours. Independent dense NumPy distances verify nearest-neighbour geometry.
No counts or clone outcomes select k, tie handling or a scientific threshold.

`neighbour_graph(xy,k,cell_ids=...)` validates aligned, unique, nonempty string IDs
before indexing. The public pipeline supplies its specimen cell IDs. Helper calls
without IDs preserve legacy behavior. Untied fixture geometry is unchanged.
Reproducibility is conditional on the recorded nearest-neighbour implementation,
not a promise of identical tie choices across arbitrary software/hardware.

The computational estimand is the identity-aligned graph under fixed coordinates
and k. Exact permutations and independent distances establish this correction.
Falsifiers were a failed old counterexample, non-invariant candidate edges, self
edges, non-nearest geometry, accepted ambiguous IDs or existing-test failure.
None occurs in the executed validation. Biological n=0 in these fixtures.

Spatial assignment remains a row-normalized smoothing heuristic conditional on
RNA likelihoods and graph choice. Tissue-state gradients and normal admixture
remain alternatives to genomic clones. Matched-region DNA at demonstrated
resolution and boundary/mixing challenges, rather than reproducible smoothing,
would address those biological explanations. Adaptive lambda is not calibrated
false-split control. The original three RNA recovery failures, failed error/power
criterion and missing biological null panel remain unresolved.

## Compute-node evidence and preserved artifacts

Probe23030369 completes13s exit0. Environment23030649 completes3m18s exit0 on
an I/O node, copying existing dependency versions into a new prefix. No new
dependencies or shared installations. Validation23030652 completes2m43s exit0
on ne1dc6-007,nslab8CPU32GB:365 Python pass,1 GPU-only skip,4 Rust pass. The full
existing suite and scientific assertions remain unchanged;14 new identity cases
are added. Pytest126.39s,time-v peak2497968KiB,scheduler peak2252632KiB.

Wheel SHA256:
`fd507a0bc21643547084e47ebb8462fc0eddd2f6f11849a6dc5840ceade66546`.
Unchanged native SHA256:
`e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.
GPU23024875 remains queued on its prior wheel; it does not certify this new path.

Source-manifest comparison permits only two existing package changes:
`python/haplodose/clones.py`,`python/haplodose/pipeline.py`. Added test:
`tests/test_spatial_identity.py`; every older test is byte-identical to baseline.
Added diagnostic:`benchmarks/hardening/spatial_identity_probe.py`.
Migration,scientific-contract,status and claim-ledger documents record the result.
`git diff --check` passes.

Build command:`slurm/hardening_validate.sbatch dev/spatial_identity_20261109 dev/spatial_identity_20261109/candidate_v1 dev/spatial_identity_20261109/environment`.
Exact allocations/commands, protocols, baseline/candidate snapshots, before/after
fixtures, dependency inventories, wheel/test/resource records are under
`dev/spatial_identity_20261109`. Frozen results and historical decisions remain
unchanged. Corrected tied-coordinate reruns need new output identities; old
analyses are not retroactively edited or declared scientifically validated.
