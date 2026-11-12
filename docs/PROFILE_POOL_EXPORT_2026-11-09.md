# Exact profile pools preserve the estimand before aggregation

New single-cell/spatial profiles retain the actual observed and expected pools
used by the final fit. This permits exact log-ratio reconstruction and aggregation
without reconstructing responsibilities from float16. It is evidence preservation,
not a demonstrated accuracy gain. Inference and scientific criteria are unchanged.

## Estimands and exported evidence

The existing Kinker comparison averages clone log effects by MAP cell counts;
pseudobulk uses pooled counts. Nonlinearity, clone-specific centering, exposure
and fractional membership make these different quantities. Neither automatically
identifies DNA dosage. Historical profiles lack the observed/expected pools,
and their saved float16 responsibilities cannot reconstruct them exactly.

New `profiles.tsv.gz` rows add `pooled_count` (full-count observations weighted by
final responsibilities), `expected_count` (conditional reference/nuisance expectation
before RNA-state effect), `log_ratio_centre` and `count_units='input_count_units'`.
Pools need not be integers. The relation is
`log_ratio = log((pooled_count + 0.5)/(expected_count + 0.5)) - log_ratio_centre`.
It applies to bin-summed `log_ratio`, not gene-median `log_ratio_gene`.

For a declared population, sum observed and expected pools before one half-count
offset and one declared centering step. Averaging clone log effects or adding a
half-count per clone targets a different quantity. Include the reference-like
component when targeting the whole fitted population; unfitted cells remain
unrepresented. QC records `profile_pool_contract`; empty profiles have matching
headers. Existing columns/readers are retained. Historical outputs are not backfilled.

These are conditional fitted RNA quantities. Reference mismatch and transcriptional
state remain alternatives to DNA dosage. Fractional pools do not acquire a calibrated
count likelihood by being exported. Corresponding independent DNA at the claimed
resolution, not count conservation alone, would establish biological accuracy.

## Independent computational gates

The new invariant test reconstructs log ratios and independently encloses pools
using original bin counts and exact rounding intervals of saved float16 weights.
Adjacent-representable-value midpoints bound original responsibilities; accumulation
allowances depend on operation count, precision and scale. Rounded weights are
never treated as exact. Full isolated-wheel validation: **385 Python tests and
four Rust tests pass; one GPU-only test skipped**. Original tests are unchanged.

Separate old/new wheel runs use identical simulated raw inputs in fitted and
empty modes. The runner captures full-precision responsibilities, observed/expected
pools and centres at the `summarise` boundary. All captured arrays match exactly;
all96 new populated rows exactly match the captured pools/centres. Independent
matrix multiplication verifies observed pools from responsibilities and bin counts.
Original profile fields and16 other artifacts per mode remain identical: exact
DataFrame comparison for tables, bytes for arrays/other files. QC/manifests acquire
new metadata/source identity and are excluded from that equality claim. Empty
headers pass. This is a computational fixture, not biological replication.

## Compute and preservation

Task: `dev/profile_pool_export_20261109`. `baseline` preserves pre-edit dirty
source; `candidate_v1` is tested; `parity_source_v1` contains the parity runner.

| Job | Account/QOS | CPUs / RAM | Elapsed | Result |
|---|---|---|---:|---|
|23072229|innovation|4 /8GB|3m29s|isolated dependency environment ready|
|23072230|nslab|8 /32GB|2m19s|385 Python +4 Rust passes;1 skip|
|23072668|innovation|4 /8GB|39s|fitted/empty parity passes|

The pending environment job moved from I/O to CPU when Slurm estimated an overnight
wait. It was neither restarted nor duplicated; success dependencies remained intact.
All execution used compute nodes; inference backend was asserted CPU. GPU validation
was not run. No new production dependencies or shared installations were introduced.

Wheel SHA256: `2ba2a00107a89bc0aa899c8afaf274134544c47f287336521cce53008280cd42`.
Native SHA256 remains `e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.
Evidence: `validation_23072230/{installed_artifact.json,pytest.txt,cargo_test.txt}`,
`analysis_23072668` and `parity_v1/report.json`. Jobs/resources are in `jobs.json`.

## Decision and evidence gaps

Keep the exports to preserve primary fitted evidence. Do not promote a new
aggregate or claim better DNA accuracy. The next scientific analysis must compare
matched estimands/operators without DNA-selected weights or tuning.

Kinker values were reconstructed from published CPM, while the current inference
contract requires declared raw counts. Original molecule/library provenance remains
unverified. They were **not relabelled raw or rerun** through that interface here.
A historical lossy reconstruction, or a precise export on another fixture, cannot
settle this measurement-scale gap. Exact historical executed source/environment
and specimen correspondence also remain incompletely established.

Original thinning recovery remains2/5; the joint/full research candidate remains
unpromoted after its false-split/power failure. Membership-confidence and interval
coverage failures remain. Frozen results and historical decisions are unchanged.
No superiority, DNA-membership or biological calibration claim is established.
