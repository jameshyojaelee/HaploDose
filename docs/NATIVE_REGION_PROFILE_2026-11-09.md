# Optional bulk RNA profile on native regions

Implementation and compute-node validation are complete:421 Python/4 Rust tests
pass,1 GPU-only test is skipped, and all68 candidate profiles, adjusted gene
effects and fitted gene expectations reproduce exactly. All11 previous native
arrays/tables remain identical.
The supported task-specific development gain and unresolved pseudobulk comparison
are reported separately in [the performance study](PROGRAMME_GENE_PROFILE_2026-11-09.md).

## Interface and meaning

Bulk API `rna_profile_smoothing='native-regions'`, or CLI
`--rna-profile-smoothing native-regions`, adds a descriptive regional RNA profile
when `bin_stat='sum'`. Median pseudocount mode fails explicitly before creating
an output directory: its fitted expectations cannot be allocated as raw gene
count expectations. Default `none` and existing `predictive` behavior are retained.

The fitted native bin expectations are allocated among genes by the corresponding
within-bin reference counts. Normals exclude themselves from that reference.
Log((gene_count+.5)/(gene_expected+.5)) effects are centred by the median of within-bin
gene medians; ALL genes contribute to each final native-region median. State
labels provide boundaries only. Chromosomes, arms and genomic gaps reset regions.
The full RNA observations supply reference selection, nuisance fitting and regions;
none is independent confirmation or a calibrated DNA copy-state probability.

| Added artifact | Axes | Units |
|---|---|---|
|rna_log_ratio_native_regions.npy|samples.tsv x bins.tsv|natural-log relative RNA rate|
|gene_rna_log_ratio_programme.npy|samples.tsv x genes.tsv.gz|natural-log relative RNA rate|
|gene_fitted_expected.npy|samples.tsv x genes.tsv.gz|input-count units|
|rna_native_regions.json|one record per sample|centres, reference membership counts and region-selection metadata|

`run.json` records formulas, axes, units and unvalidated status. Existing
`gene_rna_log_ratio.npy` remains programme-free. Outputs are float64 memmaps written
one sample at a time. No new production dependency, native kernel, clone inference,
scientific threshold or default is introduced. Reference mismatch and regional
cell-state expression remain alternatives to genomic dosage; corresponding
independent DNA at the declared resolution and biological replication are needed.

## Required computational gates

Independent Fraction allocation and scalar-median fixtures test count conservation,
zero-reference genes, gene permutation, unequal gene counts per bin, state changes,
gap resets and invalid identities. The existing end-to-end predictive test is
retained and additionally parameterized for the new mode. It compares EVERY old
native array/table exactly, checks literal sample IDs, and reconstructs effects
and medians from exported counts/expectations using a separate scalar calculation.

The isolated noneditable-wheel job passes the full Python and Rust suites from a
neutral directory. Its dependent bulk execution reproduces all68 locked research
profiles, adjusted gene effects and gene expectations exactly, as well as every
old native array/table. There is no tolerance relaxation or omitted failed
comparison. GPU parity is not run for this export-only change and no GPU claim
follows. No original evaluation cohort is reused for parameter selection.

## Reproducibility and execution

[Design and evidence](../dev/native_region_profile_20261109/) retain the dirty
baseline, source_v1/v2, exact commands in jobs.json, isolated environment setup,
source hashes, import checks, noneditable wheel/native identity and resource records.
Source_v1 was not submitted; source_v2 corrects a replay input-directory name
before any execution. There is no first-version scientific run to reinterpret.

| Job | Role | Account/QOS | CPUs / RAM / requested time | Outcome |
|---|---|---|---|---|
|23148030|isolated environment|innovation|1 /8GB /1h|complete0:0,4m05s|
|23148031|wheel and full tests|innovation|8 /32GB /2h|complete0:0,2m53s|
|23148032|68-sample exact replay|nslab|4 /16GB /3h|complete0:0,20s|

The pytest process reports peak RSS2,504,796KiB; scheduler batch RSS is2,201,540KiB.
The replay process reports578,796KiB; its scheduler RSS is0 and is not substituted
for the process measurement. Actual CPU backend, imports, dependencies and runtime
are retained in validation_23148031 and analysis_23148032. Frozen results and all
scientific criteria remain unchanged. Independent superiority and biological
calibration remain unestablished; the original three production recovery failures
remain unresolved.

Wheel SHA256: `4e6e5f42ec26911efc9d118c511563c5cf1dfa2026a7456fe39d08765230dd07`.
Native SHA256 remains `e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.
The exact output report is [replay_v1/report.json](../dev/native_region_profile_20261109/replay_v1/report.json).
