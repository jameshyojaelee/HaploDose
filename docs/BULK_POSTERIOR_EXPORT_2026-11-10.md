# Retained bulk profile now runs directly from counts through the CLI

The optional `posterior-regions` bulk export delivers the retained RNA profile
directly from raw counts. All48tumour profiles and68gene-effect/expectation
contexts match frozen predictions EXACTLY. Both raw-count control and export
preserve all5base numerical arrays and6tables EXACTLY against frozen native
outputs.480Python/4Rust pass,1GPU-only skip. All4compute jobs complete0:0.
No default, clone model, scientific criterion or frozen result changes.

## Concrete behavior and scientific meaning

Run `haplodose bulk` with `--rna-profile-smoothing posterior-regions`; the
default remains `none`. Raw counts, explicit assembly/coordinates/gene namespace,
normal sample membership and the existing sample/gene validation remain.
The mode requires sum bins:median pseudocounts are rejected before creating
output because they do not provide the corresponding count-unit gene exposures.

Existing sample fits and RNA-state/segment/arm outputs finish unchanged. The
export prepares programme-corrected gene RNA effects and fitted gene exposures
using the identified genomic gene order, then calls the validated public
posterior RNA estimator. Existing chromosome/arm/gap resets define contiguous
blocks. Query tumours receive the additional profile; normals define its
reference and remain NaN/reference-only in that additional array. This avoids
silently querying a normal against its own panel. Normal sample fits in all
preexisting outputs remain exactly as before.

Output:`rna_log_ratio_posterior_regions.npy`,float64,samples.tsv xbins.tsv,
natural-log relative RNA units; `rna_posterior_regions.json` records decisions,
working parameters,numerical diagnostics,reference identity and status. The
corresponding `gene_rna_log_ratio_programme.npy` and `gene_fitted_expected.npy`
use samples.tsv xgenes.tsv.gz; expectations have input-count units. Reference
membership is explicit, with identity tied to this schema2run manifest. The
working smoother records CPU backend independently of native inference backend.
No calibrated DNA state, SE, confidence interval or detection-limit claim is
attached. [README usage](../README.md).

Estimand:regional relative RNA abundance against the declared normal panel.
The observation model is the retained normal-centred, reference-weighted
Gaussian boundary smoother on prepared RNA effects, with robust interval gene
medians and posterior averaging over complete runs. The same bulk beta grid,
tau floor, leave-one-out normal weighting and zero-centre working fit remain.
The ported [public estimator](../python/haplodose/posterior.py) is unchanged.

Identifying assumptions:appropriate reference composition, preprocessing and
regional RNA location. Shared lineage mismatch or coordinated transcription can
produce similar profiles without a DNA change. Corresponding independently
processed DNA at matching resolution, plus independent RNA references, is the
distinguishing measurement. No new biological specimens or independent
evaluation are added. Prepared reference spreads remain working weights, not
identified technical variance or uncertainty coverage. RNA profile agreement
does not validate genomic clone membership, spatial spots, purity,WGD or ancestry.

## Frozen decisive gates, including the complete measurement chain

The [design](../dev/bulk_posterior_export_20261110/DESIGN.md) precedes these runs.
Any numerical/identity discrepancy rejects delivery; no new epsilon, filtering,
feature/reference selection or scoring-rule change is permitted.

1. Full existing Python/Rust suites plus3new integration tests. Tests exercise
   the real CLI parser/forwarding, original output equality, literal001sample
   identity, permuted source genes, an intra-arm genomic gap, independent gene
   effect reconstruction from counts/expected exposures, reference-only rows
   and pseudocount rejection before creation. All480Python tests pass;
   1GPU-only skip and8existing HiGHS forwarding warnings,148.10s. All4Rust pass.
2. Installed CLI control on the original68samples/20normals, exact input file
   hashes and original configuration. All5base arrays and6tables match frozen
   native values EXACTLY. Backend is verified CPU in the inference outputs.
   The old optional predictive-only extra array is not an output of mode none;
   its prior tests remain intact. No failing comparator disappears from a
   performance family because this is delivery parity, not a new ranking study.
3. Installed CLI export on identical raw counts/config. All5base arrays and
   6tables match BOTH the completed control and original frozen values EXACTLY.
   All68programme-corrected gene effects and fitted expectations equal the
   independently retained RNA contexts EXACTLY. All48tumour regional profiles
   equal frozen averaged_gaussian_reference_median predictions EXACTLY;
   all20reference-only rows remain NaN. Full CLI outputs lock without DNA access.

[New integration tests](../tests/test_bulk_posterior_export.py),
[Python report](../dev/bulk_posterior_export_20261110/validation_23158727/pytest.txt),
[Rust report](../dev/bulk_posterior_export_20261110/validation_23158727/cargo_test.txt),
[raw-count control](../dev/bulk_posterior_export_20261110/control_v1/report.json),
[complete count-to-export report](../dev/bulk_posterior_export_20261110/candidate_v1/report.json),
[locked RNA artifact](../dev/bulk_posterior_export_20261110/candidate_v1/RNA_LOCKED.json).

Exact native parity proves that this export does not change the tested native
fits/calls. The68context/48profile equality closes source-count→preprocessing→
reference→gene exposure→RNA profile for these retained inputs. It does not
establish a universal data contract from shape matching or an independent
biological conclusion. All selected identities, source/config hashes, actual
CLI argv and installed code/environment are in the manifests.

## Retained accuracy evidence and unchanged limitations

The delivered single-reference estimator has prior exposed-regression mean
Pearson0.7648167082 versus strongest pseudobulk0.7302613671,23informative
evaluation patients in the original48patient24/24split. Loss recall0.6027473
and FPR0.0437372 are the old conditional outcomes. Exact predictions transport
that estimator to the CLI; this is not a new scoring experiment, prospective
replication or restored uncertainty claim. The prior69-family promotion gate
fails. The separate eight-reference ensemble and its loss-recall cost remain
separate; this CLI option does not silently average reference fits or choose
them using DNA. [Prior profile evidence](POSTERIOR_REGION_AVERAGE_2026-11-09.md),
[public estimator/ensemble delivery](NATIVE_POSTERIOR_PROFILE_2026-11-10.md),
[ensemble uncertainty and recall limits](REFERENCE_BAGGED_PROFILE_2026-11-10.md).

Computational assessment:specified full-count delivery gates pass. Statistical
assessment:conditional descriptive RNA model; prior confidence/superiority gates
remain failed or unestablished. Biological assessment:no new independent DNA
or replication. No default promotion, genomic validation or SOTA completion.

## Executed code, jobs and preservation

Surgical changes extend only [bulk.py](../python/haplodose/bulk.py) and
[cli.py](../python/haplodose/cli.py); the validated posterior kernel and every
old test remain unchanged. One new test file adds3cases. The delivery runner
[bulk_posterior_export.py](../benchmarks/hardening/bulk_posterior_export.py)
calls the actual installed `cli.main(argv)`, not a replacement fitter. Its
recorded argv is the argument parser invocation; no unrecorded shell command
is asserted. Raw counts/metadata/annotations and frozen native/RNA outputs are
read only. A failed control would block the candidate through both afterok
and explicit complete-report assertions.

| Job | Role | Account/QOS | CPU / RAM / requested time | Outcome |
|---|---|---|---|---|
|23158718|private dependency copy|innovation|1 /8GB /1h|complete0:0,4m19s|
|23158727|noneditable wheel/full suites|nslab|2 /32GB /1h|complete0:0,3m10s|
|23158731|full raw-count native control|innovation|1 /16GB /2h|complete0:0,19s|
|23158732|full raw-count posterior export|nslab|1 /16GB /2h|complete0:0,3m36s|

CPU/v6&amd; original005/capturedruntime003 for exact analyses, new isolated
Conda prefix and immutable source, neutral cwd, protected frozen paths. No
arrays,interactiveQOS,GPU,new dependency or shared install. Old posterior,
clone and pending pilot environments/snapshots remain untouched.

Measured process peaks:pytest2,145,692KB, raw-count control573,108KB and
posterior export585,988KB. The environment-copy process peak was not separately
measured; scheduler accounting is retained.

Wheel SHA256 `96d01000def84cb3cd5bbb82424758bec41651691345bb5a120454d8c5284cc2`.
Native SHA256 `e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`,
unchanged. [Installed artifact](../dev/bulk_posterior_export_20261110/validation_23158727/installed_artifact.json),
[all commands](../dev/bulk_posterior_export_20261110/jobs.json),
[accounting](../dev/bulk_posterior_export_20261110/accounting_final.tsv),
[resources](../dev/bulk_posterior_export_20261110/resource_summary_final.json),
[preservation](../dev/bulk_posterior_export_20261110/preservation_final.json).

README/documents add usage and evidence; every169prior ledger values stays.
Frozen scientific artifacts, old tests, scientific cutoffs and clone method
remain unchanged. The independent guarded clone run currently passes4/5full
originals, including108and109;110passes confirmation but awaits final merging.
Fresh phase8false-split/power remains gated on complete5/5. This delivered bulk
improvement does not substitute for those remaining goals.
