# Haplodose

Copy-number clones from single-cell and spatial transcriptomes. Rust streaming
aggregation plus JAX kernels; the same code runs on CPU threads and CUDA.

Status (2026-11-01): the expression model, the allele model (phased SNPs, cnLOH) and the
spatial Potts prior are implemented. The working tree currently fails 2 tests and the
sections below describe M1 only. See `docs/STATUS.md` for the current state, decisions,
results and open work. Outputs are relative RNA dosage states and RNA-inferred clones,
not absolute DNA copy number or DNA clones.

## Model in one paragraph

Per-cell counts are summed into fine genomic bins (default 10 expressed genes)
in one streaming pass. Each specimen's pooled profile is compared with a mixture
of reference groups (Poisson-lognormal: sampling noise plus a bin-level baseline
mismatch tau). Candidate groups come from Leiden on windowed, depth-regressed log
ratios. Groups merge unless their pooled log-rate difference contains a resolved
contiguous copy-number step; that difference does not depend on the reference.
Cells are then assigned to clone profiles by negative-binomial likelihood (EM),
clones are re-segmented, and each clone is searched once more for subclones.
Proposals use even bins and difference tests use odd bins. Clone `_c0` is copy
neutral relative to the reference (normal-like). Sequencing depth shifts bin
composition (gene-dense chromosomes most), so a per-bin depth slope, estimated
within clones, enters the cell likelihood, pooled expectations and difference
tests; without it, real SMC specimens split by depth. See `docs/related_work.md`.

Known limits (M1): subclone dose is overestimated after a split (cells assigned to
the subclone are enriched for strong signal); per-cell assignment between clones
that differ by one small event is limited by depth (d' = shift x sqrt(expected
UMIs in the event)); discrete states come from a fitted dosage slope beta that
can sit at its bound.

## Install (cluster env)

```
micromamba activate haplodose     # python 3.11, jax 0.10.2 (+cuda12 plugin)
pip install --no-user --no-build-isolation --no-deps -e .   # builds the Rust extension
```
Without the extension, NumPy/SciPy fallbacks give identical results, slower.

## Run

```
haplodose run --counts counts.h5 --cells cells.tsv --positions positions.tsv \
  --out OUT [--exclude-genes mask.txt]
haplodose simulate --out SIM
```
`cells.tsv` lists cells in matrix order with `cell_id, specimen, patient, role,
ref_group, ref_split` (optional `dataset`). Roles: `tumour_epithelium` or `query`
(fitted), `reference` (`ref_split` fit/holdout), anything else ignored.
Set `HAPLODOSE_PROTECTED_PATHS` (path-separated) to refuse writes under source data.

Outputs: `cells.tsv.gz` (clone, posterior, p_normal, assignment margin, malignant
call), `posteriors.npy`, `clones.tsv`, `segments.tsv` (state, log dose, posterior,
min detectable log dose), `profiles.tsv.gz`, `differences.tsv` (resolved segments
between clone pairs), `tree.tsv`, `log.tsv` (merges/splits), `reference_check.tsv`
(held-out normals called aneuploid), `specimens.tsv`, `bins.tsv`, `genes.tsv.gz`,
`bin_counts.npy`, `qc.json`.

## Tests

`slurm/build_test.sbatch` (compute node): `cargo test`, editable build, `pytest`.
