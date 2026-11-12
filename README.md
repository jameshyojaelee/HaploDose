# HaploDose

Estimate regional RNA dosage from single-cell, spatial and bulk RNA-seq, and identify candidate clones in cellular data.

- Detect genomic regions with increased or decreased expression relative to reference cells.
- Group cells by shared RNA dosage patterns; genomic clone identity requires independent DNA evidence.
- Add phased SNP counts to detect allelic imbalance, conditional on donor, phase and mapping integrity.
- Use spatial neighbours to smooth candidate assignments in spatial data.
- Produce a regional profile for each bulk sample relative to a normal reference panel.
- Export cell assignments, genomic segments, clone profiles and QC; simulate data for testing.

Outputs are RNA-based estimates, not validated DNA copy-number calls.
Version 0.2.0 introduces explicit input declarations and schema-2 evidence records.
See the [scientific contract](docs/SCIENTIFIC_CONTRACT.md),
[migration guide](docs/MIGRATION_0.2.md), [status](docs/STATUS.md) and [model details](docs/related_work.md).
The proposed 5% specimen false-split guarantee has not been biologically established.
Development tests found substantial false grouping under correlated expression state and
undercoverage of approximate RNA-effect intervals after selection. The original thinning-mode
recovery regression is repaired, but recovery remains fragile across independent simulation draws.
See the [exact confirmation and bulk operator investigation](docs/JOINT_CONFIRMATION_2026-11-08.md),
[confirmation and merge ablation](docs/CONFIRMATION_2026-11-07.md),
[earlier membership investigation](docs/MEMBERSHIP_2026-11-07.md),
[recovery-cause investigation](docs/RECOVERY_CAUSES_2026-11-07.md),
[earlier follow-up evidence](docs/RECOVERY_2026-11-07.md) and
[preserved earlier report](docs/HARDENING_2026-11-06.md).

## Install and run

In a fresh isolated Conda environment with Python 3.11+, the project dependencies and Rust build tools:

```bash
pip wheel --no-build-isolation --no-deps --wheel-dir wheels .
pip install --no-deps wheels/haplodose-0.2.0-*.whl
haplodose run --counts counts.h5 --cells cells.tsv --positions positions.tsv --out OUT \
  --counts-kind raw --genome-build hg38 --coordinate-system one-based-closed --gene-namespace symbol
haplodose simulate --out SIM
```

Inputs: gene counts (10x H5, h5ad or CSR H5), cell metadata and genomic gene positions.
For Ensembl or custom feature IDs, include `gene_symbol` in the positions table
for default ribosomal exclusion. Symbols provide category annotation only;
count IDs remain unchanged. Missing annotation fails explicitly. Single-cell
`--keep-ribosomal` disables this filter when that is the intended analysis.
Optional: phased SNP counts (`--alleles`). Run `haplodose run --help` for options.
For spatial smoothing, supply finite `x` and `y` coordinates in a common frame per specimen.
For independent sections fitted together, supply `section` IDs to limit edges to each section.
Coordinate or section fallback is warned and recorded in `qc.json`.
The declaration must describe the source measurements: rounded normalized values are not raw counts.
Choose the coordinate system and feature namespace matching the input files. New outputs use zero-based,
half-open genomic intervals. Use a fresh isolated environment and output directory for a new version.

For prepared natural-log gene RNA contrasts, an opt-in profile API combines
smoothing scales using conditional masked-gene RNA prediction:

```python
from haplodose.stacking import stacked_rna_profile
profile, diagnostics = stacked_rna_profile(gene_log_ratios, gene_bins, genomic_blocks)
```

`gene_bins` must be complete zero-based integer IDs, and `genomic_blocks` has
one contiguous arm/gap label per bin. Preserve gene identities and order.
The output is descriptive RNA dosage, without DNA states or calibrated intervals.
It improved the exposed single-cell profile regression and failed the bulk
replacement gate; it does not change bulk or clone defaults. See the
[method, full comparisons and negative results](docs/RNA_STACKED_SMOOTHING_2026-11-10.md).

Both profile APIs also accept `gene_segments=True` to add a gene-level median
segment expert. This optional mode further improved the single-cell regression
but also failed bulk transfer. Its working penalty is not calibrated error
control; see the [method and validation](docs/GENE_LAPLACE_PROFILE_2026-11-10.md).

With `gene_segments=True`, both APIs also accept `average_segments=True` to
average fitted segmentations before stacking. This further improved pooled-query
profiles in the exposed single-cell regression. Per-cell DNA accuracy, clone
benefit and bulk/spatial transfer remain unvalidated. Existing defaults stay
unchanged; these descriptive RNA profiles have no calibrated intervals.
See the [complete installed and source-count checks](docs/NATIVE_PARTITION_AVERAGE_2026-11-10.md).

Both APIs accept `solver="dual"` for a smaller linear program with the same
masked-gene loss. The default primal solver preserves historical weight choices.
Equivalent loss optima can give different full profiles when weights are tied.

To start from a cells-by-genes count matrix, use explicit matrix identities and
an ordered feature table (`gene`, `bin`):

```python
from haplodose.stacking import stacked_count_profile
profile, gene_effects, diagnostics = stacked_count_profile(
    counts, cell_ids, gene_ids, features, genomic_blocks,
    query_cells, reference_cells, counts_kind="raw",
)
```

Query and reference cells must be disjoint. Feature-table order defines the
prediction folds; matrix axes are aligned by identity. This produces one
descriptive specimen profile, without assigning cells to genomic clones.
`counts_kind="reconstructed"` is available only for this descriptive API and
does not confer count-model calibration; inferential `run`/`bulk` still require raw counts.
See the [input contract and exact source-count validation](docs/STACKED_COUNT_PROFILE_2026-11-10.md).

For bulk data with prepared gene effects and fitted gene expectations, use the
optional normal-reference profile estimator:

```python
from haplodose.posterior import posterior_rna_profile
profile, diagnostics = posterior_rna_profile(
    query_gene_effects, query_gene_expected, gene_bins,
    reference_gene_effects, reference_gene_expected, genomic_blocks,
    reference_id="normal-panel-v1", counts_kind="raw",
)
```

All arrays must use the same ordered genes and preprocessing; reference arrays
are references by genes. At least three reference profiles are required.
This averages relative RNA profiles over working regional boundaries. Its bulk
profile regression improved over pseudobulk, while broader superiority and
uncertainty remain unestablished. It supplies no calibrated DNA states or
intervals and changes no bulk or clone default. The docstring describes the
prepared-data contract and computational limits.

`haplodose.posterior.average_rna_profiles` combines already mapped profiles using
explicit region and reference-fit identities. Map each fit to common regions
before averaging; averaging bins before taking window medians is a different
estimator. Reference-fit count is not biological replication.
See the [exact delivery checks and scientific limits](docs/NATIVE_POSTERIOR_PROFILE_2026-11-10.md).

To produce the retained profile directly from bulk counts:

```bash
haplodose bulk --counts counts.tsv --samples samples.tsv --positions positions.tsv --out BULK \
  --counts-kind raw --genome-build hg38 --coordinate-system one-based-closed --gene-namespace symbol \
  --bin-stat sum --rna-profile-smoothing posterior-regions
```

The additional `rna_log_ratio_posterior_regions.npy` uses the `samples.tsv`
and `bins.tsv` axes. Normal samples define the reference and have NaN rows in
this additional profile. Existing inference outputs and defaults remain.
Choose declarations matching your source measurements and annotation.
See the [complete count-to-export checks](docs/BULK_POSTERIOR_EXPORT_2026-11-10.md).
