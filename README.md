# HaploDose

Infer relative RNA dosage and candidate clones from single-cell and spatial RNA-seq.

- Detect genomic regions with increased or decreased expression relative to reference cells.
- Group cells by shared RNA dosage patterns; genomic clone identity requires independent DNA evidence.
- Add phased SNP counts to detect allelic imbalance, conditional on donor, phase and mapping integrity.
- Use spatial neighbours to smooth candidate assignments in spatial data.
- Export cell assignments, genomic segments, clone profiles and QC; simulate data for testing.

Outputs are RNA-based estimates, not validated DNA copy-number calls.
Version 0.2.0 introduces explicit input declarations and schema-2 evidence records.
See the [scientific contract](docs/SCIENTIFIC_CONTRACT.md),
[migration guide](docs/MIGRATION_0.2.md), [status](docs/STATUS.md) and [model details](docs/related_work.md).
The proposed 5% specimen false-split guarantee has not been biologically established.
Development tests found substantial false grouping under correlated expression state and
undercoverage of approximate RNA-effect intervals after selection. The original thinning-mode
recovery regression is repaired, but recovery remains fragile across independent simulation draws.
See the [membership and merge investigation](docs/MEMBERSHIP_2026-11-07.md),
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
Optional: phased SNP counts (`--alleles`). Run `haplodose run --help` for options.
The declaration must describe the source measurements: rounded normalized values are not raw counts.
Choose the coordinate system and feature namespace matching the input files. New outputs use zero-based,
half-open genomic intervals. Use a fresh isolated environment and output directory for a new version.
