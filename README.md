# HaploDose

Infer relative RNA dosage and candidate clones from single-cell and spatial RNA-seq.

- Detect genomic regions with increased or decreased expression relative to reference cells.
- Group cells into clones and subclones with shared dosage patterns.
- Add phased SNP counts to detect allelic imbalance, including candidate copy-neutral LOH.
- Use spatial neighbours to support clone assignments in spatial data.
- Export cell assignments, genomic segments, clone profiles and QC; simulate data for testing.

Outputs are RNA-based estimates, not validated DNA copy-number calls.
Under development: the latest documented test run has two failures;
see [status and limitations](docs/STATUS.md) and [model details](docs/related_work.md).

## Install and run

In the existing `haplodose` environment (Python 3.11+, Rust build tools):

```bash
micromamba activate haplodose
pip install --no-user --no-build-isolation --no-deps -e .
haplodose run --counts counts.h5 --cells cells.tsv --positions positions.tsv --out OUT
haplodose simulate --out SIM
```

Inputs: gene counts (10x H5, h5ad or CSR H5), cell metadata and genomic gene positions.
Optional: phased SNP counts (`--alleles`). Run `haplodose run --help` for options.
