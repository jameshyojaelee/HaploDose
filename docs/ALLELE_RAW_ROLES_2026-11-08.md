# Stable allele roles precede public raw-input discovery

A research front end now partitions whole SNP sites before the unmodified public
reader and complete raw-input pipeline. Eight public runs preserve identical
scientific outputs within each of the two existing cell models when only reserved
observations change. An independent audit verifies site roles, all row/count
conservation and actual emitted tables/posteriors. This establishes tested
computational isolation, not independent confirmation or clone accuracy.

## Measurement and identifying assumption

The estimand is dependence of discovery artifacts on designated reserved allele
observations, conditional on the same raw RNA counts, reference cells, assembly
and site metadata. The existing NB RNA and phased allele models are unchanged.
Shared ASE, donor phase, mapping bias and reference uncertainty remain possible
sources of statistical dependence between physically separate SNP observations.
Corresponding DNA allelic measurements would distinguish a genomic contrast
from transcriptional imbalance; this computation cannot.

The role hash uses assembly, canonical chromosome, positive one-based position,
fixed seed and schema version. Discovery occupies exactly half the BLAKE2b hash
space, without rebalancing observed counts or site numbers. AD, DP, phase,
coverage, gene, site order and other sites cannot alter an existing coordinate's
role. All rows at that coordinate receive the same role. Explicit biallelic
REF/ALT and map/gene metadata are required by this research interface; conflicting
annotations fail. It does not verify supplied alleles/assembly against a FASTA.
The historical public input reader is unchanged.

## Locked perturbations and exact results

The fixture contains 480 genes, 240 reference cells, 300 query cells and 520 sites
across four chromosomes. Its simulated populations are a computational challenge,
not independent biological replication. The baseline input contains 117445 rows,
151488 alternate reads and 299181 total reads. The splitter preserves 60403 rows
for discovery and 57042 for reserved evidence, with no loss or duplication.

Four contexts were fixed before fitting: baseline, reserved-only AD reversal,
reserved-only AD/DP zeroing, and one new reserved site at chr1:1. Perturbations
include reference rows. All original roles remain unchanged and all discovery
TSVs are byte-identical. A metadata-only lexical search selects the added site;
neither SNP retention nor method outcome selects it.

All eight public pipeline calls use only the discovery TSV, with declared raw
counts, hg38 and one-based closed expression positions. CloneParams retains its
defaults except fixed seed and thinning mode. This includes gene mapping/binning,
ten programmes, reference fitting, gene vetoes and the complete specimen search.
The fixture has no spatial graph. Cell evidence is separately binomial or explicit
beta-binomial concentration 2; outcomes do not select between those models.

Each run selects 480 genes, 48 bins and 257 loaded SNPs. Within each model all
16 decompressed scientific tables/confirmation records, posterior arrays and
numerical QC are identical across the four contexts. Only recorded input paths,
timing and resources are excluded from the QC identity comparison. Both baseline
models assign 300 query cells to groups of 60/120/120, including the reference-like
group. This is not cell-level DNA validation or a recovery benchmark verdict.

## Independent audit and decision

The audit recomputes role hashes using hexadecimal conversion, reconstructs the
multiset of every input row and independently totals integer AD/DP. It verifies
partition disjointness, exact conservation for all contexts, site addition/order
stability and phase-only role stability. Duplicate headers, conflicting ALT,
fractional position and identical REF/ALT are rejected. Actual emitted scientific
tables and posteriors are compared independently of the worker's report, and all
recorded public-run input/output hashes are checked.

Preparation job 23003081 completes in 11s. Binomial job 23003088 completes in 31s
with 880100 KiB process peak RSS; concentration-2 job 23003089 completes in 29s
with 876668 KiB. Audit 23003140 completes in 8s with 277516 KiB. All exit zero on
compute nodes. Immutable evidence is under
`dev/allele_raw_roles_20261108/{prepare_v1,binomial_v1,beta_v1,audit_v1}`.
The audit source_v1 was preserved unexecuted; source_v2 removes an unnecessary
conditional context manager before execution. No criterion was changed.

Retain this research implementation. Production still uses its original full-SNP
discovery path and retains the previously demonstrated local dependency. Reserved
files have no fitted significance or clone acceptance interpretation. Role
isolation alone does not establish statistical independence: shared-cell ASE
must be challenged before an independent-confirmation claim. Observed and null
search mechanisms must still match for any joint-channel error-control procedure.

No package/default, frozen artifact, scientific threshold or preexisting test
changes. The validated wheel remains at 348 Python and four Rust passes, with
one GPU-only skip; that unchanged suite was not rerun. Original thinning remains
2/5 with three unresolved failures, the broader joint error/power gate remains
failed, and certified biological null n remains zero.
