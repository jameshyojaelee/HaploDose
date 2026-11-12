# Bian shared-reference contract

This checks exact existing reference assets and autosomal GTF coordinate
compatibility before RNA count or independent DNA fragment construction. It is
reference QC only. It supplies no specimen truth, normal reference, genomic
profile, clone membership, biological calibration or performance improvement.
Bian remains secondary/descriptive. The candidate/reference/comparator/scorer
measurement freeze and all original scientific criteria remain required.

## Source and exact checks

The read-only UCSC hg19 FASTA,FAI,chromosome sizes and bisulfite index files are
listed in the existing CRC DNA reference manifest. Job23231180 checks all15
assets against their recorded bytes,MD5 and SHA256. It independently scans the
FASTA and compares every contig name and length with FAI and chromosome sizes.
No reference/index files are written or rebuilt, and no CRC workflow is executed.

The existing Ensembl GRCh37 release75 GTF has its original source URL and archive
CHECKSUMS alongside it. Those checksums use BSD sum plus1024-byte blocks, not MD5.
The job independently computes the rolling16-bit sum, verifies the rounded block
count, records SHA256 and reads gzip to EOF. Every autosomal feature must have
valid1-based inclusive coordinates within the exact hg19 contig. Names1..22 map
only to chr1..chr22. Off-autosome features are counted and excluded, with no
guessed mitochondrial,alternate-contig or sex-chromosome aliases. Gene IDs cannot
silently select a different autosome. Per-chromosome feature/exonic-gene totals
describe the annotation; they are not specimen RNA counts.

Passing establishes identity of these exact native hg19 artifacts and coordinate
compatibility of this autosomal annotation. It does not prove every GRCh37 and
hg19 sequence identical or support lifting genomic outcomes between assemblies.
Later interval construction must convert once to0-based half-open coordinates.
Normal-reference composition,DNA callability and biological identification remain
separate assumptions. Existing CRC pilot masks are not adopted as universal
sensitivity criteria by this check.

## Executed scope and dependencies

One job23231180:nslab account/nslab QOS,I/O1CPU2GB1h,no arrays. It streams existing
FASTA/index assets and parses the existing GTF in the isolated Conda/runtime,
with144 source identities before/after. Static AST and shell syntax pass;
executable fixtures/outcomes are pending at submission. Fixtures include known
BSD values,streaming update/block arithmetic,explicit aliases and invalid/exact
coordinate boundaries. Numerical gates are exact integer/cryptographic checks.

The raw staging job23230977 and its dependent independent reader23230978 continue
unchanged. Completed raw-file integrity and this reference contract are inputs to
a future measurement protocol, not permission to infer biology from one cell.
The current candidate wheel is frozen as an artifact but retains its known
recovery,error/power and uncertainty failures. No validated release or SOTA claim
follows. Biological n analyzed=0.

[Frozen design](../dev/bian_reference_contract_20261112/DESIGN.md),
[reference checker](../dev/bian_reference_contract_20261112/check.py),
[source identities](../dev/bian_reference_contract_20261112/SNAPSHOT.json),
[submission](../dev/bian_reference_contract_20261112/submission.json),
[raw-input chain](BIAN_RAW_COUNT_CHAIN_2026-11-12.md),
[primary protocol qualification](BIAN_PROTOCOL_QUALIFICATION_2026-11-12.md).

## Completed reference result

Job23231180 completed0:0 on ne1di6-002 in36s. Timed process36.33s,peak32304KiB;
Slurm sampled12984KiB. All144 source identities pass before and after; all six
fixture groups pass. All15 existing reference assets match the original exact
byte,MD5 and SHA256 values. The independent FASTA scan matches both FAI and
chromosome sizes across93 contigs and3137161264 bases.

The39344043-byte compressed GTF matches its declared BSD checksum743 and38422
blocks,passes gzip CRC, and has SHA256
`770f48d9aecd9bdecdf147fd141c23670b215bccbdbfa459048f45b02581c2ed`.
All22 autosomes have valid in-bounds features and exonic gene identifiers.
The full report preserves per-chromosome counts and all excluded off-autosome
feature counts. No sequence equivalence across assemblies or biological truth
is inferred from this coordinate check. No shared asset was repaired or replaced.

[Executed result](../dev/bian_reference_contract_20261112/run_23231180/results/REPORT.json),
[all15 asset identities](../dev/bian_reference_contract_20261112/run_23231180/results/ASSET_CHECKS.json),
[executed fixtures](../dev/bian_reference_contract_20261112/run_23231180/results/FIXTURES.json).

The raw staging job has independently verified DNA mate1 against archive bytes
and MD5 and is downloading mate2. This is3/4 verified compressed files, not a
completed DNA paired-read or independent-reader result. The same live job and
dependency remain in place; there is no restart or duplicate download.
