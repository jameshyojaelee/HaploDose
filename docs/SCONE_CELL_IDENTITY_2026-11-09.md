# Deposited numbering is consistent; raw-to-matrix biological identity remains open

All98 legacy HCT116 GEO records have unique titles001-098 and unique description
numbers097-194. Each description number exceeds its own title number by96.
This is a cohort-wide naming observation, not evidence of a specimen swap or an
adopted numeric offset. The raw RNA/DNA chain remains linked within each literal
GSM/SRR accession. Its exact biological correspondence to published matrix
columns is a separate claim.

## Question, assumptions and independent evidence

The estimand here is administrative correspondence among accession, run, title,
description and matrix header labels. The observations are literal deposited
SOFT/ENA fields and first-line labels, not RNA/DNA values. Computational
correctness means exact identities, hashes, uniqueness and parser agreement.
Statistical or biological correspondence requires the deposited aliases to refer
to the same physical cell. Competing explanations include different numbering
conventions, source mislabeling or specimen miscorrespondence. An explicit
author crosswalk or independent molecular correspondence could distinguish them;
RNA-DNA profile agreement must not be used to select the mapping being validated.

Frozen DESIGN.md selects all98 previously listed legacy HCT116 scONEseq records.
Collector23050170 retrieves their self-sample quick SOFT metadata, reuses the
preserved pilot and reads only the two deposited cell-line matrix headers.
Every metadata record and both headers are available. No expression/DNA value,
genotype or clone outcome is parsed. Gzip/HTTP buffering can fetch additional
compressed bytes; only the first decompressed line is inspected and retained.
Header hashes are not represented as full-matrix hashes.

Independent audit23050279 reconstructs the selection from ENA using a different
literal parser, checks all98 records and both header files, and verifies the
reused pilot against its original raw-chain hash. There are no duplicate
description ordinals. Both parsers reproduce the offset pattern exactly.

## Explicit candidate mappings, with no automatic adoption

RNA labels include98 HCT116 columns,One.HCT.01 throughOne.HCT.98. The DNA header
includes49 One.HCT columns with the_dedup suffix. The author pipeline's
[DNA mapping script](https://github.com/0YuLei0/scONE-seq-data-processing/blob/main/RawDataProcessing/3_DNA_UMI_mapping.sh)
explains that output suffix, but the inspected author files do not explicitly
mapOne.HCT ordinals toGSM accessions. Current code is not proof of executed code.

The independently exported candidate_crosswalk.tsv retains both title-based and
description-based hypotheses for every record, including absent labels:

| Convention | RNA header labels present | DNA header labels present |
|---|---:|---:|
| Title ordinal |98|49|
| Description ordinal |2|1|

[Inference] Header coverage favors the title naming convention. It does not
certify biological identity. For the lockedGSM5609521 pilot, bothOne.HCT.01 and
One.HCT.97 exist in both matrices, so existence alone cannot distinguish them.
Do not adopt a numeric join by choosing the better RNA-DNA profile match.
Missing DNA columns are missing measurements, not certified CNA-negative cells.

The next frozen SUBMISSION_DESIGN.md requests original submitted filenames and
library aliases without accessing those files. Job23050358 retrieves621 ENA
run metadata records and verifies all98 target run/BioSample/experiment identities
against the preserved table. All98 library_name values are theirGSM accessions;
sample titles/descriptions use the title convention. Submitted filename, format,
size and MD5 fields are empty for all98 targets. This source therefore supplies
no original-filename link toOne.HCT. An empty ENA field does not establish that
original names never existed or cannot be recovered from another primary source.

## Execution and decision

| Job | CPU/account | Seconds | Process peak RSS KiB | Exit |
|---|---|---:|---:|---:|
|23050170 metadata/headers|2/innovation|40|82176|0|
|23050279 independent audit|2/nslab|2|79104|0|
|23050358 submission metadata|1/nslab|2|82944|0|

All execute on compute nodes in the existing isolated validated environment.
Artifacts, full literal fields, candidate mappings, header hashes, immutable
source snapshots, request URLs, native execution inventories and resource logs
are under`dev/scone_cell_identity_20261109`. New diagnostic scripts are
scone_cell_identity.py,scone_cell_identity_audit.py andscone_submission_identity.py
underbenchmarks/hardening. No package source or existing test changes; full suite
not rerun for metadata-only diagnostics. Latest package evidence remains365
Python+4 Rust pass,1GPU-only skip. GPU23024875 remains pendingPriority.

Decision: the deposited dual numbering is verified, and the title-based matrix
mapping is an explicit administrative candidate. Exact raw-to-published per-cell
correspondence remains unestablished by inspected sources. No comparison using
that candidate becomes independent validation. No cohort admission, certified
biological null specimen, CNA truth, default or scientific criterion changes.
Original recovery, false-split/power and uncertainty failures remain unresolved;
frozen results and historical decisions remain intact.

## Original SRA filenames support the locked pilot naming link

Further primary metadata changes the pilot's administrative evidence. Job
23050583 resolves SRX12479606 to exactly one SRA UID, retrieves its full XML,
and verifies SRR16194948 and SAMN22044808. The XML declares two original FASTQs:

| Original submitted filename | Declared bytes | Declared MD5 |
|---|---:|---|
|One-HCT-01_read1.fastq.gz|94184843|fe9e9f4f849000a2a39b560ade9ce763|
|One-HCT-01_read2.fastq.gz|91182090|d564df20763f8bdcaa24eb57a4a0a850|

These are metadata-declared checksums, not independently measured file hashes.
No original FASTQ is downloaded here. They differ from the archive-converted
ENA FASTQ file identities already used in the raw-read chain; filename evidence
does not prove original-versus-converted sequence or quality equality.

The first collector recognized generic `FILE` elements but omitted NCBI's native
`SRAFile` elements. Its frozen XML and initial incomplete inventory remain
unchanged. The helper now recognizes both forms, including namespaced elements.
Three hand-specified schema fixtures provide independently known expected file
attributes. Audit 23050835 uses a separate DOM parser on the preserved actual
XML and matches the corrected ElementTree inventory exactly: two original
FASTQs and two primary ETL records. No new network request is used by the audit.

The original submitted stem `One-HCT-01` supplies a primary naming bridge to
the existing matrix label `One.HCT.01`, under the explicit separator conversion
in SRA_AUDIT_DESIGN.md. This supports the locked pilot's title-based administrative
mapping more strongly than header membership alone. It is independent of RNA
or DNA profile agreement. The description-based `One.HCT.97` candidate remains
in the historical table but is not the submitted filename's stem.

This is one pilot's naming evidence. It does not certify biological identity,
the cohort-wide crosswalk, original-to-ENA base/quality equality, exact executed
preprocessing or DNA truth. No cohort-wide conversion, genotype claim or method
promotion follows. The earlier decision records remain above; the pilot naming
link is now supported, while the stronger biological/whole-cohort claims remain
unestablished.

Both jobs complete with exit 0 in 2 seconds. Process peak RSS is 83720 KiB for
retrieval and 82176 KiB for the independent audit. Original source/output are
`sra_source_v1` / `sra_v1`; corrected helper/audit source is
`sra_audit_source_v1`, with new output `sra_audit_v1`. The XML is 8563 bytes,
SHA256 `04e19587fbcaf82a88fa41160442f80173284470d3b1b27f235617e2be6102d1`.
Current scripts are `benchmarks/hardening/scone_sra_identity.py` and
`benchmarks/hardening/scone_sra_identity_audit.py`. No inference source, existing
test, frozen result or scientific criterion changed. Full package suite not
rerun for this diagnostic metadata parser correction; targeted executable
fixtures and independent actual-record parity pass on compute nodes.
