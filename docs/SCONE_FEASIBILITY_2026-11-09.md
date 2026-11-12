# scONE primary metadata pass; specimen-level DNA access unresolved

The existing catalogue candidate scONE-seq remains unadmitted. Two compute-node
jobs verify primary archive metadata and five author-source identities without
opening RNA/DNA measurements, per-cell annotations, clone labels or reads. The
result narrows the access gap; it does not validate a method or supply independent
single-CNA controls. The abstract and published biological descriptions were
visible during literature review, so the literature review is not fully blinded.

## What the archive establishes

The [primary study record](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA768428)
links GSE185269 and PRJNA768428. The collected GEO self-series metadata contains
621 distinct GSM identifiers. ENA returns621 runs with621 distinct sample aliases
and621 sample accessions. Literal GEO/ENA alias sets agree exactly; this is an
administrative reconciliation, not proof of same-cell molecular correspondence.

All run titles fall into six literal patterns after masking standalone numeric
indices: HCT116 scONEseq98, extracted HCT116 nuclei98, live whole HCT116 cells92,
HUVEC+H9 mixtures67, NPC43 scONEseq170 and Smart-HCT96. None contains the recorded
brain/GBM/glioma/astro/patient terms. Generic source tags are all TRANSCRIPTOMIC,
while525 strategies are OTHER and96 RNA-Seq. Neither tag alone establishes the
actual assay layer; combined DNA/RNA libraries require protocol interpretation.

GEO lists GBM RNA, annotation, normalized scWGS depth and mBAF supplementary URLs.
Those files were not opened. [Unverified] Raw DNA access for the brain specimen
has not been established. Absence of a matching run title cannot establish that
raw data are unavailable elsewhere. This gap prevents claiming an independently
reconstructable brain DNA measurement chain from the records inspected so far.

No runs, cells, mixed libraries or sample IDs are counted as independent donors.
The GEO narrative describes one frozen brain specimen. No biological null
specimen is certified here; the whole-specimen5% bound still lacks suitable
independent replication. Current biological null n remains0.

## Consequential source assumptions

The [author DNA instructions](https://github.com/0YuLei0/scONE-seq-data-processing/blob/main/DNA_related/README.md)
describe initial assumed-diploid calls,5Mb allele-frequency information mapped
onto500kb depth bins, and a sample-specific tetraploid correction. A500kb depth
matrix therefore cannot automatically validate500kb allele resolution, and
normalized depth times two is not independent absolute-CN truth.

The inspected [RNA workflow](https://github.com/0YuLei0/scONE-seq-data-processing/blob/main/RawDataProcessing/4_RNA_UMI_quantification.sh)
uses UMI deduplication followed by Kallisto/Salmon quantification. That does not
establish deposited values as raw observed integer counts. A filename saying
gene counts, or rounding values, cannot clear the raw-count model contract.
Actual executed references, tool versions and observation scale remain necessary.

Five source contents reproduce their Git blob identities exactly, including
retained line endings. This establishes the fetched files, not the source version
executed for published measurements. External scripts were not run. A search
snippet exposing a competitor's repackaging was seen; none of its rounded-count,
normalized-depth or label construction was adopted as primary truth.

## Decision and smallest next evidence

Do not admit the brain specimen, construct truth from author clone labels, change
the inference method, or report a new validation pass. The missing evidence is
an explicit specimen/library/cell mapping to raw DNA access and executed molecular
separation/reference settings. If supplied by an archive or authorized local
resource, freeze that correspondence before bounded upstream inspection and
independent DNA calling. DNA coverage, segmentation/ploidy ambiguity and stability
must determine the supported resolution without selecting it to match RNA.

Cell-line raw libraries are a separate potential challenge, not additional brain
patients or certified single-clone nulls. Any later development/evaluation use
requires a frozen role and measurement contract before opening outcomes. Existing
DEFND/OneK1K resource questions remain unresolved; this check does not replace them.

Computational correctness is supported for source identities, metadata counts and
literal joins. Statistical identifiability of genomic effects and biological
validity of CNA membership are not established. Regional transcription and DNA
processing/ploidy ambiguity remain competing explanations. Independent same-cell
DNA at demonstrated resolution, rather than comparator agreement, would address
membership. No result changes the original RNA thinning2/5 or failed error/power
criteria; no SOTA or whole-plan completion claim follows.

## Execution and files

Metadata23029513 completed5s,processpeak86144KiB. Independent audit23029836
completed2s,processpeak81152KiB. Both exit0 on innovation CPU nodes; scheduler
memory sampling reports0K for these short jobs, so process time-v peaks are used.
No new dependencies, production changes or shared installations. No full package
suite rerun for metadata-only scripts; git diff --check passes.

Added scripts: `benchmarks/hardening/scone_feasibility_metadata.py` and
`scone_metadata_audit.py`. Protocols, source contents/identities, isolated execution
inventory, immutable snapshots, exact URLs, raw metadata, reconciled reports and
resource records are under `dev/scone_feasibility_20261109`. No measurement URL was
followed. All frozen results, historical decisions and scientific criteria remain.
