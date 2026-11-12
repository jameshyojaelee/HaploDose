# NextGen RNA/DNA evidence: metadata feasibility, not dataset admission

The scientific goal is to find legitimately accessible, corresponding RNA/DNA
measurements on independent biological units after consuming the fixed
50-unit test. Additional computation cannot create donor independence, assay
correspondence or genomic truth. No measurement matrix is downloaded or fitted
in this investigation.

## Source and two preserved attempts

The primary paper, [Neiswender et al., Nature](https://doi.org/10.1038/s41586-026-10843-7),
links the public [Figshare deposition](https://doi.org/10.6084/m9.figshare.29472362.v1).
The official [26Q3 release announcement](https://forum.depmap.org/t/announcing-the-26q3-release/4692)
also describes additional profiling and changes to relative copy-number
processing. Neither publication certifies independence or suitability for this
package's endpoint.

The official portal metadata endpoint returned a human-verification page in
the browsing tool. No verification was bypassed. Public Figshare API metadata,
its README and its model table were accessible through ordinary HTTP on a
compute node. The README says profiling matrices reside in the DepMap fileset
**NextGen Model Manuscript 2026**. Derived article tables are not substituted
for those primary measurements.

1. **Job 23226521** freezes the metadata-only protocol and fetches the public
   article inventory and README. Its initial filename whitelist does not
   recognize the actual `model_metadata.csv`; an empty candidate list therefore
   means unassessed availability, not absence of new specimens. This first
   attempt and its outputs remain unchanged.
2. **Job 23226622** freezes a separate schema correction, reads the actual
   documented model table and compares identifiers with the old 24Q4 roster.
   No RNA, DNA, segment or CRISPR measurement values are decoded.

[Original protocol](../dev/nextgen_metadata_20261111/DESIGN.md),
[first report](../dev/nextgen_metadata_20261111/run_23226521/results/REPORT.json),
[schema correction](../dev/nextgen_metadata_20261111/DESIGN_v2.md),
[corrected report](../dev/nextgen_metadata_20261111/run_23226622/results/REPORT.json).

## Actual metadata and identification boundary

The table has **1,999 rows**, with **309 unique ModelIDs absent from the
2,105-model old roster**. It supplies lineage, disease, model type, culture
attributes and assay-availability flags. It has no PatientID, RRID,
CellLineName, aliquot or matched condition/profile identifiers. An unfamiliar
ModelID is not an independently new donor or a verified RNA/DNA pair.

The independent summary reader 23226762 confirms that all 1,690 known exposed
ModelIDs lie in the old roster and none of the 309 unfamiliar identifiers
overlaps that exposure list. This check does not identify related donors or
previously exposed cultures under other identifiers. **Zero independent
biological units are qualified by this metadata alone.** This is an admission
result, not a claim that the underlying release contains zero usable units.

| Load-bearing requirement | Current evidence | Decisive next evidence |
|---|---|---|
| Biological independence | ModelID only; donor/family mapping missing | Authoritative donor/culture relationships cross-checked against exposure records |
| Measurement correspondence | Availability flags, no condition/aliquot map | RNA and DNA profile metadata establishing corresponding material and assay |
| Input scale | Matrices not obtained | Release-specific file definitions, transformations and compatible RNA reference scale |
| DNA truth and resolution | No candidate segments or callable masks inspected | Independently constructed DNA support, coordinates, sensitivity and segmentation uncertainty |
| Legitimate measurement access | Public deposition metadata accessible; portal endpoint verification page | Supported public download or authorized export of the required files |

The strongest alternative to an apparent new validation cohort is old or
related biological material under new identifiers. To falsify that alternative,
the mapping must establish independence before evaluating predictions. Likewise,
accurate processed relative-DNA prediction would not identify absolute copy,
cell membership or spot mixtures.

The 26Q3 announcement reports correction of unexpectedly low relative CN values
and low-mappability filtering. **[Unverified] Whether any 24Q4 specimen or locus
used here is affected.** Do not overwrite historical truth or select favorable
loci from this release note. A later version-specific comparison must use an
independent upstream processing/mappability specification, preserve both truth
versions and distinguish a truth correction from a method improvement.

## Identities and compute

The public model table is file 64494438, 272,338 bytes, MD5
`97b26ec52b8a05e52df7abd507b8139f`, SHA-256
`a000a5ee1ebe2b2488a146b2a413ec4b93472488e6ad31b56126f522d652557c`.
The README SHA-256 is
`cdf44320ba3ef22e487547613f2a0470a87fcf107d685de009ea65a913f4c91a`.
The old roster SHA-256 is
`b7a0c1385e6cef30132b56aff61f1261d11e3f490490b355c430d32ee0dbdcfa`.
Downloaded file MD5s agree with deposition metadata.

Both jobs complete 0:0 in three scheduler seconds on nslab CPU. The first
requests 1 CPU / 4 GB / 1 h; the second 1 CPU / 1 GB / 30 min. Submission
headers were reported beforehand. No array, new dependency, inference backend,
shared installation modification or unrelated running-job change is involved.

[First source identities](../dev/nextgen_metadata_20261111/SNAPSHOT.json),
[second source identities](../dev/nextgen_metadata_20261111/SNAPSHOT_v2.json),
[submission one](../dev/nextgen_metadata_20261111/submission.json),
[submission two](../dev/nextgen_metadata_20261111/submission_v2.json),
[independent exposure cross-check](../dev/dosage_response_reserved_20261111/summary_23226762/results/REPORT.json).

No new dataset is admitted, no parameter is selected and none of the consumed
50 units is retuned. The public recovery, false-split and uncertainty failures
remain unchanged. The full scientific goal is active and unachieved.

## Supported metadata route checked while recovery runs

A subsequent primary-source search finds an official
[metadata-inventory endpoint announcement](https://forum.depmap.org/t/provide-an-open-endpoint-for-latest-version-retrieval/4652).
The provider documents an intentionally public inventory route, with download
URLs omitted. Opening that exact documented route through the browsing tool
returns a verification page, rather than inventory JSON. No verification is
bypassed, undocumented storage path guessed, or measurement matrix opened.

The official [organoid availability discussion](https://forum.depmap.org/t/organoid-dependency-data/4671/7)
states that the temporary NextGen dataset was integrated into the 26Q3 release.
The earlier model table therefore should not be treated as the final word on
whether richer metadata exists. However, neither that announcement nor the
verification response supplies donor, culture or RNA/DNA correspondence for
the 309 unfamiliar identifiers. Qualified independent units remain zero under
the current evidence. A supported metadata export remains needed before any
new-cohort admission; existing consumed units remain regression-only.
