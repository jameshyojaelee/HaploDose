# Prospective measurement feasibility: metadata only

The goal is independent validation of RNA-derived profiles and cell memberships,
with biological replication and DNA truth at the claimed resolution. This review
identifies possible same-cell challenges; it does not accept a cohort for
confirmatory use or establish the required single-CNA null population.

## Verified archive leads

| Candidate | Verified access/measurement metadata | Potential use and unresolved gate |
|---|---|---|
| DEFND-seq, GSE224149 | ENA returns60 runs across13 sample accessions. The authors' Zenodo archive lists DNA fragment files for two GBM patients, paired with GEO entries GSM7016419 and GSM7817780. | Possible same-nucleus membership/profile challenge. The queried ENA titles identify GBM RNA libraries; corresponding raw GBM DNA read access was not established by this manifest. Fragment-level provenance, barcode correspondence and DNA resolution need verification. |
| Hybrid BAG-seq, PRJNA773107 | ENA returns38 runs:18 labelled OTHER,10 RNA-Seq and10 WGS. The paper describes five uterine-cancer patients and DNA/RNA capture from the same nuclei. | Possible same-nucleus challenge with separate DNA-only/RNA-only controls. Libraries and mixed-patient pools are not independent donors. Molecular-layer separation, collisions, independent DNA calling and callable resolution remain unresolved. |
| HIPSD&R-seq, EGAS50000000691 | The paper deposits raw data at EGA and describes modified whole-genome DNA/RNA co-assays. | Access is unverified. Requires modality-specific interpretation and sensitivity checks for sparse DNA/metacell truth. No access request or data fetch was made. |

Sources: [DEFND-seq study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12272635/),
[author-deposited GBM fragment metadata](https://zenodo.org/records/13984673),
[hybrid BAG-seq study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12476631/),
[HIPSD&R-seq study](https://link.springer.com/article/10.1186/s13059-024-03450-0).

The DEFND manifest mixes actual ATAC controls, modified DNA assays and RNA
libraries. A GENOMIC tag or use of a Multiome kit does not establish whole-genome
DNA measurement. Each library must be assigned from its experimental protocol.
The primary DEFND report also identifies possible residual accessibility bias;
DNA copy profiles must be checked independently of RNA-defined groups.

Hybrid BAG separates RNA- and DNA-origin templates using exonic versus strictly
intergenic mapping. [Inference] This makes molecule classification and the
separate DNA-only controls consequential for truth independence. Published
co-clusters must not be copied as certain clone truth. The metadata include
aliases such as1+5tdr01 and2+5tdr01; donor assignment requires explicit verification,
not interpreting each mixed library as another specimen.

## Decision and next admissible actions

A later [OneK1K archive check](ONEK1K_FEASIBILITY_2026-11-08.md) verifies public
RNA/DNA series metadata and a possible source of donor replication. It does not
certify any single-CNA null specimens: donor correspondence, blood/PBMC assay
differences and DNA detection limits remain unresolved. No measurements were opened.

The follow-up [protocol/source review](PROSPECTIVE_SOURCE_REVIEW_2026-11-08.md)
identifies the specific remaining gates: hybrid BAG's RNA-derived hotspot mask
and broad DNA bins, and DEFND's artifact-specific barcode/reference identity and
fragment multiplicity convention. Six source blobs and three repository commits
are recorded. No per-cell data were opened and neither cohort was admitted.

Neither two GBM patients nor five uterine-cancer patients can provide59 independent
single-CNA null specimens. Repeated libraries, thousands of cells and mixed-donor
pools do not close this replication gap. These studies might support a narrower
prospective membership challenge after the method and evaluation procedure are
frozen, while the population false-split guarantee remains withheld.

Before accepting either public study, establish the exact RNA/DNA barcode chain,
assembly, raw-count provenance, independent donor identities and DNA detection
limits. Construct DNA truth without RNA grouping. Assess coverage, resampling,
breakpoint uncertainty and assay-specific biases at the intended resolution.
Do not choose specimens or resolution based on HaploDose/comparator outcomes.
No raw reads, expression matrices, fragment files or clone assignments were
downloaded or inspected in this investigation. No inference was run.

Published abstracts, methods and availability snippets were read, including
some authors' biological findings. This literature review is therefore not
entirely outcome-blind. Own-method prospective outputs and per-cell data remain
unopened. Prior Schmid, CCLE, MIX-seq, CPTAC, SPECTRUM and REBC evaluations stay
regression-only. Bian retains its earlier descriptive-only decision; the user's
withdrawal of Tahoe remains in force.

## Reproducible metadata record

Compute job22964928 completed in10s. `prospective_metadata.py` requests only the
ENA read_run metadata endpoint and verifies available fields before querying.
The returned URLs are recorded but not followed. Responses are capped at10MiB.
The manifests, query URLs, checksums, unresolved eligibility fields, source snapshot
and execution environment are under `dev/prospective_metadata_20261108/metadata_v1`.
Run/sample accession counts are explicitly not converted into biological n.
