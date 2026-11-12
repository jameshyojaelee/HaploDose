# DNTR-seq correspondence: metadata audited, archive repair failed

The administrative DNA/RNA crosswalk passes its independent metadata checks.
It establishes neither molecular identity nor clone truth. The attempted repair
of six missing archive mates failed its exact integrity gate on the first file;
the partial and failure are preserved. No production method or frozen scientific
result changed. DNTR remains an exposed regression cohort.

## Question, contract and scope

The estimand is unique declared DNA/RNA library correspondence for an explicit
physical-cell identifier. Observations are primary GEO cell characteristics,
GEO experiment links, ENA run metadata and frozen download records. The identifying
assumption is that these declarations correctly describe cells and molecular
layers. Reused well identifiers, multi-cell wells, mismatched conditions or
unmatched libraries can produce misleading pairs. Independent molecular genotype
concordance, followed by DNA callable sensitivity, would distinguish administrative
agreement from a valid same-cell genomic measurement.

This investigation follows the existing primary-measurement plan. It does not
replace the three recovery failures, reopen held-out evaluation, or supply new
biological replication. Source review found that joint membership/merging had
already reached 5/5 recovery and failed the false-split/power gate; repeating that
weighting variant was not treated as new evidence. Bian remains descriptive-only.

The [frozen metadata design](../dev/dntr_correspondence_20261109/DESIGN.md)
prespecified exact explicit GEO cell IDs, molecular-source layer roles, one RNA
and one DNA library, count1 and matching nonmissing context. No identity was
inferred from filename punctuation. All ambiguous and negative records remain.

## Inputs, code and execution chain

Primary sources were [GEO GSE144296](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144296),
the [GSM4839244 example](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4839244),
the [GEO SOFT endpoint](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE144nnn/GSE144296/soft/GSE144296_family.soft.gz),
ENA PRJNA603321 and the [authors' protocol repository](https://github.com/EngeLab/DNTRseq).
The browser series page presented a challenge; compute-node SOFT retrieval worked.
[Source/exposure notes](../dev/dntr_correspondence_20261109/SOURCE_NOTES.md)
record third-party snippets encountered without accepting processed truth.

The execution chain is:

1. Frozen [run table](../benchmarks/heldout2/dntrseq_gse144296_runs.tsv) and
   [download manifest](../data/external/dntrseq_gse144296/manifest.tsv), with primary
   metadata responses and exact requests recorded in
   [metadata_v1](../dev/dntr_correspondence_20261109/metadata_v1).
2. Immutable source_v1/v2/v3 snapshots under
   [the task directory](../dev/dntr_correspondence_20261109), isolated validated
   [installed artifact](../dev/exact_profile_merge_20261109/validation_23121461/installed_artifact.json),
   and exact submission arguments in [jobs.json](../dev/dntr_correspondence_20261109/jobs.json).
3. Corrected [metadata_v2 report](../dev/dntr_correspondence_20261109/metadata_v2/report.json),
   with all joins, classifications, file metadata and run manifests retained.
4. Independent [audit report](../dev/dntr_correspondence_20261109/audit_v1/report.json),
   independently reconstructed records, fixtures and explicit missing-file list.

The ENA response SHA256 is
`185157fcf5d08185127388573c8bac6ba2fd8dd86b0710355910a6bcbe7432fc`;
the GEO compressed response SHA256 is
`245af63fbbeb3df760ad68c213fd9f76f101866cf7d13cf41b56284532698b65`.
Corrected runs reused verified cached responses rather than refetching a changing
source. Whitelisted metadata were retained; embedded measurement tables were
prohibited. No count matrix, DNA segment, clone assignment or sequencing record
was decoded. The later transfer streamed compressed bytes solely for integrity.

Code: [metadata reconciliation](../benchmarks/hardening/dntr_correspondence_metadata.py),
[independent audit](../benchmarks/hardening/dntr_correspondence_audit.py) and
[bounded archive staging](../benchmarks/hardening/dntr_missing_inputs.py).
Historical executed analysis code/reference hashes for the deposited measurements
remain unresolved; the public protocol repository does not establish them.

## Every attempt, including failures

| Attempt/job | Resources and result | Interpretation and retained evidence |
|---|---|---|
| Metadata v1 /23130012 | innovation CPU,1CPU,4GB,2h; completed25s,0:0 | First classification accepted literal `NA` as context. Retained as a superseded result. |
| Corrected v2 /23130112 | innovation CPU,1CPU,4GB,2h; completed18s,0:0 | Missing-value amendment; exact same metadata,115 pairs now explicitly context-incomplete. |
| Independent audit /23130113 | nslab CPU,1CPU,4GB,2h; completed8s,0:0 | Eleven fixtures, original NA counterexample and corrected behavior; all4869 joins/3232cell keys/9738file records pass independent reconstruction. |
| First staging submission | nslab I/O,1CPU,2GB,2h; rejected before job creation | Inherited `v6` constraint unavailable on I/O nodes. Test-only diagnostic retained. |
| Staging /23130269 | nslab I/O,1CPU,2GB,2h, explicit `amd`; failed17s,1:0 | Exact size/MD5 assertion fails on first mate. No archive accepted; failed run manifest and partial retained. |
| Dependent audit /23130270 | innovation I/O,1CPU,2GB,2h; never executed | `afterok` parent failed; observed `DependencyNeverSatisfied`, then task-owned pending job cancelled. Audit not run. |

Memory accounting reported0 for the first three batch steps and failed staging
batch, and8K for the staging extern step. These are scheduler records, not reliable
peak memory estimates. See [accounting](../dev/dntr_correspondence_20261109/documentation_accounting.tsv)
and the complete generated [run journal](DEVELOPMENT_RUNS.md). No arrays,
interactive QOS, shared installation or production code changes were used.

The [missing-value amendment](../dev/dntr_correspondence_20261109/MISSING_VALUE_AMENDMENT.md)
treats blank/NA/N/A/unknown/not available/not provided/not applicable as missing,
rather than interpreting equal unknown strings as known agreement. The independent
audit reproduces the original error and checks the corrected classification.
It also tests cardinality, context mismatch, multi-cell records, missing IDs and
record permutations. This repairs metadata interpretation, not genomic inference.

## What the primary metadata support

All4869 current ENA runs match the4869 frozen records, with unique GEO joins and
no added/missing runs or changed read metadata. There are3232 explicit cell IDs:

| Classification | Cell IDs | Meaning |
|---|---|---|
| Paired metadata only |1504| Declared one-cell RNA/DNA pair, known matching context; all HCT116. |
| Missing or discordant context |115| Declared single-cell identity pairs, insufficient context certification. Physical identity is not automatically disproved. |
| Not certified single-cell |8| HCA00102.A01/A02/A03/A04/A09/A10/A11/A12 have count2 in both layers. |
| Unresolved library cardinality |1605| Missing/ambiguous partner structure; do not silently label all as RNA-only. |

The1504 context-complete pairs comprise Control192, ETO100uM16, ETO10uM224,
ETO1uM224, ETO30uM48, ETO3uM224, X-ray1Gy192, X-ray3Gy192 and X-ray5Gy192;
the time point is48h. These are cell counts, not independent donor or null-specimen
counts. A control condition does not certify a homogeneous CNA population.

All4869 records declare hg38 and gene-level RNA/bin-level DNA count products.
This is declared assembly and aggregation metadata, not an executed reference
hash, measured callable resolution, raw UMI evidence or genotype concordance.
No layer or cell-line conflicts were found in the independent run mapping.

## Archive gap and unsuccessful repair

Of9738 expected mates, six are absent and lack corresponding frozen manifest
records. Their combined expected compressed size is830,629,012bytes. Existing
9732 files were checked for presence/size and metadata consistency; their full
current contents were not rehashed. Historical checksum records must not be
described as current checksum validation.

The six gaps are SRR12860636 mate1, SRR12860671 mate2, SRR12864211 mate1,
SRR12861527 mate2, SRR12861634 mate2 and SRR12862247 mate1.
[Exact paths, URLs, bytes and MD5s](../dev/dntr_correspondence_20261109/audit_v1/missing_files.json)
were frozen before the [repair design](../dev/dntr_correspondence_20261109/REPAIR_DESIGN.md).
Only a new staging namespace was permitted. Original raw paths, manifest and
historical outputs were not repair targets.

The first transfer, SRR12860636 mate1, produced an802-byte `.partial` instead of
6,765,711 expected bytes and failed the combined size/MD5 assertion. The expected
MD5 is `b8a80fc63b3b53b5d586e97b527ef5bc`. No verified repair mapping exists.
[Unverified] The remote response's nature and failure cause remain unresolved;
neither corruption, archive unavailability nor a network gateway is established.
No checksum criterion was weakened and no retry replaced the failed artifacts.
The independent repair audit remains not run because its prerequisite failed.

## Three separate decisions and next decisive checks

- **Computational:** metadata audit passes; literal-NA correction independently
  verified. Archive restoration fails integrity and remains incomplete.
- **Statistical identification:** only a declared administrative pair is
  identified. Molecular correspondence and DNA contrast sensitivity are unverified.
- **Biological validation:** no clone truth, added independent specimen or
  certified one-CNA null population. No evaluation admission or method promotion.

The smallest archive follow-up is a bounded compute-node response/header and
checksum diagnostic for the first failed URL, preserving the existing partial.
A fresh retry may be accepted only on exact declared size and MD5 plus an
independent checksum audit. This distinguishes transfer failure from accepted
archive content without automatically reprocessing the cohort.

Before membership scoring, freeze representative cell IDs without reading new RNA
outcomes; verify suitable matched genotype evidence and reconstruct DNA truth
without RNA grouping. Determine callable coverage/resolution and sensitivity to
the contrast under DNA resampling. Missing biological replication remains a
separate gap. Passing these checks could support a regression diagnostic, not a
new prospective or5% specimen-level validation claim.

## Subsequent transport diagnosis: immediate failure explained

Jobs23130416 and23130429 run on nslab I/O nodes, each1CPU/2GB/2h, with no
sequence decoding. Both complete0:0 in3s and17s, respectively. Their process
resource records and complete submission arguments are retained in
[the resource summary](../dev/dntr_correspondence_20261109/transport_resource_summary.json)
and [jobs.json](../dev/dntr_correspondence_20261109/jobs.json).

The [first diagnostic](../dev/dntr_correspondence_20261109/transfer_diagnostic_v1/report.json)
independently hashes the retained802-byte partial: MD5
`130fb596e3deee5a65b653e4a5e9dce6`, no gzip magic. Its bytes are an Apache HTML
directory index. HEAD and a bounded range GET redirect to a trailing-slash URL
and reproduce that content. This establishes the immediate integrity failure;
it supersedes the earlier checkpoint's unknown response type. Why the archive
path is served as a directory remains unresolved.

The [follow-up](../dev/dntr_correspondence_20261109/transfer_diagnostic_v2/report.json)
uses HEAD for the other five gaps and system curl independently for the first.
Four additional URLs also advertise802-byte HTML directory responses.
SRR12861527 mate2 instead advertises `application/x-gzip` and804,194,663bytes,
matching expected size. That header is not checksum validation or an accepted
repair. Curl independently confirms the first URL's301 redirect and HTML final
response. Its nonfatal libcurl version-information warning is preserved;
exit0 and the recorded response support transport agreement, not a claim of
library independence or biological validity.

The original design and first diagnostic's input/output hashes remain intact.
Follow-up scope is in a separate
[amendment](../dev/dntr_correspondence_20261109/TRANSFER_DIAGNOSTIC_FOLLOWUP.md).
Source_v4 executed the first diagnostic. Source_v5 is an unused intermediate
snapshot; source_v6 records the explicit follow-up design/parameter and executed
the second. A briefly appended design note was moved into the separate amendment
before submission, restoring and checking the first run's exact design hash.
No failed transfer was replaced, no production code changed and no scientific
threshold was relaxed.

Do not repeat the six-file staging unchanged: five current endpoints do not
provide advertised archives. Preserve the incomplete-input status; usable
same-cell pairs must explicitly exclude missing mates. A future bounded truth
investigation can select complete declared pairs prospectively by metadata,
without RNA outcomes. Molecular correspondence and DNA callable sensitivity
remain the next scientific gates; transport diagnostics add no null specimen,
clone truth, power or accuracy claim.
