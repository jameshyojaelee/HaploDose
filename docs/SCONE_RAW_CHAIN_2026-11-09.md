# Raw scONE reads expose library-specific processing uncertainty

Bounded primary-read inspection finds a material mismatch between one legacy
library and the archived processing settings. Independent integer reconstruction
passes. This is a development measurement-chain investigation, not validation
of HaploDose, admission of a cohort or certification of single-clone controls.
Frozen method results, original recovery assertions and error/power criteria
remain unchanged.

## Question and scientific contract

The measured quantities are literal paired-read prefixes and fixed-adapter
positions in the first250000 pairs of metadata-selected libraries. These are
deterministic sequence diagnostics, not random read sampling, molecule counts,
cellular replication or biological specimens. No confidence interval for a
population frequency follows from this sampling rule.

The identifying assumption for interpreting a DNA adapter position is the
source-specified visible5-base tag followed by a candidate5- or8-base UMI and
the full fixed19-base Tn5 ME sequence AGATGTGTATAAGAGACAG. Different chemistry,
adapter trimming, sequencing error and incidental motif occurrences are competing
explanations. The full exact ME at the predicted offsets distinguishes the two
layouts for reads retaining that sequence. It cannot establish the layout of all
reads, RNA UMI length or the commands executed for the published measurements.

Computational correctness is assessed by independent raw-file identity checks,
parsers and counting mechanisms. Statistical identification concerns the retained
sequence layout only. Biological CNA profiles and same-cell clone membership
remain unvalidated. No cells, libraries or seeds increase biological null n=0.

## Legacy whole-cell library: independent checks complete

The locked lexical GSM/run rule selects GSM5609521/SRR16194948, whose two raw
archives total171949409 bytes. ENA MD5 and byte counts match, with SHA256 identities
retained. Literal GEO title HCT116 scONEseq001 conflicts with description HCT116
single cell097. Numeric suffixes therefore cannot safely join processed rows;
literal GSM/SRR identity and an explicit crosswalk are required.

| First250000 read2 records | Exact prefix count |
|---|---:|
| DNA TCATG |125443|
| GEO RNA ATCGT |79673|
| Archived-source RNA AACTT |25|
| Other prefixes |44884|

Within the DNA-prefix stratum, the exact forward ME occurs112521 times at
zero-based offset10, predicted by5 tag+5 UMI bases, versus3 at offset13,
predicted by5 tag+8 UMI bases. All offsets0..40 and the reverse ME are retained;
no favourable offset or frequency cutoff was selected. An independent actual-read
regex calculation agrees with literal slices exactly, and synthetic5/8 controls
pass. This supports the5-base candidate layout for the retained-ME-positive
legacy DNA reads. RNA UMI identity remains unresolved and cannot be borrowed from
the DNA result.

The deposited author archive6796059 passes its declared checksum. Its barcode,
separation, driver, DNA mapping and RNA quantification scripts are byte-identical
to the five inspected current files. The archive retains AACTT as RNA barcode,
an8-base UMI extraction setting, and a downstream trimmed-5RNA filename although
the barcode record is named RNA. These facts expose unresolved source settings;
the deposited source does not prove those commands generated this library.
No external author script was executed or modified.

## One justified later-library contrast

The discrepancy justifies one further metadata-selected library, rather than
automatic study-wide reprocessing. Before opening its reads, NUCLEUS_DESIGN.md
locks GSM6683025/SRR22047451: lexical GSM/run among OTHER titles describing
HCT116 extracted nuclei. Declared paired size295064675 bytes. The same250000-pair
bound, six exact motifs, all offsets and candidate DNA layouts are retained.
AACTT additionally has its own positional stratum in this branch; the legacy
branch and all earlier artifacts are unchanged.

Jobs23036460/23036461/23036462 complete staging, independent reconstruction and
layout inspection with afterok dependencies. Both archive sizes and MD5 match;
the independent parser reproduces every count exactly. The bounded results are:

| First250000 read2 records | Legacy whole cell | Later nucleus |
|---|---:|---:|
| DNA TCATG |125443|130832|
| GEO RNA ATCGT |79673|15|
| Archived-source RNA AACTT |25|13996|
| Other prefixes |44884|105157|
| DNA ME at offset10: candidate5-base UMI |112521|27|
| DNA ME at offset13: candidate8-base UMI |3|69691|

Independent actual-read adapter matching and both synthetic controls also pass
for the nucleus library. This supports different DNA layouts in retained-ME
reads and different visible RNA tags in these two specific libraries. It does
not establish RNA UMI identity or classify the other-prefix reads. Platform,
preparation and library date differ from the legacy case; any observed difference
cannot isolate a single causal factor or generalize to all study libraries.

## Primary primer table recovered; executed-library provenance remains uncertain

The [published supplement](https://pmc.ncbi.nlm.nih.gov/articles/instance/9812385/bin/sciadv.abp8901_sm.pdf)
is the next source for the independent RNA-primer contract. Browser access
returned a checking-browser page. The one prespecified official Europe PMC
supplementaryFiles retrieval job23037137 fails after22s with
IncompleteRead(465492 bytes read). Its failed manifest and exit1 are preserved;
no complete package, PDF or primer table was obtained from that route. This is an
access failure, not a numerical test or evidence that the published primers are
wrong. No retry of that locked endpoint or alternative specimen was selected.

PRIMER_SOURCE_DESIGN.md and primer_source_v1 record the bounded request, archive
member rule and existing conda pdftotext identity. The parser/extraction path
was not reached and is not validated by this failed retrieval.

The subsequent official NCBI OA route23037689 completes3s but returns404, with
source_accessible=false. It used a retired endpoint. The
[current official documentation](https://pmc.ncbi.nlm.nih.gov/tools/pmcaws/)
confirms the legacy service was removed inAugust2026. After verifying its supported
replacement, PRIMER_CLOUD_DESIGN.md locks retrieval of the same supplement from
the metadata-identified public PMC cloud object. Job23038001 succeeds3s, with
the declared MD5ea454b5fbfccc83e77de13f7ce8c5fe2 exactly matching22583764 bytes.
SHA256 is266debbe0bebdeecf9da29ca6cb12d69f5e7b53a4ea7ba4e885b255a4fd65cad.
The source is acknowledged as NIH NLM NCBI PubMed Central Article Datasets,
accessed2026-11-09. Both earlier results and their source snapshots remain.

Table S2 on PDFpage27 specifies the DNA tag TCATG and eight variable UMI bases,
and RNA tag AACTT with eight variable bases followed by short GGGHN/TTTVN or
polyT priming sequences. This independently supports the later candidate design;
it does not explain away the legacy read mismatch or prove which chemistry and
commands were executed. Earlier patent examples describe another candidate
protocol, not exact executed-library provenance.

RNA_LAYOUT_DESIGN.md now locks a direct RNA-flank diagnostic for the two already
selected250000-read2 subsets, both5/8 candidate offsets, all motif positions and
independent exact matching. This tests RNA primers directly rather than transferring
DNA layout. Job23038319 completes9s with all nine synthetic controls and independent
actual-read counts exact. There are no too-short records in the three declared
tag strata for either candidate short-primer pattern.

| Exact RNA-prefix stratum | Reads | Candidate UMI bases | GGGHN | TTTVN | >=20T then VN |
|---|---:|---:|---:|---:|---:|
| Legacy ATCGT |79673|5|38875|10433|1990|
| Legacy ATCGT |79673|8|848|1389|533|
| Later AACTT |13996|5|100|181|11|
| Later AACTT |13996|8|6625|2450|221|

H=ACT,V=ACG,N=ACGT. Each class is counted separately; counts are not added to
claim known molecules. All overlapping triplet offsets0..40, both RNA tags and
DNA diagnostic strata are retained. Genomic low-complexity sequence can match a
short pattern incidentally; the DNA strata are diagnostic contrasts, not certified
negative controls. The source-specified priming structures directly support the
five-base candidate in legacy RNA and eight-base candidate in later RNA. This
falsifies a universal eight-base rule for the source-compatible legacy primer
structure. It does not prove every read's UMI identity, exclude mixed chemistry
or identify exact executed commands. No cutoff or processing rule was deployed.

## Decision and evidence needed next

Do not apply one archived demultiplexing/UMI setting universally, round deposited
RNA values into raw counts, reconstruct truth from RNA groups, or choose trimming
to improve CNA agreement. A justified library-specific count builder requires
library-specific RNA primer/UMI provenance as well as DNA layout, compatible
reference identity, literal cell correspondence and an independently checkable
counting definition. The bounded evidence now identifies supported candidate RNA
and DNA layouts for these two libraries. Exact executed chemistry, mapping and
deduplication, count scale and processed-cell correspondence still need independent
checks before downstream truth or calibration. These gaps are not reasons for
more CNA fits.

The brain raw-DNA correspondence gap described in SCONE_FEASIBILITY remains.
Same-cell DNA at documented depth and resolution is still required for membership
validation; biological replication is still required for specimen error bounds.
Original thinning2/5, three recovery failures and failed false-split/power criteria
remain. No method or scientific default is promoted by these checks.

## Execution and reproducibility

All numerical work runs on compute nodes using isolated validated wheels,
immutable source snapshots and fresh outputs. The linked chain is inputs and
literal library selection -> source/environment record -> output -> independent
audit. Protocols, source identities, raw metadata and reads, complete reports,
recursive hashes, commands and time-v resource records are under
`dev/scone_raw_chain_20261109`; jobs.json identifies each execution.

| Completed job | Result | Runtime | Process peak RSS KiB |
|---|---|---:|---:|
|23034188 legacy staging/prefixes|pass|58s|92104|
|23034467 independent prefix audit|pass|5s|96612|
|23035054 archived-source identity|pass|7s|92744|
|23035665 adapter layout/reference|pass|7s|82948|
|23036460 nucleus staging/prefixes|pass|191s|92940|
|23036461 nucleus independent prefix audit|pass|5s|92612|
|23036462 nucleus adapter layout/reference|pass|5s|83896|
|23037137 primary primer-source access|failed: incomplete HTTP response|22s|84168|
|23037689 retired NCBI OA route|completed; source unavailable HTTP404|3s|83272|
|23038001 current PMC cloud supplement|checksum and extraction pass|3s|108644|
|23038319 direct RNA-primer layout|exact counts and9 controls pass|9s|84780|

These add no production dependency or package changes. Full package suite was
not rerun for diagnostic-only scripts; latest isolated validation remains365
Python+4 Rust passes, with1 GPU-only skip. The separate prior-wheel GPU job
23024875 remains pending; no GPU pass is claimed. Frozen results and historical
decisions remain intact.

Diagnostic scripts added under `benchmarks/hardening`: scone_raw_chain.py,
scone_raw_chain_audit.py, scone_archived_source.py, scone_read_layout.py,
scone_primer_source.py, scone_primer_cloud.py and scone_rna_layout.py. The later
library branches preserve the earlier default branch and executed snapshots.
Current STATUS.md, claims.json and task STATE.md/jobs.json record these decisions;
none of the package or original test sources changed in this diagnostic batch.

Primary input identities are linked through the literal
[legacy GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5609521),
[later GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6683025)
and [deposited author-code archive](https://doi.org/10.5281/zenodo.6796059).
These source records support provenance; the integer findings above come from
the separately hashed raw reads and independent computations.
