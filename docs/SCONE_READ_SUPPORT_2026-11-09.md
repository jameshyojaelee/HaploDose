# Exon assignment misses gene-body signal, but more counts do not establish accuracy

The raw-read investigation now distinguishes annotation support from mapping
failure in two locked scONE development subsets. It supports a bounded comparison
of count definitions; it does not identify every intronic fragment as RNA, select
a new package default or validate genomic CNA/clones. Frozen results and the
original three recovery and false-split/power failures remain unchanged.

## Scientific question and decisive observations

The estimand is observed alignment support for annotated exons and gene bodies,
conditional on the same assembly, annotation, paired reads and mapping. The
working observation is a paired RNA-tagged fragment, with explicit candidate UMI
metadata. Gene-body-only support is compatible with nascent/intronic RNA, but
residual genomic material and incomplete annotation remain competing explanations.
NoFeatures with observed exon support would instead contradict a simple intronic
loss explanation and trigger investigation of counting/coordinates.

Computational correctness concerns exact identities, coordinates, support sets
and integer counts. Statistical identification of RNA origin and DNA dosage
remains separate; alignment/gene overlap alone cannot settle it. Biological
validation still needs suitable independent DNA measurements and replication.
Cells/read panels do not increase certified biological null n=0.

## Fixed read panel, with successful and ambiguous controls

Before opening detailed assignment results, NOFEATURES_DESIGN.md fixes lexical
first50 NoFeatures,10 Assigned and10 Ambiguity templates per library. Eligibility
requires exactly one detailed CORE record and two primary mapped NH=1 mates.
No alternative reads or libraries are selected. This purposive panel investigates
mechanisms; it estimates no population error/frequency or biological power.

| Library | NoFeatures: gene-body only | NoFeatures: no annotated gene body | NoFeatures: exon support | Assigned controls with reported exon gene | Ambiguity controls with multiple exon genes |
|---|---:|---:|---:|---:|---:|
|GSM5609521/SRR16194948|36/50|14/50|0/50|10/10|10/10|
|GSM6683025/SRR22047451|39/50|11/50|0/50|10/10|10/10|

Detailed featureCounts output reproduces the original gene counts and full summary
exactly. Seven boundary controls and a spliced/indel CIGAR fixture pass. The first
calculation uses zero-based half-open observed match blocks and an interval index;
brute-force annotation checks agree for the first four selected templates per
library. Independent audit23044269 reconstructs eligibility and all140 selected
identities from CORE/BAM, uses a separate character-wise CIGAR parser and observed
one-based positions, and scans original closed GTF intervals directly. All observed
positions, source sequence/quality, gene/exon sets, classifications and controls
agree exactly. No numerical tolerance determines these integer/coordinate results.

Additional audit23045176 closes the explicit prealignment sequence link: undoing
SAM reverse-strand orientation and STAR's originalR2-first mate ordering reproduces
all280 prepared FASTQ mate sequences and qualities exactly, including soft-clipped
bases. The prepared records already have independent raw-archive slice/quality
conservation. This links observed bases through the executed mapping to the
reported loci without claiming that each genomic alignment is uniquely correct.
The earlier BAM/annotation audit is retained alongside this stronger check.

The panel contains both gene-body-only and unsupported-by-current-gene-body
fragments. It therefore does not support a universal intronic explanation, an
annotation-failure explanation or a claim that these are all genuine RNA.

## Three fixed count definitions on the same mapped templates

The verified gene-body-only observations justify COUNT_VARIANTS_DESIGN.md, locked
before these new counts. The family is exactly exon-only, direct gene-body and
exon-priority rescue. featureCounts changes only the feature type from exon to
gene, retaining the same reads, annotation, unstranded paired-fragment settings
and multimapping/ambiguity rules. The
[primary featureCounts documentation](https://subread.sourceforge.net/featureCounts.html)
describes RNA feature summarization; this is a changed RNA observation definition,
not a cross-assay variance model or inferred genomic quantity.

| Definition or consequence | Legacy | Nucleus |
|---|---:|---:|
| Original exon assignments |22678|1473|
| Direct gene-body assignments |38243|7109|
| Exon-priority assignments |40596|7350|
| Explicit NoFeatures rescues |17918|5877|
| Exon-assigned templates not Assigned by direct gene-body mode |2378|243|
| Assigned targets changed between direct exon/body modes |88|23|

All63140 gene IDs are retained in the matched count table. Native direct counts
equal independent per-gene tallies of unique Assigned CORE QNAME/target pairs.
Exon-priority preserves every original exon assignment and rescues only a unique
Assigned gene-body target with a single original NoFeatures record. Ambiguous,
multimapped, unmapped and other categories are not rescued. Rescued identities
are retained; no template can be counted twice. Every increase reconciles exactly
with the rescued-template list and no baseline gene count decreases.

This demonstrates a tradeoff in measurement definitions: simply switching to
gene bodies increases yield while losing or changing some existing exon
assignments. Exon-priority preserves those assignments computationally, but that
does not prove their biological target or prove the new fragments' RNA origin.
No method is promoted by increased counts. These remain PCR-affected read
fragments, not UMI-deduplicated molecules or calibrated RNA abundance.

## Independent DNA chain and remaining evidence

The next independent measurement route uses the already audited DNA-tagged,
fixed-ME candidate fragments from the same literal libraries. A native BWA
reference-index job23044283 is running on a private copy of the verified FASTA.
It does not hold a shared reference lock during the long build. Its first launch
23043961 failed before analysis because the new shell wrapper attempted direct
execution of a nonexecutable launcher; explicit Bash invocation corrects that
and a clean module stack isolates BWA. The failed log and source remain.

DNA mapping23044784 and independent read/reference audit23044785 are queued with
afterok dependencies. Neither is a result yet. They will retain MAPQ/coverage
support without selecting a cutoff, deduplicating, calling DNA CN or grouping
cells. Suitable DNA sensitivity, processing stability and reference uncertainty
must precede any claim that these two cells supply profile or clone truth.
Gene-body overlap plus DNA alignment cannot by itself certify RNA origin or
solve missing biological replication. No new cohort is admitted here.

## Reproducible chain and compute

The linked chain is frozen RNA BAM/GTF/tool identities -> executed immutable code
and environment -> detailed CORE/read panels/count definitions -> independent
one-based annotation audit. All protocols, selected IDs, sequence/quality/CIGAR,
literal gene sets, rescued templates, commands, native identities, manifests and
resource records are under `dev/scone_count_chain_20261109`.

| Completed job | Checks | Seconds | Process peak RSS KiB |
|---|---|---:|---:|
|23043091|fixed read panel and original-count identity|34|847300|
|23044269|independent reconstruction of all140 templates|17|318800|
|23044457|matched count definitions and rescue conservation|14|253120|
|23045176|all280 prepared FASTQ mate sequence/quality links|20|resource record linked in task|

Added scripts under `benchmarks/hardening`: scone_nofeatures.py,
scone_nofeatures_audit.py,scone_count_variants.py,scone_dna_index.py,
scone_dna_map.py,scone_dna_map_audit.py; native module launcher is
`scripts/scone_bwa_readonly.sbatch`. Original package/tests and frozen artifacts
are unchanged by this batch. Full suite not rerun for these diagnostic-only
changes; latest package validation remains365 Python+4 Rust pass,1GPU-only skip.
Prior-wheel GPU23024875 remains pendingPriority, not passed. jobs.json records
all live dependencies, failures and resources. No SOTA or genomic calibration claim.

## Completed DNA mapping, identity audit and support census

Index23044283 completes in57m15s, retaining the private FASTA, all five BWA
index artifacts, exact native/library identities and matching194-contig FAI.
Its FASTA SHA256 remains
`c974c719f1a50f1d7fa0885c05428a98e2889a92d7dd6b061b6bf9f4a15de1c7`.
The module's actual version string is0.7.18-r1243-dirty; its binary hash is
recorded rather than assuming a pristine upstream build. Mapping23044784 and
independent audit23044785 then complete in43s and9s. No dependency failure is
ignored. All primary mates reproduce the prepared original sequence and quality
exactly after strand reversal; both assemblies match all194 reference lengths.

| DNA support | Legacy SRR16194948 | Nucleus SRR22047451 |
|---|---:|---:|
| Retained primary templates |112521|69691|
| Mapped primary mates |218406|138417|
| Unmapped primary mates |6636|965|
| Proper, both mapped, same-chromosome templates |106289|67802|
| Different-chromosome templates |344|554|
| Both mapped, same chromosome, nonproper templates |1410|832|
| Templates with at least one unmapped mate |4478|503|
| Proper templates with both MAPQ60 |83434|56507|

The full read/pair MAPQ histograms are retained without a cutoff. These are
aligner diagnostics, not calibrated mapping probabilities. The paired-library
identities link RNA and DNA within each literal accession; the earlier legacy
cell001/cell097 metadata conflict still prevents an assumed join to processed
cell labels. Matching read identities does not resolve that conflict.

The next frozen protocol,DNA_SUPPORT_DESIGN.md, asks whether raw fragment depth
can be treated as distinct-molecule support. Census23050103 and separate
SAM-parser audit23050115 complete in12s and9s. Every template is exported,
including nonproper/unmapped templates excluded from the fingerprint census.
Candidate UMIs are allACGT, with lengths5 and8 exactly in their respective
libraries. Fingerprints use ordered primary-mate observed positions, orientation
and CIGAR, with or without candidate UMI. They do not infer unclipped molecular
endpoints or mark PCR duplicates.

| Proper-template census | Legacy | Nucleus |
|---|---:|---:|
| Raw templates |106289|67802|
| Coordinate fingerprints |104043|67318|
| Coordinate-plus-UMI fingerprints |104376|67328|
| Coordinate-plus-UMI multiplicity1 fingerprints |102489|66858|
| Multiplicity2 fingerprints |1862|466|
| Multiplicity3 fingerprints |24|4|
| Multiplicity4 fingerprints |1|0|

There are repeated exact records, but these counts do not identify every repeat
as PCR or certify each unique fingerprint as an independent molecule. Candidate
UMI collisions and alignment ambiguity remain possible. Repeated coordinates
with different UMIs also occur. Do not silently replace read depth with either
fingerprint count, or infer biological independence from low repeat multiplicity.

Fixed500kb,1Mb and5Mb tables include all194 contigs, empty bins and final partial
bins. Leftmost observed template positions usefloor((SAM_POS-1)/width). All three
resolutions conserve exactly the proper-template total. The separate audit
reconstructs every QNAME/fingerprint/UMI, all10435 bin rows per library, all61
MAPQ strata per library and all multiplicity histograms from the original BAMs.
Integer equality is the criterion; no floating tolerance or scientific cutoff
was selected. Stratum-specific unique counts are not summed as global molecule n.

| Completed job | Seconds | Process peak RSS KiB | Scheduler sampled peak KiB |
|---|---:|---:|---:|
|23044283 index|3435|4543688|4602036|
|23044784 mapping|43|5793084|5845864|
|23044785 mate/reference audit|9|271180|not captured by sampler|
|23050103 support census|12|284160|not captured by sampler|
|23050115 independent census audit|9|364812|not captured by sampler|

Outputs are`dna_index_v2`,`dna_map_v1`,`dna_map_audit_v1`,`dna_support_v1` and
`dna_support_audit_v1` under the same task; each retains its source snapshot,
execution environment, commands, native identities, inputs and output hashes.
The two new diagnostic scripts are`scone_dna_support.py` and
`scone_dna_support_audit.py`. Existing source/tests and historical outputs remain
unchanged. Full package suite not rerun for these diagnostic-only additions.

The completed chain establishes record/reference integrity and conditional
support in two deterministic read prefixes. It does not establish calibrated
DNA dosage, callable sensitivity, ploidy, independent clones, random sampling or
biological replication. No GC/mappability correction, smoothing, segmentation,
resolution selection, UMI duplicate rule or CNA reference is promoted. Additional
DNA truth requires justified coverage/bias/sensitivity and specimen correspondence;
more compute cannot turn these two library prefixes into certified biological
null specimens. Biologicalnulln0, original three recovery failures, failed
false-split/power requirements and all scientific criteria remain unchanged.
