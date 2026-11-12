# Auditable scONE read preparation and reference identity

The two locked development libraries now have independently audited candidate
RNA/DNA read layers and a verified reference/tool chain. This is progress toward
independent measurements, not new CNA accuracy or clone validation. Their
biological null-specimen count remains0. Original thinning2/5, three recovery
failures and the failed false-split/power criterion remain unchanged.

## Estimand, assumptions and scope

The initial estimand is the exact read-pair subset and sequence transformation
under an explicitly versioned library layout. The observation is the first250000
paired records per already selected raw archive, not random sampling. The working
measurement model is source-compatible molecular tags and primer structure;
damaged adapters, mixed chemistry and incidental prefix matches remain competing
explanations. Direct RNA/DNA flank evidence is reported in
[the upstream investigation](SCONE_RAW_CHAIN_2026-11-09.md).

Computational correctness requires exact identity, coordinate, sequence/quality
and count conservation. Statistical identification of RNA molecules and genomic
alterations does not follow from those checks. Biological validity requires
independent DNA and supported cell correspondence, at demonstrated resolution.
Published labels and source-script agreement do not supply that evidence.

## Preparation: every emitted record independently reconstructed

Schema1 names each literal GSM/SRR and candidate layout: legacy RNA ATCGT,
DNA TCATG,5 candidate UMI bases; later RNA AACTT,DNA TCATG,8 candidate bases.
Both mates receive the same original pair ID and candidate UMI tag. Original R1
sequence/quality remain intact. R2 removes the exact tag and UMI. DNA additionally
requires/removes the full fixed19-base ME; it does not blindly discard20 bases.
RNA retains its RT primer sequence for the declared alignment treatment.

No fuzzy rescue, quality or convenient read-length cutoff, duplicate removal or
selection on gene/CNA outcomes enters preparation. All exact RNA-tag pairs with
valid candidate UMI and nonempty insert remain, without selecting the promising
primer subset. Deterministic FASTQs retain order and exact sequence/quality cuts.

| Library | Fixed input pairs | Candidate DNA | Candidate RNA | DNA tag without full fixed ME | Other prefix |
|---|---:|---:|---:|---:|---:|
|GSM5609521/SRR16194948|250000|112521|79673|12922|44884|
|GSM6683025/SRR22047451|250000|69691|13996|61141|105172|

No invalid candidate UMIs or empty post-trim inserts occur in these subsets.
Other prefix includes the alternative RNA tag; it is unresolved, not certified
contamination or biological absence. Exclusion is explicit, and these candidate
layers are not complete counts of true RNA/DNA molecules.

Audit23040016 reconstructs every emitted mate directly from raw reads using a
separate text parser and regex rather than producer byte slicing. Every ID,
candidate UMI, sequence, quality, order, output hash and exclusion count agrees
exactly. Each library has250000 distinct input pair IDs. No output leftovers are
accepted. These integer/identity checks use no numerical tolerance.

## Reference and tools: actual identities rather than filenames

Reference23040441 verifies all194 FASTA contigs against the retained parent
10X GRCh38-2020-A sequences, after only documented primary-contig name changes.
The reference is this derived FASTA, not an independently downloaded Ensembl
assembly. Byte hashes match its original manifest; every contig length also
matches FAI and STAR. All15 STAR index artifacts and reference/annotation build
hashes agree. Ensembl112 GTF records have valid one-based closed bounds, all
canonical chromosomes have exon support, and no records refer to absent contigs
or extend past present contigs. The protected setup lock is shared throughout.

Two earlier inventory failures remain preserved:23039565 treated nonzero ldd as
failure,23040177 assumed STAR was ELF. Inspection established that generic STAR
is a CPU dispatcher shell script. The successful gate explicitly inventories
the concrete dynamic STAR-plain native executable, its ELF program headers and
loaded libraries; samtools and featureCounts receive the same concrete checks.
Those failures were in the new diagnostic, not evidence of an invalid assembly.

Tool staging23040739 copies the exact27 installed conda packages into a fresh
offline environment, with no dependency updates or source-environment changes.
All selected package names/versions/builds/URLs/MD5 match;842 source conda-meta
identities remain stable. All staged regular files are hashed and escaping
symlinks rejected. Concrete versions are STAR2.7.11b,samtools1.22.1/htslib1.22.1,
featureCounts2.1.1. Conda relocation can change binary bytes; actual staged hashes
are recorded rather than claiming source-byte identity.

## Fixed RNA mapping/counting gate

RNA_MAP_DESIGN.md fixes source RNA STAR settings, originalR2-first order and
retained RT priming sequence. It explicitly uses STAR-plain and no shared-memory
index. Native STAR gene counts and featureCounts paired unstranded exonic
gene_id fragment counts will both be retained. The latter requires both -p and
--countReadPairs, consistent with the
[featureCounts primary documentation](https://subread.sourceforge.net/featureCounts.html).

These are PCR-affected integer read-fragment measurements. They are not UMI
molecule counts, calibrated RNA uncertainty or independent DNA truth. Priming,
PCR, intronic RNA, strandedness and alignment ambiguity remain consequential
measurement assumptions. No DNA reads are mapped with an RNA pipeline.

Mapping23041602 and afterok audit23041603 pass. The independent SAM parser
reconciles every primary BAM mate to the prepared QNAME/UMI:159346 primary records
for79673 legacy templates and27992 for13996 nucleus templates, exactly one first
and second mate per original candidate UMI-tagged read ID. No unexpected IDs or
repeated primary mates occur. Count table integers independently sum to Assigned.
The complete summary family, native STAR counts and flagstat are retained.

| Fixed RNA subset | Input pairs | STAR uniquely mapped templates | Exonic fragments assigned | featureCounts NoFeatures records |
|---|---:|---:|---:|---:|
| Legacy |79673|55784|22678|31467|
| Later nucleus |13996|11512|1473|9955|

Columns have deliberately explicit units. In particular, featureCounts summary
categories can include separate multialignment records; their totals are not
unique input-template totals. The audit establishes primary-pair conservation
directly from BAM records, not by summing those summary categories. Counting units
and independent output-sum checks are correct, but every genomic alignment or
gene assignment is not yet independently validated.

The later subset's substantial unique mapping with few exonic assignments
identifies the next consequential assumption. Intronic/nascent signal, intergenic
signal, residual genomic material and counting/annotation behavior are alternatives.
Low exonic counts alone do not distinguish them or establish low RNA content.
The smallest next analysis is a locked bounded panel of primary unique
NoFeatures templates, with independent CIGAR/block and gene-body/exon overlap
inspection. Any gene-body counting path needs a separately declared estimand and
matched support; more assigned fragments alone cannot justify its promotion.

No HaploDose fit, new calibration,
prospective cohort admission, spatial assignment or DNA truth is produced here.

## Execution and reproducible records

All executable checks/builds/mapping run on compute nodes. Scripts, protocols,
immutable snapshots, actual commands and environment inventories, exact source
and output identities and process/scheduler resource records reside under
`dev/scone_count_chain_20261109`. jobs.json preserves failures and dependencies.

| Job | Status | Seconds | Process peak RSS KiB |
|---|---|---:|---:|
|23039565 reference gate v1|failed on wrapper ldd|104|85256|
|23040177 reference gate v2|failed on wrapper ELF assumption|92|85248|
|23040441 reference gate v3|pass|97|85072|
|23040015 bounded preparation|pass|45|130700|
|23040016 independent preparation audit|pass|6|134228|
|23040739 isolated tool environment|pass|17|217912|
|23041602 RNA mapping/counting|pass|80|30296600|
|23041603 primary-pair/count-unit audit|pass|3|resource record linked in task|

New diagnostic scripts are scone_reference_gate.py,scone_prepare_reads.py,
scone_prepare_reads_audit.py,scone_tool_environment.py,scone_rna_map.py and
scone_rna_map_audit.py under `benchmarks/hardening`. The package and original
tests are unchanged by this batch; full suite not rerun for these diagnostic-only
changes. Latest package validation remains365 Python+4 Rust pass,1GPU-only skip.
The prior-wheel GPU23024875 remains pending; no GPU or superiority pass claimed.
Frozen results, historical decisions and scientific criteria remain unchanged.
