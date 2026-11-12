# Bian RNA preprocessing: source-primer contamination and unresolved tags

Implemented a candidate source-specific RNA preprocessor and a stricter version
that rejects complete known-primer overlaps with the intended UMI/RT region.
All executable checks pass on compute nodes. The raw reads nevertheless show
unresolved tag structure,so no tagged specimen FASTQ,molecule counts,gene matrix,
alignment or genomic outcome has been constructed. This is an upstream correctness
repair in an isolated measurement path,not a demonstrated HaploDose performance
improvement. Frozen package/results and all original scientific criteria remain.

## Scientific object and why these checks were needed

The technical object is eligibility of a read under the proposed capture layout:
read2 carries an eight-base cell barcode,eight candidate tag bases and an RT
polyT tract. The [cited primary RNA protocol](https://link.springer.com/article/10.1186/s13059-018-1416-2)
and its preserved primer table support the barcode and primer sequences. They do
not prove the exact processing/design of this deposited human library. Uniform
random eight-base tags,independent molecules and a calibrated count likelihood
are separate assumptions. Neither a matching barcode nor a successful parser
establishes them.

The competing explanations for apparent tag bias include frequency dominance,
primer-derived sequence mistaken for a tag,and source-specific synthesis/layout/
processing. Distinct-tag weighting,full known-primer matches and the expected
RT tract can distinguish parts of these explanations. They cannot determine
true gene-wise molecules or DNA dosage without further measurements.

The specimen/run remains the original metadata-selected CRC04_LN1_349 /
GSM3241543 / SRR7461687. No genomic outcome selected it. There is one technical
cell and no biological replication for an inferential claim. Former evaluations
stay regression-only; Bian remains secondary/descriptive.

## 1. Literal preprocessing and independent reconstruction

Version1 recognizes only exact anchored source TSO/PCR sequences,longest match
first,and slices qualities with bases. Unmatched residues remain. The candidate
terminal-A rule uses at least25 bases,explicitly derived from the assay's RT
tract length; it is not a recovered historical parameter. A biological A25 tract
cannot be distinguished from a technical tail by that sequence rule. A matched
untrimmed ablation remains required before adopting it in measurement analysis.
The post-trim N exclusion is the source protocol's exact integer rule
`10*N > retained length`. Homopolymer candidate tags remain valid; foreign
barcodes,invalid tags,empty reads and ambiguous identities are explicit outcomes.

Literal fixtures test primers,internal-primer retention,tail24/25,sequence-quality
alignment,identities,homopolymer tags,empty/N-rich reads and exact N boundaries.
An independent regex/index implementation reproduces all100000 candidate records
in the first deposited block,including output digests. Both full raw input
SHA256s match. There are99665 eligible candidates and335 empty reads after
candidate trimming. Exact TSO/PCR/no-primer counts are40574/6125/53301.
This agreement establishes implementation consistency,not correct biological
interpretation of the proposed tags or removed tails.

## 2. Tag bias persists beyond one read-weighted observation

The first block revealed C-poor candidate positions. The follow-up locked three
disjoint100000-pair blocks at starts0,871083,1742166,chosen from the already known
1842166-record library size. These are records in deposited order,not independent
biological specimens or random validation samples. Full raw SHA256s match.

| Block start0 | Distinct candidate8mers | Exact read1-sequence/tag pairs | Dominant candidate tag;reads |
|---|---:|---:|---|
|0|7845|97892|TGTACTCT;9262|
|871083|7799|97899|TGTACTCT;9262|
|1742166|7832|97927|TGTACTCT;9210|

Giving each observed distinct tag unit weight does not remove the C-poor pattern:
its per-position C fractions remain roughly1.7–4.2% across these panels. The
read-weighted pattern has two positions near10% C and six much lower positions.
All integer totals are independently reconstructed from the complete emitted
frequency table. A99-to1 two-tag fixture checks the different weighting estimands.

Distinct8mers are not molecules,and distinct read1/tag pairs do not establish
molecular independence. The former collides across genes; the latter can vary
with fragmentation or recur for different biological molecules. The experiment
does not prove PCR is the sole cause or identify the intended tag alphabet.
No rare C-containing tag is discarded and no4^8 capacity is assumed.

## 3. Complete source-primer overlap falsifies part of the candidate layout

The dominant8mer matches part of the source primer's reverse complement.
Version2 checks the complete23-base source PCR reverse complement inside the
region intended to hold UMI/RT sequence. Exact starts8..40 are inconsistent with
that proposed capture layout; a later primer in the RNA insert is retained.
This boundary comes from barcode8 + UMI8 + RT tract25,not an optimized cutoff.
Every decision agrees with independent exhaustive start-position comparison.

The old literal fixtures pass unchanged under version2. A new counterexample
is eligible under version1 but rejected under version2; a later-primer positive
control remains eligible. All matched real early-primer reads were formerly
eligible and are now rejected:

| Block start0 | Exact T25 structure | Early complete source primer | Other/unresolved layout | Newly rejected by guard |
|---|---:|---:|---:|---:|
|0|69811|9000|21189|9000|
|871083|70121|8976|20903|8976|
|1742166|69949|8932|21119|8932|

The guard rejects26908 formerly eligible reads across300000 inspected pairs.
Most complete-primer matches start at read2 position11. Bounded examples retain
exact read IDs and read2 prefixes,including repeated primer sequence. The
observed sequence contradicts the proposed tag/RT interpretation; the chemical
origin of those reads is not established. Some might also fail later alignment,
so these rejections are not a measured improvement in gene counts or CNA accuracy.

The exact-T25 subset still has very few Cs in its proposed tags under both
read and distinct-tag weighting. That leaves intended design,selection or
processing unresolved. Exact T25 is a descriptive structural subset,not a new
whole-specimen success criterion. Version2 also does not certify the remaining
272104 eligible candidates as independent molecules. All986 empty-read and two
invalid-tag exclusions remain. No downstream gene count or genomic outcome was
used to choose this guard.

## Compute,evidence and next decision

| Job | Purpose | nslab CPU resources | Source identities before/after | Result;process peak |
|---|---|---|---:|---|
|23231829|Version1/literal and independent100000-record check|1CPU2GB1h|150|COMPLETED0:0,7s;45900KiB|
|23231893|Three-block tag discrepancy|1CPU2GB1h|157|COMPLETED0:0,21s;116028KiB|
|23231957|Version2 structural guard and unchanged/new fixtures|1CPU2GB1h|166|COMPLETED0:0,22s;43972KiB|

All jobs use the existing isolated runtime,neutral directories and one thread.
Scheduler RSS was not usefully sampled. No arrays,new dependencies,shared writes,
public package edits or frozen-artifact changes. No full Python/Rust suite was
rerun for this isolated read-QC work; the relevant executable checks are above.

The run-specific original primer specification/read-processing source remains
an evidence gap. The official GEO text request was inaccessible; the original
Science supplement was not retrieved. Later modified scTrio protocols are not
substituted for this run's source. The original publication identity is recorded
by [PubMed](https://pubmed.ncbi.nlm.nih.gov/30498128/),which does not itself supply
the missing processing details. Additional simulation cannot resolve that gap.

Retain version2 as a tested isolated input guard; do not claim full tag extraction,
uniform capacity,calibrated molecule counts or genomic validation. Keep the
complete candidate/reference/comparator/scorer measurement freeze before gene/
DNA outcome analysis. The three production recovery failures,whole-specimen5%
false-split/power tradeoff,uncertainty and all-competitor superiority remain
unresolved. No scientific criterion has been relaxed.

[Version1 design](../dev/bian_preprocess_20261112/DESIGN.md),
[declared parameters](../dev/bian_preprocess_20261112/parameters.json),
[version1 code](../dev/bian_preprocess_20261112/preprocess.py),
[independent implementation](../dev/bian_preprocess_20261112/reference.py),
[version1 result](../dev/bian_preprocess_20261112/run_23231829/results/REPORT.json),
[tag-bias design](../dev/bian_tag_bias_20261112/DESIGN.md),
[tag-bias result](../dev/bian_tag_bias_20261112/run_23231893/results/REPORT.json),
[complete tag-frequency table](../dev/bian_tag_bias_20261112/run_23231893/results/TAG_FREQUENCIES.tsv),
[strict guard design](../dev/bian_layout_guard_20261112/DESIGN.md),
[version2 implementation](../dev/bian_layout_guard_20261112/preprocess_v2.py),
[guard result and bounded read examples](../dev/bian_layout_guard_20261112/run_23231957/results/REPORT.json),
[unchanged and new fixtures](../dev/bian_layout_guard_20261112/run_23231957/results/FIXTURES.json).
