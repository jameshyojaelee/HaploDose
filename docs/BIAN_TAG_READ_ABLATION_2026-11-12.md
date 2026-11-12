# Tag validity and read evidence need different checks

## Completed result, 2026-11-12

Job **23232643 completed 0:0**, in 3m21s on nslab CPU, two CPUs and
48 GB requested. All **201** source identities agree before and after execution.
The original fixtures, every primary read identity, the complete private-index
hashes and native featureCounts agreement across all 54,849 annotated genes
pass. Annotated zero counts remain present. STAR 2.7.11b, featureCounts 2.1.1
and the exact commands appear in the saved report.

Of the 9,000 rejected-tag records, **46 reads across 43 genes** are uniquely
assigned, including **28 CIGAR-N reads across 27 genes**. The other primary
statuses are 8,922 unmapped, 25 without features, five ambiguous and two
nonunique. Two additional secondary/supplementary records are auxiliary
alignments, not additional source reads. Native nonunique totals include these
auxiliary records; the primary-ID accounting remains exact.

The first three lexical assigned splice examples include their bounded original
R1/R2 sequence and alignment fields in `BOUNDED_SPLICE_EXAMPLES.json`. Unspliced
exon evidence also permits DNA contamination; splice-supported alignments
justify the narrower read/tag distinction but do not prove molecular validity
or genomic specificity. STAR MAPQ 255 is its alignment field, not a Phred
confidence claim; uniqueness here uses NH=1.

The measured maximum resident set from `/usr/bin/time -v` is **29,066,876 KiB**.
Slurm's sampled 49,476 KiB missed the brief native-index allocation and must not
be substituted for that maximum. All rejected tags stay excluded from molecular
counts. No CNA, clone, uncertainty or biological error-control improvement is
established. The separately versioned export tests this read conservation
correction without changing the frozen guard or original results.

[Completed report](../dev/bian_tag_read_ablation_20261112/run_23232643/results/REPORT.json),
[bounded read examples](../dev/bian_tag_read_ablation_20261112/run_23232643/results/BOUNDED_SPLICE_EXAMPLES.json),
[completion record](../dev/bian_tag_read_ablation_20261112/completion_checkpoint.json),
[versioned correction](BIAN_READ_TAG_EXPORT_2026-11-12.md).

The following design and initial pending checkpoint retain their original context.

Implemented a matched ablation to test whether the new candidate tag-layout
guard also removes usable RNA read evidence. Invalid tags must remain excluded
from molecular counts. Their R1 sequence may still be informative, so that
second claim needs a direct measurement before any correction is promoted.

## Frozen experiment

Use exactly9000 records from the first100000 SRR7461687 pairs that frozenv1
called eligible and frozenv2 rejected for complete source-primer overlap in R2.
Keep the same v1 R1 trimming. Neither candidate is retuned and the original
guard remains unchanged. All source/raw identities are required. No selection
by gene,DNA outcome,clone identity or effect size.

After successful RNA23232347, reuse its private STAR index and exon annotation,
tool binaries and mapping/assignment settings. Only the recorded thread count
changes from8 to2. Index files are completely hashed before and after. Native
STAR single-end/NH/unmapped and unstranded featureCounts processing match the
parent path. The existing index is reused, not rebuilt or modified.

The new independent read parser has hand-counted spliced/unspliced positives,
nonunique/ambiguous/unmapped negatives and duplicate-primary rejection. Every
selected primary ID must appear exactly once. Every annotated gene, including
zeros, must match native featureCounts. No exact-tag/molecule counts are
exported from this rejected stratum. Preserve all failures and status categories.

## Distinguishing observations and falsification

If the rejected material is primarily primer/adapter artifact, it may produce
little unique exon evidence. Genuine transcript fragments with unusable tags
can still align and be assigned. Unspliced exon evidence also fits DNA
contamination; it cannot by itself resolve RNA origin. Report assigned CIGAR-N
reads separately and retain up to the first3 lexical uniquely assigned spliced
records with original bounded R1/R2 and alignment fields. These are supporting
alignment observations, not proof of independent molecules or biological truth.

Zero gene assignments rejects the proposed gene-read information-loss mechanism
for this stratum. Unspliced-only assignments leave the main alternative open.
Splice-supported assignments justify testing separate read eligibility in a
subsequent isolated correction. They do not justify admitting invalid tags,
choosing a convenient mapping fraction or claiming improved CNA inference.
All outcomes retain the full family of annotated genes and rejection reasons.

## Execution and scientific scope

Job **23232643**: account/QOS nslab,CPU partition,two CPUs,48GB,one hour,
dependency `afterok:23232347`,201 identities. STAR index loading requires tens
of GB; the9000-read analysis itself is bounded. No array,GPU,interactive QOS,
new dependency,shared installation or alteration of a running job. Static AST
and launcher syntax pass. Initially pending Dependency; no executable fixture,
gene-count result or scientific decision is claimed yet.

One technical cell,biological inference n=0. Bian remains secondary/descriptive.
No DNA profile is opened or generated,no genomic/clone score is calculated,and
no uncertainty or error guarantee follows. Original recovery,false-split/power,
uncertainty,comparator and genomic-admission criteria stay unchanged. All434
prior claim entries are preserved; the separate entry records this design.
Public source and frozen results remain untouched.

[Frozen design](../dev/bian_tag_read_ablation_20261112/DESIGN.md),
[candidate parser and native workflow](../dev/bian_tag_read_ablation_20261112/run.py),
[resource header](../dev/bian_tag_read_ablation_20261112/run.sbatch),
[submission/dependency](../dev/bian_tag_read_ablation_20261112/submission.json),
[201 identities](../dev/bian_tag_read_ablation_20261112/SNAPSHOT.json),
[original preprocessing guard and evidence gap](BIAN_PREPROCESS_2026-11-12.md),
[prerequisite native RNA count pilot](BIAN_RNA_COUNTS_2026-11-12.md).
