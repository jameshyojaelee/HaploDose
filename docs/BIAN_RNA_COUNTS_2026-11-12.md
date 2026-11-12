# Corresponding RNA measurement pilot

Implemented a bounded real-input RNA alignment and counting path using the
qualified preprocessing and native counting components. It exports assigned-read
counts separately from distinct exact gene/tag counts. The latter are candidate
tag measurements, not validated molecular counts: the library-specific tag
chemistry, collisions, sequencing errors and molecule independence remain
unresolved. No HaploDose accuracy improvement is claimed from this preparation.

## Fixed measurement and competing interpretations

Use CRC04_LN1_349/GSM3241543/SRR7461687, first100000 deposited pairs, exactly as
in the earlier preprocessing check. Rehash both complete original FASTQs and
the qualified native hg19 FASTA/Ensembl GRCh37r75 GTF before use. The whole
93-contig genome remains available for competitive mapping; only the previously
qualified autosomal aliases receive exon annotation. Explicit coordinate/bounds
checks prevent silent assembly or feature reassignment.

Frozen preprocessingv2 removes exact known primers/candidate terminal A tails
and rejects the demonstrated early source-primer overlap. It retains homopolymer
tags and makes no uniform tag-complexity assumption. These restrictions are
recorded measurement choices, not original-Bian protocol validation. Attempts
to access the original sample page and Science supplement through the web tool
remain unsuccessful in this continuation; no new protocol source was obtained.
The existing source request remains an evidence gap.

STAR2.7.11b is used as a splice-aware RNA aligner, featureCounts2.1.1 for exon/gene
assignment and samtools1.22.1 for native SAM/BAM inspection. These tools are
appropriate to RNA reads. The private STAR index uses the exact qualified FASTA
and autosomal annotation, sjdbOverhang149 from original150-base reads and no
adaptive second pass. Gene assignment is unstranded, uniquely mapped and
unambiguous. Unstranded counting is a declared working choice; chemistry is not
inferred from an API accepting the reads.

Real alignable transcripts should produce uniquely assigned exonic reads.
Primer artifacts, nonunique sequence and overlapping-gene ambiguity remain
explicit rejections. Distinct gene/tag counts can fall because of true PCR
repetition, tag collisions or biased tag chemistry; counting alone cannot
distinguish these explanations. No Poisson/NB uncertainty or genomic validity
claim follows from integer-valued outputs.

## Decisive computational checks

Before real processing, reuse the unchanged literal preprocessing/counting
fixtures and early/late-primer counterexamples. After native processing,
independently reconstruct counts with tuple sets and a separate SAM parser.
Every eligible primary read ID must be present exactly once. Assigned-read
counts must equal native featureCounts for every annotated gene; exact gene/tag
counts must equal the frozen counter. Retain all annotated zero-count genes.
No minimum mapping/assignment fraction is chosen, and empty assignment is no
evidence. Failure stops execution and retains the failed artifacts.

## Execution checkpoint

Job **23232347**, account/QOS nslab, CPU partition, eight CPUs,96GB,four hours,
constraint `v6&amd`,193 source identities. It is authoritatively **PENDING
(Resources)** at the initial checkpoint, with no start-time estimate and no
downstream job. No fixture, index build or actual RNA alignment/counting result
is claimed yet. AST and shell syntax checks pass; source431 prior ledger values
remain unchanged after appending the separate technical design entry.

The initial vendor-specific CPU constraint unnecessarily limited scheduling:
this task does not invoke JAX/HaploDose inference. After verifying the same job
was still pending, `scontrol update JobId=23232347 Features=` removed that
constraint. Account/QOS,CPUs,memory,time,code and scientific checks are unchanged;
the original frozen launcher is preserved. Only this pending job was amended.
The exact scheduler records before/after are retained in the
[resource amendment](../dev/bian_rna_counts_20261112/pending_resource_amendment.json).

Running checkpoint:the same job starts on ne1dc4-017 after that amendment.
All193 initial source identities pass. Unchanged preprocessing/count fixtures
and early/late-primer counterexamples pass;complete FASTQ/FASTA/GTF hashes pass.
Real preprocessing reproduces the earlier locked first-block results exactly:
90665 eligible,9000 early-primer rejections,335 empty. The private autosomal exon
annotation retains54849 genes. STAR index construction is running;alignment,
gene counts,independent reconstruction and after-execution identities remain
pending. This is process evidence,not a completed measurement result.

Subsequent verified wait:Slurm confirms23232347still RUNNING on ne1dc4-017 after
eight minutes. The native index log records completed suffix-array chunk writes
from its27-chunk construction;sampled peak RSS26516300KiB remains within96GB.
No terminal failure exists,and the same handle is retained without restart or
resubmission. Gene counting and its independent acceptance checks remain pending.
Bounded inspection of the separate CRC truth-construction source confirms that
its existing DNA classes depend on reference selection,GC fitting,masks and run
thresholds. No other-project file/result was modified or adopted as certified
truth. This continuation adds no performance or biological claim.

This pilot extends the technical measurement path only. It does not read DNA
regional profiles, published clone labels or RNA genomic outcomes, or perform
reference fitting, copy-state/clone inference or genomic scoring. The Bian source
remains exposed and secondary/descriptive. Technical cells=1,biological n for
inference=0. The full genomic evaluation freeze and all original error,power,
uncertainty and competitor criteria remain required. Public source and frozen
results are unchanged; the overall goal remains active and unachieved.

## Completed measurement result

Job23232347 completes **0:0** in44m40s,timed process44m35.23s. All193 identities
pass before/after;the unchanged fixtures and complete input/reference hashes pass.
Slurm sampled peak RSS36391068KiB and timed-process maximum37897652KiB are
different memory summaries,not interchangeable measurements.

Of100000 deposited pairs,90665 have eligible prepared R1 records. Every primary
ID is retained exactly once. Native featureCounts and independent reconstruction
agree for all54849 annotated genes,including zeros. There are **26780 assigned
reads across4284 genes**,and **22364 distinct exact gene/tag pairs**. Independent
tag-set counting exactly agrees with the frozen counter. These tags remain
candidate labels;neither chemistry nor independent molecular counts is validated.

Primary outcomes retain44302 unmapped,4668 ambiguous-gene,10622 no-feature and
4293 nonunique records. The7632 secondary/supplementary alignments are auxiliary
records,not additional input reads. STAR reports42365 reads unmapped as too short;
this is an observed processing diagnostic,not proof of low biological RNA content
or a reason to replace the lexical cell with a better-looking one.

The completed computational gates release staging23232546and rejected-tag
ablation23232643. Both now run independently. The ablation checks whether the
9000 excluded R1 reads contain useful splice-supported gene evidence despite
unreliable tags. Reference staging preserves distinct RNA NC426/DNA NC297 IDs;
no paired NC exists in the CRC04 crosswalk. Neither follow-up changes this
producer or its original output.

The separate result entry preserves all435 preceding claim values. No genomic
profile,DNA truth,clone evaluation,false-split/power,uncertainty or superiority
claim follows. Public code and frozen results remain unchanged;the full goal
is active and unachieved.

[Full completed report](../dev/bian_rna_counts_20261112/run_23232347/results/REPORT.json),
[all annotated gene counts](../dev/bian_rna_counts_20261112/run_23232347/results/GENE_COUNTS.tsv),
[resource record](../dev/bian_rna_counts_20261112/run_23232347/resources.txt),
[completion checkpoint](../dev/bian_rna_counts_20261112/completion_checkpoint.json),
[rejected-read matched ablation](BIAN_TAG_READ_ABLATION_2026-11-12.md),
[assay-specific reference candidates](BIAN_REFERENCE_CANDIDATES_2026-11-12.md).

[Frozen design](../dev/bian_rna_counts_20261112/DESIGN.md),
[executed candidate implementation](../dev/bian_rna_counts_20261112/run.py),
[resource header](../dev/bian_rna_counts_20261112/run.sbatch),
[submission](../dev/bian_rna_counts_20261112/submission.json),
[source identities](../dev/bian_rna_counts_20261112/SNAPSHOT.json),
[qualified preprocessing and its source gaps](BIAN_PREPROCESS_2026-11-12.md),
[native synthetic counting fixture](BIAN_COUNT_FIXTURE_2026-11-12.md),
[corresponding DNA pilot](BIAN_DNA_READ_QC_2026-11-12.md).
