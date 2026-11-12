# Bounded corresponding-DNA read QC

This experiment checks whether the corresponding bisulfite DNA library supplies
usable paired alignments under a declared processing path. It does not construct
DNA copy-number truth or score HaploDose. The Bian source remains exposed,
secondary and descriptive. Biological replication for inference is zero; one
technical cell is inspected.

## Frozen experiment

Cell CRC04_LN1_349, GSM3153081/SRR7214590, was selected before new outcomes in the
earlier correspondence investigation. Use its first 100,000 deposited read pairs.
Rehash both complete input FASTQs and all 15 previously qualified hg19 reference
assets before analysis. No RNA grouping or outcome is used.

Keep both front-trimming arms, zero and nine bases. Nine bases comes from the
existing CRC adaptation, not a verified original Bian parameter. Both arms use
Cutadapt 4.4, exact full 13-base AGATCGGAAGAGC adapter matching, quality 20 and
minimum length 50. These are technical working settings, not validated biological
thresholds. Bismark 0.24.1 uses paired nondirectional mapping with Bowtie2 and
one parallel instance. Paired positional deduplication precedes name sorting.
Samtools 1.22.1 explicitly sorts lexicographically with `-N`.

The strict independent fragment counter runs at MAPQ 10, 20 and 30 before and
after deduplication. MAPQ 20 is the primary technical screen; other values are
descriptive sensitivities. All arms, empty outputs and rejection reasons remain
in the report. No mapping-rate threshold or convenient winning arm is selected.
Internally computed regional bins are discarded without inspection or export.
There is no normal selection, callable genomic mask, segmentation, copy state,
clone label or genomic outcome comparison in this task.

The smallest computational gate roundtrips the existing independent synthetic
SAM fixture through native BAM conversion and name sorting, then requires exact
equality of the complete counter result. A failed check stops real processing.
The counter's earlier compute-node fixtures independently establish reciprocal
mate-coordinate handling, CIGAR/query/reference bounds, actual name ordering,
MAPQ missingness and empty-input accounting. They do not establish molecule
independence or historical specimen errors.

## Execution checkpoint

Job **23232222** uses account/QOS nslab, CPU partition, eight CPUs, 64 GB,
two hours, constraint `v6&amd`; no array or downstream job. It runs on
ne1dc6-006. All **185** frozen source identities pass the initial check. The
native fixture and full input/reference identity checks pass. Mapping is
running at this checkpoint; completion and scientific outcomes are not claimed.
The first sampled process memory is 16,768,660 KiB. Final runtime, memory,
terminal state and the after-execution identity check remain pending.

## Completed result

Job23232222 completes **0:0** in3m55s; timed process3m52s. All185 identities
pass before and after. Slurm sampled peak RSS is16,864,972KiB; the timed parent/
child-process maximum is3,408,776KiB. These are different memory summaries and
should not be conflated. No pending or running task/dependent job remains.

| Fixed front trim | Retained pairs | Unique mapped pairs | Deduplicated mapped pairs | Autosome MAPQ10 | Autosome MAPQ20 | Autosome MAPQ30 |
|---|---:|---:|---:|---:|---:|---:|
| 0 | 97,984 | 21,104 | 20,978 | 9,266 | 9,179 | 2,770 |
| 9 | 97,393 | 28,440 | 28,252 | 22,920 | 22,653 | 16,063 |

The three final columns are after deduplication; complete before/after rejection
categories remain in the JSON. Bismark's unique-pair counts and deduplication
report agree with the strict counter's query-group totals in both arms. Native
fixture, BAM integrity, reference dictionary, name ordering and reciprocal-mate
checks pass. There are no secondary/supplementary records in these mapped BAMs.
Low mapping quality explains the remaining counted rejections in this subset.

[Observed pattern, not a guarantee] The fixed nine-base intervention increases
the number of accepted high-quality paired observations. This is consistent
with a processing-sensitive mapping problem, but does not establish the original
library's exact primer chemistry, molecular independence, genomic uniformity,
normalisation validity or DNA sensitivity at a particular CNA resolution.
Mapping yield alone cannot select a scientifically valid truth-processing
default. No regional DNA values were inspected/exported and no claim of improved
HaploDose accuracy follows. Full-library processing, suitable normal controls,
callable regions, resampling and uncertain truth construction remain unrun.

The completed technical result is appended separately to the claim ledger;
the prior design entry and all430 preceding values remain unchanged. The full
goal remains active and unachieved.

[Full executed result](../dev/bian_dna_readqc_20261112/run_23232222/results/REPORT.json),
[resource record](../dev/bian_dna_readqc_20261112/run_23232222/resources.txt),
[completion checkpoint](../dev/bian_dna_readqc_20261112/completion_checkpoint.json).

## Relation to the scientific review

The [November 10 review](scientific_code_review_2026-11-10.md) correctly requires
measurement provenance and truth at the claimed resolution. This experiment
addresses a bounded technical prerequisite to those requirements. It cannot
resolve RNA/DNA nonidentifiability, validate individual-cell membership, provide
independent null specimens or repair the false-split and uncertainty failures.

A new small-file hash reconciliation finds all **91** protected package, native,
data and test identities unchanged from the source with recorded 521-Python and
four-Rust validation. That is source correspondence, not a fresh test execution.
Those public suites were not rerun for private read-QC code. Every one of the
429 prior claim-ledger entries remains value-identical after appending the
technical design entry. Frozen results, production thresholds and scientific
criteria are unchanged.

[Prespecified design](../dev/bian_dna_readqc_20261112/DESIGN.md),
[executed implementation](../dev/bian_dna_readqc_20261112/run.py),
[resource header](../dev/bian_dna_readqc_20261112/run.sbatch),
[submission record](../dev/bian_dna_readqc_20261112/submission.json),
[source identities](../dev/bian_dna_readqc_20261112/SNAPSHOT.json),
[native fixture result](../dev/bian_dna_readqc_20261112/run_23232222/results/NATIVE_FIXTURE.json),
[protected source reconciliation](../dev/bian_dna_readqc_20261112/review_current_source_check.json),
[independent fragment-counter report](BIAN_DNA_COUNTER_2026-11-12.md).
