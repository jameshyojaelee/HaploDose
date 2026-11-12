# Assay-specific reference candidates remain distinct cells

## Completed continuation, 2026-11-12

Job **23232996 completes0:0 in26m28s**. All **223** source identities agree
before and after. Both RNA mates reverify bytes/MD5/SHA256 and complete paired
QC after exact JSON key normalization; positive/negative QC fixtures pass. The
copied DNA R1 prefix is retained independently and the exact206 range contract
passes. Both complete DNA mates match original archive bytes/MD5, receive SHA256
identities and pass full paired validation. Native md5sum and gzip checks pass
for all four complete files. Both failed parents remain untouched.

Result: **1,243,474 RNA pairs** for CRC04_NC_426/SRR8096032 and
**13,643,476 DNA pairs** for CRC04_NC_297/SRR7214625, all150bases in both mates.
Timed peak RSS is42,968KiB; Slurm sampled28,828KiB. The two NC candidates remain
distinct cells, not a same-cell pair, certified normal baseline, genomic truth
or independent biological replication. Full genomic admission remains false.

All five full-library RNA/DNA jobs23233126–23233130 have now started, each with
257 passed source identities. The four DNA native known-answer fixtures pass.
Full measurements and genomic scores remain pending. The frozen results,
original recovery/error/power/coverage and comparator requirements are unchanged.

[Completed staging report](../dev/bian_reference_resume_v2_20261112/run_23232996/results/REPORT.json),
[native four-file checks](../dev/bian_reference_resume_v2_20261112/run_23232996/results/EXTERNAL_CHECKS.json),
[completion and resource record](../dev/bian_reference_resume_v2_20261112/completion_checkpoint.json),
[full-library measurement stage](BIAN_FULL_MEASUREMENTS_2026-11-12.md).

The following submission and failure checkpoints preserve their original context.

## Continuation and serialization correction

V2 job **23232996** is submitted with account/QOS nslab, I/O partition, one CPU,
2 GB and two hours. It freezes **223** source/provenance identities. Its result
is pending at this checkpoint; no complete-file or scientific pass is inferred
from partial transfer progress. Static AST and launcher syntax checks pass.
The earlier failures and stronger biological requirements remain below.

Fresh214-identity continuation **23232945 failed1:0 in17s**, before any transfer.
Both completed RNA files pass bytes/MD5/SHA256 again. Recomputed complete QC
then differs from the archived JSON in integer versus string read-length keys.
The full-object equality correctly stops the workflow. Preserve its source,
output and resource record. V2 compares the entire recomputed QC after JSON
round-trip, with a fixture verifying key normalization and another rejecting
changed counts. No QC value or archive criterion is relaxed. The fresh v2
retains nslab I/O1CPU2GB2h and all original selection and biological limits.

[First continuation design](../dev/bian_reference_resume_20261112/DESIGN.md),
[failed continuation](../dev/bian_reference_resume_20261112/run_23232945/output.txt),
[v2 design](../dev/bian_reference_resume_v2_20261112/DESIGN.md),
[v2 exact comparison and transfer](../dev/bian_reference_resume_v2_20261112/run.py),
[v2 actual submission](../dev/bian_reference_resume_v2_20261112/submission.json).

## Terminal staging failure, 2026-11-12

Job **23232546 failed 1:0**, after 10m03s, with `TimeoutError` while reading the
DNA R1 HTTPS response. Both RNA mates had passed exact archive byte/MD5 checks,
strict full paired FASTQ validation (1,243,474 pairs) and native md5sum/gzip
checks. The retained DNA R1 partial is 872,415,232 bytes against 1,181,548,139
expected bytes. It is not a verified complete file; DNA R2 was not transferred.
The original run, traceback, completed RNA report and partial remain untouched.
The 203-entry source check passed before execution; the failure prevented the
after-check and overall report, so full staging is **not complete**.

This is a network-transfer failure, not evidence about DNA quality, diploidy,
reference suitability or biological replication. The separate-cell identities,
lexical selection and all admission criteria remain unchanged. Any continuation
must use a fresh destination and verify exact archive identity before QC.

[Failed run](../dev/bian_reference_candidates_20261112/run_23232546/output.txt),
[completed RNA-only report](../dev/bian_reference_candidates_20261112/run_23232546/results/SRR8096032_REPORT.json).

The next reference investigation requires baseline material. Inspection of the
frozen609-pair crosswalk finds57 CRC04 matched RNA/DNA pairs and **zero matched
CRC04 NC pairs**. An initial small metadata-inspection attempt indexed the empty
NC list and raised `IndexError` before any transfer or computation. Subsequent
inspection of the complete crosswalk confirmed the actual absence; it was not
a missing-file download failure. No same-cell normal pair is fabricated.

The separate frozen ENA table provides NC-labelled RNA and DNA libraries with
different physical-cell identifiers. A new design stages the lexical first
exact title in each assay among that table's available rows. It retains all
candidate identifiers and forbids outcome-based replacement after failure.
This is preparation for a later reference-sensitivity experiment, not admission
of a normal reference or biological null panel.

| Role candidate | Author cell ID | GEO sample | Raw run |
|---|---|---|---|
| RNA reference candidate | CRC04_NC_426 | GSM3443126 | SRR8096032 |
| DNA reference candidate | CRC04_NC_297 | GSM3153116 | SRR7214625 |

The RNA library is declared RNA-Seq/TRANSCRIPTOMIC and the DNA library
Bisulfite-Seq/GENOMIC; both are paired-end. Total four-file compressed bytes:
**2,666,506,201**. These are two technical cells from the same declared patient,
not one paired cell or independent biological replicates. NC is an author label;
diploidy and single-CNA normality remain unverified. Separate-cell DNA cannot
validate the RNA reference cell's individual genotype or clone membership.

## Executable checks and dependencies

Reuse the original exclusive archive downloader and strict paired-FASTQ reader,
including unchanged malformed/truncated/discordant-mate fixtures. Require exact
provider byte counts/MD5s, record SHA256 and complete FASTQ conservation/quality
summaries, then independently run native md5sum and gzip integrity checks. Keep
all partial failures, with no hidden retry, overwrite or change of candidate.
First10000-read head composition is technical description; it does not prove
tag chemistry, molecule independence or normal identity.

The job waits for successful RNA pilot23232347. It additionally requires that
pilot's final primary-ID, independent tag, native gene-read and annotated-zero
checks, positive transcript evidence, and identical before/after source records.
Positive evidence here means nonempty assignment only; no selected mapping-rate
threshold is introduced. The existing candidate software freeze must remain
present. These are computational prerequisites, not genomic validation.

Job **23232546** requests account/QOS nslab, I/O partition, one CPU,2GB,two hours,
dependency `afterok:23232347`,203 source identities. Streaming transfer/rehashing
and complete FASTQ inspection justify I/O. It is initially pending Dependency;
no transfer/read-QC result is claimed. The original RNA job remains running and
is not modified. No array,interactive QOS,new dependency or shared installation.
Static AST and launcher syntax pass; executable fixtures run on compute nodes.

## Scientific limits and next decision

This stage performs no real tag extraction, gene counting, alignment, genomic
binning, normalization, segmentation, DNA truth or HaploDose scoring. Bian stays
exposed and secondary/descriptive. Any later normalization must separately
justify reference composition, genomic callability and DNA baseline stability.
The alternative explanation remains reference-cell or tissue-state difference;
corresponding DNA and reference perturbation would help distinguish it from
genomic dosage. These two filesets alone cannot resolve that alternative.

A failed reference candidate stays failed and does not trigger a search for one
that better matches RNA/DNA outcomes. A completed download establishes input
integrity only. Full candidate/reference/comparator/scorer freezing and all
original recovery,false-split/power,uncertainty and competitor criteria remain.
No biological n for inference,genomic negative-control certification,production
promotion or frozen-result change follows. All433 preceding ledger values remain
unchanged; the separate new entry records this design and pending job.

[Frozen design](../dev/bian_reference_candidates_20261112/DESIGN.md),
[exact inputs and complete candidate IDs](../dev/bian_reference_candidates_20261112/INPUTS.json),
[executed code](../dev/bian_reference_candidates_20261112/run.py),
[resource header](../dev/bian_reference_candidates_20261112/run.sbatch),
[submission/dependency](../dev/bian_reference_candidates_20261112/submission.json),
[203 source identities](../dev/bian_reference_candidates_20261112/SNAPSHOT.json),
[complete frozen correspondence](../dev/bian_pairing_feasibility_20261111/crosswalk_23230821/results/REPORT.json),
[prerequisite RNA pilot](BIAN_RNA_COUNTS_2026-11-12.md).
