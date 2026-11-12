# From source-identified libraries to usable measurements

The scientific goal remains accurate genomic profiles and supported CNA-defined
populations across RNA modalities, with suitable independent validation. The
[completed Bian metadata gate](BIAN_PAIRING_FEASIBILITY_2026-11-11.md) resolves
609 candidate library pairs, not raw-count provenance or DNA sensitivity. This
bounded technical stage investigates those upstream requirements before a
biological performance experiment. It does not replace any original criterion.

## Fixed source and permitted scope

UseCRC04_LN1_349, already chosen lexically from the UMI-protocol metadata before
outcomes. RNA=GSM3241543/SRR7461687; bisulfite=GSM3153081/SRR7214590. Exactly
four archive-listed files contain2685628275 compressed bytes:295329590 RNA and
2390298685 DNA. The archive byte counts and MD5 values are fixed in the input
manifest. All new files live under a new task directory; no CRC-project or frozen
result is modified.

This stage is download/read-QC only, consistent with the historical unopened
source contract. Do not count genes, align reads, construct DNA truth, select
normal references or score HaploDose. One technical cell cannot establish
population performance, membership reliability or biological false-split control.
Biological n analyzed=0. Bian remains secondary/descriptive.

## Software artifact freeze

The previously validated noneditable0.2.0 wheel has SHA256
`5aeb402d0b3b3cace1000826801b36960053da33d7f5325e8f20b9548898cebb`.
All26 current package modules match it byte-for-byte. Copy that artifact into
the new task and verify the installed isolated package/native binary before
read access. Record default clone/profile parameters, pipeline signature,
environment inventory, source identities and actual backend from a neutral cwd.
No dependency, shared installation, public code or default changes.

CPU job23230976 now completes the fresh import/artifact check:26 modules,
native binary and wheel agree; actual JAX backend is CPU. The source manifest's
136 identities pass before/after. This is a software-artifact freeze, **not** a
claim that the full genomic evaluation protocol or a scientifically validated
v-next is frozen. The old521-Python/four-Rust validation is historical, not rerun.
Genomic evaluation and full measurement analysis remain gated by the existing
candidate/reference/comparator/scorer freeze.

## Decisive technical observations

Raw molecule tags may be retained, stripped or modified upstream. The GEO
protocol's reference to a UMI method cannot establish which of those holds in
the deposited run. The [cited Dong2018 primary RNA protocol](https://link.springer.com/article/10.1186/s13059-018-1416-2)
describes cell barcodes and UMIs in the RT primer and cell segregation using
read2. Its current primary methods text specifies8-base barcode and8-base UMI;
a search-result PDF snippet had a shorter N-run and is not sufficient evidence
for an extraction rule. No barcode/UMI length or read position is adopted for
Bian without run-specific evidence. The Science primary-page request returned
403; no inaccessible supplement or recalled parameter is treated as verified.
The cited RNA protocol uses a mouse reference; that alignment choice is not
transferred to this human dataset. Bian's reported hg19 build still needs an
executed reference/annotation contract before mapping.

The stager checks exact archive bytes/MD5, records SHA256, streams gzip to EOF,
and validates FASTQ framing, sequence/quality lengths, printable qualities and
paired identifiers. It records read counts, length distributions, quality ranges
and first10000-record head-position composition. These are technical quantities,
not UMI counts or independent DNA fragments. Composition sampling is lexical,
not a random biological sample. No QC-percentage threshold chooses an outcome.

Independent GNU md5sum/gzip/wc checks must reproduce each file's checksum,
compressed integrity and exactly4 times the recorded read count in decompressed
lines. Both subprocess failure codes are checked. Fixtures reject truncation,
length mismatch and discordant mates. A checksum or format failure stops the
chain. Partial downloads keep an explicit `.partial` suffix and are not overwritten
or silently admitted. Network failure does not authorize a duplicate live job.

Passing those tests establishes intact raw files only. UMI retention/extraction,
gene counting and fragment deduplication require additional source-specific
qualification; independently reconstructed DNA sensitivity is still needed for
genomic interpretation. More reads do not supply missing biological replication.

## Compute and current dependency chain

| Job | Purpose | Resources | Dependency |
|---|---|---|---|
|23230976|Candidate artifact/import check|nslab CPU1CPU4GB1h|None|
|23230977|Four-file staging and read QC|nslab I/O1CPU4GB2h|afterok23230976 plus actual frozen-artifact record|
|23230978|Independent compressed/checksum/line verification|nslab CPU1CPU2GB1h|afterok23230977 plus complete QC record|

All account/QOS values are nslab. No arrays, interactive QOS, GPU allocation,
new dependencies, shared installs or unrelated-job changes. Existing isolated
Conda/runtime, one thread, neutral directories and136 before/after identities.
CPU and I/O partition limits were checked before submission. The I/O job is
specifically for staging and streaming these2.69GB compressed read files.
Static Python AST and shell syntax pass. Raw staging and the independent reader
are pending completion at this checkpoint; successful artifact import is not
evidence of biological performance. The full goal remains active and unachieved.

[Frozen design](../dev/bian_raw_count_chain_20261112/DESIGN.md),
[input/count-scale/selection contract](../dev/bian_raw_count_chain_20261112/INPUTS.json),
[artifact checker](../dev/bian_raw_count_chain_20261112/freeze.py),
[staging and FASTQ checks](../dev/bian_raw_count_chain_20261112/readqc.py),
[independent reader](../dev/bian_raw_count_chain_20261112/reader.py),
[source identities](../dev/bian_raw_count_chain_20261112/SNAPSHOT.json),
[actual CPU artifact record](../dev/bian_raw_count_chain_20261112/freeze_23230976/results/CANDIDATE_FROZEN.json),
[exact submissions](../dev/bian_raw_count_chain_20261112/submissions.json).

## RNA producer-QC checkpoint; DNA download remains live

Both RNA files match archive bytes/MD5. The producer's complete gzip/FASTQ pass
finds1842166 paired reads, exactly150 bases in both mates, with printable quality
code ranges35–74. Its truncation, length-mismatch and discordant-mate fixtures
pass. The independent GNU reader still waits for the full staging job; do not
describe the complete chain as passed yet.

[Observed pattern, not a guarantee] In the first10000 records, read2 positions0–7
are exactly`AACAACCA` in every record. Positions8–15 vary, followed by a T-rich
region beginning at16. This is consistent with the cited RNA protocol's barcode/
UMI layout and argues against assuming that all molecule-tag sequence was removed
from this run. It does not validate an extraction rule, molecule deduplication,
independence, complete-run contamination rate or a genomic measurement. No tags
are extracted and no reads trimmed or mapped in this stage. The sampled records
remain lexical; no population error rate or biological inference is estimated.

[Completed RNA producer record](../dev/bian_raw_count_chain_20261112/readqc_23230977/results/SRR7461687_REPORT.json).

Artifact job23230976 completed0:0 in30s; timed process29.18s, peak423248KiB.
Scheduler RSS was not sampled. At the saved live checkpoint, I/O23230977 is
RUNNING, DNA R1 is still an explicit partial file, and reader23230978 is PENDING
for dependency. Live sampled staging RSS is32288KiB. The running job is neither
restarted nor duplicated. Its exact byte progress and UTC checkpoint are saved
separately, because a partial download is not a verified input or completed QC.
One running,one pending,one completed,zero errors in this chain.

## Later checkpoint: transport timeout preserved; new continuation

The original live stage subsequently failed1:0 with a socket read timeout during
DNA mate2 after830472192 bytes. Three completed files passed archive checks;
paired DNA QC and independent verification did not run. Pending reader23230978
was canceled because its failed dependency could never satisfy. Original logs,
source,files and partial remain unchanged. A separate checked range continuation
23231469 and dependent reader23231472 preserve all original admission criteria.
[Exact terminal evidence,resources and continuation](BIAN_RAW_RESUME_2026-11-12.md).
