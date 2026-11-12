# Bian raw staging: preserved timeout and checked continuation

Original staging23230977 is terminal FAILED1:0 after30:15.70. The failure is
`TimeoutError: The read operation timed out` while reading DNA mate2,after
830472192 written bytes. Three completed FASTQs retain their archive checks;
DNA paired-read QC and the independent reader did not complete. Process peak
40300KiB; Slurm sampled32484KiB. This is a transport failure,not a failed
scientific result or evidence that the archive is inaccessible. The original
code,logs,completed files and partial remain unchanged.

Reader23230978 was pending with DependencyNeverSatisfied and was canceled.
No running or unrelated job was changed. Continuation was submitted only after
authoritative Slurm records confirmed the original parent terminal.

## Exact continuation contract

Job23231469 rehashes the three completed files against original archive bytes/
MD5 and recorded SHA256 and reuses them read-only. It copies the partial into
a new exclusive output,records prefix SHA256 and requires unchanged source size/
mtime. One HTTP Range request asks for missing bytes only. Require206,exact
Content-Range start/end/total and exact remaining Content-Length. Fixtures reject
ignored ranges,wrong offsets/totals and wrong lengths. The180s socket timeout
replaces the failed60s technical allowance; no scientific criterion changes.
There is no hidden retry,overwrite or duplicate live job.

Only full archive bytes/MD5 admit the new file; SHA256 is recorded. The unchanged
strict DNA paired FASTQ/gzip checks must then pass. Original RNA QC is reused
only after exact input rehash. A complete report appears only after all checks.
On failure the new partial remains unverified. Independent reader23231472 must
then pass unchanged external md5sum,gzipCRC and decompressed line checks for all
four files. A copied prefix or HTTP response is not a passed file or QC result.

| Job | Resources | Dependency |
|---|---|---|
|23231469|nslab/nslab I/O1CPU4GB1h|Original parent terminal FAILED|
|23231472|nslab/nslab CPU1CPU2GB1h|afterok23231469 plus complete producer report|

All146 source identities are frozen; static AST/shell pass. Existing isolated
runtime,neutral directories,no arrays,new dependencies,shared installs,public
changes or frozen-output edits. Biological n=0. At the live checkpoint the range
response is206 with exact contract,fixtures pass,and the new partial grows;
full archive verification,paired QC and independent reader remain pending.

Scope remains download/read-QC only. No real trimming,UMI extraction,gene counts,
DNA fragments/truth,normal reference or genomic outcomes. Preserve the original
candidate/reference/comparator/scorer freeze,Bian secondary/descriptive role and
all recovery,false-split,power,uncertainty and competitor criteria. Full goal
active and unachieved.

[Frozen design](../dev/bian_raw_resume_20261112/DESIGN.md),
[checked continuation](../dev/bian_raw_resume_20261112/resume.py),
[source identities](../dev/bian_raw_resume_20261112/SNAPSHOT.json),
[actual submissions](../dev/bian_raw_resume_20261112/submissions.json),
[preserved failure](../dev/bian_raw_count_chain_20261112/readqc_23230977/output.txt),
[original stage and independent criteria](BIAN_RAW_COUNT_CHAIN_2026-11-12.md).

## Four archive files now verified; complete paired scan remains live

The new DNA mate2 is1224034324 bytes and passes its original archive MD5. Its
exact HTTP206 response covered830472192–1224034323 of1224034324 bytes.
All four FASTQs therefore have archive byte/MD5 verification. The retained old
partial was not changed. DNA paired FASTQ scanning remains live in23231469;
reader23231472 remains dependency-pending. The full QC chain has not passed.
Sampled live RSS is28988KiB. This checkpoint is file-integrity evidence only;
no molecule counts,DNA truth,biological result or package-performance gain.

[Full resumed-file identity](../dev/bian_raw_resume_20261112/resume_23231469/results/R2_VERIFIED.json),
[exact range response](../dev/bian_raw_resume_20261112/resume_23231469/results/RANGE_RESPONSE.json),
[retained source-prefix identity](../dev/bian_raw_resume_20261112/resume_23231469/results/RETAINED_PREFIX.json),
[reverified completed files](../dev/bian_raw_resume_20261112/resume_23231469/results/REUSED_FILES.json).

## Producer completes; independent reader starts

Job23231469 subsequently completed0:0 in6:00.61,process peak37696KiB,Slurm
sampled28988KiB. All146 source identities pass before/after. The complete DNA
scan finds12334976 paired reads,all150 bases per mate,with printable quality
codes35–74. DNA paired identifiers,sequence/quality lengths,gzip integrity and
both archive MD5s pass. Original RNA remains1842166 paired150-base reads after
exact input rehash. These read counts are not independent DNA fragments,UMIs,
genomic coverage or biological replication. No deduplication or DNA sensitivity
claim follows. Independent reader23231472 is running;full independent chain
remains unverified until that job and its external checks pass.

[Complete producer report](../dev/bian_raw_resume_20261112/resume_23231469/results/REPORT.json),
[DNA paired-read record](../dev/bian_raw_resume_20261112/resume_23231469/results/SRR7214590_REPORT.json).

## Complete independent raw-input QC passes

Reader23231472 completed0:0 in1:25; Slurm sampled7036KiB. All146 source identities
pass before and after both the producer and reader. External md5sum and gzip
integrity checks pass for all four files. Each RNA mate has7368664 decompressed
lines and each DNA mate49339904 lines,exactly four times the producer read counts.
The independent tool versions and all per-file checks are retained in the result.
No job in this replacement chain remains pending or running. The original
FAILED1:0 stage and canceled impossible reader remain part of the record.

This completes the bounded raw-file integrity task. There are1842166 RNA pairs
and12334976 bisulfite DNA pairs for the already prespecified technical cell.
These counts do not establish usable unique fragments,DNA resolution or CNA
truth. No trimming,UMI extraction,gene counts,alignment,DNA deduplication or
genomic outcome has been run on the specimen. The reference contract and synthetic
counting component passed separately; real measurement construction still needs
its complete frozen protocol. All original scientific and release criteria,
frozen artifacts and historical decisions remain unchanged.

[Independent completed result](../dev/bian_raw_resume_20261112/reader_23231472/results/REPORT.json),
[producer source checks after execution](../dev/bian_raw_resume_20261112/resume_23231469/source_after.txt),
[reader source checks after execution](../dev/bian_raw_resume_20261112/reader_23231472/source_after.txt).
