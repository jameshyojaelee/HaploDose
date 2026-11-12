# Bian RNA library-protocol qualification

This bounded follow-up checks the published barcode table against the existing
RNA technical-QC observation, before choosing any molecule extraction rule.
It makes no performance or genomic claim. The current software artifact is frozen;
the full measurement candidate/reference/comparator/scorer contract is not.
Bian remains secondary/descriptive and the prior scientific criteria are unchanged.

The [primary Dong2018 methods](https://link.springer.com/article/10.1186/s13059-018-1416-2)
place an eight-base cell barcode and eight-base UMI in the RT primer, identify
read2 as carrying barcode information, and describe gene-level UMI deduplication.
This method was developed for RNA; the mouse reference in that paper is not
appropriate for the human Bian specimen. The primary linked additional file2
contains barcode-RT primer information. A match to that table corroborates the
protocol layout but does not establish Bian-specific extraction, cell identity,
molecular independence, deduplication, or genomic validity.

The exact prespecified RNA library remains CRC04_LN1_349/GSM3241543/SRR7461687.
Its completed producer-QC report found AACAACCA at read2 positions0–7 in each of
the first10000 records. Those records are lexical, not a random biological sample.
The follow-up downloads only the publisher-linked XLSX, retains the exact bytes
and SHA256, parses workbook relationships and shared/inline strings, and reports
all matching cells with addresses and full values. No alternative protocol is
searched to obtain a match. There is no read trimming, UMI extraction, alignment,
gene counting, DNA truth construction, reference selection or outcome scoring.

Job23231112:nslab account/nslab QOS,CPU1CPU1GB1h,no array,existing isolated
environment,neutral directory,142 source identities before/after. Static AST and
shell syntax pass; executable fixtures and retrieval remain pending at submission.
The fixtures check shared strings, inline strings, rich strings and invalid shared
indices. ZIP CRC failure, an empty table or inaccessible attachment must be
recorded as a failure rather than interpreted as a compatible protocol.

The existing I/O job23230977 continues unchanged; reader23230978 waits for it.
No unrelated or running jobs were modified. Biological n=0. Full benchmark and
biological validation remain unachieved.

[Frozen design](../dev/bian_protocol_qualification_20261112/DESIGN.md),
[executable metadata inspection](../dev/bian_protocol_qualification_20261112/qualify.py),
[source identities](../dev/bian_protocol_qualification_20261112/SNAPSHOT.json),
[submission](../dev/bian_protocol_qualification_20261112/submission.json),
[existing raw-input chain](BIAN_RAW_COUNT_CHAIN_2026-11-12.md).

## Completed protocol-table result

Job23231112 completed0:0 on ne1dc6-008 in2s; timed process1.63s and peak21036KiB.
Scheduler RSS was not sampled. All142 source identities pass before and after.
The four executable string-parser fixtures pass. The publisher attachment is
18588 bytes, SHA256
`2e24f8f104b2a2afc9091c4d453c350bc5d5927f9425bd4ee58578ddb6bd1f2b`.
The extracted workbook has494 cells. Sheet `Barcode-RT`,cell `B14` contains the
observed AACAACCA sequence immediately before eight Ns and the polyT tract.
No standalone eight-base cells were found; the match occurs inside the complete
RT-primer sequence. Keeping that distinction avoids interpreting an empty
standalone-barcode list as absence of the barcode from the source.

This independently sourced primer record corroborates the previously observed
read2 structure. It does not validate its full-run error rate or establish the
RNA/DNA specimen's genomic state. It also cannot identify how many distinct
eight-base strings represent independent molecules without gene assignment and
an explicit collision/error model. There is no new recovery, false-split,
coverage, genomic or competitor-performance result. The primary Science
supplement and Bian-specific extraction implementation remain unverified.

[Exact result](../dev/bian_protocol_qualification_20261112/run_23231112/results/REPORT.json),
[complete extracted table](../dev/bian_protocol_qualification_20261112/run_23231112/results/TABLE.json),
[preserved publisher workbook](../dev/bian_protocol_qualification_20261112/run_23231112/results/13059_2018_1416_MOESM2_ESM.xlsx),
[executed fixtures](../dev/bian_protocol_qualification_20261112/run_23231112/results/FIXTURES.json).
