# DEFND bounded measurement identity: real provenance gates remain closed

The scientific goal remains accurate genomic RNA dosage profiles and supported
CNA clones/memberships, with independent superiority over the fixed competitors.
This work implements and executes the next measurement-chain prerequisite from
the deep review. It finds a real unresolved specimen-name discrepancy and
barcode-orientation uncertainty. It delivers no new accuracy improvement or
biological validation. Original recovery2/5,continuous-state errors14/200,
matched-generator power195/200 and selected coverage54/1000 remain unchanged.

## What changed and what the observations establish

Implemented private bounded readers, suffix exclusion, strict source-supported
header parsing, explicit whitelist diagnostics and a fail-closed specimen-alias
check. No public package source/default changed. These are technical admission
tools, not an integrated count-matrix loader or a new inference method.

The lexical diagnostic remains GSM7016419/TB5471, fixed before inspection.
GSM7817780 remains the previously declared second specimen and is unopened here.
The RNA header contains two annotation columns and1821 cell columns. All1821
sanitized IDs are distinct16nt sequences, separated from opaque author suffixes
by underscores; all are literal forward ARC-v1 GEX whitelist members. Author
cluster-number values are never interpreted, exported or saved. No expression
row is decompressed.

The exact archived TB5471 DNA file's header says **id=PTO026**. It reports
pipeline_name=dna10x,pipeline_version=dna10x-1.0.0 and194 primary contigs, with206
comment lines occupying5548 bytes. Its reference_path is
/opt/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa. Reference FASTA/GTF
hashes,reference_version and mkref_version are empty. A declared path cannot
certify reference bytes or prove the historically executed command.

The first50 complete DNA fragment rows have50 distinct literal16nt IDs. The
whole prespecified orientation/reference table is:

| Input | Literal forward ATAC | Literal forward GEX | Reverse-complement ATAC | Reverse-complement GEX |
|---|---:|---:|---:|---:|
|1821 RNA header IDs|0|1821|0|0|
|50 DNA prefix IDs|0|0|45|0|

[Inference] The DNA pattern is consistent with a reverse-complement ATAC
namespace for45 inspected sequences. It does not establish the executed
orientation/correction settings; five sequences remain outside those references.
No stored IDs are transformed, direction selected, barcode corrected, bead
translation applied or RNA/DNA intersection performed. A coordinate-sorted
50-row prefix cannot estimate library-wide barcode matching or clone diversity.

All13 cached GEO sample records were scanned for exact PTO026/TB5471 tokens in
identity fields, excluding shared processing/protocol text. The only match is
the TB5471 RNA supplementary filename in GSM7016419. No explicit alias binding
was found in those fields. The linked
[BioSample record](https://www.ncbi.nlm.nih.gov/biosample/SAMN32979122) confirms
the accession and generic glioblastoma source; its inspected page does not supply
that alias. The direct SRA browser request was inaccessible, not evidence that
metadata are absent. This is **not proof of a specimen mix-up**. The alternative
is a harmless historical processing alias, which requires an independent binding.
The new literal-name/alias gate fails for this specimen; no binding was supplied.

| Assessment | Estimand/model and current status | Decisive independent evidence |
|---|---|---|
|Computational correctness|Exact molecular prefixes,archive ranges,writer fields and vendor membership; source-supported checks pass; initial inspection failure preserved|Independent fixtures,identical request hashes,before/after manifests and vendor references|
|Statistical identification|Library-qualified physical-cell correspondence remains unresolved; matching chemistry sequences alone is insufficient|Original specimen/library manifest,executed whitelist/orientation,barcode correction and translation history|
|Biological validity|No DNA states,resolution,memberships or performance qualified; biological validation n=0|Corresponding-cell independently constructed DNA evidence with documented sensitivity,plus suitable biological replication|

## Every attempt and the explicit protocol changes

1. [Technical protocol version2](../dev/defnd_bounded_identity_20261112/DESIGN_V2.md)
   was frozen before access. It explicitly changes the earlier ordering that
   required an actual grammar and full evaluation design before any header
   inspection. A bounded,label-excluding schema investigation now precedes the
   evaluation freeze. The old protocols remain unchanged. Full candidate,
   comparator,endpoint,reference,callability and failure design must still freeze
   before RNA values,DNA profiles,predictions or scores are read. This does not
   admit a cohort or weaken a biological error/power requirement.
2. Qualification23238705 passes1006 suffix perturbations and13 rejection cases,
   independent compression/fragment fixtures and exact range handling. Invalid
   response headers cause zero body reads. All56 source identities match.
3. Real inspection23238706 fails1:0: an unrecognized RNA prefix and the fixed
   100-comment DNA cap. Both servers honor206 and exact Content-Range. Preserve
   that code,exit status,request records and opaque compressed DNA prefix. No
   raw RNA bytes/decorated fields are saved. No silent parser retry occurs.
4. Separate [diagnostic](../dev/defnd_bounded_identity_20261112/DIAGNOSTIC_DESIGN.md)
   23238729 passes1000 suffix perturbations, reports only structural failure
   classes, and repeats the exact RNA range because decorated source was not
   retained. Its prefix hash matches. It examines only the original first100
   DNA comments; those are not recognized by its closed key vocabulary. No
   additional DNA request or fragment acceptance follows. All79 identities match.
5. Inspect the author-declared dna10x repository's README,driver,address,
   error-correction and count/writer sources. All five returned blobs match
   independently derived Git blob identities. Current writer code documents ten
   key=value fields and one primary_contig comment per reference sequence.
   Current address code applies Tn5 endpoint shifts; correction can retain
   unknown IDs, and the driver permits reverse-complementing the whitelist.
   None proves the actual historical version/options. Do not execute,install,
   fix or silently transfer that pipeline's observation model into HaploDose.
6. Separate [source-supported version3](../dev/defnd_bounded_identity_20261112/STRUCTURE_DESIGN_V3.md)
   23238753 passes. It fixes the independently declared two-annotation contract
   and requires every cell column to satisfy the observed exact grammar. Missing
   annotations,malformed/shortened IDs and collisions reject. Replace only this
   version's100-comment resource cap with1MiB comment bytes,accepting only the
   documented keys. The1MiB compressed DNA/16MiB inflated caps and50-row target
   stay. This resource-policy change is explicit; the old failure stays failed.
   200 suffix perturbations,12 rejection cases,an independently specified
   200-contig/50-row fixture and the original complete qualification all pass.
   All95 source identities and both vendor hashes match before/after. RNA's
   final declared repeat prefix is identical; DNA uses the original cached bytes.
7. Separate [cached reconciliation](../dev/defnd_bounded_identity_20261112/RECONCILIATION_DESIGN.md)
   23238773 passes independent complement examples,exact metadata-token
   boundaries,shared-protocol exclusion and five invalid-alias cases. No
   network requests. All112 source identities and vendor hashes match. It
   reports the full orientation table without choosing a direction and rejects
   the unsupported alias. No more barcode rescue/search is justified here.

## Exact execution and inputs

Every job uses nslab account AND QOS,cpu partition,1 CPU,4GiB and10min,no arrays.
No unrelated running or pending jobs are changed. The existing isolated Conda
environment/runtime is used read-only with neutral output directories,Python-I,
one thread and no installs. The explicit execution backend is Python standard
library,not inference/CUDA. Peak memory below is process /usr/bin/time-v RSS;
Slurm's zero sampled RSS is unavailable evidence, not zero memory use.

| Job | Terminal progress | State/exit | Node | Slurm/timed runtime | Peak KiB |
|---|---|---|---|---|---:|
|23238705 qualification|1/1 ████████|COMPLETED0:0|ne1dc4-046|5s/1.93s|16140|
|23238706 first inspection|1/1 ████████|FAILED1:0,preserved|ne1dc4-046|5s/3.88s|233456|
|23238729 diagnosis|1/1 ████████|COMPLETED0:0|ne1dc4-002|1s/1.32s|27264|
|23238753 source reader|1/1 ████████|COMPLETED0:0|ne1dc5-006|3s/2.53s|333472|
|23238773 reconciliation|1/1 ████████|COMPLETED0:0|ne1dc4-063|3s/2.07s|229952|

Running0,pending0,completed4,failed1. The dependent inspection ran only after
qualification succeeded and repeated qualification before requesting data.
No downstream inference/evaluation job is admitted.

RNA source: same public NCBI FTP host/path as the locked provider URL,with ftp
normalized to HTTPS. Each of three explicitly declared requests is exactly
bytes0–524287 of an8922600-byte artifact; every prefix SHA256 is
cf1655abbe1c79d8a9656c9b7cc6f9cd5b5f8487aba6b3ff3ccb124e31dc9afa.
DNA source: exact locked Zenodo content URL,bytes0–1048575 of13318383087 once;
prefix SHA256 bba059bff1eb1c880d470be1f38e5543bcd713c2850d7f303f41b0c55fb1568d.
These are partial-byte identities,not full-file checksum or compression-integrity
verification. RNA compressed body may contain later compressed data; expression
rows are never inflated or analyzed. No whole13.3GB download or staging follows.

All source,design,request,qualification,sanitized identity,resource and environment
records are under [the execution directory](../dev/defnd_bounded_identity_20261112/).
[Frozen author source](../dev/defnd_bounded_identity_20261112/PUBLIC_DNA_SOURCE.json),
[final structural result](../dev/defnd_bounded_identity_20261112/structure_23238753/results/REPORT.json),
[reconciliation result](../dev/defnd_bounded_identity_20261112/reconcile_23238773/results/REPORT.json).
The current-source observations trace to the
[author's dna10x repository](https://github.com/simslab/dna10x); exact blob
identities/content are retained locally. The official SAM specification was
inspected as a possible header explanation; that possibility was not adopted.

Executed launcher commands are the recorded sbatch files with their absolute
task path; inspection used afterok:23238705 and kill-on-invalid-dep=yes. Static
Python AST and bash-n checks pass. Public pytest/Rust suites are not rerun:
public source/native paths did not change. Their existing recorded validation
is not relabeled as new execution.

## Decision and next evidence

Keep this specimen out of same-cell evaluation until an independent manifest
binds PTO026 to TB5471/GSM7016419 and the original whitelist/orientation,
correction and translation history are supported. The async user request asks
for an existing local archive path, not permission to weaken a gate. Five unknown
prefix barcodes cannot be repaired by selecting their nearest reference sequence
or choosing the direction with maximum RNA/DNA agreement.

After that technical gate, independent DNA normalization/callability,copy-state
truth,assay resolution and biological replication still need qualification.
More computation cannot substitute for those measurements. No accuracy score,
clone-recovery success,5% biological guarantee,uncertainty calibration or SOTA
claim follows. Frozen results and all previous decisions remain unchanged.
