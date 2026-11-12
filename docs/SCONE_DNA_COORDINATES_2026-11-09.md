# Published DNA has variable intervals and limited chromosome support

The scONE cell-line normalized DNA table contains 5,393 variable-length intervals,
with no interval exactly 500 kb. Its coordinate representation is consistent
with one-based closed intervals. All autosomes reach their audited GRCh38 FAI
ends, but the X block stops at 114,820,832; Y and mitochondrial DNA are absent.
These are properties of published-table support, not observed neutral DNA or a
validated biological deletion.

The source is the exact previously header-locked
GSE185269_scONEseq_Cellline_scWGS_normalized_counts_matrix.txt.gz. The fixed pilot
label remains One.HCT.01_dedup, supported by original SRA filenames before any
assay comparison. No mapping is selected from profile agreement.

## Scientific contract

The estimand is the deposited physical interval footprint and row identity. The
observations are CHR, START, END and feature fields. DNA normalized values and
RNA values are not parsed. The retained compressed archive includes other assay
columns; its full hash is recorded without treating those values as inspected.

Correct interval computation, statistical interpretation of normalized DNA and
biological validation are separate. This analysis establishes the first. It does
not establish the original reference sequence, unique-mappability mask, counting
units, normalization configuration, baseline ploidy or DNA detection sensitivity.
Published-table coordinates cannot certify that every included base was observed.

The strongest competing explanation for an RNA–DNA profile discrepancy at this
stage is incompatible physical support or upstream processing, rather than an
RNA model defect. A decisive comparison must first allocate independently
processed DNA to the actual intervals, retain outside-support observations, and
justify callable support and normalization. A fixed-bin array index comparison
would not test the intended regional estimand.

## Primary source and independent coordinate checks

Four current Ginkgo genome-bin source files have exact Git blob checks. Its
[documentation](https://github.com/robertaboukhalil/ginkgo/blob/master/genomes/scripts/README.md)
supports both fixed bins and bins based on uniquely mappable
sequence. Three additional interval-source files also pass exact blob checks.
The [process.R writer](https://github.com/robertaboukhalil/ginkgo/blob/master/scripts/process.R)
constructs START using previous boundary plus one and
END using bin boundaries. The current source is not proof of the author's exact
executed historical version. The actual table independently supplies the same
closed-interval adjacency: each included chromosome begins at 1 and every
subsequent START equals the preceding END plus one.

Collector 23051701 verifies the exact header and inventories all literal fields,
widths, adjacency and FAI bounds. Four independent-answer extent fixtures cover
zero, one, terminal and singleton cases. There are no duplicate feature IDs,
duplicate intervals, or extents incompatible with both assessed conventions.

Independent audit 23052386 reparses the full preserved archive using a separate
CSV reader and reproduces every row, identity, width histogram, adjacency and
extent summary. It exports `start0=START-1, end0=END`, preserving source row and
feature identifiers. Integer equality is required; no floating tolerance or
scientific cutoff decides these coordinate checks.

The first audit attempt, 23052098, stopped because it iterated Counter keys when
checking multiplicities. The correction inspects Counter values and retains the
same requirement that every interval occur once. Unique and duplicate fixtures
pass in the corrected run. Failed source, output and traceback remain unchanged.

Examples demonstrate why physical support matters:

| Published interval | START | END |
|---|---:|---:|
| First chromosome 1 bin |1|1191004|
| Second chromosome 1 bin |1191005|1776825|
| First chromosome 15 bin |1|20369064|
| Final included X bin |114248236|114820832|

These wide extents are not a uniform 500 kb resolution guarantee. The original
mappability/masking provenance is still needed to interpret physical support.

| Outside published-table support | Zero-based start | Exclusive end |
|---|---:|---:|
| X tail |114820832|156040895|
| Y, entirely absent |0|57227415|
| MT, entirely absent |0|16569|

The reason for omitted territory is unestablished. Do not interpret a missing
row as a neutral observation, deletion, zero count or lack of raw-read support.

## Execution, artifacts and decision

| Job | Stage | Seconds | Process peak RSS KiB | Exit |
|---|---|---:|---:|---:|
|23051701|coordinate inventory|17|90524|0|
|23052098|initial independent audit|10|89088|1|
|23052386|corrected independent audit|10|92400|0|

All jobs run on compute nodes in the existing isolated validated environment.
Evidence is under `dev/scone_dna_coordinates_20261109`: original archive and
coordinate inventory in `coordinates_v1`, successful independent conversion in
`coordinate_audit_v2`, and the retained failed attempt in `coordinate_audit_v1`.
Each has an immutable source snapshot, command/environment manifest, source
identities and input/output hashes. The new scripts are
`benchmarks/hardening/scone_dna_coordinates.py` and
`benchmarks/hardening/scone_dna_coordinates_audit.py`.

Decision: use the audited literal variable intervals for any subsequent
resolution-matched comparison and retain outside-support status explicitly.
The earlier fixed-bin raw DNA census remains a separate descriptive result;
it is not aligned to this published table. No count comparison, segmentation,
DNA truth, calibration, cohort admission or scientific default is promoted.
Original reference sequence, callable sensitivity and normalization remain open.
The biological null count stays zero. Original recovery, false-split/power and
uncertainty failures remain unresolved. Frozen results and criteria are unchanged.
The full package suite was not rerun for these diagnostic-only additions;
targeted compute-node fixtures and independent actual-table checks pass.

## Raw starts occupy omitted X territory

The next locked experiment allocates only the existing SRR16194948 pilot DNA
BAM to the actual published intervals. No new alignment, assay-value inspection,
MAPQ threshold, normalization or segmentation is introduced. These observations
are aligned primary-mate starts and leftmost proper-pair starts, separately;
neither measures physical base coverage or independent molecules.

| Unit | All observed | Inside published intervals | Outside |
|---|---:|---:|---:|
| Mapped primary-mate starts |218406|210746|7660|
| Leftmost proper same-chromosome pair starts |106289|102827|3462|

All 112,521 templates retain their two primary-mate records. There are 104,376
distinct coordinate-plus-candidate-UMI fingerprints among eligible pairs, exactly
matching the prior descriptive census; this does not reproduce original Ginkgo
deduplication. Twenty-four proper pairs have different mate-start targets. Their
two starts remain distinct observations rather than being silently coassigned.

The omitted X tail contains 1,648 mapped primary-mate starts, including 1,372 at
MAPQ 60, and 802 eligible proper-pair left starts. Every MAPQ stratum is retained.
All terminal outside-table starts occur on X. Other outside starts include Y,
MT and noncanonical contigs, each explicitly retained. A lexical, quality-blind
panel of ten mate records beyond the X boundary and ten within its final included
bin preserves sequence, quality, flags and CIGAR for bounded inspection. Existing
all-mate sequence/quality audits link the BAM to the prepared raw-read subset.

The competing assertion that the omitted X territory has no mapped observations
is falsified by these retained reads. Their presence does not prove reliable
dosage, unique locus assignment, copy number or why the published table omits
that territory. Missing rows must remain outside support rather than become
neutral observations or deletions. Original masks, historical preprocessing,
count units, normalization and assay sensitivity remain unresolved.

Count job 23053246 (nslab, CPU2/4GB/1h) completed in 13 seconds with process peak
264,204 KiB and exit 0. Independent audit 23053775 (innovation, CPU2/4GB/1h)
completed in 6 seconds with process peak 228,864 KiB and exit 0. Both ran on
compute nodes using the existing isolated wheel/native environment. Each passes
seven boundary/gap fixtures. The independent audit uses original one-based closed
coordinates and NumPy inclusive-END search, without importing the counter's
zero-based assignment. It exactly reconstructs all 225,042 mates, 112,521
templates, 5,393 interval rows, 33,866 quality strata, 20 panel records and all
reported integer summaries. No floating tolerance is relevant to these counts.

New scripts are `benchmarks/hardening/scone_dna_interval_counts.py` and
`benchmarks/hardening/scone_dna_interval_counts_audit.py`. Frozen designs,
source snapshots, outputs, input/output hashes and execution manifests are in
`dev/scone_dna_coordinates_20261109`, under `interval_counts_v1` and
`interval_count_audit_v1`. The first counter's frozen report correctly records
the audit as not yet run at its execution; the later audit is a separate linked
record. No original output is overwritten. This is a computational support
result, not a correction of the three recovery failures, uncertainty or the
false-split/power failure. No biological validation claim or default changes.

## Deposited truth construction carries explicit assumptions

Source-only compute job 23053989 completes in 3 seconds, process peak 87,552 KiB,
exit 0 (nslab, CPU1/2GB/1h). It verifies the previously metadata-checksummed
Zenodo archive SHA256 and retains exactly eight named DNA source/dependency
members. Each retained file matches its archive member byte for byte. No author
code executes, no network request occurs and no assay values are parsed.
`dna_author_source_inventory_v1/report.json` records all member hashes; the
corresponding design and immutable snapshot are retained alongside it.

The primary code distinguishes several layers that cannot be called interchangeable
DNA truth:

| Source, line | Deposited operation | Consequence for validation |
|---|---|---|
| `RawDataProcessing/2_DNA_BWA_mapping.sh`, 9–15 | Reverse mate order for first alignment; retain alignments with MAPQ at least 59 | Our unfiltered original-order diagnostic is not a reproduction of that workflow |
| Previously inspected `3_DNA_UMI_mapping.sh`, 23–46 | Select flag64 mate, UMI-tools deduplicate, recover its pair, realign and export MAPQ-filtered BED records | Coordinate-plus-UMI fingerprints are not the deposited deduplication/count unit |
| `DNA_related/process_segnorm.R`, 7–15 | Exclude Y, rescale by a density-mode estimate, segment, multiply by 2 and round | Integer CN includes a baseline assumption; it is not directly measured absolute DNA copy number |
| `DNA_related/mergeLevels_multi_segmentation.R`, 18–26 | Replace nonpositive values with 1e-7, median-normalize/log2, winsorize and jointly segment using multipcf | Published segments are processed measurements with shared estimation, not independent bin observations |
| `DNA_related/README.md`, integer-CNV sections | Describes diploid count-based scaling and explicitly restricts allele-based correction to its tetraploid sample | Do not port that correction to other baseline ploidies or treat it as general independent truth |
| `DNA_related/Transfer_CHISEL_mbaf_to_ginkgo.R`, 4–6 and 23–30 | Restrict to autosomes, find overlaps and keep only the first overlap per destination feature | Replicated coarse BAF does not establish finer-resolution allelic truth; multiple-overlap evidence is discarded |
| `DNA_related/extract_baf_from_chisel.R`, 6–18 | Use the first barcode's bin annotation and assign other cells' pulled values positionally | Exact coordinate/order correspondence must be checked before adopting these matrices |

[Inference] The allele-transfer operation cannot distinguish two allele profiles
that agree on the retained overlap and differ only on discarded overlaps. This
is an information limit, even if every downstream computation executes. Nominal
500 kb labels do not convert the stated 5 Mb CHISEL input into 500 kb independent
allele evidence. Our observed variable physical widths remain the relevant
destination coordinates.

The integer-CN script also emits messages rather than stopping for the listed
identity checks, and uses `all(names1 != names2)` rather than exact identity.
Consequently, the guard can miss a partial permutation. Its local `seg2cluster`
helper compares adjacent values without chromosome input and accumulates short
runs before emitting a group (lines 65–99). These are source-level hazards to
check before any reuse, not executed failures demonstrated for the pilot. The
archive is not proof of the actual historical code, parameters or intermediates.
No deposited source is edited and no published measurement is invalidated solely
by this inspection. None of these files explains the exact terminal X omission;
its origin remains open.

Decision: retain published normalized depth as an uncertain relative measurement,
integer CN as model-dependent, and coarse allele evidence at its original
resolution. Do not admit it as absolute-CN or fine-scale allele truth. Before a
comparison, require the original Ginkgo configuration/bin-mask identities,
matching count units and depth/normalization provenance. For allele validation,
require explicit keyed cell/coordinate joins, all overlaps with declared weights,
and original coarse support rather than first-overlap pseudo-resolution. The
smallest decisive checks are a partially permuted ID fixture, a destination bin
straddling two discordant source bins, and a same-valued chromosome boundary.
If the historical executed workflow and suitable intermediates remain missing,
do not substitute those reconstructed examples for biological validation.

### Minimal source counterexamples executed

Job 23054448 completes on an innovation CPU node in 7 seconds, process peak
77,568 KiB, exit 0. The installed module supplies R 4.3.3; its executable path,
hash, version and loaded module inventory are retained. No packages are installed.
The job verifies the source-inventory hashes, parses the deposited integer-CN
script without sourcing it, and executes only its first identity guard and
isolated local grouping helper on tiny synthetic matrices.

The original guard produces no message for sample labels `(a,b,c)` versus
`(a,c,b)`. For an entirely different label vector it emits a message and continues.
Thus the stated source guard neither detects every permutation nor enforces
failure. These are executed counterexamples to the guard's sufficiency; whether
the pilot ever encountered such inputs remains unverified.

For twelve equal values spanning two declared six-bin chromosomes, the grouping
helper returns one group. Its interface has no chromosome input. This establishes
that chromosome-contained runs cannot be reconstructed from that interface
alone. It does not establish that the author intended chromosome-contained
grouping, or that any actual pilot inference is wrong.

An independent inclusive-coordinate calculation uses a destination `[4,7]` and
two source bins `[1,5]`, `[6,10]` with values 0.1 and 0.4. Exact overlap lengths
are 2 and 2. First-overlap projection returns 0.1 or 0.4 when source order changes;
a change confined to the discarded source leaves its output unchanged. The
length-weighted descriptive average is 0.25. That average is a reference
calculation, not a validated allele model or proposed scientific default. The
actual author GenomicRanges transfer function is not executed in this test.

The design, source snapshot, bounded R program, stdout/stderr, report and execution
manifests are in `source_counterexamples_v1` and its linked source directory.
Script: `benchmarks/hardening/scone_dna_source_counterexamples.py`; module launcher:
`scripts/scone_r_source_readonly.sbatch`. Original deposited files are unchanged.
The source-inventory job's earlier no-execution statement remains accurate for
that job; this separate protocol explicitly permits only these bounded calls.
No assay values, historical workflow reconstruction, DNA truth or clone accuracy
claim is added. These checks justify requiring explicit identity and resolution
validation before using these deposited outputs as a reference.
