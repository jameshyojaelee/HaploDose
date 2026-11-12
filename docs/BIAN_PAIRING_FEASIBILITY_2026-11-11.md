# Bian same-cell correspondence and measurement admission

The scientific goal requires cell-level genomic evidence, not pooled RNA/DNA
agreement alone. The preceding [exact-path experiment](EXACT_PATH_MOVES_2026-11-11.md)
rejects its proposed fitting correction. This metadata-only check returns to an
existing candidate source before opening any measurement or fitting another model.
It preserves Bian's historical secondary/descriptive role and candidate-freeze
requirement. It does not admit an independent evaluation panel.

## Primary measurement constraints

The [GEO RNA sample protocol](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2696969)
distinguishes Tang RNA libraries for CRC01/CRC02 from multiplexed UMI libraries
for the other patients. The processed files contain FPKM or TPM; their RefSeq
gene summaries and reported hg19 build do not establish raw-count provenance or
the actual executed preprocessing. Earlier catalogue text calling the entire
cohort full-length needs this correction. Do not convert normalized expression
to integer pseudo-counts for a calibration claim.

The existing CRC source
`/gpfs/commons/home/jameslee/CRC/src/49_bian_fragment_profiles/readout.py`
concerns one patient,CRC01. Other project source files perform DNA-group/RNA
expression joins for additional patients. Those source descriptions are not an
unopened HaploDose benchmark, independent genomic truth qualification or proof
of ten usable paired specimens. No measurement values from those outputs are
read in this experiment. All other-project files and jobs remain untouched.

## Exact feasibility gate

Parse only the existing GEO filename/size catalogue and local asset presence.
Keep all clinical RNA/methylome files, exact cell stems, declared RNA units,
missing layers and ambiguous assets. A filename match is candidate correspondence,
not certified same-cell identity. RNA/DNA specimen joins still require provenance;
DNA truth needs independent reconstruction, coverage and resolution qualification.

An independent prefix/suffix parser must reconstruct every candidate join
exactly. Fixtures cover missing layers, duplicate assets, mixed RNA units,
malformed names, nonclinical controls and layer/unit contradictions. Normalized
RNA fails the raw-count admission gate. Ambiguous pairs abstain. A descriptive
profile test requires a separate frozen method/reference/comparator protocol;
raw reprocessing or a new scientific claim does not follow automatically.

Biological n analyzed=0. Patient/cell counts from this exercise are metadata
counts, not biological validation replicates. Cells, regions and libraries cannot
fill the missing59-independent-single-CNA-null requirement. Preserve every
original recovery, error/power, uncertainty and all-comparator criterion.

## Execution

Job23230680 requests nslab/nslab CPU,1CPU1GB1h, constraintv6&amd. Existing
read-only isolated Conda/runtime, neutral cwd, one thread and99 before/after
source identities. Static Python AST and shell syntax pass. The program imports
only the standard library, does not open a measurement matrix or decompress
an expression/DNA file, and records no genomic inference. No downloads, arrays,
new dependencies, shared installations or changes to frozen results. Outcome pending.

[Frozen design](../dev/bian_pairing_feasibility_20261111/DESIGN.md),
[metadata implementation and independent parser](../dev/bian_pairing_feasibility_20261111/study.py),
[input and primary-source record](../dev/bian_pairing_feasibility_20261111/INPUTS.json),
[source identities](../dev/bian_pairing_feasibility_20261111/SNAPSHOT.json),
[submission](../dev/bian_pairing_feasibility_20261111/submission.json).

## Completed filename gate and explicit quarantine

Jobs23230680/23230716 complete0:0. Two independent canonical parsers agree;
all missing-layer, duplicate, mixed-unit, malformed-name and layer/unit fixtures
pass. Of2504 catalogue records,2464 have canonical clinical filenames,19 have
noncanonical clinical RNA filenames,20 describe nonclinical cell-line controls
and one is the archive. All2504 are accounted for. The second job retains the
19 records rather than treating them as nonclinical or silently renaming them.

| Canonical patient identifier | Candidate pairs | RNA-only | DNA-only | Canonical RNA unit |
|---|---:|---:|---:|---|
|CRC01|241|125|168|FPKM|
|CRC02|51|83|8|FPKM|
|CRC03|0|131|0|TPM|
|CRC04|57|34|36|TPM|
|CRC06|0|71|0|TPM|
|CRC09|14|85|13|TPM|
|CRC10|70|15|54|TPM|
|CRC11|176|16|78|TPM|
|CRC12|0|0|23|None|
|CRC13|0|0|210|None|
|CRC14|0|0|36|None|
|CRC15|0|0|60|None|

There are609 canonical candidate pairs across six patient identifiers,560
RNA-only cells and686 DNA-only cells. These are filename-derived candidates,
not six certified independent specimens. The twelve identifiers and historical
ten-patient description are not interchangeable sample-size estimates. All609
paired RNA assets and all19 noncanonical RNA assets are present locally.
No paired DNA asset is in the expression directory; this does **not** mean DNA
is absent elsewhere in the CRC project. Those other locations are not inventoried
by this filename check. No RNA or DNA measurement is read.

The19 noncanonical names repeat the `scTrioSeq2Rna` prefix and contain NC-coded
stems. They remain quarantined for identity and reference selection. Bounded
GEO page checks onGSM3443113/GSM3443126 returned tool access errors and
GSM3443131 returned a browser-verification page. No author-validated aliases
were established. No filename-normalization heuristic may select them as normals.
NC-coded identity also does not establish biological diploidy.

[Canonical complete report](../dev/bian_pairing_feasibility_20261111/run_23230680/results/REPORT.json),
[complete accounting and quarantine](../dev/bian_pairing_feasibility_20261111/qualification_23230716/results/REPORT.json),
[qualification source](../dev/bian_pairing_feasibility_20261111/qualification.py),
[qualification source identities](../dev/bian_pairing_feasibility_20261111/qualification_snapshot.json),
[qualification submission](../dev/bian_pairing_feasibility_20261111/qualification_submission.json).

The original admission rule fails for inferential raw-count execution: every
catalogued clinical RNA input is normalized FPKM/TPM. No independent cohort or
descriptive measurement analysis is admitted. Original candidate-freeze and
secondary/descriptive requirements remain. More metadata does not establish
the biological5% false-split guarantee or resolve the existing power failure.

## Bounded official raw-library metadata follow-up

The archived series metadata links GSE97693 toPRJNA382695/SRP103839.
The [official ENA file-report API](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/file-reports.html)
supports study run metadata. One fixed request23230763 fetches the API field
inventory and the study's run/sample/title/alias/layout/source/FASTQ-metadata
columns. Targets are all19 quarantined GEO accessions and RNA/DNA accessions
for the first lexical paired cell in each protocol,CRC01_LN1_162 and
CRC04_LN1_349. Exact GSM tokens only; no measurement or FASTQ download.
Maximum20MiB per response and60s timeout; no automatic alias resolution,
count construction, candidate tuning or cohort admission.

The job requests nslab/nslab CPU1CPU1GB1h, with109 before/after identities.
Static AST/shell checks pass. Outcome pending; missing metadata remains a gap.
[Frozen raw-metadata design](../dev/bian_pairing_feasibility_20261111/RAW_METADATA_DESIGN.md),
[request implementation](../dev/bian_pairing_feasibility_20261111/raw_metadata.py),
[source identities](../dev/bian_pairing_feasibility_20261111/raw_metadata_snapshot.json),
[submission](../dev/bian_pairing_feasibility_20261111/raw_metadata_submission.json).

## Completed archive verification

Official metadata job23230763 completes0:0 and retrieves2492 run records.
All23 prespecified targets match exact GSM accession tokens; every requested
title/alias field is available. The two metadata responses are19766 and1104575
bytes, with saved request URLs and SHA256 identities. No FASTQ is downloaded.

Crosswalk job23230821 completes0:0. Every one of the609 canonical RNA/DNA pairs
matches exact archive accessions, author-supplied cell titles, RNA versus
bisulfite strategy and transcriptomic versus genomic source. Independent exact
accession joins agree and wrong-cell/accession/assay/missing-metadata fixtures
pass. This advances library correspondence, not same-cell integrity or truth
sensitivity. No clinical/genomic call is made.

All19 repeated-prefix RNA filenames now have explicit author-supported identities.
The original quarantine and failed GEO-page attempts remain preserved. The new
identity map supersedes the file-identity uncertainty only: none of these cells
is automatically admitted as a normal reference, and no original asset is renamed.
For example,GSM3443113's archive title specifiesCRC03_NC_292 and RNA runSRR8096019.
The two lexical pairs have separate RNA and bisulfite raw-library accessions:

| Cell | RNA run | Bisulfite run |
|---|---|---|
|CRC01_LN1_162|SRR5824391|SRR5825324|
|CRC04_LN1_349|SRR7461687|SRR7214590|

Raw-library metadata availability does not establish retained UMIs, read/molecule
deduplication, fragment independence or assay sensitivity. Full count construction
and gene/mapping contracts remain unverified. Reported protocol/build is not an
executed alignment inventory. Do not convert deposited FPKM/TPM into raw counts.

[Official metadata and exact requests](../dev/bian_pairing_feasibility_20261111/raw_metadata_23230763/results/REPORT.json),
[all609 archive-qualified candidate pairs](../dev/bian_pairing_feasibility_20261111/crosswalk_23230821/results/REPORT.json),
[explicit19-file identity map](../dev/bian_pairing_feasibility_20261111/crosswalk_23230821/results/FILE_IDENTITY_MAP.tsv),
[crosswalk source](../dev/bian_pairing_feasibility_20261111/crosswalk.py),
[frozen crosswalk design](../dev/bian_pairing_feasibility_20261111/CROSSWALK_DESIGN.md),
[crosswalk submission](../dev/bian_pairing_feasibility_20261111/crosswalk_submission.json).

Final metadata job23230843 includes every supported alias in the correspondence
table, so any additional paired libraries are retained rather than excluded by
the canonical filename rule. It checks the resulting pairs against the same
archive source and independently verifies the pair set by RNA/DNA set intersection.
Same nslab CPU1CPU1GB1h,120 source identities, no measurement access or download.
Outcome pending. The existing freeze, descriptive-role and scientific gates remain.
[Final correspondence code](../dev/bian_pairing_feasibility_20261111/finalize.py),
[final identities](../dev/bian_pairing_feasibility_20261111/final_snapshot.json),
[final submission](../dev/bian_pairing_feasibility_20261111/final_submission.json).

## Final decision and preserved scientific limits

Final23230843 completes0:0. The full correspondence table retains1874 clinical
cell stems, all19 supported aliases and all609 candidate RNA/DNA pairs. Every
candidate pair passes the exact archive accession/title/assay check; independent
RNA/DNA set intersection agrees. The19 aliases add no DNA-paired cells. There
are no unresolved metadata joins in this fixed candidate pair set. This is a
library-correspondence milestone, not609 validated genomic memberships or six
independent biological validation replicates.

[Final complete correspondence](../dev/bian_pairing_feasibility_20261111/final_23230843/results/CANDIDATE_CORRESPONDENCE.tsv),
[final authoritative report](../dev/bian_pairing_feasibility_20261111/final_23230843/results/REPORT.json).

| Job | Work | Completion | Process wall time | Peak process RSS |
|---|---|---|---|---|
|23230680|Canonical filename join|1/1 ████████;0:0|1.05s|20960KiB|
|23230716|Noncanonical clinical accounting|1/1 ████████;0:0|0.35s|18776KiB|
|23230763|Official raw-library metadata|1/1 ████████;0:0|8.92s|31144KiB|
|23230821|Author title/assay crosswalk|1/1 ████████;0:0|0.08s|26944KiB|
|23230843|Complete alias-aware correspondence|1/1 ████████;0:0|0.10s|28528KiB|

All jobs use nslab account/QOS CPU1CPU1GB1h. Scheduler runtimes are2s,0s,10s,
1s,0s; zero RSS reports reflect unsampled short jobs, not zero memory use.
Source checks pass before/after at99/104/109/115/120 identities respectively.
Running0,pending0,done5,errors0; no inference or downstream sequencing job.
The91 protected package/test/native-source identities remain unchanged. The
prior521-Python/four-Rust suite is not rerun for metadata-only changes.

**Passed:** filename accounting; exact source-identified library correspondence;
explicit resolution of the19 file identities; assay-scale distinction; all
independent metadata joins and rejection fixtures.

**Still unverified:** executed raw-RNA preprocessing and retained molecule tags;
gene/coordinate contracts; specimen integrity beyond archive titles; usable DNA
coverage, independent truth reconstruction and resolution; biological normality;
genomic clone labels and suitable independent replication. Count-model execution
on FPKM/TPM remains inadmissible. The archival raw-run locations identify a path
to a bounded measurement investigation, not permission to reprocess the cohort
automatically or bypass candidate freezing.

The next measurement work must first freeze an appropriate candidate and
reference/comparator/scorer contract, then inspect the implicated RNA-library
provenance and corresponding DNA sensitivity. Keep the source secondary and
descriptive. Measurement analysis remains subject to the existing candidate-freeze
requirement. Preserve all original recovery,
false-split/power, coverage, all-competitor and biological-validation criteria.
Public recovery2/5, state errors14/200, power195/200 and selected correlated
coverage54/1000 remain unchanged. No performance upgrade or SOTA claim is made.
