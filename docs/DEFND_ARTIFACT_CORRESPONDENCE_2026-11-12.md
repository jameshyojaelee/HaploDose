# DEFND artifact correspondence and a metadata schema correction

Two existing prospective leads now have exact provider-declared RNA/DNA file
joins. A private parser omission is corrected in a separate cached replay,with
all prior attributes preserved. The metadata also exposes a consequential
leakage hazard: RNA column headers contain author cluster numbers. No measurements
or cluster values are opened,no inference is run,and neither cohort is admitted.
This is preparation for independent performance evaluation,not an accuracy gain.

## Scientific contract and fixed selection

The estimand is archive artifact correspondence. Metadata are declarations,not
physical-cell or genotype measurements. A matching accession could still link
the wrong layer,barcode convention or assembly. Artifact-specific reference and
execution provenance plus a library-qualified barcode crosswalk would distinguish
those alternatives. DNA sensitivity and biological replication remain separate.

The lexically first diagnostic accession is GSM7016419; GSM7817780 is retained
in the same table. Both were identified in the November8 prospective plan,
before any HaploDose outcomes. This follows that existing gate rather than
adding another discovery cohort. No candidate,comparator,endpoint or scientific
threshold is changed. Former evaluations remain regression-only.

## Exact new archive evidence

The public [Zenodo record13984673](https://zenodo.org/records/13984673) declares
the DNA fragments. The complete [GEO family SOFT metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE224nnn/GSE224149/soft/GSE224149_family.soft.gz)
supplies the RNA links and processing declarations. No linked measurement is fetched.

| Accession | RNA filename | DNA filename | DNA catalogue bytes / MD5 |
|---|---|---|---|
| GSM7016419 | GSM7016419_TB5471.matrix.txt.gz | GSM7016419_TB5471.fragments.tsv.bz2 | 13318383087 / 5bd43ff77caeb21189651502d1ea9f7b |
| GSM7817780 | GSM7817780_PTO037.matrix.txt.gz | GSM7817780_PTO037.fragments.tsv.bz2 | 6838939966 / 956978248393b529ec82caa960ec352d |

These are exact matching accession/patient filename stems and provider catalogue
identities. Checksums of undownloaded DNA are not verified. Both RNA entries
declare a gene-by-cell count matrix whose column names include barcode and
cluster number; human-only processing declares GRCh38/Gencodev24. A declaration
is not an artifact-specific reference fingerprint or proof of the exact executed
source. No independent biological n is qualified from filenames alone.

Author cluster labels must be excluded before any future preprocessing,reference
selection,discovery or truth construction. A future frozen reader must extract
library-qualified barcodes without exposing or using those labels. No actual
headers or cluster values are read here. Published clusters cannot substitute
for independently constructed DNA profiles or memberships.

## Every attempt and the exact correction

Normal GEO web pages return a verification page. No challenge is automated or
bypassed. The web tool also cannot retrieve the Zenodo API or public FTP index.
The compute job uses the exact standard public metadata-distribution endpoints,
ordinary requests and redirects restricted to the archive hosts. It requests no
controlled data,credentials or measurement download. Both metadata requests
succeed within the frozen1/4MiB response and16MiB inflation bounds. Their full
responses,headers and SHA256 identities remain.

**23238269 completes0:0** with10 exact before/after source identities. Original
synthetic fixtures verify repeated attributes,target isolation,duplicate rejection,
missing targets and internal accession mismatch. Zenodo artifacts and both GEO
sample records resolve. The parser retains only an unnumbered optional
supplementary-file key; actual GEO uses supplementary_file_1,so RNA file links
are omitted. Its first success report has this limited scope and remains unchanged;
it admits no cohort and does not prove complete artifact correspondence.

**23238318 completes0:0** in a separate schema version,using only the cached
responses. The correction normalizes precisely numbered supplementary-file
keys before calling the original parser. New fixtures retain numbered and
unnumbered fields. An independent original-text scan reconstructs both ordered
file lists. All previous retained attributes and cached response hashes remain
exact; all23 source/cached identities match before/after. Original source,
responses and results remain untouched. No second network request or criterion
relaxation follows the repair.

| Job | Progress | State / exit | CPUs / memory / limit | Runtime / process peak |
|---|---|---|---|---|
| 23238269 metadata | 2/2 ██████████ | Complete0:0 | 1 / 2GiB / 30min | Slurm2s,timed2.12s / 20280KiB |
| 23238318 cached schema | 2/2 ██████████ | Complete0:0 | 1 / 1GiB / 10min | Slurm1s,timed0.99s / 18496KiB |

Both use nslab account/QOS and CPU,existing isolated read-only Conda/runtime,
one thread and fresh neutral outputs. Short-job Slurm RSS0 is an unavailable
sample,not zero memory. Zero running/pending jobs and zero failed Slurm jobs in
this chain; the optional-field omission is preserved and corrected separately.
No arrays,installs,shared changes,running-job modifications or public-suite rerun.

## Admission decision and next decisive gate

Provider file correspondence passes. Barcode identity,DNA reference/command
provenance,fragment multiplicity convention and DNA sensitivity remain unverified.
The assay's accessibility bias also needs independent controls. **Do not admit
genomic truth or compute accuracy from these metadata.** No new reference,cells,
clone labels or biological validation units are established here.

Before measurement access,freeze the candidate/comparator/endpoint design and
a reader that excludes embedded author labels. Then the existing bounded
header/coordinate/barcode gate can test unique library-qualified identities
and declared assembly/multiplicity against the relevant artifact provenance.
A contradiction must stop the path; absent provenance remains explicit. No
full20GB staging,clustering,method tuning or resolution selection is justified
by the filename match alone. A two-patient lead cannot replace59 independent
single-CNA null specimens. Original recovery2/5,state errors14/200,power195/200,
coverage54/1000,fixed competitors and frozen artifacts remain unchanged.

[Original fixed design](../dev/defnd_artifact_metadata_20261112/DESIGN.md),
[first driver](../dev/defnd_artifact_metadata_20261112/run.py),
[first complete report](../dev/defnd_artifact_metadata_20261112/run_23238269/results/REPORT.json),
[schema correction contract](../dev/defnd_artifact_metadata_20261112/SCHEMA_V2.md),
[cached schema driver](../dev/defnd_artifact_metadata_20261112/schema_v2.py),
[corrected joins and all checks](../dev/defnd_artifact_metadata_20261112/schema_v2_23238318/results/REPORT.json),
[original source inventory](../dev/defnd_artifact_metadata_20261112/SNAPSHOT.json),
[schema source inventory](../dev/defnd_artifact_metadata_20261112/SNAPSHOT_V2.json),
[first submission](../dev/defnd_artifact_metadata_20261112/SUBMISSION.json),
[schema submission](../dev/defnd_artifact_metadata_20261112/SUBMISSION_V2.json).
