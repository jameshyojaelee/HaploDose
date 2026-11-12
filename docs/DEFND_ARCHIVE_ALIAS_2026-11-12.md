# Independent archive metadata do not resolve the DEFND alias

The preceding bounded inspection made progress by finding an actual specimen
name discrepancy,not by improving performance. This follow-up checks the
already-linked RNA experiment/sample against independent ENA metadata. It does
not resolve TB5471 versus the DNA header's PTO026. The existing async request
for an original local manifest remains pending;no repeat question was sent.

The estimand is accession-qualified naming correspondence. A harmless processing
alias and a wrong specimen remain competing explanations. Explicit identity
metadata can support an alias;absence here cannot prove a mix-up. Neither archive
accession agreement nor filename similarity establishes same-cell DNA/RNA truth.

The first public ENA file report returns402 bytes and only file-catalogue
columns. Assuming this endpoint supplied all record metadata was incorrect;
the driver rejects the incomplete identity schema and exits1:0. Preserve that
request,body,SHA256 b8435e3769eac427674aca7e0c187bd50c2536fe1fd13aa8b547e8273581fe68
and failure. No FASTQ link is followed. The report identifies SRR23292060 but
does not supply an explicit PTO026/TB5471 alias. Do not weaken identity checks.

A separate driver uses the official
[ENA Browser API XML interface](https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html)
for exactly SRX19235312 and SRS16639724,one request each,1MiB cap/45s socket
timeout,HTTPS/same-host redirects. Exact accession and sample cross-links pass.
The experiment alias is GSM7016419_r1,library_name GSM7016419;its descriptor
points to SRS16639724. The sample primary accession SAMN32979122 has secondary
SRS16639724 and submitter alias GSM7016419. Neither record's inspected identity
fields contains PTO026 or TB5471. Raw XML/headers/hashes are preserved. No raw
reads,genomic profiles,cell memberships or biological outcomes were opened.

| Job | Terminal progress | State/exit | Node | Slurm/timed runtime | Process peak KiB |
|---|---|---|---|---|---:|
|23238847 file-report check|1/1 ████████|FAILED1:0,preserved|ne1dc3-005|3s/2.65s|19688|
|23238869 linked XML check|1/1 ████████|COMPLETED0:0|ne1dc4-046|38s/37.58s|20172|

Both jobs use nslab account AND QOS,cpu partition,1CPU/1GiB/10min,no array,
existing read-only isolated Conda/runtime,neutral outputs and Python-I. Exact
source checks10/24 pass before/after. File-report fixtures reject four malformed
correspondences;XML fixtures independently check expected alias/library fields,
description exclusion and three invalid record/link cases. Static AST and bash-n
checks pass. Public pytest/Rust suites are not rerun because public source did
not change. Running0,pending0,completed1,failed1;downstream inference/evaluation
remains unsubmitted. Unrelated running/pending jobs remain untouched.

[Initial design](../dev/defnd_archive_alias_20261112/DESIGN.md),
[source-supported XML design](../dev/defnd_archive_alias_20261112/XML_DESIGN.md),
[complete actual results/resources](../dev/defnd_archive_alias_20261112/CHECKPOINT.json),
[first failed report](../dev/defnd_archive_alias_20261112/run_23238847/results/REPORT.json),
[actual XML result](../dev/defnd_archive_alias_20261112/xml_23238869/results/REPORT.json).

This closes the additional public-metadata check without supplying the required
alias or executed barcode history. No new performance gain,cohort admission,
DNA truth/resolution qualification,biological replication,production/default
change or scientific-criterion pass follows. Frozen results remain unchanged.
The5/5 research recovery candidate still fails its joint error/power gate;
production recovery2/5,state errors14/200 and selected coverage54/1000 remain.
The full genomic-profile,clone,spatial and independent-superiority goal is not
achieved. Another run of a rejected model would not resolve this evidence gap.
