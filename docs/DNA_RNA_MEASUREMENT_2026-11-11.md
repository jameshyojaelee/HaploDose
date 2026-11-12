# DNA/RNA development measurements:availability and identity resolved in part

The scientific goal remains accurate genomic gain/loss profiles and supported
CNA clones from RNA. This investigation checks whether new DNA/RNA development
measurements could supply information missing from RNA-only fitting. It finds
usable metadata correspondence and consequential annotation ambiguities;it does
not establish a new performance improvement,clear a training cohort or validate
genomic membership. Existing recovery,false-split,power and uncertainty criteria
remain unchanged.

## What the source measurements are

The official24Q4 README,download manifest and release-specific notes distinguish:

| File | Measurement and processing | Restriction |
|---|---|---|
| OmicsExpressionProteinCodingGenesTPMLogp1.csv | RSEM unstranded log2(TPM+1),ModelID rows | Transformed abundance,not raw-count uncertainty;ModelID presence alone does not identify the contributing ProfileID |
| OmicsExpressionRNASeQCGeneCountProfile.csv | RNA-SeQC2 read counts,processed in unstranded mode,Gencode v29,ProfileID rows | Counts are reads,not UMIs;library protocol and processing mode remain separate |
| OmicsCNGene.csv | Linear relative DNA copy ratios from WGS/WES segments;BioMart November2020 to ENSG,then Gene.csv to identifiers | Not absolute copies or log2(CN+1);mapping and DNA masking limit interpretation |
| OmicsCNSegmentsProfile.csv | ProfileID-indexed relative copy segments | Coordinate conventions,callability and sensitivity still require qualification |
| OmicsAbsoluteCNSegmentsProfile.csv | PureCN absolute segment estimates | Uncertain model-based reference,not ABSOLUTE-method output or certified clone truth |

The generic current GitHub README retains a legacy transformation note,so the
release-specific README and notes govern these inputs. Four CCLE benchmark files
now correctly name PureCN in comments/docstrings and future scorer provenance.
Numerical ASTs match after removing only module documentation and the declared
provenance string;executable shell lines match exactly. No frozen score,truth
file,threshold or historical metadata is rewritten.

[Official24Q4 release](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248),
[24Q2 scale and PureCN changes](https://forum.depmap.org/t/announcing-the-24q2-release/3312),
[24Q4 identifier and method changes](https://forum.depmap.org/t/announcing-the-24q4-release/3564),
[downloaded release README](../dev/dosage_response_measurement_20261111/diagnose_v2_23224672/results/README.txt),
[exact attribution-only comparison](../dev/dosage_response_measurement_20261111/attribution_parity.json).

## Fixed candidates and measurement availability

All182 candidates from the preceding exposure/PatientID exclusions remain
provisional. All have model-level TPM/DNA rows,default RNA count-profile rows,
and default relative-DNA segment profiles. Their default RNA/DNA ModelCondition
IDs agree for171;the other11 are explicitly retained in the inventory rather
than silently treated as corresponding material. Identical metadata IDs do not
prove same aliquots or same-cell truth. ModelCondition and condition-default
tables are staged but their full consistency has not yet been assessed.

DNA assay metadata comprise122WGS and60WES profiles. RNA library metadata comprise
83unstranded and99stranded profiles;the available counts were processed with
the unstranded operator. No candidate has its default DNA ProfileID in the
absolute-segment file. This follows the conservative earlier prepared-truth
exclusions and does not establish that the candidates are high-ploidy,neutral
or otherwise biologically homogeneous. Relative-DNA association is the possible
target;an absolute-DNA validation claim has no reference here.

All five original measurement files match their manifest MD5 identities;SHA256
identities are recorded. File streaming decodes only headers and identifier
prefixes. No expression value,copy-ratio value,segment coordinate or count
payload is converted or used for selection. HepG2/C3A resolves through its
complete row's ModelID/RRID and was already excluded. Nine names remain
unresolved;no identities are guessed and no newly matched model is admitted.

## Annotation failures and restricted identity projection

The strict original admission check correctly rejects repeated NCBI feature IDs.
The complete inventory shows:

| Matrix | Rows | Columns | Consequential header issue |
|---|---:|---:|---|
| RSEM TPM |1673|19193|17NCBI IDs each associated with two distinct labels |
| Relative DNA |1929|38590|Stable feature IDs unique;229symbol disagreements with unique common RNA IDs |
| RNA-SeQC counts |1690|57975|44canonical/suffixed Ensembl duplicate pairs and92ERCC external controls |

**Correction to the early diagnostic description:**136labels failed its canonical
gene-label pattern. These comprise44suffixed Ensembl labels plus92ERCC labels.
The frozen second design's statement that all136 were suffixed is incorrect;
its source and failure remain intact. The final namespace-aware inventory and
this correction supersede that interpretation.

The diagnostic crosswalk requires one canonical header column and one Approved
HGNC row per Ensembl/NCBI ID,an autosomal protein-coding locus,and exact release
Gene.csv correspondence. Every alias of an ambiguous ID is quarantined.
Neither choosing the first column,summing duplicates nor guessing a symbol is
allowed. External controls and unknown features cannot enter genomic joins.
All exclusions remain visible. The resulting identity-only projection contains
17766count/DNA,17633TPM/DNA and17623three-way genes. These are identity-qualified
features,not a validated DNA-callable mask or an approved analysis input.

Nine hand-calculable checks pass:all-alias quarantine,suffixed-pair quarantine,
literal-duplicate rejection,permutation identity,cross-namespace column alignment,
ambiguous HGNC rejection,no numerical-payload decoding,control exclusion without
index shifts,and unknown-label non-admission. The original strict rejection
remains;successful diagnostic completion is not successful admission.

## Executed attempts and resources

Every job uses account/QOSnslab,CPU1CPU4GB1h,the existing isolated Conda/runtime,
one BLAS thread and a new output directory. No shared installation changes.

| Job | Result | Scheduler elapsed | Analysis peak process RSS,KiB |
|---|---|---:|---:|
|23224462|Failed1:0;strict duplicated-stable-ID gate|2s|36864|
|23224633|Failed1:0;diagnostic parser encounters suffixed label;four release metadata files staged and verified|5s|54664|
|23224672|Complete0:0;header-only inventory,including noncanonical labels|1s|60828|
|23224790|Failed1:0 after seven fixtures pass;ERCC namespace not yet represented|7s|89620|
|23224819|Complete0:0;nine fixtures and full metadata inventory|8s|235420|

All three failed attempts remain as input/parser evidence;none is erased.
The final full inventory does not change the scientific decision to withhold
admission. Scheduler sampling misses these short-process peaks;process resource
logs and scheduler records are kept separately. There are no pending downstream
jobs. No package suite is rerun because protected production/test source remains
unchanged. Benchmark attribution edits receive static AST/shell checks.

[Original design and strict gate](../dev/dosage_response_measurement_20261111/DESIGN.md),
[preserved original source](../dev/dosage_response_measurement_20261111/SNAPSHOT.json),
[failed strict check](../dev/dosage_response_measurement_20261111/run_23224462/output.txt),
[header diagnostic result](../dev/dosage_response_measurement_20261111/diagnose_v2_23224672/results/REPORT.json),
[restricted projection design](../dev/dosage_response_measurement_20261111/DESIGN_v2.md),
[control/unknown-label correction](../dev/dosage_response_measurement_20261111/DESIGN_v3.md),
[final executable inventory](../dev/dosage_response_measurement_20261111/metadata_v3.py),
[hand-calculable fixtures](../dev/dosage_response_measurement_20261111/qualify_v3.py),
[executed fixture results](../dev/dosage_response_measurement_20261111/run_v3_23224819/qualification.json),
[complete result](../dev/dosage_response_measurement_20261111/run_v3_23224819/results/REPORT.json),
[all candidate correspondence and availability](../dev/dosage_response_measurement_20261111/run_v3_23224819/results/measurement_metadata.tsv),
[identity-qualified feature projection](../dev/dosage_response_measurement_20261111/run_v3_23224819/results/identity_qualified_features.tsv),
[final execution snapshot](../dev/dosage_response_measurement_20261111/SNAPSHOT_v3.json),
[all job accounting](../dev/dosage_response_measurement_20261111/accounting_complete.tsv).

## Scientific decision and next boundary

Do not train or allocate prospective roles from this inventory. First finish
condition/identity and historical-exposure qualification,coordinate/gene mapping,
and the resolution/callability assessment of the relative-DNA reference. The
available processed reference cannot establish its own biological accuracy.
Any later RNA-DNA learning must treat it as an uncertain measurement,separate
development from evaluation before values are opened,and test within-lineage
variation rather than learning lineage-CNA associations uncritically.

Regional transcription,culture differences and assay processing remain competing
explanations of RNA-DNA association. Corresponding independent DNA or replication,
matched preprocessing and donor-heldout prediction must distinguish them. There
is no new same-cell/spot membership truth,no59-independent-specimen single-CNA
panel,and no absolute-DNA reference in these candidates. Training on this panel
could not by itself validate clones or establish biological false-split control.
Former evaluation specimens remain regression-only. The overall genomic/SOTA
goal and every prior negative result remain open.

The final preservation check confirms91protected source/data/test identities,
all368prior ledger values,five immutable execution snapshots,and the four
attribution-only benchmark comparisons. Seven Python sources parse,five new
launchers pass shell syntax,and833local report links resolve. Both successful
jobs and all three failed attempts are terminal;no dependent analysis remains.
The original521-Python/four-Rust validation is retained as historical evidence,
not rerun here. [Exact preservation record](../dev/dosage_response_measurement_20261111/preservation_checkpoint_v1.json).
