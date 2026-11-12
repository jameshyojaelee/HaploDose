# Independent genomic-response development:metadata candidates only

The scientific goal remains genomic gain/loss profiles and supported CNA clones
from RNA,with independent evidence at the claimed resolution. The corrected
[ordered RNA model](ORDERED_AMPLITUDE_JOINT_2026-11-11.md) repairs high-depth
optimization/recovery but retains seven low-depth failures. Earlier normal-RNA
prediction gains can worsen DNA compatibility. Before adding another RNA-only
latent model,check whether disjoint paired DNA/RNA observations could support
learning genomic response. No such model is trained or cleared here.

The completed [measurement inventory](DNA_RNA_MEASUREMENT_2026-11-11.md) now
verifies actual row/profile presence and original file identities. It finds171
matching default ModelCondition pairs,zero absolute-segment profiles,and an
explicit identity-only autosomal gene crosswalk. It also identifies17repeated
NCBI IDs,44Ensembl duplicate pairs,92ERCC controls and229symbol disagreements.
Numerical payloads remain undecoded. Candidates are still not cleared for
training,independent validation or biological error-control claims.

## Scientific contract

If justified later,the estimand would be predictive RNA-DNA association at
explicit regional resolution. It would not establish a causal dosage
coefficient,absolute copy number,genomic membership probability or biological
false-split guarantee. Lineage,transcriptional state,culture differences and
shared preprocessing are strong alternatives to a dosage mechanism.
Independent paired measurements and donor-heldout prediction must distinguish
useful genomic information from those alternatives. Metadata pairing and
successful computation cannot do that.

This gate reads only identifiers and availability metadata. No RNA expression,
DNA gene/segment measurement matrix or ploidy payload is decoded. The existing
DepMap24Q4 Model and default-profile tables establish only default-profile IDs.
They do not prove same-culture/aliquot correspondence or suitable truth.

## Exact exposure exclusions

The captured106-file results catalogue contains samples,lines and eligible-line
metadata. Add the nine-line Andor manifest. Resolve exact ModelIDs and normalized
published names against Model.csv;exclude matched models and every model with
the same nonempty PatientID. Exclude missing patient identity,ambiguous default
profiles,non-cancerous and engineered entries requiring separate review.
No favourable outcome,gene expression or DNA value selects a candidate.

The first pass reveals an important extra exposure source:CCLE preparation
wrote DNA truth identifiers beyond its scored tissue groups. The second pass
therefore includes every Sample identifier in that prepared truth sidecar,
reading only its identifier prefix. It also excludes the entire withdrawn
panel's cell-line metadata. Tahoe remains withdrawn;its expression is unused.
Known CCLE underscore-lineage suffixes receive the same stripped-name matching.
All original exclusions and PatientID closure are reapplied.

| Metadata result | Initial pass | Stricter exposure pass |
|---|---:|---:|
| Model.csv entries |2105|2105|
| Models with paired default-profile metadata |1597|1597|
| Exposure metadata tables |107|109|
| Matched exposed ModelIDs |1142|1690|
| Models sharing exposed PatientIDs |1229|1764|
| Provisional candidate models |402|182|
| Distinct candidate PatientID labels |390|180|

The180 labels are not certification of180 independent biological donors or
single-CNA populations. Candidate models span22 metadata lineages. Ten
unresolved cell-line aliases remain in the stricter report;some other fields
may already resolve those rows,but no unsupported identity correspondence is
assumed. Candidates are explicitly **not cleared for training**.

## Genuine evidence gaps and next admission boundary

Before opening candidate measurement values or training,complete the historical
exposure/alias check,verify measurement presence and assay/units/assembly/gene
identity,establish the available RNA-DNA correspondence and DNA sensitivity,
and freeze development/evaluation roles by biological identity. Within-lineage
DNA variation and response identifiability have not been examined. Do not
consume the entire provisional panel as training or treat it as independent
clinical validation.

The paired default metadata do not supply the59 independent biological
single-CNA controls needed for the stated zero-error5% upper-bound target.
They also do not validate individual-cell genomic membership or spatial spots.
More computation cannot replace those measurements. Former evaluation models
and related donors remain excluded from any subsequent development use.

## Executed code and reproducibility

Both jobs run on nslab account/QOS CPU,1CPU4GB1h. The first uses the existing
isolated analysis Conda environment. The second uses existing read-only
scrna_cnv PyArrow only to export the three withdrawn-panel identifier columns,
then the isolated analysis environment for exclusions. No installation or
shared-environment change occurs. The export records its environment/source
hash and explicitly records that ploidy values are not decoded. Both passes,
all matched source IDs,exclusions and unresolved labels remain preserved.

Jobs23223829/23224054 complete0:0 in3s/14s scheduler elapsed,respectively.
The strict set-analysis process records26236KiB maximum RSS;short-job scheduler
sampling does not provide an informative peak. These are metadata feasibility
checks,not numerical/biological validation or new package-suite executions.
[Complete resource accounting](../dev/dosage_response_feasibility_20261111/accounting_complete.tsv).

[Initial fixed design](../dev/dosage_response_feasibility_20261111/DESIGN.md),
[initial executable metadata gate](../dev/dosage_response_feasibility_20261111/metadata.py),
[initial result23223829](../dev/dosage_response_feasibility_20261111/run_23223829/results/REPORT.json),
[stricter exclusion design](../dev/dosage_response_feasibility_20261111/DESIGN_v2.md),
[identifier-only export](../dev/dosage_response_feasibility_20261111/export_withdrawn.py),
[stricter executable gate](../dev/dosage_response_feasibility_20261111/metadata_v2.py),
[stricter result23224054](../dev/dosage_response_feasibility_20261111/run_23224054/results/REPORT.json),
[provisional candidate identifiers](../dev/dosage_response_feasibility_20261111/run_23224054/results/candidate_metadata.tsv),
[withdrawn/prepared-truth export provenance](../dev/dosage_response_feasibility_20261111/run_23224054/withdrawn_metadata_provenance.json),
[immutable v1](../dev/dosage_response_feasibility_20261111/SNAPSHOT.json),
[immutable v2](../dev/dosage_response_feasibility_20261111/SNAPSHOT_v2.json).
