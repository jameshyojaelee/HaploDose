# OneK1K: possible donor replication, unresolved CNA truth

OneK1K is a potentially useful biological-control resource, not an accepted
single-CNA null panel. The archive describes 982 healthy RNA donors and blood
DNA genotyped with Illumina Global Screening Array-24 v2.0, including 18
replicates. GSE196830 links RNA series GSE196735 and DNA series GSE196829. The
DNA record lists a raw archive and a processed genotype matrix. Neither was
downloaded. Array rows, repeated assays and pooled RNA libraries must not be
counted as independent donors.
[DNA series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196829),
[RNA series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE196735).

An author methods report describes reconstructing DNA mosaic-alteration calls
from deposited IDAT pairs using gtc2vcf and MoChA. Its cell-level Y-loss labels
are based on RNA reads, so they cannot serve as independent DNA membership
truth for an RNA method. The later report describes 981 matched RNA/array donors;
this does not establish the eligible donor count for HaploDose.
[Author methods](https://eprints.soton.ac.uk/507098/1/2025.06.01.25328624v1.full.pdf),
[published report](https://pmc.ncbi.nlm.nih.gov/articles/PMC13135724/).

[Inference] The smallest decisive admission work is a donor/specimen/assay
crosswalk plus an independently frozen DNA measurement protocol. It must resolve
replicates, array manifest and cluster/reference versions, collection timing,
whole-blood versus PBMC correspondence, RNA demultiplexing, and contamination.
Genotype agreement establishes identity, not absence of somatic CNA. Whole-blood
array sensitivity cannot automatically be transferred to a PBMC subpopulation.

After protocol freezing, use DNA-specific intensity and phased-allele measurements
to define callable regions and effect/prevalence limits, with assay-appropriate
QC and sensitivity controls. Do not select donors, thresholds or genomic regions
using HaploDose outputs. A negative array call without a supported detection
limit is insufficient to certify the primary single-CNA event. A regional RNA
programme and a low-prevalence DNA alteration remain competing explanations.
Suitable independent DNA at that resolution, rather than more RNA simulations,
must distinguish them. Same-cell membership remains outside this assay's truth.

Decision: retain the lead for a metadata-first admission protocol; no biological
null specimen is certified, and the 59-independent-specimen requirement remains
unmet. No new cohort is admitted for confirmatory evaluation. Literature findings
were visible during this review; per-donor measurements and HaploDose outcomes
remain unopened. The bounded local search found no OneK1K mention in the inspected
docs, benchmark JSON/Python and package sources; this is not proof that no prior
exposure exists elsewhere.

Compute job22995678 completed in16s. Three self-series SOFT records, query URLs,
SHA256 hashes, isolated wheel identity and environment are preserved under
`dev/onek1k_feasibility_20261108/metadata_v1`; source is `candidate_v1` and script
`benchmarks/hardening/onek1k_metadata.py`. Each request was capped at2MiB and
validated as series metadata without sample data tables. Supplementary links were
recorded but not followed. Process peak RSS was82456KiB; scheduler MaxRSS was
unavailable. No inference, DNA calling, installation or cohort download was run.

## Locked six-sample metadata pilot

Job23001177 retrieves three lexical RNA and three lexical DNA self-sample records
in6s, with a2MiB cap per response and no measurement tables. All six succeed.
RNA samples GSM5899873–GSM5899875 are pools of12,14 and14 people. Their metadata
explicitly link `GenotypeSamples.txt.gz` donor lists and individual barcode maps.
DNA samples GSM5901579–GSM5901581 expose individual IDs1_1,2_2,3_3 and IDAT-pair
URLs. No supplementary links were followed. Pool/sample order is not used as a
donor join. This establishes a concrete source to investigate, not actual matched
donor counts or CNA-negative controls.

Donor-list format job23001231 completes in3s. Its three bounded supplementary
requests yield12,14 and14 literal IDs, with no duplicates within or across pools.
None matches the three pilot array IDs. This does not establish a sample mismatch:
accession order is not donor identity. The next gate records a full exact array
metadata lookup, retaining unmatched/ambiguous IDs without numerical-offset joins.
No counts, intensities or CNA calls have been inspected.

## Array version discrepancy before DNA processing

The preserved GSE196829 series metadata describes GSA-24v2.0, while all three
pilot array sample records describe GSA-24v1.0. Their linked platform is GPL28939.
This discrepancy may be a metadata error or assay heterogeneity; neither is
established. No manifest/cluster resource is selected from these descriptions.

Compute job23001408 retrieves the platform's self-metadata prefix successfully
in1s. Its generic title is a customizable GSA-24 BeadChip and does not resolve the
version conflict. The quick response contains no annotation table; no feature
table or specimen measurements are consumed. Exact response/hash/source records
are under `dev/onek1k_crosswalk_20261108/platform_v1`. Correct manifest and cluster
identity remain necessary before array calling, independently of RNA agreement.

## Bounded pilot IDAT metadata

Job23001585 completes in8s. Four binary parser fixtures pass, including malformed
directories and multi-byte string lengths. Both locked GSM5901579 channels have
IDAT version3, barcode203007520013, positionR01C01, generic chip type
`BeadChip 24x1x4`, and780509 probe entries. This agrees with their deposited
file/sample metadata. Both manifest-name fields are empty. Intensity/probe arrays
are not decoded and no genotypes or CNA calls are computed. Compressed files are
retained with restricted permissions, hashes and exact URLs under
`dev/onek1k_crosswalk_20261108/idat_v1`; no cohort is admitted.

[Inference] The primary format-reference source has a candidate lookup entry
`GSAMD-24v2-0_20024620` with this chip type and probe count. Its lookup also uses a
mid-block count, which this probe did not extract. Matching counts are not unique
manifest identity, and no EGT identity follows. Do not substitute this hint for
independently verified manifest/probe correspondence or original cluster resources.
The v1/v2 discrepancy remains open; no conversion resource is chosen.

The preserved-file extension23001659 completes in3s, retains all four parser
fixtures, and verifies original input hashes. Both channels also have780509
mid-blocks; beadset, description and VLN strings are empty. This completes the
bounded resource-header trace without recovering a manifest or EGT name. The
native lookup's chip/probe/mid-block combination points to the same GSAMD-v2
candidate, but exact resource identity still requires independent manifest/probe
correspondence and a justified cluster resource. More header parsing cannot
supply identifiers that are absent. Original IDATs and first-probe results remain
unchanged; extension outputs are in `dev/onek1k_crosswalk_20261108/idat_v2`.

Selection, responses, execution/source identities and the next bounded donor-list
format gate are under `dev/onek1k_crosswalk_20261108`. Counts, genotype intensities,
barcode maps and CNA calls remain unopened. No cohort is admitted.

## Full administrative crosswalk and assay metadata audit

Job23001291 completes all1104 self-sample metadata retrievals in18m24s without
retrieval failures. Independent audit23001354 passes in6s: all40 literal IDs in
the three locked RNA donor lists have exactly one array accession. The audit
reconstructs the joins from preserved SOFT records and donor lists and verifies
the recorded hashes. These are40 administrative matches, not40 eligible null
specimens, genotype-confirmed identities or independent biological replicates.

Nineteen array records, GSM5902664–GSM5902682, lack the explicit individual-ID
field. They remain unresolved. Inspected endpoint titles/descriptions contain
suffixes such as155_2 and920_2; suffix stripping is not an authorized donor join
and the count must not be substituted for the archive's stated replicate count.
No administrative duplicate appears among the explicit individual-ID fields.

Assay metadata audit23001857 completes in2s. All1104 sample records, including
all40 administrative matches, stateGSA-24v1.0 and platformGPL28939. Every record
has a filename-consistent green/red IDAT pair, with no repeated chip/position
identifier. These metadata checks do not verify biological independence or the
actual array resource. The series-v2/sample-v1 conflict is cohort-wide in the
inspected descriptions; it remains unresolved by the pilot's generic IDAT fields.
Peak process RSS is96000KiB for the join audit and84480KiB for the assay audit.
Outputs, execution manifests and independent sources remain under
`dev/onek1k_crosswalk_20261108/{array_metadata_v1,audit_v1,assay_v1}`.

The [manufacturer's public GSA resource page](https://support.illumina.com/array/array_kits/infinium-global-screening-array/downloads.html)
lists parent-arrayv2 CSV/BPM/EGT resources. Their availability does not establish
identity with the candidate custom GSAMD array. The
[primary MoChA-WDL documentation](https://github.com/freeseek/mochawdl/blob/master/README.md)
describes parent-array manifests as a possible restricted route for custom GSAMD
data, omitting custom probes, while requiring correct manifest/cluster resources.
That route is not adopted: probe/address/allele correspondence, a justified EGT,
callable support and sensitivity must first be established independently of RNA.
No manifest checks are disabled, resources selected or array conversion run.

Decision: the administrative pilot join gate passes; assay identity, specimen
correspondence and DNA sensitivity gates remain open. Preserve biological nulln=0
and withhold the5% biological false-split guarantee. More simulation cannot fill
the missing resource provenance or establish CNA absence. No production method,
frozen result, scientific threshold or regression-only dataset role changes.

### Parent-resource address gate,2026-11-09

The [OneK1K LOY author's conversion script](https://github.com/ad2n15/LOY_Onek1k/blob/main/scripts_for_calling_mCA_from_SNParrays/step_1_idat_to_gtc.sh)
explicitly namesGSA-24v2-0_A1.bpm andGSA-24v2-0_A1_ClusterFile.egt;
the [VCF script](https://github.com/ad2n15/LOY_Onek1k/blob/main/scripts_for_calling_mCA_from_SNParrays/step_2_gtc_to_vcf.sh)
also namesGSA-24v2-0_A1.csv andGRCh37. Preserved source responses identify
this candidate restricted route; filenames do not prove the author's exact
binary versions or compatibility with the deposited custom array.

Frozen`dev/onek1k_crosswalk_20261108/PARENT_ADDRESS_DESIGN.md` permits only
the lockedGSM5901579 IDAT identifier field102 and manufacturer resource metadata.
Compute job23047526 completes in25s; two hand-encoded checks pass. Manufacturer
GRCh37 CSV/BPM/EGT archives and extracted members have recorded hashes. The CSV
contains665608 unique assay names requiring680357 distinct nonzero addresses.
Both pilot channels contain the same780509 identifiers. Of the parent-required
addresses,680191 are present and166 absent;665442 probes have complete address
support and166 lack an address. Another100318 pilot addresses lie outside the
parent assay address set. They are not automatically all custom biological probes.

Independent job23047779 completes in10s. A separate integer/NumPy IDAT parser and
streamed CSV comparison reproduce every665608 exported row, all counts and all
166 missing-probe identities. Full parent-address compatibility is therefore
falsified for this pilot. No manifest-check bypass, missing-address imputation,
intensity decoding, genotype conversion or cohort admission follows. Complete
address overlap itself would not establish allele/probe sequence identity,
cluster validity, normalization or DNA sensitivity. A restricted manifest or
the exact custom resources still need justified binary/probe correspondence
before conversion. Independently certified biological nulln remains0.

Outputs:`parent_addresses_v1` and`parent_address_audit_v1`; immutable sources
`parent_address_source_v1` and`parent_address_audit_source_v1` in the same task.
Process peak RSS3510408/407292KiB respectively; both exits0, CPU backend recorded
in the isolated validated environment. No package source, existing tests,
scientific criteria or frozen historical outputs changed.

One prespecified public custom-CSV retrieval23048243 fails with HTTP404 before
archive decoding. The [primary researcher resource page](https://jennysjaarda.github.io/PSYMETAB/data_sources.html)
links`https://glimdna.org/assets/gsamd-24v2-0_20024620_a1.zip`; its public
redirected resource is unavailable in the executed request. Frozen protocol
`CUSTOM_RESOURCE_DESIGN.md`, source`custom_resource_source_v1`, failed output
`custom_resource_v1/run_manifest.json` and traceback are retained. No guessed
mirror, other cohort's cluster resource or intensity conversion is substituted.
The exact custom-resource gap remains open.
