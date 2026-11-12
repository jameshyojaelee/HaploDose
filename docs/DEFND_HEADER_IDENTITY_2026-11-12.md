# DEFND identity preparation: exclude labels and translate bead barcodes explicitly

The goal remains accurate genomic RNA profiles and supported CNA memberships
with independent performance evidence. New private identity helpers now pass
compute-node qualification: author-label changes leave identities unchanged,
and explicit ARC bead translation agrees exactly with the vendor reference.
Actual DEFND headers,correction maps and genomic truth remain unverified.
No patient data,method scores or biological validation units are added.

## Why literal matching would be unsafe

The previous metadata task establishes that deposited RNA headers contain author
cluster numbers. The new helper requires a declared molecular prefix grammar,
delimiter,orientation and whitelist. It exports only library-qualified molecular
identity. It does not search suffix labels for a plausible barcode,guess an
orientation,drop GEM groups,aggregate duplicates or choose a best matching rule.
An unknown or ambiguous identity stops before indexing. Synthetic prefixes and
delimiters are fixtures,not assertions about the real deposit format.

There is a separate chemistry requirement. [Official10x documentation](https://www.10xgenomics.com/support/jp/software/cell-ranger-arc/latest/analysis/atac-barcoded-bam)
describes distinct paired ATAC and GEX barcode sequences on Multiome beads.
The relevant rule is bead identity for the declared kit. No ATAC count/noise
model is transferred to DEFND genomic DNA. The installed ARC2.0.2 source loads
two ordered vendor lists to build the paired translation. Merely taking a
reverse complement is not a substitute for an independently specified pairing.

Full reference qualification finds **736320 ordered pairs**,zero sequence
intersection between the ATAC/GEX lists and zero identity pairs. Thus literal
sequence matching cannot join these two raw reference namespaces. This result
does not establish whether the deposited DEFND DNA IDs are raw ATAC,already
translated GEX,or separately corrected. Applying the map to the wrong namespace
would create another error; the helper explicitly rejects repeated translation.

## Scientific contract and independent tests

The estimands are molecular header identity and vendor-declared bead pairing.
The observations are formatted strings and ordered reference lists. Required
identifying information is the exact library,grammar,orientation,chemistry,
whitelist versions,pairing direction and correction history. A matching-looking
name could instead reflect an author label,different molecular namespace or
upstream merged barcode. Artifact-specific provenance and a validated crosswalk
would distinguish these alternatives. Neither estimand identifies a DNA CNA.

The prefix helper passes1006 fixed/seeded label substitutions for two independently
specified IDs,including suffixes containing another barcode. Library namespaces,
GEM1/2,an independent reverse-complement example and ID/value permutations pass.
All17 invalid cases reject,including duplicate/unknown IDs,short/long sequences,
implicit GEM removal,automatic orientation,empty columns and decorated whitelists.
Error messages do not expose the raw author labels. The seed is202611120031.

The separate paired helper requires a one-to-one ordered map,distinct source/
target namespaces and a provenance identity. It preserves GEM suffixes and
rejects unknown IDs,wrong namespaces,repeated translation,unequal/empty lists,
duplicate source/target IDs,shortened/decorated sequences and noncanonical groups.
All11 counterexamples reject. A hand-declared mapping differs from identity and
reverse complement. Applying the same row permutation to both lists preserves it.

On all vendor reference entries,**4417920 native comparisons and4417920 inverse
roundtrips are exact**,covering both directions and bare/GEM1/GEM2 IDs. Candidate
file reading is independent of the vendor loaders. The reference extracts only
the installed loader functions and TranslationTable AST; it does not import or
run a Cell Ranger pipeline. All six actual vendor file hashes match before/after.
No numerical tolerance is needed for string identities. Reference barcodes are
not independent cells or biological specimens.

| Vendor reference | Compressed bytes | SHA256 |
|---|---:|---|
| ATAC737K-arc-v1 | 2465078 | 9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26 |
| GEX737K-arc-v1 | 2142553 | 1279bb7ba28ddcd805e2f16ebd5d6933b7b1e7f128e18aabfe6ae8cd7b7bf710 |

## RNA barcode correction is another unresolved upstream measurement

The GEO-declared RNA repository is inspected through the connected GitHub tool.
Four exact public blobs are recorded; their Git blob hashes and byte sizes
rederive independently. No external RNA source is executed or installed.

| File | Bytes | Git blob SHA1 |
|---|---:|---|
| README.md | 1033 | 23bfea57cc3e6ca2b4cf99547f88f3d914cf6608 |
| DropSeqPipeline8.py | 13224 | 46b3a04e9e462ea994d998430dff02670b7ffe4d |
| collapse.py | 5964 | 6e9d479128289a3cba3e57f7fa5b3295e8aa43ca |
| collapse2.py | 5700 | c04c4a032ca8ae5c964b34b362d76916d66bf7de |

The inspected [collapse source](https://github.com/simslab/DropSeqPipeline8/blob/main/collapse.py)
can shorten16-base technology barcodes during truncation repair. The inspected
[second correction](https://github.com/simslab/DropSeqPipeline8/blob/main/collapse2.py)
can merge barcodes using RNA UMI/gene overlap and writes old/new tracking maps.
**[Unverified] Whether these exact versions and settings produced the deposits.**
If such corrections occurred, padding/truncation or an apparent exact match
cannot recover a justified physical-cell mapping. The archived correction maps
and run settings,or decisive bounded upstream evidence,remain required where
the identities depend on those operations. These inspected files also do not
establish the released matrix's author-label delimiter. Current source capability
is not evidence that a particular specimen contains shortened or merged IDs.

## Every execution and preservation record

| Job | Progress | State / exit | Allocation and node | Runtime / process peak |
|---|---|---|---|---|
| 23238413 header contract | 1006/1006 ██████████ | Complete0:0 | CPU1/1GiB/10min,ne1dc4-068 | Slurm1s,timed1.12s / 11584KiB |
| 23238599 paired reference | 4417920/4417920 ██████████ | Complete0:0 | CPU1/4GiB/10min,ne1dc3-010 | Slurm1m18s,timed1m13.53s / 473916KiB |

Both use nslab account/QOS,the existing isolated read-only Python3.11.18
environment and fresh neutral outputs. The first32 and second46 source/cached
identities match before/after. The paired job loads the existing system
cellranger-arc/2.0.2 module; source and both barcode assets are locked on compute.
Slurm sampled peak for pairing is372452KiB; the header job's0 is an unavailable
short-job sample. CPU standard-library execution,no inference kernel or GPU.
Two complete,zero running/pending,zero failed. All headers were announced.
No arrays,installs,production dependencies,shared/running-job changes or public
Python/Rust rerun. Original header code and completed source inventory remain.

## Admission and next decisive observation

This is a qualified private identity layer,not an integrated matrix loader or
same-cell genomic validation. Before patient measurement access,the actual
grammar and correction/translation status plus complete evaluation design must
be frozen. A bounded header/fragment inspection must distinguish raw versus
translated namespaces and check unique library-qualified pairs without exposing
author labels or using CNA agreement to choose a mapping. A partial coordinate-
sorted prefix cannot establish whole-file identity or population matching rates.
An unknown format,collision,unresolved correction or assembly contradiction stops
the path. Broader staging and scoring do not follow automatically.

DNA reference/provenance,multiplicity convention,assay accessibility controls,
resolution and biological replication remain separate gates. Cohort admission
remains false. Original recovery2/5,state errors14/200,power195/200,coverage54/1000,
fixed competitors and all frozen results remain unchanged. No accuracy gain or
all-competitor superiority is established by barcode-reference agreement.

[Header contract](../dev/defnd_header_contract_20261112/identity.py),
[qualification design](../dev/defnd_header_contract_20261112/DESIGN.md),
[header result](../dev/defnd_header_contract_20261112/run_23238413/results/QUALIFICATION.json),
[explicit paired implementation](../dev/defnd_header_contract_20261112/paired.py),
[pairing design](../dev/defnd_header_contract_20261112/PAIRING_DESIGN.md),
[complete paired result](../dev/defnd_header_contract_20261112/pairing_23238599/results/QUALIFICATION.json),
[actual locked vendor identities](../dev/defnd_header_contract_20261112/pairing_23238599/results/ASSETS_LOCKED.json),
[exact public RNA source snapshots](../dev/defnd_header_contract_20261112/PUBLIC_RNA_SOURCE.json),
[first submission](../dev/defnd_header_contract_20261112/SUBMISSION.json),
[paired submission](../dev/defnd_header_contract_20261112/PAIRING_SUBMISSION.json).
