# Benchmark datasets with DNA truth

Compiled 2026-11-01. In-house sets come from an inventory of `~/CRC`; public sets come from a literature
search, with accessions checked against the papers' data-availability statements and the GEO/ENA records.
UNVERIFIED marks anything that was not confirmed. Nothing listed here has been downloaded yet unless it says so.

## In-house (CRC), in use

| Set | RNA | DNA truth | n | Status |
|---|---|---|---|---|
| U-CAN Visium HD | `share/Spaceranger/novaseq/U*` (probe-based, so no SNP alleles) | WGS FACETS: TCN, LCN, CF, purity, ploidy (`~/CRC/data/external/eva_cnv/vcf`) | 20 patients (3 IGF2 TD; 8 with ploidy >= 2.7) | `benchmarks/ucan_spatial/` |
| Moorman (HTAN HTA8) | 10x 3' scRNA | MSK-IMPACT log2 segments, hg19 | 9 specimens / 8 patients | `benchmarks/moorman/` |
| GSE149224 cell lines | Drop-seq-like, HCT116/SW480/RKO, 5-FU doses | DepMap 24Q4 gene CN | 3 lines | `benchmarks/cell_lines/`. The GEO SW480 labels do not match DepMap SW480, so this set only tests batch and dose robustness |
| YCC (GSE273567) | 10x, reads available | none (native Numbat only) | 5 | comparator only |
| Bian scTrio-seq2 | same-cell full-length RNA + methylome | CNV from bisulfite reads | 10 patients | not started |
| TCGA, U-CAN bulk, CPTAC, PCAWG | bulk RNA | ASCAT / FACETS / GISTIC / consensus CN | ~620 / 1,063 / 105 / 2,778 | planned: gene dosage response and event library |

## Public, ranked by how decisive a test they give

| # | Dataset | Accession (access) | Cancer | RNA | DNA truth | Why use it |
|---|---|---|---|---|---|---|
| 1 | Andor 2020 NAR Genom Bioinform | scDNA PRJNA498809; scRNA PRJNA598203 / GSE142750 (open) | 9 gastric cell lines | 10x, FASTQ | 10x scCNV, 2-12 clones per line | Only set with open reads for both scDNA and scRNA. Used by Schmid, Chang, Numbat and Clonalscope, so results are directly comparable |
| 2 | McPherson 2025 Nature (MSK SPECTRUM) | RNA GSE180661 (open); DNA dbGaP phs002857 (controlled); processed Synapse syn66366960 (openness UNVERIFIED) | HGSOC | 10x 3', FASTQ | DLP+ scWGS, per-cell WGD | Strongest test of WGD (27/41 patients WGD-high) |
| 3 | Kinker 2020 Nat Genet | GSE157220 (open) | 198 pooled lines | 10x | CCLE/DepMap CN (on disk: DepMap 24Q4) | Line identity gives clone truth at scale |
| 4 | DNTR-seq, Zachariadis 2020 Mol Cell | GSE144296 (open; ALL patients restricted) | HCT116 (CRC line), A375 | plate-based, same cell | scWGS, same cell | Per-cell assignment accuracy |
| 5 | Funnell 2022 Nature (SIGNALS) | processed Zenodo 6998936 (open); raw EGAS00001006343 (controlled) | TNBC, HGSC, 184-hTERT | 10x | DLP+ haplotype-specific CN | LOH and allele-specific truth |
| 6 | XClone GX109-T1c | RNA GSE232733 (open); WES/scDNA EGAS00001007854 (controlled) | gastric | 10x 5' | WES + scDNA | LOH |
| 7 | COLO320-HSR | WGS PRJNA506071; multiome PRJNA672109 (open) | CRC line | multiome | WGS | CRC line with ecDNA; used by Schmid |
| 8 | Slide-DNA-seq, Zhao 2022 Nature | PRJNA768453 (open) | CRC | Slide-seq v2 | Slide-DNA-seq on an adjacent section | Spatial CRC with spatial DNA |
| 9 | scONE-seq, Yu 2023 Sci Adv | GSE185269 (open) | glioma, NPC43 | same cell | scWGS | Per-cell |
| 10 | Ji lab GI (Clonalscope; Chang 2026) | dbGaP phs001818, phs001711 (controlled) | gastric; CRC liver mets | 10x | scDNA / WGS | CRC with WGS; needs a dbGaP application |
| 11 | HTAN WUSTL (CalicoST) | dbGaP phs002371 (controlled) | CRC, breast | Visium | WES | Spatial CRC with WES; needs a dbGaP application |

Gap: the search found no open-access CRC 10x scRNA set with matched bulk WGS/WES. The U-CAN Visium HD
sections (in-house) are therefore the main CRC test with DNA truth.
