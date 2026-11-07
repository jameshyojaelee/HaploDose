# Claim 5 run plan: HaploDose in the Schmid et al. 2025 benchmark pipeline

Status 2026-11-02: wrapper written and tested end to end on Kinker development data. No Schmid raw or truth
file has been opened (listings and the truth manifest only). Nothing below is installed. Claim definition:
`docs/STATUS.md`, item 5.

## Pinned inputs

- Pipeline: `external/benchmark_scrnaseq_cnv_callers`, commit `319ea34ddfc46bdd4c0cb1ed534f4ece0f5c44ca`
  (see `external/README.md`).
- Gene positions for inferCNV and HaploDose: `external/annotations/hg38_gencode_v27.txt`, md5
  `0c08b487706b2b8362470570d5b9eeb8`; copy to their `data/annotations/`.
- Truth: `data/external/schmid2025_truth/` (Zenodo 20260649, md5 checked against the API at download).
- HaploDose: frozen copy given as `--code` (PYTHONPATH), today
  `results/scaling_frozen1424_20261102_1541/code` (a symlink to `results/heldout_20261102_1424/code`). Replace
  with the v-next frozen path once v-next is frozen. The wrapper refuses to run if `haplodose` imports from
  anywhere else. Check the frozen copy with `sha256sum -c code_sha256.txt` in the job, as `test/run_test.sbatch`
  does.

## Files in benchmarks/schmid

| File | Role |
|---|---|
| `run_haplodose.py` | Wrapper. Reads their three inputs plus the gene order file, writes CSR H5, cells and positions tables, runs `python -s -m haplodose.cli run` on the frozen code (CPU by default, `JAX_PLATFORMS=cpu`), and writes `<dataset>_haplodose_gene_by_clone.tsv`, `<dataset>_haplodose_cell_clone.tsv` and `run_haplodose.json` (args, code SHA-256, timings). |
| `lib_haplodose.R` | `read_haplodose()` and `read_haplodose_individual_cells()`, written to mirror `read_copykat()`. They return a per-gene GRanges where the score = cell-weighted mean clone log ratio over non-reference cells + 2, plus the cell count. |
| `haplodose.smk` | Two rules in their style. `run_haplodose` takes the same inputs as `run_copykat`, with the code path from `--config haplodose_code=...`. `evaluate_cnv_prediction_wgs_haplodose` is their WGS evaluation rule plus the two HaploDose inputs, writing to `evaluation_haplodose/`. |
| `make_patched_scorers.py` | Writes `evaluate_results_haplodose.R` and `evaluate_results_permethod_haplodose.R` from the pinned originals with literal, count-checked edits, plus `evaluate_results_haplodose.patch`. Three edits: source `lib_haplodose.R`, add an optional HaploDose block, add `haplodose = "HaploDose"` to `method_names`. |
| `snakemake_shim.R` | A stand-in `snakemake` S4 object, so their scripts run with Rscript where snakemake is not installed. |
| `test/` | Kinker reshaping (`build_kinker_input.py`), `write_count_matrix.R` (their `write.table` call), `score_test.R`, and `run_test.sbatch`. |

Setup in their tree: copy `run_haplodose.py`, `lib_haplodose.R` and the two `*_haplodose.R` scorers into
`snakemake_pipeline/scripts/`, then add `include: "haplodose.smk"` to `workflow.sm`. Do this in a working copy
of the clone, not the pinned clone.

## Test result (Kinker development data, not Schmid)

Job 22619319 ran on ne1dc6-008 (AMD EPYC 9254, a v6 node) in 52 s, state COMPLETED. Output is in
`results/schmid_wrapper_test_20261102_2058/`.

- Input: Kinker pool18, rotation r0. 554 C32_SKIN query cells, plus 1,653 cells from the 12 reference lines,
  each line its own reference group. 22,722 genes. `count_matrix.txt` was written with R `write.table`, as their
  preprocessing does.
- Wrapper: the count parse round trip is identical to the source counts (re-derived from `counts.h5`, sum
  38,102,774 UMI). HaploDose imported from the frozen path. Reference split was 1,158 fit / 495 holdout. All 554
  query cells were assigned (clone c1 511, clone c2 43). 11,605 genes went out. HaploDose took 19.5 s and
  1.5 GB peak RSS.
- Reader checks pass. The pseudobulk equals the mean of the per-cell matrix (to 1e-9), and no reference cell
  enters the pseudobulk.
- Their functions, per-method path, 1,657 of 1,751 truth tiles: Pearson 0.514, AUC gains 0.844, AUC losses
  0.783, truncated AUC gains 0.815 and losses 0.695. There were 187 gain tiles and 369 loss tiles.
- These numbers are a plumbing check only. The truth is DepMap gene CN relative to the reference mix, put on
  1 Mb tiles with thresholds of ±0.15 that I set for the test. It is not calibrated and not comparable to any
  Schmid number.
- Not run:
  - F1 (`crfsuite` is missing).
  - The two full patched scorers. They need all comparator outputs and `library(numbat)` / `library(crfsuite)`.
  - The Snakemake rules (snakemake is not installed).

## Decisions to fix before the claim 5 run (each changes numbers)

1. **Bin set.** Their combined scorer applies `combine_range_objects` in sequence, so each method's bins are
   the intersection across all methods. Frozen HaploDose returns autosomes only (the test output has no chrX
   gene), so adding it removes chrX tiles for every tool. It also changes the comparators' bins relative to
   the paper.
   - Recommendation: the primary analysis is their combined scorer with HaploDose included (one common bin
     set, which is the paper's headline approach).
   - Report two sensitivity analyses: (a) the comparators on the unpatched scorer (paper bins), and
     (b) `evaluate_results_permethod_haplodose.R` (each method on its own bins).
2. **Which outputs count as "the six tools".** The paper's summary ranks 9 outputs: inferCNV, SCEVAN and Numbat
   each in CNV and expression form.
   - Recommendation: rank HaploDose against all 9 outputs as the paper lists them. Note that the paper's SCEVAN
     (CNV) key is `scevan_vega` in `visualize_overall_method_performance.R` but `scevan_cnv` in
     `evaluate_results.R`; resolve this when wiring the WES and CNVkit rules.
3. **Rank rule.** The paper's overall order is the mean of normalised Pearson (negative set to 0), AUC gains and
   losses (rescaled as (AUC - 0.5) x 2, floored at 0) and max F1. STATUS item 5 says "Pearson and gain/loss
   AUC" and "best mean rank on both metrics".
   - Recommendation: three per-dataset ranks (Pearson, AUC gains, AUC losses), averaged over the 4 datasets
     and the 3 metrics. HaploDose is SOTA only if it is first on that mean.
   - Also report the paper's normalised overall score with F1, descriptively.
   - Write this into STATUS before any Schmid run.
4. **HaploDose settings** (wrapper defaults, chosen before any Schmid data):
   - Value column `log_ratio_gene`, which matches claim 1.
   - Reference holdout 0.3, as in the Kinker and MIX-seq runs.
   - Seed 20261103. CLI defaults otherwise. CPU.
   - Each Schmid reference group (e.g. Endothelial, Melanocytes, Myofibroblasts in BCC) becomes one HaploDose
     reference type.
5. **Comparator versions.**
   - Paper versions: inferCNV 1.10.0, copyKat 1.1.0, SCEVAN 1.0.1, Numbat 1.4.0, CaSpER 0.2.0, CONICSmat
     0.0.0.1.
   - Installed versions here: copyKat 1.2.5 and SCEVAN 1.0.3 (cnv_field), Numbat 1.5.2 (numbat_env).
   - Recommendation: use the paper's versions. As a gate, reproduce their published comparator values on one
     primary dataset before adding HaploDose to the scoring. A large gap means our rerun is not their benchmark.

## Environments per tool

| Tool | Needs | Covered by an existing env? |
|---|---|---|
| inferCNV | R, infercnv, gene order file | Partly. infercnv is in the R/4.4.1 module user library (version not recorded; paper used 1.10.0). |
| copyKat | R, copykat | Partly. cnv_field (R 4.4.3) has copykat 1.2.5; paper used 1.1.0. |
| SCEVAN | R, SCEVAN, yaGST | Partly. cnv_field has SCEVAN 1.0.3 and yaGST; paper used 1.0.1. |
| Numbat | numbat; `pileup_and_phase.R` with cellsnp-lite, Eagle 2.4.1, 1000G hg38 panel, genetic_map_hg38_withX; a **BAM** per dataset | Partly. numbat_env has numbat 1.5.2 (paper 1.4.0), cellsnp-lite, eagle and samtools. The panel and map are in `CRC/data/external/gse273567/`. BAMs do not exist yet. |
| CaSpER | CaSpER, BAFExtract, samtools, biomaRt annotation (**internet** in `casper_gene_annotations.R`), a **BAM** | No. |
| CONICSmat | CONICSmat, `getGenePositions` from ensembl.org (**internet**), `chromosome_arm_positions_grch38.txt` (in their repo) | No. |
| HaploDose | haplodose env (Python 3.11) plus the frozen code | Yes. |
| Scorer | data.table, dplyr, GenomicRanges, reshape2, ggplot2, ggpubr, ROCR, pROC, viridis, numbat, crfsuite | Partly. rnaseq env (R 4.4.3) has all but numbat and crfsuite. |
| Workflow | snakemake 7.32.4 | No. |
| scrna_cnv (micromamba) | infercnvpy only (Python) | Not used by their pipeline. |

## New dependencies (listed, not installed)

Recommended route: one environment built from their `conda_env_R.yml`, with R 4.1.0, bioconductor-infercnv
1.10.0, snakemake 7.32.4 and cellsnp-lite 1.2.3. This about 500-package env is the new production dependency.
Then:

1. GitHub R packages, pinned to the paper's versions by tag or commit. Their README installs unpinned HEAD.
   - `akdess/CaSpER` (0.2.0)
   - `navinlabcode/copykat` (1.1.0)
   - `miccec/yaGST`
   - `AntonioDeFalco/SCEVAN` (1.0.1)
   - `diazlab/CONICS/CONICSmat`
2. CRAN packages: `numbat` 1.4.0 (archive), `Seurat` (only for preprocessing), `clevr`, `crfsuite`.
3. `ROCR`. `lib.R` loads it, but I found no `r-rocr` line in `conda_env_R.yml`, so it was probably installed by
   hand [Inference].
4. BAFExtract (C++ build from `akdess/BAFExtract`), plus its hg38 genome list and pileup inputs, for CaSpER.
5. Eagle 2.4.1 for Numbat phasing. The tarball is at `CRC/data/external/gse273567/Eagle_v2.4.1.tar.gz`
   (md5 7e517c616ae7dd281b878de6e5f7fa26). The env build copies it to `external/schmid_env/tools/` and
   unpacks it there.
6. CellRanger 7.0.0 with `refdata-gex-GRCh38-2020-A`, as in their `preprocessing_datasets/cellranger_*.sh`.
   COLO320 also uses `--chemistry=ARC-v1`. Only CellRanger 8.0.1 is available as a module, and the 2020-A
   reference has not been located. Using 8.0.1 is a version deviation to record. [Unverified] whether it
   changes counts for these runs.
6b. `zellkonverter` and `Seurat`. Their reference preparation calls `readH5AD`, and neither package is in
   `conda_env_R.yml`. Not needed if references are subset as in "Open design point" below.
7. Network access from the job for CaSpER (biomaRt) and CONICSmat (ensembl.org). Otherwise their scripts need a
   local-annotation patch, which then becomes a deviation from the paper.

Minimal fallback, if the paper env is rejected: add `crfsuite` and `numbat` to the rnaseq env for the scorer.
Run each tool in the closest existing env and record the version differences above as a limitation.

## Data gaps (from listings only)

- **MM199** (primary set): the SRA run and the GSE141299 tar are being fetched after the ENA step of data-fetch
  job 22612843.
  - The annotation is rebuilt inside `prepare_MM_sciCNV_scran.Rmd`, lines 250-302 (clusters and markers, then
    `input_MM/sample_annotation.txt`). This was reported by data-fetch; I have not re-read those lines.
  - `input_MM_old` is used only for plots and cannot be obtained.
  - Plan: run their Rmd as-is minus the old-annotation plotting lines, and record this as a deviation.
- **COLO320-HSR**: `SRR12900668_2.fastq.gz` was still downloading when listed, not failed (data-fetch). Wait
  for its `.ok` marker before alignment. The reference file `Adult_healthy_large_intestine.h5ad` has no
  stated source; see "Open design point".
- **BCC06, BCC06post**: their count matrices come from GEO GSE123813 processed files. data-fetch is downloading
  them in job 22619347 to `data/external/yost2019_gse123813/suppl/`. The FASTQs (SRR8315756, SRR8315759) are needed only for BAM-based tools (Numbat, CaSpER).
- **A375, HCT116** (line-overlap set): DNTR-seq is plate-based full-length RNA. HaploDose and its count model
  were calibrated on droplet UMI counts, and `prepare_A375_DNTRseq.R` also builds a TPM matrix. The wrapper
  rejects non-integer input. Running HaploDose here is a modality port and must be labelled as such, or these
  sets dropped. `DNTRseq_A375_rna` holds a single SRR of 42 MB, which looks incomplete [Inference].
- **MCF7**: FASTQ present (SRR10018149).

## Open design point: external normal references (fix in writing before any Schmid run)

The authors' scripts give only local reference filenames. data-fetch (job 22619248) fetched the best-match
public sources. I opened only their annotations (obs) and checked that counts are integers. This is
normal-cell reference data, not benchmark data. The reference choice changes every tool's result, so record
the rule below, or a replacement, in STATUS before any Schmid file is opened.

**COLO320.** `prepare_COLO320_rep1_withRef.R` does the following:
- Reads `Adult_healthy_large_intestine.h5ad` and uses `X` as raw counts.
- Intersects genes with COLO320 by symbol.
- Keeps `category == "Epithelial"`.
- Runs `set.seed(1)` then `sample(1:ncol, 5000)`, with group label `Epithelial`.

There is no age, health or region filter in the script, so the file name implies the file was already subset.

Candidates, both in `data/external/gutcellatlas_elmentaite2021/`. Both have integer `X` and gene symbols.
- `Full_obj_raw_counts_nosoupx_v2.h5ad`: 428,469 cells. It has a `category` column with the same labels the
  script filters on.
- `epi_raw_counts02_v2.h5ad`: 142,113 cells. All of them are also in the full object.

Healthy adult epithelial cells (`Diagnosis == "Healthy adult"`) by `Region`:
- LargeInt 38,612 (19,104 3' and 19,508 5' chemistry)
- REC 17,348
- SmallInt 16,692
- APD 4,689

Proposed rule:
- Use `Full_obj_raw_counts_nosoupx_v2.h5ad`. The script filters on `category`, which implies a multi-lineage
  object [Inference].
- Keep `Diagnosis == "Healthy adult"` and `Region == "LargeInt"`; rectum (REC) and appendix (APD) are coded
  separately in the atlas.
- Keep `category == "Epithelial"`, in atlas row order. That gives n = 38,612.
- Intersect genes by symbol.
- In R (4.x, default `sample.kind = "Rejection"`), run `set.seed(1); idx <- sample(1:n, 5000)`. The indices
  depend only on n, so they can be computed in R and applied in Python without zellkonverter.

This reproduces the authors' 5,000 cells only if their file held the same cells in the same order. That cannot
be verified, so the rule is a documented best match, not a replication [Inference].

Fallback, fixed now: if the files differ in a way that matters, keep the same rule (any n, same seed and draw).
Report the reference as a deviation from the paper. Write the 5,000 cell IDs to the run directory so every
tool gets the identical reference.

Do not pick between LargeInt and LargeInt + REC after seeing results. If a sensitivity analysis is wanted,
name it now.

**MCF7** (line-overlap set only). `prepare_MCF7_dataset_withRef.R` does the following:
- Reads `TS_Mammary.h5ad` and uses `layers/raw_counts` (integer in our copy).
- Keeps `method == "10X"` and `compartment == "epithelial"`. There is no subsampling.
- Uses `cell_ontology_class` with spaces replaced by `_` as reference groups.

In `data/external/tabula_sapiens_v1_figshare/TS_Mammary.h5ad` (figshare file 34701997, donor TSP4) this gives
5,630 cells:
- `luminal_epithelial_cell_of_mammary_gland`: 5,160
- `basal_cell`: 470

The rule is deterministic, so it matches exactly if the figshare v1 file is the one the authors used.
[Unverified]: their path `data/Tabular_sapiens/` does not name a version. Fallback: use the same rule on this
file and record the version.

**Selected cell IDs (saved 2026-11-02, job 22619450).** They are in
`results/schmid_reference_cells_20261102_2119/`: `COLO320_reference_cells.tsv` (5,000 cells drawn from 38,612),
`MCF7_reference_cells.tsv` (5,630 cells) and `summary.json`.
- Scripts are in `benchmarks/schmid/reference/`.
- The COLO320 draw ran in schmid_r41 with R 4.1.0 (Mersenne-Twister, Rejection) and is saved as
  `colo320_r_indices.txt`. It was applied to the atlas obs in row order. Only annotations were read, no counts.
- Every tool's COLO320 and MCF7 input must use exactly these IDs.

**For all tools.** The external references come from other donors and other 10x chemistries (the gut atlas
mixes 3' and 5'), while COLO320 is ARC-v1 multiome GEX. This reference mismatch is the same for all tools.
For HaploDose, each reference group is one reference type, and the 30% holdout check reports how well the
reference fits.

## Resources per dataset

| Step | Ask per dataset | Basis |
|---|---|---|
| CellRanger count (COLO320, BCC06, BCC06post, MCF7, MM) | 16 CPUs, 128 GB, 90 h limit, `cpu`, `--prefer=v6` | [Estimate] 3-8 h for BCC/MCF7, 15-30 h for COLO320 (43 GB FASTQ); see "CellRanger plan" |
| HaploDose wrapper | 8 CPUs, 32 GB, 90 h limit | Test: 2,207 cells in 20 s at 1.5 GB. Schmid sets are larger; [Estimate] under 1 h and under 16 GB |
| inferCNV | 1 CPU, 32 GB, 48 h | their rule |
| copyKat, SCEVAN, CONICSmat | 1 CPU, 32 GB, 24 h each | their rules |
| Numbat preprocessing + run | 1 CPU, 32 GB, 24 h + 48 h | their rules (`--ncores {threads}`, so more CPUs would help) |
| CaSpER (BAFExtract + run) | 24 GB 24 h, then 64 GB 48 h | their rules |
| Scorers (combined + per-method) | 1 CPU, 24 GB, 24 h | their rules |

Per primary dataset, the comparator total is about 300 GB-hours of memory reservation [Estimate]. Run the
datasets sequentially or as separate jobs (no arrays). Wall time limits default to 90 h, per CLAUDE.md.

## Comparator env and gate (2026-11-02, user-approved)

**Env `schmid_r41`.** It is built by `external/schmid_env/build.sbatch`, with outputs in
`external/schmid_env/build_<ts>/`: env export, `versions.tsv`, `sessionInfo.txt`, tool versions.
- The conda part is their `conda_env_R.yml`, exact builds first. The fallback drops build strings, and the job
  records which one was used.
- GitHub packages are installed from local clones at the commits in `external/schmid_env/refs.tsv`.
- CRAN packages come from the Posit snapshot of 2024-02-29, the date numbat 1.4.0 was the current release.
  Bioconductor is 3.14.
- BAFExtract is built with make in `external/schmid_env/tools/`, and Eagle 2.4.1 is unpacked there.

Pinning caveats:
- SCEVAN 1.0.1 spans several commits from 2021-12 to 2025-02, and some of them change the code. I chose the
  last 1.0.1 commit (92761f6, 2024-06-11) because the paper's env export postdates 2024-12-31 [Inference].
  The v1.0.1 tag (90b45b9, 2024-01-16) is the alternative. This is not an exact pin.
- copyKat 1.1.0 is pinned at f34eff7 (2022-09-23). The 2026 tag commit differs only in its description text.
- zellkonverter is not installed. The reference subsets are made without it, per amendment (f).
- Build result (jobs 22619382, 22619434, 22619440):
  - All seven pinned versions load in the env: inferCNV 1.10.0, copyKat 1.1.0, SCEVAN 1.0.1, CaSpER 0.2.0,
    CONICSmat 0.0.0.1, Numbat 1.4.0, yaGST 2017.8.25.
  - The conda solve used the exact build strings.
  - Five conda R packages were replaced by newer CRAN versions: ggplot2, igraph (1.3.0 -> 2.0.2), Matrix,
    scales, vctrs.
  - tidyverse is not installed. Details are in `external/schmid_env/README.md`.

**Gate** (`benchmarks/schmid/gate/`, prepared, not submitted). It waits for team-lead's go on data completeness
and the v-next freeze.
- `run_gate.sbatch OUT DATASET INPUT_DIR TRUTH BAM_DIR` runs the six comparators and their unpatched scorer for
  one primary dataset, in a working copy of the pinned pipeline. It uses the schmid_r41 env and requests 16 CPUs,
  160 GB and 90 h; snakemake gets `--cores 8` and a 150 GB memory budget. It fails early if an input or
  annotation is missing.
- `compare_published.py` writes ours vs published per output and metric. The published values must be
  transcribed from the paper's Source Data into a TSV; none are typed into code.
- Pass rule: STATUS claim 5 amendment (g).
- `published_values.tsv` (183 rows) is made by `transcribe_published.py` from the paper's Source Data,
  `41467_2025_62359_MOESM11_ESM.zip:figure_source_data/suppfigure_4_metrics_extended.tsv` (md5 of the zip in
  `external/schmid_env/ref/MANIFEST.tsv`). The zip holds TSVs, not sheets, so `sheet` is `NA (tsv)`.
  - Datasets: COLO320, BCC06, BCC06post, MM (primary) and A375, MCF7, HCT116 (line overlap). Metrics: Pearson,
    AUC gain and AUC loss, as raw `value` (not `norm_value`).
  - The AUC is the ROCR full AUC our scorer writes as `auc_gains`/`auc_losses` (lib.R:1309, 1325). The source's
    `norm_value` equals max(0, (AUC - 0.5) * 2) on all 122 AUC rows, which matches the paper's summary script.
    Main figure 2A shows the pROC partial AUC instead, a different quantity, so it is not used.
  - Pearson is rounded to 3 decimals in the source. Two Pearson values are negative (BCC06 scevan_cnv -0.021,
    scevan -0.030); their `norm_value` is floored at 0, and the gate uses the raw value.
  - Absent from the source: A375 casper and HCT116 numbat_cnv. A375 and HCT116 are plate-based.

**Gate reference inputs** are in `external/schmid_env/ref_20261102_2215/`, linked as `external/schmid_env/ref`
(io jobs 22619466 and 22619469; scripts `external/schmid_env/fetch_refs.sbatch` and `fetch_refs2.sbatch`). md5,
bytes and source for every file are in its `MANIFEST.tsv`. `run_gate.sbatch` links from the pinned dir and
checks that each file exists.
- `chromosome_arm_positions_grch38.txt`: copied from the CONICS repo at the pinned commit 053cf267, where it ships
  (md5 182816974240dfbcc622dd7739a2dd59). It is not derived. Cross-check against UCSC hg38 cytoBand: chromosome
  ends match exactly, and centromere edges differ by up to 2.97 Mb (a different centromere model). It has 44
  arms (1p-22q), with no X, as in the original.
- `genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz`: the Numbat-docs SourceForge URL. BGZF, passes `gzip -t`.
- `1000G_hg38` (Eagle panel): **deviation.** The Numbat-docs URL (pklab.med.harvard.edu) returned an HTML
  bot-check page, not the zip (kept as `dl_1000G/1000G_hg38.zip.blocked.html`). I did not try to get around it.
  The panel used is a copy of the CRC panel: chr1-22 `genotypes.bcf` and `sites.vcf.gz` with indexes, built
  from `/gpfs/commons/datasets/1000genomes/hg38` (1000G phase 3 GRCh38 lift; hashes in
  `source_inputs.sha256`).
  - The missing chrX does not matter here: their `pileup_and_phase.R` phases chr1-22 only (lines 233 and 269).
  - The remaining difference is how the panel was built from 1000G, which could shift Numbat's phasing slightly
    [Inference].
  - To use the exact file, someone needs to download `1000G_hg38.zip` by hand in a browser and drop it into a
    new ref dir.
- `genetic_map_hg38_withX.txt.gz`: copied from CRC. Its md5 is the same as the Eagle 2.4.1 tarball's copy
  (70f0be81e1ad1314ccc88e48091bcba1).
- `genome_input_BAFExtract_hg38/`: the README Dropbox files (`hg38.list` as `hg38_genome_list`; `hg38.zip`
  unpacked to `dl/hg38`, with `hg38 -> dl/hg38`). The zip has 27 entries: 1-22, X, Y, M `.bin` files plus
  `chr_ids.txt`. No `-preprocess_FASTA` was needed.
- Still missing: the CellRanger outputs per dataset (BAM, BAI, `barcodes.tsv.gz`). See below.

## CellRanger plan (planned, NOT submitted; hold for team-lead)

The script is `benchmarks/schmid/cellranger/run_cellranger.sbatch OUT DS`, with one job per dataset and no
arrays.
- It mirrors the authors' `preprocessing_datasets/cellranger_*.sh`: `--sample=<run>`, the bcl2fastq names they
  used (made as symlinks), and COLO320 `--chemistry=ARC-v1`.
- Module: `module load CellRanger/8.0.1`, with a capital C. The lowercase `cellranger/8.0.1` module is an
  incomplete install: its `count --help` fails on the login node. The capital-C install prints the 8.0.1 usage.
- Reference: `refdata-gex-GRCh38-2020-A` is used read-only from
  `/gpfs/commons/groups/sanjana_lab/jmorris/211215_BeeSTINGseq/downloads/` (`reference.json` version 2020-A,
  GENCODE v32). md5 of `genome.fa`, `genes.gtf` and `reference.json` are in the ref `MANIFEST.tsv`.
- Recorded deviations from their 7.0.0 runs:
  - CellRanger 8.0.1 (amendment f).
  - `--create-bam=true`: 8.0.1 requires the flag (verified in its `--help`); 7.0.0 always wrote a BAM.
  - 16 cores and 120 GB instead of 6 and 64. This changes speed only [Inference].
  - `--include-introns` keeps its default, which is true in both versions.
- The staging step was tested without running CellRanger. BCC06, BCC06post and MCF7 stage correctly. COLO320
  stops because R2 has no `.ok` marker yet (download in progress). MM stops because its FASTQs are absent.

| DS | Runs | FASTQ size | Notes | Est. wall time |
|---|---|---|---|---|
| BCC06 | SRR8315756 | 9.1 GB | R1/R2 complete | 3-8 h [Estimate] |
| BCC06post | SRR8315759 | 3.8 GB | R1/R2 complete | 2-5 h [Estimate] |
| MCF7 | SRR10018149 | 7.5 GB | R1/R2 complete; line-overlap, not primary | 3-8 h [Estimate] |
| COLO320 | SRR12900668 | about 43 GB | ARC-v1 (multiome GEX); R2 still downloading | 15-30 h [Estimate] |
| MM | SRR10574405-12 (8 runs, R1/R2/I1) as RTBM199T L001-L008 | not fetched | needs a FASTQ fetch first | unknown |

- Header for each job: `--account=innovation --partition=cpu --prefer=v6 --cpus-per-task=16 --mem=128G
  --time=90:00:00`.
- Total ask with all 5 running: 80 CPUs and 640 GB, under the 7 TB QOS cap.
- Disk: BAMs are about 1-1.5x the FASTQ size, so roughly 70-100 GB in total [Estimate].
- Submit command: `sbatch -J cr_<DS> --output=<OUT>/logs/%x_%j.out
  benchmarks/schmid/cellranger/run_cellranger.sbatch <OUT> <DS>`, with OUT a new
  `results/schmid_cellranger_<stamp>`.

## Order of work

1. Fix decisions 1-5 and the external reference rule in STATUS, after v-next is frozen and before any Schmid
   file is opened.
2. Build the paper env and resolve the data gaps (MM source, COLO320 file integrity, BCC GEO matrices,
   10x reference).
3. Prepare inputs with their preprocessing scripts. Run the comparators and reproduce their published values
   on one dataset.
4. Run the wrapper on the frozen v-next code. Score with the patched combined scorer, then the two sensitivity
   variants. Report the result whatever it is.
