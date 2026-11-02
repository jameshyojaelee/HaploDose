# Haplodose: status, decisions and open work

Written 2026-11-01. Covers work from 2026-10-31 to 2026-11-01. Plan file:
`~/.claude-teams/plans/t3-infercnv-like-caller-handoff-dapper-cat.md` (milestones M0-M4).
Git: one commit `d6ed692` (2026-11-01 08:35), made while 2 tests were failing.
Current plan: `~/.claude-teams/plans/inspect-the-codebase-clearly-vast-lampson.md`.

## Manuscript direction (user decision, 2026-11-02)

HaploDose is intended to become a full scientific manuscript: a methods paper
with a focused benchmark and a clear account of what RNA can reliably establish
about genomic dosage and tumour clones. The reference for presentation is
Wang et al., *Categorization and analysis of 16 computational methods for
estimating cell potency from single-cell RNA-seq data* (Cell Systems, 2026;
[DOI](https://doi.org/10.1016/j.cels.2026.101762)). That paper is a Synthesis with
targeted reanalyses; borrow its treatment of assumptions, mathematical notation,
implementation and validation limits. This direction does not establish current
readiness, acceptance prospects or a target journal.

The manuscript should distinguish relative RNA dosage, DNA copy-number
alterations and clone identity as separate estimands. Bulk DNA agreement supports
aggregate profile fidelity, not cell-level clone accuracy. Explain which
HaploDose components improve which task through ablations, including cases where
simple gene-ratio or pseudobulk baselines perform as well or better. Existing CNV
benchmarks are documented in [related work](related_work.md); the contribution
must go beyond another broad tool comparison.

Evidence requirements before manuscript claims are finalized:

- Freeze the method, parameters, scoring and primary endpoints before independent
  evaluation. Record datasets already used to select outputs or tune the method
  as development evidence, even when no model parameters were fitted on them.
- Validate at the resolution claimed, using matched DNA where available and
  documenting truth quality, sample correspondence and uncertainty. Synthetic
  implants complement biological validation; they do not replace it.
- Use patients or independently derived cell lines as biological replication
  units, account for repeated samples and shared pools, and specify covariates,
  effect sizes, uncertainty and the multiple-testing correction and family.
- Include simple baselines, fair comparator inputs and post-processing, component
  ablations, reference mismatch, expression-program confounding, rare-clone
  sensitivity and false splits. Report failures and negative results.
- Produce reproducible figures and tables from versioned code, inputs, seeds and
  environment inventories. Resolve material numerical and node-dependent
  behavior before making stability or reproducibility claims; compare resource
  use under controlled conditions.

Current results below are development evidence from multiple code versions, not
a frozen manuscript result set. Broad superiority and reliable subclone recovery
are not established. Bulk and scRNA-seq remain the primary benchmarks; spatial
RNA-seq provides a secondary evaluation. The package and manuscript scope remain
general-purpose; the CRC/IGF2 motivation is a research question, not an established
mechanism.

## Current state (read this first)

- **2026-11-02: 61 Python tests pass on v6 (jobs 22601385, 22601594; 3 Rust in 22601012); Potts fails on v4** (node-dependent,
  see the Kinker subclone section). Repo `clones.py` = ablation variant C (shared tumour doses, Jaccard event
  dedup, all-pairs split test).
- Kinker vs the field: subclone design v3 Haplodose 0.167-0.169 vs infercnvpy 0.165 ARI (tie; region found
  0.36-0.40 vs 0.67); single-line rho 0.854 (log_ratio_gene 0.863) vs pseudobulk 0.863. CopyKAT, SCEVAN, R inferCNV and Numbat runs
  in progress.
- Published CRC results (`~/CRC/results/haplodose_benchmark/{SMC_full,Moorman_full}/`)
  are from jobs 22589677 and 22589678 (2026-10-31). They predate ribosomal-gene
  exclusion, alleles, spatial mode, depth-matched tests, vetoes and the
  cell-type mismatch filter. All later real-data runs are debug runs under
  `logs/runs/` and are not published.

## Validation against DNA (2026-11-01, code after the event_mixtures NaN fix)

Scorers: `benchmarks/ucan_spatial/`, `benchmarks/moorman/`, `benchmarks/cell_lines/`, `benchmarks/ycc/`.
Comparator runners (HaploDose2): `benchmarks/comparators/`. Datasets: `docs/benchmark_datasets.md`.

**U-CAN Visium HD vs WGS FACETS** (20 patients, one section each; `results/ucan_spatial_20261101_0911/`):
- Spearman RNA vs DNA log2(TCN/ploidy) across arms: median 0.62 (0.19 to 0.93). U0396 and U0991 have no
  tumour clone (all cells neutral). U0438 (IGF2 TD) 0.20 with 0/27 altered arms called.
- Discrete arm calls: 147/371 DNA-altered arms in the right direction; 34/423 neutral arms called altered.
- Continuous profile: right sign on 287/336 altered arms. Arm-level RNA SD on DNA-neutral arms 0.085 vs
  median +0.10 on gained arms. The noise floor is per-gene baseline mismatch, not sampling.
- Subclones: 109/166 subclonal DNA segments (CF < 0.8 x purity, >= 10 Mb) overlap a tumour clone split,
  vs 95.5 by chance (random placement). Subclone detection is not better than chance. 21 to 36 clones in
  U0007, U0046, U2955 (over-splitting).
**Moorman vs MSK-IMPACT** (9 same-sample pairs, 8 patients; `results/moorman_20261101_0911/`): median rho 0.76;
  75/99 altered arms right; 77/252 neutral arms called altered (DNA threshold |log2| 0.2, not purity adjusted).
**Cell lines GSE149224**: GEO SW480 labels do not match DepMap SW480 (expression Spearman 0.07; HCT116 and RKO
  0.42, 0.43), so only HCT116 (rho 0.25, same as DepMap's own bulk RNA ceiling 0.29), batch/dose ARI (<= 0.003)
  and line separation (median ARI 0.89) are usable.
**Fix**: `event_mixtures` gave 0 x log 0 = NaN for cells with no expected UMIs in a region, silently dropping
  region split candidates (Moorman). Now `xlogy`; tests green (job 22594015). Effect on all benchmarks: none
  material (scale subclone ARI still 0.19).

## Head-to-head and reference experiments (2026-11-01)

U-CAN Visium HD vs FACETS, 20 patients (comparators by HaploDose2, `results/comparators_20261101_0915/`,
table `collected_20261101_102632.tsv`; malignant/clone calls for both inferCNVs are tool + our rule):

| Tool | median rho | no result | altered right | neutral called | subclonal splits (chance) | median time / RSS |
|---|---|---|---|---|---|---|
| Haplodose | 0.62 | 2 | 147/371 | 34/423 | 109 (95.5)/166 | 84 s / 1.9 GB |
| infercnvpy | 0.62 | 1 | 168/355 | 32/416 | 121 (111.3) | 51 s / 1.0 GB |
| T3 prototype | 0.61 | 0 | 21/331 | 3/386 | 57 (32.7), cn_state runs | 40 s |
| R inferCNV | 0.64 | 13 (16/27 jobs failed) | 66/129 | 62/187 | 16 (15.4) | 678 s / 4.5 GB |

Paired per-patient rho, Haplodose vs infercnvpy / prototype / R inferCNV: higher in 12/18, 14/18, 6/6;
Wilcoxon p 0.23, 0.021, 0.031; Holm over these 3: none below 0.05. Moorman same-sample median rho 0.76-0.79
for all tools. Numbat on YCC: 10-34 min per sample vs Haplodose with alleles 1-6 min.

Experiments on the 20 U-CAN patients (none used for fitting anything):
- **External normal colon epithelium reference** (bin2cell cells from JGCRC4NT, JGCRC14NT mucosa, P3_NAT,
  P5_NAT; 40 um tiles; `benchmarks/ucan_spatial/external_ref/`, `results/ucan_spatial_extref_20261101_1156/`):
  tumour clone found in 20/20 (was 18); median rho 0.71 vs 0.62 on 18 common patients; altered arms right
  200/364 vs 147/371; neutral called 46/416 vs 34/423; tau 0.227 vs 0.232; paired Wilcoxon p 0.23 (10/18
  higher). Caveats: different segmentation (bin2cell), platform/sequencer batch, 5-14x deeper tiles.
- **infercnvpy with the same external reference** (`results/comparators_extref_20261101_1228/`): median rho
  0.64 (0.62 with its baseline reference; higher in 14/19, Wilcoxon p 0.029, nominal). Haplodose vs infercnvpy,
  both with the external reference: median 0.71 vs 0.64, Haplodose higher in 14/20, p 0.18; altered arms right
  200/364 vs 187/357; neutral called 46/416 vs 44/413. Haplodose still worse on U0247 (0.09 vs 0.41) and U0438
  (IGF2 TD; 0.35 vs 0.75).
- **U0438 diagnosis: in-section immune/stromal references leak tumour RNA.** A plain pool of all query
  bins vs the external epithelium scores median rho 0.72 (U0438 0.88 vs Haplodose 0.35;
  `results/ucan_pool_diagnostic_20261101_1329/`), so cell selection was not the loss. The in-section
  immune and stromal bins' own profiles track the patient's tumour DNA (median rho 0.43 and 0.54; U0438
  immune 0.86), and U0438's fitted baseline was 98% immune, which cancels the tumour signal
  (`results/ucan_reference_leak_20261101_1341/`). Also: the depth slope learned on query bins tracks DNA
  (median rho 0.39 vs -0.16 on reference bins; `results/ucan_depth_slope_20261101_1331/`), since deeper
  Visium HD bins hold more tumour.
- **Epithelium-only reference** (in-section immune/stromal set to role other;
  `results/ucan_spatial_epionly_20261101_1342/`, infercnvpy on the same inputs
  `results/comparators_epionly_20261101_1342/`): Haplodose median rho 0.75 (U0438 0.89, U2964 0.91),
  altered right 237/364, neutral called 60/416; infercnvpy 0.63, 213/357, 77/413. Haplodose higher in 15/20,
  p 0.033 nominal (Holm over the 4 comparisons run on these sets: 0.13). Versus Haplodose with in-section
  references: 10/20 higher, p 0.87. Adding `--depth-source reference` (new flag, default query): 0.76,
  11/20 vs without, p 0.67; not made default. Recommendation for Visium HD: no in-section reference groups.
- **Excluding 1,682 low-dosage genes** (U-CAN/TCGA bulk slope < 0.25; `results/dosage_response_20261101_0931/`,
  slopes replicate U-CAN vs TCGA r 0.83): worse, median rho 0.60, altered right 128/371. Not adopted.

## Bulk RNA mode and U-CAN bulk benchmark (2026-11-01)

User direction (2026-11-01): bulk RNA and scRNA-seq are the main benchmarks; Visium HD is secondary; datasets
outside the in-house CRC sets are allowed.

`haplodose bulk` (`python/haplodose/bulk.py`, test `tests/test_bulk.py`): one profile per sample against a
normal panel. Pooled-normal rate, per-bin normal variance, normal-sample programmes high-passed within arms,
beta start chosen by HMM likelihood on a grid (fit_pooled alone never moved beta off 0.6 in low-purity
samples), arm calls as z against leave-one-out normals scaled by sample tau (|z| > 3). Reported log ratios
are centred by the plain median over bins.

U-CAN (`benchmarks/ucan_bulk/`): 1,062 tumours + 120 normals (exon counts, 18,522 protein-coding genes),
truth WGS FACETS per arm log2(TCN/ploidy); 1,034 tumours scored (28 FACETS purity/ploidy NA). Unit =
patient; Holm over method pairs. `results/ucan_bulk_20261101_1412/`:

| Method | median per-tumour rho | median arm AUC gain | median arm AUC loss |
|---|---|---|---|
| Haplodose segmented dose | 0.649 | 0.851 | 0.846 |
| Haplodose log ratio (median-centred) | 0.627 | 0.887 | 0.904 |
| naive gene log2 CPM ratio, arm median | 0.630 | 0.887 | 0.907 |
| infercnvpy (samples as cells) | 0.565 | 0.810 | 0.773 |

Dose beats naive per tumour (589/1,034 higher, Holm p < 1e-4); both beat infercnvpy. Diagnostics: per-gene
median bins, no programmes, mean-of-logs reference and 1-CPM pseudocount changed nothing material
(`results/ucan_bulk_variants_*`, `results/ucan_bulk_diag_ref_*`); the AUC gap (0.854 loss before) was the
HMM-fit centre (`results/ucan_bulk_diag_centre_*`), now replaced by the plain median for reported values.
TCGA pan-cancer (HaploDose2, `benchmarks/tcga_bulk/`, `results/tcga_bulk_20261101_1405/score_centred/`):
16 types with >= 10 solid normals, ABSOLUTE truth, 6,243 paired tumours. Median per-tumour rho / arm AUC
gain / loss: Haplodose log ratio 0.637 / 0.909 / 0.890; Haplodose dose 0.638 / 0.889 / 0.845; naive 0.670 /
0.926 / 0.915; infercnvpy 0.532 / 0.823 / 0.753. Naive > Haplodose pooled (Holm p < 1e-4) and in 12/16 types
(Holm < 0.05 within type); not significant in KICH, KIRP, READ, THCA. Raw STAR counts vs Xena pseudo-counts
(COAD+READ, n 527) change nothing (0.732 vs 0.737). Conclusion: at arm level, Haplodose bulk does not beat
a gene-level median log ratio.
Sub-arm (HaploDose2, `benchmarks/tcga_bulk/subarm.py`, `results/tcga_bulk_20261101_1405/subarm_20261101_1431/`):
gene-anchored ~10 Mb windows, DNA at hg19 gene positions (Ensembl GRCh37 r75). Pooled medians, dose / HD
log ratio / naive / naive smoothed (101 genes) / infercnvpy: within-arm rho (arm medians removed, n 5,055)
0.335 / 0.257 / 0.253 / 0.243 / 0.181; focal-loss recall (mean, n 3,013) 0.53 / 0.46 / 0.47 / 0.48 / 0.27;
focal-gain recall 0.58 / 0.59 / 0.59 / 0.51 / 0.50. Dose vs each method Holm < 1e-4 on within-arm rho and
loss recall; within-arm rho dose > naive in 12/16 types, Holm < 1e-4 in 9 (BRCA, COAD, HNSC, KIRC, KIRP,
LIHC, LUAD, READ, UCEC); naive higher in STAD (Holm < 1e-4). Thresholds at the 5th/95th percentile of DNA-neutral
windows. FPR-matched rerun (seeded random tie-break; `score_fprmatched_random_20261101_1442/`, pooled realised
FPR equal to target for all methods): focal-loss recall dose vs naive 0.254 vs 0.174 at 1% FPR, 0.379 vs
0.276 at 2%, 0.531 vs 0.472 at 5% (Holm < 1e-4, n 3,013, dose higher in 1,349, naive in 871), 0.603 vs 0.630
at 10%, 0.689 vs 0.781 at 20% (10-20% descriptive, not tested). Gains: no difference at any FPR (0.580 vs
0.589 at 5%). Segmentation gives cleaner strict-threshold loss calls; at loose thresholds the HMM's zeroed
weak events cap its recall.

Andor 2020 SNU-601, full-depth scDNA truth (HaploDose2, `results/andor_20261101_1224/scdna_fulldepth/SNU-601`):
1,527 cells, 15 DNA clones, 19/39 arms differ between some clone pair (chance precision 0.49). RNA clones /
difference AUROC / arms called / precision / recall: Haplodose 5 / 0.66 / 13 / 0.77 / 0.53; infercnvpy 6 /
0.77 / 24 / 0.50 / 0.63; R inferCNV 5 / 0.84 / 16 / 0.56 / 0.47; prototype 4 / 0.65 / 6 / 0.67 / 0.21.
One line, 39 arms: no ranking is supported. Other 8 lines (10 specimens) have only low-depth cell-tagged scDNA (median 4.1k fragments per cell). Simulated depth check on SNU-601 (`benchmarks/andor/simulate_lowdepth.py`, seed 20261101; `depth_check_SNU-601/compare/`): low-depth truth gives 2 clones vs 15, ARI 0.12, difference-arm recall 2/19. The 8-line table (`scores_combined/lowdepth_pooled.tsv`; 11 positives; precision 0.06-0.10 for all tools, a lower bound) cannot rank tools; no tests run. Andor result = SNU-601 only; 8-line table is a supplement sanity check.

Andor rerun on repo code 2026-11-02 (clones.py variant C + HaploDose2 fix b; frozen `dev/hd2_20261101/repo_0042`;
`results/andor_rerun_20261102_0042/`), vs the 2026-11-01 12:24 code (all changes since then, not one fix).
SNU-601 full depth (1 specimen, no test): RNA clones 5 -> 12 (DNA 15), difference AUROC 0.66 -> 0.82, arms called
13 -> 22, true 10 -> 17, precision 0.77 -> 0.77, recall 10/19 = 0.53 -> 17/19 = 0.89, profile rho 0.77 -> 0.78.
Comparators (AUROC / precision / recall): infercnvpy 0.77 / 0.50 / 0.63, R inferCNV 0.84 / 0.56 / 0.47; Haplodose
now has the best precision and recall, R inferCNV a slightly higher AUROC. Low depth, 8 lines (supplement):
pooled arms called 54 -> 97, true 4 -> 10 of 11 positives, precision 0.07 -> 0.10, median AUROC 0.89 (n 7) ->
0.97 (n 8). Cost: over-splitting of DNA-clonal lines (MKN-45 1 -> 3 RNA clones, 2 false arms; NUGC-4 1 -> 4, 6
false arms; HGC-27 4-6 clones vs DNA 2; SNU-16 6-7 vs DNA 2). Open: split calibration (see Kinker subclone).

Allele work (WGD, LOH from bulk RNA): 246,275 common exonic SNPs (1000G GRCh38, EUR AF 0.05-0.95,
biallelic; `results/ucan_allele_sites_20261101_1403/`), pileup `benchmarks/ucan_bulk/alleles/pileup.sh`.
Streaming mode (`-T`) 3.6 min / 0.5 GB per sample vs > 29 min indexed (`-R`); full cohort 1,182 samples,
`results/ucan_allele_pileup_20261101_1518/` (jobs 22596652-57). Arm minor-allele fraction without germline
genotypes: `alleles/arm_baf.py` (HWE genotype mixture from EUR AF, beta-binomial het sites, per-arm b on a
grid, per-sample overdispersion). Two test tumours: per-tumour rho of b vs FACETS purity-mixed expected minor
fraction 0.91 (U0001) and 0.86 (UM052); every FACETS LOH arm imbalanced. Scoring (`alleles/score_alleles.py`:
LOH and imbalance AUC per arm, WGD classifier allele vs dosage, IGF2 TD vs non-TD), job 22596673).
Results (`results/ucan_alleles_20261101_1523/score_20261101_185027/`; 1,034 tumours, 28 skipped; unit = patient):
- Per arm vs FACETS: median per-tumour rho 0.73; median per-arm AUC 0.93 for LOH, 0.88 for any imbalance. At a
  call threshold set to the 99th percentile of normal-sample arms, 289/23,740 DNA-balanced arms are called
  (1.2%) and 3,600/16,555 DNA-imbalanced arms (22%): strict calls catch only clear imbalance.
- WGD (out-of-fold AUC, logistic regression, 5-fold x 10 repeats, seed 0): label FACETS ploidy >= 2.7 (295
  WGD): dosage features 0.81, allele features 0.88, both 0.89; label major copy >= 2 over > 50% of genome (328):
  0.83 / 0.90 / 0.92. Adding allele features to dosage: +0.06 to +0.10 (95% paired bootstrap over patients)
  and +0.07 to +0.11. Alleles carry the WGD signal that dosage alone misses.
- IGF2 TD (78) vs non-TD (956), Holm within each source: RNA broad gain chr2 9.0% vs 1.8% (Holm p 0.008);
  RNA arm log ratio chr2 and chr20 higher in TD (Holm p < 1e-4 each); RNA WGD score 0.21 vs 0.13 (nominal p
  0.016, Holm p 0.078, not significant). DNA for comparison: gain chr2 27% vs 10% and chr20 64% vs 38% (Holm
  p < 1e-3 each), FACETS ploidy Holm p 0.13. RNA reproduces the slide-14 chr2 and chr20 gain enrichment; the
  WGD difference is in the same direction but not significant in RNA or in DNA ploidy at n = 78 TD.

GSE149224 HCT116 arm rho (Haplodose 0.25, lowest of four) is truth-limited, not a method ranking. Truth is
DepMap line minus RKO: only 7/39 arms differ for HCT116 (|diff| >= 0.2), so rho is mostly rank noise on
neutral arms; Haplodose gets 5 of the 7 (9q, 8q, 7q, 10p down, 17q up) and misses chr20. For SW480 (22
informative arms), all four tools show no 7p, 3q, 13q gain (DepMap +0.6 to +0.74 vs RKO) and none sees 18q
loss; the tools agree with each other and not with DepMap, so the GSE149224 stocks likely differ from DepMap's
[Inference]. Not used for ranking from here on.

Kinker 2020 (GSE157220; `benchmarks/kinker/`): 53,513 cells, 198 lines in 8 pools + an 8-line custom pool.
GEO has linear CPM and no cell labels: counts recovered exactly (library = 1e6 / min nonzero CPM, max
deviation 2e-11); lines assigned per cell by correlation with DepMap expression of the pool's lines
(`results/kinker_input_20261101_1535/`). 181/198 lines kept (`prepare.eligible`): 7 lines lack DepMap
expression and their cells went to neighbours (custom pool dropped; pool 10 NCIH2077 -> NCIH1048, RD), so
lines with assigned cells off the published count by > 10% are dropped. Design (`results/kinker_input_v2_20261101_1551/`):
two rotations per pool, half the lines as reference, each other line a single-line specimen; truth = DepMap CN
minus the UMI-weighted reference-line CN. Jobs 22596804 (Haplodose) and 22596806 (infercnvpy per pool);
`results/kinker_20261101_1536/` and `kinker_infercnvpy_20261101_1537/` are from cancelled runs on the
uncleaned design and should not be used.
Results (`results/kinker_compare_20261101_*/`, `benchmarks/kinker/compare.py`; n = 181 lines, one query each,
paired Wilcoxon on arm-level Spearman rho, Holm over 3 tool pairs): median rho Haplodose 0.854, pseudobulk
0.863, infercnvpy 0.813; sign on informative arms (|truth| >= 0.2, median 15 per line) 0.964, 0.967, 0.931.
Haplodose > infercnvpy on 130/181 lines (Holm p 2e-11); pseudobulk > Haplodose on 119/181 (Holm p 2e-6, median
gap 0.01, same in single- and multi-clone specimens). A single-line specimen is all one population, so its
pseudobulk is the best possible estimate; Haplodose cannot beat it here, only match it. Haplodose gives 1 tumour
clone in 115/181 specimens, 2-8 in 66 (Kinker report CNA subpopulations in some lines [Unverified], so not all
splits are errors).

Kinker subclone design (`benchmarks/kinker/implant.py`, `score_implant.py`, `run_implant.sbatch`): in each
single-line query (166 lines, >= 100 cells) one event is implanted in a random share of cells (binomial thinning;
size 50/150/500 expressed genes x loss/gain x 15%/40%, 13-14 lines per setting). Same expression state, so only
the event separates the two groups. infercnvpy runs with `--rule all` (every query cell clustered; the
held-out-reference malignancy null is meaningless when the reference is other cancer lines). Mean ARI vs the
implant split, unit = line, paired Wilcoxon, Holm over 1 pair: before the review fixes Haplodose 0.126 vs
infercnvpy 0.157 (p 0.086; `results/kinker_implant_20261101_1619/`), after 0.135 vs 0.157 (p 0.19;
`results/kinker_implant_v2code_20261101_1717/`). Haplodose leads only on 500-gene losses (0.83 vs 0.67 at 40%,
0.47 vs 0.30 at 15%); it finds almost nothing at 50-150 genes or for gains, where infercnvpy reaches 0.07-0.33.
[Inference, later found incomplete] Event proposals are whole chromosomes or halves, so a 15-bin event is diluted in its region.

Diagnosis (2026-11-01 evening). (1) Region ranking: the same chr11 (FTH1, MALAT1) and chr16 (MT2A) regions
topped the event proposals in every specimen; region totals with one median dispersion let one overdispersed
gene dominate. `event_mixtures` now uses per-bin NB likelihoods with a shared rate ratio, a per-cell NB size
factor, sliding windows (`CloneParams.event_windows` = (8, 24) proposal bins; CLI `--event-windows`) and drops
candidates whose cell split agrees >= 90% with a better one; `tests/test_numerics.py` has a test that fails on
the old code. Result (`results/kinker_implant_windows_20261101_1854/`): proposals now hit the implanted region
(TOV21G chr12, SKMEL30 chr3), but ARI fell to 0.116 vs infercnvpy 0.157 (Holm p 0.003); 500-gene losses lost
0.115 ARI on average from more accepted false splits that then starve a small true subclone below min_cells
(CAOV3: chr12 and chr10 splits at score ~10 before the true chr11 split at 48). Single-line: rho 0.854 and
115/181 one-clone lines (116 before). Scale 100k: 404 s (was 148 s), ARI unchanged 0.19. (2) Gain design flaw:
v2 gains thinned implanted cells outside the region to 2/3, so they had a third fewer UMIs; Haplodose models
depth as a covariate and reads that as depth. Design v3 (`implant.py` GAIN_MODE depth_neutral, same regions and
cells as v2) thins the whole specimen to 2/3 except the region in implanted cells. Reruns with and without
windows: jobs 22598244 / 22598245, `results/kinker_implant_v3[_nowin]_20261101_2021/`. Open: the split
acceptance test is not calibrated (review: even/odd reuse), which more candidates make worse.

v3 results (n = 166 lines, paired Wilcoxon, Holm over 1 pair): windows on 0.147 vs infercnvpy 0.165 (p 0.087),
off 0.153 (p 0.18); windows are now off by default (`event_windows = ()`). Cause of the remaining loss (MCF7,
15-bin loss in 57/143 cells): EM drained the small clone because each clone was segmented on its own and a
small clone's profile is more shrunk. Fix (`SpecimenModel.shared_profiles`, used in `em`): tumour clones take
the pooled tumour dose plus their weighted mean deviation from the pool on each run where the pooled path, own
path and resolved pairwise differences are constant. Kept version (deviation only on resolved differences):
0.169 vs 0.165 (p 0.93; `results/kinker_implant_v3_shared_20261101_2109/`), region found 0.40 vs 0.67. Also kept:
event dedup is Jaccard >= 0.7 on minority cell sets (raw agreement >= 0.9 dropped disjoint rare groups), and
`try_splits` requires every pair of parts to differ, not only the first two. Ablation on frozen copies, all on
v4 (`dev/hd1_ablate_20261101_2350/`; 166 lines implant, 181 single-line):

| variant | implant ARI | region found | single-line rho | one-clone lines | tests |
|---|---|---|---|---|---|
| A shared doses | 0.169 | 0.37 | 0.854 | 70 | 58/59 |
| B A + deviation on every run | 0.147 | 0.39 | 0.854 | 75 | 57/59 |
| C A + Jaccard + all-pairs (repo) | 0.167 | 0.36 | 0.854 | 70 | 58/59 (v4), 59/59 (v6) |
| D B + C | 0.149 | 0.38 | 0.854 | 76 | 57/59 |

infercnvpy: 0.165, region found 0.67. Where implanted events are lost (variant C, `trace_region.py`): an event
proposal overlaps the implanted region in 49% of specimens, and most of those are merged back by the odd-bin
agglomeration; 22% accepted. Tried, not adopted: variant E, a targeted held-out test of each event proposal's
region (robust z vs same-size windows, |z| >= 4, gene-consistency veto), also guarding merge_unsupported:
ARI 0.173 (noise-level vs C), one-clone lines 70 -> 61, rejected most proposals that hit real implants
(150-gene 40% gains: 0.64), and failed `test_expression_programme_is_not_a_clone` (no programme correction; a
co-regulated cluster shifts even and odd bins alike). Reading: the proposed cell groups are too mixed for any
held-out test to have power; proposal quality, not the acceptance test, limits small events.
Proposal diagnostic (`dev/hd1_proposals_20261102_0600/`, r0, 150- and 500-gene settings, 51 specimens, up to 50
event candidates kept): Leiden proposals have ARI ~0 with the implanted cells; an event proposal overlapping a
150-gene region exists in nearly every specimen, but the best one has ARI 0.15-0.29 (windows on: 0.26-0.51) and
ranks 3-9 (5 are kept). Fate follows quality: ARI > 0.6 accepted 11/13, ARI < 0.1 mostly merged. Oracle bound
(`oracle.py`; region known, 2-component binomial mixture on each cell's region UMIs, fitted without labels, real
counts): ARI 0.32-0.63 at 50 genes, 0.58-0.84 at 150, 0.80-0.95 at 500, vs Haplodose 0.00-0.34 and infercnvpy
0.00-0.32 at 50-150. So there is large headroom, and it is in localising the region. Tried, not adopted:
variant F, windows (8, 24) on and candidates ranked by gain as a robust z within their kind (half, whole,
window size): ARI 0.171 (noise vs C 0.167), better rank for 150-gene 40% gains (1.3), worse for losses (10-17),
500-gene 15% gains 0.094 -> 0.016. Same oracle on Haplodose's proposal genes only (used genes in even bins,
`oracle_even.tsv`): 150 genes 0.40-0.62 (median region UMIs ~60 vs ~190 with all genes), 500 genes 0.52-0.89.
So the even/odd split plus the gene filter costs ~0.2 ARI at 150 genes, and localisation another ~0.1-0.2.
The even/odd design is what makes the odd-bin test valid (labels chosen on even bins); using all bins needs a
different calibration, e.g. a null for the best window's mixture gain from a parametric bootstrap.
One-clone lines were 115/181 before shared doses (truth unknown: Kinker
lines can carry real subclones). `test_spatial::test_potts_improves_coherent_clones` fails on v4 for the old
own-profile code too (accuracy 0.75 < 0.85, job 22600052) and passes on v6 (22597830, 22601012): the fixture's
per-cell signal is weak by design (plain run: the copy-neutral clone takes 334 tumour cells; the chr12 split
proposal separates the subclones barely above chance, `dev/hd1_potts_20261101_2250/dbg_chr12_22600993.out`), so
the outcome flips with CPU float differences. Open. A/B noise: identical code on v4 vs v5 differed by 0.004 ARI
(HaploDose2); `--constraint` on the sbatch line does not override `--prefer` in the script, so pin in the script.

Proposal floor (#5). Field rescore by floor stratum (HaploDose2): at/above 30 implanted cells (137) Haplodose
0.192 vs infercnvpy 0.180 ARI; below (29) 0.051 vs 0.093, so the gap to infercnvpy is in the below-floor
specimens. `--min-cells` experiment on repo code (clones H, v4, `dev/hd1_ablate_20261101_2350/out_R1_*_mc{15,20}`):
min-cells 30 / 20 / 15: below-floor ARI 0.031 / 0.081 / 0.150, at/above 0.195 / 0.183 / 0.184, all 166 0.167 /
0.165 / 0.178 (mc15 vs 30: 58/64, p 0.80); single-line one-clone lines 70 / 46 / 38 of 181, gene rho 0.860 in
all. Split scores of small (< 30 cell) splits do not separate helpful ones (implant below-floor, ARI +0.1:
scores 5-25) from single-line small splits (median 11.6, 90th percentile 21.4), so a size-dependent score
threshold is not justified without a calibrated null. Default stays 30; `--min-cells 15` is the option when rare
subclones matter more than extra splits.

Single-line profile gap to pseudobulk (HaploDose2, `results/singleline_decomp_20261101_2354/`): summed-count
bins let a few highly expressed genes dominate a bin (about 0.01 rho); the median of per-gene ratios within the
same bins recovers it (0.862 vs pseudobulk 0.863). Added (HaploDose2, add-only): `profiles.tsv` log_ratio_gene
and `arms.tsv` median_log_ratio_gene (per-gene log ratio, genes with E_g >= 1 in the clone, plain median in bin
or arm; bulk.py not touched); existing columns unchanged. Kinker single-line rho 0.854 -> 0.863 (181 lines,
109/72, Holm p 0.001 over 2 gates), but the column was chosen after the Kinker decomposition, so that is not
independent evidence. Independent check, Andor SNU-601 full depth (1 specimen, no test;
`results/andor_gene_20261102_0113/`): profile rho 0.765 -> 0.815, difference AUROC 0.76 -> 0.87; low depth
profile rho median 0.71 -> 0.77 (9/10 specimens higher). Costs: U-CAN spatial rho 0.64 -> 0.63 (20 patients,
p 0.52), small clones' region contrast 0.419 -> 0.411 (65 clones < 100 cells, p 0.009). Tests 61/61 on v6.
Andor numbers move with node type (SNU-601 v4: 11 clones, AUROC 0.76; v6: 12, 0.82, same code): quote ranges.
Remaining gap to pseudobulk (field run, 181 lines; d = rho gene - rho pseudobulk): no clone-count effect
(pseudobulk ahead in 39/71 one-clone vs 67/110 multi-clone lines, Fisher p 0.44; median d -0.004 in both); a
small shift over all lines, about 0.004 median, of which about 0.0015 is bin coarsening (pb_gene 0.8634 vs
pb_binmed 0.8619). Not pursued.

Kinker field comparison 2026-11-02 (`results/kinker_field_20261102_0529/`; Haplodose = repo code, clones.py
variant C; unit = cell line; paired Wilcoxon, Holm over the pairs against Haplodose (single-line) or over all
tool pairs (subclone)). Single-line, 181 lines, median rho vs DepMap arms / mean sign on informative arms:
Haplodose log_ratio_gene 0.861 / 0.966, Haplodose bins 0.854 / 0.964, pseudobulk 0.863 / 0.967, infercnvpy
0.813 / 0.931, CopyKAT 0.765 / 0.915, SCEVAN 0.710 / 0.896. Haplodose beats infercnvpy (138/43 lines, Holm p
9e-15), CopyKAT (169/12, 7e-29) and SCEVAN (174/7, 5e-30); pseudobulk is slightly ahead (106/74, Holm p 0.03).
Subclone v3, 166 lines, mean ARI / region found: Haplodose 0.167 / 0.36 (v4 run `out_C_implant`), infercnvpy
0.165 / 0.67, CopyKAT 0.089 / 0.48, SCEVAN 0.027 / 0.17. Haplodose vs infercnvpy p 0.91 (tie); Haplodose beats
CopyKAT (Holm p 0.0015) and SCEVAN (7.5e-6). Caveats: CopyKAT has no subclone output, so its clones are our
Leiden post-processing (`cnv_post`); SCEVAN's clones are its own (`--clones native`), fitted per pool, and it
called a median of 1 clone per line. All tools got other cancer lines as the normal reference (none was built
for that). infercnvpy single-line scores are from the 2026-11-01 16:09 evaluate run on the same input.
R inferCNV (job 22598581, `results/kinker_infercnv_r_20261101_2110/`; comparison with all tools in `results/kinker_field_20261102_0529/single_line_with_infercnv_r/`): single-line median rho 0.835, sign 0.954; Haplodose higher in 118/181 lines (Holm p 1.4e-6 over 5 pairs). Subclone design v3, all five tools with the 2026-11-02 scorer (`results/kinker_field_20261102_0529/implant_all5/`, ari_all, 166 lines, paired Wilcoxon, Holm over 10 tool pairs): Haplodose 0.167, infercnvpy 0.165, R inferCNV 0.163, CopyKAT 0.089, SCEVAN 0.027. Haplodose vs infercnvpy p 0.91, vs R inferCNV 0.51 (ties); Haplodose beats CopyKAT (Holm 0.0025) and SCEVAN (1.2e-5). By floor stratum (>= 30 / < 30 implanted cells; 137 / 29): Haplodose 0.192 / 0.051, infercnvpy 0.180 / 0.093, R inferCNV 0.175 / 0.107. Effect sizes for single-line (median paired rho difference, 95% bootstrap CI): vs R inferCNV +0.018 [0.009, 0.026], infercnvpy +0.033 [0.025, 0.048], CopyKAT +0.095 [0.083, 0.110], SCEVAN +0.144 [0.132, 0.161], pseudobulk -0.004 [-0.009, -0.001] (`single_line_ci/`). These are development-benchmark numbers; confirmatory claims follow `docs/heldout_protocol.md`.

Compute (2026-11-02). Rough timings from the Kinker runs, not a controlled comparison (different nodes, CPU
counts and batching): Haplodose 46 min for both rotations (2 runs of ~53k cells, 8 CPUs, v3 node), 3.1 GiB;
infercnvpy 13 min (16 per-pool runs), 4.6 GiB; SCEVAN 6.3 h, 16 GiB; CopyKAT 9.2 h, 22 GiB (both 16 CPUs, v4).
Controlled scaling benchmark: `results/scaling_20261102_0819/` (`benchmarks/scaling/`: prepare.py, run_one.sbatch,
summarise.py). Real Kinker r0 query and reference cells (45,977 available): 5k / 20k without replacement, 50k
(30,431 unique) / 100k / 1M with replacement (duplicates renamed); 3 repeats up to 50k, 1 at 100k, 1M for
Haplodose and infercnvpy only; every run on cpu,
`--constraint=v6`, 16 CPUs, 72 h cap, 200G memory cap (SCEVAN 5k and 20k first ran with 32G / 96G and ran out of
memory; rerun at 200G as repeats 11-13). Wall time from /usr/bin/time -v; peak memory = max(sacct MaxRSS of the
batch step, /usr/bin/time), since /usr/bin/time misses forked R workers and sacct samples every 30 s. Haplodose
code frozen in `code/` (08:19). Jobs in `jobs.txt`.

Field comparators (env `cnv_field`, `slurm/build_cnv_field.sbatch`: CopyKAT 1.2.5 @01924c0, SCEVAN 1.0.3
@5a49b88, R 4.4; runner `benchmarks/comparators/run_r_field.py`, `benchmarks/kinker/run_field.sbatch`). Both use
the fit reference lines as normal cells, as for inferCNV. SCEVAN's `plotConsensusCNA` fails with >= 10
subclones (grep matches subclone1 and subclone10); `scevan_run.R` replaces it with a no-op. Runs:
`results/kinker_{copykat,scevan}_20261101_2246/`, `results/kinker_implant_v3_{copykat,scevan}_20261101_2246/`.
R inferCNV: `benchmarks/kinker/run_infercnv_r.sbatch` (pinned `--constraint=v6`: SIGILL on a v4 node).
Numbat on Kinker needs the ENA BAMs (`benchmarks/kinker/alleles/fetch_bams.sbatch`; HTTPS gave 403 pages and
truncations, now FTP with size-checked resume, job 22601406) and an hg19 1000G panel (`build_hg19_panel.sbatch`,
done: `data/external/numbat_hg19_panel/`). Cell IDs are `<barcode>-<channel>-<pool>`; ENA titles are
`Pool<channel>_<pool>`, so each line's barcodes per channel BAM are known. Pileup and phasing must be per line
(each line is a different person). Implanted subclones exist only in thinned counts, not in reads, so Numbat and
allele mode can be scored on the single-line design (vs DepMap) only, not on the subclone test.

## Efficiency review response (2026-11-01, `docs/efficiency_review_2026-11-01.md`)

All five numerical findings reproduced and were fixed; each has a test in `tests/test_numerics.py` that fails on
the old code. PLN log pmf rewritten in deviance form (was +16 at y = mu = 1e7); HMM messages renormalised per step
in scan, associative scan and Viterbi (posteriors moved up to 0.12 under a constant emission offset); programme
score gradient includes the softmax term and runs 5 Fisher steps (3 did not converge; off by 0.29 with unequal
dispersions); event mixtures solve the weighted NB score equation (`clones.nb_rate`) with a convergence stop;
pooled expectation uses UMI-weighted depth mean and variance (`SpecimenModel.depth_log_rate`; worst bin error 2%
vs 24% for the mean alone). Efficiency: Viterbi skipped where paths are unused, sparse depth-matched weights,
agglomeration pair cache, randomized SVD for programmes, power-of-two JAX chunk shapes, sparse bulk binning,
block-aligned gene reads, compact allele reading (new test). Bulk normals are now fully leave-one-out (rate,
variance, programmes); U-CAN bulk not rerun yet. Not done: even/odd independence (needs null simulation),
regional dispersion calibration, neighbour search, bulk beta-grid batching, Rust kernels. The review lists
scikit-learn 1.9.1; the env has 1.7.2. Code snapshot `results/code_snapshot_20261101_1717.tar.gz`.
Effect (new code, jobs 22597228-30): Kinker single-line median rho 0.854 (unchanged); subclone ARI 0.126 -> 0.135;
scale 100k: 148 s, 3.45 GiB, ARI 0.19, normals called malignant 0.20 (unchanged). The fixes are correct but do
not explain the subclone weakness.

## Review response 2026-11-02 (`docs/efficiency_review_2026-11-02.md`)

Split: HaploDose1 #1, #2, #7 (clones.py) and #5 (reporting strata); HaploDose2 #3, #6, #8, #9; #4 a decision (see
Decisions). HaploDose2 part, in the repo; full suite 65/66 on v4 (Potts, v4 only; job 22604409):
- #9 `kernels/estep.py`: cell log-likelihoods are computed relative to clone 0 per bin from float64 log-rate
  differences (`cell_loglik` now returns values up to a per-cell constant; every consumer uses within-cell
  differences). The review's case (+0.01 in 20 of 2,000 bins) returned 0.0; now within 2% of float64 (new tests).
- #3 `pipeline.py`: log_ratio_gene splits each bin's expectation by gene rates mixed with the specimen's fitted
  reference weights (and the in-sample normal type), not raw pooled reference counts. The review's pure-type-A
  query: max |log_ratio_gene| 0.02, was -1.0 (new test). No effect where a pool has one reference group (Kinker).
- #6 `benchmarks/kinker/score_implant.py`: ari is NaN below 2 assigned cells (was 1.0 for an empty specimen);
  primary endpoint ari_all (unassigned cells one extra label, missing specimen 0), coverage, n_no_output. The field
  implant rescore reproduces every tool's mean ARI exactly (coverage 1.0 for all; `results/review2_gates_20261102_0844/implant_rescore`).
- #8 `bulk.py`: median mode rebuilds pseudo counts without each held-out normal; fewer than 10 normals warns and
  sets run.json `arm_calls_calibrated` false (threshold unchanged). `benchmarks/ucan_spatial/dna_truth.py` also
  reports rho_log_ratio_gene.
Gate, U-CAN spatial, 20 patients (18 with tumour clones), v4, base = code of 01:13 (clones C, old E-step, pooled
gene shares) vs new = 08:44 repo (`results/review2_gates_20261102_0844/`; paired Wilcoxon, Holm over 3): bins
rho 0.635 vs 0.635 (1/4 patients higher/lower, Holm 0.16; #9 alone does not move U-CAN); gene rho 0.636 -> 0.649
(10/7, p 0.05, Holm 0.15; U0354 0.47 -> 0.64, U2761 0.62 -> 0.70); new gene vs new bins 0.649 vs 0.635 (10/8,
p 0.50). Neutral arms called 31/423 and altered right 158 vs 156/371 in both. The earlier U-CAN cost of the gene
column (7/11 against) is gone after #3. Kinker effect of #9 (+#3, a no-op with one reference group per pool),
v4, `dev/hd1_ablate_20261101_2350/out_{Hbase,C}_*` (Hbase = 08:40 repo with clones.py C): v3 implant ARI
0.1655 vs 0.1672 (166 lines, 24/17 higher/lower, p 0.92); single-line bins rho 0.8540 vs 0.8540 (13/18, p 0.38),
one-clone lines 72 vs 70. No effect on Kinker.
Also (HaploDose2, on request): `profiles.tsv` assignment_dose (the dose that set the final cell likelihoods,
`fit_specimen` doses; normal-like clone 0) next to log_dose (the clone's own HMM fit, as segments.tsv) for #2;
test in test_pipeline.py (15/15). #5 strata in score_implant.py (`below_floor` = n_implanted < 30, the default
min_cells; all 166 stay in the primary endpoint). Field implant runs (`results/review2_gates_20261102_0844/
implant_rescore_floor`), mean ari_all / region found, descriptive: 137 at or above the floor: Haplodose 0.192 /
0.42, infercnvpy 0.180 / 0.72, CopyKAT 0.102 / 0.51, SCEVAN 0.032 / 0.20; 29 below: Haplodose 0.051 / 0.07,
infercnvpy 0.093 / 0.48, CopyKAT 0.025 / 0.34, SCEVAN 0.004 / 0.03. Haplodose's deficit against infercnvpy sits
in the below-floor specimens; there it almost never reports the region (2/29).

HaploDose1 part, in the repo (clones.py = variant H; new `tests/test_review_20261102.py`):
- #1 `shared_profiles`: deviations are averaged within atomic intervals, cut at every accepted segment boundary of
  every tumour pair and at every pooled-path change. The review's probe (tumour 1 vs 2 opposite events, tumour 3
  spanning both): contrast 0.005 before, now at least 0.9 x the independent fit's 0.099 with opposite signs, and
  invariant to clone order. Test fails on the old code (job 22604390).
- #2 the deviation has the fitted programme contrast removed (as in `differences`), and segments pass the
  gene-consistency veto with MAP labels (`supported(..., log=False)`); alleles are not used there because the
  allele pools are updated after this step in `em`. Reporting: `fit_specimen` returns both `doses` (assignment)
  and `prof` (independent fit); profiles.tsv still writes the independent fit. Open: write the assignment dose
  as its own column.
- #7 `neighbour_graph` removes the cell by index and orders ties by (distance, index): 0 self-edges at
  triplicated coordinates (was 20 in the review's probe; test fails on old code). The directed row-normalized
  graph is documented as a smoothing heuristic, not a Potts mean-field gradient. Spatial event proposals use the
  moment-matched NB size (sum m)^2 / (phi sum m^2) for neighbourhood sums (was n x one cell's size).
Gate, v4, H vs Hbase: implant ARI 0.1668 vs 0.1655 (50/48, p 0.68), region found 0.39 vs 0.37; single-line bins
rho +0.0003 (47/46, p 0.43); gene rho -0.0013 (36/64 lines, p 0.019 nominal, Holm 0.056 over the 3 H tests);
one-clone lines 70 vs 72. Adopted for correctness; the benchmark effect is at most -0.001 rho.

## Package layout

| Path | Content |
|---|---|
| `crates/core/src/` | Rust (PyO3): streaming per-cell bin sums, grouped column sums (`aggregate.rs`) |
| `python/haplodose/io.py` | count readers (10x H5, h5ad, CSR H5), output guard (`HAPLODOSE_PROTECTED_PATHS`), gene and site lists |
| `python/haplodose/bins.py` | fine bins (10 genes), ribosomal exclusion, chromosome and arm boundaries, arm labels |
| `python/haplodose/data/centromeres_hg38.tsv` | p-arm ends from UCSC cytoBand (CRC copy, md5 in file header) |
| `python/haplodose/profile.py` | pooled Poisson-lognormal HMM segmentation, specimen fit (centre, beta, tau), clone difference test |
| `python/haplodose/programmes.py` | normal-cell expression programmes as per-cell nuisance terms |
| `python/haplodose/alleles.py` | phased SNP reading, normal-cell SNP filter, allelic-imbalance HMM, cell allele likelihood, vetoes |
| `python/haplodose/clones.py` | per-specimen clone model: proposals, merges, EM, splits, Potts prior |
| `python/haplodose/pipeline.py` | end-to-end run, reference check, outputs |
| `python/haplodose/kernels/` | JAX: HMM (scan on CPU, associative scan on GPU, per-position transitions), NB/PLN likelihoods, cell E-step, programme scores |
| `python/haplodose/simulate.py` | simulator: clones, depth effects, programmes, alleles with phase switches, spatial bands |
| `python/haplodose/cli.py` | `haplodose run`, `haplodose simulate` |
| `benchmarks/window_grid.py` | programme high-pass grid (done) |
| `benchmarks/scale.py` | streaming large simulation + run, records wall time and peak memory |
| `slurm/build_test.sbatch` | cargo test, editable install, pytest (compute node) |
| `slurm/window_grid.sbatch`, `slurm/scale.sbatch` | benchmarks |
| `slurm/gpu_smoke.sbatch` | prepared, **never submitted** (each GPU run needs explicit approval) |
| `logs/dbg/` | debug scripts and their logs (see last section) |
| `logs/runs/` | debug real-data runs (not published) |
| `results/` | window grid and scale benchmark outputs |

CRC side (outside this repo): `~/CRC/src/44_haplodose/` (adapter: `masks.txt`,
`allele_exclude_genes.txt`, `allele_exclude_sites.txt` with chr11:2130588),
`~/CRC/slurm/haplodose_cpu.sbatch` (staged publish with archive). Journal entries
H.2 and H.3 appended to `~/CRC/docs/notes/2026-10-30_meeting_followups.md`.

## What the model does now

1. **Binning.** Genes are summed into bins of 10 expressed genes in one streaming pass.
   Never features: user mask (CRC: IGF2 and aliases), ribosomal protein genes
   (`bins.RIBOSOMAL`), and genes no reference type matches within log 5
   (`pipeline.type_mismatch_genes`, listed in `type_mismatch_genes.tsv`).
2. **Baseline.** Mixture of reference types fitted to the pooled specimen. Per-bin
   mismatch variance from reference heterogeneity (`SpecimenModel.v`) is added to tau.
   With alleles, cells that are confidently balanced in the specimen's strong
   imbalance runs become an extra in-sample normal reference (`alleles.balanced_cells`).
3. **Programmes.** Top 10 SVD components of normal-cell residuals. Arm means are removed
   from the loadings (default `programme_window=0`), so arm-level copy number cannot be absorbed.
4. **Clone search.**
   - Leiden proposals (direct igraph Leiden; scanpy was removed from this path because
     of a scanpy/scipy `getnnz` bug), merged unless their pooled profiles differ.
   - EM over clones with NB cell likelihood plus depth slopes, programme scores and
     allele evidence.
   - Divisive splits proposed by Leiden, per-region two-component mixtures and per-run
     allele mixtures.
   - Even bins propose and odd bins test; even SNPs propose and odd SNPs test.
5. **Split and merge tests.**
   - Depth-matched pooling (coarsened exact matching on log depth) for split tests,
     merges after EM and reported differences.
   - Gene-consistency veto: one-sided sign test that most genes in a segment move.
   - Allele veto: expression differences with well-powered, unchanged allelic imbalance.
6. **Alleles.**
   - Numbat-format phased counts.
   - HMM over balanced and imbalanced (4 levels × 2 orientations) states.
   - Phase-switch probability (1 − e^(−2νd))/2 with d in cM; ν and the beta-binomial
     κ are chosen by likelihood.
   - SNPs biased in same-patient reference cells are dropped.
   - Joint states per segment, including cnLOH.
7. **Spatial.**
   - Potts prior over an 8-nearest-neighbour graph, mean-field.
   - λ chosen by how well neighbours' data-only posteriors predict each cell's own data.
   - The mean-field state is carried across E-steps and splits.
   - Proposals use neighbourhood-summed features. Cell likelihoods always use the
     cell's own counts.
8. **Outputs.** `cells.tsv.gz`, `posteriors.npy`, `clones.tsv`, `segments.tsv` (with
   allele columns and `joint_state` when alleles are given), `allele_segments.tsv`,
   `differences.tsv`, `profiles.tsv.gz`, `arms.tsv` (new, untested), `tree.tsv`,
   `log.tsv`, `reference_check.tsv`, `specimens.tsv`, `bins.tsv`, `genes.tsv.gz`,
   `type_mismatch_genes.tsv`, `bin_counts.npy`, `qc.json`.

## Decisions and the evidence for each

| Decision | Evidence |
|---|---|
| Expression programmes as nuisance terms | SMC held-out normals called aneuploid 0.121 observed vs 0.063 expected before; 0.050 vs 0.043 after (job 22577293) |
| Arm-mean removal instead of a 25-bin window | Synthetic grid (job 22589519, `results/window_grid_20261031_203023/`, 20 datasets per setting). Mean ARI / tumour clones (truth 2) / false programme differences / normals obs vs exp: off 0.14 / 3.85 / 126 / 0.013 vs 0.004; w25 0.17 / 2.45 / 25 / 0.010 vs 0.004; arm 0.27 / 1.35 / 0 / 0.004 vs 0.005. On SMC, every high-pass option could absorb at most about 5% of any real clone difference (job 22589503) |
| Exclude ribosomal protein genes | Whole-chr19 splits (SMC14/23/25, KG146M2, KG74P, KG136P) were driven by RP genes: SMC14 chr19 UMI-weighted shift +0.28, without RP genes +0.08; SMC18 chr8 control unchanged (job 22593201). Moorman clone ARI vs infercnvpy rose from 0.111 to 0.195 |
| Depth-matched difference tests | YCC155 clones were depth strata (median UMIs 7,700 vs 2,400 vs 1,770) separated by whole chromosomes |
| Exclude cell-type mismatched genes | Malignant agreement with Numbat: YCC155 0.56 → 0.985, YCC230 0.81 → 0.977; YCC230 clones 20 → 11; SMC normals 0.047/0.027 → 0.027/0.020 (jobs 22593864-74) |
| Potts λ from neighbour-predicted likelihood | λ from pseudo-likelihood on smoothed labels ran away (λ = 16 with shuffled coordinates). The current rule passes both spatial tests |
| Per-run allele splits with phase-free test | Summing all runs gave 457/54 splits driven by truncal events; requiring "clearly balanced" in the other group rejected true cnLOH splits |
| HMM transitions keep equal exit rates (no change; review 2026-11-02 #4) | `sticky_transition` uses T_ii = 1-p, T_ij = p q_j/(1-q_i), so the long-run state occupancy is proportional to q_i(1-q_i), not q_i: neutral 0.379 for the copy-number HMM (q 0.2/0.5/0.2/0.1) and 0.452 for the difference HMM (q 0.15/0.7/0.15), against the 0.5 and 0.7 written in STATE_PRIOR/DIFF_PRIOR (recomputed by HaploDose2). q sets where a state change goes, p how often; every state is left at the same rate. Valid as is; switching to an occupancy-preserving redraw would change state exit rates and need recalibration on null false calls, so it is not done to gain sensitivity |

## Tried and reverted

- `cell_doses='resolved'` (only resolved segments in the cell likelihood): lost the
  expression subclone in S1.
- Seven-state (fractional step) difference HMM: more false splits on real data (Moorman
  12 of 19 multi-clone, U0438 24 clones, depth test failed) and did not recover the
  scale subclone.
- Proposal refinement by even-bin reassignment: broke the S1 subclone tests and did not
  change the scale result.

## Latest debug results (unpublished; code as of jobs 22593864-74, last green tests 22593863)

**YCC with alleles vs native Numbat** (`logs/runs/ycc6_1101_0815/`):

| Sample | Malignant agreement | Numbat CNV segments matched | Tumour clones | Clone ARI vs Numbat |
|---|---|---|---|---|
| YCC230 | 0.977 | 15 of 19 | 11 | 0.04 |
| YCC061 | 0.513 | 20 of 26 | 3 | 0 |
| YCC097 | 0.549 | 0 of 2 | 2 | NA |
| YCC155 | 0.985 | 11 of 14 | 5 | 0 |
| YCC162 | 0.957 | 10 of 21 | 1 | 0 |

**SMC and Moorman, expression only** (`logs/runs/*_v6_1101_0815/`):

| Dataset | Multi-clone specimens | Held-out normals, observed vs expected |
|---|---|---|
| SMC | 8 of 22 | 0.027 vs 0.020 |
| Moorman | 9 of 19 | 0.007 vs 0.002 |

**Visium HD, spatial mode** (`logs/runs/spatial6_1101_0815/`). Spearman of arm-level
RNA against patient DNA log ratio, Haplodose arm median vs prototype:

| Sample | Haplodose | Prototype |
|---|---|---|
| U0354 | 0.565 | 0.564 |
| U0440 | 0.610 | 0.638 |
| U0247 | 0.230 | 0.227 |
| U0438 | 0.064 | 0.638 |

**Simulations:**
- **cnLOH clone found only with alleles.** Accuracy 0.83 with about 3 SNP reads per cell
  on the event [Estimate: the best achievable is about 0.89].
- **Spatial.** At 800 UMIs Potts reaches about 0.99 accuracy vs 0.66 without
  (debug job 22593332, earlier code).
- **Scale, 100k cells.** 11,000 genes, 4 specimens, depth 1,500 (job 22593859):
  145 s, 3.35 GiB peak. The subclone (a 250-gene half-chromosome gain) was not found
  (median ARI 0.19). 20% of normals are called malignant [Inference: information limit
  at d′ ≈ 3 with a 90:10 prior; the reference check agrees].

## Open problems (not done)

1. **Failing tests** (see the top section).
2. **Scale subclone missed.** Region mixtures rank large regions above the real event
   (top-5 candidates were chr1, chr11, chr12, chr16; the event was on chr3).
   The per-region moment fix broke S1. Needs a calibrated ranking or more candidates.
3. **Over-splitting where references do not match the query** (YCC230 11 clones vs
   Numbat 4; clone ARI vs Numbat about 0 everywhere).
4. **YCC061 and YCC097 malignant calls disagree with Numbat** (0.51, 0.55).
   No in-sample normals were found there.
5. **U0438 spatial profile does not track DNA** (0.06 vs prototype 0.64).
   The baseline is 0.997 immune and β sits at its bound.
6. **Borderline small clones are unstable across CPU types** (SMC19 on v4 Intel vs
   v6 AMD; floating-point differences of about 1e-16). A clone-stability report is
   planned, not built.
7. **`programme_overlap` is computed per clone vs normal**, so truncal events dominate it.
   It should be per clone difference.
8. **Known model limits.**
   - Events under about 60 genes cannot pass the 3-test-bin minimum.
   - Balanced 2:2 gains escape the allele veto.
   - A programme's arm-mean part is left to dosage on short arms.
   - Potts absorbs randomly scattered normal cells.

## Planned next steps

1. Restore green tests (revert or fix the region NB size). Run the closest tests, then
   the full suite.
2. Fix candidate ranking for region mixtures (scale subclone), then rerun the 100k and
   1M scale benchmarks (`slurm/scale.sbatch`, target ≤ 8 GiB peak). Record wall time
   and memory per 100k cells.
3. Rerun all real data (YCC ×5, SMC, Moorman, Visium HD ×4) on green code. Publish SMC
   and Moorman through `~/CRC/slurm/haplodose_cpu.sbatch`. Extend the CRC adapter
   for allele runs (YCC) and spatial runs (Visium HD, U0247 positions for samples
   without `positions.tsv`).
4. Clone-stability output (refit under perturbed seeds and subsamples; report how
   often each clone recurs).
5. Per-difference `programme_overlap`; move `arms.tsv` from untested to tested.
6. Benchmarks still owed from the plan:
   - comparison with the frozen prototype and infercnvpy on the synthetic grid;
   - patient-mixture benchmark;
   - Moorman arms vs MSK-IMPACT.
7. GPU (M4): submit `slurm/gpu_smoke.sbatch` only after explicit approval
   (`--prefer=b6k` after a smoke test).
8. Rewrite `README.md` (still describes M1 only) and add journal entries for this work.

## Rules this work follows

- Show the sbatch header before submitting.
- Account `innovation`, partition `cpu`, `--prefer=v6`, 90 h wall time.
- No arrays.
- Tests on compute nodes.
- CRC source data read only.
- Outputs to new directories.
- Never edit the frozen prototype (`~/CRC/src/43_cnv_caller/`).
- Each GPU run needs approval.
- RNA clone groups are not validated DNA clones; patient is the biological unit.

## Debug scripts (`logs/dbg/`)

| Script | Purpose |
|---|---|
| `dsrun.sbatch` + `summarise.py` | run a CRC prepared dataset, summary table |
| `yccrun.sbatch` + `numbat_compare.py` | YCC with alleles vs native Numbat |
| `spatrun.sbatch` + `dna_compare.py`, `dna_only.sbatch` | Visium HD run and DNA arm comparison |
| `agreement.py/.sbatch` | malignant and clone agreement with CRC infercnvpy calls |
| `chr19_genes.py` | gene-level view of whole-chromosome splits |
| `smc_window_check.py`, `window_shapes.py` | programme high-pass checks |
| `cnloh_dbg.py`, `spatial_dbg.py`, `ycc_alleles.py`, `prog_test_seeds.py`, `smc_subset.py` | simulation and allele debugging |
