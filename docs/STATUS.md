# Haplodose: status, decisions and open work

Written 2026-11-01. Covers work from 2026-10-31 to 2026-11-01. Plan file:
`~/.claude-teams/plans/t3-infercnv-like-caller-handoff-dapper-cat.md` (milestones M0-M4).
There is no git repository in `~/Haplodose`; nothing is committed.

## Current state (read this first)

- **The working tree is not green.** The last change (negative-binomial size of
  each region's totals estimated by moments in `SpecimenModel.event_mixtures`,
  `python/haplodose/clones.py`) fails 2 of 40 Python tests (job 22593885):
  `test_subclone_assignment_above_chance` and `test_subclonal_difference_reported`
  (`tests/test_pipeline.py`). The last fully passing state was job 22593863
  (40 passed + 3 Rust), before that change and before the `arms.tsv` output.
  Next step: revert the moment-based `r_nb` (restore
  `r_nb = len(part) / max(median(phi[cols[part]]), 1e-3)`, times the neighbourhood
  size for spatial data) or fix it so the S1 subclone is still found.
- `arms.tsv` output (added last, in `pipeline.summarise`) has not been run yet.
- Scale run 22593886 (100k cells, 11,000 genes) was still running when this was
  written; its output goes to `results/scale_100k_20261101_0828/`.
- Published CRC results (`~/CRC/results/haplodose_benchmark/{SMC_full,Moorman_full}/`)
  are from jobs 22589677 and 22589678 (2026-10-31). They predate ribosomal-gene
  exclusion, alleles, spatial mode, depth-matched tests, vetoes and the
  cell-type mismatch filter. All later real-data runs are debug runs under
  `logs/runs/` and are not published.

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
