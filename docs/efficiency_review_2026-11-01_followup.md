# HaploDose follow-up: algorithms, efficiency, and competitive performance

Reviewed 2026-11-01, America/New_York. Source fingerprints recorded at 2026-11-02T01:57:31Z.

## Assessment

HaploDose has improved substantially since the [first review](efficiency_review_2026-11-01.md). The remaining obstacles include structural losses of biological information, inconsistent nuisance and variance models, incomplete split validation, and expensive candidate searches. These deserve priority over a language rewrite.

**The current evidence does not establish overall superiority.** There are promising strengths, especially strict-threshold focal-loss detection and throughput, but rare/focal clone recovery remains weak. The appropriate objective is to demonstrate improvements on prespecified tasks with independent DNA truth and matched false-positive rates. Winning every task, assay, and information regime is not something the current model or evidence can guarantee.

This is an inspection and diagnostic report. Implementation files and existing results were preserved. Controlled probes establish specific code behavior; their contribution to real-cohort errors remains unmeasured. Recommendations below are proposed changes, not completed fixes.

## What changed since the first review

The working tree now contains stable PLN evaluation, normalized HMM messages, the corrected softmax score derivative, NB rate estimation, depth-variance correction, sparse bulk binning and matched weights, agglomeration caching, randomized programme SVD, fewer unnecessary Viterbi passes, and full leave-one-out bulk reference fits. These are substantive improvements.

The depth correction remains an approximation, and five programme iterations are not a convergence guarantee. Current clone proposals also include per-bin NB region mixtures, cell size factors, sliding windows, and duplicate-split filtering. These introduce new accuracy and efficiency issues examined below.

## Current benchmark evidence

These figures were read directly from saved artifacts during this review. They are existing runs, not cohort reruns performed for this audit.

| Evaluation | HaploDose | Comparator or earlier version | Interpretation |
|---|---:|---:|---|
| Depth-adjusted Kinker v3 implants, 166 cell-line specimens, mean clone ARI | 0.1473 | infercnvpy 0.1651 | Paired Wilcoxon p = 0.0868; Holm over the one tool pair. No demonstrated overall advantage. |
| Same experiment, 500-gene gains in 15% of cells, 13 lines | 0.0076 | infercnvpy 0.2040 | Descriptive subgroup result; substantial rare-gain weakness. |
| Same experiment, 500-gene losses in 40% of cells, 14 lines | 0.7549 | infercnvpy 0.6527 | Descriptive subgroup result; a promising operating regime. |
| Kinker single-line arm correlation, 181 lines, median | 0.8543 | pseudobulk 0.8634; infercnvpy 0.8126 | Beats infercnvpy but trails the simple pooled baseline. |
| 100,000-cell simulation, four specimens | 404.52 s; 3.45 GiB; median ARI 0.19377 | Before sliding windows: 147.94 s; 3.45 GiB; ARI 0.19377 | 2.73× recorded runtime with identical reported accuracy, sensitivity, and normal false-call metrics. |
| TCGA focal-loss recall, patient mean at 1% per-direction FPR, 3,013 evaluable patients | 0.2541 | naive 0.1745 | Existing segmentation advantage at a strict threshold. |
| Same endpoint at 5% / 10% per-direction FPR | 0.5305 / 0.6026 | naive 0.4720 / 0.6298 | The advantage changes with the operating threshold. |

Sources: [Kinker v3 summary](../results/kinker_implant_v3_20261101_2021/score_20261101_204602/summary.json), [settings](../results/kinker_implant_v3_20261101_2021/score_20261101_204602/by_setting.tsv), [single-line comparisons](../results/kinker_compare_20261101_1616/paired_tests.tsv), [window scale run](../results/scale_100k_windows_20261101_1854/scale_summary.json), [earlier scale run](../results/scale_100k_v2code_20261101_1717/scale_summary.json), [TCGA FPR curves](../results/tcga_bulk_20261101_1405/subarm_20261101_1431/score_fprmatched_random_20261101_1442/curve_pooled.tsv).

The v3 implant results are newer than the pending-run description in `STATUS.md`. The timing ratio is an artifact comparison, not an isolated causal benchmark with repeated timings. Region-overlap recall alone is insufficient: the implant scorer accepts any overlapping reported difference, which can reward broad or excessive calls. Use clone ARI alongside event precision, breakpoint accuracy, and null false splits. A single-line pseudobulk is a strong baseline, but not a mathematical ceiling: robust gene weighting or denoising could improve it.

## Highest-priority findings

### 1. Distinct rare-clone proposals are discarded as duplicates

**Confirmed by direct evaluation of the selection rule.** [clones.py](../python/haplodose/clones.py), `event_mixtures`, lines 812–817.

The duplicate test is `max(mean(g == h), mean(g != h)) >= 0.9`. This is dominated by cells belonging to neither small population. Two entirely disjoint groups of 40 cells among 1,000 have 92% label agreement and zero minority-set overlap. The second proposal is discarded. Both groups exceed the default `min_cells=30`.

In general, two disjoint minority groups with fractions f and g are suppressed whenever f + g <= 0.1. This creates an explicit rare-clone blind spot independently of sequencing depth or event strength. It does not explain the 15% implant result by itself.

**Improvement:** compare minority memberships using overlap/Jaccard or a balanced contingency-table measure, allowing label reversal. Preserve disjoint minorities. Verify identical and complemented splits collapse, while distinct small populations survive.

### 2. The final validation of a multi-child split checks only its first pair

**Confirmed with an isolated control-flow probe.** [clones.py](../python/haplodose/clones.py), `try_splits`, lines 831–875.

After agglomeration, only `members[:2]` enter the depth-matched comparison and gene/allele checks. If that pair passes, all remaining children are instantiated. A stub supplying three surviving children recorded exactly one validated pair, `(0, 1)`, and accepted all three children.

The other children have undergone the earlier agglomeration comparisons. They have not undergone the same final depth matching and vetoes as the first pair. This makes acceptance depend on arbitrary group ordering and weakens protection against nuisance-driven splits.

**Improvement:** accept validated binary splits, or validate every added child against its parent remainder under the same calibrated rule. Require invariance to relabeling the proposal groups. This is especially important when increasing proposal diversity.

### 3. Allele comparisons cannot distinguish opposite-haplotype alterations

**Confirmed analytically and by a direct call.** [alleles.py](../python/haplodose/alleles.py), `_theta2`, `differ`, and `veto`, lines 252–265 and 345–394; [clones.py](../python/haplodose/clones.py), `agglomerate` and `merge_unsupported`.

The comparison squares `(A/D - 0.5)`. Thus BAF 0.8 and 0.2 have the same statistic. With 40 SNPs, depth 1,000 per SNP, and kappa 10,000, the actual functions returned:

| Check | Result |
|---|---:|
| Allele difference score between the two groups | 0.0 |
| Allele veto of an expression difference | True |

Opposite-haplotype events with equal total dosage can therefore merge despite strong allelic evidence. A model that pools balanced mixtures first can also lose the initial imbalance needed to propose such clones.

**Improvement:** retain a shared phase/orientation model across clones and compare signed haplotype fractions, integrating phase uncertainty. A shared switch error should affect both clones together, preserving their relative direction. Add mirrored-CNA and opposite-haplotype LOH cases. This is necessary to compete on allele-specific clone reconstruction: [CalicoST explicitly targets mirrored subclonal CNAs](https://www.nature.com/articles/s41592-024-02438-9).

### 4. Baseline-offset updates omit nuisance effects that the cell likelihood includes

**Confirmed in a deterministic model-consistency probe.** [clones.py](../python/haplodose/clones.py), `_update_offset`, lines 666–673; compare `loglik` and `_update_phi`.

The offset update computes expected counts from library totals and clone doses. It omits the current per-cell depth effect and programme scores. The next cell likelihood adds those effects to the newly estimated offset. Consequently, an effect already explained by the model can be reintroduced through the offset.

Probe: three bins, baseline `(0.6, 0.3, 0.1)`, slopes `(-0.8, 0.2, 1.2)` centered against the baseline, 100 cells at depth 100 and 100 at depth 10,000. Set counts exactly to `L * softmax(log(base) + gamma*x)`, with no CNA or residual baseline error. Starting at the known correct parameters, `_update_offset` returned `(-2.10955, 0, 2.30045)`. Applying that offset plus the existing depth term moved a cell's bin probability by up to 0.3054; the maximum pooled compositional error was 0.2334.

This is a deliberately strong stress case, not an estimate of typical patient error. It establishes that the true model parameters are not a fixed point of this update. [Inference] Similar inconsistency can redistribute programme signal into the baseline and dispersion.

**Improvement:** accumulate expectations using the full normalized cell model in bounded chunks, then optimize a constrained/regularized offset conditional on depth, programmes, and clone assignments. Monitor the same objective before and after each update. Separate reference mismatch from biological dosage through explicit constraints.

### 5. Known reference variance is added after already being absorbed into tau

**Confirmed in a controlled variance experiment.** [profile.py](../python/haplodose/profile.py), `fit_pooled`, lines 97–112, and `group_tau`; [clones.py](../python/haplodose/clones.py), `profiles` and `bin_tau`.

`estimate_tau` subtracts Poisson sampling variance but does not subtract the supplied reference variance `v`. Segmentation then uses `tau**2 + v`. When `v` genuinely describes bin-level reference noise, that component is counted twice.

Probe: RNG seed 81; 12,000 bins; reset every 500 bins; expected count 1,000,000/bin; independent normal log offsets with SD 0.3; Poisson counts; supplied `v=0.09`; no additional residual tau. `fit_pooled(..., rounds=0, v=v)` returned tau 0.299675, and the combined SD became 0.424034 instead of approximately 0.3. Further `fit_pooled` rounds do not update tau.

**Improvement:** estimate residual tau conditional on known per-bin variance and sampling noise, with a robust heteroskedastic likelihood or calibrated moment estimator. Do not simply subtract a median variance in every setting: varying variances and correlated offsets require validation. Verify recovery of known variance components and calibration before measuring sensitivity gains.

### 6. Expression-only search cannot recover after everything becomes normal-like

**Confirmed by a direct function call.** [clones.py](../python/haplodose/clones.py), `split_candidates`, lines 700–716.

Expression proposals iterate over `range(1, len(doses))`; clone 0 is searched only for allele-based proposals. With one normal-like clone and no alleles, `split_candidates` returns `[]`. A diluted tumour event missed during initialization has no expression-based route back into the search.

**Improvement:** permit calibrated residual-event proposals within clone 0, with stricter normal-null control. Preserve the neutral option while allowing a detectable minority to emerge. Include mixtures whose pooled profile looks neutral but whose minority is identifiable.

The fixed 20,000-cell proposal cap also limits discovery: a 0.1% population contributes 20 cells in expectation, below `min_cells=30`, however many cells the full specimen contains. [Estimate] Prefer multiscale sketches, event-based enrichment, or adaptive sampling with explicit sensitivity accounting.

### 7. Reported clone differences bypass the final vetoes

**Established by source inspection.** [pipeline.py](../python/haplodose/pipeline.py), `summarise`, lines 340–357.

The final `differences.tsv` and tree distances use all segments returned by `matched_differences`; they do not call `model.supported`. A pair can survive merging because of one supported event while other vetoed events reappear in its reported differences. The tree also uses expression differences alone, so an allele-only separation need not contribute to its distance.

**Improvement:** share a single evidence-evaluation path between acceptance and reporting. Report candidate, supported, and rejected evidence explicitly with rejection reasons. Include validated allele differences in a clearly defined clone-distance model. An average-linkage tree should be labeled as a similarity tree unless a biological phylogeny model is fitted.

### 8. Data-dependent gene exclusion can remove true amplifications and affect unrelated datasets

**Confirmed by source inspection and a small input probe.** [pipeline.py](../python/haplodose/pipeline.py), `type_mismatch_genes`, lines 211–233, and `run`, lines 74–76.

Genes differing by more than fivefold from every reference type are classified as cell-type mismatches. The function returns a union across datasets, and that union defines one shared binning. A gene mismatched in dataset A is removed from dataset B even when B's reference matches it. A two-dataset input confirmed this union behavior.

The assumption that sufficiently large differences must be cell type is also unsafe for high-copy focal amplification. The four-state model and permanent exclusion can remove precisely the signal needed for an amplification caller.

**Improvement:** preserve the full genomic feature map; use dataset-specific reliability weights and robust gene likelihoods. Evaluate exclusions against validated high-amplification events. Separate poor-reference warnings from irreversible feature removal within a run.

## Further mathematical limitations

### Five Fisher steps can still leave programme estimates far from the optimum

[kernels/estep.py](../python/haplodose/kernels/estep.py), `_chunk_scores`, lines 27–53, starts at zero on every call and caps each step's maximum bin log-rate movement at 0.5. Five steps permit at most 2.5 total movement, regardless of how far the optimum lies away.

With `W=(-1,1)`, baseline `(0.5,0.5)`, counts `(1,999)`, library 1,000, phi 0.001, and prior precision 1e-6, the function returned u = 2.48085. The unpenalized equal-dispersion optimum is `0.5*log(999) = 3.45338`, and the tiny prior barely changes it. The existing convergence test covers a milder case.

Use warm starts, a score or objective convergence check, and a safeguarded step. The corrected derivative is valuable, but the `R @ log_rate` surrogate still differs from the responsibility-weighted clone likelihood when assignments are uncertain. Any speed optimization should preserve or explicitly bound that approximation.

### A single dosage slope does not represent purity mixing exactly

[profile.py](../python/haplodose/profile.py), `CN`, `LOG_DOSE`, and `segment`; [bulk.py](../python/haplodose/bulk.py), `BULK_PROFILE` and `start_beta`.

The model uses `beta * log(C/2)`. Even in a favorable linear RNA-dosage model with matched normal and tumour baselines, mixing tumour fraction p with diploid normal cells gives `log((1-p) + p*C/2)`. At p = 0.2, the effective beta is 0.1520 for C=1, 0.2351 for C=3, and 0.2630 for C=4. No single beta matches all three. RNA contribution can also differ from cell fraction.

Current expression states are C=1,2,3,4, with no explicit homozygous deletion or higher-copy state. Global centering further limits absolute interpretation. Spatial posteriors represent uncertainty over clone labels; they do not explicitly represent within-spot mixtures of normal and tumour cells.

**Improvement:** retain a robust continuous relative-dosage output, and evaluate a separate purity-aware allele-specific model where the data identify it. Treat purity and dosage compensation as distinct quantities. Report ambiguity between ploidy/purity solutions. A uniform balanced doubling is not identifiable from library-normalized expression alone.

### Current cnLOH labels exceed what the state model identifies

[alleles.py](../python/haplodose/alleles.py), `joint_state`, lines 290–302, returns `cnLOH` whenever relative expression state 2 is confidently allelically imbalanced. It does not establish loss of one parental allele or absolute total copy number two. A weak gain, an allele-specific amplification, expression compensation, or an alternative ploidy baseline can produce relative-neutral expression with imbalance.

Report “relative-neutral dosage with allelic imbalance” as the directly supported observation. Reserve a stronger cnLOH interpretation for an identified joint allele-copy/purity model or orthogonal evidence. This matters when comparing with DNA-defined cnLOH truth.

### The phase-free imbalance moment is biased at low depth

For beta-binomial A with depth d, fraction b, concentration kappa, define `h=(1+(d-1)/(kappa+1))/d` and theta=b-0.5. Then

`E[(A/d-0.5)^2] = theta^2 + (0.25-theta^2)*h`.

The subtraction in `_theta2` yields expectation `theta^2*(1-h)`, not theta². At d=2, kappa=50, b=0.8, exact enumeration of the beta-binomial gives 0.04412 rather than 0.09. Unequal depths can therefore confound comparisons. Its variance formula also uses a Gaussian squared-error approximation.

Prefer a direct beta-binomial comparison or a corrected moment with an appropriate variance. Include unequal-depth nulls. Correction cannot restore haplotype orientation after squaring; that requires the separate change in finding 3.

### Segment uncertainty omits heterogeneity, correlation, and selection

[profile.py](../python/haplodose/profile.py), `segments_table`, lines 137–155, reports `sqrt(tau²/n_bins + 1/U)`. Single-cell reporting passes group tau without the per-bin reference variance used by the HMM. Bulk reporting substitutes a genome-wide median SD. Neither propagates nuisance-fit uncertainty, correlated genomic residuals, or adaptive segment selection.

For an illustrative equicorrelated residual with n=20 and correlation 0.3, the residual component's true SE is 2.59× its independence-based SE. This is a mathematical scenario, not a measured HaploDose correlation. Also, mean marginal state posterior across a run is not the probability that the entire run is correct.

Calibrate segment-level uncertainty and minimum-detectable effects using held-out specimens or genomic block resampling. State whether uncertainty conditions on selected boundaries and fitted clones.

## Computational efficiency

### Region proposals: avoid fitting every window from scratch

The new per-bin mixture addresses real dispersion heterogeneity, but each overlapping window runs up to 200 mixture steps, each with two `nb_rate` calls of up to 100 fixed-point iterations. A small CPU profile with 240 cells, 240 bins, three chromosomes, and one implanted gain took 0.299 s under cProfile: 879 `nb_rate` calls used 0.167 s; 969 `logp` calls used 0.092 s. This is a local diagnostic, not a cohort speed estimate.

An efficient RNA-count proposal stage can use NB scores and information:

`score_cj = (Y_cj - m_cj)/(1 + phi_j*m_cj)`

`information_cj = m_cj/(1 + phi_j*m_cj)`.

After handling fitted size factors and composition, chromosome prefix sums give inexpensive multiscale region statistics. Screen for coherent between-cell heterogeneity, then fit exact per-bin mixtures only on promising regions. Prefix sums do not replace the exact mixture fit; they provide a candidate screen whose sensitivity must be measured.

Warm-start neighboring windows, cache denominators where valid, and stop against objective improvement relative to the null gain or a parameter/score tolerance. The present relative stopping rule scales by the absolute likelihood after dropping arbitrary constants and also stops on a likelihood decrease. Record convergence and reject or investigate failed updates. Benchmark safeguarded Newton/root solving for the scalar NB rate before moving this kernel to Rust.

### Remaining repeated work and memory costs

| Location | Remaining cost | Proposed improvement |
|---|---|---|
| `try_splits`, clones.py:825–829 | Reloads the same proposal cells for each candidate; reconstructs full-specimen depth moments inside the candidate loop | Cache per-parent counts/moments or group sufficient statistics for the candidate batch. |
| `pooled_depth_slopes`, clones.py:239–243 | Dense cells × (5 × clones) one-hot matrix and float64 conversions | Direct grouped aggregation or sparse one-hot weights. At 1M cells and 40 clones, one float32 matrix is 0.75 GiB; one float64 copy adds 1.49 GiB, excluding other allocations. [Estimate] |
| `matched_differences`, clones.py:439–451 | Sparse weights improve memory, but every cell still appears in up to K-1 pairs; lists and concatenated COO arrays coexist | Batch pairs, omit zero weights, reuse group membership/depth order. At 1M cells and K=40 there can be 39M stored entries before eliminating zeros. [Estimate] |
| `merge_unsupported`, clones.py:501–507 | Recalculates all pairs after every merge despite the separate cache in `agglomerate` | Cache only when inputs are unchanged; invalidate pairs whose MAP memberships or nuisance inputs change. Soft responsibilities mean this requires more care than the initial agglomeration cache. |
| `leiden_groups`, clones.py:200–205 | Repeated PCA and exact nearest-neighbor searches on up to 20,000 proposal cells | Reuse embeddings/graphs where valid; evaluate approximate candidate graphs while keeping exact evidence checks. |
| `start_beta`, bulk.py:59–65 | Unused initial tau calculation; nine beta candidates, each with two full posterior fits and another forward/backward evaluation | Reuse initial summaries, batch the grid, and add a forward-only likelihood function for final ranking. |
| `reference_model`, bulk.py:50–52 | Full SVD repeated for leave-one-out normal fits | Cache only quantities whose exclusion can be exact; evaluate truncated SVD or explicit cross-fitting without leaking held-out samples. |
| E-step and HMM boundaries | Repeated host/device conversion and compiled shapes varying with clone/group count remain | Keep reusable parameters on device; bucket shapes and mask padding; profile compile time separately from warmed execution. |
| `read_alleles` | Compact tables still accumulate in memory before concatenation and sparse assembly | Partition by specimen or use two-pass sparse assembly for very large allele files. Compact input is an improvement, not bounded streaming memory. |

The memory estimates refer to one specimen at the stated size. They are not measured peak RSS. The power-of-two comment in `_chunks` is only partly true: when the element cap binds, it can produce an arbitrary row count, and K remains a changing compiled dimension.

## What would support a competitive claim

An independent [2025 scRNA-seq benchmark](https://www.nature.com/articles/s41467-025-62359-9) found no single caller best across all tasks. That is a reason to define a rigorous target, not to lower the ambition. HaploDose should aim to improve the accuracy/runtime frontier across explicit regimes.

1. **Freeze development and evaluation cohorts.** Use bulk patients and single-cell biological specimens as the units. Keep patients, related samples, and repeated perturbations together in folds. Repeated implant seeds within a line are not independent biological replicates. Account for shared pools/references in Kinker through pool-aware sensitivity analyses. Record modality, reference design, depth, purity, batch, truth quality, and event prevalence.
2. **Use appropriate competitors for each assay.** For scRNA counts, include CopyKAT, SCEVAN, inferCNV, infercnvpy, and a robust pseudobulk baseline; with alleles, include Numbat and XClone. For spatial mixtures, include CalicoST and an appropriate spatial baseline. For bulk, include the current robust gene/arm baseline and bulk-compatible RNA-CNA methods such as CaSpER after verifying input requirements. DNA callers provide truth; they are not interchangeable RNA-input competitors. The comparator-build script is preparation, not evidence of completed CopyKAT/SCEVAN comparisons.
3. **Measure the complete task.** Report continuous dosage error/correlation, event precision and recall, breakpoint error, gain/loss separately, allele-specific events, malignant-call calibration, rare-clone recall, ARI, false splits on null specimens, failures, runtime, and peak memory. Match inputs and reference access. Include phasing/pileup/preprocessing costs separately and end-to-end for allele-aware workflows.
4. **Fix operating points before testing.** Calibrate thresholds and clone-search rules on development/normal data. Report precision-recall and partial ROC curves on untouched test specimens. The current TCGA scorer sets thresholds using DNA-neutral windows from the scored cancer type: useful retrospective discrimination evidence, but not an independently calibrated deployment threshold. Its nominal 5% FPR is per direction; the union of upper and lower tails can approach 10%.
5. **Calibrate the whole search.** Even/odd bins are reused through all-bin EM, candidate selection, and repeated refinement. Use null runs preserving depth, gene covariance, programme structure, and spatial dependence. Calibration must include every proposal, merge, and selection step. Name the family of hypotheses: for example, added clones per specimen and event calls across the evaluated genome. A 0.95 posterior is not automatically 5% FDR for either family.
6. **Use paired biological-unit inference.** Estimate paired patient/specimen effects with confidence intervals; use an appropriate clustered bootstrap when observations share patients or reference pools. Prespecify primary endpoints and apply Holm across the planned method/endpoint family, or a declared FDR procedure for exploratory screens. Report missing outputs and failure rates rather than silently dropping difficult specimens.

Agreement with Numbat is useful for diagnosis but is not independent DNA validation. [Numbat's primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC10289836/) also makes clear why allele information and pooled clone composition matter. Compare allele-aware tools under matched information access.

## Recommended development order

| Priority | Work | Acceptance evidence |
|---|---|---|
| 1 | Repair rare-proposal deduplication, multi-child validation, normal-clone recovery, and reporting consistency | Controlled cases pass; relabeling proposals does not change accepted evidence; full-null false splits do not increase. |
| 2 | Make offsets, depth/programme effects, dispersion, and reference variance internally consistent | Known-parameter simulations recover variance and nuisance effects; conditional objectives improve; held-out calibration improves. |
| 3 | Add shared-phase clone comparisons and honest joint-state interpretation | Mirrored alterations and allele-only clones remain separable; uncertain ploidy/purity cases are labeled appropriately. |
| 4 | Screen multiscale proposals cheaply, then refine exact mixtures | Candidate recall at least matches current search on development data, with measured lower runtime and bounded memory. |
| 5 | Combine robust broad dosage with local segmentation and adaptive information-based bins | Preserve arm-level baseline performance while improving focal precision/recall at fixed FPR on untouched specimens. |
| 6 | Evaluate a purity-aware, richer allele-copy model | Improved low-purity/LOH/high-amplification accuracy without unsupported absolute-copy claims. |

Adaptive bins should respect chromosome arms, centromeres, gaps, expression information, and SNP coverage. Current ten-gene bins plus at least three odd test bins impose a substantial floor on focal-event resolution. A fixed transition probability per gene bin also imposes different physical-length priors in gene-rich and gene-poor regions. Evaluate physical-distance or information-aware transitions and local gene-resolution refinement rather than simply reducing every bin size.

Keep Python for the API, scientific workflow, and most orchestration. NumPy/SciPy/JAX/Rust already perform substantial compiled work. After removing redundant work, benchmark a compiled region-mixture kernel and small-state CPU HMM if profiles still justify it. GPU benefits are unmeasured here. A language rewrite would preserve the scientific defects unless the model and search also change.

## Verification and reproducibility

Existing compute allocation: job 22345309, host `ne1dc5-003.nygenome.org`, CPU partition, eight allocated CPUs and 90G. Diagnostic processes used bounded CPU affinity and one or two library threads. No jobs were submitted, no packages installed, and no production source edited.

Executed test group A:

```bash
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  MKL_NUM_THREADS=2 RAYON_NUM_THREADS=2 JAX_PLATFORMS=cpu taskset -c 4,5 \
  /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
  -m pytest -q -p no:cacheprovider \
  tests/test_numerics.py tests/test_hmm.py tests/test_profile.py \
  tests/test_alleles.py tests/test_programmes.py tests/test_bulk.py
```

Result: **32 passed in 108.72 s**. A second, nonoverlapping group covered `test_pipeline.py`, `test_spatial.py`, `test_native.py`, `test_bins.py`, and `test_likelihood.py` with the same command prefix: **27 passed in 107.71 s**. Together these covered all 11 current Python test files, with **59 passing tests**. Passing these tests does not establish calibration or exclude the independently demonstrated failures above.

Additional probes used the same environment with `taskset -c 6` and one thread. They evaluated the exact deduplication predicate, the actual offset and allele functions, the all-normal proposal branch, global gene-mask behavior, programme scoring, variance estimation, purity-mixture arithmetic, and an isolated multi-child validation stub. Numerical constructions and parameters are specified in the findings. The candidate profiler used RNG seed 81 after the variance probe, 240 cells at library size 4,000, 240 uniform bins, phi 0.1, three 80-bin chromosomes, an 80-cell 1.5× event in bins 80–111, default windows, and `min_cells=15`. Its timings include profiling overhead.

Installed package versions observed in this audit:

| Package | Version |
|---|---|
| HaploDose | 0.1.0 |
| NumPy | 2.4.6 |
| SciPy | 1.17.1 |
| JAX / jaxlib | 0.10.2 / 0.10.2 |
| scikit-learn | 1.9.1 |

`STATUS.md` mentions scikit-learn 1.7.2; the interpreter used for this audit reports 1.9.1. Record the environment of each benchmark separately.

HEAD was `27dbbe9613d3a10560ca7800802b063cbdeba55e`, with substantial pre-existing modified/untracked source, tests, and benchmarks. HEAD alone does not identify the reviewed code. SHA-256 fingerprints:

```text
dd8d09dcae1179e96537706b74acb7b63e13b525468a8fcd7894d4a9ab7b8e38  python/haplodose/clones.py
c681e451d608427f587a7fcd552946f875a54738c7ff73b592b7e930a0c16529  python/haplodose/alleles.py
ab13f84d57c72bd2de04e13c108e85ab4b5ffd0e82b89ecdadcec094b508fe8c  python/haplodose/profile.py
726198a0cadf32ab2476ca0d9eb53d83d8f3fb8cf93f9c519672583513f67f45  python/haplodose/kernels/estep.py
8d69022370335060d4eeff3a7f938a1d8904b8c72506aa1f023e7d6bd650e33f  python/haplodose/pipeline.py
6e4fcb34d7c91c8f6c284176200ae41044928608364476ff43445363f54cc223  python/haplodose/bulk.py
```

Report verification: all 22 local links resolve; Markdown fences are balanced; no trailing whitespace. The six source fingerprints above were unchanged on the final verification read. The only file added by this audit is this report.

Not run: real-cohort refits, independent rescoring from raw DNA/RNA, full clone-search null calibration, proposed algorithm ablations, GPU benchmarks, or a Rust rebuild. Existing native-extension tests exercise the installed extension. No universal-superiority claim follows from this audit.
