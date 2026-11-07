# HaploDose scientific and codebase audit

Reviewed 2026-11-02. Scope: the current working tree, including uncommitted changes, rather than only HEAD or the frozen benchmark package. Review and documentation only; no implementation, result, protocol, or readiness changes were made.

## Assessment

HaploDose has a coherent computational design and a substantially stronger implementation than the earlier reviews describe. The existing Python suite passed **89 tests, with 1 GPU test skipped**, on the current CPU environment. Pooled profiling, count-based cell assignment, reference-only expression programmes, sparse aggregation, and explicit relative-dosage outputs are useful foundations.

**[Inference] The main obstacle to a rigorous manuscript is validation of the scientific claims, followed by a few input and scoring weaknesses.** Passing implementation tests does not establish independent clone validation, calibrated false-split rates, absolute DNA copy number, or superiority across methods. The most useful next work is to make those claims testable and verify them at the appropriate biological resolution. Additional datasets or faster kernels alone will not resolve these issues.

This audit examined the package, Rust aggregation source, tests, prior reviews and responses, manuscript/status documents, and selected benchmark preparation, scoring, comparator, and scaling code. It did not reanalyse biological cohorts or inspect unopened SPECTRUM/Bian data. Existing benchmark numbers in status documents were not independently reproduced and are not endorsed here as current results.

## Prioritized findings

P1 means a correctness or interpretation issue to resolve before relying on the affected claim. P2 means a material limitation or robustness improvement. A source-confirmed dependency establishes that the claimed guarantee does not follow; it does not quantify the resulting bias. Real-data effects are marked where unmeasured.

### 1. P1: Count splitting does not make the complete clone-validation procedure independent

Evidence: [clones.py](../python/haplodose/clones.py), initialization at lines 393–409, proposal-half merge labels at line 711, programme updates at line 825, EM at lines 838–887, and split proposals at lines 913–941.

The newer thin-mode merge check forms labels from proposal-half counts, an improvement over the earlier implementation. However, those likelihoods still use clone doses, offsets, dispersions, depth effects, programme scores, mixture proportions, and library totals estimated using the full counts. Clone doses are repeatedly rebuilt from full-count responsibilities. Split candidates also start from full-count MAP membership. The same test half is reused for adaptive merging, candidate selection, refinement, and reporting. Default even/odd mode retains the earlier full-data assignment dependency and cannot assume neighbouring genomic bins are independent.

Thus, proposal-half assignment alone does not establish independence from the test half. The opening docstring's claim that a split cannot be supported by the noise that created it is stronger than the implementation supports. A posterior threshold of 0.95 is also not a 5% family-wise false-split rate.

The thinning distribution itself is appropriate for negative-binomial counts under its assumptions. Its independence result requires the relevant NB model and size parameters. The code acknowledges estimated sizes; additionally, sums of heterogeneous gene-level NB variables need not be exactly NB, and biological correlations survive a marginal count model. The existing moment test verifies a correctly specified, one-bin synthetic case, not this complete adaptive procedure. [Primary NB splitting paper](https://arxiv.org/abs/2307.12985), [general data-thinning paper](https://arxiv.org/abs/2301.07276).

Recommended action: describe current cross-fitting as an evidence-separation heuristic. For a formal independence claim, fit all assignment-defining quantities on training data or independent references and reserve validation data for a prespecified decision. In either design, measure specimen-level false splits for the complete search across depth separation, expression programmes, reference mismatch, rare populations, and spatial structure. **[Unverified] The magnitude of false-split inflation from the remaining dependencies has not been measured in this audit.**

### 2. P1: Allele-balanced tumour cells can become a whole-genome normal reference

Evidence: [alleles.py](../python/haplodose/alleles.py), `balanced_cells` at line 453; [pipeline.py](../python/haplodose/pipeline.py), lines 145–158; `CloneParams.in_sample_reference` defaults to true.

The selection establishes balance at the specimen's selected imbalance loci. It cannot establish that a cell is diploid or nonmalignant elsewhere. A tumour clone that lacks another clone's LOH event can pass this selection while sharing truncal CNAs with that clone. Its entire expression profile is then added as an in-sample reference component.

Reproduced selection counterexample: two groups of 50 cells, 80 SNPs, depth 4 per cell/SNP, with allele counts 2/4 versus 4/4. `balanced_cells` selected **all 50 balanced cells and none of the other group**. Designating both groups as neoplastic does not change that calculation; the function has no independent evidence of normality. This was an allele-selection probe, not an end-to-end measurement of lost CNA sensitivity.

**[Inference] In a tumour-only sample, the added component can absorb shared tumour dosage and change normal-like assignments or truncal profiles.** Require independent normal-cell evidence before treating these cells as a normal reference, or explicitly model them as an uncertain query-derived baseline. Test two neoplastic clones with shared truncal CNAs and a private allele event, with no true normal cells. Preserve the existing option as an ablation rather than assuming it is always beneficial.

### 3. P1: Missing RNA arm values can improve the clone-difference AUROC

Evidence: [Andor scorer](../benchmarks/andor/score.py), `max_pair_difference` at line 44; [SPECTRUM scorer](../benchmarks/spectrum/score.py), lines 133–152, which imports that function.

Missing clone-arm contrasts are replaced by zero. This assigns absent measurements the biological meaning “no difference.” Reproduced with four arms: truth `[1, 1, 0, 0]`, observed pair contrasts `[1, 1, NaN, NaN]`. The scorer returns `[1, 1, 0, 0]` and **AUROC 1.0**, although the second clone has measurements on only two arms.

Recommended action: retain an explicit measurement mask and report arm coverage per clone, specimen, and patient. Prespecify a common evaluation domain, minimum coverage, and a failure-aware rule; do not select a different favourable denominator for each tool. A completely observed profile with fewer than two inferred clones is a different case from missing profile measurements. **[Unverified] Whether current biological runs contain consequential missing arms is not established here.**

### 4. P1: The SPECTRUM “SOTA” flag can be true without the principal comparators

Evidence: [paired_test.py](../benchmarks/spectrum/paired_test.py), lines 71–83.

The comparison family contains only tools with at least one valid patient score. `sota` requires superiority only over that family. A synthetic fixture with 12 patients, HaploDose AUROC 0.8, infercnvpy AUROC 0.5, and every other comparator absent produced **`sota: true`**, with Numbat, CopyKAT, and SCEVAN among `not_run`.

This follows the documented availability rule; it is not a discrepancy between that rule and the code. It does mean the flag establishes superiority only over the available subset. A comparator with one successful patient also enters with all its other missing patients imputed to 0.5, mixing discrimination and operational success.

Recommended action: distinguish “superior to the evaluated subset” from a field-wide leadership decision. For future confirmatory evaluation, fix required comparators and acceptable coverage before outcomes are examined, retain failures in the accounting, and report discrimination and reliability separately. Preserve existing frozen protocols and disclose any later amendment or sensitivity analysis.

### 5. P1: Shared references and reused biological lines limit independent superiority claims

Evidence: [held-out analysis](../benchmarks/heldout/analyse.py), lines 40–64; [sensitivity analysis](../benchmarks/heldout/sensitivity.py), `cluster_signflip`, `mixed_model`, and `claim1_decision`; [MIX-seq preparation](../benchmarks/heldout/prepare_mixseq.py), lines 78–98; [CCLE scoring](../benchmarks/ccle_bulk/score.py), `primary`.

The MIX-seq primary Wilcoxon tests and bootstrap resample lines despite shared references. The protocol now discloses this, adds common-gene scoring, and supplies dependence sensitivities. These are meaningful improvements. The primary claim decision nevertheless remains conditional on line independence. Holm correction addresses the comparison family, not correlated observations.

The six pool/rotation reference sets are also not six separately collected experiments: reciprocal rotations reuse the same pool's cells as queries and references. **[Inference] Treating those six sets as independent clusters can still miss dependence between reciprocal rotations.** Report the original pool as an additional grouping level; with only three compositions, broad across-pool inference is necessarily limited. A mixed-model optimizer returning a finite standard error does not supply additional independent replication.

CCLE lines share GTEx panels within lineage groups and share lineage/protocol characteristics. In addition, some lines' DNA truth was already used in development. The exclusion sensitivity is appropriate and should accompany any independence claim; new bulk RNA from a previously used line supplies a new assay, not a new biological line.

Recommended action: keep conditional line-level results labelled, show pool/lineage effects and reference-panel sensitivity, and obtain independent biological evaluation for generalization claims. For a frozen analysis, report the prespecified result and these limitations together rather than replacing an inconvenient primary test.

### 6. P1: Fractional inputs pass validation and can create negative thinned counts

Evidence: [io.py](../python/haplodose/io.py), `CountSource.read` at line 104; [clones.py](../python/haplodose/clones.py), `thin_counts` at lines 106–124 and `Rows._fold` at line 175. Bulk input has no equivalent comprehensive count-value validation.

The count reader checks finiteness and nonnegativity but accepts fractional h5ad/CSR values. Thinning rounds each value before drawing an integer proposal count, then subtracts that draw from the original fractional value. For 100 cells with count 0.6, size 5, seed 1, **44 test counts were negative**, with minimum approximately −0.4. A temporary CSR H5 containing 0.6 was accepted by `CountSource`.

Recommended action: enforce the raw-count contract before binning/thinning, with a clear error for normalized or log-transformed matrices. Validate finite, nonnegative bulk counts and positive sample totals. If estimated fractional counts are to be supported, explicitly define their likelihood and disable procedures that require integer sampling rather than silently rounding them.

### 7. P2: Metadata permits silently combining patients and datasets into one specimen

Evidence: [io.py](../python/haplodose/io.py), `read_cells` at line 134; [pipeline.py](../python/haplodose/pipeline.py), lines 112–134.

`read_cells` validates columns, cell-ID uniqueness, and matrix order. It accepts the same specimen ID under two patients and two datasets. The pipeline pools query rows by specimen and then chooses the dataset and patient from its first row, using that choice for the reference and SNP filter.

Reproduced: two query rows with specimen `S`, patients `P1/P2`, and datasets `D1/D2` pass the reader. **[Inference] Real metadata with locally reused specimen names can therefore fit unrelated cells against the wrong reference or genotype context.** Validate the required one-patient/one-dataset mapping for each fitted specimen, or use an explicit composite key. Also validate role/split values and nonempty identifying fields. Do not require global specimen uniqueness for unrelated reference rows unless the fitting contract needs it.

### 8. P2: Bulk calibration flags and segment precision exceed their statistical support

Evidence: [bulk.py](../python/haplodose/bulk.py), lines 172–197; [profile.py](../python/haplodose/profile.py), `segments_table` at line 137; [pipeline.py](../python/haplodose/pipeline.py), lines 362–363; [clones.py](../python/haplodose/clones.py), `bin_tau` at line 607.

Bulk `arm_calls_calibrated` becomes true solely when the panel contains at least ten normals. The accompanying text correctly says this is panel-size-based, not empirical calibration. Sample size alone cannot establish calibration under reference mismatch, batch effects, or a changed assay. The reported normal call rate also uses normals that define the median/MAD null and is explicitly optimistic. `|z| > 3` is neither a multiple-testing correction nor a genome-wide false-call guarantee.

Separately, the HMM includes per-bin reference variance `v`, but single-cell segment standard errors and minimum detectable effects receive only the scalar group `tau`. Bulk passes a median standard deviation, losing local heterogeneity. For three bins with one million counts each and tau 0.02, `segments_table` reports SE **0.01156**; including independent per-bin variance 1 in the same approximation gives **0.57747**. This demonstrates omitted variance, not an independently calibrated confidence interval. Genomic correlation and segmentation selection create further uncertainty.

Recommended action: use separate fields for panel-size adequacy and measured predictive calibration. Pass the relevant local variance into precision summaries and label them conditional/model-based. Assess held-out biological normals and the declared arm/segment hypothesis family before attaching error-control language.

### 9. P2: Bulk focal-event scoring accepts windows with very little DNA truth

Evidence: [subarm.py](../benchmarks/tcga_bulk/subarm.py), `window_flags` at lines 159–169 and `per_patient` at line 253; imported by the CCLE scorer.

The truth denominator is the number of genes with any DNA value, not all genes in the window. A ten-gene window with one known neutral gene and nine missing values is labelled **DNA-neutral**, reproduced directly with `window_flags`. One observed focal-loss gene can similarly define a focal-loss window. The docstring's description that all genes are neutral is therefore too strong without a coverage qualification.

Recommended action: prespecify minimum DNA coverage and minimum observed genes, preserve missingness, and report excluded windows consistently across methods. Ensure RNA window summaries also meet a common coverage rule.

The fixed-FPR endpoint is useful as a discrimination benchmark: thresholds are estimated from DNA-neutral windows in the same evaluation cohort. It does not demonstrate a deployable threshold giving 5% FPR on new samples. Thresholds, normal panels, and lineage sharing also induce dependencies that a bootstrap of already-computed per-line recalls does not propagate. Distinguish per-direction FPR from `fpr_neutral_*`, which combines both tails.

The CCLE rule combines a signed-rank p-value with a positive mean difference, while its mean-effect CI is descriptive. This can legitimately satisfy the declared rule even if that CI crosses zero, as the status document acknowledges. It does not establish a confidently positive mean improvement. Signed-rank tests concern the difference distribution, with symmetry assumptions, rather than directly testing its arithmetic mean. Match the effect estimand, interval, and decision rule in future protocols. [SciPy's test definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html).

### 10. P2: Genome assembly and interval conventions are not enforced by the public interface

Evidence: [io.py](../python/haplodose/io.py), `read_positions` at line 148; [bins.py](../python/haplodose/bins.py), lines 27–98; [CLI](../python/haplodose/cli.py).

Position tables carry no assembly or coordinate-origin contract. Arm labels and default programme boundaries always use bundled hg38 centromeres, while the package accepts other coordinates without an assembly check. Consecutive-gene bins can span centromeres or large physical gaps, and the HMM transition probability is per bin rather than per base pair.

**[Inference] Non-hg38 inputs can receive wrong arm summaries or programme boundaries, and gene-poor regions can receive misleading physical-resolution interpretations.** This audit did not establish an assembly mismatch in an existing result. Make the supported assembly explicit, record it in outputs, expose matching centromeres, and validate interval conventions. Treat gap-aware transitions/bins as a separately tested method change, not a routine annotation fix.

### 11. P2: Allele-level uncertainty is not carried consistently into cell assignments

Evidence: [alleles.py](../python/haplodose/alleles.py), `normal_snp_filter` at line 140, beta-binomial emissions at line 203, cell likelihood at line 297, and concentration selection at line 353; [clones.py](../python/haplodose/clones.py), `loglik` at line 804.

Pooled states use a beta-binomial model, but cell evidence is a plug-in binomial likelihood with fixed shrinkage and no concentration parameter. Bias screening is binomial, and the reference concentration is estimated after filtering imbalanced sites, then can be increased using tumour-pool likelihood. Multiple SNPs from one transcript and allele reads already included in expression counts also require a stated joint/conditional model; marginal likelihood addition alone does not prove independent evidence.

These choices can be useful approximations, but pooled calibration does not automatically validate cell posterior confidence. **[Unverified] This audit did not measure assignment miscalibration from them.** Test unequal allele depth, cell-specific ASE, correlated SNPs, mapping bias, phase errors, and allele/expression dependence. Derive a compatible per-cell model or calibrate the composite evidence weight. The previously repaired mirrored-imbalance contrast should remain in those checks.

## Model reasoning and manuscript interpretation

The count-based models fit the RNA modality. Their consequential assumptions should be explicit: reference adequacy, compositional normalization, the four relative dosage states, common dosage slope, pooled residual variance, and the extent to which normal-derived programmes transfer to tumour cells. Avoid interpreting their successful execution as evidence that each assumption holds.

Three additional points matter:

- **RNA dosage needs an RNA estimand as well as DNA concordance.** The manuscript outline currently puts DNA truth under relative RNA dosage. DNA concordance validates use of RNA as a DNA proxy; it does not directly validate estimation of relative RNA expression. Use separate RNA recovery/predictive checks and DNA concordance, then clone-membership truth at cell resolution.
- **Balanced WGD is unidentifiable from normalized dosage plus allele fractions alone.** Doubling every haplotype copy count leaves normalized composition and B-allele fractions unchanged. The outline's statement that WGD is therefore claimed through allele features needs qualification. The [bulk allele benchmark](../benchmarks/ucan_bulk/alleles/score_alleles.py) is a supervised predictor of DNA-derived WGD/ploidy proxies, not direct identification of genome doubling. Preserve its out-of-fold design, validate externally, and distinguish classifier performance from mechanistic identifiability. Its patient bootstrap holds fitted predictions fixed and does not include refitting uncertainty.
- **An arm-difference AUROC is not clone recovery.** It can validate localization of between-population heterogeneity without validating the number of clones, their genotypes, ancestry, or cell membership. SPECTRUM's aggregate DNA/RNA correspondence should be described at that resolution. A linkage tree of RNA differences is likewise a similarity tree unless independent evidence supports ancestry.

The HMM prior also deserves an explicit design decision. For the current sticky transition, stationary occupancy is proportional to `q * (1 - q)`, not the supplied initial/jump weight `q`. Recalculation gives neutral occupancy **0.3788**, rather than 0.5, for dosage states; difference-null occupancy **0.4516**, rather than 0.7, for difference states. The calculated distributions satisfy the implemented transition to approximately 1.2e-8. The rows are normalized correctly; this is not a row-sum bug. Decide whether initial occupancy, stationary occupancy, or common state duration is intended, and evaluate changes through calibration rather than silently replacing a valid transition construction.

Gene filtering also constrains the attainable result. `type_mismatch_genes` assumes extreme query/reference expression cannot plausibly be dosage, and the four-state HMM truncates available relative states. **[Inference] These choices can lose extreme focal amplifications or absorb them as mismatch.** Record exclusions, state saturation, and resolution limits; challenge the assumptions using DNA-backed high-amplitude events. Do not tune those filters around the CRC motivation or a favourable benchmark result.

## Computational efficiency

The current hybrid design is appropriate: Rust aggregates sparse counts; NumPy/SciPy organize the fit; JAX compiles likelihood and HMM kernels. There is no evidence here supporting a language rewrite. Optimize measured end-to-end work, including input layout and clone count.

The table gives source-based costs and proposals. No new scaling benchmark or heap/I/O profile was run.

| Path | Cost or risk | Improvement and acceptance check |
|---|---|---|
| `shared_profiles`, matching, and split validation | All tumour pairs generate HMM tracks each update. Forty tumour clones create 780 pairs; proposal `max_groups` is not a global bound on accumulated clones. | Batch exact pair tests under a memory budget, cache valid unchanged results, and stage candidates. Measure preservation of accepted differences and false splits. |
| Thin-mode `Rows` access | Each request redraws deterministic beta-binomial halves using a Python loop per cell. Repeated proposal/test pooling repeats random generation and reads. | Profile thinning separately; compare a bounded fold cache or fused grouped accumulation against current memory/I/O. Preserve row-key determinism and exact complement counts. |
| Raw-gene loading and reporting | A specimen's subsample and clone reporting reread source storage blocks. Interleaved specimens amplify I/O. | Aggregate final gene pools for several specimens in a streaming pass, or cache a bounded used-gene representation. Measure bytes read as well as wall time. |
| Output retention | `post_store` retains all responsibilities; profile rows accumulate as Python dictionaries. Pooled and E-step work create additional arrays. | Stream specimen outputs and consolidate deterministically. **[Estimate]** One million cells × 40 float32 responsibilities alone is 153 MiB, before copies, likelihoods, and tables. |
| Rust grouped aggregation | Each Rayon accumulator has groups × genes float64 entries. Memory depends on group count and concurrency despite sparse input. | Budget accumulator size and group batches. **[Estimate]** Sixteen such buffers for 500 groups × 30,000 genes total about 1.79 GiB, before reduction/output overhead; actual Rayon buffer lifetime needs measurement. |
| Bulk reference/programme fitting | Each leave-one-out normal refits a full SVD; beta-grid ranking computes posterior passes even when only likelihood is required. Bulk tables are loaded densely. | Consider seeded truncated SVD, forward-only likelihood ranking, and bounded sample input. Preserve fold-specific transformations and compare scientific outputs within declared tolerances. |
| Allele loading | Chunked TSV reads are retained in `parts`, concatenated, and converted to sparse matrices; chunking does not bound total ingestion memory. | Use a site-index pass followed by bounded sparse assembly, or a documented on-disk intermediate. Measure peak memory against nonzero allele rows. |
| Iteration stopping | EM stops on hard-label changes below 0.5%, despite moving probabilities and nuisance estimates. Programme scores restart from zero for five steps. | Monitor a compatible objective/residual and parameter changes, with a cap; consider warm starts. First define the objective for the alternating approximations. |

For resource claims, retain the existing matched-hardware, fixed-input, repeated-run design and report failures and preprocessing costs. Scaling inputs above the source cell count contain duplicated real cells: useful for throughput, not evidence about independent biological complexity or rare-clone recovery. Slurm measurements are sampled and should not be presented as exact instantaneous peaks. The local `sacct --noconvert` byte conversion was spot-checked and agreed with `--units=G`; no unit-conversion bug was found.

## What has improved since earlier reviews

Current source and passing regression tests support the following repairs; they should not be repeated as outstanding bugs:

- Zero-UMI/depth-disjoint comparisons return an explicit untestable result rather than NaN difference segments; zero-SNP specimens run without an allele-channel crash.
- Mirrored allelic imbalance is compared using a signed per-SNP contrast, and allele vetoes include ratio-preserving exceptions. Joint-state wording no longer equates relative-neutral imbalance with established cnLOH.
- Shared profiles preserve atomic event boundaries, gene reference shares use fitted mixture weights, and `assignment_dose` distinguishes assignment profiles from diagnostic HMM doses.
- Cell likelihoods retain small contrasts using likelihood ratios; HMM/Poisson-lognormal numerical regression checks pass on this host.
- Duplicate-coordinate spatial graphs exclude self-edges. Median-mode bulk leave-one-out processing now rebuilds the reference-dependent pseudo-count transformation inside each fold.
- Held-out sensitivity code now addresses common genes, native cell presence, and some reference-set dependence. These additions reduce earlier fairness problems without eliminating every inferential limitation.

The README's statement that the latest tests have two failures is stale relative to this audit. Passing here is one-host evidence, not proof that previously documented node-dependent thresholds or GPU results are resolved.

## Recommended sequence toward the manuscript

| Stage | Concrete deliverable | Evidence needed to proceed |
|---|---|---|
| 1. Scientific contracts | Define supported counts, assembly, specimen keys, relative dosage, candidate clones, and calibration terminology. Repair missing-measurement scoring and separate required-comparator decisions. | Adversarial count/metadata/coverage cases cannot silently yield invalid inputs or unsupported success flags. Preserve frozen result provenance. |
| 2. Clone calibration | Define what training/test separation guarantees and calibrate the complete adaptive procedure. Challenge the in-sample baseline. | Biological null specimens, tumour-only balanced/imbalanced mixtures, depth-confounded and programme-confounded negatives; patient/line-level false splits with intervals. |
| 3. Focused method contribution | Compare pooled sharing, programmes, depth correction, thinning, alleles, and spatial support against simple alternatives. | Ablations improve the endpoint they target without unacceptable false splits, lost rare clones, or reduced coverage. Include negative results and detectable-effect curves. |
| 4. Independent evaluation | Freeze code, input identities, endpoints, scoring domains, failure rules, comparator set, and inference unit. | Independent patients/lines; matched DNA where claimed; cell-level DNA/RNA for membership. Report effect sizes, compatible uncertainty, and the fixed correction family. |
| 5. Efficiency and reproducibility | Profile and optimize accepted methods; generate versioned figures/tables and a compact final evidence ledger. | Matched inputs/hardware, repeated wall time and memory, I/O accounting, environment inventories, deterministic provenance, and declared numerical tolerances. |

For the internal CRC/IGF2 question, retain a separate patient-level analysis. Establish IGF2 tandem-duplication status and genome-wide gains from DNA; evaluate clone membership with matched evidence where possible. State the biological n, repeated samples, gain-burden estimand, and prespecified confounders such as purity, ploidy/WGD, subtype, and assay. Decide whether ploidy/WGD is a confounder or a mediator before adjustment, because the choice changes the question. RNA-derived clone groups can support an association analysis; they cannot by themselves establish that the duplication caused genome-wide gains. Keep this application out of the general package's public motivation.

## Verification and provenance

Executed on compute host `ne1dc5-003.nygenome.org`, inside existing Slurm allocation **22345309** (8 CPUs, 90 GB). Tests used CPUs 4–5 with two numerical-library threads; probes used CPU 6 with one thread. No job was submitted, no dependency was installed, and no build or analysis result was overwritten.

| Check | Outcome |
|---|---|
| Existing Python suite | **89 passed, 1 skipped in 351.68 seconds**; GPU equivalence skipped because no GPU was visible. Native Rust-versus-SciPy Python tests ran with the existing extension. |
| Fractional-count ingestion/thinning | CSR reader accepted 0.6; 44/100 negative complementary test counts at seed 1, size 5. |
| Mixed specimen metadata | Reader accepted specimen S spanning patients P1/P2 and datasets D1/D2. |
| Allele-balanced reference selection | Selected 50/50 balanced cells and 0/50 imbalanced cells in the 80-SNP counterexample. |
| Missing-arm scoring | Two observed positive arms plus two missing negative arms produced AUROC 1.0. |
| Available-comparator decision | Synthetic 12-patient fixture returned `sota=true` with only infercnvpy present. |
| DNA-window coverage | One neutral DNA gene out of ten labelled the window neutral. |
| Segment precision and HMM occupancy | Recomputed the variance example and stationary distributions reported above. |
| Working-tree/source stability | Initial dirty state recorded; 20 fingerprinted package/scorer files remained unchanged throughout inspection and tests. |
| Report validation | All 25 local links resolved; code fences balanced; no trailing whitespace; all 8 embedded source fingerprints matched. |
| Git whitespace check | `git diff --check` passed before and after report creation. |

Test command from the repository root:

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=python \
  OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 RAYON_NUM_THREADS=2 \
  JAX_PLATFORMS=cpu taskset -c 4,5 \
  /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
  -m pytest -q -p no:cacheprovider
```

Environment: Python 3.11.18; HaploDose 0.1.0; NumPy 2.4.6; SciPy 1.17.1; JAX/JAXlib 0.10.2; scikit-learn 1.9.1; pytest 9.2.0. The existing native extension was available. No environment changes were made.

Not run: Rust rebuild or standalone Cargo tests, new cohort inference, comparator inference, GPU execution, multi-host reproducibility, complete-search false-split calibration, new large-scale timing/memory/I/O measurements, or re-estimation of biological benchmark effect sizes. Probe results establish the specified constructed cases, not their prevalence in biological data.

HEAD: `54b0cf2a4b94c1dfa42ae4a5e792e88bd0a4b2e2`. The working tree already contained 13 modified tracked files and multiple untracked code/data-development directories and reports. HEAD alone does not identify this reviewed implementation. Selected SHA-256 fingerprints:

```text
2d6101d59dd3bb68eb7d7a09a07163e6f100addcf1d5f597c117b67df627b942  python/haplodose/clones.py
bce796c10b98e9b33e9f456fc4f25c8ec69fb6b9f9316e826c7156c938485316  python/haplodose/pipeline.py
fb48bc9f40549a9aeb06fdb49e44c62851af444eda6c47f435ae236f5b83bc0e  python/haplodose/alleles.py
466944e8d4137592f8ed0f7f95d35a5e7087ad72398a836634ff82269fbc03e9  python/haplodose/bulk.py
cc2c5818203fb734a3c68af80cc0060f53d0ce8f268eeea02b24cbc385705a63  python/haplodose/profile.py
d86465cf8aafc1e32c8f2bf66d631d376873e5b7a7748ec2bf5a12f8a3f5c7d2  python/haplodose/kernels/estep.py
7435c4763683a8d27ad769d025422ae3e861a6c0dbd89d452bc275491a2d0781  benchmarks/spectrum/score.py
ee87286bb412538904f4a4c23d6914a0261452e8ffaa192c904ea3f9f373a2a2  benchmarks/ccle_bulk/score.py
```

Only this Markdown document was authored by this audit. Existing reviews and user changes were preserved. Recommendations above are proposed work, not completed fixes or permission to alter frozen evaluations.
