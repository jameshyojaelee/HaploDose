# HaploDose: scientific rigor, correctness, efficiency, and manuscript review

Reviewed 2026-11-03 against the current working tree, including uncommitted changes. This is a review, not an implementation or a change to benchmark decisions. Package code, frozen runs, protocols, and existing reports were preserved.

**Assessment.** HaploDose has a useful computational foundation: pooled genomic profiles, count-based cell assignment, reference-derived expression programmes, sparse aggregation, bounded likelihood batches, and explicit relative-dosage outputs. **[Inference] The present evidence supports a focused RNA-dosage methods paper more strongly than a claim of calibrated cell-level clone identification.** The highest priorities are a reproduced EM numerical defect, calibration of the complete clone search, and benchmark coverage and inference rules. More datasets alone will not resolve those problems.

The existing Python suite passed **95 tests, with one GPU test skipped**. Separate synthetic probes reproduced defects not covered by that suite, including nonfinite cell assignments after removing a small clone.

This review extends the [previous audit](scientific_codebase_audit_2026-11-02.md) to the newer covariance split null and the current manuscript evidence. Findings below distinguish executed counterexamples, source-confirmed behavior, and effects on biological datasets that remain unmeasured. P1 means resolve before relying on the affected output or claim; P2 means a material limitation or improvement.

## 1. Prioritized findings

### F1. P1: Removing a small clone can produce NaN cell assignments

Evidence: [clones.py](../python/haplodose/clones.py), lines 915–934, especially the division at line 934; [estep.py](../python/haplodose/kernels/estep.py), lines 36–61.

After dropping a clone below `min_cells`, EM normalizes the responsibilities of the surviving clones. A cell confidently assigned to the removed clone can have exactly zero floating-point probability in every survivor. The denominator is then zero. The resulting NaN responsibilities enter programme-score estimation and contaminate subsequent likelihoods and posteriors.

**Reproduced using the actual `SpecimenModel.em` implementation.** A 40-cell, 40-bin synthetic input with ten high-depth cells assigned to a component below the 30-cell minimum produced ten rows with nonfinite programme scores, likelihoods, and final responsibilities. The log recorded `drop_small` followed by `n_eff: nan`. No model methods were mocked. This checks an EM starting state, not the frequency of this failure in full biological runs.

Recommended correction: after changing the component set, recompute responsibilities from surviving log-likelihoods and mixture weights in log space. Define the fallback for a cell with no finite likelihood, and check finiteness before accepting the update. Add a regression for removing a high-confidence rare component with programmes enabled. **[Unverified] The prevalence of this failure in existing cohort outputs is unknown.**

### F2. P1: The split null does not establish a 5% specimen-level false-split rate

Evidence: [clones.py](../python/haplodose/clones.py), lines 447–484, 762–787, 898–970, 973–1007, 1128–1278, and 1381–1404; [CLI](../python/haplodose/cli.py), lines 48–55.

The new null simulates a clone's proposal cells and repeats its proposal-and-score search. This addresses selection among the expression proposals in that search. It does not replay the complete procedure that selected those cells, estimated their nuisance quantities, created the initial groups, and repeatedly split and merged clones.

Several qualifications are consequential:

- Initial groups are not null-gated unless `split_null_initial=True`. The frozen claim-4 configuration uses the covariance null without this flag.
- The procedure tests multiple clones and cycles without an explicit specimen-level correction. Beating 19 null draws is a single-search Monte Carlo rule, not genome-wide or specimen-wide error control.
- Null parameters are estimated from selected observed cells and reused without a full refit. This is an approximate fitted-null simulation, not an exact exchangeability argument.
- Thin-mode merge assignments use proposal counts but retain doses, mixture weights, offsets, programme scores, depth quantities, and dispersions informed by full counts. Final reported differences use full-data MAP assignments. The same validation data are reused adaptively.

The thinning distribution itself is appropriate for negative-binomial counts with the required size parameter. Independence depends on that parameter and the count model; estimated sizes and correlated biological effects require assessment. A successful marginal moment test does not validate the complete adaptive analysis. [Neufeld et al., negative-binomial count splitting](https://arxiv.org/html/2307.12985v1).

The [current gate results](STATUS.md), lines 912–950, already show the tradeoff: the covariance variant has clone-count agreement across five seeds in 0.519 of 181 lines; pooling the initial groups raises agreement to 0.917 but lowers implant mean ARI from 0.150 to 0.074. These are documented results, not rerun here. Stable under-splitting is not accurate recovery. Extra RNA groups in MKN-45 also do not conclusively prove false splits: the documentation itself limits what its DNA data can exclude.

Recommended next step: define the error event as an extra tumour clone in an independently justified one-clone specimen, then evaluate the **entire frozen workflow**, including initialization, nuisance fitting, selection, and stopping. Use independent specimens across depth, reference mismatch, cell-state structure, and spatial structure. Report false-split frequency and binomial uncertainty alongside power. Until then, describe the null as a model-based split filter and avoid promising calibrated false-split control or independent validation.

### F3. P1: With alleles, the observed split score and simulated null score differ

Evidence: [clones.py](../python/haplodose/clones.py), lines 988–992, 1180–1190, 1248–1251, and 1275–1277.

Observed candidates can come from allele-defined splits, and their scores add positive allele-difference evidence. Simulated candidates are expression-only: `_null_model` removes alleles and the null invokes `_clone_candidates`, not the allele proposal path. Consequently, an observed expression-plus-allele score is compared with an expression-only null distribution.

The comment that omitting gene and allele vetoes makes null scores conservatively higher does not cover the omitted **positive** allele score or omitted allele proposals. The direction and size of miscalibration require measurement, but there is no established stochastic ordering between these different statistics.

Recommended correction: simulate and score both channels under a compatible joint null, or restrict this gate to expression evidence and assess allele splits separately under a declared rule. Reject unsupported option combinations or clearly label their scope. The documented SPECTRUM primary is expression-only, so this finding does **not** invalidate that primary through an allele dependency; it affects the public allele-plus-null combination and descriptive allele analyses.

### F4. P1: Allele balance is insufficient evidence for a whole-genome normal reference

Evidence: [alleles.py](../python/haplodose/alleles.py), lines 453–471; [pipeline.py](../python/haplodose/pipeline.py), lines 145–158; [clones.py](../python/haplodose/clones.py), line 97.

`balanced_cells` selects cells balanced at the specimen's strong imbalance loci. The pipeline adds their entire expression profile as an `in_sample_normal` reference, enabled by default when alleles are supplied. A neoplastic clone lacking another clone's private LOH can satisfy this selection while retaining shared truncal CNAs. Allele balance also cannot distinguish 1:1 from 2:2 copy states.

**[Inference] This baseline can absorb shared tumour dosage or alter normal-like assignments in tumour-only samples.** The earlier audit's selection counterexample remains consistent with the current code; this review did not repeat the end-to-end sensitivity experiment needed to measure the effect.

Require independent normal-cell evidence, or explicitly treat this as an uncertain query-derived reference and expose its influence. Compare profiles and assignments with and without it. The CLI currently has no switch for `in_sample_reference`, although the Python parameters do; an explicit user-facing option would make this assumption reviewable.

### F5. P1: Missing predictions can improve the clone-difference benchmark

Evidence: [Andor scorer](../benchmarks/andor/score.py), lines 44–56; [SPECTRUM scorer](../benchmarks/spectrum/score.py), lines 129–155; [R comparator wrapper](../benchmarks/comparators/run_r_field.py), lines 108–125.

Missing clone-arm contrasts are replaced with zero, meaning no inferred difference. **Reproduced:** for truth `[1, 1, 0, 0]`, complete contrasts `[1, 1, 2, 2]` score AUROC 0.0. Making the last two contrasts missing changes the returned score to `[1, 1, 0, 0]` and AUROC to 1.0. Missing measurements therefore can improve the ranking.

The comparator wrapper also replaces unmapped/NaN gene values with zero. With `--rule all`, absent cells receive zero profiles and are analyzed. Native coverage and successful biological prediction must remain distinct.

Recommended correction for future scoring: retain measurement masks, require a prespecified common evaluation domain and coverage threshold, and report native cell/gene/arm coverage. A successfully measured one-clone sample and a missing profile require distinct statuses. Preserve frozen primary analyses; disclose any coverage sensitivity as such. **[Unverified] The effect on existing biological rankings was not measured.**

### F6. P1: SPECTRUM can declare success with an incomplete required comparator family

Evidence: [paired_test.py](../benchmarks/spectrum/paired_test.py), lines 64–85; [current run rules](STATUS.md), lines 966–973.

The scorer includes only comparators with at least one valid patient score. Its `sota` flag then requires superiority only over that available subset. **Reproduced:** a synthetic 21-patient fixture with HaploDose AUROC 0.9 and infercnvpy 0.5 returns `sota=true` while Numbat, CopyKAT, and SCEVAN have no output. This availability rule is documented in the script, but the newer run rules specify four comparisons and Holm over four.

A comparator with one valid patient enters the family while all its remaining failures become AUROC 0.5. This mixes discrimination with operational reliability. It can be a declared composite endpoint, but should not be presented as accuracy among successfully analyzed samples.

Recommended action: reconcile the implementation and dated protocol before interpreting the decision; identify required versus optional comparators independently of their outcomes. Report completion and discrimination separately, retain failures, and keep an incomplete required comparison from supporting an unrestricted superiority label. Do not silently change a frozen endpoint.

### F7. P1: Fractional inputs can generate negative thinned counts

Evidence: [io.py](../python/haplodose/io.py), lines 104–112; [clones.py](../python/haplodose/clones.py), lines 166–184 and 235–245; [bulk.py](../python/haplodose/bulk.py), lines 116–137.

The sparse reader requires finite nonnegative values but accepts fractional normalized inputs. Thinning rounds those values to integers before sampling and subtracts the sampled integer from the original fractional value. **Reproduced:** the reader accepts 0.6; with 100 such counts, size 5 and seed 1, 44 complementary counts are negative, with minimum approximately −0.4.

Validate the raw-count contract before binning and thinning. If fractional estimated counts are intentionally supported, specify their statistical treatment and disable incompatible integer-sampling procedures. Bulk input also needs finite/nonnegative value checks, positive library totals, and sample-ID validation. Fractional allele AD/DP values are currently cast to integers rather than explicitly validated.

### F8. P1: A reused specimen name can silently combine patients and reference contexts

Evidence: [io.py](../python/haplodose/io.py), lines 134–145; [pipeline.py](../python/haplodose/pipeline.py), lines 98–140.

The reader checks cell IDs and matrix order, but the pipeline groups query cells by `specimen` alone and chooses patient and dataset from the first row. **Reproduced:** a table containing one specimen name across two patients and two datasets passes validation.

Require one patient and one dataset per fitted specimen, or define and propagate an explicit composite specimen key. Validate nonempty identifiers and role/ref-split vocabulary before fitting. This is a scientific identity issue: the current behavior can pool unrelated patients and use the wrong reference or SNP filter. The counterexample establishes acceptance and the downstream code dependency; no real cohort was found to contain this error during this review.

### F9. P2: Bulk focal truth accepts windows with almost no observed DNA

Evidence: [subarm.py](../benchmarks/tcga_bulk/subarm.py), lines 159–169 and 255–305, reused by the [CCLE](../benchmarks/ccle_bulk/score.py) and [CPTAC](../benchmarks/cptac_bulk/score.py) scorers.

Window truth uses the number of genes with a DNA value as its denominator, with only `n > 0` required. **Reproduced:** one neutral gene among ten genes labels the window neutral despite nine missing values; one observed focal-loss gene similarly labels it a focal-loss window. This is weaker than the docstring's description of all genes being neutral or half the window's genes carrying an event.

Require a fixed minimum number and fraction of genes with DNA truth, and analogous RNA coverage for every method. Record retained and excluded windows. Maintain common targets across methods; the Kinker/MIX-seq primary scorer also derives arm truth and eligibility from each tool's used genes, so its existing common-gene sensitivity remains important.

The 10-Mb endpoint measures window-level focal-loss discrimination. It does not establish breakpoint accuracy, recovery of each focal event, or sensitivity to very small tandem duplications.

### F10. P2: Uncertainty and decision rules need to match the scientific estimand

Evidence: [held-out analysis](../benchmarks/heldout/analyse.py), lines 43–63 and 133–165; [CCLE primary](../benchmarks/ccle_bulk/score.py), lines 51–68; [bulk window scoring](../benchmarks/tcga_bulk/subarm.py), lines 265–305.

Using lines or patients instead of cells as the inferential unit is appropriate. Remaining limits are:

- MIX-seq lines share reference pools. Reciprocal rotations reuse cells, so six pool/rotation labels are not automatically six independent biological experiments. The protocol discloses dependence and adds sensitivities; those do not create new replication.
- CCLE lines and CPTAC cases share reference panels within groups. The fixed-FPR threshold is estimated from DNA-neutral windows in the same evaluation group. Bootstrapping already-computed recall values conditions on the fitted reference and threshold rather than propagating their estimation.
- Signed-rank significance and an arithmetic-mean confidence interval answer different questions. A signed-rank rejection plus positive observed mean can meet the stated rule while uncertainty about the mean spans zero. [SciPy's Wilcoxon definition](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html).

I independently recomputed sample counts, mean effects, signed-rank p-values, and Holm adjustment from the **already reported per-case recalls**. These agree with the saved primary JSON files:

| Dataset and biological unit | Comparison | n with focal loss | Mean recall difference | Holm p, family of 2 | Saved 95% bootstrap CI of mean |
|---|---|---:|---:|---:|---|
| CCLE, cell line | HaploDose − naive | 822 | +0.02431 | 0.003158 | −0.00117 to +0.04933 |
| CCLE, cell line | HaploDose − infercnvpy | 822 | +0.08516 | 1.643 × 10⁻¹⁴ | +0.06227 to +0.10760 |
| CPTAC-3, case | HaploDose − naive | 618 | +0.02839 | 0.02522 | −0.00287 to +0.06030 |
| CPTAC-3, case | HaploDose − infercnvpy | 618 | +0.13559 | 4.294 × 10⁻²¹ | +0.10903 to +0.16325 |

Intervals above were read from saved primary outputs, not independently bootstrapped here. The RNA/DNA-to-window scoring was not rerun. The CPTAC status table reports approximately 2 × 10⁻²¹ under a Holm heading; that is the raw p-value, whereas the saved corrected value is approximately 4.29 × 10⁻²¹. This transcription difference does not change the decision.

Retain the prespecified bulk decisions, but describe the advantage over the naive baseline as small with mean-effect intervals spanning zero. For future confirmatory work, align the effect, interval, and test; report cohort/lineage/protocol heterogeneity; and propagate the shared fitting steps when estimating unconditional generalization uncertainty. A same-cohort 5% per-direction FPR comparison is a discrimination benchmark, not a deployable 5% threshold for new patients. `fpr_neutral_*` combines both tails and is a different quantity.

### F11. P2: Calibration and precision labels exceed the uncertainty represented

Evidence: [bulk.py](../python/haplodose/bulk.py), lines 172–197; [profile.py](../python/haplodose/profile.py), lines 137–155; [pipeline.py](../python/haplodose/pipeline.py), lines 362–363; [clones.py](../python/haplodose/clones.py), lines 667–669; [alleles.py](../python/haplodose/alleles.py), lines 203–208 and 297–301.

`arm_calls_calibrated` is true when there are at least ten normals. Its accompanying text correctly admits that the criterion is panel size, not calibration on these data. Replace this semantic mismatch with separate fields for panel adequacy and observed predictive calibration. Normal-panel self-assessment and `|z| > 3` do not establish family-wise error control across arms or samples.

The HMM includes per-bin reference variance, but single-cell segment precision receives only scalar `tau`. Bulk substitutes a median SD across bins. **Reproduced:** three bins with one million counts each and `tau=0.02` give reported SE 0.01156; adding independent per-bin variance 1 to the same approximation gives 0.57747. This is an omitted-variance demonstration, not a calibrated alternative confidence interval. Genomic dependence and selection of the segment add further uncertainty.

The allele HMM uses beta-binomial emissions, while cell assignment uses a shrunken binomial log-likelihood. Fixed shrinkage does not propagate cell-level overdispersion or clone-fraction uncertainty. **[Inference] High allele depth can give overly sharp assignment evidence under allele-specific expression.** Validate posterior reliability by allele depth and assay; do not transfer clone-level calibration to cell-level probabilities without checking it.

### F12. P2: Assembly, physical resolution, and modality assumptions need explicit contracts

Evidence: [bins.py](../python/haplodose/bins.py), lines 27–98; [io.py](../python/haplodose/io.py), lines 148–156; [CLI](../python/haplodose/cli.py), lines 15–17 and 39–40; [bulk median mode](../python/haplodose/bulk.py), lines 98–106 and 136–151.

The package accepts arbitrary gene positions but defaults to hg38 centromeres without an assembly declaration. Bins can cross centromeres or large unexpressed gaps, and expression-HMM transitions depend on bin count, not physical distance. Make genome build, interval convention, gene identifiers, duplicate handling, and maximum genomic gaps explicit. Report resolution in both expressed genes and base pairs.

The core NB/PLN design is suited to count data; modality assumptions still differ between UMI scRNA, full-length amplified read counts, mixed spatial spots, and bulk RNA. A spot containing several populations is not a single tumour cell. Optional bulk median mode converts median log ratios into fractional pseudo-counts and then applies a count likelihood; its posterior scale is an approximation requiring separate validation. Validate per modality instead of interpreting successful execution as evidence of transferability.

## 2. Logical interpretation of the model and benchmarks

| Component or output | Supported interpretation | Limitation to preserve |
|---|---|---|
| Normalized cell likelihood and median-centred profiles | Relative RNA composition against the fitted reference | Uniform scaling of every genomic bin is unidentifiable; neither absolute ploidy nor whole-genome doubling follows from expression alone. |
| `beta` and states 1–4 | A constrained dosage-response scale and relative state representation | Purity, dosage response, baseline choice, and centring are not separately identified by this fit. Avoid interpreting beta as a purity estimate. |
| Reference mixture fitted to the specimen | A useful adjustment for known reference types | Tumour dosage aligned with a reference-type contrast can be absorbed by the mixture. Challenge this with truncal events and mismatched references. |
| High-passed programmes | Suppression of some localized expression structure | Whole-arm programme effects remain available to the dosage model by design. Focal events can compete with programmes. Neither direction is guaranteed to be biological CN. |
| Alternating clone fitting | A practical combination of NB assignment, pooled PLN profiling, and regularized updates | The heterogeneous objectives and surrogate score updates do not establish exact monotone EM for a single stated likelihood. Document the approximations and stopping criterion. |
| `normal_like`, `malignant_call`, `tree.tsv` | Reference-like membership, an operational RNA call, and an average-linkage similarity tree | These are not independently validated normal identity, histopathology, or an evolutionary phylogeny. |
| SPECTRUM clone-difference AUROC | Agreement about which arms vary among populations within a patient | Maximum pairwise differences discard which alterations co-occur in which clone. Different clone partitions can have the same score. This cannot validate per-cell clone identity. |
| DNA-derived benchmark clones | A stronger source of structural evidence than RNA-to-RNA agreement | Inferred DNA truth has depth, ploidy, passage, sample/site correspondence, and clustering uncertainty. Separate-cell DNA/RNA does not supply same-cell labels. |
| Seed stability | Repeatability under a fixed dataset and implementation | Identical clone counts do not imply identical membership, accurate boundaries, or calibrated uncertainty. Report label agreement and uncertainty in addition to counts. |

These distinctions are especially important for the internal CRC/IGF2 question. Establish tandem duplication and genome-wide gain burden from DNA; use independent patients as biological n, account for repeated samples, and state whether purity, ploidy/WGD, subtype, and assay are confounders or part of the mechanism being studied. RNA grouping can support investigation of an association, but does not establish that IGF2 duplication caused the gains. Keep the package's public scope general.

## 3. Computational efficiency

The existing implementation already addresses several major bottlenecks: sparse native aggregation, disk-backed bin counts, chunked likelihoods, bounded proposal subsampling, sparse depth/pair weights, cached agglomeration results, randomized programme SVD, and a GPU-oriented HMM scan. Preserve those choices unless profiling identifies a better tradeoff.

| Priority | Current cost and source | Focused improvement |
|---|---|---|
| High | `gene_loader` repeatedly opens HDF5 and reads enclosing blocks; proposal checks and final gene summaries add reads after the advertised two passes. Interleaved specimens can repeatedly touch the same blocks. `pipeline.py:280–299,343–356`. | Instrument bytes/blocks read and time by stage. Consider one streaming pass for gene-level summaries or a bounded block cache keyed by input identity. Preserve row order and pooled-count invariants. |
| High for null-enabled runs | The null materializes cell-by-bin means, simulated counts, and proposal arrays up to 20,000 cells; repeating the search multiplies compute. `clones.py:1199–1278`. | Set an explicit memory budget for simultaneous arrays; cache deterministic reusable quantities; distinguish proposal, simulation, HMM, and I/O time before optimizing. Preserve the tested decision rule. |
| Medium | All-pair clone comparison scales quadratically in clone count; sparse cell-by-pair weights still scale approximately with cells times participating clones. `clones.py:610–612,691–713`. | Batch pairs when necessary and drop stored zero weights. Balance bounded working memory against rereading counts. Preserve the pair-result cache already present in agglomeration. |
| Medium | Thin mode reconstructs deterministic beta-binomial draws repeatedly, using a Python loop per cell. `clones.py:166–184,248–274`. | Profile first; consider caching the active proposal subset or one on-disk half. Preserve row-keyed reproducibility and exact nonnegative count conservation. |
| Medium | `post_store` retains all specimen responsibilities, and profile tables accumulate Python dictionaries before writing. `pipeline.py:105–107,182–198,341,391–395`. | Stream completed specimen outputs or write bounded columnar blocks; preserve final column-to-clone mapping. |
| Lower until measured | Bulk refits reference programmes for each leave-one-out normal. `bulk.py:145–151`. | Profile this against repeated HMM/grid fits before replacing SVD. The matrices have normals by bins dimensions, so SVD is not automatically the dominant cost. |

**[Estimate]** At 20,000 proposal cells and 3,000 bins, one float64 array occupies 0.447 GiB; several coexist in null/proposal operations. At one million cells and 3,000 bins, the float32 bin-count file occupies 11.18 GiB on disk. A memmap controls allocation, but page residency, I/O, and temporary copies still matter. Kernel `max_elements` budgets a nominal tensor, not the sum of all live temporaries.

For efficiency claims, retain identical inputs, pinned code and environments, controlled CPU generation and threads, failure/OOM accounting, and repeated timings. Larger scaling inputs contain resampled duplicate cells ([prepare.py](../benchmarks/scaling/prepare.py)); this tests computational scaling, not recovery of additional biological diversity. Include specimen count, bins, clones, null draws, and programme count alongside cell count.

The scaling summary takes the larger of Slurm MaxRSS and `/usr/bin/time` memory. Verify what the site's accounting plugin includes before equating this with total process-tree peak memory; Slurm defines MaxRSS across tasks and its measurements depend on accounting configuration. The installed `sacct --noconvert` returned integer bytes in the checked example, so I did not find a unit-parser defect there. GPU comparisons also need device-memory measurements. [Slurm accounting documentation](https://slurm.schedmd.com/sacct.html).

## 4. Path toward the manuscript

**[Inference]** The strongest contribution would explain when separating pooled dosage estimation from cell assignment is useful, where it fails, and which components improve that tradeoff. A quota of benchmark wins is less informative than a frozen, reproducible claim tied to a biological estimand.

| Order | Work | Evidence needed before making the corresponding claim |
|---|---|---|
| 1 | Repair finite-value handling and validate input identities/counts; make scoring coverage and required comparators explicit. | Targeted regression cases for F1 and F5–F9; no silent NaN assignments or unsupported success flags. Preserve the original frozen results and label corrected sensitivities. |
| 2 | Evaluate the complete clone procedure on justified nulls and realistic alternatives. | Specimen-level false splits with intervals; rare-clone power by cell count, depth, effect size, expression state, reference mismatch, allele coverage, and spatial mixing. Keep truth quality visible. |
| 3 | Establish a focused method contribution through matched ablations. | Paired biological-unit comparisons of pooled sharing, depth terms, programmes, reference choices, alleles, spatial support, and null gating against simple baselines. Include power lost when false splits are reduced. |
| 4 | Freeze a final method and evaluation protocol for the remaining independent question. | Exact data and code identities, required comparator family, common scoring domain, failure rules, primary effect and uncertainty, multiplicity family, and biological n. Previously exposed lines remain development/overlap evidence even with a new assay. |
| 5 | Validate at the resolution claimed and generate manuscript artifacts. | Aggregate DNA for aggregate profiles; event-level truth for events; same-cell or independently anchored evidence for membership. Reproducible figure/table commands and environment inventories. |

The current records already preserve valuable negative results: MIX-seq did not meet the prespecified profile-superiority claim, the memory claim failed at 20k cells, and the newer clone gate did not satisfy all its acceptance criteria. Preserve those outcomes. The two bulk datasets provide a repeated task-specific signal, but a second cohort for the same endpoint is replication of that endpoint, not a new biological capability.

Documentation should be reconciled after scientific decisions, not used to silently change them. The README's old failure count and the manuscript outline's pending results lag the newer STATUS entries. The outline also associates separate-cell Andor evidence with cell-level clone identity too closely. Keep a compact evidence ledger with frozen version, dataset exposure, unit, endpoint, result, uncertainty, and limitation. Public bulk outputs should record the package/environment and source/input identities as consistently as single-cell outputs. Benchmark scripts still contain absolute workspace paths and runtime imports from shared scorers; portable configuration and hashes of those transitive dependencies are needed for an external reproduction.

## 5. Verification and provenance

Review coverage: all package Python modules, Rust aggregation source and wrappers, the test suite, status/protocol/manuscript documents, and selected Kinker/MIX-seq, Andor/SPECTRUM, CCLE/CPTAC, comparator, scaling, and Schmid integration code. External comparator implementations and every download/launch script were not exhaustively audited. The Wang manuscript exemplar was referenced in the project instructions; its DOI did not resolve through the browsing tool, and this review does not claim to have evaluated its full text.

Checks ran on compute host `ne1dc5-003.nygenome.org` inside existing allocation `22345309` (8 CPUs, 90 GB). The Python suite was restricted to four CPUs and four numerical-library threads; synthetic probes to one separate allowed CPU. No job was submitted, no dependency installed, and no production source changed.

| Check | Outcome |
|---|---|
| Python source syntax | 116 package, test, and benchmark files parsed successfully with `ast.parse`. |
| Existing Python suite | **95 passed, 1 skipped in 368.77 seconds**, exit 0. The GPU comparison was skipped because no GPU was visible; native-versus-SciPy checks ran against the existing extension. |
| EM component removal | Ten of 40 synthetic cells acquired nonfinite programme scores, likelihoods, and responsibilities. |
| Count, metadata, and scoring probes | Fractional counts, mixed specimen identities, missing-arm AUROC, one-gene window truth, and incomplete comparator-family success reproduced. |
| Segment precision | Omitted-reference-variance example independently calculated. |
| Reported bulk endpoints | n, means, signed-rank p-values, and Holm adjustment recomputed from saved CCLE/CPTAC per-case recalls; matched primary JSON. |
| Raw biological analysis | Not run. No unopened SPECTRUM/Bian inputs or new held-out outcomes inspected. |
| GPU execution and Rust rebuild | Not run. The Python suite uses the existing native extension where available. |
| Source preservation | All 133 fingerprinted baseline files remained unchanged through the review. |
| Report checks | Local links resolved, code fences balanced, no trailing whitespace; `git diff --check` passed. |

Evidence is saved under [dev/scientific_review_20261103](../dev/scientific_review_20261103/):

- [probes.py](../dev/scientific_review_20261103/probes.py) and [observations.json](../dev/scientific_review_20261103/probe_outputs/observations.json): reproducible synthetic counterexamples.
- [pytest.log](../dev/scientific_review_20261103/pytest.log) and [pytest_result.json](../dev/scientific_review_20261103/pytest_result.json): suite output, command, exit status, and elapsed time.
- [environment.json](../dev/scientific_review_20261103/environment.json), [source_sha256.json](../dev/scientific_review_20261103/source_sha256.json), and [git_status_before.txt](../dev/scientific_review_20261103/git_status_before.txt): environment, 133 initial source/document fingerprints, and the user's dirty-tree baseline.
- [published_endpoint_recheck.json](../dev/scientific_review_20261103/published_endpoint_recheck.json): independent arithmetic checks on reported bulk endpoints.

Environment: Python 3.11.18, NumPy 2.4.6, SciPy 1.17.1, pandas 2.3.3, JAX/JAXlib 0.10.2, scikit-learn 1.9.1, pytest 9.2.0. HEAD was `54b0cf2a4b94c1dfa42ae4a5e792e88bd0a4b2e2`; hashes identify the reviewed working files beyond that commit.

Reproduce the bounded checks on an allocated compute node, with an unused output directory for probes:

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=python \
  OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 RAYON_NUM_THREADS=4 \
  JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false taskset -c 0,2,4,5 \
  /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
  -m pytest -q -p no:cacheprovider tests

env PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=python \
  /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
  dev/scientific_review_20261103/probes.py /path/to/new/review_probe_output
```

The CPU IDs above are those used in this allocation; select allowed CPUs in a different allocation. The review runner imposed a 600-second suite time limit. Suggested scientific calibration and ablation experiments are recommendations, not completed checks.
