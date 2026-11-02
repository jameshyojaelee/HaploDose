# HaploDose: fresh mathematical, algorithmic and efficiency audit

Reviewed 2026-11-02. This report inspects the working tree after the November 1 reviews and the subsequent shared-profile, proposal-validation and gene-level-output changes. It does not assume those earlier findings are still unfixed. Source fingerprints and executed checks are below. Review only: no package or benchmark implementation was changed.

## Main conclusions

HaploDose has improved. The newest field comparison supports strong single-line RNA dosage profiles and an advantage over CopyKAT and SCEVAN on the current implant design. It does **not** establish superiority over infercnvpy for subclone recovery, or over all packages across assays and endpoints.

The highest-priority new problem is in `shared_profiles`: collapsing all pairwise event intervals to a Boolean mask can erase the boundary between distinct events. An actual-function numerical probe reduced a well-supported contrast between two tumour clones by approximately 95% when a third clone's event bridged their boundaries. The new gene-level reporting path also uses a different reference composition from the fitted model and can generate a large apparent loss in a perfectly matched, copy-neutral construction.

Keep Python as the public interface and orchestration layer. Existing JAX and Rust already move substantial numerical work outside interpreted Python. The main opportunities are preserving information, improving local-event proposals, making model components agree, bounding pairwise work and reducing repeated reads. A whole-language rewrite would preserve the current mathematical problems.

Real-cohort consequences of the new failure cases are **[Unverified]**. The report distinguishes reproduced behavior, analytical consequences, existing benchmark artifacts and proposed improvements. No new cohort inference or competitor run was performed.

## What changed since the previous audit

Sources: [current clones.py](../python/haplodose/clones.py), [pipeline.py](../python/haplodose/pipeline.py), [previous follow-up](efficiency_review_2026-11-01_followup.md), and its [subsequent critique and experiments](efficiency_review_2026-11-01_followup_critique.md).

| Earlier finding | Current status |
|---|---|
| Disjoint rare groups discarded by raw agreement | Fixed: minority-set Jaccard replaces raw agreement. |
| Only the first two children checked before accepting a multi-child split | Fixed: every pair is checked; the weakest pair determines the score. |
| Sliding windows add cost without demonstrated gain | Default is now off. Do not repeat the earlier recommendation as unfinished work. |
| Reported differences bypass acceptance vetoes | Fixed: `summarise` now calls `supported`. |
| Small clones lose shared background dosage during EM | New `shared_profiles` addresses this; findings 1–2 below concern its implementation. |
| Highly expressed genes dominate summed-bin reports | New `log_ratio_gene` output addresses this; finding 3 concerns its reference construction. |
| Reference variance absorbed into tau and then added again | Still present. The critique records an ablation without a clear benefit and with increased spatial neutral-arm calls. Do not assume correcting this alone improves accuracy. |
| Global mismatch-gene exclusion across datasets | Still present. Dataset-specific/no-mask variants did not improve the recorded Kinker experiment. This remains a generalization issue, not a proven immediate performance win. |
| Mirrored allelic events, offset inconsistency, expression search excluding clone 0 | Still present. These are carried-forward issues, not new discoveries. |

The earlier suggestion that mirrored-allele comparison necessarily needs a new shared phase model was too broad. Counts already share the input SNP orientation across groups. A signed per-SNP comparison can exploit that orientation. The parked heterogeneity implementation still needs low-depth calibration and proposals that can discover opposing subclones when their pooled allele fraction is balanced.

## New findings

### 1. Shared-profile masking can erase strongly supported event boundaries

**Priority: high. Reproduced using the actual implementation.**

Source: [clones.py](../python/haplodose/clones.py), `shared_profiles`, lines 384–429, especially 412–426.

The function unions all pairwise difference segments into one Boolean `mark[k]` per clone. It then averages the clone-to-pool deviation over each contiguous marked run. The mask retains neither the identity of the pair nor the internal event boundaries. A third clone can connect two intervals that distinguish the first two clones in opposite directions.

Probe: 240 bins, one chromosome, no programmes or depth effects, and four pooled groups including normal-like group 0. All counts start at 100,000 per bin. In bins 40–79 and 80–119, respectively:

| Group | First interval | Second interval |
|---|---:|---:|
| Normal-like | 100,000 | 100,000 |
| Tumour 1 | 180,000 | 120,000 |
| Tumour 2 | 120,000 | 180,000 |
| Tumour 3 | 60,000 | 60,000 |

The actual difference HMM detects tumour 1 versus tumour 2 as two opposite events, each with posterior 1.0, mean log difference ±0.40546 and |z| ≈ 89.13. Each also differs from tumour 3 across the whole interval. Those longer differences join the masks.

Maximum tumour-1 versus tumour-2 contrast is 0.09916 in the independently segmented doses, but only 0.00503 after shared-profile replacement. The residual contrast mainly reflects boundary/weighting details; most of the biologically distinct pattern is averaged away. This is a deliberately high-depth structural probe, not a patient effect estimate.

**Improvement:** form a common partition from every relevant event boundary, including changes in direction, pooled-profile boundaries and clone-profile boundaries. Estimate a deviation within each atomic interval. Adding a third clone must not remove an already resolved contrast between two existing clones. Test adjacent events, overlapping events, opposite directions and clone-order permutations. Longer term, fit a pooled background plus regularized clone-specific deviations on a shared partition, rather than editing doses after separate HMM fits.

### 2. Shared doses bypass some evidence filters and differ from reported doses

**Priority: high for model consistency. Confirmed by source inspection; cohort impact unmeasured.**

Sources: [clones.py](../python/haplodose/clones.py), lines 409–426, 468–517, 665–714 and 1008; [pipeline.py](../python/haplodose/pipeline.py), lines 314–379.

Three paths now have different meanings:

1. Shared-dose intervals come from `differences` on soft pools, without the depth-matched pools and gene/allele vetoes used for final split acceptance.
2. `differences` can subtract programme effects, but the deviation added back in `shared_profiles` is calculated from raw pooled log ratios and depth expectations. It does not subtract the fitted programme contrast. **[Inference]** nuisance expression can therefore enter the dose inside an interval selected for a real event.
3. Cell assignments use shared doses, while `fit_specimen` finishes with a fresh call to ordinary `profiles`. The `profiles.tsv` and `segments.tsv` outputs describe the independent fits, not the shared dose that generated the last expression likelihoods. Even before the change, group-0 profiling and its fixed zero assignment dose had different semantics; shared profiles extend the discrepancy to tumour clones.

Distinct estimators can be useful, but their roles must be explicit. A reported profile should not silently appear to be the profile used to classify cells when it is not.

**Improvement:** define one model-result contract containing the dose used for assignment, the independent diagnostic profile, programme correction, accepted event partition and final responsibilities. Reuse an evidence evaluator with explicit modes for soft-pool fitting and hard-group validation. Apply a programme-corrected contrast when constructing shared doses. Add an invariant that stored assignment doses reproduce the reported cell likelihoods, conditional on the stored nuisance parameters. Assess a stable objective or held-out predictive score; the current EM stopping rule only checks hard-label changes below 0.5%, so unchanged labels do not establish parameter convergence.

### 3. The new gene-level output does not use the specimen's fitted reference mixture

**Priority: high because the new output is used in the headline profile comparison. Reproduced helper behavior and confirmed call path.**

Source: [pipeline.py](../python/haplodose/pipeline.py), lines 80–82, 142–165 and 242–266.

The bin model fits specimen-specific reference weights. However, `ref_gene[dataset]` sums raw gene counts over all fit reference groups, and `gene_shares` uses that pooled composition to split each bin expectation among genes. It ignores the model's mixture weights and any added in-sample reference. Consequently, `E_bin` and its within-bin gene shares can describe different reference populations.

Probe: five bins, three genes per bin. Query counts exactly match reference type A. Four bins have gene counts (100,100,100); the remaining bin has (240,30,30). Type B matches the other bins but has (30,135,135) in the changed-composition bin. Both types have exactly the same bin totals. With equal raw-count pooling of reference types, the gene-level helper reports **−1.00164 natural-log units** in the changed-composition bin despite a zero bin-level log ratio. Using type A's gene shares gives −0.01135, the small pseudocount effect. There is no CNA in the construction.

The identical bin totals in that first probe do not identify the reference mixture from bins alone. A second probe makes type A identifiable: place the changed-composition bin third, and set type B's other per-gene counts to 10, 20, 300 and 70 in bins 1, 2, 4 and 5. Both types still have total 1,500. The actual `mixture_weights` function assigns **0.99818 to type A**. Using that fitted bin expectation, centered bin ratios range only from −0.00418 to +0.00109, yet the gene-level output still reports **−1.00218** in the third bin. Thus this is not merely an unidentified-mixture example.

**Improvement:** retain per-reference-group gene rates and mix them with the specimen's fitted weights before calculating conditional within-bin shares. Handle the in-sample reference explicitly. Verify that a pure reference-type query yields a near-zero gene profile regardless of how many unrelated reference cells were supplied. Re-evaluate mixed-reference cohorts; single-line gains from a pooled reference do not validate this case.

The current helper also aligns the medians only on bins containing eligible genes. Its docstring promises equality of the overall bin median. With fallback bins, those statements need not agree. Define and test the intended centering convention.

### 4. The HMM's initial prior is not its long-run occupancy prior

**Priority: medium. Exact mathematical consequence; intended behavior needs an explicit decision.**

Sources: [kernels/hmm.py](../python/haplodose/kernels/hmm.py), `sticky_transition`, lines 148–153; [profile.py](../python/haplodose/profile.py), `STATE_PRIOR` and `DIFF_PRIOR`.

The transition uses `T_ii = 1-p` and `T_ij = p*q_j/(1-q_i)` for i ≠ j. Its stationary distribution is proportional to `q_i*(1-q_i)`, not `q_i`. This follows directly from detailed balance and was checked numerically against the actual matrix.

| Model | Supplied neutral probability | Stationary neutral probability |
|---|---:|---:|
| Copy-number HMM, q=(0.2,0.5,0.2,0.1) | 0.5000 | 0.3788 |
| Difference HMM, q=(0.15,0.7,0.15) | 0.7000 | 0.4516 |

These are long-chain, uninformative-emission limits, not measured false-call rates. The current construction deliberately gives states equal exit rates, so it is mathematically valid if that is the desired prior. It is not a model with 70% long-run neutral occupancy in difference tracks.

**Improvement:** state separately the intended occupancy distribution and expected segment lengths. If occupancy q is intended, a redraw transition `(1-p)I + p*1*qᵀ` preserves q, but changes state-specific exit probabilities and must be recalibrated. A distance-dependent transition can also avoid giving sparse and gene-dense genomic intervals different physical-length priors. Validate occupancy, segment lengths and null event calls together; do not switch the formula merely to increase sensitivity.

### 5. The current implant benchmark contains minorities below the proposal size floor

**Priority: high for interpreting and improving rare-clone detection. Verified from current score tables.**

Sources: [clones.py](../python/haplodose/clones.py), `min_cells`, `consolidate`, `event_mixtures` lines 858–860; [implant score table](../results/kinker_field_20261102_0529/implant/per_specimen.tsv).

**29 of 166 implanted specimens contain fewer than 30 implanted cells**, all in the 15% settings. The smallest minority has 15 cells. In the 500-gene, 15%-gain setting, 8 of 13 specimens fall below the default floor; the median minority is 29 cells.

An exact true binary event proposal with fewer than 30 minority cells fails the explicit size check, even if its region is known. This is a structural restriction before asking whether the HMM can resolve the event. It is not a proof that the final soft-assignment model can never produce a smaller MAP group; responsibilities and later assignments can change membership.

**Improvement:** separate minimum proposal membership from the evidence requirement for retaining a clone. Consider effective information or informative UMIs, with a stricter calibrated false-split rule for small groups. Report performance stratified above and below the configured floor. Preserve all 166 specimens in the primary benchmark; changing the evaluation population after inspecting failures would obscure the limitation. The previous 20,000-cell proposal cap remains a separate obstacle for very rare clones in large specimens.

### 6. Entirely unassigned specimens receive perfect ARI in the implant scorer

**Priority: high for benchmark reliability. Reproduced. Not affecting the inspected field table.**

Source: [benchmarks/kinker/score_implant.py](../benchmarks/kinker/score_implant.py), lines 24–35.

The scorer removes unassigned cells and passes the remaining arrays to `adjusted_rand_score`. With no remaining cells, the installed sklearn returns 1.0. Calling the actual scorer with a two-cell specimen whose labels are both `NA` produced `n_scored=0`, `clones=0`, `ari=1.0`.

All four tools in the inspected 166-specimen field table have `n_scored/n_cells = 1.0`; there are no empty scored specimens. This bug does **not** explain their current ARI ordering.

**Improvement:** fail or explicitly mark insufficient coverage before scoring, retain unassigned-cell counts and publish a prespecified failure-aware endpoint. The same common-cell denominator must be used across tools. Do not silently exclude missing specimens from aggregate comparisons. A regression should include empty, partially assigned and completely missing outputs.

### 7. Spatial graph construction can leak each cell into its own prior

**Priority: medium; spatial mode is secondary to bulk and scRNA. Reproduced.**

Source: [clones.py](../python/haplodose/clones.py), `neighbour_graph`, lines 1012–1023; `choose_lambda`, lines 326–341.

The graph builder discards the first returned nearest neighbour, assuming that row is the cell itself. At duplicated coordinates, tied neighbours can appear in another order. The graph then contains self-edges, violating the stated exclusion of a cell's own data from its prior.

A deterministic probe with 60 cells at 20 triplicated coordinates and k=1 produced 20 self-edges. With random binary cell labels (seed 81), likelihood separation 2 and equal priors, `choose_lambda` selected 1.0. Removing the diagonal in that diagnostic selected 0.0. This isolates the leakage; it is not a full spatial accuracy benchmark.

Two additional mathematical issues deserve explicit treatment:

- The kNN graph is directed and row-normalized. For a conventional symmetric Potts energy, the mean-field gradient uses symmetric edge weights. `G @ R` on an asymmetric graph is not generally that gradient. Either construct a symmetric weighted energy or document this as a smoothing heuristic without a variational-objective guarantee.
- Spatial event proposals add NB sizes in proportion to neighbour count. Independent NB counts with unequal means and common dispersion do not, in general, sum to that NB distribution. For means 10 and 100 and phi=0.2, the true summed variance is 2,130; the implemented equal-mean approximation gives 1,320. The reduction in uncertainty is overstated. With independent cells, a moment-matched effective size is `(sum m)^2 / (phi*sum m^2)`. Biological/spatial covariance would require an additional term.

**Improvement:** remove self indices explicitly, handle coordinate ties deterministically, decide the graph-energy contract and match pooled variance to heterogeneous depths. Test disconnected slices, duplicate coordinates, variable density and null spatial labels.

### 8. Three-normal bulk panels cannot support the implied stability of a fixed z threshold

**Priority: medium for bulk calibration. Analytical concern with a controlled simulation; full-pipeline rates not measured here.**

Source: [bulk.py](../python/haplodose/bulk.py), lines 98–100 and 149–158.

Bulk accepts as few as three normal samples, estimates an arm-specific center and MAD scale, and calls at |z|>3. The 1.4826 MAD factor is a population-normal consistency factor, not a finite-panel predictive calibration. Estimated centers and scales are uncertain, especially with three normals.

Controlled null: independent standard-normal panel observations and one independent standard-normal target, 50,000 repetitions per panel size, seed 20261102, code's median/MAD and 0.001 SD floor, and inflation factor held at 1. Observed fractions exceeding |z|=3 were:

| Normal panel size | Null exceedance |
|---|---:|
| 3 | 30.634% |
| 10 | 5.252% |
| 30 | 1.294% |
| 100 | 0.496% |

These are Monte Carlo estimates for the isolated threshold, not HaploDose cohort false-positive rates. Actual leave-one-out normal profiles are dependent and sample-tau inflation can change the result. Even under an ideal Gaussian null, treating the estimated scale as known is unsafe. The package does not currently promise that z=3 gives a particular FDR; it should not acquire such a claim without calibration.

**Improvement:** use predictive calibration that includes reference estimation, shrink arm scales across an appropriate normal panel, and return uncertainty when calibration is underpowered. Calibrate the declared family of arm/segment hypotheses, rather than interpreting each z threshold as a genome-wide error rate. Do not lower thresholds simply to improve recall.

A separate residual leave-one-out issue exists in optional `bin_stat='median'`: gene-level pseudo-counts are constructed using **all** normals at lines 114–118 before the per-normal exclusion at line 127. Refitting `reference_model` on those precomputed pseudo-counts does not remove the held-out normal from the transformation. Default summed-count mode does not have this particular leakage. Recompute reference-dependent transformations within each fold.

### 9. Cell likelihoods still lose small contrasts during float32 summation

**Priority: medium for numerical robustness. Reproduced; not established as the cause of CPU-dependent cohort results.**

Source: [kernels/estep.py](../python/haplodose/kernels/estep.py), `_chunk_loglik` and `posteriors`.

The earlier audit fixed HMM message cancellation. Cell likelihoods still sum large common terms separately for each clone and then subtract them during posterior normalization.

Probe: 2,000 uniform bins, phi=0.1, two clone rates differing by +0.01 log units in the first 20 bins before normalization, counts exactly at clone 0's mean. At library sizes 20,000 and 200,000, the float64 reference log-likelihood differences were −0.004950 and −0.008976. The current function returned **0.0** for both differences. These are small individual effects, but they demonstrate loss of an identifiable contrast.

**Improvement:** accumulate likelihood ratios relative to a reference clone or stable per-bin baseline before the large reduction. For a clone/reference log-mean difference delta, use `log1p((m_ref/(r+m_ref))*expm1(delta))` for the NB denominator ratio. Selective float64 reduction is another option; measure its cost. Add invariance tests over library depth, bin count and clone permutation. Hardware-repeatability checks should report margin distributions and discrete model-selection changes, not just whether a seed was set.

## Computational inefficiencies in the current architecture

These are source-based scaling observations unless a measured result is explicitly given. No new full-scale timing benchmark was run.

| Path | Current cost or duplication | Improvement and validation |
|---|---|---|
| `shared_profiles` every EM round | Reprofiles all groups, profiles the tumour pool, then evaluates every tumour pair. With 40 tumour clones this adds 780 difference-HMM tracks per update, up to 6,240 across eight updates, before merge/report tests. | Use a shared interval partition and regularized clone-to-pool contrasts, reserving exact pairwise tests for ambiguous neighbours. Batch pair computations under a memory budget. Verify pairwise event preservation first. |
| `shared_profiles` and `fit_amplitudes` | Each pair can fit programme amplitudes separately; overlapping comparisons repeat matrix work. | Cache reference-only matrices; batch weighted solves when weights allow. Pool-dependent results must be invalidated when responsibilities or nuisance parameters change. |
| New gene-level reporting | `summarise` rereads raw gene counts per specimen. `gene_loader` reads enclosing storage blocks, so interleaved specimens can cause repeated reads of mostly irrelevant rows. | After fitting, aggregate genes for all specimens in one streaming pass, or use a reusable sparse used-gene cache. Keep per-group reference statistics from the original pass. Measure bytes read and read amplification, not only CPU time. |
| Retained cohort results | `post_store` holds every specimen's responsibilities until the final output; profile tables accumulate as Python dictionaries. | Stream tables and per-specimen posteriors, then consolidate deterministically. **[Estimate]** 1M cells × 40 float32 responsibilities alone is about 153 MiB; the out-of-core count matrix does not make all outputs bounded-memory. |
| Proposal region search | Each region fits up to 200 mixture iterations, each using repeated NB rate solves. Extra overlapping windows already failed to pay for themselves in recorded 100k tests. | Screen windows with cached residual/information prefix sums, then fit the exact count model on a short list. Record scan, fit and validation costs separately. |
| `pooled_depth_slopes` | Dense cells × depth-stratum/clone one-hot arrays remain. | Direct grouped accumulation or sparse labels, including first and second depth moments. This earlier recommendation remains applicable. |
| `merge_unsupported` | Recomputes all depth-matched pair tests after each merge; `agglomerate` already has a pair cache. | Cache unchanged hard-group pairs within a merge pass and invalidate merged groups and nuisance-dependent entries. Avoid rebuilding all sparse weights. |
| Bulk beta search and reference fits | Grid search repeatedly computes posterior passes even when only likelihood is needed; each leave-one-out reference uses a full SVD. | Forward-only likelihood ranking, batched beta candidates, truncated SVD with verified exclusion. Preserve fold-specific feature transformation for median mode. |
| EM stopping and nuisance updates | Hard-label stability can stop with moving probabilities/parameters; five programme-score steps restart from zero. | Warm-start nuisance scores and use a residual/objective tolerance with an iteration cap. Measure scientific fit and total work together. |

The new all-pairs split validation is an important correctness improvement. Its cost is a reason to bound or stage proposals, not to restore first-pair-only validation. `max_groups` caps a Leiden proposal, not obviously the total number of clones accumulated through repeated splitting; enforce and report a global resource budget separately from a biological evidence threshold.

**Python decision:** retain the current hybrid architecture. Compile or fuse a hotspot only after profiling a representative workload with counts, clone count, programmes and input layout recorded. Native prefix-sum/grouped-aggregation kernels may be useful after the algorithm is settled. JAX already compiles the cell and HMM kernels; GPU benefit remains unmeasured here. A 0.3-second microprofile from the previous review is insufficient evidence for a language migration.

## What the latest comparisons actually establish

Sources: [single-line field summary](../results/kinker_field_20261102_0529/single_line/summary.json) and [implant field summary](../results/kinker_field_20261102_0529/implant/summary.json), read directly in this review.

| Method/output | Single-line median arm rho, 181 lines | Implant mean ARI, 166 lines |
|---|---:|---:|
| HaploDose gene-level output | 0.8613 | 0.1672, using current clone model |
| HaploDose original bin output | 0.8545 | Same clone model |
| Pseudobulk baseline | 0.8634 | Not a clone-discovery method |
| infercnvpy | 0.8126 | 0.1651 |
| CopyKAT plus wrapper clustering | 0.7651 | 0.0886 |
| SCEVAN, native subclone labels | 0.7098 | 0.0273 |

The single-line artifact applies paired Wilcoxon tests with Holm correction over five comparisons against HaploDose. HaploDose exceeds infercnvpy, CopyKAT and SCEVAN on that endpoint; pseudobulk is slightly ahead (Holm p=0.0317). The implant artifact uses Holm over six tool pairs. HaploDose versus infercnvpy has p=0.910, while its comparisons with CopyKAT and SCEVAN have adjusted p≈0.00147 and 0.00000753. These are artifact results, not independently rerun inference. Repeated use of Kinker to choose changes makes it a development benchmark.

Several evaluation details matter before treating this as field-wide superiority:

- All tools receive other cancer lines as reference cells. This is a useful relative-profile stress test, but it differs from the diploid-normal setting expected by several tools. Include an independent normal-reference benchmark.
- `evaluate.py` forms DNA arm truth and RNA arm scores from **each run's used genes**. `compare_tools.py` intersects specimen IDs, not the underlying genomic evaluation units. Freeze common genomic coordinates and truth for a fair primary comparison; retain native-resolution results as a secondary view.
- `run_r_field.py` maps comparator rows to the nearest genomic start, then averages 50-gene bins. For window outputs, interval overlap and native boundaries are preferable to nearest-start interpolation. Record assembly and interval conventions explicitly. This review did not establish an assembly mismatch.
- With `--rule all`, missing comparator cells are filled with zero profiles and given labels. Therefore 100% scored-cell coverage does not establish 100% native inference coverage. QC files confirm small native omissions, for example CopyKAT r0/pool9 returned 6,030 of 6,032 sent cells. Report both forms of coverage; do not imply these small omissions explain the large profile gaps.
- `region_found` is any overlap with any reported clone difference. It does not require the correct cells, direction or event length, and the wrappers generate differences using a different rule from HaploDose. The 0.36 versus 0.67 result is a useful diagnostic, not a comparable event-recall estimate by itself.
- Pseudobulk is a strong baseline, not a mathematical ceiling on DNA prediction: a well-validated denoising model could improve on noisy pooled expression.

An independent multi-omics benchmark found different strengths for malignant classification, CNA profiles, breakpoints and cnLOH, including strong Numbat performance. That makes an allele-aware comparator and matched DNA truth essential for a broad claim. [Primary benchmark study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11879432/).

Numbat jointly uses expression, allelic information and haplotype phasing; its benefit cannot be assessed by RNA-count implants that leave the underlying reads unchanged. [Numbat paper](https://www.nature.com/articles/s41587-022-01468-y). CalicoST addresses allele-specific spatial CNA inference, including mirrored events and tumour admixture, which are different endpoints from expression-only clone clustering. [CalicoST paper](https://www.nature.com/articles/s41592-024-02438-9).

## The most promising algorithmic direction

The existing known-region diagnostics point to substantial headroom. I recalculated their setting-level means from [oracle.tsv](../dev/hd1_proposals_20261102_0600/oracle.tsv) and [oracle_even.tsv](../dev/hd1_proposals_20261102_0600/oracle_even.tsv):

| Implanted size | Known-region mixture ARI, all genes, range across settings | Known-region mixture ARI, retained proposal genes |
|---|---:|---:|
| 50 expressed genes | 0.322–0.627 | 0.139–0.418 |
| 150 expressed genes | 0.579–0.842 | 0.401–0.621 |
| 500 expressed genes | 0.801–0.953 | 0.522–0.888 |

These are existing real-count implant diagnostics, not a formal upper bound. The region is supplied, the binomial model ignores overdispersion, and the proposal-gene comparison changes both filtering and the even-bin restriction. It cannot attribute the entire difference specifically to even/odd splitting. Nevertheless, it identifies information discarded before final segmentation.

**[Inference]** The most useful development experiment is a count-aware multiscale proposal system with an explicit false-split calibration:

1. Estimate per-cell/bin residuals under the current depth/programme baseline. Cache NB score contributions `(Y-m)/(1+phi*m)` and information `m/(1+phi*m)`, preserving chromosome boundaries.
2. Use prefix sums or an equivalent multiscale representation to screen local regions cheaply. Rank regions using comparable evidence, not absolute fit gain across arbitrarily different widths. Preserve candidates supporting distinct rare cell sets.
3. Refit shortlisted regions with a composition-aware count model and a background complement. For an ideal conditional-binomial formulation, a relative regional multiplier theta changes expected fraction p to `theta*p/(1-p+theta*p)`. Actual RNA data require overdispersion/covariance treatment; the oracle's binomial assumption is not automatically valid for deployment.
4. Calibrate the complete adaptive search, including region choice, candidate ranking, repeated splitting and nuisance estimation. If using all bins to regain power, ordinary Poisson thinning does not automatically create independent validation data for overdispersed or correlated RNA counts. Use an appropriately derived count-splitting model or null resampling that preserves the relevant structure, and verify calibration on biological null samples.
5. Pass accepted intervals to a shared-partition clone fit that preserves their boundaries. A better proposal is wasted if a later update averages the distinguishing event away.

For bulk, keep a robust broad-scale gene-level profile and assess whether segmentation improves focal-event accuracy at fixed false-call rates without degrading arm profiles. The current four-state model, single dosage slope and normalized RNA counts do not identify arbitrary purity/ploidy solutions. Do not force a single point copy number where multiple DNA configurations are observationally compatible. These are scientific constraints on what a package can win, not Python limitations.

## Recommended order and acceptance criteria

| Order | Work | Required evidence before adoption |
|---|---|---|
| 1 | Preserve shared-profile boundaries; align fitted and reported profile semantics; match gene shares to reference weights. | Reproduced cases above pass, clone-order invariance, no loss of existing programme/normal protections; paired cohort ablations. |
| 2 | Repair missing-output scoring and freeze common genomic truth/coverage. | Empty output cannot score perfectly; identical genomic denominator across tools; native failures reported. |
| 3 | Local proposal screening plus exact shortlist fitting; separate proposal-size and retained-clone evidence requirements. | Better rare/focal-event recovery at the same full-search false-split rate; report the 29 below-floor specimens separately and in the overall result. |
| 4 | Decide HMM occupancy/duration priors, finite-normal-panel calibration, stable likelihood contrasts and convergence checks. | Controlled null calibration, parameter recovery, CPU-repeatability diagnostics, no unsupported FDR interpretation. |
| 5 | Optimize repeated pair fits, raw-gene reads and output retention. | Matched hardware/input layout, wall time and peak RSS, I/O measured, identical accepted scientific behavior within declared tolerances. |
| 6 | Improve mirrored-allele discovery/comparison and spatial graph/variance handling. | DNA-backed allele events, unequal-depth nulls, mirrored clones whose pooled BAF cancels, duplicate-coordinate invariants. |

To support a leadership claim, predeclare separate primary endpoints for bulk profiling, expression-only scRNA clone discovery and allele-aware inference. Evaluate an untouched patient/line set after development. Use biological specimen/patient as the unit, cluster repeated implants within line and account for shared reference pools; cells, bins and repeated seeds are not independent patients. Report paired effect sizes and confidence intervals, with Holm correction over the prespecified method/endpoint family. Use a declared exploratory FDR family when screening many variants. Include neutral-sample false splits, failed samples, breakpoint/event precision-recall, runtime, memory and preprocessing costs.

## Checks executed and limits

- Read current source, prior reports, critique, status, relevant tests and benchmark wrappers. Compared six current source hashes with the preceding review; `clones.py` and `pipeline.py` changed while `alleles.py`, `profile.py`, `bulk.py` and `kernels/estep.py` matched its fingerprints.
- Ran the entire existing Python test suite: **61 passed in 166.61 seconds**.
- Ran the actual-function shared-profile and gene-reference probes, graph/self-prior probes, scorer empty-output probe, HMM stationary-distribution check, float32 versus float64 cell-likelihood calculation and finite-panel null simulation described above.
- Read existing field scores and known-region diagnostic tables; recomputed coverage, minority-size counts and oracle setting means. These checks validate the inspected artifacts, not the original biological data or entire comparator pipelines.
- Executed on compute host `ne1dc5-003.nygenome.org`, within existing SLURM allocation 22345309 (8 CPUs, 90G). No new job was submitted. Test suite restricted to CPUs 4–5; numerical probes used CPU 6 and one numerical-library thread.
- No new cohort refits, full-search null simulations, source-fix ablations, GPU runs, comparator inference or Rust rebuild/tests. Those checks are **not run** in this audit.
- Existing unrelated modified/untracked files were preserved. A new unrelated `benchmarks/figures/` directory appeared during the review; it was not created or edited by this audit. Audited package hashes remained unchanged through inspection.
- Report validation: all 20 local links resolved, code fences were balanced, no trailing whitespace was found, and `git diff --check` passed. The only file authored by this audit is this Markdown report.

Test command, from the repository root:

```bash
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
  MKL_NUM_THREADS=2 RAYON_NUM_THREADS=2 JAX_PLATFORMS=cpu taskset -c 4,5 \
  /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
  -m pytest -q -p no:cacheprovider
```

Environment observed: Python 3.11.18; HaploDose 0.1.0; NumPy 2.4.6; SciPy 1.17.1; JAX/JAXlib 0.10.2; scikit-learn 1.9.1; pytest 9.2.0. The local status document's scikit-learn 1.7.2 statement does not describe this interpreter.

HEAD was `27dbbe9613d3a10560ca7800802b063cbdeba55e`, with a pre-existing dirty working tree. HEAD alone does not identify the audited code. SHA-256:

```text
071f6496b466c2e42284e98ad4b7cc8aeb9e7f02cf9fdbe063d88894918cddf6  python/haplodose/clones.py
2d04a761ad6f0722225859ead5bf1f64bc13294455e02e0d1db3ff51ca25ffdd  python/haplodose/pipeline.py
c681e451d608427f587a7fcd552946f875a54738c7ff73b592b7e930a0c16529  python/haplodose/alleles.py
ab13f84d57c72bd2de04e13c108e85ab4b5ffd0e82b89ecdadcec094b508fe8c  python/haplodose/profile.py
6e4fcb34d7c91c8f6c284176200ae41044928608364476ff43445363f54cc223  python/haplodose/bulk.py
726198a0cadf32ab2476ca0d9eb53d83d8f3fb8cf93f9c519672583513f67f45  python/haplodose/kernels/estep.py
b4de6517c2293e7eae0d77226319603852979ef9283e061e22acdea3d940352d  python/haplodose/kernels/hmm.py
```

### Reproducing the most consequential new structural case

Run on an allocated compute node with the same interpreter and thread limits. This constructs pooled profiles directly; it does not simulate patients or estimate biological prevalence.

```python
import numpy as np
from haplodose.clones import CloneParams, Rows, SpecimenModel

J = 240
reset = np.zeros(J, bool)
reset[0] = True
Y = np.full((4, J), 100000.)
Y[1, 40:80], Y[1, 80:120] = 180000, 120000
Y[2, 40:80], Y[2, 80:120] = 120000, 180000
Y[3, 40:120] = 60000
L = Y.sum(axis=1)
XL = np.zeros((4, 2))
m = SpecimenModel(
    Rows(Y, np.arange(4)), L, np.full((J, 1), 1/J), reset,
    np.full(J, .1), CloneParams(n_programmes=0), gamma=np.zeros(J),
)
print(m.differences(Y[1], L[1], XL[1], Y[2], L[2], XL[2])[0])
p = m.shared_profiles(Y, L, XL)
for key in ('dose_own', 'dose'):
    print(key, np.max(np.abs(p[key][1] - p[key][2])))
# Audited version: dose_own 0.09916258; dose 0.005027741.
```
