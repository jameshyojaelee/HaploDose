# HaploDose: deep scientific code review

Review date: 2026-11-10. Reviewed implementation: working tree at Git HEAD `58157f2a3b4e76228c444664c36448571ebc6578`, package version `0.2.0`, including uncommitted package changes. Source identities appear at the end. This is a fresh assessment of the current implementation and newer evidence; the [2026-11-09 review](scientific_code_review_2026-11-09.md) remains a historical record.

**Assessment.** HaploDose has a defensible exploratory use as a regional RNA profiling and candidate grouping package. Its newer profile estimators have measured improvements on exposed, pooled-query regression datasets. The implementation also now resolves the earlier ribosomal-namespace and final allele-snapshot defects. These are substantive advances. The reviewed evidence still does not establish calibrated DNA copy states, genomic clone membership probabilities, specimen-wide biological false-split control, or calibrated uncertainty for the routine segment outputs.

The central scientific obstacle is separating genomic dosage from other mechanisms producing the same RNA observations. Additional obstacles are adaptive selection, nuisance-model dependence, inconsistent objectives across fitting stages, and the gap between pooled-profile benchmarks and cell-level claims. Numerical precision and exact regression reproduction strengthen computational credibility but do not resolve those scientific questions.

Only this document was added by this review. No package code, tests, defaults, historical results, readiness decisions, or claims were changed. No jobs were submitted, packages installed, or current tests executed.

**Scope and evidence discipline.** The review covers the public cellular and bulk workflows; RNA and allele likelihoods; reference, depth and programme fitting; adaptive splitting and merging; confirmation; uncertainty; input identities and native aggregation; the new smoothing, stacking and boundary-averaging modules; selected tests; and selected benchmark/scoring code and saved results. It does not claim exhaustive inspection of every research script, raw sequencing file, development snapshot, or external dataset. Large matrices and DNA truth were not reprocessed.

The initial tree contained 24 modified tracked paths and 831 untracked entries. These are the user's baseline. Relative to the preceding source inventory, four existing package files changed and five Python modules were added. All 26 current Python modules and two native Rust source files match the source snapshot associated with the saved 508-test validation exactly. That correspondence is stronger than relying on the unchanged version number, but it is not a fresh execution in the current environment. Evidence includes bounded standard-library arithmetic on saved tables, explicit derivations, source inspection and primary literature. Reported simulations remain simulations; no new biological specimens were analysed. Concurrent work added a shared-profile merge research script and report during this review; their pending experiment is not evidence of a completed correction, and the package/native source remained unchanged.

**Priority map.** “Blocking” means blocking the specified scientific claim, not exploratory execution. “High” identifies a material inferential or modelling limitation. “Moderate” identifies a bounded implementation, resource, interface or reporting issue. Several limitations are already accurately acknowledged in the README and scientific contract.

| ID | Priority | Finding | Main evidence |
|---|---|---|---|
| 1 | Blocking | Relative RNA observations do not uniquely identify genomic dosage or clones | Algebra and model source |
| 2 | Blocking | The complete adaptive procedure lacks established biological false-split control | Source and conditional simulation evidence |
| 3 | Blocking | Pooled-profile agreement does not validate individual-cell genomic membership | Benchmark units, truth resolution and transfer failure |
| 4 | High | Refinement combines objectives without a common convergence guarantee | Source and preserved transition investigations |
| 5 | High | Reference mixtures and nuisance factors can absorb or preserve the wrong signal | Explicit counterexamples and source |
| 6 | High | Coarse depth matching is not exact adjustment; gene vetoes use another population | Source and analytical counterexample |
| 7 | High | Heterogeneous gene sums and estimated noise do not justify exact NB thinning | Distributional derivation |
| 8 | High | Routine segment scales omit dependence, fitting and selection uncertainty | Source and recounted coverage |
| 9 | High | Frozen-cell confirmation supplies nominal RNA diagnostics, not genomic confirmation | Source and role contracts |
| 10 | High | Allele evidence remains conditional on expression state, noise and measurement integrity | Source and retained falsification |
| 11 | High | Spatial smoothing does not model within-spot mixtures or validate genotype | Source and assay mismatch |
| 12 | High | New RNA prediction gates do not establish independent predictive or DNA accuracy | New source and transfer evidence |
| 13 | High | Bulk amplitudes and boundary probabilities are working RNA quantities | Algebra and new source |
| 14 | High | Exposed benchmark gains remain conditional; broad superiority gates do not pass | Independent recounts and scorer design |
| 15 | Moderate | Boundary averaging has cubic-time steps despite quadratic-time documentation | Direct loop-count derivation |
| 16 | Moderate | Legacy names and some explanations still imply stronger scientific objects | Direct source/documentation mismatch |
| 17 | Moderate | Genomic resolution and measurement provenance need narrower interpretation | Binning and input contracts |

**What is being estimated.** The cell assignment model is approximately

`m_ckj = L_c × softmax_j(log(base_j) + offset_j + dose_kj + x_c gamma_j + U_c W_j)`

with a product of NB scores and working variance `m_ckj + phi_j m_ckj²`. Here L is the observed total over retained features, x is centred log L, k indexes fitted RNA components, and component zero has dose fixed to zero. Reference composition, query-specific reference weights, offsets, depth slopes, programmes and dispersion all influence the interpretation of a dose. See [model initialization](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:703), [log rates](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:936) and [assignment likelihood](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1206).

Pooled clone segmentation instead uses a four-state Poisson-lognormal HMM with log-intensity location `log(E_j) + centre + beta log(C_j/2)`, where the working grid is C = 1, 2, 3, 4. The PLN unconditional mean includes its lognormal correction. E is not automatically that unconditional mean. These are RNA count model families, appropriate modalities to investigate here; the questions are whether their aggregation, dependence, nuisance and selection assumptions fit this application. See [profile.py](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/profile.py:95).

The newer public profile APIs operate on prepared RNA contrasts and use Gaussian or Laplace working smoothing losses. They deliberately do not supply count-model calibration. Their output is separate from the default clone controller.

| Quantity | Supported interpretation | Interpretation requiring additional evidence |
|---|---|---|
| `log_ratio`, gene contrasts | Centred relative RNA contrasts against a chosen baseline | Absolute DNA copy number |
| `rna_log_effect` | Conditional HMM posterior mean of a working RNA state effect | Unbiased genomic dosage |
| `assignment_dose` | Profile used in cell assignment | Necessarily the same as an exported resegmented profile |
| `clone_posterior`, `p_normal` | Plug-in responsibilities over fitted RNA components | Calibrated genomic membership or normal-cell probability |
| `n_eff` in clone summaries | Sum of component responsibilities | Independent biological n or universally a Kish effective sample size |
| New stacked/posterior-region profiles | Descriptive regional RNA estimates | Genomic states, uncertainty intervals or clone assignments |
| `tree.tsv` | Similarity dendrogram | Tumour evolutionary phylogeny |

**1. Blocking: RNA dosage and genomic identity remain different estimands.**

Softmax is invariant to adding a constant to every bin. Balanced genome-wide scaling therefore disappears from this compositional model. A regional transcriptional effect can also produce the same RNA intensity vector as a DNA dosage effect. No classifier can distinguish mechanisms with identical observable distributions without additional information or restrictions.

A second ambiguity involves the shared offset. When the normal-like component has zero population weight, adding a regional vector h to the offset and subtracting h from all tumour doses leaves the tumour means unchanged. Keeping a small numerical prior weight on component zero does not produce observed normal cells or validate a normal baseline. The convention `dose_0 = 0` fixes a parameterisation; it does not establish biological diploidy.

Pairwise contrasts can remove shared reference terms under a common nuisance model. They do not remove clone-specific transcription, lineage differences or residual depth composition. The same warning applies to gene-coherent changes: coordinated expression can span many genes, while DNA dosage compensation means not every gene responds equally to an actual CNA.

The appropriate claim is relative RNA dosage and candidate RNA populations. Genomic identification requires independently justified assumptions or corresponding DNA evidence at the claimed resolution. Keep the package general-purpose; a project-specific hypothesis about one alteration cannot supply those assumptions for other datasets.

**2. Blocking: local support rules are not specimen-wide false-split control.**

The public default remains even/odd-bin splitting, eight refinement rounds, three search cycles and no fitted-parent simulation gate. The search repeatedly chooses candidate groups and intervals, updates memberships and nuisance parameters, and merges or retains groups. A posterior threshold within that process is not a frequentist error bound over the whole procedure. See [defaults](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:85), [merge logic](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1066), [split selection](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1527) and [final search](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1725).

In thinning mode, a merge may assign cells using the proposal fold, but its profiles and nuisance fits have already used full counts. That local reassignment does not undo selection. The optional fitted-parent null is narrower still: it conditions on selected cells and fitted nuisance quantities. The public pipeline correctly rejects its combination with gene-aware support or positive allele evidence when the simulator cannot reproduce those sources. Preserve those restrictions; see [pipeline.py](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:119).

The latest completed production support regression is encouraging but small and exposed. I aggregated all 30 case records and independently maximized their saved confusion-table assignments by exhaustive rectangular matching:

| Fixed simulated condition | Before | After final-support correction | Correct assignments before → after |
|---|---:|---:|---:|
| No added regional-state null, additional groups | 0/10 | 0/10 | 5,994/6,000 → 5,994/6,000 |
| Continuous-state null, additional groups | 0/10 | 0/10 | 5,998/6,000 → 5,998/6,000 |
| Two populations, exact two-profile chr12 recovery | 9/10 | 10/10 | 4,249/6,000 → 4,271/6,000 |

The denominator includes every simulated query cell, including normal-like/unassigned errors absent from the tumour confusion table. These are not 6,000 independent specimens. The evidence is the [complete paired report](/gpfs/commons/home/jameslee/Haplodose/dev/final_membership_support_20261110/production_regression_summary_v1/report.json). This recount checks stored outcomes and assignment tables, not regenerated counts or independently reconstructed raw cell labels.

Zero errors in ten independent null draws gives a one-sided 95% binomial upper bound of 25.8866%, before accounting for multiple null strata. Thus 0/10 cannot establish a 5% target. The separately frozen fresh protocol uses 200 seeds per condition and three jointly protected endpoints. Independently inverting the binomial CDF at alpha = 1/60 gives upper bounds 2.02636%, 4.59656% and 5.33450% for 0, 3 and 4 errors respectively. Those arithmetic requirements are coherent. No final fresh-study report existed when checked; partial outcomes were not scored in this review. See the [frozen protocol and measurement qualification](PRODUCTION_SUPPORT_FRESH_2026-11-10.md).

Even a passing fresh experiment would establish performance only under its declared generators. The legacy “matched” gene-count generator is not necessarily an exact bin-level NB law. It cannot isolate pure selection error from model misspecification or establish biological calibration. Define the specimen-level error event, biological null, power endpoint and entire allowed procedure before independent evaluation.

**3. Blocking: the strongest new profile result is pooled, not cell-specific.**

The 181-line profile regression pools query cells and compares regional RNA profiles with relative cell-line DNA. The replication units are query lines, with shared reference-group dependence. Cells, genes, windows and repeated reference fits do not multiply the biological sample size. Cross-culture bulk DNA is also not same-cell truth.

The distinction is empirically relevant. I recomputed means from the saved [individual-cell score table](/gpfs/commons/home/jameslee/Haplodose/dev/single_cell_profile_transfer_20261110/panel_v1/scores.tsv). For three selected cells from each of two exposed lines, partition averaging changes mean DNA Pearson correlation by +0.00282910 in `2313287_STOMACH` and −0.00622930 in `42MGBA_CENTRAL_NERVOUS_SYSTEM`. The latter fails the prespecified both-line transfer gate. All six candidates were selected by their RNA fitting criterion. This is a negative transfer result against bulk-line DNA, not proof that the declining cell estimates are genomically wrong. Actual cell-specific CNA, transcriptional variation and reconstructed-count effects remain alternatives.

Before a genomic membership claim, use corresponding-cell DNA or a justified independent assignment standard, explicit DNA callability, uncertain truth labels, and a specimen-level evaluation. Report assignment accuracy, false distinctions, abstention and clone abundance together. Pooled correlation cannot substitute for these endpoints. The broader primary [single-cell CNV benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC12491403/) likewise treats profile recovery, normal-cell behavior and subclonal structure as different tasks; its findings are context, not validation of HaploDose.

**4. High: the controller is a hybrid algorithm, not demonstrated EM for one objective.**

Cell assignments use compositional NB scores. Pooled profiles use PLN segmentation, state shrinkage and shared-profile restrictions. Dispersion is updated from hard-MAP residual moments. The offset uses a retained-population residual that excludes depth and programme terms. Programme fitting has another conditional objective. These steps are not established M-steps for a single joint likelihood. See [shared profiles](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:838), [refinement](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1240), [offset](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1347) and [dispersion](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1363).

The current code correctly states that stable MAP labels or a fixed iteration budget does not certify convergence. The saved [transition investigation](PROFILE_TRANSITION_2026-11-09.md) records decreases in the conditional NB objective. Such decreases disprove monotonicity of that objective; they do not automatically make a regularised HMM step erroneous. A paper should describe the actual hybrid procedure and its tradeoffs rather than imply a common EM guarantee.

The newer combined research controller also illustrates why independently helpful interventions need combined evaluation. Its [saved four-case reader](/gpfs/commons/home/jameslee/Haplodose/dev/final_membership_support_20261110/saved_combined_v1/report.json) retains a one-group failure, a three-group failure and a state-null false split. The case-four trace shows a refinement preview with two supported groups rejected by a constraint requiring an increase in group count, restoring three groups. This is a research-path interaction, not evidence that the simpler corrected production path has those same outcomes. See [combined evidence](COMBINED_RECOVERY_REGRESSION_2026-11-10.md).

The final-support correction is nevertheless well motivated: it enforces the existing mass threshold on the returned responsibilities rather than the previous iteration. It does not establish optimal model complexity or genomic support. Future controller changes need one explicit selection rationale plus joint recovery/error evaluation, not a collection of case-specific rescue rules.

**5. High: reference adaptation and nuisance fitting do not identify the cause of an RNA shift.**

Query counts select the mixture of reference compositions. This can absorb biological signal as a reference change. For a two-region reference A = (0.5, 0.5), multiplying region one by 1.5 and renormalising gives (0.6, 0.4). That is exactly indistinguishable from choosing an unperturbed reference B = (0.6, 0.4). Accurate optimisation of the reference mixture cannot resolve this ambiguity. See [mixture fitting](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:462).

High-pass programme construction removes whole-arm means from the programme basis. It therefore preserves whole-arm CNA-shaped signal, but also preserves whole-arm transcriptional state or reference mismatch. “Protected from nuisance removal” is a design property, not demonstrated genomic specificity. Conversely, focal dosage signals can overlap retained programme directions and be attenuated. See [programme construction](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/programmes.py:61).

Depth slopes ordinarily come from the query and depend on current grouping. They may absorb genuine abundance-associated biology or fail to remove nonlinear state effects. Fixed reference slopes require transportability across cell type, assay and depth range. Reference sensitivity should therefore vary donor/line, lineage, library chemistry and fit composition while preserving a fixed evaluation estimand. Counts pooled across references weight high-RNA contributors more heavily; this is not automatically donor-balanced reference construction.

**6. High: coarse depth matching does not guarantee cancellation, and the gene veto evaluates another population.**

The [matching docstring](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:955) says depth-driven composition cancels regardless of functional form. The implementation balances mass in eight coarsened log-depth strata. Equal stratum frequencies do not imply equal depth distributions within strata.

A direct counterexample uses group A log depths `0, 0.2, ..., 1.4` and group B log depths `0.1, 0.3, ..., 1.5`. Each combined quantile stratum contains one observation from each group, so matching retains equal weights everywhere. Mean log depths remain 0.7 versus 0.8. Replicating each value increases support without changing the imbalance. Even a linear depth effect need not cancel; nonlinear effects need not either. This is an analytical counterexample to the guarantee, not a newly run HaploDose simulation. Report residual depth balance and restrict the claim to coarse adjustment.

There is a second mismatch. Regional contrasts use depth-matched weighted counts, while [gene_consistent](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1037) sums unweighted gene counts over a proposal subsample with hard labels. It can therefore evaluate another population. It also returns permissive support when too few cells or informative genes are available. A passing support flag may mean that the veto was uninformative, rather than that concordance was demonstrated.

The sign test treats concordant gene signs as binomial evidence after interval and direction selection. Correlated genes and transcriptional programmes do not meet the simple independent-sign null. The source explanation that CNA shifts every gene while cell state shifts only a minority is not generally valid. The earlier [gene-pool investigation](PRODUCTION_GENE_POOLS_2026-11-09.md) found all observed vetoes passing in its traced cases, so this concern must not be presented as a proven cause of their recovery failures. Preserve separate flags for support, insufficient evidence and model-based veto decisions.

**7. High: aggregation and thinning require the actual count law, not moment agreement alone.**

For independent gene counts with means m_g and NB dispersions phi_g, the bin sum has

`E[Y_bin] = sum_g m_g`, and `Var[Y_bin] = sum_g m_g + sum_g phi_g m_g²`.

A moment-matched bin dispersion reproduces those two moments but does not make the sum NB. Independent NB variables are closed under addition under a common success probability, equivalently a common Gamma scale in the Poisson-Gamma representation. Heterogeneous gene means and dispersions generally violate that condition. Shared cell-state factors add covariance terms as well.

The beta-binomial split in [thin_counts](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:169) can provide independent parts under the requisite NB law and correctly supplied size. A plug-in size learned from the same data, an arbitrary moment-matched bin, or adaptive reuse of both parts does not inherit that theorem automatically. The primary [data-thinning paper](https://www.jmlr.org/papers/volume25/23-0446/23-0446.pdf) supplies a distributional construction, not a blanket guarantee for every fitted pipeline using its randomisation.

Likewise, multiplying marginal NB terms while fixing L to the observed sum of those same bins produces a working/composite score. It is not automatically the exact joint distribution conditional on L. Some dependence is introduced by compositional normalisation even before biological correlation. The current confirmation docstring acknowledges the composite-likelihood interpretation. Distinguish count conservation, moment matching, fold independence and valid post-selection inference as separate properties.

**8. High: the routine segment diagnostic is not a calibrated confidence interval.**

The [segment exporter](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/profile.py:152) computes an independent-bin scale proportional to `sqrt(sum_j(tau_j² + v_j + 1/(Y_j+0.5))) / n_bins`. It omits cross-bin covariance, reference and nuisance fitting, membership estimation, and selected interval boundaries. Its centre is an observed pseudocount log-ratio statistic, distinct from the HMM posterior-mean state effect. The current metadata appropriately marks the quantity unvalidated and does not assign it a validated parameter-interval estimand.

I independently recomputed coverage directly from `lower <= true_mu_I <= upper` in the preserved [interval table](/gpfs/commons/home/jameslee/Haplodose/dev/recovery_20261107/uncertainty_simultaneous_v1/intervals.tsv), rather than trusting the stored coverage flag:

| Saved generator | Legacy fixed-interval coverage | Legacy selected-interval coverage |
|---|---:|---:|
| Matched homogeneous | 934/1,000 | 818/1,000 |
| Matched heterogeneous | 985/1,000 | 913/1,000 |
| Correlated heterogeneous | 703/1,000 | 54/1,000 |

These are conditional simulation results for the preserved diagnostic construction, not coverage estimates on biological cohorts. They nevertheless rule out casually reading these scales as universally nominal 95% intervals.

The standalone [systematic interval helper](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/uncertainty.py:12) covered all 1,000 intervals in each of the same six saved settings. Its contract assumes independently fixed expected intensities, known Gaussian latent covariance and a frozen genomic universe. It does not run as the default pipeline uncertainty calculation. Model-conditional simultaneous coverage does not establish practical width, power or validity after estimating those quantities from the query. A future uncertainty claim must name the parameter, conditioning set, selection family, dependence model and biological replication unit.

**9. High: confirmation isolates roles but does not resolve shared state or calibrate nominal p-values.**

[FrozenDiscovery and score_confirmation](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/confirmation.py:128) have useful safeguards: cell identities cannot overlap discovery, candidates are fixed, scoring uses proposal features, and exact candidate-score ties abstain. Test contrasts are then regional-versus-background RNA differences. Their delta-method variance omits cross-feature covariance and their one-sided p-values are explicitly nominal and unadjusted.

Disjoint cells from the same specimen share biological state, batch and reference mechanisms. Disjoint features can also be correlated. Confirmation of the same RNA pattern therefore does not decide whether it is genomic. Candidate-specific programme fits on proposal features do not automatically remove the corresponding test-feature state signal. Exact-tie abstention is computationally appropriate, but near ties still have no calibrated genomic confidence threshold.

If these diagnostics are later used for formal selection, prespecify the family of contrasts and the entire decision rule before calibration. Freezing one RNA-derived candidate list is valuable process isolation; the scientific null still needs to describe the state alternatives that should not count as distinct genomic clones.

**10. High: allele information helps only under a defensible measurement and biological model.**

Pooled alleles use a beta-binomial HMM with reference-derived overdispersion and tumour-fitted parameters. Cell scoring uses a binomial default or an explicitly supplied cell-level concentration. Pool-level and cell-level concentration parameters are different quantities; using the former for the latter would not be justified by a shared parameter name. See [filtering](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:166), [cell likelihood](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:323) and [pool fitting](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:428).

Reference SNP filtering is itself selected on observed imbalance before dispersion estimation. Allele-specific expression, imprinting, mapping bias, phase error and cell-state-specific allelic effects remain alternative explanations. The preserved [shared-state experiment](ALLELE_SHARED_STATE_2026-11-08.md) directly shows why reserved SNP evidence alone cannot identify a genomic split when both discovery and reserved observations share the same allele-expression state.

The [allele reader](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:79) aligns cells and sites and accumulates counts, but its table schema cannot by itself establish molecule deduplication, consistent REF/ALT identity, patient-specific phasing or absence of overlapping-read dependence. Those properties belong to the upstream measurement chain. Locus-wise balance also does not establish normal identity, balanced gain absence or whole-genome ploidy.

The final allele-snapshot defect is now resolved at source level and in saved targeted tests, as documented below. That correction ensures internally consistent scores; it does not calibrate allele probabilities or make the final iterative state a fixed point.

**11. High: spatial smoothing and mixed-spot inference need separate scientific contracts.**

The [spatial posterior update](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:763) uses a directed row-normalised neighbour graph and finite iterative smoothing. The strength is chosen using neighbour-predicted cell likelihoods. This is a spatial regularisation heuristic, not an established symmetric Potts posterior sampler or energy optimiser. Positive spatial structure can reflect histology, state, technical effects or mixed composition, not only genotype.

For multi-cell spots, a responsibility over competing single-component profiles is not automatically an estimate of within-spot clone fractions. A genuine mixture likelihood combines expected contributions before scoring counts and must address RNA contribution versus cell proportion. Segmented cellular spatial assays and multi-cell capture spots therefore require different validation. Spatially coherent RNA error can look biologically plausible and improve neighbour prediction.

Cell-ID-aware graph construction is a valuable integrity safeguard. Section isolation is a separate caller responsibility: the pipeline checks patient/dataset consistency within a specimen, but [graph construction](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:302) has no separate section identifier. [Inference] Combining unregistered sections under one specimen can therefore connect unrelated locations. In addition, one nonfinite coordinate prevents graph construction for the entire fitted specimen; the resulting `spatial_graph=False` is recorded in QC, but the requested spatial analysis has changed. Validate section-specific specimen IDs or explicitly registered coordinates and make this fallback visible to users.

Evaluate genotype/state discordance, boundaries, spot admixture and section registration with corresponding truth at the intended resolution. The primary [spatial CNA benchmark](https://pmc.ncbi.nlm.nih.gov/articles/PMC13642305/) separates profile, classification and subclone-related tasks; its results do not validate transfer of HaploDose's cell model to mixed spots.

**12. High: new profile-selection gates establish conditional RNA fit, not independent prediction or DNA truth.**

The [stacking API](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/stacking.py:103) excludes each target gene from its base predictors, then fits nonnegative simplex weights to those target values. The base predictions are out of their own training subsets. The selected ensemble's reported loss is still used to fit and choose the ensemble. It is therefore not an unbiased held-out estimate of the selected method's generalisation error. The source correctly says this.

Prepared contrasts have already undergone shared normalisation, feature selection, reference fitting and centring. Within-bin rank modulo ten determines folds, so nearby genes and biological programmes remain correlated across them. The optimiser gives each eligible gene a target loss contribution; gene-rich genomic regions can contribute more than gene-poor regions. None of this is inherently inappropriate for a descriptive smoother. It must not be described as independently validating a DNA profile or selecting the biologically best resolution.

The [gene-Laplace expert](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/gene_segments.py:17) and [partition average](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/average_segments.py:17) use fitted working scales and segment penalties. The latter normalises Gibbs weights over partitions; it does not integrate parameter priors into calibrated Bayesian model evidence. Its source appropriately distinguishes these weights from boundary probabilities with empirical coverage.

The L1 primal/dual bound is a useful optimisation check. It bounds objective suboptimality conditional on the supplied prediction matrix. It does not bound statistical prediction error, fitted-weight uniqueness, or biological uncertainty. If two experts have the same masked predictions but different full-data profiles, multiple equally optimal weights can produce different exported profiles. That is an algebraic possibility, not a reproduced failure in these saved results. If stability of the exported profile becomes a claim, check it over equivalent optima or perturbations, rather than using the objective gap alone.

The failed cell-transfer example in finding 3 supplies an actual limit to interpreting an improved RNA selection loss. Preserve these modules as explicit descriptive options and freeze an independent evaluation before a general accuracy claim.

**13. High: bulk beta and posterior-region weights are not purity or calibrated DNA probabilities.**

For purity p and DNA copy state C, a simple linear mixture has relative signal `1 − p + p C/2`. HaploDose's working effect is instead `beta log(C/2)`. These forms are not equivalent for a single beta across states. With p = 0.5, exact matching would require beta = 0.41504 for C = 1, 0.55034 for C = 3 and 0.58496 for C = 4. This independently checked counterexample suffices to reject an interpretation of fitted beta as a purity estimate. Expression response, centring and reference composition further complicate the relationship. See [bulk fitting](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:74).

The median-bin mode constructs `E_bin × exp(median gene log ratio)`. This is a transformed pseudo-observation, not a count generated by summing independent reads. A count likelihood applied to it is a working approximation. The new native/posterior-region export correctly requires sum bins when allocating fitted bin exposures to genes; preserve that restriction. See [median construction](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:108) and [export validation](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:156).

The [new posterior-region estimator](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/posterior.py:301) estimates normal-reference locations/spreads, uses count-related working weights, chooses a beta from the existing grid, and averages gene summaries across constant-state runs. These are conditional smoothing probabilities. Grid selection, reference estimation, upstream correction and state misspecification remain outside a calibrated posterior interpretation.

Normal effects supplied to this API may themselves have been prepared against overlapping leave-one-out reference panels. Holding out one already prepared normal row is not equivalent to excluding that donor from every upstream transformation of the remaining rows. The minimum of three reference profiles is a computational requirement, not evidence of calibrated spread estimation or independent biological n. Arm z scores likewise remain descriptive. The new bulk export's NaN/reference-only rows avoid querying normals against their own final panel, which is a useful safeguard.

**14. High: recent benchmark improvements are real within their scope, but do not establish broad superiority.**

From all 181 unique rows in [paired_profiles.tsv](/gpfs/commons/home/jameslee/Haplodose/dev/partition_averaged_stack_20261110/summary_v1/paired_profiles.tsv), I independently obtain mean 1-Mb Pearson correlation 0.71353557 for the averaged stack versus 0.71135089 for the preceding gene-Laplace stack: mean gain 0.00218468, with 133 wins, 36 losses and 12 ties. This is a small, consistent exposed regression gain. It is not an effect on clone membership or event sensitivity.

The saved [comparison report](/gpfs/commons/home/jameslee/Haplodose/dev/partition_averaged_stack_20261110/summary_v1/report.json) uses eight reference-group resampling units and 20,000 bootstrap draws, Holm correction across eight contrasts, and basic simultaneous intervals. All eight conditional gates pass. I inspected those declarations and decisions; I did not rerun their bootstrap. Eight reference groups are not 20,000 biological replicates. The quoted interval and p-value are conditional approximate bootstrap outputs, not established frequentist coverage under repeated method development.

I also recounted all 5,792 entries in the [16-method field table](/gpfs/commons/home/jameslee/Haplodose/dev/partition_averaged_stack_20261110/field_v1/scores.tsv): 181 cases, two endpoints, and 167 cases complete across the family. CopyKAT lacks 13 populations and SCEVAN 14 for each endpoint. On the common 167 cases, the new stack's mean 1-Mb Pearson is 0.70931634 and mean arm Spearman is 0.85084330. The [fixed 30-comparison report](/gpfs/commons/home/jameslee/Haplodose/dev/partition_averaged_stack_20261110/field_v1/report.json) marks every formal comparison untestable under its declared missing-family rule. Common-complete point means cannot replace that rule.

The 181-case and 167-case means use different scoring populations and should not be combined into one performance narrative. RNA outputs locked before a particular DNA read protect that execution against direct outcome-dependent adjustment. They do not erase earlier exposure of the cohort during method development. Likewise, correlation is invariant to positive affine rescaling: strong correlation can coexist with wrong amplitude, baseline or event thresholds.

For a manuscript, freeze the estimator and endpoint family before independent specimens, report all comparator failures, and preserve simple pseudobulk/smoothing baselines. Include event direction, amplitude, callability, failure/abstention and membership endpoints where claimed. The older [Kinker scorer](/gpfs/commons/home/jameslee/Haplodose/benchmarks/kinker/evaluate.py:63) uses paired line-level Wilcoxon tests with Holm correction; newer clustered analyses do not retroactively validate independence assumptions of older scores. Shared reference uncertainty needs to match the actual contrast and resampling unit.

**15. Moderate, direct resource/documentation defect: posterior boundary averaging is not quadratic-time as documented.**

The public [posterior_rna_profile docstring](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/posterior.py:323) claims quadratic memory/time within a genomic block. Quadratic storage is a reasonable description of the dense interval arrays, but several executed loops are cubic-time in block length B:

- [_run_probabilities](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/posterior.py:199) sums all covering intervals separately for every bin.
- [_average](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/posterior.py:292) again sums every covering interval for each bin and each of three summaries.
- [_interval_summaries](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/posterior.py:274) recomputes precision-weighted numerator and denominator over all selected genes for each interval, even though the public estimator ultimately returns the weighted-median branch.

The number of interval/bin incidences is `sum_j (j+1)(B-j) = B(B+1)(B+2)/6`. At B = 1,000 this is 167,167,000 terms for one covering-interval pass, before the repeated passes and gene-level work. This is an exact loop-count derivation, not a runtime measurement. Actual time and memory at that size were not benchmarked in this review.

Correct the complexity statement before users plan large-block allocations. An eventual implementation improvement could avoid unused summaries, accumulate interval contributions with range additions, and maintain prefix weighted sums. That would require numerical and scientific-output equivalence checks. No such optimisation was implemented here. Until then, retain genomic blocks and report measured scaling across block size and genes per bin; do not infer genome-wide scalability from small-block regression success.

**16. Moderate: output labels and some explanatory text still encourage incorrect interpretation.**

The README and scientific contract are substantially more careful than some lower-level names and docstrings. Examples include `malignant_call`, `clone_posterior`, `cn_state`, and the bulk module introduction's “copy-number profile” and purity-response description. The [gene veto explanation](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1040) and [coarse matching guarantee](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:959) are also stronger than the underlying biology or mathematics.

The clone summary's [n_eff](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:670) is responsibility mass. Elsewhere, a weighted support calculation actually uses `(sum w)²/sum w²`. Those quantities serve different purposes, and neither is independent patient n. A tumour component can exceed the responsibility-mass minimum while having no full-count MAP-assigned cells, as the combined research trace demonstrates. That is not by itself a deletion criterion, but it matters when reporting clone abundance or interpreting group counts.

Keep backward-compatible fields if required, but make their estimands and calibration status unavoidable in tables and manuscript notation. Separate candidate RNA grouping, biological tumour classification and genomic lineage throughout. An average-linkage tree should not be used as evolutionary history without an evolutionary model and corresponding evidence.

**17. Moderate: genomic resolution and raw-count declarations have limited guarantees.**

[Binning](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bins.py:29) groups consecutive expressed genes, with arm/gap boundaries. Gene density, expression filtering, namespace annotation and missing reference support therefore change physical resolution. A fixed HMM transition per bin is not a fixed transition rate per base pair. A minimum number of bins or genes does not establish a universal physical detection limit. The newer smoothers also operate on supplied gene/bin/block identities rather than inferring assay resolution from coordinates.

The explicit build, coordinate-system and namespace contracts, canonical tie handling and exact count aggregation are valuable. Still, count integrality does not prove raw measurement provenance; rounded normalised counts can pass an integer check. The [descriptive count API](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/stacking.py:201) properly distinguishes raw and reconstructed inputs and validates identities before pooling. Its reconstructed option should remain confined to descriptive interpretation.

Benchmark and manuscript reports should specify retained genes, physical coverage, gap handling, count origin and reference construction. Mapping each fitted profile to common evaluation regions before averaging is necessary when the region operator is nonlinear. The [reference-ensemble helper](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/posterior.py:372) appropriately enforces aligned region identities and identical finite support; an arbitrary reference label still cannot prove biological independence or matching preprocessing.

**Defects resolved since the preceding review.** These should not remain listed as open source defects.

| Earlier issue | Current behavior and evidence | Remaining limit |
|---|---|---|
| Ribosomal regex applied to Ensembl/custom IDs | [bins.py](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bins.py:56) now requires explicit `gene_symbol` annotation for non-symbol category filtering, preserves IDs, and rejects missing required annotations. Saved old-wheel fixtures show two intended failures; candidate checks pass. See [namespace report](RIBOSOMAL_NAMESPACE_2026-11-10.md). | Annotation correctness remains external. The change can alter selected genes or reject formerly accepted inputs; historical analyses need their original source. |
| Allele parameters refitted after returned LL/R | [em](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1240) now freezes the last refinement-time allele state through final scoring. Saved old-wheel fixtures show three intended failures and the candidate passes them. See [snapshot report](ALLELE_FINAL_SNAPSHOT_2026-11-10.md). | Consistency is not convergence or allele calibration; fit pools need not equal pools from final responsibilities. |
| Final responsibilities could fall below the existing support minimum | [final support checks](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1299) remove unsupported tumour components, align associated arrays and recompute responsibilities/pools without another nuisance fit. The 30-case exposed production comparison records a recovery benefit. | This is threshold enforcement, not an optimal-complexity or false-split guarantee. |

**Strengths worth preserving.** The package separates fit and held-out references, requires explicit measurement declarations, records content-addressed provenance and environments, guards fresh output destinations, and maintains identity alignment. The newer APIs state their descriptive scope, preserve defaults, and expose numerical diagnostics. Exact native count checks, weighted pooling corrections, independent small-model comparisons and saved negative results are useful scientific infrastructure. The new profile delivery reproduces installed outputs and source-count examples, rather than relying only on research-script outputs. These strengths make remaining limitations easier to investigate.

**Recommended next scientific steps, in order.**

1. Freeze a narrow contribution: pooled regional RNA estimation, candidate RNA grouping, or genomic membership. Assign distinct validation endpoints to each retained claim. Use notation that separates counts, relative RNA effects, DNA states and latent memberships.
2. Finish and report the already frozen fresh production study at its declared final horizon. Preserve its generator qualification, both error strata and power interval. Do not use a favourable conditional result to imply a biological error guarantee or alter the method mid-study.
3. Evaluate the frozen profile estimator on independent biological units with matching DNA and suitable reference variation. Keep raw versus reconstructed inputs, pooled versus individual-cell endpoints, and common versus missing comparator populations explicit.
4. For genomic clone claims, obtain truth at the membership resolution and a credible biological single-genotype null that includes transcriptional and spatial state. Model truth uncertainty and correlated samples; do not count cells or adjacent sections as independent patients.
5. Investigate hybrid-controller failures using a declared common criterion and joint error/power assessment. Preserve simpler production baselines and rejected research interactions. Report support mass, MAP counts and uncertainty separately.
6. Correct the coarse-matching and complexity explanations, then address computational scaling only with output-preservation checks. Define any future interval parameter and simultaneous family before calibration. None of these recommendations authorises new implementation or submissions within this review task.

**Verification actually performed.**

- Recorded the initial dirty working tree and package/native source hashes; compared all 28 current source files with the saved installed-validation source snapshot, with exact agreement.
- Parsed all 26 Python package modules with `ast.parse`, without importing HaploDose or invoking inference.
- Recounted the complete 30-case production support report; independently verified saved confusion-table assignment maxima and used all-cell denominators.
- Recomputed 181-case paired profile means/wins/ties, the complete-family population and missing counts in all 5,792 field score records, and the two individual-cell transfer mean differences.
- Recomputed the 12 coverage cells from saved interval endpoints and truth values. Derived the NB aggregation, depth-matching counterexample, purity mismatch, binomial limits and cubic interval-count formula independently of the package implementation.
- Inspected saved validation logs: 508 Python tests passed, one GPU-only test skipped, ten warnings; four Rust tests passed. All 28 current package/native source files match that test snapshot. These are historical executed checks, not tests run by this review and not biological calibration.
- Checked report links, source-line targets and whitespace after writing; rechecked reviewed package/native source identities for concurrent changes.

**Not run:** current-environment `python -m pytest -q`; `cargo test --release --manifest-path crates/core/Cargo.toml`; GPU tests; new simulations; full benchmark reruns; raw-count or DNA reprocessing; runtime/memory scaling experiments. Executable checks belong on an authorised compute allocation. No active fresh-study partial outcomes were scored, and no result readiness was changed.

**Source identity.** The full reviewed implementation consists of 26 Python files and two Rust files under `python/haplodose` and `crates/core/src`, excluding generated build output. Selected SHA-256 identities for the load-bearing files follow; they identify the dirty working tree more precisely than the unchanged package version or Git HEAD.

| File | SHA-256 |
|---|---|
| `python/haplodose/clones.py` | `2bad993a565af2e3ed8950081485405ca60c443d5602786f68dd1e46abddbb13` |
| `python/haplodose/pipeline.py` | `7f57ac4c63d7b2ca37b58a823464eedfd3d9b596487bf0c7aa11710eb69fb7c5` |
| `python/haplodose/profile.py` | `7db040cbb9503bb6f48b64465af23fc39eb7190f8e79c8c02f539dd7f5cc3612` |
| `python/haplodose/alleles.py` | `112109a2066158b069c4e3fb87ca5417ff0a1de7164625d6a3211a609deb958e` |
| `python/haplodose/confirmation.py` | `d7dcc9bd9ce09b7e65c45a7adc0f462ee18a64e154caf2df3b840ea8d274b1ac` |
| `python/haplodose/bins.py` | `3b604339287ecf14a13cce84cca5ba81be55fb1a092b89fbd42ad662950d188d` |
| `python/haplodose/bulk.py` | `623b2fbb3ed90855f47ede5d1ba5fc7f0702d500412a268e9b82f5a511d333fd` |
| `python/haplodose/stacking.py` | `32566f396c0e9c5d39af2add8fdad42dafa7e0fdbc5206db6fd7f40f8c7870e1` |
| `python/haplodose/gene_segments.py` | `68fbddfed2829a14bb9de8608b23bf381548b5cb769a33533c956905c58b46e0` |
| `python/haplodose/average_segments.py` | `b0078b0a6f880f1a10e4451bcde4f7b2509651186bcc33b7472b46cadd18f8d7` |
| `python/haplodose/posterior.py` | `9440f0960b792d03491114d2df108f7f440e0246ef6ab10ec9a30730de663946` |
