# HaploDose: deep scientific code review

Review date: 2026-11-09. Scope: the current working tree, including existing uncommitted package changes. Git HEAD: `58157f2a3b4e76228c444664c36448571ebc6578`; package version: `0.2.0`. HEAD alone does not identify the reviewed implementation. Selected source hashes are recorded below.

**Assessment.** HaploDose has a useful architecture for exploratory regional RNA profiling and candidate RNA grouping, substantial numerical safeguards, and an unusually explicit record of negative results. The present evidence does not support interpreting its groups as validated genomic clones, its responsibilities as calibrated genomic membership probabilities, or its segment diagnostics as confidence intervals. The principal obstacles are identification and statistical calibration, with additional inconsistencies between the models used at different fitting stages. Faster optimization or a larger test suite cannot by themselves resolve those obstacles.

The current README and scientific contract acknowledge much of this. This review does not reclassify acknowledged limitations as newly discovered software bugs. It distinguishes structural limitations, observed failures in specific saved experiments, and narrower source-level defects. It also identifies concrete issues that deserve focused checks: namespace-dependent ribosomal exclusion, final allele-state/likelihood freshness, and differences between the populations used by regional tests and gene vetoes.

No implementation, defaults, readiness records, historical results, or claim decisions were changed. No jobs were submitted. The sole repository artifact added by this review is this document.

**Evidence and review boundaries.** I inspected the public pipeline, count/identity handling, binning, reference construction, programme learning, clone controller, pooled and cell likelihoods, HMM kernels, allele model, confirmation module, bulk mode, uncertainty helper, native aggregation, provenance, representative tests, and selected benchmark generators/scorers. I examined the scientific contract, historical validation protocol, recent mechanism reports, and selected saved JSON/TSV evidence. This was not a line-by-line audit of every development snapshot or every benchmark script. Raw sequencing, large count matrices, and corresponding DNA were not reprocessed.

The initial working tree contained 23 modified tracked paths and 614 untracked entries. These were preserved. Concurrent development added other reports, a script, and scheduler logs during the review; they were not created or edited by this review. Package/native source hashes remained stable. Evidence inspection closed at approximately 16:50 America/New_York. This review uses static inspection and bounded reading/arithmetic on existing evidence, not new scientific inference runs. Historical test logs are identified as historical. No current-wheel Python/Rust tests or CPU/GPU inference were run.

**Priority map.** “Blocking” below means blocking the indicated scientific claim, not blocking exploratory execution. “High” indicates a material modeling or inferential limitation. “Moderate” indicates a bounded correctness, interface, or reporting problem.

| ID | Priority | Finding | Evidence type |
|---|---|---|---|
| 1 | Blocking | RNA observations do not uniquely identify DNA dosage or genomic clones | Algebra and source |
| 2 | Blocking | Adaptive clone search has no established specimen-wide error control | Source and recounted saved experiments |
| 3 | Blocking | Independent biological validation at the claimed resolution is still missing | Validation contracts and cohort roles |
| 4 | High | The clone controller combines different objectives and does not certify a common optimum | Source and saved transition diagnostics |
| 5 | High | Reference fitting, depth correction, and programmes can absorb or preserve the wrong biological signal | Source and explicit counterexamples |
| 6 | High | Gene aggregation and estimated noise do not establish independent NB thinning | Distributional derivation and source |
| 7 | High | Regional matching and the gene veto can test different cell populations | Source; impact on current cases is limited by negative evidence |
| 8 | High | Segment scales omit dependence, fitting, and selection uncertainty | Source and independently recounted coverage |
| 9 | High | Frozen-cell confirmation is useful isolation, but not calibrated genomic confirmation | Source |
| 10 | High | Allele evidence remains conditional on noise, phase, and shared expression state | Source, primary literature, existing falsifications |
| 11 | High | Spatial smoothing does not model mixed spots or establish genotype confidence | Source |
| 12 | High | Bulk state amplitudes, purity, and pseudocounts have different statistical meanings | Algebra, source, recounted comparison decisions |
| 13 | Moderate | Ribosomal exclusion silently changes with the feature namespace | Direct source defect |
| 14 | Moderate | Final allele parameters can be newer than returned likelihoods/responsibilities | Direct call-order concern; numerical impact unmeasured |
| 15 | Moderate | Several documentation and output names still overstate the scientific object | Direct documentation/source mismatch |

**What the implementation estimates.** The measurement model is for RNA counts. Negative-binomial cell models and Poisson-lognormal pooled working models are appropriate families to investigate for RNA; no methylation or other assay-specific variance method is being silently transplanted here. Their appropriateness for a particular aggregation, nuisance fit, or selection procedure still requires separate evidence.

For cell or spatial row c, bin j, and component k, the expression assignment model has the form

`m_ckj = L_c * softmax_j(log(base_j) + offset_j + dose_kj + x_c*gamma_j + U_c W_j)`.

It scores bins using NB variance `m + phi_j*m^2`. L is the observed total over retained genes; x is centered log L. Component zero has assignment dose fixed to zero. References determine the feature universe and reference compositions; query data choose reference mixture weights, and ordinarily update depth slopes, dispersion, offsets, memberships, and profiles. See [clones.py:703](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:703), [clones.py:1206](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1206), and [estep.py](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/kernels/estep.py).

Pooled segmentation instead uses a four-state PLN HMM with log-intensity location

`log(E_kj) + centre_k + beta*log(C_j/2)`, with `C_j in {1,2,3,4}`,

and latent variance `tau_k^2 + reference_variance_j`. Here E is a baseline exposure before the state effect, not the unconditional mean of the PLN: the latter includes the lognormal mean correction. The states are a working grid of relative RNA effects. They do not independently measure DNA C. See [profile.py:95](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/profile.py:95).

| Export or quantity | Defensible interpretation | Interpretation not established |
|---|---|---|
| `log_ratio` | Centered, pseudocount-adjusted bin RNA contrast against a fitted baseline | Absolute DNA copy number |
| `log_ratio_gene` | Median gene contrast within each bin, with separate centering | The same estimand as a summed-count contrast |
| `rna_log_effect` / `log_dose` | Conditional HMM posterior mean RNA state effect | An unbiased DNA dosage estimate |
| `assignment_dose` | Profile actually used in the expression assignment likelihood | Necessarily the same profile as the separately resegmented export |
| `clone_posterior` | Plug-in responsibility over fitted RNA components | Validated probability of genomic membership |
| `p_normal` | Responsibility of the reference-like component | Probability that a cell is biologically normal |
| `resolved` | Passed specified working-model support rules | Controlled specimen-level false-positive probability |
| `tree.tsv` | Average-linkage similarity dendrogram | Tumour evolutionary history |

**1. Blocking: DNA state and clone identity are not identified by the RNA model.**

The softmax model is invariant to adding a constant across all bins. Balanced genome-wide scaling therefore disappears. More generally, a regional transcriptional change can give exactly the same RNA intensity vector as a regional dosage change. A high likelihood or a precise RNA contrast cannot distinguish two mechanisms with the same observation distribution.

There is also a shared-offset ambiguity. If the reference-like component has zero population weight, adding any regional vector h to `offset` and subtracting h from every tumour `dose` leaves every tumour mean unchanged. The numerical floor on the reference-like mixing weight does not establish an independently observed normal anchor. A fixed zero dose for component zero is a parameter convention, not biological truth. Relevant code is [clones.py:936](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:936) and [clones.py:1323](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1323).

Pairwise RNA differences can cancel a shared reference term under the assumed common nuisance model. They do not cancel clone-specific transcription, cell-type mixture, or all depth-related composition differences. The claim that a pairwise contrast is reference-independent is consequently conditional; it is not a route to genomic identification without further assumptions.

**Required response:** retain RNA-based terminology, define the genomic resolution of any proposed validation, and use independently constructed DNA truth for genomic claims. Report observationally equivalent examples as limits, not as optimization failures. For downstream associations, the replication unit must be donor/specimen or independently derived line, with repeated samples nested. Cells or spots assigned by this model cannot become independent patient replicates.

**2. Blocking: the adaptive search does not control the probability of an unsupported extra distinction.**

The public default is `split_mode='evenodd'`, with no split-null filter. Proposal and difference features are separated, but initial nuisance fits and later membership/profile updates use full counts. The same test features contribute repeatedly to agglomeration, splitting, merging, and final reporting. In thinning mode, proposal-only reassignment before merging still uses profiles and nuisance parameters learned from full counts. Selection is not removed by that local reassignment. See [clones.py:85](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:85), [clones.py:1066](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1066), [clones.py:1503](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1503), and [clones.py:1722](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1722).

A 0.95 HMM posterior cutoff is not a 5% familywise error guarantee over proposed groups, intervals, pairs, iterations, specimens, and tuning choices. The optional fitted-parent simulation filter does not establish such a guarantee either. It conditions on selected parent cells and fitted nuisance quantities. The public pipeline correctly rejects its use with the gene-aware workflow because the bin-only null cannot reproduce the gene veto; it also rejects positive allele evidence with that expression-only null. Those rejection checks should be preserved. See [pipeline.py:119](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:119).

Saved full-search research experiments make the practical concern concrete. I recounted all 600 evaluation records in [evaluation.jsonl](/gpfs/commons/home/jameslee/Haplodose/dev/recovery_20261107/stage2_global_rank_v3/evaluation.jsonl):

| Frozen research procedure | Matched NB single-population null | Correlated-state single-population null | Two-population detections |
|---|---:|---:|---:|
| Raw groups | 0/200 | 184/200 | 200/200 |
| NB-calibrated global gate | 0/200 | 184/200 | 200/200 |
| State-calibrated global gate | 0/200 | 8/200 | 65/200 |

These are binned, known-reference simulation experiments, not error rates for today's complete public raw-data workflow and not biological replication. The state-calibrated gate's 4% point estimate has an 8.9639% one-sided upper bound in the predeclared nine-member Bonferroni family. I independently reproduced that bound by binomial-tail inversion. Its detection rate is 32.5%, and the unchanged raw two-population memberships have median ARI 0.07642. Rejecting a global null and correctly assigning cells are different endpoints.

The more recent exposed gene-workflow research replay after exact-profile merging has 12/100 state-null errors in each weight-unit arm, versus 41/100 in its saved predecessor. I recounted the 600 records across the twelve `panel_*_v1/results.json` files under [the exact-merge experiment](/gpfs/commons/home/jameslee/Haplodose/dev/exact_profile_merge_20261109). The six-member simultaneous upper bound is 21.8909%, independently reproduced here. This improvement is real within that replay but does not establish 5% control, and it is not a current production error-rate estimate.

**Required response:** define the entire procedure and its specimen-level error event before calibration; include the complete allowed search and preprocessing; report power and membership accuracy alongside error. Any fitted-state null needs an explicit applicability argument. More draws from the wrong null cannot repair misspecification.

**3. Blocking: available validation does not establish biological clone accuracy at the required resolution.**

The historical held-out protocol is explicit about frozen source, comparator families, missing output, common genes, and shared references. It also records negative conclusions. Its older MIX-seq results cannot validate a subsequently revised method. The current [scientific contract](/gpfs/commons/home/jameslee/Haplodose/docs/SCIENTIFIC_CONTRACT.md:225) correctly treats previously exposed cohorts as regression evidence and states that a sufficient unopened biological clone panel is unavailable.

Important distinctions remain:

- Bulk or line-level RNA/DNA correlation validates an aggregate profile endpoint, not cell-level assignment.
- Cell-line identity is not proof that all cells share one CNA state. Cross-culture DepMap profiles also need not equal the sampled culture's DNA state.
- Synthetic implants establish recovery of an imposed RNA perturbation under that simulator. They do not independently show genomic specificity or realistic allele, transcription, and purity behavior.
- Repeated seeds, libraries, sections, subsamples, and cells do not increase the number of independent patients. Shared reference pools induce dependence between otherwise distinct query lines.

The current Kinker scorer explicitly labels its revision retrospective and requires complete profile support. However, it still performs line-level Wilcoxon comparisons, while reference-pool dependence needs a separate analysis. Its truth is a contrast of gene-level DepMap values against a UMI-weighted reference, and its HaploDose endpoint is a cell-weighted mean of clone log effects. Neither should silently be called a pooled-count effect. See [evaluate.py:38](/gpfs/commons/home/jameslee/Haplodose/benchmarks/kinker/evaluate.py:38), [evaluate.py:63](/gpfs/commons/home/jameslee/Haplodose/benchmarks/kinker/evaluate.py:63), and [evaluate.py:127](/gpfs/commons/home/jameslee/Haplodose/benchmarks/kinker/evaluate.py:127).

For a single prespecified null population and zero errors, the one-sided 95% binomial upper bound is `1 - 0.05^(1/n)`. It first falls below 5% at n=59: 5.0339% at 58, 4.9508% at 59. I independently checked this arithmetic. It assumes independent, correctly established null biological specimens; multiple strata/families require a larger sample or a different declared design. It is not a general prescription to collect 59 arbitrary specimens.

**Required response:** freeze a narrower estimand, method, comparator setup, scoring universe, failure policy, and multiplicity family before the next independent evaluation. For genomic membership, corresponding-cell or suitably justified independent DNA labels are essential. Include DNA callability and uncertain truth rather than treating absent DNA events as confirmed negatives.

**4. High: the clone controller is a hybrid refinement procedure, not established EM for one joint model.**

Cell assignments use a compositional NB score. Profile updates use pooled PLN segmentation, median centering, discrete state shrinkage, and shared-profile restrictions. Dispersion is updated by hard-MAP residual moments; the offset uses a retained-population pooled residual that omits depth and programme terms. Programme scores optimize a separate conditional responsibility-weighted NB objective. These components can each be useful without jointly constituting an EM M-step. See [clones.py:838](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:838), [clones.py:1240](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1240), and [clones.py:1339](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1339).

Stopping when fewer than 0.5% of MAP labels change, or when the iteration budget is exhausted, does not imply stable responsibilities or a stationary objective. The code now records that convergence is not certified. That is appropriate.

The saved applied-transition investigation records NB objective decreases in 31 of 67 transitions. This is evidence against describing those transitions as monotone updates of that unpenalized conditional NB objective, not proof that every decrease is a software bug: HMM regularization and sharing optimize different criteria. The later [production profile trace](/gpfs/commons/home/jameslee/Haplodose/docs/PRODUCTION_PROFILE_TRACE_2026-11-09.md) also finds distinct failure locations across the three failed recovery cases. It does not support a universal claim that sharing destroys the true contrast.

The preserved five-case production thinning-mode recovery table contains two passes and three assertion failures, with no execution errors; I recounted it directly. This is a fixed development regression panel, not a 40% biological sensitivity estimate or a measurement of all default even/odd behavior. The reference-estimated gene research candidate reports 3/5 recovery, versus production 2/5; its pilot does not establish error control. See [recovery table](/gpfs/commons/home/jameslee/Haplodose/dev/recovery_20261107/recovery_panel_v1/summary.tsv) and [reference-noise report](/gpfs/commons/home/jameslee/Haplodose/docs/REFERENCE_GENE_NOISE_2026-11-09.md).

Two reports added during this review extend that development sequence: a depth-stratified candidate reports 4/5 original recovery but 13/100 state-null false splits versus 1/100 for its production comparator, with no power gain; a subsequent continuous-state veto reports 0/5 original recovery. Both reports reject promotion. I inspected their stated results but did not independently recount these later experiments. They reinforce the need to evaluate recovery and false distinctions together; they must not be pooled with the different historical null experiments above. See [depth-stratified noise](/gpfs/commons/home/jameslee/Haplodose/docs/DEPTH_STRATIFIED_GENE_NOISE_2026-11-09.md) and [continuous-state safeguard](/gpfs/commons/home/jameslee/Haplodose/docs/CONTINUOUS_STATE_GUARD_2026-11-09.md).

**Required response:** either specify a coherent joint objective and valid updates, or present the algorithm explicitly as a staged estimator with separate objectives and validation. Assess profile, membership, and nuisance stability separately. Preserve the current failure cases; do not replace them with a convenient passing seed or weaken their assertions.

**5. High: reference and nuisance corrections are not independently anchored biological corrections.**

The reference mixture is optimized against the query's pooled counts under a Poisson composition objective before the clone interpretation is established. A genuine dosage-shaped signal can therefore be represented by a different reference mixture. A simple exact example is sufficient: reference A has composition `(0.5, 0.5)`. Multiplying its first region by 1.5 and renormalizing gives `(0.6, 0.4)`. If reference B has that composition, the query can be fit perfectly by B with no dosage effect. Moreover, the reference-heterogeneity term is zero at a pure reference weight, so a concentrated mixture does not establish a well-matched baseline. See [clones.py:462](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:462) and [clones.py:713](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:713).

Reference rates pool counts within dataset/ref_group, not donor-balanced independent units. Deeper cells and more represented donors contribute more. That may be the intended abundance-weighted estimand, but it is not an equal-donor population baseline. Reference sampling uncertainty and donor/batch transport are not integrated into cell responsibilities. Held-out reference cells from the same groups provide a useful diagnostic, not independent-donor validation. See [pipeline.py:572](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:572) and [pipeline.py:715](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:715).

Depth correction also has an origin and transport problem. Even reference-derived slopes are applied around query-centered log depth. Slopes alone do not identify the reference population's intercept or depth distribution. In the default query-derived mode, depth can track state, CNA burden, or spatial cell mixture, so regression may remove biological signal. The quintile estimator returns zero when it lacks enough usable groups; zero is then an estimation limitation, not evidence that depth has no effect. See [clones.py:654](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:654) and [clones.py:725](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:725).

Finally, high-passing normal programmes intentionally prevents them from explaining a whole-arm shift. This protects arm-scale dosage under that assumption, but also forces arm-scale transcriptional state outside the nuisance span. It cannot discriminate the two causes. Shorter genuine events may instead compete with programme loadings. See [programmes.py:61](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/programmes.py:61).

**Required response:** define the reference population and weighting estimand; record donor/batch composition and depth support; evaluate normal-to-query transport. Report sensitivity to scientifically justified references without selecting one by DNA agreement on evaluation data. A statistically better fit to one conditional objective does not establish improved dosage recovery.

**6. High: a sum of heterogeneous gene counts need not have the NB law used for thinning.**

If independent gene counts have NB means m_g and dispersions phi_g, their sum has

`E(Y) = sum_g m_g` and `Var(Y) = sum_g m_g + sum_g phi_g*m_g^2`.

A moment-matched bin dispersion is `sum(phi_g*m_g^2)/(sum m_g)^2`. Matching these two moments does not make the distribution NB. Exact NB closure requires the common success-probability condition, equivalently a shared Gamma scale. Equal gene dispersions with unequal gene means generally do not satisfy it. This matters because the simulator generates heterogeneous gene means and the production procedure aggregates those genes before applying beta-binomial thinning. See [simulate.py:108](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/simulate.py:108), [bins.py:29](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bins.py:29), and [clones.py:169](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:169).

The exact thinning theorem is distributional and requires the relevant nuisance parameter. It does not follow from count conservation, approximate moments, or near-zero empirical correlation. These conditions are explicit in [Neufeld et al., Data Thinning for Convolution-Closed Distributions](https://www.jmlr.org/papers/volume25/23-0446/23-0446.pdf), including the NB closure condition and nuisance-parameter discussion. In HaploDose, the reference-estimated bin size is frozen for splitting, while cell-likelihood dispersion can later change. The updated likelihood is not automatically the true law of those previously generated folds.

There is a second qualification: L is computed from the same observed bins being scored. A product of independent NBs with means summing to L is a working/composite likelihood when L is the observed total; it is not automatically the normalized joint conditional distribution given that total. This distinction is material for posterior calibration and null simulation, even if ranking performance is useful.

**Required response:** test the observation law and nuisance estimation jointly with the full workflow. The current conservation and NB-moment tests are valuable computational checks but insufficient statistical validation. Do not conclude that all recovery failures are caused by binning: existing oracle gene-versus-bin comparisons retain substantial bin-level assignment information. A gene-level replacement needs its own evidence, especially with estimated rather than known noise.

**7. High: regional matching and gene consistency can evaluate different populations.**

`matched_differences` constructs pair-specific depth-matching weights and pools the corresponding test counts. `gene_consistent` instead takes unweighted gene sums from MAP-labeled cells in a bounded proposal subsample. It does not receive those matching weights. Thus a depth/state imbalance can affect the gene veto even when the regional contrast is based on a different, matched population. Coarse depth strata also do not ensure exact equality of the continuous depth distributions within each stratum. See [clones.py:955](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:955) and [clones.py:1037](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1037).

The one-sided sign test is not calibrated after selecting intervals and genes from these same observations. Coexpressed genes are not independent Bernoulli signs, and real CNA responses need not shift every gene. The current code also bypasses the veto when gene/cell support is insufficient. “No veto” then means unavailable evidence, not positive independent confirmation.

This concern should not be misidentified as the known cause of the current three production recovery failures. The [production gene-pool diagnostic](/gpfs/commons/home/jameslee/Haplodose/docs/PRODUCTION_GENE_POOLS_2026-11-09.md) reports that all 54 observed production gene checks passed. Research fixes that aligned weighted populations improved some exposed recoveries but subsequently failed state-null or weight-unit checks. Those outcomes constrain the recommendation.

**Required response:** expose the population, weights, depth overlap, and missing gene support for each comparison. Align the estimands before evaluating a replacement veto, then assess selection and dependence. Preserve negative results rather than promoting a weighted sign test merely because it recovers five familiar fixtures.

**8. High: routine segment uncertainty is a conditional scale diagnostic with demonstrated coverage failure.**

The exported approximate SE is based on an unweighted mean of bin log ratios, using independent-bin variance terms. It omits fitted references, memberships, nuisance uncertainty, spatial/genomic correlation, and interval selection. It is also a different object from the posterior mean state effect. Current metadata correctly state that no validated parameter-interval estimand is supplied. See [profile.py:152](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/profile.py:152).

I independently recounted coverage from the numeric lower/upper endpoints and true effects in the saved 4.2 MB [interval table](/gpfs/commons/home/jameslee/Haplodose/dev/recovery_20261107/uncertainty_simultaneous_v1/intervals.tsv), rather than trusting its `covered` field:

| Conditional generator | Fixed interval, legacy scale | Selected interval, legacy scale | Explicit simultaneous helper, fixed / selected |
|---|---:|---:|---:|
| Independent, homogeneous | 934/1000 | 818/1000 | 1000/1000 / 1000/1000 |
| Independent, heterogeneous | 985/1000 | 913/1000 | 1000/1000 / 1000/1000 |
| Correlated, heterogeneous | 703/1000 | 54/1000 | 1000/1000 / 1000/1000 |

The preserved analysis uses Holm correction over six one-sided undercoverage tests; five reject nominal 95% coverage in this replication. These are conditional simulation endpoints, not clinical or biological coverage estimates. Selection includes the largest absolute mean among HMM runs, so these numbers should not be described as coverage for every ordinary emitted segment.

The optional [systematic_rna_intervals helper](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/uncertainty.py:12) has a coherent and much narrower contract: conditional Poisson counts, independently fixed expectations, known Gaussian covariance, and a frozen eligible bin universe. It allocates error to simultaneous count and all-interval Gaussian bounds. Its conservatism in the saved study does not make fitted residual covariance, median pseudocounts, or fractional responsibility pools eligible inputs. Its substantial width/power cost must remain visible.

A further issue applies upstream: fractional weighted pools are passed to a PLN count working likelihood. For fixed weights, `Var(sum wY)` involves squared weights, not the sum of weighted means alone; fitted weights add further dependence. The depth-matched difference helper partly accounts for squared weights, but this does not turn every pooled profile into an exact count likelihood.

**Required response:** keep scale diagnostics separate from parameter intervals. Do not attach the approximate SE to the HMM effect or reinterpret three SEs as a minimum detectable DNA change. A deployable uncertainty method must cover its actual estimand under estimated references, selection, and dependence, with useful power.

**9. High: frozen confirmation isolates execution roles but does not establish genomic confirmation.**

The confirmation implementation has real strengths: disjoint cell IDs, frozen discovery arrays/contrasts, proposal-only assignment and nuisance scoring, exact-tie abstention, and no feedback into discovery. Those are meaningful protections. See [confirmation.py:128](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/confirmation.py:128) and [pipeline.py:435](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:435).

Nevertheless, proposal and test bins can share cell state, donor effects, and transcriptomic dependence. Test summaries use observed count-based delta-method NB variances without cross-feature covariance. Nominal one-sided p-values are explicitly unadjusted, with no specimen-wide decision rule. The test contrast is region-versus-background pooled RNA odds; it is not necessarily the same depth/programme-adjusted estimand that selected the discovery contrast. Candidate-specific programme scores assist assignment, but are not a full correction of the test-group contrast.

Finally, discovery fits can use all discovery features, and neighboring cells in spatial data are not independent biological replicates merely because their IDs differ. Frozen execution should not be equated with a confirmatory genomic experiment.

**Required response:** retain this module as an explicitly labeled candidate RNA replication diagnostic. If formal testing is pursued, specify the test estimand, assignment-induced selection, covariance assumptions, multiplicity family, and sampling unit. Treat genotype validation as a separate requirement. There is no evidence here that simply applying BH to the existing nominal p-values would make them calibrated.

**10. High: allele evidence adds information but does not remove genomic ambiguity or model-selection effects.**

The allele model distinguishes balance from oriented imbalance and handles phase switches. Current joint labels correctly avoid declaring copy-neutral LOH from relative-neutral RNA plus imbalance. These are substantive strengths. The use of pooled PLN expression, beta-binomial alleles, phase modeling, and per-cell allele scoring has precedent in [the Numbat primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10289836/); that precedent does not validate HaploDose's particular search, priors, thresholds, or data handling.

Four unresolved conditions matter:

- The default cell allele likelihood remains binomial. An independently supplied cell-level beta-binomial concentration is available, but the pooled concentration is not automatically a valid cell-level parameter. Stochastic ASE and linked molecules/sites can create more variation than the default assumes. See [alleles.py:323](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:323).
- Normal bias filtering and concentration estimation reuse the retained SNP set; the first query-pool fit additionally selects concentration multiples and phase-switch rate by likelihood. Reported posteriors condition on those selected quantities. See [alleles.py:166](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:166) and [alleles.py:428](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:428).
- Even/odd SNP roles do not make the full allele workflow independent: imbalance intervals and oriented fractions come from pooled HMM evidence, and cell memberships are fitted using allele data. Shared allele-expression state can persist across reserved SNPs without a DNA clone difference. The scientific contract already records such a falsification.
- `read_alleles` identifies sites by chromosome/position, orients each row by its GT, sums duplicate cell/site rows, and takes first metadata values. It does not establish donor-specific phase consistency, REF/ALT identity, deduplicated molecules, or that duplicated rows are disjoint observations. Those require an upstream contract and integrity evidence. This is an interface limitation, not a claim that the inspected datasets contain such errors. See [alleles.py:79](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/alleles.py:79).

The allele veto also compares summaries under assumed one-copy imbalance behavior. Dividing an RNA contrast by fitted beta does not identify the underlying allele-specific copy states. The existing exemptions for LOH-like states and large changes reduce particular errors but do not supply a full joint DNA model.

**Required response:** report actual reference/allele support, assay-specific molecule counting, phase identity, and the selected noise model. Evaluate shared ASE, mapping bias, sparse support, and mirrored imbalance as separate failure modes. Independent RNA SNPs alone cannot establish genomic identity when the competing state mechanism has the same allele-expression law.

**11. High: spatial rows are assigned exclusive components even when their biological content is mixed.**

`neighbour_graph` constructs a directed, row-normalized nearest-neighbor graph. The posterior update adds `lambda * G R`; adaptive lambda maximizes an in-sample neighbor-prediction criterion. This is a smoothing heuristic rather than optimization of a symmetric Potts energy, as the current graph docstring correctly states. A positive selected lambda can occur under a spatially random finite sample. See [clones.py:763](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:763) and [clones.py:1825](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1825).

More fundamentally, the likelihood uses one latent component per row. Soft assignment probabilities represent uncertainty over that exclusive label; they are not within-spot clone fractions. For multicellular spots, changing proportions of tumour, stromal, and immune cells can generate smooth regional RNA contrasts. Generic x/y metadata do not establish one-cell observation units, spatial barriers, or appropriate neighbor distances. Spatial smoothing can reinforce coherent cell state as readily as genotype.

**Required response:** distinguish segmented single-cell spatial data from mixed spots in the supported scientific scope. Test genotype-state discordance, mixed spots, boundaries, and coordinate/section integrity. If mixed-spot clone fractions are a target, they require a measurement model for mixtures and truth at that resolution. Keep cell-ID tie handling, but do not treat graph determinism as biological validation.

**12. High: bulk beta is a working amplitude, and median pseudocounts do not inherit count likelihood guarantees.**

The bulk module describes beta as absorbing purity times dosage response. Even in a simple linear RNA mixture with tumour fraction p, equal per-copy expression, and diploid normal cells, the regional effect is

`log((1-p) + p*C/2)`.

This is generally not `beta*log(C/2)` for one shared beta across all C. At p=0.5, matching the two formulas requires beta=0.4150 for C=1, 0.5503 for C=3, and 0.5850 for C=4. These values were independently calculated here. Different RNA yield per cell, compensation, ploidy, and compositional normalization add further ambiguity. Beta should therefore be called an effective state-amplitude parameter, not a purity estimate. See [bulk.py:1](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:1), [bulk.py:74](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:74), and [profile.py:112](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/profile.py:112).

With `bin_stat='median'`, per-gene median log ratios are exponentiated into pseudocounts and passed through the PLN path. This can be a useful robust estimator, but its sampling distribution is not Poisson just because it has count-like units. The current run metadata labels the approximation; downstream interpretation must preserve that label. Arm calls similarly use a normal-panel median/MAD scale and a fixed absolute z threshold, without calibrated multiplicity or predictive error control. Leave-one-out fitting does not remove all dependence because feature selection and the final arm reference scale use the shared normal panel. See [bulk.py:108](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:108), [bulk.py:266](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:266), and [bulk.py:297](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bulk.py:297).

The recent joint bulk uncertainty evidence does not justify a strong-baseline superiority claim. I recounted 80 comparison records: 64 testable, 16 untestable, and six passing conditional gates. All six concern Pearson correlation for the two continuous-amplitude candidates against `hd_lr`, `naive`, or `naive_smooth`. None establishes an advantage over the stronger posterior HMM or naive segmentation controls. The analysis has 20 original normal donors, 24 calibration patients, and 24 comparison patients; one million bootstrap draws and eight resampled reference fits do not increase biological n. Its Holm family is 80, with separate Monte Carlo uncertainty accounting. See [saved summary](/gpfs/commons/home/jameslee/Haplodose/dev/bulk_amplitude_uncertainty_20261109/summary_v2/summary.json) and [the full report](/gpfs/commons/home/jameslee/Haplodose/docs/BULK_AMPLITUDE_UNCERTAINTY_2026-11-09.md).

**Required response:** separate continuous RNA estimation from discrete state classification and purity. Compare on the same patients, genes, callable DNA regions, and fixed false-positive operating points. Preserve undefined endpoints and reference sensitivity. Do not replace an unresolved biological comparison with more bootstrap precision.

**13. Moderate, direct defect: the advertised default ribosomal exclusion depends on gene-ID namespace.**

The CLI advertises ribosomal protein genes as excluded by default. `CountSource` selects Ensembl IDs for `gene_namespace='ensembl'`, and `make_bins` applies the ribosomal symbol regex directly to those identifiers. Valid Ensembl IDs start with ENS; they cannot match the RPL/RPS symbol expression. Neither `read_positions` nor `make_bins` maps those IDs to symbols for this filter. Equivalent counts supplied under Ensembl IDs can therefore retain genes removed under symbols. See [cli.py:24](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/cli.py:24), [io.py:130](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/io.py:130), and [bins.py:41](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/bins.py:41).

This is scientifically relevant because the exclusion was introduced to prevent coordinated, highly expressed ribosomal clusters from driving groups. The filtering difference follows directly from the source; no inference run is needed to establish it. [Unverified] Its effect on any particular dataset's assignments has not been measured in this review.

**Recommended correction:** apply biological category filtering through a declared annotation mapping, or explicitly require namespace-matched exclusion IDs and report when automatic category filtering cannot be applied. A focused check should compare the retained biological genes for equivalent symbol and Ensembl inputs. Do not silently guess mappings or strip version suffixes without identity checks.

**14. Moderate, source concern: final allele parameters can be updated after the returned assignment scores.**

At the end of `SpecimenModel.em`, the code computes LL, computes R, pools those responsibilities, and then calls `set_allele_pools`. That call replaces `B` and `allele_post`. The method returns the earlier LL and R without evaluating the likelihood under the newly stored B. `loglik` explicitly uses `self.B`, and final segment export uses the updated allele posterior. Therefore, when that last allele update changes B, the returned assignment scores and the stored allele model do not describe the same parameter snapshot. See [clones.py:1294](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1294), [clones.py:798](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:798), and [pipeline.py:640](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:640).

This is distinct from the intentional `log_dose` versus `assignment_dose` separation already documented for RNA. [Unverified] The size and frequency of the allele mismatch, including whether it changes a MAP label in realistic data, were not measured here. A finite-iteration estimator need not converge to a fixed point, but it should identify which parameter state generated each exported responsibility.

**Recommended correction/check:** on a nonspatial allele-enabled fixture, compare the returned relative LL with a fresh `model.loglik(res['doses'])` before downstream mutation. Then either preserve/export the assignment-time allele snapshot or finalize scores against an explicitly frozen final snapshot. Do not describe a single extra iteration as proof of convergence. This dynamic regression was not run.

**15. Moderate: documentation and legacy field names still contradict the narrower contract.**

The README links `related_work.md` as “model details,” but that file still says each segment reports a minimum detectable log dosage, describes even/odd separation as cross-fitting without the adaptive qualification, and lists now-implemented allele/spatial functions as planned. The manuscript draft and outline retain the assertion that spatially random clones give lambda=0. Several source docstrings also use stronger language than their implementation supports, including reference depth slopes that “cannot absorb” tumour copy number and gene consistency premised on every gene shifting. See [related_work.md:42](/gpfs/commons/home/jameslee/Haplodose/docs/related_work.md:42), [manuscript_outline.md:108](/gpfs/commons/home/jameslee/Haplodose/docs/manuscript_outline.md:108), [manuscript draft:187](/gpfs/commons/home/jameslee/Haplodose/docs/manuscript_draft_20261103.md:187), and [clones.py:683](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:683).

Runtime exports retain `cn_state`, `malignant_call`, `normal_like`, and `p_normal`. Their metadata contain appropriate caveats, but a standalone table can still be misread. `n_eff` in the clone table is summed responsibility mass, not a Kish effective sample size or biological replication count. `summarise` assigns every eligible discovery row to a MAP component, including low-confidence or exact-tie cases; unlike confirmation scoring, it has no general ambiguity abstention. Component retention/merging also relies on support rules: failure to resolve a difference, including missing depth overlap, is not evidence of genomic equivalence. The logs acknowledge this, but the main labels are easier to consume than the logs. See [pipeline.py:595](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:595), [pipeline.py:346](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/pipeline.py:346), and [clones.py:1104](/gpfs/commons/home/jameslee/Haplodose/python/haplodose/clones.py:1104).

**Recommended correction:** make one current methods specification authoritative and label historical documents clearly. Carry candidate/uncalibrated status into consumable tables; distinguish assignment confidence, lack of information, and evidence of a difference. Preserve backward-compatible fields if needed, but do not rely on a distant JSON caveat to define their scientific meaning.

**Strengths that should be preserved.**

- Explicit raw-count, assembly, coordinate, and feature-namespace declarations; rejection of fractional/nonfinite counts and ambiguous feature IDs; exact representability checks before float32 aggregation.
- Reference-only feature selection and independent dataset subruns; held-out cells do not select the gene universe in the single-cell pipeline.
- Stable log-likelihood-ratio arithmetic, normalized HMM messages, explicit mixture optimizer checks, and retention of numerical ambiguity rather than arbitrary threshold decisions.
- Corrected squared-weight sampling terms for matched pools; separation of assignment profiles, diagnostic profiles, observed pools, expected pools, and centering constants.
- Fresh-output guards, package/native hashes, environment inventories, failure manifests, seeds, and preservation of historical failures. These make scientific conclusions more auditable.
- Rejection of unsupported joint null configurations, explicit unvalidated uncertainty metadata, exact-profile redundancy repair, and preservation of simple baselines and failed research candidates.

These are evidence of careful computation and reporting. They are not substitutes for the validation requirements above. The historical [403-pass Python test log](/gpfs/commons/home/jameslee/Haplodose/dev/exact_profile_merge_20261109/validation_23121461/pytest.txt) and [four-pass Rust test log](/gpfs/commons/home/jameslee/Haplodose/dev/exact_profile_merge_20261109/validation_23121461/cargo_test.txt) were inspected. They belong to that recorded build; one GPU-only Python test was skipped. I did not rerun them or establish that their wheel is byte-identical to every file in this working tree.

**A focused scientific development sequence.** These are recommendations, not jobs or implementation authorized by this review.

1. Fix the bounded namespace and output-snapshot issues, and reconcile the current methods text. Verify those behaviors with small identity/consistency fixtures on a compute node. Preserve all existing thresholds and frozen outputs.
2. Choose the primary object: regional relative RNA estimation, candidate RNA grouping, or genomic clone inference under additional evidence. Write the observation unit, reference population, covariates, effect units, uncertainty target, and allowed biological interpretation in one specification.
3. Resolve or explicitly delimit the observation law. Separate raw gene counts, binned counts, fixed weighted pools, fitted weighted pools, and median pseudocounts. State which likelihood is exact, approximate, or composite at each stage. Do not infer independent thinning from moment agreement.
4. Use the existing failed cases to compare a small number of justified mechanism changes. Keep a strong reference-normalized pseudobulk/gene baseline and a simple segmented baseline. Measure error, genuine contrast recovery, all-cell assignment accuracy, and abstention together. An improvement on one endpoint must not hide a decline on another.
5. Freeze the complete revised method before independent biological evaluation. For aggregate dosage, use independently measured corresponding DNA at matching resolution. For clone identity, use corresponding-cell labels or a justified external clone linkage. For mixed spatial spots, evaluate fractions or aggregate regions rather than asserting single-cell identity.
6. Treat donors/specimens/independent lines as biological replication, account for shared references and repeated samples, and freeze the multiple-testing family. Keep missing calls, uncertain DNA, and failed methods visible. Estimate uncertainty for the actual paired population and scientific estimand.

A manuscript can be organized around these estimands, assumptions, mechanisms, and failure modes. The requested presentation model, [Wang et al.'s 16-method cell-potency review](https://pubmed.ncbi.nlm.nih.gov/42956748/), explicitly motivates unified notation and careful sensitivity analysis; it is a presentation reference, not a CNA comparator or validation dataset. A new methods claim for HaploDose still needs a demonstrated contribution beyond its combination of existing likelihood, HMM, pooling, and grouping components. The most defensible current material is the computational architecture and the reproducible evidence delimiting what its RNA outputs can and cannot support.

**Verification performed for this review.**

| Check | Outcome |
|---|---|
| Initial `git status --porcelain=v1` and `git rev-parse HEAD` | Baseline recorded; 23 modified tracked paths, 614 untracked entries |
| Bounded `rg`, numbered source reads, selected JSON/TSV inspection | Core execution path, representative tests, benchmark/scientific evidence reviewed |
| Standard-library `ast.parse` of all package Python modules | 21 modules parsed successfully; syntax only |
| Saved five-seed thinning recovery recount | 2/5 passed; no recorded execution errors |
| Saved complete-search evaluation recount | 600 records; table in finding 2 reproduced; two-clone median ARI reproduced |
| Saved exact-merge research-panel recount | 600 records; matched null 0/100, state null 12/100 in each arm; two-population exact recovery 99/100 and 98/100 |
| Coverage recomputed from numeric interval endpoints | All twelve 1000-interval coverage counts reproduced |
| Independent binomial-tail inversion | 8/200 family-9 upper 0.0896385514; 12/100 family-6 upper 0.2189094284; zero-error n=59 boundary reproduced |
| Bulk comparison-record recount | 80 contrasts, 64 testable, six conditional gates; passing comparisons identified |
| Historical build logs | 403 Python passed, one skipped; four Rust passed, for the recorded historical build |
| Source SHA-256 recheck | Reviewed package/native sources unchanged across the review checks |
| New inference, scientific reruns, current-wheel tests, GPU parity | **Not run**; review performed without a compute allocation or new job submission |

The standard-library evidence scripts only read saved tables/JSON and performed bounded arithmetic; they did not import HaploDose, fit a model, or execute a test suite. Initial `python` invocation was unavailable; inspection used `python3`. The review did not install dependencies, alter environments, create branches, commit, or publish.

**Selected source identity.** SHA-256 values below identify the actual inspected files, including uncommitted changes. Remaining source hashes and the initial Git listing were also recorded in temporary review files; this table keeps the load-bearing identities in the durable report.

| Source | SHA-256 |
|---|---|
| `python/haplodose/clones.py` | `343bb77b7dc0811bdeabccaaa78ffb49cecaee5eca002c8d085abc39d5f709de` |
| `python/haplodose/pipeline.py` | `7f57ac4c63d7b2ca37b58a823464eedfd3d9b596487bf0c7aa11710eb69fb7c5` |
| `python/haplodose/profile.py` | `7db040cbb9503bb6f48b64465af23fc39eb7190f8e79c8c02f539dd7f5cc3612` |
| `python/haplodose/alleles.py` | `112109a2066158b069c4e3fb87ca5417ff0a1de7164625d6a3211a609deb958e` |
| `python/haplodose/bulk.py` | `b4d44361d059239fb363f77ebac8b0e5b21e0909830291be1e7adf45939de79b` |
| `python/haplodose/confirmation.py` | `d7dcc9bd9ce09b7e65c45a7adc0f462ee18a64e154caf2df3b840ea8d274b1ac` |
| `python/haplodose/bins.py` | `0c3acfd3817dd98a4e3c9bdd8573fb986b679ca4aaac4591cf6e2fca39d5eabb` |
| `python/haplodose/kernels/estep.py` | `95c9fcfd95991dab105f102d133e5848f2a113149475dc1a3fe7174781d94aab` |
| `crates/core/src/aggregate.rs` | `e35f1172c9a2118f815e555ee36fe454385ca2ee5b1bcdafd5e709cb632ce525` |
