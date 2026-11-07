# Recovery failures and the false-split/power tradeoff

The goal is to distinguish computational defects from failures of inference and
identification, repair supported defects, and retain the original recovery and
error/power criteria. This follow-up uses frozen simulation inputs as development
and regression evidence. It does not create new biological validation.

## Outcome

One numerical defect is repaired: a failed L-BFGS line search no longer immediately
removes a potentially valid proposal. The identical bounded objective is retried
with SLSQP, then must pass all the existing optimizer-success, projected-score,
curvature and objective-error checks. Neither numerical certification nor scientific
cutoffs were relaxed. Compute job 22846715 passed **311 Python tests, one GPU-only
skip, and four Rust tests** from an isolated noneditable wheel.

The three complete recovery assertions still fail. The repair restores an omitted
proposal in one case, but does not solve its later merging or normal-cell errors.
Two nuisance ablations and a more accurate thinning-law diagnostic do not justify
additional default changes. The paired 998-fit false-split/power replay has the
same endpoint counts as before; the 5% error criterion remains unmet.

The original positive-control assertion passes in the full suite, and the other
previously passing seed (20261107) passes replay job 22848498. The five-case panel
therefore remains **2/5**, with all three previously failing cases still failing.
These repetitions check regressions; they are not independent new evaluation.

## Where each case fails

All three source matrices and fitting metadata are the exact saved files from the
previous five-case panel. Truth columns were removed before fitting and joined only
afterward. Return-event instrumentation records existing intermediates without
changing calculations. The original recovery assertion remains unchanged.

| Seed | Traced mechanism | Independent distinguishing evidence | Decision |
|---|---|---|---|
| 20261108 | L-BFGS reported `ABNORMAL` and the first-half chr12 proposal was discarded, despite small scores and positive local curvature. After retry, this proposal passes initial confirmation, but the groups are subsequently merged. Normal false calls rise from 9/90 to 10/90. | Retry gain 24.237872389437065 agrees with full-NB-PMF/Nelder–Mead optimization within 2e-12. Final posterior arithmetic agrees independently within 1.8e-15. | Retain the numerical correction, including the worse normal-error result. Do not claim full recovery. |
| 20261109 | The relevant proposal exists, but its test-half pooled contrast is only -0.0682 against an HMM alternative of magnitude 0.2189. It is rejected during initial agglomeration, before depth matching or gene veto. | Maximum nonneutral posterior in the target region is 0.00889. Supplying true group labels only after fitting yields mean -0.2173 and minimum posterior 0.9930, passing the unchanged HMM and gene checks. | The inferred partition is the immediate bottleneck; lowering the posterior cutoff is not justified. |
| 20261110 | The corresponding proposal exists, but its pooled contrast is -0.0752 against an alternative of magnitude 0.2243. It fails the same initial HMM gate. | Maximum nonneutral posterior is 0.000878. The post-fit true-group contrast is -0.2264, with minimum posterior 0.999998 and passing gene support. | Retain the failed assertion. This is not an indexing error or a missing supported HMM core. |

The estimated quantity in the HMM comparison is a contrast between **assigned
cell pools**. That can be much smaller than the latent component contrast when
assignments mix populations. These are different estimands. The post-fit true-label
probe establishes that the binned test data can support the contrast when grouping
is known; it does not make those labels available to the method. A replacement
confirmation procedure would need to account for uncertain membership and retain
whole-workflow error/power assessment, not reinterpret the existing posterior.

For the numerical counterexample, the original local gain-error estimate was
1.64e-16, below its unchanged 4.73e-8 budget. The independent likelihood calculation
has its own subtraction-error budget. The retry performs actual optimization;
it does not simply override the failed status flag. A regression test exercises
this failure mode while retaining the original numerical requirements.

## Information in the measurements versus information used by the method

Job 22847127 regenerated each simulator realization and verified exact equality of
CSR data, indices, row pointers and cell identifiers against the frozen counts.
It then computed an independent gene-level NB likelihood using the **true** gene
rates, dispersion, per-cell depth and population proportions. These quantities are
not supplied to HaploDose.

| Seed | Oracle normal false calls | Oracle tumour detection | Oracle tumour assignment accuracy, all genes | Retained genes |
|---|---:|---:|---:|---:|
| 20261108 | 0/90 | 510/510 | 86.47% | 86.47% |
| 20261109 | 0/90 | 510/510 | 80.59% | 80.59% |
| 20261110 | 0/90 | 510/510 | 82.16% | 82.35% |

This rules out simple absence of information in the source counts as a sufficient
explanation. Retaining the pipeline's genes does not account for the observed
loss. The oracle changes several unavailable nuisance inputs together, so it does
**not** identify a single deployable correction or prove that estimated nuisance
parameters can attain this performance. Binning, nuisance estimation and latent
membership remain distinct mechanisms to assess.

The posterior/prior decomposition was checked separately after inference. With
independently recomputed `softmax(loglik + log_prior)`, the maximum discrepancy is
1.4–1.8e-15; the smallest decision margin exceeds 0.0047. Numerical rounding does
not explain the classifications. Of 10 false normal-cell calls in repaired case
20261108, seven have likelihood evidence favoring normal but are assigned to the
more prevalent tumour population after the prior is applied. Its fitted normal
proportion is 0.1444, near the generating 0.15. This is a classification/error-cost
issue under the fitted model, not evidence that the prior was coded incorrectly.
Choosing a prior or posterior threshold from the known labels is not a correction.

## Mechanism checks that did not support a default change

Two one-component controls reused all five frozen inputs and all six original
scientific criteria: reference-derived depth slopes, and removal of programme
correction. Each passed only **1/5** cases, losing the original passing seed
20261031 and recovering none of the three failures. S3 remained a single tumour
group. These results do not justify either setting as a new default. They do not
exclude every possible interaction between nuisance estimation and grouping.

A separate saved-proposal calculation integrated the actual thinning construction
under the assumed full-count NB law:

`H | P ~ NB(full_mean * P, full_size)` and
`P ~ Beta(frozen_size/2, frozen_size/2)`.

Gauss–Jacobi integration at 32, 64 and 128 nodes gave stable likelihoods and labels.
It accounts for the mismatch between the frozen splitting size and subsequently
refitted dispersion. It did not improve proposal-label accuracy: case 20261109
changed from 69.4% to 68.8%, and case 20261110 from 71.71% to 71.12%. Thus this
more elaborate likelihood is not supported as the recovery solution. This result
does not establish that the existing working-NB marginal is exact, or that the
halves are independent after dispersion refitting. Those assumptions remain
unverified; routine confirmation does not receive an independence claim.

## Error, power and biological identification

The paired replay preserves the two independent 199-specimen calibration banks,
200 evaluations per condition, the original generators and the nine-member
simultaneous one-sided confidence family. These previously examined simulation
seeds are **regression evidence**, not a new independent evaluation.

| Procedure | Matched NB null: errors / upper bound | State null: errors / upper bound | Two-clone condition: detections / lower bound |
|---|---|---|---|
| Raw groups | 0/200; 2.56% | 184/200; 96.10% | 200/200; 97.44% |
| NB-calibrated gate | 0/200; 2.56% | 184/200; 96.10% | 200/200; 97.44% |
| State-calibrated gate | 0/200; 2.56% | 8/200; 8.96% | 65/200; 24.34% |

Raw membership median ARI in the two-clone condition remains 0.0764. The gate does
not improve memberships and continues to abstain from supporting them. The
original event counts eligible tumour-group distinctions, excluding the
reference-like group; it is not the entire biological false-distinction contract.
No error cutoff, confidence family, sample size or comparator was changed to
obtain a success verdict.

Job 22848584 independently recounted the decisions and inverted binomial tail
probabilities for the confidence bounds. All nine clone bounds agree to displayed
precision. It also rechecked the preserved uncertainty experiment without rerunning
or changing it; its undercoverage conclusions remain unchanged. The report is
`dev/recovery_next_20261107/evidence_v1/report.json`.

An independent identification check explains why this cannot be solved solely by
optimizing the clustering score. In the generator, a regional CNA effect
`delta = 0.6 log(1.5) = 0.2433` is added to Gaussian state variation of SD 0.25.
The resulting two-component marginal density is continuous and unimodal: with
`a = delta/2`, the centered density derivative has the sign of
`-x + a tanh(a*x/sigma²)`, whose derivative is bounded above by
`-1 + a²/sigma² = -0.7633`. Reinterpreting that same total regional effect as
transcriptional state leaves the normalized RNA rates and NB emissions identical.
RNA observations alone cannot select its biological label.

A separate IID diagnostic compares that latent mixture with a moment-matched
Gaussian state law. Numerical quadrature converges to KL divergence 0.0001249822
per cell. Pinsker's inequality bounds level-5% power at 600 IID cells by about
24.36%, even if the latent effect were directly observed. **This is a different,
matched-moment null and IID alternative, not the original fixed-balanced
benchmark.** It does not invalidate or relabel the original 32.5% detection result.
The stronger exact observational-equivalence example does not depend on this
bound. Genomic interpretation needs independently justified restrictions on state
variation or corresponding DNA evidence.

## Reproducibility and remaining limits

Evidence is under `dev/recovery_next_20261107`: `baseline`, `trace_20261108/109/110`,
`after_retry`, `followup_v1`, `postfit_20261108/109/110`, `oracle_v1`, `law_v1`, the
two `ablation_*_v1` directories, and `rank_v1`. Each computation uses an immutable
source snapshot and a recorded isolated wheel. Initial sources, patches and all
previous `dev/recovery_20261107` results are preserved. New results never replace
old outcomes. The production edit is in `python/haplodose/clones.py`; its new
regression test is in `tests/test_thin_optimization.py`. Diagnostic runners are
under `benchmarks/hardening`.

The bounded preservation audit verified 331 initial source files, 331 files in
the previous final source, nine recorded historical-result hashes, the original
four criteria/generator files, and three wheel/native pairs. All 22 live package
source/data files match the candidate tested in job 22846715. This is not an
exhaustive rehash of historical raw data. `preservation_audit.json` records the
checks; `jobs_final.tsv` records resources, runtime, memory and exit status for
all 21 jobs in this follow-up. All completed with exit code zero. A completed
diagnostic job is compatible with a failed scientific assertion reported inside it.

The numerical repair is supported. The three full recovery failures, the 5%
false-split guarantee and reliable CNA membership are **not resolved**. The new
analyses locate the failures more precisely and reject several tempting fixes;
they do not establish biological validity or justify releasing a newly calibrated
clone caller. No biological replication has been added by these computations.
