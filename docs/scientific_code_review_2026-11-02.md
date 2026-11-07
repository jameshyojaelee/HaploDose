# HaploDose scientific and code review

Reviewed 2026-11-02. Saved from the preceding conversational review. Findings describe the source inspected during that review; subsequent or concurrent edits have not been reassessed here. No source files were modified by the review.

The review found material correctness and validation issues that should be resolved before using HaploDose’s clone-confidence or superiority results as manuscript evidence. The foundation is useful: RNA count models, bounded count processing, stable likelihood calculations, and targeted regression tests.

## Prioritized findings

### 1. P1: Comparisons with no depth overlap can produce spurious difference segments

`matched_differences` ([clones.py](../python/haplodose/clones.py), reviewed line 498) can give both groups zero matched weight. `difference_segments` ([profile.py](../python/haplodose/profile.py), reviewed line 173) then divides by zero, and its posterior threshold does not reject NaN.

**Reproduced:** two groups with disjoint depths returned a difference segment with `mean_diff=NaN`, `z=NaN`, and `posterior=NaN`. This nonempty segment list can preserve an unsupported clone distinction. Validate matched depth and finite statistics explicitly; report comparisons without adequate overlap as untestable.

### 2. P1: Allele logic can discard real genomic differences

The allele veto ([alleles.py](../python/haplodose/alleles.py), reviewed line 375) treats sufficiently precise, unchanged allelic imbalance as evidence against a dosage difference. Balanced gains, such as local 1:1 → 2:2 copies, retain the same allele fraction. The implementation does not implement the exception mentioned in its docstring.

**Reproduced:** identical balanced counts at 100 SNPs triggered the veto. Separately, allele comparison (reviewed line 345) returned zero difference for allele fractions 0.9 versus 0.1 because it compares squared imbalance magnitude. Replace the universal veto with state compatibility and retain signed haplotype contrasts.

### 3. P1: The claimed split-validation independence does not survive EM

The EM loop ([clones.py](../python/haplodose/clones.py), reviewed line 686) assigns cells using all bins. Those assignments subsequently define groups tested on odd bins, which have already influenced membership. Odd bins are also reused across adaptive split decisions.

Consequently, the even/odd partition does not establish an independent validation set throughout the algorithm. A 0.95 HMM posterior is not demonstrated false-split control. Calibrate the **entire adaptive procedure** on specimens with no true subclone, including depth, expression-state, reference-mismatch, and spatial confounding.

### 4. P1: Benchmark methods are not consistently evaluated against identical targets

[evaluate.py](../benchmarks/kinker/evaluate.py) (reviewed line 42) calculates DNA arm truth using each run’s own retained genes. Different gene filters therefore change the truth and potentially the eligible arms between methods.

Additionally, the [CopyKAT/SCEVAN wrapper](../benchmarks/comparators/run_r_field.py) (reviewed line 107) replaces unmapped values with zero. Under `--rule all`, cells absent from native output are also clustered from zero profiles, so reported coverage can include cells the tool did not return.

Use a prespecified scoring universe, preserve missingness, and distinguish native outputs from wrapper-derived assignments. **[Unverified]** The effect on existing rankings requires rescoring. Preserve frozen results and disclose any corrected analysis separately.

### 5. P1: Superiority tests do not account for shared reference pools

[analyse.py](../benchmarks/heldout/analyse.py) (reviewed line 43) bootstraps individual lines and applies Wilcoxon tests while the protocol acknowledges shared pools and references. The `superior` flag nevertheless treats these as confirmatory results.

Holm correction handles multiple comparisons; it does not correct dependence between observations. Wilcoxon assumes independent paired differences. [SciPy documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.wilcoxon.html)

**[Inference]** Uncertainty may be understated. With only three pool compositions and reciprocal reference rotations, stronger evidence requires additional independent datasets or an explicitly justified dependence model.

### 6. P2: Specimens with zero usable SNPs crash allele inference

After coverage and bias filtering, a specimen can have no retained SNPs. `switch_probs` ([alleles.py](../python/haplodose/alleles.py), reviewed line 164) accesses `cm[0]`.

**Reproduced:** an empty SNP set raises `IndexError`. Continue with expression-only inference for that specimen and record why the allele channel was unavailable.

### 7. P2: Some output labels assert more than the model establishes

`joint_state` ([alleles.py](../python/haplodose/alleles.py), reviewed line 290) labels relative-neutral RNA dosage plus imbalance as `cnLOH`. That evidence does not itself establish absolute DNA copy neutrality or loss of a haplotype.

Likewise, bulk output metadata ([bulk.py](../python/haplodose/bulk.py), reviewed line 188) sets `arm_calls_calibrated=True` solely when the panel contains at least ten normals. Panel size alone does not demonstrate calibration.

Prefer explicit evidence labels, candidate genomic interpretations, and separately documented calibration results. Treat `normal_like` as similarity to the reference, not proof that a cell is biologically normal.

### 8. P2: GPU initialization errors are not handled consistently

`device.select` ([device.py](../python/haplodose/device.py), reviewed line 37) catches `RuntimeError`, but the installed JAX backend raises `AssertionError` when CUDA initialization fails. The existing clear-error test fails on this environment. Handle the supported backend failure paths without obscuring unrelated programming errors.

## Computational efficiency

The largest efficiency opportunities are structural:

| Location | Improvement |
|---|---|
| Depth pooling, `clones.py`, reviewed line 232 | Replace dense cell-by-stratum one-hot matrices with grouped accumulation. **[Estimate]** At one million cells and 40 groups, the float32 matrix plus a float64 copy require approximately **2.24 GiB**. |
| Shared profiles, `clones.py`, reviewed line 384 | Bound pairwise batches and avoid repeated unchanged comparisons. Forty tumour clones require **780 pair comparisons per update**. Invalidate caches whenever pools or nuisance parameters change. |
| Gene-level reporting, `pipeline.py`, reviewed line 343 | Pool genes for completed specimens in one streaming pass; avoid repeatedly reading enclosing blocks for interleaved specimens. |
| Bulk fitting, `bulk.py`, reviewed line 48 | Use truncated SVD for retained programmes and forward-only HMM calculations when ranking beta candidates needs only likelihood. |

## Direction for the manuscript and biological goal

Prioritize scientific calibration over further architectural expansion. A defensible manuscript contribution would demonstrate reliable clone discovery, quantified failure conditions, and bounded resource use against simple baselines. Pooling and haplotype-aware HMMs already have substantial precedent in [Numbat](https://pmc.ncbi.nlm.nih.gov/articles/PMC10289836/).

For the CRC/IGF2 question, establish tandem duplication and genome-wide gain burden from independent DNA measurements, then test their relationship to RNA-derived candidate clones. Use patients as biological replicates and prespecify ploidy, purity, batch, and other relevant covariates. Relative RNA composition cannot identify uniform genome doubling or establish why additional gains occurred.

## Verification completed

Checks ran on the existing compute allocation, without a new job submission:

- Existing Python suite: **75 passed, 1 failed, 1 skipped**, in 274 seconds. Failure: `tests/test_device.py::test_gpu_without_visible_gpu_fails_clearly`; skip: GPU comparison.
- Four targeted probes reproduced the depth-overlap, balanced-veto, mirrored-allele, and empty-SNP cases.
- **91 Python files** passed AST parsing.
- Real-cohort reruns, corrected benchmark rankings, fresh Rust build tests, and GPU accuracy validation: **not run**.

The Python suite was invoked through `pytest.main` with `-q -p no:cacheprovider`, a fresh scratch `--basetemp`, and `tests`, using `/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python`, `PYTHONPATH=python`, `JAX_PLATFORMS=cpu`, and affinity restricted to two allocated CPUs. Scratch directory: `/scratch/haplodose-review-20261102-cv457c3j`.

Observed environment: Python 3.11.18; NumPy 2.4.6; SciPy 1.17.1; JAX/JAXlib 0.10.2; scikit-learn 1.9.1; pytest 9.2.0.
