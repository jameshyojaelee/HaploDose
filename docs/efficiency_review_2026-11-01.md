# HaploDose efficiency and algorithm review

Review date: 2026-11-01.

## Recommendation and scope

Keep Python as the public interface and orchestration layer. HaploDose already
uses Rust for sparse aggregation and compiled JAX for numerical kernels. The
largest opportunities identified in this review concern algorithms, data layout,
and mathematical implementation. Address these before considering a language
rewrite.

The review covered the current working tree, including bulk analysis, clone
fitting, alleles, HMMs, likelihoods, I/O, Rust kernels, tests, and benchmarks.
Existing uncommitted changes were included in the inspection and preserved. No
implementation files were changed. This document records the review and its
verification limits.

## Mathematical and numerical findings

### 1. Pooled-count likelihood loses accuracy at high counts

Source: [`pln_logpmf`](../python/haplodose/kernels/likelihood.py).

The function forces float32 and subtracts large, nearly equal terms. With the
observed count equal to the expected count and `tau=0.15`, controlled probes gave:

| Count | Current log probability | Float64 evaluation of the same quadrature |
|---:|---:|---:|
| 100,000 | -10.500 | -10.535 |
| 1,000,000 | -13.000 | -12.837 |
| 10,000,000 | **+16.000** | **-15.140** |

A positive log probability is impossible for this discrete probability mass
function. Independent, scaled SciPy integration confirmed the last reference
value: approximately `-15.13991642`.

**Improvement:** formulate emissions using stable Poisson deviance expressions,
remove state-independent constants where appropriate, and selectively use
float64. Merely enabling JAX float64 would not fix the explicit casts.

### 2. HMM message accumulation damages posterior accuracy

Source: [`_forward_scan`, `_backward_scan`, and `forward_backward`](../python/haplodose/kernels/hmm.py).

The implementation carries unnormalized forward/backward values across chromosome
restarts. Adding the same constant to every state's emission must leave
posteriors unchanged. Controlled tests violated this invariant:

| Positions | Common emission offset per position | Maximum posterior change | Maximum posterior row-sum error |
|---:|---:|---:|---:|
| 1,000 | -10 | 0.00127 | 0.00049 |
| 1,000 | -100 | 0.01051 | 0.00390 |
| 10,000 | -10 | 0.01804 | 0.00391 |
| 10,000 | -100 | **0.14471** | **0.03174** |

These probes used four states and chromosome restarts every 250 positions.

**Improvement:** normalize messages at each step, track likelihood scale
separately, and restart chromosome messages independently. The standard scan
algorithm is suitable; its numerical implementation needs strengthening.

### 3. Programme-score updates omit part of the normalization derivative

Source: [`_chunk_scores`](../python/haplodose/kernels/estep.py).

The code normalizes expected counts with a softmax, then calculates the gradient
as though bin means had independent log links. For one controlled example, the
implemented likelihood gradient was **-7.500**, while finite differences gave
**-14.041**. Independently deriving the softmax correction also gave **-14.041**.

**Improvement:** derive the gradient and curvature from the actual normalized
negative-binomial objective. The programme loadings need composition-weighted
centering in that derivative.

Also, `R @ log_rate` is a surrogate for the responsibility-weighted objective.
Document and evaluate that approximation.

### 4. Regional negative-binomial mixtures use a Poisson-style parameter update

Source: [`SpecimenModel.event_mixtures`](../python/haplodose/clones.py).

The rate multiplier is updated using a ratio of pooled observed and expected
counts. With heterogeneous expected counts, this generally does not maximize the
stated negative-binomial likelihood.

In a small example, the code's estimate was **0.251**, while numerical likelihood
maximization gave **0.511**. The likelihood derivative at the code's estimate was
clearly nonzero.

**Improvement:** solve the one-dimensional weighted NB score equation and stop on
convergence. The regional dispersion approximation also needs calibration:
summing unequal NB means does not generally produce an NB with size
`number_of_bins / median(phi)`.

[Inference] These approximations can affect which regions enter the five
highest-ranked split proposals. Their contribution to existing missed subclones
remains unmeasured.

### 5. Pooled depth correction approximates a nonlinear expectation

Source: [`SpecimenModel.expected`](../python/haplodose/clones.py).

The pooled model evaluates a normalized rate at the pooled mean depth covariate.
The cell model instead implies pooling each cell's normalized expected rate:

`average(softmax(depth effect))` generally differs from
`softmax(average depth effect)`.

A two-depth numerical example produced a **10.3% relative error in one bin's
expectation**.

**Improvement:** accumulate expected counts alongside observed counts in bounded
chunks, or quantify and constrain this approximation over realistic depth
distributions.

These are controlled numerical findings. **Their effects on biological
conclusions have not been measured.**

## Avoidable computation and memory use

| Area | Current inefficiency | Recommended improvement |
|---|---|---|
| [Depth-matched comparisons](../python/haplodose/clones.py), `matched_differences` | Builds a dense cells × twice-the-number-of-pairs weight matrix, then casts it to float64. Pooling performs work on many zeros. | Use sparse weights or direct accumulation over participating groups. Cache reusable depth-stratified summaries where possible. |
| [Agglomeration](../python/haplodose/clones.py), `agglomerate` | Recomputes every pair after each merge, including unchanged pairs. | Cache pair results and invalidate only comparisons involving the merged group. |
| [Bulk aggregation](../python/haplodose/bulk.py), `run_bulk` | Constructs a dense genes × bins one-hot matrix and multiplies through its zeros. | Use a sparse assignment matrix or grouped reduction with existing SciPy/Rust capabilities. |
| [Bulk fitting](../python/haplodose/bulk.py), `start_beta` and `fit_sample` | Performs **32 full segmentations per sample without programmes and 42 with programmes** under the defaults. Most calls discard Viterbi paths. | Separate likelihood-only, posterior, and decoding operations. Batch beta-grid candidates and samples. Decode final paths when needed. |
| [Regional mixtures](../python/haplodose/clones.py), `event_mixtures` | Always runs 200 iterations per region and repeatedly calculates gamma-function terms that cancel from comparisons. | Remove those constants, add convergence criteria, and batch region calculations. |
| [Programme learning](../python/haplodose/programmes.py), `learn_programmes` | Computes the entire reduced SVD, then retains roughly ten components. | Use a seeded truncated/randomized SVD, already available through installed dependencies. |
| [Proposal neighbours](../python/haplodose/clones.py), `leiden_groups` | Default 20-dimensional PCA coordinates trigger brute-force neighbour search. This repeats during refinement. | Measure graph reuse and approximate neighbour proposals, retaining downstream biological checks. |
| [Gene-level checks](../python/haplodose/pipeline.py), `gene_loader` | Random subsamples become many tiny HDF5 reads because only consecutive selected rows are combined. | Read larger storage-aligned blocks and select requested rows in memory. |
| [Allele input](../python/haplodose/alleles.py), `read_alleles` | Reads chunks but retains all filtered tables before concatenation and sparse-matrix construction. | Build compact arrays incrementally or partition by specimen. Input memory is not actually bounded by `chunksize`. |

### Scaling calculations

- **Pairwise pooling:** [Estimate] with one million cells **in one specimen** and
  40 compared groups, the dense weights plus one float64 copy alone require
  approximately **17.4 GiB**. There are `40 × 39 / 2 = 780` pairs and 1,560 weight
  columns. This estimate excludes count matrices and other working memory.
- **Agglomeration:** collapsing 40 groups to one can evaluate **10,660 pairs**.
  Caching unchanged comparisons reduces this to **1,521**, preserving the same
  merge decisions. These counts assume merging continues to one group.

### Measured regional-mixture optimization

For 1,200 cells and 660 bins, removing cancelling gamma-function terms reduced
median runtime from **2.97 to 2.08 seconds**, approximately **1.43×**. Candidate
labels and ordering were identical.

The comparison used two unprofiled runs per implementation, alternating their
execution order. Original timings were `3.026` and `2.905` seconds; timings after
removing constants were `2.070` and `2.081` seconds. The modified function existed
only in the diagnostic process; repository code was not edited.

This is a kernel measurement, **not an end-to-end speedup claim**.

## JAX compilation and execution boundaries

Source: [`cell_loglik`, `_chunks`, and `_pad`](../python/haplodose/kernels/estep.py),
and [batched HMM functions](../python/haplodose/kernels/hmm.py).

Current chunking only partly stabilizes shapes. Trailing outer chunks, changing
clone counts, and shrinking pair batches still produce new compiled shapes.
Controlled probes observed separate cached executables as these dimensions
changed. JAX documents this shape specialization behavior:
[JAX compilation documentation](https://docs.jax.dev/en/latest/201/jit.html).

**Improvement:** use fixed batch sizes with masks, sensible shape buckets,
reusable device arrays, and fewer NumPy/JAX round trips.

The neighbour-search concern was confirmed by the installed scikit-learn source
and its documented dimensionality heuristic:
[scikit-learn neighbour algorithms](https://scikit-learn.org/stable/modules/neighbors.html).

## Statistical design issues

### Even/odd bins do not remain an independent validation split

Sources: [`SpecimenModel.em`, `split_candidates`, and `try_splits`](../python/haplodose/clones.py).

EM assignments use all bins, and subsequent adaptive proposals reuse the test
bins. A posterior threshold of 0.95 does not establish a 5% false discovery rate
across those searches. Validate the complete selection procedure with calibrated
null simulations and independent evidence.

### Bulk normal calibration is only partly leave-one-out

Source: [`reference_model` and `run_bulk`](../python/haplodose/bulk.py).

The tested normal is excluded from the mean reference, while variance and
programme loadings still include it.

[Inference] This can make normal-based calibration optimistic, especially with
small panels. Cross-fit all reference components by biological sample.

## Python, Rust, and analysis resolution

1. **Keep Python** for interfaces, metadata, workflow control, and scientific
   iteration.
2. **Fix correctness and remove unnecessary work first.**
3. **Compare selected CPU kernels with Rust after profiling.** The small-state
   HMM and regional mixtures are reasonable candidates. JAX remains useful where
   large, stable batches justify compilation and accelerator execution.

Bulk analysis also warrants evaluating an inexpensive arm-level path. The
existing [benchmark notes](STATUS.md) report strong simple arm-median results and
a more specific benefit from segmentation for focal losses. That supports making
computational effort depend on the requested resolution. Those cohorts were not
independently rescored during this review.

## Verification and limits

**21 existing tests passed** using the following command in the existing compute
allocation:

```bash
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 \
    MKL_NUM_THREADS=2 RAYON_NUM_THREADS=2 JAX_PLATFORMS=cpu \
    taskset -c 4,5 \
    /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s \
    -m pytest -q -p no:cacheprovider \
    tests/test_hmm.py tests/test_likelihood.py \
    tests/test_programmes.py tests/test_native.py
```

The test run completed in 19.71 seconds. Passing these tests did not establish
the numerical invariants examined by the additional probes.

The numerical probes covered finite-difference gradients, scalar NB likelihood
maximization, nonlinear pooled expectations, HMM emission-offset invariance,
float64 and independent quadrature comparisons, JAX compiled-shape caches, bulk
segmentation call counts, and local aggregation/mixture comparisons. Probes ran
in temporary Python processes; they were not added as permanent regression
tests.

Environment inspected:

| Component | Version |
|---|---|
| Python | 3.11.18 |
| NumPy | 2.4.6 |
| SciPy | 1.17.1 |
| JAX | 0.10.2 |
| jaxlib | 0.10.2 |
| scikit-learn | 1.9.1 |
| HaploDose | 0.1.0 |

Full integration tests and cohort reruns were **not run**. No jobs were
submitted, no dependencies were installed, and no implementation changes were
made. Proposed changes that affect the statistical model require scientific
revalidation; local speed measurements do not establish end-to-end gains.
