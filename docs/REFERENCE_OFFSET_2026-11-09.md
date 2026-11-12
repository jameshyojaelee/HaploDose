# A normal anchor is useful only for effects it shares with the query

Removing the free shared offset is not a recovery repair: all five original
cases fail, versus2/5 under current production. Independently labelled normal
counts can separate a shared technical effect from a common query effect under
the stated observation model. Mismatched normals can reverse a gain, and even
matched normals cannot distinguish a tumour transcriptional state with the same
RNA observation law. No package default is changed.

## Estimands and identifying information

The primary mechanistic target is the DNA-linked RNA response, a regional
natural-log effect0.6*log(3/2) in the gain conditions. It is not absolute DNA copy
number. Secondary targets are RNA abundance relative to a fixed baseline and
RNA abundance relative to calibration normals. These are different estimands.
Small error for the latter does not establish accurate DNA inference when the
normal/query relation or the DNA-to-RNA link is wrong.

The conditional model is independent NB RNA counts with known nominal exposures,
baseline, depth slopes and dispersion. Exposures need not equal realized count
row sums. Normal calibration cells are labelled by construction independently of
query inference. This is a known-model, binned-count experiment, not validation
of whole-workflow preprocessing, unknown library exposures or patient biology.
Reference transport and absence of a tumour-specific expression programme are
load-bearing assumptions. Independently corresponding DNA and normal identity
would distinguish their biological interpretations.

## Whole-workflow restriction fails

Task `dev/reference_offset_20261109` preserves the initial dirty source and uses
the prior validated wheel/environment read-only. `offset_restriction.py` changes
only shared-offset updates to zero. All original genes, references, nuisance
models, proposals, merging, confirmation, parameters and six criteria remain.
All77 recorded updates have zero previous offset; there is no partial fallback.
Job23096638 finishes all five cases, with0/5 original recovery passes. This
rejects zero offset as a sufficient correction before any new error/power claim.

The frozen inference metadata identify external reference cells. They do not
independently identify normal cells inside each query specimen. Scoring truth
is not used to manufacture such labels, and held-out reference observations are
not used to fit a specimen-specific anchor. An oracle normal-labelled workflow
is therefore not presented as an available-data package comparison.

## Independent normals: benefit, noise and failed transport

Job23097090 generates16 independent simulation draws under the frozen
`ANCHOR_DESIGN.md`:40 bins,60 calibration normals and240 query cells, keys
[20261109,903,replicate]. Conditions share observations where specified and are
dependent within each draw. Biological n=0. Known depth slopes enter all methods;
all40 bins are common, without smoothing or segmentation.

Two conditional NB profiles and two simple pooled-ratio baselines differ only
in using the fixed reference or subtracting the independent-normal profile.
Each pooled observed/expected pair receives one half-count. Fit functions accept
counts and known model covariates, never condition or truth labels. Estimating a
free RNA profile on normals and query separately, then taking their centred
difference, is the normal-anchored conditional contrast. It is not the full
HaploDose clone algorithm or a separately calibrated genomic caller.

Mean squared error against the DNA-linked RNA-response target:

| Condition | NB fixed baseline | NB normal anchor | Pooled fixed baseline | Pooled normal anchor |
|---|---:|---:|---:|---:|
| No effect |0.000912|0.004847|0.001057|0.005326|
| Gain, transportable baseline |0.001074|0.006140|0.001232|0.006365|
| Gain plus shared technical effect |0.029589|0.006118|0.029562|0.007080|
| Gain with mismatched normals |0.001074|0.032521|0.001232|0.032997|
| Tumour state, no DNA gain |0.014398|0.016773|0.014348|0.017613|

These are descriptive paired simulation means, not patient effect estimates or
superiority tests. NB anchoring improves the shared-effect case beyond the
matched pooled anchor but adds normal-sampling noise in the no-offset cases.
Neither result justifies an unconditional normal-anchor default.

In the shared-effect condition, NB anchoring gives mean regional contrast0.240379
against the true DNA-linked response of about0.243279. With mismatched normals
the mean becomes−0.110272, and the contrast is negative in all16 draws. The
normal-anchored own-relative-RNA MSE remains0.005879 there: the result can describe
RNA relative to the wrong normals reasonably while misrepresenting the genomic
response. This is a target/transport failure, not numerical failure.

Gain and state conditions have exactly identical complete observation arrays
and fitted profiles, with different DNA-linked targets. The NB anchored state
contrast is positive in all16 draws, averaging0.242292 despite no DNA gain.
No RNA optimizer or normal-based uncertainty calculation distinguishes those
observationally equivalent biological explanations. This complements the prior
[shared-offset identification counterexample](OFFSET_OBJECTIVE_2026-11-09.md).

## Independent audit and reproducibility

Job23097562 passes an independent audit: all30 original recovery criteria,
all77 zero-offset updates, exact replay of960000 count entries, input hashes,
320 case/method metric records and20 summaries. Scalar pooled arithmetic uses
integer count sums and separate normalization. Thirty multi-bin gradient checks
compare refined finite differences of full SciPy NB logpmfs with the analytic
gradient. All160 normal/query fits have saved stationarity diagnostics; the
panel reports no optimizer failures. Complete gain/state count and fit equality
is checked, without treating those duplicate conditions as replication.

All three CPU jobs finish without execution errors, using innovation and nslab.
The recovery process takes72.11s with1323884KiB peak RSS; count panel7.86s with
198376KiB; audit7.55s with197844KiB. Job headers, actual allocation, backend,
runtime, memory and dependencies are recorded in `jobs.json`,
`accounting_final.tsv` and each `analysis_JOB/execution.json`. The wheel is pinned
to SHA256 `711325c2020268d339822afc896efbd36c1f1370751317e32368e7eabb55b93f`.
No build, new dependency, shared installation or existing package/test change.

## Decision

Reject zero offset as a full recovery fix. Retain the normal-anchor estimator as
research code with an explicit transport assumption; do not infer a query-normal
anchor from scoring truth or choose normals by DNA agreement. The added normal
sampling noise motivates testing normal-only shrinkage, with every selection
step restricted to normals and its target kept separate from DNA validity.
That is a future experiment, not a demonstrated correction. Any candidate must
still pass the original whole-workflow recovery and false-split/power gates.
Those earlier failures, frozen outputs and the5% specimen criterion remain.
