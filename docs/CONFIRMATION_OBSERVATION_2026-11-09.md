# Confirmation pools are correct; removing nuisance or extra tau does not rescue them

All27 captured binary chr12 confirmation comparisons pass independent count,
exposure, depth-moment, nuisance, variance and HMM checks. Ten complete workflows
retain exact public-output parity and all60 original recovery criteria. The
production recovery result remains2/5. No production correction is justified by
these checks; the remaining failure concerns what the working model supports,
not a demonstrated pool-construction or HMM arithmetic defect.

## Scientific contract

The estimand is the conditional log RNA-rate difference between the two assigned
cell pools. It is not the difference between pure genomic clones. The executed
observation model uses smoothed pooled log rates, a conditional Poisson sampling
term, estimated residual tau, depth/programme shifts and a three-state Gaussian
HMM with RNA effects minus/zero/plus beta*log(1.5).

The identifying assumptions are that the nuisance model preserves the relevant
contrast and this variance/state family describes the selected pools. The
strongest alternatives are mixed/selected expression groups, regional cell state,
reference mismatch or an inappropriate state family. A numerical calculation
cannot distinguish genomic dosage from an observationally equivalent RNA state.
Independent corresponding DNA and suitable biological/state replication remain
necessary. Biological n=0; these are exposed regression mechanisms.

The [frozen design](../dev/confirmation_observation_20261109/DESIGN.md) predates
execution. Both ordinary bin and known-size10 gene paths use hard confirmation
and merging, the same validated wheel, all five canonical matrices, and unchanged
posterior, core-length, group-size and false-split criteria. No truth is opened
by the observer. All actual-chr12 binary candidate comparisons are retained,
including unsuccessful whole-chromosome and second-half alternatives.

## What was measured and independently reconstructed

The observer copies the actual difference-function inputs and intermediates:
Ya/La/XLa, Yb/Lb/XLb, base, gamma, programmes, applied shift, centering, sampling
variance, tau, emissions, transition/initial probabilities, posteriors, paths,
reset/mask, model parameters and returned segments. It does not replace a fit.
The source frame is identified before existing experimental wrappers are entered.

The independent audit:

- Joins exact specimen/global rows to the previously frozen canonical count
  halves and saved proposal labels. Scalar integer sums reproduce all27 pooled
  count comparisons exactly; exposures and depth moments agree. No inferred
  label purity or regenerated count matrix substitutes for those observations.
- Reconstructs rates and conditional sampling variance directly. It independently
  computes the depth formula and the programme ridge prediction with augmented
  least squares, rather than the production normal-equation fit. Guards include
  arithmetic magnitude, operation count and matrix conditioning.
- Reconstructs weighted-median centering and the robust first-difference tau,
  including the existing floor. It independently reconstructs credible cores from
  the saved posterior/path under the unchanged every-bin/minimum-length rule.
- Runs separate NumPy long-double forward/backward and Viterbi recursions,
  splitting at every reset. A243-path exhaustive fixture agrees to3.33e-16;
  float64/extended-precision comparisons and operation/magnitude guards pass.
  JAX posterior discrepancies are recorded for every comparison.
- Independently recounts all60 original scientific criteria from cell identity
  and difference tables. Computational agreement includes the original failures.

The float32 guard allows both HMM passes, normalization, logsumexp and exponentiation
over each reset block. It is deliberately conservative and is a checked numerical
error estimate, not a formal libm enclosure or statistical confidence interval.
Some other windows overlap this broad guard around0.95 and remain numerically
unresolved for interpretation. The rejected first-half production109/110 cases
are far from that boundary; numerical precision cannot account for their loss.

## Prespecified conditional ablations

These operate only on saved pools, not through a new complete pipeline. Depth-only
and uncorrected tracks retain the executed sampling variance and tau, and use the
same weighted-median centering. The tau-free ablation retains the executed
contrast and removes only extra residual variance. All retain the state means,
transition/initial probabilities,0.95 posterior criterion and minimum core length.
Removing tau is deliberately not asserted to be a valid uncertainty correction.

For the first-half chr12 candidates, the table gives the **largest nonzero-state
probability among bins starting before10Mb**. This is a diagnostic bound on support,
not a replacement for the required contiguous credible core. The territory retains
boundary-crossing bins rather than silently trimming unfavorable observations.

| Production source case | Executed | Remove programme shift | Remove depth and programme shifts | Remove extra tau | Required chr12 core |
|---|---:|---:|---:|---:|---|
|20261108|0.997886|0.993188|0.999751|0.999955|Present in all four:bins212–216|
|20261109|0.008886|0.003628|0.004381|0.015703|Absent in all four|
|20261110|0.000878|0.000434|0.001027|0.001229|Absent in all four|

These are conditional, dependent calculations on one realized pool comparison
per case. They are not error-rate or power estimates. The 20261108 core is not a
new recovery: its accepted distinction is still lost later in the full workflow.

The known-size gene path is also retained completely. For20261110 its first-half
maximum rises from0.757297 to0.926667 when tau is removed, still without a core.
Removing programme or both nuisance shifts likewise gives no core. Gene109 can
pass this initial stage but fails later; some of its conditional boundary decisions
overlap the conservative numerical guard. No canonical gene-path result is promoted
as an independently calibrated method.

## Disposition and remaining hypotheses

1. Pool identity/count arithmetic and the HMM recursion pass their independent
   checks. This rejects a computation defect of those kinds as the explanation
   for these observed production comparisons.
2. Removing programme correction, both nuisance corrections or extra tau is
   insufficient for production109/110's first-half comparisons. Do not disable
   these terms or launch a large null/power study as though rescue had occurred.
3. Very low posterior support under a fixed three-state working model does **not**
   prove that the RNA measurement contains no clone information. The earlier
   known-parameter gene/bin oracle and the retained proposals establish different
   conditional questions. Mixed-pool attenuation, the state family and predictive
   evidence after selection remain distinct mechanisms to assess.
4. A finer state grid is not an automatic fix. The historical fractional-state
   experiment increased false splitting and remains a failed attempt. A replacement
   observation/statistic would need an explicit estimand, numerical reference,
   recovery benefit and the complete-workflow state-null/power experiment before
   promotion. Current false-split and biological-calibration claims remain failed
   or unestablished.

## Jobs, preservation and reproducibility

| Job | Role | Allocation | Scheduler result |
|---|---|---|---|
|23124265|Ten-workflow executed-payload capture|innovation CPU,4CPU/8GB/2h|Completed2m21s,exit0:0|
|23124266|Independent reconstruction/ablation audit v1|nslab CPU,1CPU/8GB/1h|Completed8s,exit0:0|
|23124319|Audit v2 adds scalar count-pool/exposure/moment reconstruction and core recount|nslab CPU,1CPU/8GB/1h|Completed8s,exit0:0|

Audit v2 is an additional independent check on the same trace; no fits were
repeated for it and no earlier audit was weakened. Source v1/v2 and both audits
remain intact. Both source snapshots preserve all741 recorded hashes unchanged.
`jobs.json`, `accounting_final.tsv`, `analysis_JOB/execution.json` and process
`resources.txt` retain exact commands, dependencies, imports, backend, exits and
memory. Scheduler peak trace RSS was1286796KiB; process records must also be read,
as scheduler samples can understate peaks and report0 for short audit jobs.
Process peaks were1622412KiB for capture and88684/101084KiB for audits v1/v2;
process elapsed times were136.84,5.41 and6.52s, respectively. None exceeded8GB.

The retained artifact still has403 Python/4 Rust test passes; the full production
suite was not rerun because only research observers/auditors changed. Frozen
scientific outputs, original tests, package defaults, dependencies and scientific
criteria are unchanged. All payloads and negative results are under
[confirmation_observation_20261109](../dev/confirmation_observation_20261109).
