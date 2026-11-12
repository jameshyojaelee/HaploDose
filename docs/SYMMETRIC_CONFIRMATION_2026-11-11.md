# Symmetric prediction from one frozen count split

The scientific goal remains independently supported genomic RNA-derived profiles,
CNA-defined populations and spatial assignments,with the original correctness,
error and power requirements. This experiment tests whether the arbitrary roles
of an unchanged count allocation discard useful membership information. It
changes neither the observations nor the biological criteria and does not turn
RNA confirmation into genomic identification.

Completed result,2026-11-11: all nine conditional mechanism probes pass. All five
required contrasts resolve positive with component mass at least30 in both
directions;all four state controls resolve negative. The independent reader
reproduces every decision. This permits a separately frozen full-workflow test;
it does not establish production recovery or specimen error control. Original
production recovery2/5,continuous-state false splits14/200 and correlated
selected-interval coverage54/1000 remain unchanged.

## Completed results and independent verification

| Case | Forward | Reverse | Symmetric | Numerical diagnostic | Decision |
|---|---:|---:|---:|---:|---|
| Original20261031 | -0.776271 | 3.228360 | 1.226044 | 0.000713 | Positive,supported |
| Original20261107 | 6.179664 | 2.919889 | 4.549777 | 0.000876 | Positive,supported |
| Original20261108 | 2.675432 | 1.641792 | 2.158612 | 0.000833 | Positive,supported |
| Original20261109 | 1.860061 | 0.554473 | 1.207267 | 0.000514 | Positive,supported |
| Original20261110 | 1.493710 | 2.486849 | 1.990280 | 0.000706 | Positive,supported |
| State202611096001 | 0.000000055 | -3.323827 | -1.661913 | 0.001818 | Negative |
| State202611096002 | -0.190350 | -0.936875 | -0.563613 | 0.000784 | Negative |
| State202611096004 | -0.815211 | -1.602175 | -1.208693 | 0.001452 | Negative |
| State202611096006 | -1.225231 | -0.372439 | -0.798835 | 0.001189 | Negative |

The independent reader23207992 completes0:0 on ne1dc4-016 in16s. It reconstructs
81 exact rational score-array sums,27 averaging calculations,nine reverse
exposure/initialization records and36 independent adaptive full-PMF Gaussian
integrals. Maximum probability-score discrepancy is1.46369e-13,within the
recorded adaptive/tail and rounding diagnostics. Counts,roles and original split
identities are conserved. It imports no fitting implementation.

All12 jobs complete0:0;peak process RSS299904KiB. All91 protected source/test/data
identities and all290 preceding claim-ledger values remain exact. Source and
reader snapshots preserve their repository-relative paths. No public source,
default,scientific threshold or frozen result changed. Original checkpoint
records remain unedited in the claim ledger.

[Exact results](../dev/symmetric_confirmation_20261111/summary_v1/REPORT.json),
[independent reader](../dev/symmetric_confirmation_20261111/reader_23207992/results/REPORT.json),
[all job accounting](../dev/symmetric_confirmation_20261111/accounting_final.txt),
[process resources](../dev/symmetric_confirmation_20261111/resources_final.json),
[preservation](../dev/symmetric_confirmation_20261111/preservation_final.json).

## Evidence motivating a different prediction rule

The earlier redraw diagnostic varied technical splits while holding original
source counts fixed. Neither ordinary nor known-size gene thinning was a general
rescue;several cases changed success under redraws. Selecting a favourable draw
is prohibited. The latest shared-noise/compositional screen remains7/9:one
required original contrast is negative and one state control unresolved. It uses
a fixed half-share binomial allocation. There is no reason to privilege one
half's role when the same fitting procedure can be applied to both.

[Earlier complete redraw evidence](SPLIT_DRAW_STABILITY_2026-11-09.md),
[latest unchanged conditional failures](COMPOSITIONAL_STATE_2026-11-11.md),
[frozen new design](../dev/symmetric_confirmation_20261111/DESIGN.md).

For D|Y~Binomial(Y,1/2) and C=Y-D,the assumed allocation law is exchangeable.
Define F as the C|D alternative-versus-null score with parameters fitted to D.
Define R by the same procedure with the roles reversed. The new primary is
`S=(F+R)/2`. Averaging removes the antisymmetric role component. Under exact
exchangeability its variance cannot exceed the variance of F. It does not add
data,average every possible allocation or create two independent replicates.
This is not a full Rao-Blackwell calculation or a specimen-level error guarantee.
Equal weights are fixed by role symmetry,not selected using recovery outcomes.

## Estimand,observation model and unresolved assumptions

The estimand is predictive RNA separation of discrete populations from the fixed
continuous-state competitor in nine exposed mechanism probes. The unchanged
working model has gene NB sampling,shared gene Gamma intensity across halves,
a shared Gaussian cell state,linear outside-event depth response,a two-component
alternative with depth-independent prevalence,and bounded compositional regional
response. Fitted reference/noise,selected cells/features and original proposal
history remain conditioned upon. This is not independent post-selection inference.

Binomial halves retain the full gene NB shape and share technical intensity.
The conditional model comparison uses full-count marginal likelihood minus
fitting-half marginal likelihood;the allocation mass cancels. Treating the two
halves as independent NB observations would use another law. Null and alternative
use the same share,depth and observation model within each direction.

Reverse fitting uses archived C,its own regional median initialization,and outside
C counts obtained by subtracting archived outside D from exact full counts.
The selected event-bin union is excluded in both directions. Standardize depth
and construct A/(A+B) from the fitting half. Do not use old q or target D in the
reverse fits. Three reverse fits lock before the target file opens. The forward
fits and all27 per-cell arrays remain byte-identical,hashed and exactly summed.

Gene dependence,noise/reference misspecification and non-Gaussian transcriptional
state compete with a genomic explanation. A transcriptional population with the
same RNA distribution as a CNA population is not distinguishable by this rule.
Independent technical measurements test the count law;corresponding-cell DNA
and appropriate biological nulls test genomic membership/error. Biological n=0
here;two dependent directions and nine reused cases are not new specimens.

## Fixed gates and numerical qualification

No redraw,new seed,share change,optimizer/model/basis adjustment or threshold
selection. Keep allfive original cases andfour selected state probes,all counts,
features,populations,IDs,original source expectations and gene shapes. Fits retain
24/48/96 quadrature,variance boundary,two null/seven alternative starts and the
same optimizer settings. Required true contrasts must have component
responsibility masses at least30 in BOTH directions and a positive resolved
symmetric score. Every state score must resolve negative. Allnine required;
zero or unresolved scores fail. Reverse-only success cannot rescue the primary.

The numerical diagnostic is at least half the sum of both directional diagnostics,
as well as the direct symmetric convergence/rounding calculation. Dependence does
not justify a sqrt(2) sampling-SE reduction. These are numerical diagnostics,not
biological confidence intervals. A failed primary rejects unchanged symmetrization
without direction weights,extra redraws,allocation shares or model/threshold
sweeps. A pass only permits a separately frozen full adaptive workflow plus fresh
false-split/power assessment,keeping the5% specimen upper-bound requirement.

Qualification23207973 completes0:0 in24s. All16 exact rational NB/binomial
conditional pairs,symmetric-score identities,role exchange,cancellation,identity
permutation and identical-model zero checks pass. The unchanged bounded-model
gradient/variance checks pass:maximum rational error3.0531e-16,independent adaptive
full-PMF error1.6653e-16. All27 archived forward score sums are exact. Allnine
reverse inputs conserve source/event/outside counts. No new fit or result is
inferred from numerical qualification alone.

[Numerical checks](../dev/symmetric_confirmation_20261111/prepared_v1/NUMERICS.json),
[locked inputs and exact forward records](../dev/symmetric_confirmation_20261111/prepared_v1/INPUTS_LOCKED.json),
[executed orchestration](../dev/symmetric_confirmation_20261111/study.py),
[unchanged observation/optimizer](../dev/symmetric_confirmation_20261111/frozen_model.py),
[unchanged shared-noise conditional likelihood](../dev/symmetric_confirmation_20261111/frozen_joint.py),
[independent rational fixtures](../dev/symmetric_confirmation_20261111/checks.py).

## Independent reader and execution

The independent reader imports no candidate fitting code. It recomputes exact
rational forward,reverse and symmetric sums;checks the ideal arithmetic mean
against per-cell floating evaluation;recounts allnine unchanged decisions;and
reconstructs every source/target/outside count,covariate and initialization.
Independent full NB PMF Gaussian integration probes reverse fits using the same
RNA-only cell-selection rule as the earlier validated reader. Reusing a verified
forward result is distinguished from rerunning it. All failed and unresolved
cases remain in the primary denominator.

[Independent reader](../dev/symmetric_confirmation_20261111/reader.py),
[reader source identities](../dev/symmetric_confirmation_20261111/reader_manifest.json),
[fit/qualification source identities](../dev/symmetric_confirmation_20261111/source_manifest.json).

All jobs use nslab account/QOS,CPU partition,one CPU,no arrays. Fits and reader
use the existing read-only validated review_depth_balance environment,neutral
fresh output directories,-I,one BLAS thread and source/resource/environment
records. No installation,public/default/test change,frozen-result write or edit
to historical ledger values. Relative-path source_v1 retains allseven initial
source identities without basename collisions.

| Job | Work | RAM/time | Checkpoint status |
|---|---|---|---|
|23207973|Numerical qualification and input preparation|8GB/1h|Complete0:0,24s|
|23207982|Original20261031 reverse fit|8GB/3h|Running|
|23207983|Original20261107 reverse fit|8GB/3h|Running|
|23207984|Original20261108 reverse fit|8GB/3h|Running|
|23207985|Original20261109 reverse fit|8GB/3h|Running|
|23207986|Original20261110 reverse fit|8GB/3h|Running|
|23207987|State202611096001 reverse fit|8GB/3h|Running|
|23207988|State202611096002 reverse fit|8GB/3h|Running|
|23207989|State202611096004 reverse fit|8GB/3h|Running|
|23207990|State202611096006 reverse fit|8GB/3h|Running|
|23207991|Allnine summary|4GB/1h|Dependent on all fits|
|23207992|Independent reader|8GB/1h|Dependent on summary|

[Initial working tree](../dev/symmetric_confirmation_20261111/initial_git_status.txt),
[289 previous ledger entries](../dev/symmetric_confirmation_20261111/claims_before.json),
[91 protected source/data/test identities](../dev/symmetric_confirmation_20261111/prior_source_identities.json),
[qualification submission](../dev/symmetric_confirmation_20261111/prepare_submission.json),
[all eleven fit/summary/reader commands and resources](../dev/symmetric_confirmation_20261111/case_chain_submission.json).

The complete genomic/clone/spatial/SOTA goal remains active and unachieved.
Passing arithmetic,an exposed conditional screen or repeated compute cannot
replace suitable independent DNA truth and biological replication.
