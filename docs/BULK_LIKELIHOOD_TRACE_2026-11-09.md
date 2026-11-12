# Bulk loss evidence changes before segment averaging

All13 implicated patient fits reproduce exactly from integer counts. Nine of10
raw-only loss windows are locally loss-favoured but cease to be so after HMM
context. Context also strengthens the15 posterior-only detections. This identifies
a real tradeoff, not a numerical defect or evidence that removing context is a
valid correction. Frozen results and scientific criteria remain unchanged.

## Scientific contract

The estimand is regional natural-log RNA abundance relative to the fixed normal
reference, represented by a conditional four-state RNA-HMM posterior mean. The
observation model is pooled Poisson-lognormal with estimated reference variance,
programme effects, centre and slope. A local posterior uses the identical
emissions and prior with dependence removed; its probabilities are not calibrated
DNA-state probabilities. Programme estimates and other nuisance parameters still
depend on HMM fitting, so a local-stage ablation is not a complete alternative.

Interpretation as DNA loss requires corresponding independent DNA at the declared
resolution and sufficient sensitivity, plus a sufficiently controlled RNA
reference and biological state. Regional transcription/reference mismatch is the
strongest competing explanation. Existing ABSOLUTE segments are uncertain
GRCh37 reference measurements; RNA hg38 positions are gene anchored, not lifted
over. These49 windows in13 exposed patients are mechanism investigations, not49
independent specimens or prospective validation. No new biological replication.

## Measurement chain and computational checks

The source raw-STAR panel has48 tumour patients and20 normal donors, with the
existing24/24 calibration/comparison split. Keep every callable comparison loss
window:19 detected by both raw ratio and posterior,10 raw-only,15 posterior-only,
and5 missed by both. Saved genes reconstruct exact counts/library totals; identical
reference and defaults reproduce posteriors, effects, Viterbi paths, centres,
slopes, variance and raw-ratio arrays. Recover original genes, bins, windows,
truth flags and support through the frozen scorer. Gene-weighted medians match
the scorer; count totals include each bin once.

Job23060452 completes in25s, process peak853392KiB. The independently written
audit23060883 completes in6s and verifies14703 scalar integer bin sums. On222
distinct patient/bin pairs it evaluates1776 infinite-domain SciPy integrals,
using Brent roots and two quadrature refinement targets rather than production
Newton/Hermite quadrature. Actual float32 count/log-mean/sd inputs are matched.
Maximum state log-Bayes-factor discrepancy is1.15070e-5; maximum refinement
change5.32907e-15, estimated log integration error2.74233e-9 and local probability
discrepancy1.14558e-6. Quadrature error estimates are not formal bounds.

Independent float64 block-wise HMM recursions differ by at most1.21801e-6 in
probability, satisfying the unchanged existing1e-4 HMM criterion. Failed audit
23060869 is retained: its independent implementation initially omitted declared
genomic-block resets. Correcting that audit restores the same chain definition;
the package and tolerance did not change.

## Mechanism, identification and falsification

| Existing call category | Windows | Local median loss probability | Smoothed median loss probability |
|---|---:|---:|---:|
| Both |19|0.729667|0.999979|
| Raw-only |10|0.722132|0.051097|
| Posterior-only |15|0.406695|0.999164|
| Neither |5|0.237053|0.094489|

Independent quadrature confirms all10 raw-only windows have gene-weighted median
local loss mass above one half;9 fall below after context. Of15 posterior-only
windows,5 have local loss mass above half and all15 do after context. One half
only describes mass on loss versus all other RNA states; it is not a calling
threshold. Numerical discrepancies are far smaller than these changes.

Two of the5 neither-detected windows have fitted loss-state amplitudes shallower
than the fixed posterior-score threshold. Every posterior effect is a convex
combination of the state effects, so that threshold is unreachable for those
fits. None of the10 raw-only windows has this constraint. Thus state amplitude
and HMM context are distinct mechanisms; neither justifies retrospectively
altering thresholds to match DNA. The independently reserved biological evidence
needed to distinguish CNA from transcription remains absent.

Exact replay failure would falsify the computational trace. A claimed benefit
from removing context must survive matched whole-panel false-positive, recall
and precision assessment, followed by a complete-workflow study. A gain confined
to selected loss windows is insufficient. Neither read inspection nor broader
reprocessing is triggered: there is no observed count/identity inconsistency.

## Whole-panel conditional ablation: do not promote

Job23061233 reconstructs all48 tumour fits exactly before loading DNA values.
The final local score and HMM score use identical nuisance fits, genes, support
and window medians. Each retains the original patient-equal neutral calibration
rule on the24 calibration patients; independent exact Fraction ranks agree.
The local-score threshold is-0.191632435, versus-0.356426865 for the differently
scaled posterior score. All historical raw/HMM patient scores reproduce exactly.

| Exposed comparison patients | Mean loss recall (13 patients) | Mean neutral FPR (24 patients) | Mean loss precision |
|---|---:|---:|---:|
| Raw ratio |0.545971|0.050508|0.185764|
| HMM posterior |0.500183|0.057043|0.411818|
| Fixed-fit local ablation |0.632280|0.063186|0.223608|

These are patient-average descriptive effects, not independent significance,
calibration guarantees or superiority evidence. Precision uses the existing
confident gain/neutral negatives; missing support is not silently treated as
neutral. The increased recall trades against worse FPR and substantially reduced
precision relative to HMM posteriors. **Do not promote removal of context.**
Nuisance estimates remain HMM-derived, so no conclusion about an independently
fitted local model follows. No further physical-distance or richer model default
is justified by this conditional ablation alone. Retain simple baseline evidence.

Immutable designs, snapshots, counts-to-bin links, per-patient fitted arrays,
hashes, environment/backend identities and resource records are under
`dev/bulk_likelihood_trace_20261109`. Scripts are
`benchmarks/hardening/bulk_likelihood_trace.py`, `bulk_likelihood_audit.py` and
`bulk_local_ablation.py`. Production code, original tests and frozen artifacts
are unchanged. Original thinning recovery2/5, false-split/power failures,
uncertainty failures and missing biological null panel remain unresolved.
