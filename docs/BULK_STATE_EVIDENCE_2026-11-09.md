# Bulk runs retain per-bin RNA-model evidence before segment averaging

The historical locked bulk loss diagnostic could not separate HMM projection
from export averaging because its per-bin posteriors were not retained. New bulk
runs now preserve those actual fitted quantities, alongside all existing output.
No inference, reference, threshold or existing numeric result changes.

| Artifact | Axes | Meaning |
|---|---|---|
| `rna_state_posterior.npy` | samples, bins, four RNA states | Conditional state probabilities under the fitted RNA HMM |
| `rna_log_effect.npy` | samples, bins | Per-bin posterior mean RNA state effect, natural-log relative rate |
| `rna_state_index.npy` | samples, bins | Zero-based Viterbi RNA-state index |

Literal sample and bin order are in `samples.tsv` and `bins.tsv`. `run.json`
declares labels, dtypes, the working relative-copy states and the state-effect
formula `beta*log(state/2)`. Working states are not independently measured DNA
copy numbers. Probabilities condition on fitted nuisance quantities; they are
not validated biological state probabilities, interval coverage or specimen
error control. Each normal's reference fit excludes that sample, but feature
selection still uses the common all-normal universe and is explicitly recorded.
Median mode retains its declared reconstructed-pseudocount observation model.

This separates computational preservation, statistical interpretation and
biological validation. In particular, segment means can erase posterior variation
inside an unchanged Viterbi run. The new arrays expose that operation without
declaring that removing it improves DNA inference. RNA state and reference
composition remain competing explanations for apparent genomic signal.

## Exact preservation and numerical evidence

Two synthetic tests provide distinct per-bin posterior effects inside one neutral
run, for both sum and median aggregation. They require exact capture of arrays,
literal axes, dtype and provenance, while segment summaries retain their original
means. Every original test remains unchanged.

Isolated noneditable wheel validation23056837 completes in163s:374 Python tests
and4 Rust tests pass,1 GPU-only test skipped. Python process peak2512456KiB,
scheduler peak2264376KiB. Wheel SHA256 is
`2b285f46ed96ff1197952a1bec4246b1f7b4975ac3dc5c6d45505d36ad2fc365`;
native SHA256 remains
`e29d8e06561723441147c63b37782b9773367f175593770b3cb33e7a0a03e7a7`.
Dependency-copy23056835 completes in162s on an I/O node. No shared installation
or new production dependency changes.

Independent parity23057231 completes in27s, process peak511620KiB, exit0. Both
bulk modes use the same raw integer input and literal sample identities. The
original fit function is observed without changing its arguments or result.
Every exported state/effect/path value equals the old wheel's internal fit
exactly. All8 historical artifacts per mode are identical, including direct-gene
exports; gzip tables are compared after decompression. Original `run.json`
metadata is identical after removing only the declared additive export block.
All597 frozen baseline source files and all original tests remain unchanged.

Independent scalar RNA-effect calculations have maximum discrepancy3.09452e-8,
within the operation-derived float32 coefficient/dot-product allowances, whose
maximum is1.41614e-7. This budget accounts for stored log-state precision, beta
rounding and four-term accumulation; it is not an epsilon for changing a
scientific decision. Segment means and indices reproduce exactly in all40
captured segments. A varying posterior is retained rather than being replaced
by its segment average. No tolerance substitutes for historical artifact equality.

## Locked stage trace

Job23058249 uses only reference-refit000 of the existing exposed raw-STAR
COAD/READ panel:48 tumour patients and20 normal donors. The same genes, references,
model settings, original five comparators, callable support,24 calibration and24
comparison roles remain fixed. It adds the already fitted per-bin posterior stage
to descriptive scores, not a new tuned method. RNA fits and projections finish
before DNA values are used for scoring. The trace requires exact reconstruction
of all segment means, original DNA arrays/support and every new window median.
It completes exit0 in22s, process peak952648KiB. All2836 segment means/indices
and12480 scalar window medians reproduce exactly; original tumour IDs, DNA
arrays and callable support are identical. All five original comparators remain.

Comparison-role point summaries are:

| Stage | Mean Pearson,23 patients | Loss recall,13 patients | Achieved loss-side FPR |
|---|---:|---:|---:|
| Raw bin ratio |0.545858|0.545971|0.050508|
| Fixed-smoothed bin ratio |0.692597|0.545925|0.042894|
| Per-bin HMM posterior |0.720291|0.500183|0.057043|
| Segment-averaged posterior |0.716664|0.500183|0.056146|

All24 comparison patients remain in the population. One lacks an evaluable
correlation and11 lack a loss-recall denominator. Loss FPR is a patient-average
descriptive rate at thresholds learned from the original24 calibration patients;
it is not a certified5% population bound or the separate clone false-split rate.
No endpoint family, inferential p-value or historical153-test result changes.

Independent audit23058889 completes exit0 in3s, process peak80552KiB. Exact
Fraction-valued patient probability masses reproduce all three stage calibration
thresholds, separately from the scorer's integer-LCM search. It verifies288
patient recall/FPR values,138 extended/scalar Pearson references,138 independent
average-rank correlations and36 descriptive means. Constant/undefined cases stay
missing. Correlation arithmetic uses the existing operation-derived enclosure;
it is not a statistical confidence interval. Frozen masks and sample IDs agree.

Paired-window audit23059383 completes exit0 in3s, process peak80864KiB, repeating
the same checks and retaining the first audit. At their separately calibrated
thresholds, the posterior and segment stages call the same34 of49 callable loss
windows; both miss the same15. Thus equal recall here does not hide exchanged
loss-window identities. Neutral calls do differ:197 in common,18 posterior-only,
14 segment-only among3309 callable neutral windows. Windows are nested within
patients and do not add biological replication.

Raw ratio versus posterior calls overlap on19 loss windows;10 are raw-only,
15 posterior-only and5 called by neither. Posterior therefore detects more
pooled loss windows while having lower patient-average recall. This illustrates
why pooled window totals must not replace the prespecified patient estimand.
The next mechanistic question concerns those asymmetric changes and patient
heterogeneity, not simply the total number of detected windows.

Decision: keep per-bin evidence available, but reject removal of segment averaging
as a demonstrated bulk loss-recovery fix. The posterior's mean Pearson point
improves0.003627 over the averaged export while loss recall is unchanged and
achieved FPR increases0.000897. These descriptive changes do not establish a
benefit or equivalence. The same15 missed loss-window calls are present before
segment averaging at these thresholds. [Inference] Further investigation must
explain the10 raw-only and15 posterior-only loss calls, retaining patients as
the replication unit, rather than assuming export averaging caused all loss.
This does not establish the cause in transformed-input LUAD,
nor invalidate the separate Kinker case's larger averaging loss. The earlier
matched RNA-only simple baselines and their retained negative comparisons remain
mandatory; this new stage trace does not replace them or support superiority.

The new trace may distinguish HMM compression from segment averaging. It cannot
validate raw-count uncertainty using the historical transformed-RSEM LUAD input,
and does not do so. Existing uncertain DNA windows retain their resolution and
correspondence limits. No inferential test, superiority decision or biological
replication is added. The original clone recovery, false-split/power and routine
uncertainty failures remain unresolved.

Implementation is `python/haplodose/bulk.py`; new tests are
`tests/test_bulk_state_exports.py`. Diagnostic scripts are
`benchmarks/hardening/bulk_state_export_parity.py` and
`bulk_posterior_stage_trace.py`. Designs, immutable source snapshots, environments,
hashes, old/new captures and execution records are in
`dev/bulk_posterior_export_20261109`. Historical outputs are never retrofitted.
