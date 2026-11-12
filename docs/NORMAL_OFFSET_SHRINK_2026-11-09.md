# Normal-only shrinkage reduces noise but biases the shared-effect response

Learning anchor shrinkage from normals alone reduces the sampling-noise cost of
a fully fitted anchor in fresh conditional RNA simulations. It also leaves part
of a shared technical effect in the query estimate: regional response bias grows
even when average per-bin MSE falls. This is a conditional research benefit,
not a validated default, clone-recovery fix or genomic-error guarantee.

## Contract and frozen procedure

Task `dev/normal_offset_shrink_20261109` follows the
[normal-anchor failure analysis](REFERENCE_OFFSET_2026-11-09.md). The primary
estimand remains the DNA-linked RNA response in natural-log units. Secondary
RNA-relative targets are fixed by the actual baseline or full normal population;
they are not redefined to the selected fraction of the normal effect.

Observation model: independent NB binned RNA counts with known nominal exposures,
baseline, depth nuisance and dispersion. Nominal exposures need not equal count
row sums. There are16 fresh independent simulation draws, seed keys
[20261109,904,replicate],40 bins,60 normals and240 query cells. Five conditions
within each draw are paired and dependent. Biological n=0. Prior903 draws remain
development/regression evidence. This does not validate unknown library
exposures, preprocessing, learned nuisance terms or population DNA inference.

Normal calibration rows are split even/odd before counts are inspected. NB and
pooled predictors are fitted separately to each half. Each profile, multiplied
by t in[0,1], predicts the opposite half under the same native NB model. The sum
of heldout scores selects t, comparing endpoints and a scalar optimizer.
These bounds define interpolation from no correction to full correction; they
are not biological cutoffs. Two numerical precisions and independent refinement
are recorded. Neither query counts nor scoring truth enters this selection.

After selection, all60 normals are refitted and the selected t remains fixed.
Original half-fit predictions and heldout scores remain separate from this
descriptive refit; they are not overwritten as independent confirmation.
Copying a fraction chosen on half-sample fits to a full-sample fit may over-shrink
the less noisy full anchor. That possibility was prespecified, not added after
observing the shared-effect bias.

Six fixed comparators share all bins, known nuisance quantities and query fits:
NB and pooled profiles with no anchor, full anchor, or selected anchor shrinkage.
There is no smoothing, segmentation, adaptive clone search or new uncertainty
calibration. No comparator or failed condition is removed.

## Observed benefit and costs

Mean per-bin MSE against the DNA-linked RNA-response target, over16 draws:

| Condition | NB zero | NB anchor | NB shrink | Pooled zero | Pooled anchor | Pooled shrink |
|---|---:|---:|---:|---:|---:|---:|
| No effect |0.000860|0.005175|0.000900|0.000982|0.005663|0.001012|
| Gain |0.001103|0.005956|0.001110|0.001155|0.006650|0.001208|
| Gain plus shared effect |0.029123|0.005597|0.004924|0.028921|0.006623|0.005664|
| Mismatched normals |0.001103|0.030582|0.018416|0.001155|0.030960|0.017846|
| Tumour state, no DNA gain |0.014153|0.016157|0.014144|0.014131|0.016640|0.014011|

These are descriptive conditional simulation means, not biological-unit
superiority tests. Both NB and pooled methods benefit from normal-only shrinkage;
the idea's benefit is not exclusive to the NB optimizer. NB shrinkage has lower
MSE than pooled shrinkage under no effect, gain and shared effect; pooled shrinkage
has lower MSE under mismatched normals and the DNA-labelled state counterexample.
No universal ranking or package superiority follows.

In the shared-effect condition the true DNA-linked RNA response is approximately
0.243279. The mean regional contrast rises from0.249500 with full NB anchoring
to0.342675 with NB shrinkage. MSE improves because the variance reduction outweighs
the distributed error under this metric; the regional amplitude becomes less
accurate. A lower MSE alone therefore does not satisfy the scientific task.
The half-to-full refit mismatch is one candidate explanation [Inference]; an
intentional bias/variance tradeoff under normal prediction is another. Neither
has been isolated as the sole cause.

With mismatched normals, the NB-shrink mean contrast is−0.007130, positive in
only4/16 draws, despite a genuine gain. Normal predictive fit does not establish
normal-to-query transport. Gain/state count matrices and fitted profiles are
exactly identical. All16 state contrasts remain positive despite no DNA gain;
these are RNA contrasts, not a measured FPR of a genomic caller.

## Independent checks and retained numerical failure

Job23099596 freezes observations and truth separately. Job23099620 fits the
procedure;23099621 changes all query count values, and23099622 changes all truth
effects. In both challenges, all480 normal-profile arrays are bit-identical and
the complete normal-selection JSON is byte-identical to the original. All480
fitted vectors remain bit-identical under truth-only changes. The query challenge
does alter query estimates, and the truth challenge alters scoring, so these
are live information-leakage tests rather than unconsumed inputs.

The first independent audit23100365 fails because its coarse grid treats an
endpoint winner as settled without refining the adjacent interval. The corrected
auditor refines edge intervals too; the estimator and fitted results remain
unchanged. For example, null draw5 has selected t=0.01505245 and heldout gain
about0.0097464; the coarse archived reference chooses0 with gain0, while its finer
mesh resolves the interior optimum. Its original arithmetic bound is about
0.00000126, so enlarging the bound would have hidden a real search error.

Corrected audit23100566 reproduces the archived failures and passes all new
checks. It independently replays960000 count entries,240 scalar pooled normal
fits,480 metric records and30 summaries. It checks160 selection records,
representing64 distinct normal-input/predictor combinations, using full SciPy NB
logpmfs,33/65-point searches and local refinement. Duplicate normal inputs are
not counted as replication. Recomputing estimates with independently selected
fractions preserves all15 prespecified mean-comparison orderings. Gain/state
equality and all leakage challenges also pass. No numerical or scientific
tolerance was loosened.

## Exact descriptive decomposition

After observing the mean contrast shift, job23101499 performs a post hoc
diagnostic on the frozen estimates; it does not retune or select a method.
Exact rational arithmetic on the stored binary floats verifies1230 empirical
per-bin/regional identities:MSE=squared mean error+between-draw variation.
Between-draw variation includes the randomized baseline/exposures and counts;
it is not an isolated within-specimen sampling variance or a biological interval.

For the shared-effect NB comparison, per-bin squared mean error increases from
0.00111940 to0.00203284, while between-draw variation falls from0.00447769 to
0.00289093. Their sum explains the lower per-bin MSE. Regional-contrast MSE,
however, increases from0.000618764 for full anchoring to0.01071035 for shrinkage.
Both regional squared mean error and between-draw variation increase. The
pooled anchor shows the same direction:regional MSE0.000689229 becomes0.01276440.
This is direct evidence against promoting the method solely for per-bin MSE.
The new regional-loss analysis is explicitly descriptive and post hoc, not a
replacement primary endpoint chosen to declare a different winner.

## Scope, resources and disposition

All seven CPU jobs are terminal: six complete, one archived auditor failure.
The final preparation/fitting/leakage/audit/decomposition chain is6/6 complete. The main fit
takes13.13s with194656KiB peak process RSS; corrected audit6.07s with126520KiB.
Both innovation and nslab are used, each job1CPU/4GB/1h. Actual headers, runtime,
memory, node and dependency status are in `jobs.json`, `accounting_final.tsv`
and `analysis_JOB/execution.json` under the task directory.

Pinned read-only wheel SHA256 is
`711325c2020268d339822afc896efbd36c1f1370751317e32368e7eabb55b93f`.
`source_integrity.json` verifies all76 existing package/native/test/README/build
files unchanged. New research sources are `normal_offset_prepare.py`,
`normal_offset_shrink.py`, `normal_offset_audit.py` and
`normal_offset_error_decomposition.py` in `benchmarks/hardening`.
No dependency, shared installation, package default, existing test or frozen
result changes. New fitting code uses NumPy/SciPy on CPU; no GPU was needed.

Retain the conditional variance reduction as evidence, but withhold an
unqualified reference-correction default. The next decisive test should account
for the amount of information in the full normal refit, rather than copying a
half-sample shrinkage fraction, and assess regional amplitude as well as MSE.
Any gain must survive fresh conditional noise and the full-workflow recovery and
false-split/power gates. The original three recovery failures, correlated-state
splitting and absent biological calibration panel remain unresolved. Neither
normal-only tuning nor additional compute supplies missing reference identity
or independent DNA evidence.
