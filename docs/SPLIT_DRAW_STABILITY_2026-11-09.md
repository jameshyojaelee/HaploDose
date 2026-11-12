# Count-splitting redraw stability: separate a law change from one realized split

This is a development diagnostic on the five existing raw-count matrices.
The known-gene-size research path changes both the splitting law and the realized
count halves. Its canonical recovery improvement therefore does not, by itself,
identify the law as the causal improvement. The experiment asks whether recovery
persists under other count-splitting draws with all other random settings fixed.
It does not select a best seed or create a new clinical error-control panel.

The target is the probability, over technical count redraws conditional on each
fixed input, of satisfying all six original recovery criteria. Biological n=0.
The alternative to a robust method improvement is a favorable canonical draw.
Repeated technical redraws can distinguish these in the fixed inputs; they do
not distinguish genomic clones from transcriptional state in biological samples.

## Frozen design

`dev/split_draw_stability_20261109/DESIGN.md` and `seed_bank.json` were written
before the primary outcomes. Ten effective count-split seeds come from the first
four BLAKE2b bytes of `HaploDose|split-draw-probe|20261109|rep`, interpreted little
endian. They are distinct and exclude the canonical seed20261031.

`CloneParams.seed` remains20261031 for preprocessing, PCA, Leiden and search.
Only `Rows.set_thinning` and its corresponding nested gene-count allocation
receive the alternate effective seed. Both requested and effective seed are
recorded explicitly. Ordinary bin thinning and known-size10 gene-native thinning
both use hard confirmation and merge stages. No margin weighting, aligned gene
veto, new threshold or true clone label enters fitting.

All5 source cases, both methods and all10 draws give100 complete public workflows,
each including reference preparation and all3 query specimens. These are300
specimen-condition fits on5 reused source matrices, not100 independent specimens.
Ten draws are a bounded first stability probe, not a biological power sample size.
All outcomes remain visible. Conditional exact intervals use the original six-part
success event and Bonferroni95% simultaneous coverage over10 method/case
probabilities (.05/20 per tail). Wide intervals remain wide. The5% whole-specimen
criterion is not replaced by this redraw analysis.

## Canonical gate and failure journal

| Job | Role | Result |
|---|---|---|
|23123022|Initial canonical gate|Failed after the5 unwrapped production workflows: QC contains preprocessing fields beyond `CloneParams`, so the metadata equality assertion was incomplete.|
|23123023|Dependent initial audit|Cancelled before execution; its dependency could not succeed. It ran0s.|
|23123140|Corrected canonical gate|Passed20 workflows: wrapped versus unwrapped, both methods, all5 cases.|
|23123141|Independent gate audit|Passed all120 original criterion checks, canonical fold conservation/identity and fixed features/inputs.|

The corrected check requires every declared model parameter and exact agreement
with the frozen complete preprocessing configuration. No model or scientific
criterion changed. The failed output/source remains in `gate_v1`/`source_v1`;
the corrected evidence is `gate_v2`, `gate_audit_v2` and `source_v2`.
All public tables and posterior arrays match in both methods; ordinary production
also matches its frozen prior replay. The gate independently checks30 bin-fold
models and30 gene-fold models.

## Primary batch and independently audited results

Ten CPU jobs23123398–23123407 run the two methods in five blocks of two draws,
4 CPUs/8GB/3h each, alternating innovation/nslab. Audit23123408 depends on all
ten and requests1 CPU/8GB/2h. No arrays, interactive QOS or shared-install changes.
The wheel is the previously validated403-Python/4-Rust-test artifact; production
source is unchanged. Commands and resource/status records are in `jobs.json` and
`analysis_JOB/`. All eleven jobs completed with exit0:0. No primary job failed.

The independent audit recounted every original criterion, checked source/input,
feature and initial-dispersion identities, conserved full/proposal/test counts,
verified effective seeds and known-gene aggregation, and independently inverted
binomial tails. All600 criterion checks,300 bin-fold models and150 gene-fold
models passed. This means the calculations agree with the unchanged criteria;
it does not mean all recovery outcomes passed those criteria.

| Source case | Ordinary bin: pass/draws | Simultaneous interval | Known-size gene: pass/draws | Simultaneous interval | Bin-only / gene-only successes |
|---|---:|---|---:|---|---:|
|20261031|9/10|[0.4194, 0.9997]|4/10|[0.0635, 0.8325]|5 / 0|
|20261107|0/10|[0, 0.4507]|8/10|[0.3193, 0.9924]|0 / 8|
|20261108|4/10|[0.0635, 0.8325]|9/10|[0.4194, 0.9997]|0 / 5|
|20261109|0/10|[0, 0.4507]|0/10|[0, 0.4507]|0 / 0|
|20261110|1/10|[0.0003, 0.5806]|0/10|[0, 0.4507]|1 / 0|

Intervals are rounded for display; full precision, per-draw failures and paired
outcomes are in `primary_audit_v1/report.json`. They are conditional on these
fixed source matrices, with simultaneous coverage over the ten declared
method/case probabilities. They are not biological specimen error bounds.

[Observed pattern, not a guarantee] The gene-native law helps two cases but harms
the original successful case; neither law recovers case20261109 in these draws.
The known-gene law is not a general rescue. A successful alternative draw falsifies deterministic
information absence in that fixed measurement. It does not authorize choosing
that draw retrospectively or promoting a multi-draw search without recalibration.

No new method was promoted. Production, all scientific criteria and frozen
historical outputs remain unchanged. This experiment used the hard-stage
gene-native path; it is not a repeat of the separate weighted-stage/effective-veto
research variant that reached canonical5/5. There is no identified deployable
gene-dispersion estimator or calibrated multi-draw procedure in this experiment.

The canonical three production failures, failed state-null error bound, uncertainty
failures and missing independent biological truth remain unchanged.

## Resource and source preservation record

All15 submitted job records remain in `jobs.json`:13 completed,1 initial gate
failed,1 dependent audit cancelled before execution. The ten primary jobs plus
final audit are11/11 completed,0 running,0 pending,0 failed. Primary peak recorded
batch RSS was1393160KiB (about1.33GiB); the largest primary elapsed time was2m26s.
Step RSS is not a measurement of simultaneous total allocation usage.
`accounting_final.tsv` preserves the scheduler rows. Both immutable source
snapshots retain all730 recorded file hashes unchanged; see
`source_integrity_checkpoint.json`. Documentation was updated outside these
snapshots. No frozen output or shared installation was edited.
