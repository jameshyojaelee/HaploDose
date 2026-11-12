# Optional expression state improves null profiles but loses minority RNA

This is one fixed eligibility ablation of the preceding count-predictive
selector:three already fitted no-factor structures compete with its five
factor-enabled structures. No refits,ranks,starts,penalties or thresholds change.
Every preceding five-member score and prediction must reproduce exactly before
the new choices are locked. All290 cases and original32/64 criteria remain.
The estimands are expected RNA log odds and conditional state-zero log odds;
neither identifies DNA dosage or genomic membership. Biological n=0.

[Frozen scientific contract](../dev/predictive_state_presence_20261111/DESIGN.md),
[eligibility implementation](../dev/predictive_state_presence_20261111/study.py),
[independent reader](../dev/predictive_state_presence_20261111/reader.py),
[source inventory](../dev/predictive_state_presence_20261111/SNAPSHOT.json).

Both nslab CPU jobs complete0:0. Writer23221029 takes2m12.49s,71496KiB process
RSS;reader23221030 takes1m16.79s,96400KiB. All290 choices lock before truth.
All2320 independent full-PMF score comparisons,exact fold/order sums,numerical
choices,weighted points and every predecessor five-model replay pass.

| Condition,whole RNA MSE | Preceding selector | Optional state | Pseudobulk |
|---|---:|---:|---:|
| Homogeneous/300 |.00198067|.000971058|.000916010|
| Homogeneous/3000 |.000177963|.0000884867|.0000794514|
| Clean two-population/300 |.00348635|.00359691|.00413259|
| Clean two-population/3000 |.000202234|.000136896|.00348458|
| State plus populations/300 |.00652925|.00652925|.01619652|
| State plus populations/3000 |.000500701|.000500701|.01554463|

No-factor mean weight is.9583 in each homogeneous stratum,.5417/.9063 in the
two clean strata,andzero in both state-containing strata. The observed pattern
supports unnecessary factor estimation as a homogeneous variance cost. State
outputs remain unchanged. But low-depth minority RNA worsens,so this does not
repair rare-population information or qualify as a general correction.

The full gate fails29/32 original RNA and56/64 simultaneous tests,below the
preceding31/32+58/64. All32 high-depth comparisons remain passed. The specific
clean low-depth whole-RNA pseudobulk comparison remains passed43/48,simultaneous
lower win probability.6906934. Every failed simultaneous comparison remains:

| Target | Condition/depth | Endpoint | Comparator | Resolved wins/48 |
|---|---|---|---|---:|
| RNA |Clean/300|Minority|Global|34|
| RNA |Clean/300|Minority|Beta expected logit|31|
| RNA |State/300|Minority|Beta expected logit|33|
| State-free |Clean/300|Minority|Global|19|
| State-free |Clean/300|Minority|Beta expected logit|11|
| State-free |State/300|Minority|Pseudobulk|14|
| State-free |State/300|Minority|Global|14|
| State-free |State/300|Minority|Beta expected logit|18|

Reject general promotion. The next diagnostic must separate unknown-parameter
estimation from uncertainty about which cells carry a profile. Use actual observed
counts under the known generating distributions,without treating noiseless latent
rates as count-level power. A whole-risk Bayes action is not automatically an
optimal minority-only action,so oracle minority failure cannot establish a
universal impossibility bound. No new penalty or threshold search is justified.

All346 prior ledger values,91 protected files,public defaults,frozen results,
false-split/power and uncertainty criteria remain. Original production2/5 and
research5/5 recovery remain separate;this experiment changes neither.

[Complete independent report](../dev/predictive_state_presence_20261111/reader_23221030/results/REPORT.json),
[exact job commands and dependencies](../dev/predictive_state_presence_20261111/submissions.json).

## Bounded competitor inspection

The recorded field coverage remains168/181 for CopyKAT and167/181 for SCEVAN.
Each missing record is `incomplete_population_or_profile`. One inspected completed
CopyKAT run,r0_pool6,records7906 sent cells but7901 native matrix cells. Its
wrapper restricts query/holdout to native-present IDs under the recorded cluster
rule,leaving absent cells without profiles. This is not evidence of an OOM or
unfinished job. It does not establish the upstream filtering cause for all27
missing method/population records. No imputation,drop,reweighting or generic
competitor rerun follows;the fixed30-comparison family remains untestable.

[Recorded population coverage](../dev/native_stacked_profile_20261110/field_v2/population_coverage.tsv),
[inspected native-run QC](../results/kinker_copykat_20261101_2246/r0_pool6/copykat_20261101_235644/qc.json),
[native identity adapter](../benchmarks/comparators/run_r_field.py).
