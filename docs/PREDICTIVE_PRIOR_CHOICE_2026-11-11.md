# Predictive simplification of the regional empirical prior

The scientific goal remains accurate genomic profiles and supported clones,
with RNA estimation and genomic validation kept distinct. The119-cell empirical
prior preserves all16 minority comparisons but fails low-depth whole-profile
accuracy. More training information did not solve that tradeoff. Test one
structural simplification:allow the common regional rate to win when it predicts
eligible counts better than a fitted distribution.

For each target,remove it from ALL prior fitting and model selection. Retain
the original twofold identities among the remaining119 cells. Fit empirical
and pooled common-rate models in each direction,predict the other direction,
and sum both roles so every eligible cell is predicted once. Choose empirical
only when the gain is positive beyond its numerical diagnostic;otherwise use
the exact strong pseudobulk estimate. Both final predictors retain all119
eligible cells. No free weight,prior-strength grid,scientific cutoff,depth rule,
favourable fold or action change. Numerically unresolved choice contributes
the entire prediction difference to the error allowance.

The observation model is binomial regional count margins. Estimand:individual
regional RNA log odds against the same external reference. Exchangeability and
adequate prediction within this specimen are working assumptions. Regional
transcription/cell state can produce dosage-shaped RNA;choosing the better RNA
predictor cannot identify DNA or validate genomic membership. Independent DNA
and biological replication remain required. Biological n=0;exposed development.

Same288 cases,all120cells/24regions,48 draws/scenario,both depths,all original
baseline populations and32 simultaneous strict-win comparisons. Preserve
all2016 previous score rows and all old prediction arrays. All new RNA outputs
lock before truth scoring. A failed gate rejects the fixed criterion without
adding a penalty/cutoff or tuning the folds. Positive subgroup results cannot
replace the complete gate.

Independent evidence:exact rational predictive ratios,prior optimisation
references,target perturbation leaving its inner fits/selection unchanged,
every inner count/identity/support/pooled mean,complete-binomial-PMF predictions
for every eligible target/region,all829440 final empirical posterior actions,
exact historical rows/arrays,independent fsum endpoints and integer decision
tails. Arrays are cached once per case. Reuse fits only for exactly identical
training identity tuples and regions. No omitted scientific or numerical checks.

[Frozen design](../dev/predictive_prior_choice_20261111/DESIGN.md),
[selection implementation](../dev/predictive_prior_choice_20261111/model.py),
[study](../dev/predictive_prior_choice_20261111/study.py),
[independent reader](../dev/predictive_prior_choice_20261111/reader.py),
[source snapshot](../dev/predictive_prior_choice_20261111/SNAPSHOT.json).

Jobs23215251/23215252 use nslab/nslab CPU AMD-v6:fit4CPU8GB3h and afterok
reader1CPU4GB2h. Existing Conda/runtime003,four fixed case workers,one BLAS
thread each,neutral directories and source hashes/resources. No arrays,
dependencies,shared installations,public/default edits or frozen changes.
Before submission:AST,bash syntax,all91 protected identities and immutable
source snapshot passed. Results pending at this checkpoint.
[Exact commands](../dev/predictive_prior_choice_20261111/submissions.json).

## Completed result: mean RNA improvement, failed general gate

Fit23215251 and independent23215252 complete0:0 in4m05s/9m18s. All2016 previous
score rows and original prediction arrays reproduce exactly. The reader checks
98703360 inner predictive observations,all target-excluded inner supports/means,
all829440 final empirical posterior actions,independent fsum errors and all32
exact integer decisions. No numerical check or tolerance is removed.

| Condition | Depth | Pseudobulk whole RNA MSE | Empirical prior | Predictive choice |
|---|---:|---:|---:|---:|
| Two CNA populations | 300 | 0.004132589 | 0.005449287 | 0.004359660 |
| Two CNA populations | 3000 | 0.003484575 | 0.000688085 | 0.000593258 |
| Two CNA populations plus state | 300 | 0.016196525 | 0.011538307 | 0.010675680 |
| Two CNA populations plus state | 3000 | 0.015544627 | 0.002036112 | 0.001961058 |

[Observed pattern, not a guarantee] The fixed choice reduces mean whole RNA
error versus the empirical control in all four heterogeneous strata. It passes
29/32 original comparisons but fails three:

| Condition/depth | Endpoint/comparator | Wins/48 | Simultaneous lower strict-win bound |
|---|---|---:|---:|
| Two CNA populations/300 | Whole RNA/pseudobulk | 16/48 | 0.153698 |
| Two CNA populations/300 | Minority RNA/original Beta | 12/48 | 0.095912 |
| Two CNA populations plus state/300 | Minority RNA/original Beta | 14/48 | 0.123908 |

There is one numerical boundary overlap in the first comparison;the scientific
failure persists independently of that observation. The hard selection improves
whole error while sacrificing low-depth minority accuracy relative to Beta.
The full gate fails. Do not promote it,add a threshold or choose a favourable
depth/region. Point means are not a replacement for the simultaneous criteria.
[Full independent report](../dev/predictive_prior_choice_20261111/reader_23215252/results/REPORT.json),
[all decisions](../dev/predictive_prior_choice_20261111/run_23215251/results/comparisons.tsv).

## The state-free target changes the development direction

A separate compute-node reader independently reconstructs all2304 saved
state-free response errors,with counts/references/folds/hashes and fsum checks.
The legacy column name `genomic_response_mse` means the simulated regional RNA
log-odds response with the added state removed. It is not measured DNA dosage
or absolute copy-number truth. Biological n=0.

| Condition | Depth | Pseudobulk state-free whole MSE | Predictive-choice state-free whole MSE |
|---|---:|---:|---:|
| Two CNA populations | 300 | 0.004132589 | 0.004359660 |
| Two CNA populations | 3000 | 0.003484575 | 0.000593258 |
| Two CNA populations plus state | 300 | 0.004401166 | 0.010081169 |
| Two CNA populations plus state | 3000 | 0.003756438 | 0.010666733 |

In the state condition,lower RNA error does not correspond to better state-free
dosage response. At depth3000,the selector has RNA error.001961 versus pseudobulk
.015545,yet state-free error.010667 versus.003756. This is a direct falsification
of the interpretation that the RNA improvement establishes genomic progress.
State-free minority errors remain better than pseudobulk(.036082 versus.078806
at depth3000);the whole/minority tradeoff must both remain visible.

All40 prespecified matched state-free comparisons are descriptive,with48
independent base simulations per condition and paired depths. They do not add
a new favourable significance test or replace the failed32-way primary family.
These results reject promoting marginal RNA predictive selection as the route
to genomic profiles. Do not continue a marginal-prior/weight/strength sweep.
The next mechanistic candidate must separate shared continuous state from
discrete profiles using their joint count structure,with explicit identifying
assumptions and both RNA/state-free targets. Independent DNA remains essential.

Response-reader23216318 completes0:0 in17s on nslab/nslab CPU AMD-v6,
1CPU4GB1h,afterok23215252. It writes separate results and changes no earlier
source/output/decision. All five jobs in this turn complete0:0;the older
inefficient23213185 remains running and untouched,with no downstream work.

[Independent state-free response report](../dev/profile_response_reader_20261111/reader_23216318/results/REPORT.json),
[exact methods and estimands](../dev/profile_response_reader_20261111/DESIGN.md),
[response-reader source](../dev/profile_response_reader_20261111/reader.py),
[response-reader submission](../dev/profile_response_reader_20261111/submission.json),
[fit resources](../dev/predictive_prior_choice_20261111/run_23215251/resources.txt),
[predictive-reader resources](../dev/predictive_prior_choice_20261111/reader_23215252/resources.txt),
[response-reader resources](../dev/profile_response_reader_20261111/reader_23216318/resources.txt).

Research recovery5/5 remains established. Production2/5,state14/200 and
selected correlated coverage54/1000 remain. No production/default change,
scientific criterion weakening,frozen artifact modification,biological
validation or SOTA claim. Full goal active and unachieved.

[Final source/claim preservation checkpoint](../dev/profile_response_reader_20261111/preservation_checkpoint_v1.json),
[scheduler checkpoint](../dev/profile_response_reader_20261111/accounting_checkpoint.tsv).
