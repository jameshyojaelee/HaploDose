# Target-excluded prior estimation:real improvement,incomplete scientific gate

Using all119 eligible cells instead of60 improves whole-profile RNA MSE over
the preceding two-fold moment-Beta estimator in all four heterogeneous settings.
The new estimator nevertheless fails its full gate:26/32 comparisons pass.
Clean low-depth pseudobulk wins increase to30/48,but the simultaneous lower
strict-win probability is0.402995,below the required0.5. No default is promoted.

## Cause,estimand and implementation

The original point priors used one60-cell training fold for each target,while
the leave-one-cell-out pseudobulk baseline used119. That baseline is legitimate
and remains unchanged. The fixed candidate removes the target from integer
population sufficient statistics and uses every remaining cell. It is a unique
target-exclusion policy,not selection of a favourable fold realization.

Letm=119,N=library total,S=sum of other cells' regional counts,andQ=sum of their
squared counts. The factorial moment is(Q−S)/(mN(N−1));the distinct-cell moment
is(S²−Q)/(m(m−1)N²). Their difference estimates between-cell regional-probability
variance. Divide by the unregularised mean's Bernoulli variance,clip to[0,1],
and retain the original half-count mean and moment-Beta expected-logit action.
Integer arithmetic computes the combined variance numerator exactly;an explicit
int64 bound rejects unsupported allocation sizes. No variance threshold or
alternative prior strength is selected. A pre-execution fixture revision adds
interior-rho cases to the exact Fraction check;both prepared sources are retained.

The estimand is cell-specific regional RNA log odds relative to the specified
reference. The observation model uses binomial margins and a fitted Beta working
prior. It does not integrate uncertainty in estimated prior parameters or jointly
model compositional dependence. RNA state can produce the same profile as CNA.
Genomic-response error is reported separately;corresponding independent DNA
would be required for genomic validity. Biological n=0,48 base draws per scenario,
paired depths300/3000,120cells and24regions. These288 conditions are exposed
development regressions;no former evaluation cohort becomes prospective again.

All original1152 baseline records and the288 preceding Beta controls reproduce
exactly. All predictions lock before truth scoring. The32-way family,48-draw
denominators,exact one-sided95% simultaneous bounds,Bonferroni adjustment,tie/
numerical-failure rules and four comparators remain unchanged.

## Complete result

| Condition | Depth | Pseudobulk whole MSE | Prior two-fold Beta | Leave-one-cell-out Beta |
|---|---:|---:|---:|---:|
| Two populations |300|.004132589|.005169724|.004071228|
| Two populations |3000|.003484575|.001268426|.001153843|
| Two populations plus state |300|.016196525|.011178183|.010006824|
| Two populations plus state |3000|.015544627|.002313576|.002201251|

All six failed comparisons remain:

| Condition | Depth | Endpoint | Comparator | Wins/48 | Simultaneous lower win probability |
|---|---:|---|---|---:|---:|
| Two populations |300|Whole RNA|Pseudobulk|30|.402995|
| Two populations |300|Minority RNA|Prior Beta|26|.324923|
| Two populations |3000|Minority RNA|Individual|4|.010367|
| Two populations |3000|Minority RNA|Prior Beta|26|.324923|
| Two populations plus state |3000|Minority RNA|Individual|3|.004731|
| Two populations plus state |3000|Minority RNA|Prior Beta|29|.382954|

The high-depth minority losses expose retained shrinkage bias. More eligible
training data reduces prior-estimation error but is not a sufficient general
repair. A small positive mean against low-depth pseudobulk does not satisfy the
strict repeated-specimen criterion. Reject the unchanged estimator's general
performance claim;no fold,action,prior-strength or threshold sweep follows.

## Independent checks,resources and evidence

The independent reader verifies34560 naive target-exclusion calculations and
829440 point actions,reconstructs every whole/minority/genomic-response error,
and matches1440 preceding baseline/control rows. Maximum prior-rho discrepancy
is2.03457e-17. Exact integer binomial tails reproduce all32 decisions. Exact
Fraction fixtures,target perturbation,permutation and zero-variance limits pass.
The action is unchanged from the previously qualified expected-logit estimator;
these checks validate the changed exclusion/moment arithmetic,not calibration.

Study23212017 completes0:0 in13.42s,peak81632KiB;reader23212018 completes0:0 in
3.15s,peak84672KiB,on ne1dc6-006. Both request nslab/nslab CPU AMD-v6,1CPU4GB1h,
with original read-only Conda/runtime003,neutral directories,one thread and
before/after source hashes. No installation,production edit,frozen output,
scientific criterion or historical decision changes.

[Frozen design](../dev/leave_one_cell_prior_20261111/DESIGN.md),
[exact moments and action](../dev/leave_one_cell_prior_20261111/model.py),
[qualification](../dev/leave_one_cell_prior_20261111/run_23212017/results/QUALIFIED.json),
[complete scientific comparisons](../dev/leave_one_cell_prior_20261111/run_23212017/results/REPORT.json),
[independent result and all means](../dev/leave_one_cell_prior_20261111/reader_23212018/results/REPORT.json),
[all score records](../dev/leave_one_cell_prior_20261111/run_23212017/results/scores.tsv),
[exact submission commands](../dev/leave_one_cell_prior_20261111/SUBMISSIONS.json).

The full genomic-profile/clone/spatial/SOTA goal remains active. This improvement
does not repair the separate false-split,power,genomic membership or uncertainty
failures,and cannot replace suitable biological replication or DNA truth.
