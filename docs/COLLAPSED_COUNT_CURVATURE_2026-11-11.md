# Remaining collapsed count fits: a no-truth objective test

The [shared-incumbent correction](SHARED_INCUMBENT_RECOVERY_2026-11-11.md)
resolves one fixed witness, with better minority RNA profiles and105/120 matched
assignments. It leaves two other low-depth fits collapsed. The next experiment
asks whether those fits miss an observable likelihood direction, while keeping
the genomic identification and full false-split/power criteria unchanged.

Use the original low-depth clean/state failures, preflight cases2/4, and the
original homogeneous/state-only controls, cases0/6. Source counts/models are
locked before new computation. First apply the qualified optimizer unchanged.
No outcome labels, latent rates, generating parameters or truth functions enter
this investigation. Biological n=0; these are exposed mechanism cases.

For a fitted single-population Gaussian-state count density, form a symmetric
two-component perturbation of its mean. Its log-likelihood curvature uses
posterior count-score second moments minus conditional multinomial information.
This directly conditions on the fitted expression state; raw RNA covariance
can instead emphasize the state already explained by the model. Check the
formula on an exactly calculable two-bin example, complete-PMF finite
differences and96/192/384 integration. Finite-difference cancellation/truncation
diagnostics remain numerical checks, not confidence bounds.

Only the leading curvature direction and its existing exact three-means
projection are examined. If numerically positive, one symmetric ray is traced
with a curvature-derived scale and the existing SciPy scalar optimizer. All
failures are retained. Independently compare its best feasible point with the
current K2 point using full384-node multinomial PMFs and the same path prior.
A resolved better objective proves a remaining search miss. No witness does
not establish global optimality or an information limit. A control witness
also prevents interpreting the direction as biologically specific.

One job23228439 requests nslab CPU,1CPU8GB1h. Forty-six frozen input/source
identities are checked before/after. Static Python/shell syntax passes. Final
outcomes are pending at this checkpoint. There is no public change, new
production dependency, array, scientific threshold change, automatic full-study
expansion or frozen-result overwrite. Component support30 and all original
recovery,error/power,coverage and independent-DNA requirements remain.

[Complete prespecified design](../dev/collapsed_count_curvature_20261111/DESIGN.md),
[implementation](../dev/collapsed_count_curvature_20261111/study.py),
[source/input identities](../dev/collapsed_count_curvature_20261111/SNAPSHOT.json),
[submission](../dev/collapsed_count_curvature_20261111/submission.json).

## Completed result: no admissible optimizer addition

Job23228439 completes0:0 in38 scheduler seconds on nslab CPU. The exact
two-bin fixture gives matrix[[8,-8],[-8,8]] and anchored-direction curvature2.
Independent complete-PMF finite differences and96/192/384 integration checks
pass. All46 before/after identities match. No truth is opened.

| Original case | Projected log-likelihood curvature | Ray gain in unchanged full MAP target | Numerical allowance |
|---|---:|---:|---:|
|0 homogeneous control|231.26166|-14.865411|0.039868|
|2 low-depth clean failure|296.65974|-12.899175|0.043422|
|4 low-depth state failure|287.91201|-12.382476|0.075024|
|6 state-only control|333.16733|-10.811538|0.071694|

All four scalar searches terminate successfully. All candidate objectives are
independently worse beyond their original numerical budgets; no better feasible
same-objective point is established. **Reject adding or expanding this proposed
direction.** Keep the path prior and all scientific criteria. Do not reinterpret
likelihood curvature as a supported CNA distinction or add more starts merely
because the scalar count likelihood initially increases.

The finite-difference check is coarse: smallest-step derivatives remain2.58–5.16
below the matrix formula,and recorded truncation diagnostics are51.89–95.14.
It shows qualitative convergence,not a precision certificate for the Hessian.
The decisive rejection relies on independently evaluated complete-PMF objective
differences,whose recorded allowances are0.040–0.075 versus losses10.81–14.87.
That conclusion does not require treating the curvature calculation as exact.

The positive curvature appears in controls as well as true-population cases.
Its leading eigenvector,even after three-level projection,need not respect
genomic contiguity; the existing path prior penalizes its resulting profile.
This establishes failure of this particular unconstrained-curvature proposal,
not that all curvature approaches fail. It does not prove global optimality,
identify the prior as incorrect,or establish an information limit. A subsequent
proposal would need its own demonstrable improvement under the same complete
objective and controls. Neither a convenient new penalty nor a threshold on
this scalar curvature is justified by these outcomes.

All400 prior ledger values and91 protected package/native/data/test identities
remain exact; the ledger now has402 entries. Python and shell syntax pass. The
unchanged public suite is not rerun. Running0,pending0,done1,execution errors0
in this chain,with no downstream job or default change. Production recovery2/5,
state errors14/200,selected coverage54/1000 and the small-component support30
boundary remain unchanged. The full goal remains active and unachieved.

[Complete negative result](../dev/collapsed_count_curvature_20261111/run_23228439/results/REPORT.json),
[exact fixture](../dev/collapsed_count_curvature_20261111/run_23228439/results/FIXTURE.json),
[scheduler/resource accounting](../dev/collapsed_count_curvature_20261111/accounting_complete.tsv),
[process resources](../dev/collapsed_count_curvature_20261111/run_23228439/resources.txt),
[preservation checkpoint](../dev/collapsed_count_curvature_20261111/preservation_checkpoint_v1.json).
