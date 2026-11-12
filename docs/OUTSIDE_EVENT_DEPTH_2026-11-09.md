# Matching-covariate and margin ablations do not rescue the remaining core

Excluding the candidate event from depth matching fails case109's unchanged
30-effective-cell criterion. Integrating membership while removing the extra
classification margin restores31 effective minority cells and raises maximum
nonzero support to.131252, but still produces no qualifying core. Reject all
three as sufficient recovery fixes. No full recovery or larger false-split/power
run follows. The original production three failures remain unresolved.

## Scientific contract and two successive designs

The estimand is a conditional regional relative RNA contrast between proposed
populations after matching nuisance depth. It is not DNA alteration or certified
clone identity. Matching on total RNA includes the proposed event's own RNA;
conditioning on it could remove genomic signal. Cell-state/depth association,
reference mismatch and altered RNA composition outside the event are competing
explanations. Independent corresponding DNA and complete-workflow null/power
evidence, not a successful conditional calculation, would establish validity.

The [first design](../dev/outside_event_depth_20261109/DESIGN.md) fixes two
covariates: log totals of ALL genes outside the candidate bins from full counts
or proposal counts only. Keep original max(count,1),full-model centering,eight
pooled quantile strata,weight-mass matching,q and abs(2q-1) base margins. The
proposal covariate avoids confirmation values in that covariate; it does not
make full nuisance estimation or adaptive selection independent. Outside genes
are not certified copy-neutral exposure. No truth chooses genes or strata.

Both fail109 because31 hard-labelled minority cells have22.895 effective cells
under the margins, below30. The [second design](../dev/outside_event_depth_20261109/MARGIN_DESIGN.md)
therefore tests a different mechanism: the cell-shared likelihood already
integrates q, so additionally applying a classification margin may discard
information twice. Use ORIGINAL full-count depth and the same eight strata,
but unit base weights before matching. This is specific to the generative
membership representation, not justified for old hard-label pooled tests.
q remains in every likelihood term; absent depth overlap still gets zero mass.

Hold the previously failed exact averaged beta-binomial LR evidence fixed.
Change only matching weights, not gene/expected-count arrays,reference noise,
component normalization,effects,priors,boundaries,every-bin.95/minimum5-bin
cores,30-effective-cell requirement or gene veto. There is no favourable redraw,
oracle membership or noise constant. Biological n=0. Transcriptional states and
genomic differences can still be observationally equivalent in RNA.

## Computational checks and exact-control limitation

[Implementation](../benchmarks/hardening/outside_event_depth.py) reconstructs
source/proposal CSR folds using global cell and gene IDs. On EVERY one of the
five actual matched comparisons it independently reproduces the old executed
weights EXACTLY from full count totals,q labels/margins and the original
coarsened matching. This checks the raw covariate and weight chain before
changing it. All67 identified initial comparisons,21 gene abstentions and
other original confirmation/merge stages remain recorded.

Independent Fraction stratum mass, scalar quantile interpolation, unit-weight
mass balance, cell permutation and confirmation-only perturbation of the
proposal-depth calculation pass. The previous independent Fraction/Decimal/
complete count-label-path kernel checks remain unchanged. Effective cells use
the same squared-weight formula and exact constant-weight support as production.

There is an additional failed control: maximum initial posterior replay error
is2.168404344971009e-19. It passes the pre-existing operation-derived numerical
consistency allowance, but it does NOT satisfy the first design's literal exact
initial-array wording. Preserve and report that failure; do not call it exact,
relax its criterion, replace a historical control or select another node.
The later margin design explicitly retains this failed gate. None of these
conditional results is a full-workflow recovery verdict or eligible promotion.

## All alternatives on the actual109 chr12 comparison

The source proposal is event_chrom12_bins209-217,502 positive-cell records,471
hard-majority/31 hard-minority cells. Original matched weight mass is22.710700.
Original single-draw maximum nonzero support is.000954; the preceding averaged
same-weight calculation is.003890. Here the unchanged minimum remains30 cells.

| Matching alternative | Weight mass | Effective majority/minority | Cell gate | Maximum nonzero support | Core |
|---|---:|---|---|---:|---|
| Full counts outside event, margins |22.710700|265.078 /22.895|Fails|Not evaluated after support failure|None|
| Proposal counts outside event, margins |22.710700|268.148 /22.895|Fails|Not evaluated after support failure|None|
| Original full depth, uniform base |62|349.206 /31|Passes|.131251708|None|

The first covariate changes300 cell weights, the second471; neither resolves
minority information support. Uniform base weights change all502 weights,
retain identical original depth edges and recover the required effective-cell
count. It fails the still unchanged genomic-core rule. Increased conditional
evidence alone is not a clone-recovery success. No new full-workflow merge,
false-split, power or biological uncertainty result is claimed.

All original comparisons and abstentions remain in the
[full-outside](../dev/outside_event_depth_20261109/full_outside_v1/report.json),
[proposal-outside](../dev/outside_event_depth_20261109/proposal_outside_v1/report.json)
and [uniform-base](../dev/outside_event_depth_20261109/uniform_full_v1/report.json)
reports. Each failed continuation gate stops its larger job chain.

## Compute and preservation

| Job | Role | Account/QOS; CPUs / RAM / requested time | Outcome |
|---|---|---|---|
|23151065|original numerics|nslab;1 /4GB /1h|complete0:0,5s|
|23151066|full-outside matching|innovation;1 /8GB /1h|complete0:0,67s; scientific gates fail|
|23151067|proposal-outside matching|innovation;1 /8GB /1h|complete0:0; scientific gates fail|
|23151251|unit-weight numerics|nslab;1 /4GB /1h|complete0:0,13s|
|23151252|unit-weight matching|innovation;1 /8GB /1h|complete0:0,73s; core gate fails|

All five execute on CPU compute node006 in the existing read-only validated
noneditable wheel/environment, with a neutral cwd and import/backend assertions.
Peak process RSS489672KiB; scheduler RSS is retained separately. Commands,
versions,resources,dependencies,snapshots and shared pre-edit baseline are in
[the task directory](../dev/outside_event_depth_20261109/). No arrays,GPU,
installation,production/native/original-test or frozen-output changes.
No full package test suite is rerun for research-only helpers.

The next justified mechanism must address more than thinning noise, total-depth
self-conditioning or duplicated classification margins. These attempts do not
authorize weaker support or a claim that more compute supplies missing genomic
identification. Production stays2/5; the earlier4/5 candidate's state-null
false-split failure and the simultaneous95% upper<=5% specimen/power target remain.
