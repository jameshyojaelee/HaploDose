# Conditioning profiles on the shared offset does not rescue recovery

The research intervention matches a known-offset conditional mean but fails the
full recovery task: 0/5 original cases with the production offset update and
1/5 with the conditional NB offset optimizer. The frozen comparator arms recover
2/5 and 1/5 respectively. Reject this intervention as a recovery correction or
package default. Preserve the forward-model discrepancy and negative outcomes.

## Scientific question and locked intervention

Task `dev/profile_offset_20261109` saves `DESIGN.md` before new fitting. Historical
profile expected counts use reference composition and depth moments, whereas
cell likelihoods additionally use the current shared regional offset. This can
make profiles absorb an already fitted RNA effect rather than estimate a
conditional residual. The competing explanation is that the historical procedure
is a useful alternating heuristic whose profile/offset steps do not share a
single objective. Neither interpretation alone establishes a production bug.

The estimand in the numerical gate is expected pooled RNA count given baseline,
offset and depth moments. The full-workflow endpoints remain all six original
recovery assertions per seed. Biological n=0. There are no new biological
replicates, prospective cohorts or false-split guarantees.

Inside `profiles` only, the intervention tilts reference probabilities by the
current offset and renormalizes. Depth slopes are centered under that tilted
baseline: subtracting a common slope is an exact softmax gauge for cell means
and is required by the depth approximation's centering assumption. Original
baseline and slopes are restored before every other operation and on exceptions.
Constant offsets preserve the original baseline/slopes exactly.

Reference heterogeneity variance, programme fitting, segmentation, search,
confirmation and all scientific criteria stay fixed. Pairwise shared-profile
differences still use the historical baseline outside `profiles`. This is a
bounded component intervention, not a claimed coherent replacement of every
working model. The offset remains fitted from RNA; it is not an independent
normal or DNA anchor. Offset/common-dose nonidentifiability is unchanged.

## Independent forward gate and retained fixture failure

The first numerical job23107595 correctly rejects the fixture's unnormalized
reference columns before running the intervention. The fixture is corrected
to columns summing to one; no contract or tolerance is changed. Failed source,
output, exception and resource record remain in `source_v1/`, `numerics_v1/`
and `analysis_23107595/`.

Corrected gate23107787 passes36 cases spanning flat/heterogeneous baselines,
zero/constant/regional offsets, zero/nonzero slopes, and zero/nonzero depth
variance. Actual native-profile expected counts agree with an independent
scalar Decimal80 implementation of the depth-moment approximation. This verifies
the formula, not the accuracy of the approximation to arbitrary depth mixtures.
At zero depth variance the forward normalized mean has no depth-mixture
approximation. Zero-offset profile outputs replay bit-identically. Count
conservation, offset-gauge invariance and restoration after a native-profile
exception pass, using operation/magnitude-based float64 budgets.

The known-offset/no-dose example has four identical integer count rows:100 in
each of20 bins and200 in each of20 bins, exposure6000 per row, uniform reference,
zero depth slopes and offset0/log2 in the two regions. Corrected pooled expected
counts match400/800, while historical expectations are600 in every bin. This
establishes a conditional mean discrepancy when offset is known; it does not
establish how an offset estimated from the same RNA should be interpreted.

## Original recovery cases and conditional likelihood

Both new arms use the same pinned production wheel and all five original frozen
inputs, without truth in inference metadata. Frozen production and NB-only arms
from `offset_objective_20261109` remain fixed development comparators.

| Original seed | Conditioned profiles + production offset | Conditioned profiles + NB offset |
|---|---:|---:|
| 20261031 | Fail | Fail |
| 20261107 | Fail | Fail |
| 20261108 | Fail | Pass |
| 20261109 | Fail | Fail |
| 20261110 | Fail | Fail |

With the production offset update every S1 tumour partition has one group.
With NB offset optimization only20261108 recovers two groups. The surviving
case does not demonstrate a new interaction benefit: the frozen NB-only arm
already passes that case. All60 original criterion values across the new arms
are independently recounted, including strict normal-FPR failures; no case or
criterion is dropped. Complete criteria and fitted outputs remain saved.

| New arm | Applied profile Q decreases | Offset Q decreases | Recovery |
|---|---:|---:|---:|
| Conditioned profiles + production offset |21/53|19/53|0/5|
| Conditioned profiles + NB offset |46/62|0/62|1/5|

These are dependent update trajectories, not error rates. The NB offset gains
remain positive(minimum0.264113), with no reported optimizer failure, while
applied profile updates can decrease the cell NB objective substantially.
Pooled PLN segmentation, shared restrictions and nuisance estimation do not
constitute an NB M-step merely because their baseline is offset-conditioned.
The observation rejects this conditional-mean change as a sufficient recovery
solution. It does not isolate which remaining operation causes each failure.

## Audit, execution and decision

Independent audit23108032 passes60 original criteria,60 bounded expected-count
contexts using Decimal80 scalar calculations, and90 worst conditional-Q contexts
using independent full SciPy NB logpmfs and bounded Decimal80 ratios. Contexts
are chosen without truth: first constant/nonconstant offset within each model,
and minimum Q within model/kind. Retained arrays exactly match their pre-filter
subsets. All102 production-offset and121 NB-offset profile calls restore the
original baseline/slopes and preserve the offset. Actual inference QC records
CPU backend and device `cpu:0` for every specimen run.

Authoritative artifacts under `dev/profile_offset_20261109/`:

- `source_v2/`, `numerics_v2/report.json`: executed source and numerical gate.
- `profile_legacy_v1/`, `profile_NB_v1/`: all fits, criteria, raw/applied/offset
  traces, profile-restoration traces and bounded expected-count contexts.
- `audit_v1/report.json`: independent means, likelihoods and original recovery.
- `source_integrity.json`: all76 existing package/native/test/README/build files
  remain identical to the initial dirty source.
- `jobs.json`, `accounting_final.tsv`, `analysis_JOB/`: commands, source/wheel/
  import identities, resource headers, runtimes, process peaks, exits and failures.

Five CPU job attempts are terminal: four successful jobs and the preserved
invalid-fixture failure. Gate:1CPU/4GB/1h innovation. Each five-case arm:4CPU/8GB/
2h on innovation or nslab. Audit:1CPU/8GB/1h nslab. No arrays. Actual fit runtimes
68.04s/97.83s and process peaks1,572,776/1,540,840 KiB; audit38.35s/230,864 KiB.
AST checks pass. The full production suite was not rerun because production
code and tests are unchanged. Shared installations and frozen outputs are intact.

Research code is `profile_offset.py`, `profile_offset_numerics.py`,
`profile_offset_audit.py`, and a default-off option in `offset_objective_panel.py`.
No production option or dependency is added. Computational status: the stated
intervention and independent arithmetic checks pass. Identification: the fitted
offset is not an independent anchor. Biological validation: absent(n=0).
Reject promotion. Additional error/power simulations are not used to excuse a
failed recovery gate; no new null/power claim is made. Original production
recovery remains2/5, and false-split/power and uncertainty failures remain
unresolved under the unchanged criteria.
