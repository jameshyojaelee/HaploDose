# Joint membership recovery and a terminal search defect

## Terminal check confirms stopping-phase recovery

Job22988135 completes all five cases with unchanged merge criteria and finite
merge/refit finalization. Three cases still pass. Seeds20261109 and20261110 end
with one tumour group, accuracies0.601202 and0.598, and no qualifying chr12
distinction. Both fail the original two-group, accuracy and difference criteria.
Exactly two of15 specimen fits change, both S1 from three total groups to two
(normal plus one tumour group). Each terminal merge decreases group count;
a repeated final check changes nothing. The failed recovery results are retained.

This establishes that the earlier five-case pass depends on the unvalidated
last-cycle split. Production source now finalizes an exhausted search by EM,
merge validation and refitting until no merge remains, with a strict decreasing-
group guard and explicit stopping provenance. Four hand-calculable counterexample
tests cover stable output, one/repeated merges and a nondecreasing invalid merge.
New isolated wheel validation22991908 passes325 Python tests and four Rust tests,
with one GPU-only skip. The wheel SHA256 is
`6712fa766754432dbcc48305c8dd1be5b76aceec3cc7277a1b64d2b03b1e9000`.
Three matched recovery replays and independent audit22992925 now complete:

| Procedure | Before terminal correction | After terminal correction |
|---|---:|---:|
| Original thinning path | 2/5 | 2/5 |
| Pooled/exact conditional path | 4/5 | 2/5 |
| Joint unordered path | 5/5 | 3/5 |

The audit independently recounts90 criteria across15 cases, verifies435 frozen
baseline files and67 current core/test/native files, and reproduces the separate
terminal probe exactly in all five joint posterior arrays and exported
cell/profile/difference/segment tables. The original thinning failures108/109/110
remain. The pooled/exact path loses108/110 after terminal validation; the joint
path retains108 but loses109/110. These regressions are required by the existing
merge criterion, not grounds to weaken it. The
older wheels, numerical screen and terminal probe remain frozen. No larger joint
calibration experiment is justified by the unstable five-case result.

The goal is to preserve identifiable RNA-profile and membership information
through proposal, confirmation, refitting and merging, while controlling false
distinctions. The joint-count prototype passes all six numeric recovery criteria
in five locked development cases. A subsequent trace finds an unvalidated final
split after the search budget is exhausted. The numeric pass is therefore not
yet evidence of stable recovery, and no new large calibration sweep is launched.

## Estimand and comparison

The candidate uses an unordered pair of regional RNA-rate trajectories and cell
memberships. It retains the NB working likelihood and integrates a shared cell
membership across the entire short interval. One global label exchange preserves
relative signs across loci. Counts enter the likelihood once, with a common
population prior fitted during proposal search; fitted per-cell probabilities
are not fed back as priors on the same counts.

This changes the statistical procedure. Full-count likelihoods and memberships
are not independent thinning confirmation. QC records that distinction before
output hashes are finalized. All previous ordered/thinning failures remain.
Depth-matched subsets and nuisance estimates remain selected; the posterior is
not a calibrated genomic confidence statement. Continuous expression state and
reference mismatch remain competing explanations requiring complete-workflow
null experiments and, ultimately, corresponding DNA evidence.

All five original inputs, gene vetoes, minimum cell counts,0.95 posterior at
every retained bin and five-bin support requirement remain fixed. Two arms
differ only in ordinary versus previously implemented margin-weighted merging.
Long events beyond the declared exact-enumeration cap and non-event proposals
retain the prior procedure. This is a bounded research hybrid, not a new default.

| Seed | Ordinary merge accuracy | Margin merge accuracy | Six numeric criteria |
|---|---:|---:|---|
| 20261031 | 0.846465 | 0.846465 | pass in both |
| 20261107 | 0.738095 | 0.738095 | pass in both |
| 20261108 | 0.760870 | 0.764822 | pass in both |
| 20261109 | 0.771717 | 0.771717 | pass in both |
| 20261110 | 0.812000 | 0.812000 | pass in both |

These accuracies use the original criterion conditional on cells called tumour,
alongside its separate exactly-two-groups, normal-FPR, tumour-TPR, one-group-S2
and chr12-only distinction requirements. They are not estimates of biological
population accuracy. No inferential comparison is claimed from five exposed
simulation cases.

## Where membership information disappears

For the first chr12 proposal in seed20261109, scoring-only truth gives accuracy
0.694 for the original proposal labels and0.788 for the joint labels on the same
500 cells. At the subsequent merge, accuracy over all510 true tumour cells,
penalizing unassigned cells, is0.7490 using full-count memberships but0.6667
using proposal-only memberships. These denominators differ from the final-table
conditional accuracy above and must not be conflated.

Both arms then log **split, merge back, split again**. The final split occurs on
the last of three allowed search cycles. `fit_specimen` performs a final EM
refit after budget exhaustion but does not run `merge_unsupported` on that final
model. [Inference] The reported recovery may depend on where the oscillation is
stopped. The raw traces establish the action sequence; a terminal replay is
needed to establish whether the returned group distinction actually survives.

The next experiment applies the unchanged merge rule, refits, and repeats until
no unsupported group remains. Each merge strictly reduces the group count, so
termination is bounded without choosing a favorable split/merge phase. Any loss
of recovery will be retained. A production correction needs its own minimal
counterexample and rebuilt isolated wheel; no larger joint calibration run is
justified until this correctness issue is resolved.

## Independent checks and retained failures

Jobs22983512/22983513 complete in6m19s/5m52s, with peak RSS1524288/1476372KiB.
Audit22986640 independently recounts all60 numeric criteria, aligns99 membership
records using count-file cell IDs, verifies output/missingness correspondence,
and preserves428 baseline files and66 core/test/native files. Full-PMF numerical
fixtures run in both arms. The installed package and original tests are unchanged.

The first attempts22981710/22981711 abort before any case score on a nearly
neutral event with numerically unresolved arbitrary memberships. Version2 rejects
unsupported cores first and explicitly abstains from numerically unresolved
candidates; neither tolerance nor evidence threshold is relaxed. Failed outputs
remain frozen. Audit22983514 fails only while serializing a NumPy boolean;
22986105 flags identical missing object values as unequal. The final audit
compares missingness masks and all observed values exactly. These failed audits
remain preserved; no inference rerun follows either reporting/checking correction.

Evidence: `dev/joint_unordered_20261108/{hard_v2,margin_v2,audit_v4}`. The terminal
probe is job22988135 under `dev/terminal_merge_20261108`. Biological validity,
specimen error control, stable recovery and the original thinning failures remain
unresolved. A future candidate/control comparison should use fresh simulation
draws after freezing, retaining the exposed draws as regression evidence.
