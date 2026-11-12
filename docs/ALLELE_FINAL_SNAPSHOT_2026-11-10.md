# Keep final assignments and stored allele parameters on the same snapshot

COMPLETE:the candidate wheel passes483Python and4Rust tests,with1GPU-only
skip and8existing warnings. All three new consistency fixtures now pass;
the old wheel had failed all three on the intended assertions. All3jobs
complete0:0. Candidate validation23160542 takes3m28s,including164.69s pytest;
private environment23160541 takes4m30s. Native binary SHA is unchanged.
This establishes the bounded snapshot-consistency repair,not convergence or
biological calibration. No existing test or scientific criterion is relaxed.
[Complete Python output](../dev/allele_final_snapshot_20261110/validation_23160542/pytest.txt),
[Rust output](../dev/allele_final_snapshot_20261110/validation_23160542/cargo_test.txt),
[installed artifact identity](../dev/allele_final_snapshot_20261110/validation_23160542/installed_artifact.json),
[accounting](../dev/allele_final_snapshot_20261110/accounting_completed.tsv),
[resource evidence](../dev/allele_final_snapshot_20261110/resources_completed.json),
[preservation](../dev/allele_final_snapshot_20261110/preservation_completed.json).

Earlier preparation and baseline checkpoint below retain their original context.

Implemented the bounded correctness correction from finding14 of the
[scientific review](scientific_code_review_2026-11-09.md). The existing installed
wheel fails all THREE new consistency fixtures on a compute node:job23160540
exits0:0 because the expected test failures are confirmed,not because those
tests pass. The inner pytest exits1 with exactly three assertion failures and
zero execution errors. New-wheel validation remains pending.

Previously `SpecimenModel.em` computed final LL and responsibilities R,then
refitted allele fractions B and allele posteriors from R before returning the
earlier scores. The stored allele state could thus differ from the one used
to assign cells. This is separate from finite-iteration convergence and RNA
assignment-versus-export profile conventions.

The smallest repair removes that final post-scoring allele-only update. Keep
the last refinement-time fitted allele parameters,their fit pools and their
posterior unchanged through final scoring. Returned LL/R use that state; final
RNA pools are computed from final R. The allele fit pools need not equal pools
from the final responsibilities:they identify the finite fit state actually
scored. No extra refinement or fixed-point guarantee is implied. The source
docstring states this explicitly.

An AST comparison of complete pre/post `clones.py` proves that,apart from that
docstring,the sole executable change removes the final `if self.alleles is not
None` refit. All other source statements remain identical. Consequently the
allele-free path has the same operations; this is source equivalence,not a new
historical-output rehash or a fresh biological-performance claim. No existing
test,scientific cutoff,default,recovery criterion or frozen result is changed.

Independent numerical check:four cells with AD9/8/2/1 and DP10,one balanced
component and one fitted-binomial component. The fit is weighted AD/DP and
assignment uses a directly calculable binomial likelihood ratio. Compare the
returned LL to a fresh stored-B score EXACTLY,returned R to its softmax EXACTLY,
and RNA pools to returned R EXACTLY. Independent SciPy full binomial PMFs agree
within64*machine-epsilon*totalDP,an arithmetic allowance for this fixture.
Run one refinement,eight budgeted refinements and MAP-heuristic stopping.
These are computational checks; biological n=0 and calibration is not tested.

| Job | Purpose | Account/QOS | CPU | RAM | Wall | Status |
|---|---|---|---:|---:|---|---|
|23160540|old-wheel counterexamples|nslab|1|4GB|30min|COMPLETE;3 intended failures confirmed|
|23160541|isolated conda prefix|innovation|1|8GB|1h|RUNNING at preparation checkpoint|
|23160542|new wheel,all Python/Rust tests|nslab|2|32GB|1h|PENDING after environment|

No shared installation is modified. Both source snapshots include the dirty
working tree and contain922files; only `clones.py` differs between them. The
new fixtures exist identically in both snapshots. Existing power/false-split
and whole-cell experiments retain their original immutable source/wheels and
scientific decisions. This repair supplies no new FPR,power or DNA-validation
claim.

[Changed controller](../python/haplodose/clones.py),
[new fixtures](../tests/test_allele_final_snapshot.py),
[locked design](../dev/allele_final_snapshot_20261110/DESIGN.md),
[old-wheel failure record](../dev/allele_final_snapshot_20261110/baseline_23160540/report.json),
[full failure details](../dev/allele_final_snapshot_20261110/baseline_23160540/pytest.xml),
[source difference](../dev/allele_final_snapshot_20261110/source_difference.json),
[all commands and handles](../dev/allele_final_snapshot_20261110/jobs.json).
