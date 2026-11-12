# Where clone subdivisions enter the workflow

The objective is to recover supported CNA distinctions while meeting the unchanged
specimen false-split criterion. A correction that merely removes all subdivisions
does not satisfy that objective. This development investigation locates the origin
of extra groups before testing an initialization-only intervention.

## Exact replay trace

The preselected cases are replicate 0 from matched-NB, correlated-state and
two-clone simulations, each under merge_control and exact_matched. The selected
six are mechanism probes, not estimates of error or power. Frozen generator,
counts, parameters, genomic resets, model wheel and scientific criteria are reused.

Job22963745 reproduces every saved row field (apart from output-file provenance)
and every assignment array exactly. Read-only instrumentation records initial
groups, EM transitions, split transitions and existing merge logs. It leaves the
original merge/agglomeration method objects intact because experimental hooks
identify those caller code objects. No oracle enters inference.

| Condition | Initial Leiden groups | After initial agglomeration | Final eligible tumour groups | Later accepted splits |
|---|---:|---:|---:|---:|
| Matched NB | 7 | 1 | 1 | 0 |
| Correlated state, one CNA population | 5 | 3 | 3 | 0 |
| Two CNA populations with state variation | 6 | 4 | 3 | 0 |

Both paths produce these counts. The fourth alternative component has only one
MAP-assigned cell and does not meet the unchanged eligible-group support rule;
its soft membership mass is about 40.7. The reference component is excluded from
the table. No matched merge occurs in these replays.

[Inference] In these cases, the unsupported subdivisions survive from
initialization rather than being introduced by later accepted splits. Modifying
split confirmation therefore cannot remove them. This explains a specific route
through the observed workflow, not all 200 specimens or biological tumours.

The first attempt, job22963669, failed a provenance comparison because JSON
converts integer dictionary keys to strings. Job22963745 normalizes both record
representations to JSON and retains exact array equality; no numerical or
scientific tolerance changes. The failed output and source remain preserved.

## Prespecified causal intervention

The smallest intervention pools initial tumour-group counts, exposures and depth
moments before the first EM, refits their single pooled profile, and retains
reference/tumour mixture mass. It reuses the existing initial-pooling operation
without enabling the separate selected-parent filter. Subsequent initialization
of later proposals, confirmation, merge criteria and stopping remain unchanged.
This is a research hook, not a production change or calibration procedure.

Run all three mechanism cases and all five original recovery cases. The matched
NB case must remain bitwise identical because it already starts with one tumour
group. Improved behavior requires retaining genuine contrasts while reducing
unsupported distinctions. Loss of true distinctions falsifies this proposed repair
even if the state-null output becomes simpler. A promising small-panel result
would still require complete-workflow error and power assessment.

Jobs22963866 and22963868 completed. The matched NB control remains exactly
unchanged. The state-null example changes from three eligible groups to one;
the two-clone example retains three, so correct clone number is still not recovered.
The original five-case recovery panel remains **4/5**, with seed20261109 still
returning one S1 tumour group, accuracy0.602 and no chr12 distinction. Every
other seed passes all original assertions. Seed20261108 accuracy is0.7569,
normal false calls8/90 and tumour sensitivity506/510. No criterion is weakened.

Audit22963958 verifies142 recorded input hashes, six exact replays, independent
pooled count/depth/moment conservation, every original recovery criterion,
386 preserved baseline files and65 unchanged package/test/native files. The
depth-moment reference uses scalar sums and a float64 accumulation bound.
Mechanism/recovery jobs took2m11s and5m56s; process peak RSS was882MiB and1.38GiB.
The audit took20s. No unchanged full package suite was rerun.

This warrants a complete calibration/power experiment, not default promotion.
The new experiment under `dev/pooled_initial_null_20261108` retains199 calibration
draws per null,200 evaluation draws per condition, all old seed keys and the20-bin
reset. Five jobs22964075/22964077/22964078/22964079/22964080 run998 fits. Summary
22964088 and audit22964195 require successful upstream completion. Nine new
endpoints join the45 historical ones in54 simultaneous one-sided bounds; original
5% operating and scientific criteria remain. No outcomes are claimed while running.

Sources, designs, outputs, manifests, wheel identities and failed attempts are
under `dev/initialization_trace_20261108`. The experimental pooling trace and outer
manifest identify this intervention separately from production selected-null flags.
It does not enable the selected-parent filter or confer independent confirmation.
