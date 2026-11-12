# Proposal-region retention:bounded explanation falsified

The downstream state-versus-population comparison rejects both a state error and
its paired true contrast. This investigation tests one specific upstream cause:
native event proposals deduplicate overlapping cell partitions before the new
score evaluates their genomic measurements. Such deduplication could discard
useful intervals,but that possibility must be checked on the actual execution.

The estimand remains discrete regional RNA predictive support against the same
continuous-state competitor. The measurement chain is the first native event
proposal on the same exposed seed202611099001 under each condition. Corresponding
DNA would be required for genomic membership;these simulations provide biological
n=0 and only one independent base seed. No error-rate estimate follows.

Observation-only profiler23210672 captures the native positive-gain found list
before Jaccard deduplication,all returned sources and source arrays. Counts,
expected bin/gene intensities,gene shapes and global cell identities match the
already saved first comparison EXACTLY. Every native inference output manifest
hash verifies. The procedure itself is unchanged.

| Condition | Eligible chr12 intervals before deduplication | Retained after deduplication | Source |
|---|---:|---:|---|
| State | 1 | 1 | event_chrom12_bins211-219 |
| Two CNA populations | 1 | 1 | event_chrom12_bins211-219 |

There are no discarded eligible windows in this pair. Region retention therefore
cannot repair its observed power loss. Do not add windows,lower native gain/MAP
requirements or move to a favourable later parent. The prepared downstream
rescoring implementation was NOT RUN and is explicitly retired for this pair:
there are no additional eligible measurements to score. The existing negative
score reports already cover the complete captured family. This observation does
not prove that deduplication is harmless for all datasets.

## Execution and artifacts

First attempt23210646 fails1:0 in29s at the existing hook's specimen-path
assertion;the new launcher omitted seed_202611099001 from its output path.
Its source/output remain. Corrected snapshotv2 preserves the hook and allocation,
adds that directory component,and runs to0:0 in1m22s on ne1dc4-015. Requested
nslab/nslab CPU,Intel-v4,1CPU8GB1h,existing read-only Conda environment,one thread.
No scientific criterion,public code,default or frozen artifact changes.

[Frozen contract](../dev/event_region_retention_20261111/DESIGN.md),
[exact source correspondence and found lists](../dev/event_region_retention_20261111/captured_v2/REPORT.json),
[both execution attempts](../dev/event_region_retention_20261111/EXECUTION.md),
[submission commands](../dev/event_region_retention_20261111/submissions.json),
[resources](../dev/event_region_retention_20261111/job_23210672/resources.txt),
[prepared rescoring contract](../dev/event_region_score_20261111/DESIGN.md).
