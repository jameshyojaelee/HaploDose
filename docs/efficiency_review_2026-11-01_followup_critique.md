# Critique of `efficiency_review_2026-11-01_followup.md`

Written 2026-11-01 by HaploDose2. Each code claim was checked against the current source, not the review's line
numbers: `clones.py` and `profile.py` changed after the review fingerprints (clones.py sha256 now 169a460b...;
alleles.py, pipeline.py, bulk.py, kernels/estep.py unchanged). No code was run except reading result JSONs.

## Verdict

The review is accurate on code behaviour. All eight main findings hold in the current source. Its arithmetic
re-derives (purity effective beta 0.152/0.235/0.263; beta-binomial bias 0.04412; equicorrelated SE ratio
sqrt(6.7) = 2.59; 1M x 200 float32 = 0.75 GiB). Its weaknesses: it does not rank fixes by cost against the one
benchmark currently gating work (Kinker implants), it overstates the allele fix, and it misses that the
no-window v3 run already exists.

## Finding by finding

| # | Claim | Checked at | Holds? | My priority and why |
|---|---|---|---|---|
| 1 | Dedup keeps disjoint rare groups out (f + g <= 0.1) | clones.py:844 | Yes. 920/1000 agreement for two disjoint 40-cell groups | **High, cheap.** Same function the split-proposal work touches. Compare minority sets (Jaccard on the smaller side, allow label flip). |
| 2 | Only `members[:2]` validated, all children accepted | clones.py:864 | Yes | **High, cheap.** Accept binary splits only; extra children get proposed again next round. STATUS records false splits starving a true subclone (CAOV3), consistent with weak validation. |
| 3 | Squared allele stat cannot see opposite-haplotype events | alleles.py:252-265, 345-394 | Yes | **Medium.** The review's fix (a shared phase model) is heavier than needed. Counts are already oriented to one input phasing for all cells (alleles.py:4), so two groups share orientation per SNP. A signed per-SNP two-group beta-binomial test needs no new phase model; switch errors affect both groups equally. |
| 4 | `_update_offset` omits depth and programme terms | clones.py:695-702 | Yes | **Medium, cheap.** `_update_phi` (clones.py:704) already builds the full per-cell mean in chunks; reuse it. The 100x depth probe is extreme; real effect unmeasured. |
| 5 | Reference variance `v` counted twice in tau | profile.py:107, 110, 114 | Yes. Analytically sqrt(0.09 + 0.09) = 0.424 vs 0.3 | **High, one line, but gate it.** `estimate_tau` already takes `extra_var` (profile.py:182 uses it). Inflated tau makes calls conservative; may explain part of region_found 0.29-0.31 vs infercnvpy 0.67. Must recheck neutral false calls (U-CAN 34/423) before adopting. |
| 6 | Clone 0 gets no expression proposals | clones.py:729 | Yes | **Medium.** Deliberate null control; relevant to U-CAN U0396/U0991 (no tumour clone). Allow with a stricter threshold, test on all-normal null. |
| 7 | `differences.tsv` and tree skip `supported()` vetoes | pipeline.py ~340-357 | Yes | **Medium, cheap.** Changes benchmark numbers: Andor SNU-601 Haplodose precision (0.77) and called arms come from this table. Rescore after the fix. |
| 8 | Type-mismatch gene mask is a union over datasets | pipeline.py:211-233, 74 | Yes | **Low-medium.** log 5 at beta 0.6 needs ~29 copies to trip from dosage alone, so amplicon loss is rare; the cross-dataset union is the real bug. Make the mask per dataset. |

## What the review missed or got wrong

- **Sliding windows: the no-window v3 run exists and settles it.** `results/kinker_implant_v3_nowin_20261101_2021`:
  Haplodose ARI 0.153, region_found 0.307. With windows: 0.147, 0.289. infercnvpy 0.165, 0.675 (n = 166 lines,
  paired Wilcoxon, Holm over 1 pair: p 0.18 no-window, 0.087 windows). Windows cost 2.7x runtime at 100k cells
  for no accuracy gain. Recommend `event_windows` off by default until a proposal change beats no-window.
- **The gap is detection, not clustering.** region_found 0.31 vs 0.67 says Haplodose rarely reports the implanted
  region at all. Findings 1, 2 and 5 act on that; the efficiency section (prefix-sum screens, Rust) does not.
- **Allele fix scope** (finding 3): see table. No phase model needed.
- **TCGA FPR:** correct that 5% is per side (subarm.py:269, documented at lines 23-29) and set in-sample on the
  scored type's DNA-neutral windows. Both methods get the same rule, so the comparison is fair; only the absolute
  threshold is not deployment-calibrated.
- **0.3 s cProfile run** is too small to direct optimisation. The measured 100k timing is the useful number.
- **Competitor list** (CopyKAT, SCEVAN, XClone, CalicoST, CaSpER): fine as a goal, but none is built. Do not widen
  the comparison until the Kinker detection gap closes.

## Proposed order

1. In clones.py (owner HaploDose1, during current split-proposal work): finding 1 dedup, finding 2 binary-only
   acceptance, windows off by default.
2. profile.py finding 5 (tau with `extra_var`), gated on U-CAN neutral false calls and Kinker null false splits.
3. clones.py finding 4 (offset reuses `_update_phi` mean).
4. pipeline.py finding 7 (report through `supported()`), then rescore Andor and Kinker.
5. alleles.py finding 3 (signed per-SNP test) with a mirrored-event unit test; pipeline.py finding 8 per-dataset mask.
6. Defer: purity-aware model, adaptive bins, multiscale screen, Rust kernels, new competitors.

Acceptance for 1-4: Kinker v3 implant ARI and region_found on the same 166 lines vs the no-window v3 run,
all-normal null false splits not higher, U-CAN neutral arms called not higher than 34/423.

## Outcomes (2026-11-01, HaploDose2; `results/tau_fix_20261101_2117/`)

Code trees were frozen copies under `dev/hd2_20261101/` (base = repo snapshot about 21:15, with shared_profiles).

- **Finding 5 (tau counted twice): real in code, no measurable effect. Not adopted.** Kinker v3: every specimen's tau is
  identical, so v is about 0 there. U-CAN spatial, 20 patients: tau drops, rho vs base 9 higher / 9 lower (p 0.44),
  neutral arms called 32 -> 36/423. U-CAN bulk, 1,034 tumours: dose rho 0.649 vs 0.653, p 0.62. Paired Wilcoxon, unit =
  patient or line. The bulk arm-z inflation also had to change (tau / normal tau -> median sqrt(tau^2 + v)), or a
  test failed. The fix stays only in the dev copy.
- **Finding 7 (reports skip vetoes): applied** to `python/haplodose/pipeline.py` (summarise reports only
  `model.supported()` segments). Repo tests: 58 pass. The 1 failure, test_spatial Potts, is in clones.py work in
  progress.
- **Finding 8 (mismatch mask union): not adopted.** Kinker v3, 166 lines, all on v4 nodes, ARI vs base (0.1686,
  region_found 0.404): per-dataset mask 0.1586 / 0.368 (54 higher / 71 lower, Holm p 0.21); no mask 0.1589 / 0.361
  (48 / 69, Holm p 0.21). Holm is over the two variants. The union mask removes cell-line-specific noise and helps here.
- **CPU generation changes results.** The same code on v4 and on v5 nodes gave Kinker ARI 0.1686 vs 0.1644. Run A/B
  arms on one node feature, with `--constraint` and no `--prefer` in the script. A command-line `--constraint` with
  `--prefer=v6` in the script still landed on v6.
- **Finding 3 (squared allele statistic): fix written, parked, not applied.** A per-SNP two-group heterogeneity test,
  `sum (x_a - x_b)^2 / v` under a pooled beta-binomial null, is in `dev/hd2_20261101/allele/` (alleles.py
  `_heterogeneity`, used by `differ` and `veto`; two new tests). It separates mirrored 0.8/0.2 imbalance (z 134;
  the old rule calls it no difference). Null false rate of the differ rule (z >= 4 and gap >= 0.1; 2,000 simulations,
  60 SNPs, kappa 50): depth 2 vs 200 per SNP, 0.95% vs old 0.20%; 5 vs 200, 0.80% vs 0.15%; 10 vs 100, 0.25% vs 0%.
  At equal depth 30 it is calibrated (z mean -0.01, sd 0.91). Reasons to park: heavier null tail at low pooled SNP
  depth, which is common in scRNA; no allele benchmark with DNA truth (YCC is scored against Numbat only); and
  allele split proposals cannot raise mirrored subclones, since the pooled clone fraction is about 0.5. Revisit
  together with proposals and a variance term for q at small depth (binomial kurtosis).
