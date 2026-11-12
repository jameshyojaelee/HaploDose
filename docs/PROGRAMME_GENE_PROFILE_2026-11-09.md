# Native RNA regions retain a useful bulk profile improvement

The programme-adjusted gene median on native RNA regions improves mean comparison
Pearson from .680245 to .738056. The conditional 18-comparison simultaneous interval
for the paired gain is [.012679,.104520]. Its point correlation exceeds strong
pseudobulk .730261, but superiority over that comparator remains unresolved.
Retain this as a bulk-specific opt-in candidate; do not change defaults or claim
general superiority. The programme-adjusted predictive smoother is rejected.
Original evaluations remain regression-only and frozen artifacts are unchanged.

## Method and scientific interpretation

The estimand is a regional natural-log relative RNA abundance contrast against
the specified normal panel. Existing pooled PLN programme/bin fitting supplies
expected bin counts. Allocate each fitted bin expectation to its genes in
proportion to their pooled within-bin reference counts, excluding each normal
from its own reference. Compute log((count+.5)/(expected+.5)), centre by the
median of within-bin gene-effect medians, then take the median of ALL gene effects
in each contiguous native RNA-HMM region. Reset at chromosomes, arms and genomic
gaps. The native path defines boundaries, not the exported effect magnitude.

The reference feature universe, native programme fitting, nuisance estimates and
region selection condition on the full RNA observations. This is descriptive
estimation, not held-out confirmation or an interval construction. Within-bin
reference composition must transport to the sample after nuisance correction.
Regional transcriptional state, reference mismatch, purity and expression response
remain alternatives to DNA dosage. Corresponding independently processed DNA at
the claimed resolution can support profile compatibility; it cannot validate
cell membership. Independent biological replication and an unopened evaluation
are absent. Absolute DNA copy number, purity and genome doubling remain unidentified.

## Locked experiment and comparator controls

[Design](../dev/programme_gene_profile_20261109/DESIGN.md) freezes the two
programme-adjusted candidates and a raw-gene median on the SAME native regions.
All12 prior comparators plus retained gene predictive and width1/5 medians remain,
18 methods total. Genes, sample identities, callable support and resolution-matched
DNA windows are fixed. RNA predictions/hashes are written before DNA scores open.
The source counts and original20-normal panel come from raw_refit_000 in
hardening_20261106_2119/profile_v6;24 calibration and24 comparison tumour patients
are distinct. There are23 finite comparison correlations and13 loss-bearing
comparison patients. Cells/windows are not biological replicates.

Native bin ratios, state paths, direct raw-gene effects and retained predictors
reproduce exactly for all68 samples. Independent Fraction expectation allocation,
scalar region medians, gene permutation, count conservation and gap resets pass.
These establish implementation identities, not biological calibration.

| Method | Mean Pearson | Loss recall | Loss precision (calling patients) | Loss FPR |
|---|---:|---:|---:|---:|
|Retained gene predictive|.680245|.311493|.284772 (16)|.042193|
|Programme gene / native regions|.738056|.495375|.348000 (15)|.037392|
|Raw gene / same native regions|.737790|.500183|.265789 (18)|.045332|
|Programme gene predictive|.687919|.525595|.335287 (18)|.053511|
|Strong naive / own partition|.730261|.399679|.358571 (14)|.041244|
|HaploDose precision-weighted path mean|.733069|.577106|.458667 (15)|.046896|

Precision is conditional on emitting calls and has unequal denominators. Do not
interpret it as an unconditional paired success rate. Against strong naive, the
native-region candidate has12 correlation wins/11 losses, recall4 wins/0 losses/9
ties, and paired precision change -.009524 on14 jointly calling patients. It
does not dominate the existing precision-weighted HaploDose path: recall is lower
by .081731. All18 methods, every patient and all metrics are retained in
bulk_v1/patient_scores.tsv, summary.tsv and report.json.

The matched raw-region control changes correlation by only .000267 relative to
the programme candidate. Most correlation benefit comes from the native region
boundaries, not proof of a major programme-correction benefit. Programme correction
reduces FPR in this matched comparison. This component attribution is descriptive.

## Reference, threshold and patient uncertainty

[Post-result uncertainty design](../dev/programme_gene_profile_20261109/UNCERTAINTY_DESIGN.md)
freezes the characterization before new refits/resampling, after the primary
point results. It is not a prospective confirmation analysis. Refit all eight
previously fixed normal-reference bootstrap panels and repeat panel0 exactly.
These are resamples of existing normals, not eight independent biological panels.

Use20,000 nested draws, seed2026110914: choose one of eight reference fits, resample24
calibration patients and24 comparison patients, and re-estimate each method's SAME
exact patient-balanced5% neutral loss threshold. Evaluate the same paired-patient
estimand. The family is2 candidates x3 baselines x3 metrics =18. Report approximate
basic bootstrap Bonferroni simultaneous intervals and Holm-adjusted centered
bootstrap p-values, including missing/failed draws and negative results.
This family covers the present post-result characterization, not the cumulative
adaptive method-development history. It supplies no prospective superiority claim.
This bulk neutral-window calibration target is distinct from the5% whole-specimen
clone false-split criterion; passing either would not establish the other.

Empirical-CDF acceleration proposes quantiles only: exact integer/Fraction mass at
each returned threshold and predecessor checks the original strict rule. All1845
threshold checks and4320 lookup-score comparisons pass. Independent fsum reductions
and primary prediction/score identity pass. Monte Carlo endpoint brackets use
binomial order statistics across36 endpoints. No threshold is retuned.

For native regions versus retained gene predictive, paired Pearson gain .057811
has simultaneous interval [.012679,.104520], Holm p=.005400. Against strong naive,
gain .007795 has interval [-.017219,.029730]; against HaploDose precision-weighted
path mean, gain .004988 has interval [-.021866,.021223]. Both remain unresolved.
Loss recall versus retained has interval [-.192891,.531759], FPR improvement
[-.036988,.051696]. Loss advantages are unresolved too. The programme predictive
combination is worse than strong naive in correlation, interval
[-.076519,-.005622], Holm p=.007650; reject that addition.

All18 contrasts and reference-specific results are in uncertainty_v1/report.json
and its accompanying tables/draws. Eight reference draws, post-result selection
and the finite patient panel limit interpretation. Approximate bootstrap intervals
do not establish repeated-workflow coverage, prospective FPR or biological CNA
validity. Additional draws cannot substitute for missing independent specimens.

## Execution and decision

All12 jobs complete exit0:0 on CPU compute nodes. Numerical23145825 uses nslab
1CPU/4GB/1h (27s); primary23145826 innovation4CPU/16GB/3h (73s). Refits23146991–999
use alternating innovation/nslab4CPU/16GB/3h,59–83s each, recorded batch RSS about
365–375MB. Nested uncertainty23147000 uses nslab1CPU/8GB/3h (24s); process resource
records are retained even though its scheduler RSS reports0. No GPU is used.

Source_v1/v2, jobs.json, installed wheel/native identity, neutral-cwd import checks,
environment inventories, protected paths and per-job resource records are under
[the evidence directory](../dev/programme_gene_profile_20261109/). All failures
and old results elsewhere remain intact. No package change follows from these
research jobs themselves. Native integration requires separate installed-artifact
tests and exact replay. Clone failures and the specimen false-split/power criterion
are unaffected and remain unmet.
