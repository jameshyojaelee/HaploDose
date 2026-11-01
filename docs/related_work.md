# Related tools and what Haplodose takes from them

Checked 2026-10-31 against papers, docs or source code (PubMed/PMC full text where
marked F; docs or source D; abstract only A). Details not confirmed are marked
UNVERIFIED. Supplementary files were not read.

| Tool | Expression model | Alleles and phasing | Clones | Spatial | Implementation |
|---|---|---|---|---|---|
| Numbat (F, D) | Poisson-lognormal vs reference; reference mixture by constrained least squares | Eagle2 + panel; pseudobulk joint HMM, 15 states (mirrored haplotype pairs); switch probability (1 - exp(-2 nu d))/2, d in cM; per-cell binomial on alleles phased by pseudobulk posteriors | Iterative: subtree pseudobulk HMMs, per-cell posteriors, ML phylogeny (modified ScisTree) | Visium supported (D) | R, CPU (mclapply) |
| inferCNV (D) | Moving-average residuals; i6 HMM calibrated by spike-in; Bayesian mixture filter | none | Leiden subclusters, HMM per subcluster | none | R |
| infercnvpy (D) | 100-gene running window, log fold-change clip; no HMM documented | none | UNVERIFIED | none | Python |
| CopyKAT (F) | Freeman-Tukey + DLM smoothing; Poisson-gamma MCMC per window | none | Ward clustering | none | R |
| CaSpER (F) | Multiscale smoothing, 5-state HMM | per SNP, unphased, pooled | UNVERIFIED | none | R + C++ |
| HoneyBADGER (F) | Gaussian HMM; per-cell Bayesian model (JAGS) | 2-state binomial HMM on pooled lesser allele; no switch model | recursive cluster-pool-HMM | none | R |
| SCEVAN (F) | Total-variation smoothing; Mumford-Shah segmentation | none | Louvain, re-segment, phylogeny | none | R |
| XClone (F, D) | NB GLM, dispersions from reference cells | phased VCF, then EM phasing in 100-gene bins, then arm-level phasing | kNN smoothing + HMM | none | Python |
| CalicoST (F, D) | NB expression + beta-binomial BAF on variable bins | Eagle2, then pseudobulk HMM over (phase x BAF) corrects switches | joint HMM + HMRF; phylogeny from LOH events | HMRF, per-spot purity | Python (Gurobi) |
| STARCH (F) | 3-state Gaussian HMM | none | K clones by Viterbi + ICM | HMRF (Potts-type) | — |
| Alleloscope (F) | 4-state Gaussian HMM (scDNA/ATAC) | in-sample EM: per-SNP phase + per-cell haplotype fraction | hierarchical | none | — |
| CHISEL (F) | scDNA depth ratio + BAF | Eagle2 50-kb blocks; EM phases blocks across all cells | hierarchical | none | — |
| SIGNALS (F, D) | beta-binomial HMM on haplotype blocks (scDNA) | phases blocks by shared gains/losses across cells | UMAP + HDBSCAN | none | R |
| CONGAS / CONGAS+ (A) | Bayesian clone mixture, variational inference | — | mixture | — | GPU backend claimed |

Benchmarks:
- Schmid 2025 (Nat Commun; 6 tools, 21 datasets): no single best tool; reference
  choice matters, external references hurt; losses harder than gains.
- Song 2025 (Brief Bioinform): Numbat best overall; CopyKAT best expression-only.
- Chang 2026 (bioRxiv preprint; 12 tools, 28 datasets): Numbat leads at SNP depth >= 5;
  subclone reconstruction weak for all tools.
- Han 2026 (Nat Commun; spatial, 9 tools): continuous outputs recover more signal,
  discrete calls give fewer false positives; cnLOH, focal events, mixed spots hard.

What Haplodose does differently (design claims, to be benchmarked, not yet shown):
1. No per-cell HMM. Segmentation on pooled clones; cells assigned by NB likelihood
   to piecewise-constant clone dosage profiles.
2. Two noise terms: sampling (shrinks with depth) and bin-level baseline mismatch
   tau (shrinks with genes). Each segment reports its minimum detectable log dosage.
3. Clone differences are tested on the log-rate difference between clones, which
   does not depend on the reference baseline.
4. Cross-fitting: proposals use even bins, difference tests use odd bins.
5. One code path for CPU and GPU (JAX; associative-scan HMM on GPU); Rust streaming
   aggregation for bounded memory.
6. Planned (M2-M3): phase-switch HMM with cross-clone block phasing and per-SNP
   allelic bias from in-sample normals; Potts prior on clone labels for spatial data.

Citations (DOIs as retrieved):
- Gao T 2023 Nat Biotechnol, Numbat, 10.1038/s41587-022-01468-y
- Gao R 2021 Nat Biotechnol, CopyKAT, 10.1038/s41587-020-00795-2
- Serin Harmanci A 2020 Nat Commun, CaSpER, 10.1038/s41467-019-13779-x
- Fan J 2018 Genome Res, HoneyBADGER, 10.1101/gr.228080.117
- De Falco A 2023 Nat Commun, SCEVAN, 10.1038/s41467-023-36790-9
- Huang R 2024 Nat Commun, XClone, 10.1038/s41467-024-51026-0
- Ma C 2024 Nat Methods, CalicoST, 10.1038/s41592-024-02438-9
- Elyanow R 2021 Phys Biol, STARCH, 10.1088/1478-3975/abbe99
- Wu CY 2021 Nat Biotechnol, Alleloscope, 10.1038/s41587-021-00911-w
- Zaccaria S 2021 Nat Biotechnol, CHISEL, 10.1038/s41587-020-0661-6
- Funnell T 2022 Nature, SIGNALS, 10.1038/s41586-022-05249-0
- Wu CY 2025 Nat Methods, Clonalscope, 10.1038/s41592-025-02773-5
- Milite S 2022 Bioinformatics, CONGAS, 10.1093/bioinformatics/btac143
- Patruno L 2023 PLoS Comput Biol, CONGAS+, 10.1371/journal.pcbi.1011557
- Schmid KT 2025 Nat Commun, benchmark, 10.1038/s41467-025-62359-9
- Song M 2025 Brief Bioinform, benchmark, 10.1093/bib/bbaf076
- Chen X 2025 Precis Clin Med, benchmark, 10.1093/pcmedi/pbaf011
- Chang HC 2026 bioRxiv, benchmark, 10.64898/2026.04.12.718050
- Han S 2026 Nat Commun, spatial benchmark, 10.1038/s41467-026-77500-5
