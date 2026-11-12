"""Haplodose bulk inputs for the CCLE benchmark (prespecified SOTA claim 3, docs/STATUS.md), one directory per
lineage group, with GTEx v10 tissue-matched normals (user decision 2026-11-02).

Tumours: DepMap 24Q4 cancer cell lines with a default RNA profile (RNA-SeQC gene counts) and a default DNA profile
that has PureCN absolute segments; OncotreePrimaryDisease 'Non-Cancerous' dropped. Sample name = ModelID.
Normals: GTEx v10 RNA-SeQC 2.4.2 gene reads of the tissue(s) in GROUPS, RNA-seq freeze samples only, at most
MAX_NORMALS per group drawn with a fixed seed. Both matrices are RNA-SeQC counts, matched by unversioned Ensembl
gene ID to protein-coding chr1-22 genes of Ensembl 112 (benchmarks/ucan_bulk/make_input.py gene table).
Groups: the lineage map below, fixed before any CCLE or GTEx expression value was read. A lineage with no GTEx
organ match is excluded; a group with < MIN_LINES lines is excluded (recorded in qc.tsv).
Also writes truth for scoring: absolute_segments.tsv.gz (Sample, chromosome, Start, End, Modal_Total_CN; hg38)
and absolute_ploidy.tsv (ploidy = length-weighted mean PureCN total CN over chr1-22 segments), plus
dev_lines.tsv (ModelIDs whose DepMap CN was used as truth in development: Kinker, MIX-seq).

Usage: make_input.py OUTDIR
Writes OUTDIR/<GROUP>/{counts.tsv, samples.tsv (sample, role, patient), positions.tsv}, OUTDIR/qc.tsv, group_map.tsv.
"""
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H = Path('/gpfs/commons/home/jameslee/Haplodose')
DM = H / 'data/external/depmap_24q4_ccle'
GT = H / 'data/external/gtex_v10'
GTEX_COUNTS = GT / 'GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz'
GTEX_ATTR = GT / 'GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt'
KINKER_LINES = H / 'results/kinker_input_v2_20261101_1551/eligible_lines.tsv'
MIXSEQ_LINES = H / 'results/heldout_mixseq_input_20261102_1159/lines.tsv'
MAX_NORMALS = 150
MIN_LINES = 10
SEED = 20261102

# group: (OncotreeLineage, OncotreePrimaryDisease or None for all, GTEx SMTSD tissues)
GROUPS = {
    'LUNG': ('Lung', None, ['Lung']),
    'LYMPHOID': ('Lymphoid', None, ['Cells - EBV-transformed lymphocytes']),
    'MYELOID': ('Myeloid', None, ['Whole Blood']),
    'SKIN': ('Skin', None, ['Skin - Sun Exposed (Lower leg)', 'Skin - Not Sun Exposed (Suprapubic)']),
    'BRAIN': ('CNS/Brain', None, ['Brain - Cortex', 'Brain - Frontal Cortex (BA9)']),
    'HEADNECK': ('Head and Neck', None, ['Minor Salivary Gland']),
    'BOWEL': ('Bowel', None, ['Colon - Transverse']),
    'ESCC': ('Esophagus/Stomach', 'Esophageal Squamous Cell Carcinoma', ['Esophagus - Mucosa']),
    'EGA': ('Esophagus/Stomach', 'Esophagogastric Adenocarcinoma', ['Stomach']),
    'BREAST': ('Breast', None, ['Breast - Mammary Tissue']),
    'OVARY': ('Ovary/Fallopian Tube', None, ['Ovary', 'Fallopian Tube']),
    'PANCREAS': ('Pancreas', None, ['Pancreas']),
    'KIDNEY': ('Kidney', None, ['Kidney - Cortex']),
    'BLADDER': ('Bladder/Urinary Tract', None, ['Bladder']),
    'UTERUS': ('Uterus', None, ['Uterus']),
    'LIVER': ('Liver', None, ['Liver']),
    'THYROID': ('Thyroid', None, ['Thyroid']),
    'CERVIX': ('Cervix', None, ['Cervix - Ectocervix', 'Cervix - Endocervix']),
    'PROSTATE': ('Prostate', None, ['Prostate']),
    'TESTIS': ('Testis', None, ['Testis']),
}
# Excluded lineages (no GTEx organ match or not cancer): Soft Tissue, Bone, Peripheral Nervous System, Fibroblast,
# Biliary Tract, Eye, Pleura, Ampulla of Vater, Adrenal Gland, Vulva/Vagina, Other.

_spec = importlib.util.spec_from_file_location('ucan_make_input', H / 'benchmarks/ucan_bulk/make_input.py')
_ucan = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ucan)


def truth(dna_profiles):
    seg = pd.read_csv(DM / 'OmicsAbsoluteCNSegmentsProfile.csv', dtype={'Chromosome': str})
    seg = seg[seg.ProfileID.isin(dna_profiles.index) & seg.Chromosome.isin([str(i) for i in range(1, 23)])]
    seg = seg[seg.SegmentAbsoluteCN.notna()].copy()
    seg['Sample'] = dna_profiles.loc[seg.ProfileID].to_numpy()
    seg['chromosome'] = 'chr' + seg.Chromosome
    seg = seg.rename(columns={'SegmentAbsoluteCN': 'Modal_Total_CN'})
    seg = seg[['Sample', 'chromosome', 'Start', 'End', 'Modal_Total_CN']].sort_values(['Sample', 'chromosome', 'Start'])
    w = (seg.End - seg.Start).astype(float)
    ploidy = (seg.Modal_Total_CN * w).groupby(seg.Sample).sum() / w.groupby(seg.Sample).sum()
    return seg, ploidy.rename('ploidy')


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    g = _ucan.genes()
    g['gene_id'] = g.gene_id.str.split('.').str[0]
    g = g[~g.gene_id.duplicated(keep=False)].set_index('gene_id')
    pos = g.reset_index()[['gene', 'chromosome', 'start', 'end']]

    model = pd.read_csv(DM / 'Model.csv').set_index('ModelID')
    dp = pd.read_csv(DM / 'OmicsDefaultModelProfiles.csv')
    rna = dp[dp.ProfileType.eq('rna')].drop_duplicates('ModelID').set_index('ProfileID').ModelID
    dna_all = dp[dp.ProfileType.eq('dna')].drop_duplicates('ModelID').set_index('ProfileID').ModelID
    seg, ploidy = truth(dna_all)
    with_cn = set(ploidy.index)

    X = pd.read_csv(DM / 'OmicsExpressionRNASeQCGeneCountProfile.csv', index_col=0)
    X = X.loc[X.index.isin(rna.index)]
    X.index = rna.loc[X.index].to_numpy()
    ens = X.columns.str.extract(r'\((ENSG\d+)\)$')[0].to_numpy()
    X.columns = ens
    X = X.loc[:, pd.notna(ens)]
    X = X.T.groupby(level=0).sum()          # genes x lines; a repeated ID would be summed (none expected)
    n_ccle_genes = len(X)
    X = X.loc[X.index.isin(g.index)]

    attr = pd.read_csv(GTEX_ATTR, sep='\t', usecols=['SAMPID', 'SMTSD', 'SMAFRZE'])
    attr = attr[attr.SMAFRZE.eq('RNASEQ')]
    rng = np.random.default_rng(SEED)
    chosen = {}
    for grp, (_, _, tissues) in GROUPS.items():
        ids = np.sort(attr.SAMPID[attr.SMTSD.isin(tissues)].to_numpy())
        if len(ids) > MAX_NORMALS:
            ids = np.sort(rng.choice(ids, MAX_NORMALS, replace=False))
        chosen[grp] = ids
    with gzip.open(GTEX_COUNTS, 'rt') as fh:
        fh.readline(); fh.readline()
        header = fh.readline().rstrip('\n').split('\t')
    want = set(np.concatenate(list(chosen.values())))
    cols = [c for c in header if c in want]
    N = pd.read_csv(GTEX_COUNTS, sep='\t', skiprows=2, usecols=['Name'] + cols, index_col=0)
    N.index = N.index.str.split('.').str[0]
    N = N[~N.index.str.endswith('_PAR_Y')].groupby(level=0).sum()
    common = X.index.intersection(N.index)
    sym = g.gene.loc[common]

    dev = set(pd.read_csv(KINKER_LINES, sep='\t').query('eligible').ModelID.dropna())
    mix = pd.read_csv(MIXSEQ_LINES, sep='\t')
    dev |= set(model.index[model.CCLEName.isin(mix.line)])
    pd.DataFrame({'ModelID': sorted(dev)}).to_csv(out / 'dev_lines.tsv', sep='\t', index=False)

    qc, gmap = [], []
    for grp, (lineage, disease, tissues) in GROUPS.items():
        m = model[model.OncotreeLineage.eq(lineage) & model.OncotreePrimaryDisease.ne('Non-Cancerous')]
        if disease:
            m = m[m.OncotreePrimaryDisease.eq(disease)]
        lines = sorted(set(m.index) & set(X.columns) & with_cn)
        normals = [c for c in chosen[grp] if c in N.columns]
        gmap.append(dict(group=grp, lineage=lineage, disease=disease or 'all', gtex_tissues='; '.join(tissues),
                         n_lines=len(lines), n_normals=len(normals), included=len(lines) >= MIN_LINES))
        if len(lines) < MIN_LINES:
            continue
        C = pd.concat([X.loc[common, lines], N.loc[common, normals]], axis=1)
        C.index = sym.to_numpy()
        C.index.name = 'gene'
        s = pd.DataFrame({'sample': lines + normals, 'role': ['tumour'] * len(lines) + ['normal'] * len(normals)})
        s['patient'] = s['sample']
        d = out / grp
        d.mkdir()
        C.round(0).astype(np.int64).to_csv(d / 'counts.tsv', sep='\t')
        s.to_csv(d / 'samples.tsv', sep='\t', index=False)
        pos[pos.gene.isin(C.index)].to_csv(d / 'positions.tsv', sep='\t', index=False)
        tot = C.sum(axis=0)
        qc.append(dict(group=grp, n_tumour=len(lines), n_normal=len(normals), n_genes=len(C),
                       n_dev_lines=len(set(lines) & dev), median_total_tumour=float(tot[lines].median()),
                       median_total_normal=float(tot[normals].median()), min_total=float(tot.min())))
        print(qc[-1], flush=True)
    pd.DataFrame(gmap).to_csv(out / 'group_map.tsv', sep='\t', index=False)
    pd.DataFrame(qc).to_csv(out / 'qc.tsv', sep='\t', index=False)
    seg[seg.Sample.isin(X.columns)].to_csv(out / 'absolute_segments.tsv.gz', sep='\t', index=False)
    ploidy.to_frame().to_csv(out / 'absolute_ploidy.tsv', sep='\t')
    (out / 'params.json').write_text(json.dumps(dict(max_normals=MAX_NORMALS, min_lines=MIN_LINES, seed=SEED,
                                                     ccle_genes=n_ccle_genes, genes_common=len(common)), indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
