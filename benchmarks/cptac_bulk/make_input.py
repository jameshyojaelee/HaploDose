"""Haplodose bulk inputs for the CPTAC-3 benchmark (prespecified claim 6, docs/STATUS.md "Claim 6 prespecified
2026-11-03"), one directory per cohort, with CPTAC normal-adjacent tissue as the normal panel.

Data: DATA/files.tsv (GDC file metadata written by gdc_manifest.py before download) and the files under
DATA/<kind>/<file_id>/<file_name> (kind = rna | ascat), md5-checked by fetch.sbatch.

Rules (as written in the claim):
  cohort      from GDC primary_site + disease_type (cohort() below; fixed from the metadata survey,
              dev/cptac_survey_20261103_0612, before any data file was read). AML and CPTAC-2 are not cohorts.
  truth       AscatNGS "Allele-specific Copy Number Segment" (GRCh38), Copy_Number = total CN; ploidy =
              length-weighted mean total CN over chr1-22 (the claim 3 rule). One ASCAT file per case; if a case
              has several, the smallest file_id [not stated in the claim; deterministic choice, recorded in qc].
              ASCAT tumour samples = the file's samples with tissue_type Tumor.
  tumour RNA  STAR - Counts unstranded column, files whose samples are all Primary Tumor. Per case: the files whose
              sample list contains an ASCAT tumour sample; several -> smallest file_id; none -> case excluded.
              Pooled files (2-5 tumour samples) are kept and flagged (samples.tsv column pooled; sensitivity).
  normals     every Solid Tissue Normal STAR - Counts file of the cohort, one per case (smallest file_id),
              matched and normal-only cases alike. Cohort excluded if < MIN_NORMALS normals.
Genes: protein-coding chr1-22 genes of Ensembl 112 (benchmarks/ucan_bulk/make_input.py gene table), matched by
unversioned Ensembl gene ID to the GENCODE v36 STAR table; PAR_Y rows dropped.
Sample names: tumour = case submitter_id; normal = case submitter_id + '_N'; patient = case submitter_id.

Usage: make_input.py DATA OUTDIR
Writes OUTDIR/<COHORT>/{counts.tsv, samples.tsv (sample, role, patient, file_id, pooled), positions.tsv},
OUTDIR/{qc.tsv, cohort_map.tsv, absolute_segments.tsv.gz, absolute_ploidy.tsv, selection.tsv, params.json}.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

H = Path('/gpfs/commons/home/jameslee/Haplodose')
MIN_NORMALS = 10
COHORTS = ['CCRCC', 'UCEC', 'LUAD', 'GBM', 'CM', 'STAD', 'PDAC', 'HNSCC', 'LSCC']
HNSCC_SITES = {'Larynx', 'Other and unspecified parts of tongue', 'Floor of mouth', 'Base of tongue',
               'Other and ill-defined sites in lip, oral cavity and pharynx', 'Oropharynx', 'Lip', 'Gum', 'Tonsil',
               'Other and unspecified parts of mouth', 'Other and ill-defined sites'}

_spec = importlib.util.spec_from_file_location('ucan_make_input', H / 'benchmarks/ucan_bulk/make_input.py')
_ucan = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ucan)


def cohort(site, disease):
    """CPTAC-3 cohort from GDC case primary_site and disease_type; None = not in claim 6."""
    d = (disease or '').lower()
    if site == 'Kidney' and d == 'adenomas and adenocarcinomas':
        return 'CCRCC'
    if site == 'Uterus, NOS' and d == 'adenomas and adenocarcinomas':
        return 'UCEC'
    if site == 'Bronchus and lung' and d == 'adenomas and adenocarcinomas':
        return 'LUAD'
    if site == 'Bronchus and lung' and d == 'squamous cell neoplasms':
        return 'LSCC'
    if site == 'Brain' and d in ('gliomas', 'not applicable'):
        return 'GBM'
    if site == 'Skin' and d == 'nevi and melanomas':
        return 'CM'
    if site == 'Stomach' and d == 'adenomas and adenocarcinomas':
        return 'STAD'
    if site == 'Pancreas' and d in ('ductal and lobular neoplasms', 'not applicable'):
        return 'PDAC'
    if site in HNSCC_SITES and d == 'squamous cell neoplasms':
        return 'HNSCC'
    return None


def read_star(path):
    """GDC augmented STAR gene counts -> unstranded counts by unversioned Ensembl ID (N_* rows, PAR_Y dropped)."""
    t = pd.read_csv(path, sep='\t', comment='#', usecols=['gene_id', 'unstranded'])
    t = t[t.gene_id.str.startswith('ENSG') & ~t.gene_id.str.endswith('_PAR_Y')]
    t['gene_id'] = t.gene_id.str.split('.').str[0]
    return t.groupby('gene_id').unstranded.sum()


def read_ascat(path, sample):
    t = pd.read_csv(path, sep='\t')
    t['chromosome'] = 'chr' + t.Chromosome.astype(str).str.replace('^chr', '', regex=True)
    t = t[t.chromosome.isin([f'chr{i}' for i in range(1, 23)]) & t.Copy_Number.notna()]
    return pd.DataFrame({'Sample': sample, 'chromosome': t.chromosome, 'Start': t.Start, 'End': t.End,
                         'Modal_Total_CN': t.Copy_Number.astype(float)})


def select(meta):
    """Apply the claim 6 file rules to metadata only. Returns (tumour rows, normal rows, case-level qc)."""
    meta = meta.copy()
    meta['cohort'] = [cohort(s, d) for s, d in zip(meta.primary_site, meta.disease_type)]
    meta = meta[meta.cohort.notna()]
    asc = meta[meta.kind.eq('ascat')].sort_values('file_id')
    n_ascat = asc.groupby('case').size()
    asc = asc.drop_duplicates('case')
    rna = meta[meta.kind.eq('rna')]
    types = rna.sample_types.str.split(';')
    tum_rna = rna[types.map(lambda x: all(v == 'Primary Tumor' for v in x))]
    nor_rna = rna[types.map(lambda x: all(v == 'Solid Tissue Normal' for v in x))]
    rows, qc = [], []
    for _, a in asc.iterrows():
        tum_samples = set(a.tumour_samples.split(';')) if a.tumour_samples else set()
        c = tum_rna[tum_rna.case.eq(a.case)]
        hit = np.array([bool(tum_samples & set(x.split(';'))) for x in c.sample_ids], dtype=bool)
        c = c[hit].sort_values('file_id')
        qc.append(dict(case=a.case, cohort=a.cohort, n_ascat_files=int(n_ascat[a.case]), n_rna_candidates=len(c),
                       included=len(c) > 0))
        if len(c):
            r = c.iloc[0]
            rows.append(dict(case=a.case, cohort=a.cohort, rna_file_id=r.file_id, rna_path=r.path,
                             pooled=len(r.sample_ids.split(';')) > 1, ascat_file_id=a.file_id, ascat_path=a.path))
    tum = pd.DataFrame(rows)
    nor = nor_rna.sort_values('file_id').drop_duplicates('case')[['case', 'cohort', 'file_id', 'path']]
    return tum, nor, pd.DataFrame(qc)


def main(data, out):
    data, out = Path(data), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    meta = pd.read_csv(data / 'files.tsv', sep='\t', dtype=str, keep_default_na=False)
    meta['path'] = [str(data / k / f / n) for k, f, n in zip(meta.kind, meta.file_id, meta.file_name)]
    tum, nor, caseqc = select(meta)
    caseqc.to_csv(out / 'selection.tsv', sep='\t', index=False)

    g = _ucan.genes()
    g['gene_id'] = g.gene_id.str.split('.').str[0]
    g = g[~g.gene_id.duplicated(keep=False)].set_index('gene_id')
    pos = g.reset_index()[['gene', 'chromosome', 'start', 'end']]

    segs = [read_ascat(p, c) for c, p in zip(tum.case, tum.ascat_path)]
    seg = pd.concat(segs).sort_values(['Sample', 'chromosome', 'Start']) if segs else pd.DataFrame()
    w = (seg.End - seg.Start).astype(float)
    ploidy = ((seg.Modal_Total_CN * w).groupby(seg.Sample).sum() / w.groupby(seg.Sample).sum()).rename('ploidy')

    qc, cmap = [], []
    for co in COHORTS:
        t = tum[tum.cohort.eq(co)].sort_values('case')
        n = nor[nor.cohort.eq(co)].sort_values('case')
        cmap.append(dict(cohort=co, n_tumour=len(t), n_normal=len(n), n_pooled=int(t.pooled.sum()) if len(t) else 0,
                         included=len(n) >= MIN_NORMALS and len(t) > 0))
        if not cmap[-1]['included']:
            continue
        cols = {c: read_star(p) for c, p in zip(t.case, t.rna_path)}
        cols.update({f'{c}_N': read_star(p) for c, p in zip(n.case, n.path)})
        C = pd.DataFrame(cols).fillna(0)
        C = C.loc[C.index.isin(g.index)]
        C.index = g.gene.loc[C.index].to_numpy()
        C.index.name = 'gene'
        s = pd.DataFrame({'sample': list(t.case) + [f'{c}_N' for c in n.case],
                          'role': ['tumour'] * len(t) + ['normal'] * len(n),
                          'patient': list(t.case) + list(n.case),
                          'file_id': list(t.rna_file_id) + list(n.file_id),
                          'pooled': list(t.pooled) + [False] * len(n)})
        d = out / co
        d.mkdir()
        C.round(0).astype(np.int64).to_csv(d / 'counts.tsv', sep='\t')
        s.to_csv(d / 'samples.tsv', sep='\t', index=False)
        pos[pos.gene.isin(C.index)].to_csv(d / 'positions.tsv', sep='\t', index=False)
        tot = C.sum(axis=0)
        qc.append(dict(cohort=co, n_tumour=len(t), n_normal=len(n), n_matched_normal=int(n.case.isin(t.case).sum()),
                       n_pooled=int(t.pooled.sum()), n_genes=len(C),
                       median_total_tumour=float(tot[list(t.case)].median()),
                       median_total_normal=float(tot[[f'{c}_N' for c in n.case]].median()), min_total=float(tot.min())))
        print(qc[-1], flush=True)
    pd.DataFrame(cmap).to_csv(out / 'cohort_map.tsv', sep='\t', index=False)
    pd.DataFrame(qc).to_csv(out / 'qc.tsv', sep='\t', index=False)
    seg.to_csv(out / 'absolute_segments.tsv.gz', sep='\t', index=False)
    ploidy.to_frame().to_csv(out / 'absolute_ploidy.tsv', sep='\t')
    pd.DataFrame({'case': tum.case[tum.pooled]}).to_csv(out / 'pooled_cases.tsv', sep='\t', index=False)
    (out / 'params.json').write_text(json.dumps(dict(min_normals=MIN_NORMALS, cohorts=COHORTS), indent=1))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
