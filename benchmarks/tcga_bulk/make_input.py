"""Haplodose bulk inputs per TCGA cancer type from the PanCanAtlas RNA matrix (read only from ~/CRC).

Expression: EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz. These are NOT raw counts: RSEM
normalised gene values (upper-quartile scaled), batch adjusted across the PanCanAtlas, stored by Xena as
log2(x + 1). Pseudo-counts here are max(2^v - 1, 0). Batch adjustment leaves some values below 0 (clipped to
0, fraction recorded in qc.tsv). Genes with any NA in a cancer type are dropped for that type.
Samples: primary tumours (sample type 01) and solid tissue normals (11) from TCGA_phenotype_denseDataOnlyDownload,
one sample per patient and role. A type is kept when it has >= MIN_NORMALS normals in the matrix; its own
normals are its reference.
Genes: Xena symbols matched to protein-coding chr1-22 genes of Ensembl 112 (GRCh38 positions; names mapping
to more than one gene ID dropped, as in benchmarks/ucan_bulk/make_input.py). Older symbols that no longer
match are dropped (count in qc.tsv).

Raw-count check (--coadread-raw): GDC STAR augmented gene counts (unstranded, GENCODE v36) for COAD and READ,
written as coadread_raw/ beside coadread_xena/ (same samples, Xena pseudo-counts), so the two count sources
can be compared on identical tumours and normals.

Usage: make_input.py OUTDIR [--coadread-raw]
Writes OUTDIR/<TYPE>/{counts.tsv, samples.tsv (sample, role, patient), positions.tsv}, OUTDIR/qc.tsv.
"""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

# Same gene table as the U-CAN bulk benchmark (module shares this file's name, so load it by path).
_spec = importlib.util.spec_from_file_location('ucan_make_input', Path(__file__).resolve().parents[1] / 'ucan_bulk/make_input.py')
_ucan = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ucan)
ensembl_genes = _ucan.genes

EXT = Path('/gpfs/commons/home/jameslee/CRC/data/external')
EXPR = EXT / 'tcga_pancan_xena/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz'
PHENO = EXT / 'tcga_pancan_xena/TCGA_phenotype_denseDataOnlyDownload.tsv.gz'
COADREAD = EXT / 'tcga_coadread'
MIN_NORMALS = 10
# TCGA project codes for the _primary_disease names in the phenotype file.
CODE = {
    'kidney clear cell carcinoma': 'KIRC', 'breast invasive carcinoma': 'BRCA', 'lung squamous cell carcinoma': 'LUSC',
    'lung adenocarcinoma': 'LUAD', 'stomach adenocarcinoma': 'STAD', 'colon adenocarcinoma': 'COAD',
    'head & neck squamous cell carcinoma': 'HNSC', 'prostate adenocarcinoma': 'PRAD', 'thyroid carcinoma': 'THCA',
    'kidney papillary cell carcinoma': 'KIRP', 'liver hepatocellular carcinoma': 'LIHC',
    'uterine corpus endometrioid carcinoma': 'UCEC', 'kidney chromophobe': 'KICH',
    'bladder urothelial carcinoma': 'BLCA', 'esophageal carcinoma': 'ESCA', 'rectum adenocarcinoma': 'READ',
    'ovarian serous cystadenocarcinoma': 'OV', 'pancreatic adenocarcinoma': 'PAAD', 'cholangiocarcinoma': 'CHOL',
}


def write_type(d, X, samples, pos):
    d.mkdir()
    X.index.name = 'gene'
    X.round(3).to_csv(d / 'counts.tsv', sep='\t')
    samples.to_csv(d / 'samples.tsv', sep='\t', index=False)
    pos.loc[pos.gene.isin(X.index)].to_csv(d / 'positions.tsv', sep='\t', index=False)


def star_counts(samples):
    """GDC STAR unstranded counts for the given TCGA sample barcodes (first 15 characters)."""
    fmap = pd.read_csv(COADREAD / 'manifests/file_case_map.tsv', sep='\t')
    fmap = fmap[fmap.data_type.eq('Gene Expression Quantification')].copy()
    fmap['sample'] = fmap.sample_id.str[:15]
    fmap = fmap[fmap['sample'].isin(samples)].drop_duplicates('sample')
    cols = {}
    for r in fmap.itertuples():
        path = COADREAD / 'raw/gene_expression_quantification' / r.file_id / r.file_name
        t = pd.read_csv(path, sep='\t', comment='#', usecols=['gene_id', 'gene_name', 'gene_type', 'unstranded'])
        t = t[t.gene_type.eq('protein_coding')]
        cols[r.sample] = t.groupby('gene_name').unstranded.sum()    # PAR_Y copies share a name with chrX
    return pd.DataFrame(cols).fillna(0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('--coadread-raw', action='store_true')
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)

    g = ensembl_genes()
    pos = g[['gene', 'chromosome', 'start', 'end']]
    ph = pd.read_csv(PHENO, sep='\t', dtype=str)
    V = pd.read_csv(EXPR, sep='\t', index_col=0)
    n_symbols = len(V)
    V = V[~V.index.duplicated(keep=False)]
    V = V.loc[V.index.isin(pos.gene)]
    ph = ph[ph['sample'].isin(V.columns) & ph.sample_type_id.isin(['01', '11'])].copy()
    ph['role'] = ph.sample_type_id.map({'01': 'tumour', '11': 'normal'})
    ph['patient'] = ph['sample'].str[:12]
    ph = ph.drop_duplicates(['patient', 'role'])
    ph['type'] = ph._primary_disease.map(CODE)
    n_norm = ph[ph.role.eq('normal')].groupby('type').size()
    types = sorted(n_norm[n_norm >= MIN_NORMALS].index)
    print('normals per type:', n_norm.to_dict())

    qc = []
    for t in types:
        s = ph[ph.type.eq(t)].sort_values(['role', 'sample'], ascending=[False, True])[['sample', 'role', 'patient']]
        v = V[s['sample']]
        v = v[v.notna().all(axis=1)]
        raw = np.exp2(v) - 1
        neg = float((raw < 0).to_numpy().mean())
        X = raw.clip(lower=0)
        write_type(out / t, X, s, pos)
        tot = X.sum(axis=0)
        qc.append(dict(type=t, n_tumour=int(s.role.eq('tumour').sum()), n_normal=int(s.role.eq('normal').sum()),
                       genes_in_matrix=n_symbols, genes_matched_ensembl=len(V), genes_no_na=len(X),
                       fraction_values_below_zero=neg, median_total_pseudocounts=float(tot.median()),
                       min_total_pseudocounts=float(tot.min())))
        print(qc[-1], flush=True)

    if args.coadread_raw:
        s = ph[ph.type.isin(['COAD', 'READ'])][['sample', 'role', 'patient']]
        R = star_counts(set(s['sample']))
        s = s[s['sample'].isin(R.columns)].sort_values(['role', 'sample'], ascending=[False, True])
        v = V[s['sample']]
        v = v[v.notna().all(axis=1)]
        X = (np.exp2(v) - 1).clip(lower=0)
        R = R.loc[R.index.isin(pos.gene), s['sample']]
        common = X.index.intersection(R.index)
        write_type(out / 'coadread_xena', X.loc[common], s, pos)
        write_type(out / 'coadread_raw', R.loc[common], s, pos)
        # How the pseudo-counts relate to raw counts: per-sample scale and gene-level agreement.
        ratio = (X.loc[common].sum(axis=0) / R.loc[common].sum(axis=0))
        lx, lr = np.log2(X.loc[common] + 1), np.log2(R.loc[common] / R.loc[common].sum(axis=0) * X.loc[common].sum(axis=0) + 1)
        rho = [float(np.corrcoef(lx[c], lr[c])[0, 1]) for c in s['sample']]
        qc.append(dict(type='coadread_raw_vs_xena', n_tumour=int(s.role.eq('tumour').sum()),
                       n_normal=int(s.role.eq('normal').sum()), genes_no_na=len(common),
                       median_pseudo_over_raw_total=float(ratio.median()),
                       median_gene_pearson_log2=float(np.median(rho)), min_gene_pearson_log2=float(np.min(rho))))
        print(qc[-1], flush=True)
    pd.DataFrame(qc).to_csv(out / 'qc.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
