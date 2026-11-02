"""Build Haplodose bulk inputs from U-CAN tumour and normal RNA (read only from ~/CRC).

Counts: ~/CRC/data/interim/full_transcriptome_counts/exon_counts.tsv.gz (featureCounts exon counts, Ensembl
gene IDs, samples CRC.SW.<patient>.<T|N>). Genes: protein-coding genes on chr1-22 from Ensembl 112, keyed by
gene name; names that map to more than one gene ID are dropped. T2 (second tumour) samples are dropped. Patients are U#### and UM####.
Usage: make_input.py OUTDIR
Writes counts.tsv (gene names x samples), samples.tsv (sample, role, patient), positions.tsv.
"""
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

CRC = Path('/gpfs/commons/home/jameslee/CRC')
COUNTS = CRC / 'data/interim/full_transcriptome_counts/exon_counts.tsv.gz'
GTF = CRC / 'data/external/ensembl112/Homo_sapiens.GRCh38.112.gtf.gz'


def genes():
    rows = []
    with gzip.open(GTF, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.split('\t', 9)
            if f[2] != 'gene':
                continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
            if a.get('gene_biotype') != 'protein_coding' or 'gene_name' not in a:
                continue
            rows.append(dict(gene_id=a['gene_id'], gene=a['gene_name'], chromosome=f[0], start=int(f[3]), end=int(f[4])))
    g = pd.DataFrame(rows)
    g = g[g.chromosome.isin([str(c) for c in range(1, 23)])]
    return g[~g.gene.duplicated(keep=False)]


def main(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    g = genes()
    X = pd.read_csv(COUNTS, sep='\t', index_col=0)
    X = X.drop(columns=[c for c in X.columns if not re.fullmatch(r'CRC\.SW\.UM?\d+\.[TN]', c)])
    X = X.loc[X.index.intersection(g.gene_id)]
    name = g.set_index('gene_id').gene
    X.index = name.loc[X.index].to_numpy()
    X.index.name = 'gene'
    X.to_csv(out / 'counts.tsv', sep='\t')
    s = pd.DataFrame({'sample': X.columns})
    s['role'] = s['sample'].str[-1].map({'T': 'tumour', 'N': 'normal'})
    s['patient'] = s['sample'].str.extract(r'(UM?\d+)')[0]
    s.to_csv(out / 'samples.tsv', sep='\t', index=False)
    g[['gene', 'chromosome', 'start', 'end']].to_csv(out / 'positions.tsv', sep='\t', index=False)
    print(len(X), 'genes', s.role.value_counts().to_dict())


if __name__ == '__main__':
    main(sys.argv[1])
