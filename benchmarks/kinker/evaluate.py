"""Score Haplodose on the Kinker 2020 rotations against DepMap 24Q4 CN, with a pseudobulk baseline.

Per query line (one single-line specimen in one rotation):
- truth per gene: DepMap OmicsCNGene of the line minus the UMI-weighted mean over that pool's fit reference
  lines (the pooled reference is UMI-weighted); arm truth = median over the run's used genes in the arm.
- Haplodose: clone profiles (log_ratio, log_dose and, when the run has it, log_ratio_gene) weighted by the
  specimen's cells per clone.
- pseudobulk: line counts summed vs the pool's fit reference counts summed, per-gene log((y+0.5)/(E+0.5)),
  same genes as Haplodose, median-centred, arm median. It is the scRNA analogue of the bulk naive method.
Per line: Spearman across arms (arms with >= 5 genes); on DNA-informative arms (|truth| >= 0.2) the share
with the right sign. Methods compared per line (unit = cell line) by paired Wilcoxon, Holm over pairs.
Also per specimen: tumour clones (1 expected) and the fraction of cells called malignant.
Usage: evaluate.py [--tool NAME] INPUT_DIR DEST RUN_DIR...  (NAME labels the run columns; default haplodose)
"""
import itertools
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests

CEN = Path(__file__).resolve().parents[2] / 'python/haplodose/data/centromeres_hg38.tsv'
INFORMATIVE = 0.2
MIN_GENES = 5


def load_counts(path):
    with h5py.File(path) as h:
        X = sparse.csr_matrix((h['X/data'][:], h['X/indices'][:], h['X/indptr'][:]))
        genes = h['genes'][:].astype(str)
        cells = h['cell_ids'][:].astype(str)
    return X, genes, cells


def score_run(inp, run, X, genes, cells, cn, cen, tool):
    meta = pd.read_csv(run / 'cells.tsv.gz', sep='\t')
    used = pd.read_csv(run / 'genes.tsv.gz', sep='\t')
    bins = pd.read_csv(run / 'bins.tsv', sep='\t')
    prof = pd.read_csv(run / 'profiles.tsv.gz', sep='\t')
    clones = pd.read_csv(run / 'clones.tsv', sep='\t')
    pos = pd.read_csv(inp / 'positions.tsv', sep='\t').drop_duplicates('gene').set_index('gene')
    used = used[used.gene.isin(pos.index) & used.gene.isin(cn.index)].copy()
    used['arm'] = ('chr' + pos.loc[used.gene, 'chromosome'].astype(str).str.removeprefix('chr')).to_numpy()
    used['arm'] = used.arm + np.where(pos.loc[used.gene, 'start'].to_numpy() < used.arm.map(cen).to_numpy(), 'p', 'q')
    bins_arm = used.groupby('bin').arm.agg(lambda a: a.mode().iloc[0])
    gidx = pd.Index(genes).get_indexer(used.gene)
    row = pd.Series(np.arange(len(cells)), index=cells)
    rows, spec_rows = [], []
    for spec, c in meta[meta.design.eq('single_line')].groupby('specimen'):
        line, ds = c.truth_line.iloc[0], c.dataset.iloc[0]
        ref = meta[meta.dataset.eq(ds) & meta.role.eq('reference') & meta.ref_split.eq('fit')]
        R = X[row[ref.cell_id].to_numpy()][:, gidx]
        w = ref.groupby('line').umi.sum()
        truth_g = cn.loc[used.gene, line].to_numpy() - (cn.loc[used.gene, w.index].to_numpy() @ (w / w.sum()).to_numpy())
        # Pseudobulk baseline.
        y = np.asarray(X[row[c.cell_id].to_numpy()][:, gidx].sum(axis=0)).ravel()
        r = np.asarray(R.sum(axis=0)).ravel()
        E = r / r.sum() * y.sum()
        lr = np.log((y + 0.5) / (E + 0.5))
        g = pd.DataFrame({'arm': used.arm.to_numpy(), 'truth': truth_g, 'pseudobulk': lr - np.median(lr)})
        # Haplodose: clone profiles weighted by cells per clone.
        n = c.clone.value_counts()
        p = prof[prof.specimen.eq(spec) & prof.clone.isin(n.index)]
        for col in ['log_ratio', 'log_dose'] + (['log_ratio_gene'] if 'log_ratio_gene' in p.columns else []):
            pv = p.pivot(index='bin', columns='clone', values=col)
            nn = n[n.index.isin(pv.columns)]                 # clones without a profile are left out
            b = (pv[nn.index] * nn.to_numpy()).sum(axis=1) / nn.sum()
            g[f'{tool}_{col}'] = used.bin.map(b).to_numpy()
        a = g.groupby('arm').agg(['median', 'size'])
        keep = a[('truth', 'size')] >= MIN_GENES
        a = a.loc[keep].xs('median', axis=1, level=1)
        d = dict(specimen=spec, line=line, pool=ds, n_cells=len(c), n_arms=len(a),
                 informative_arms=int((a.truth.abs() >= INFORMATIVE).sum()))
        inf = a.truth.abs() >= INFORMATIVE
        for m in [f'{tool}_{x}' for x in ('log_ratio', 'log_dose', 'log_ratio_gene') if f'{tool}_{x}' in g] + ['pseudobulk']:
            d[f'rho_{m}'] = spearmanr(a[m], a.truth)[0]
            d[f'sign_{m}'] = float((np.sign(a.loc[inf, m]) == np.sign(a.loc[inf, 'truth'])).mean()) if inf.any() else np.nan
        rows.append(d)
        tum = clones[clones.specimen.eq(spec) & ~clones.normal_like]
        spec_rows.append(dict(specimen=spec, line=line, tumour_clones=len(tum),
                              malignant_fraction=float(c.malignant_call.mean())))
    return pd.DataFrame(rows), pd.DataFrame(spec_rows)


def main(inp, dest, *runs, tool='haplodose'):
    inp, dest = Path(inp), Path(dest)
    dest.mkdir(parents=True, exist_ok=False)
    cen = pd.read_csv(CEN, sep='\t', comment='#').set_index('chromosome').p_arm_end
    cn = pd.read_csv(inp / 'depmap_cn.tsv.gz', sep='\t', index_col=0)
    cn = cn[~cn.index.duplicated()]
    X, genes, cells = load_counts(inp / 'counts.h5')
    parts = [score_run(inp, Path(r), X, genes, cells, cn, cen, tool) for r in runs]
    lines = pd.concat([p[0] for p in parts]).reset_index(drop=True)
    specs = pd.concat([p[1] for p in parts]).reset_index(drop=True)
    lines.to_csv(dest / 'line_scores.tsv', sep='\t', index=False)
    specs.to_csv(dest / 'specimens.tsv', sep='\t', index=False)
    methods = [m for m in (f'{tool}_log_ratio', f'{tool}_log_dose', f'{tool}_log_ratio_gene', 'pseudobulk')
               if f'rho_{m}' in lines]
    tests = []
    for a, b in itertools.combinations(methods, 2):
        x = lines[[f'rho_{a}', f'rho_{b}']].dropna()
        tests.append(dict(a=a, b=b, n=len(x), median_a=x.iloc[:, 0].median(), median_b=x.iloc[:, 1].median(),
                          a_higher=int((x.iloc[:, 0] > x.iloc[:, 1]).sum()), p=wilcoxon(x.iloc[:, 0], x.iloc[:, 1]).pvalue))
    tests = pd.DataFrame(tests)
    tests['p_holm'] = multipletests(tests.p, method='holm')[1]
    tests.to_csv(dest / 'paired_tests.tsv', sep='\t', index=False)
    summary = dict(n_lines=int(len(lines)), median_informative_arms=float(lines.informative_arms.median()),
                   median_rho={m: float(lines[f'rho_{m}'].median()) for m in methods},
                   mean_sign_informative={m: float(lines[f'sign_{m}'].mean()) for m in methods},
                   tumour_clones_distribution=specs.tumour_clones.value_counts().sort_index().to_dict(),
                   median_malignant_fraction=float(specs.malignant_fraction.median()),
                   paired_tests=tests.to_dict(orient='records'),
                   correction='Holm over method pairs; unit = cell line')
    (dest / 'summary.json').write_text(json.dumps(summary, indent=1, default=float))
    print(json.dumps(summary, indent=1, default=float))


if __name__ == '__main__':
    args = sys.argv[1:]
    tool = 'haplodose'
    if args[0] == '--tool':
        tool, args = args[1], args[2:]
    main(*args, tool=tool)
