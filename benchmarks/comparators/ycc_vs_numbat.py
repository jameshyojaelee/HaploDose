#!/usr/bin/env python3
"""Agreement of one tool run (Haplodose layout) with native Numbat on a YCC sample. No DNA truth
exists for YCC, so this measures agreement with Numbat, not accuracy.

Usage: ycc_vs_numbat.py TOOL_OUT NUMBAT_DIR SAMPLE DEST
Numbat side: last clone_post_<i>.tsv (compartment_opt, clone_opt) and bulk_clones_final.tsv.gz.
- malignant: cells analysed by both; agreement and Cohen's kappa (tool malignant_call vs tumor).
- clones: ARI of tool clone vs Numbat clone_opt over cells malignant in both.
- arms (tumour clones only, weighted by cells): tool median_log_ratio vs Numbat median gene lnFC,
  Spearman; arm state agreement, where an arm is gain/loss if > 50% of weighted units are
  gain/loss (tool cn_state; Numbat cnv_state_post amp|bamp / del|bdel, loh counted neutral).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, cohen_kappa_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cmp_common as C  # noqa: E402


def arm_state(gain, loss):
    return np.select([gain > 0.5, loss > 0.5], ['gain', 'loss'], 'neutral')


def main(tool_out, numbat_dir, sample, dest):
    tool_out, nb, dest = Path(tool_out), Path(numbat_dir), C.guard(dest)
    if dest.exists() and any(dest.iterdir()):
        raise FileExistsError(f'Output directory is not empty: {dest}')
    dest.mkdir(parents=True, exist_ok=True)       # an empty one left by a failed run is reused
    post_file = sorted(nb.glob('clone_post_*.tsv'), key=lambda p: int(p.stem.split('_')[-1]))[-1]
    post = pd.read_csv(post_file, sep='\t', dtype={'clone_opt': str}).set_index('cell')
    cells = pd.read_csv(tool_out / 'cells.tsv.gz', sep='\t', dtype={'clone': str}).set_index('cell_id')
    cells = cells[cells.clone.notna() & (cells.clone != 'NA')]
    both = cells.index.intersection(post.index)
    t_mal = cells.loc[both, 'malignant_call'].astype(bool).to_numpy()
    n_mal = (post.loc[both, 'compartment_opt'] == 'tumor').to_numpy()
    mb = both[t_mal & n_mal]
    summary = dict(sample=sample, numbat=str(post_file), n_cells_both=len(both),
                   tool_malignant=int(t_mal.sum()), numbat_tumour=int(n_mal.sum()),
                   malignant_agreement=float((t_mal == n_mal).mean()) if len(both) else None,
                   malignant_kappa=float(cohen_kappa_score(t_mal, n_mal)) if len(both) and len(set(t_mal) | set(n_mal)) > 1 else None,
                   n_malignant_both=len(mb),
                   clone_ari=float(adjusted_rand_score(post.loc[mb, 'clone_opt'], cells.loc[mb, 'clone'])) if len(mb) > 1 else None,
                   tool_tumour_clones=int(cells.loc[cells.malignant_call.astype(bool), 'clone'].nunique()),
                   numbat_tumour_clones=int(post.loc[post.compartment_opt == 'tumor', 'clone_opt'].nunique()))

    # Tool arms: tumour clones, weighted by MAP cells.
    clones = pd.read_csv(tool_out / 'clones.tsv', sep='\t')
    tum = clones[~clones.normal_like].set_index('clone').n_cells_map
    arms = pd.read_csv(tool_out / 'arms.tsv', sep='\t')
    arms = arms[arms.clone.isin(tum.index)].assign(w=lambda a: a.clone.map(tum).astype(float))
    if tum.empty:   # no tumour clone called: malignant metrics stand, arm metrics are NA
        summary.update(n_arms=0, rho_arms=None, arm_state_agreement=None, note='tool called no tumour clone')
        (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        print(json.dumps(summary))
        return
    tool =arms.groupby('arm').apply(lambda g: np.average(g.median_log_ratio, weights=g.w), include_groups=False).rename('tool')
    prof = pd.read_csv(tool_out / 'profiles.tsv.gz', sep='\t')
    bins = pd.read_csv(tool_out / 'bins.tsv', sep='\t')
    bins['arm'] = C.arm_of(bins.chromosome, bins.start)   # Haplodose's bins.tsv has no arm column
    prof = prof[prof.clone.isin(tum.index)].merge(bins[['bin', 'arm']], on='bin')
    prof['w'] = prof.clone.map(tum).astype(float)
    f = prof.assign(g=prof.cn_state.ge(3) * prof.w, l=prof.cn_state.eq(1) * prof.w).groupby('arm')[['g', 'l', 'w']].sum()
    tool_state = pd.Series(arm_state(f.g / f.w, f.l / f.w), index=f.index, name='tool_state')

    # Numbat arms: tumour clones from clone_post, weighted by cells; gene-level lnFC and states.
    ncl = post[post.compartment_opt == 'tumor'].clone_opt.value_counts()
    bulk = pd.read_csv(nb / 'bulk_clones_final.tsv.gz', sep='\t', low_memory=False)
    g = bulk[bulk.gene.notna() & (bulk.gene.astype(str) != '') & bulk.lnFC.notna()].drop_duplicates(['sample', 'gene'])
    g = g[g['sample'].astype(str).isin(ncl.index)]
    start = g.gene_start.fillna(g.POS) if 'gene_start' in g else g.POS
    g = g.assign(arm=C.arm_of('chr' + g.CHROM.astype(str).str.replace('chr', ''), start),
                 w=g['sample'].astype(str).map(ncl).astype(float),
                 gain=g.cnv_state_post.isin(['amp', 'bamp']), loss=g.cnv_state_post.isin(['del', 'bdel']))
    per_clone = g.groupby(['sample', 'arm']).agg(lnfc=('lnFC', 'median'), w=('w', 'first')).reset_index()
    numbat = per_clone.groupby('arm').apply(lambda x: np.average(x.lnfc, weights=x.w), include_groups=False).rename('numbat')
    ns = g.assign(gw=g.gain * g.w, lw=g.loss * g.w).groupby('arm')[['gw', 'lw', 'w']].sum()
    numbat_state = pd.Series(arm_state(ns.gw / ns.w, ns.lw / ns.w), index=ns.index, name='numbat_state')

    table = pd.concat([tool, tool_state, numbat, numbat_state], axis=1, join='inner').reset_index(names='arm')
    table.to_csv(dest / 'arms.tsv', sep='\t', index=False)
    altered = table.numbat_state.ne('neutral')
    summary.update(n_arms=len(table),
                   rho_arms=float(spearmanr(table.tool, table.numbat).statistic) if len(table) > 3 else None,
                   arm_state_agreement=float((table.tool_state == table.numbat_state).mean()) if len(table) else None,
                   numbat_altered_arms=int(altered.sum()),
                   altered_arms_same_state=int((table.tool_state.eq(table.numbat_state) & altered).sum()),
                   numbat_neutral_arms_called_altered=int((table.tool_state.ne('neutral') & ~altered).sum()))
    (dest / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary))


if __name__ == '__main__':
    main(*sys.argv[1:])
