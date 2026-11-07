"""Tahoe-100M plate1 DMSO_TF input in the Kinker / MIX-seq benchmark format (claim 8, docs/STATUS.md).

Cells: obs_metadata rows of plate1 with drug == 'DMSO_TF' and the pass_filter value recorded in
EXTRACT_DIR/selection_rule.txt (benchmarks/tahoe/fetch_tahoe.py). Line = obs cell_name (Tahoe SNP demultiplexing);
DepMap ID from cell_line_metadata Cell_ID_DepMap. One pool, 'plate1' (both DMSO_TF samples together).
Lines: >= 50 cells and DepMap 24Q4 OmicsCNGene. Each such line keeps at most 300 cells, drawn without replacement
with numpy default_rng(20261106), lines in sorted order, cells sorted by cell_id. Only these cells are written.
Design as benchmarks/heldout/prepare_mixseq.py (benchmarks/kinker/prepare.py write_design): two rotations, half the
eligible lines (permutation seed 20261106 + 0) as reference, 70% fit / 30% holdout (default_rng(20261106)).
Counts: token IDs -> gene_metadata gene_symbol; tokens without a gene_metadata entry are dropped and counted;
counts of the kept tokens must be non-negative integers (else stop); duplicate symbols are summed.
lines.tsv flags: in_kinker, in_mixseq (any line present in those inputs), in_andor, in_schmid (MCF7, HCT116, A375,
COLO320), matched by DepMap ModelID. primary = chosen and none of the four flags.
Only per-cell UMI totals, cells per line and gene counts are reported (claim 8 blinding).
Usage: python -I prepare_tahoe.py EXTRACT_DIR OUT
"""
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'benchmarks/kinker'))
from prepare import DEPMAP, POSITIONS, depmap_rows  # noqa: E402

SEED = 20261106
POOL = 'plate1'
MIN_CELLS = 50
MAX_CELLS = 300
KEY = ['plate', 'sample', 'BARCODE_SUB_LIB_ID']
KINKER_CELLS = ROOT / 'results/kinker_input_v2_20261101_1551/cells_r0.tsv'
MIXSEQ_LINES = ROOT / 'results/heldout_mixseq_input_20261102_1159/lines.tsv'
ANDOR_LINES = ROOT / 'results/andor_20261101_1224/rna_input_v2/manifest.tsv'
SCHMID_LINES = ['MCF7', 'HCT116', 'A375', 'COLO320']


def norm(s):
    return re.sub('[^A-Z0-9]', '', str(s).upper())


def prior_ids(model):
    by_ccle = model.dropna(subset=['CCLEName']).set_index('CCLEName').ModelID
    by_name = model.assign(n=model.StrippedCellLineName.map(norm)).drop_duplicates('n').set_index('n').ModelID

    def ids(names, table, what):
        names = sorted(set(names) - {''})
        miss = [n for n in names if n not in table.index]
        print(f'{what}: {len(names)} lines, {len(miss)} without a DepMap ModelID {miss[:10]}', flush=True)
        return set(table.reindex(names).dropna())

    kinker = pd.read_csv(KINKER_CELLS, sep='\t', usecols=['line'], dtype=str).line.dropna()
    return dict(
        in_kinker=ids(kinker, by_ccle, 'Kinker'),
        in_mixseq=ids(pd.read_csv(MIXSEQ_LINES, sep='\t', dtype=str).line, by_ccle, 'MIX-seq'),
        in_andor=ids(pd.read_csv(ANDOR_LINES, sep='\t', dtype=str).line.map(norm), by_name, 'Andor'),
        in_schmid=ids([norm(x) for x in SCHMID_LINES], by_name, 'Schmid'))


def read_cells(src, keys):
    """Rows of plate1_dmso_tf.parquet whose cell_id is in keys, as (cell_id, tokens, values) arrays."""
    pf = pq.ParquetFile(src / 'plate1_dmso_tf.parquet')
    want = pa.array(sorted(keys))
    ids, toks, vals, lens = [], [], [], []
    for b in pf.iter_batches(batch_size=20000, columns=KEY + ['genes', 'expressions']):
        cid = pc.binary_join_element_wise(b['sample'], b['BARCODE_SUB_LIB_ID'], '_')
        m = pc.is_in(cid, value_set=want)
        if not pc.any(m).as_py():
            continue
        b = b.filter(m)
        ids.append(pc.binary_join_element_wise(b['sample'], b['BARCODE_SUB_LIB_ID'], '_').to_numpy(zero_copy_only=False))
        g, e = b['genes'], b['expressions']
        assert np.array_equal(pc.list_value_length(g).to_numpy(), pc.list_value_length(e).to_numpy())
        lens.append(pc.list_value_length(g).to_numpy())
        toks.append(pc.list_flatten(g).to_numpy())
        vals.append(pc.list_flatten(e).to_numpy())
    return np.concatenate(ids), np.concatenate(toks), np.concatenate(vals), np.concatenate(lens)


def main(src, out):
    src, out = Path(src), Path(out)
    out.mkdir(parents=True, exist_ok=False)
    rule = (src / 'selection_rule.txt').read_text()
    pass_value = re.search(r"pass_filter == '([^']+)'", rule).group(1)
    sel = pd.read_csv(src / 'selected_cells.tsv.gz', sep='\t', dtype=str)
    sel = sel[sel.pass_filter.eq(pass_value)].copy()
    sel['cell_id'] = sel['sample'] + '_' + sel.BARCODE_SUB_LIB_ID
    assert sel.cell_id.is_unique
    print(f'pass_filter == {pass_value!r}: {len(sel)} cells, samples {sorted(set(sel["sample"]))}', flush=True)

    clm = pd.read_parquet(src / 'metadata/cell_line_metadata.parquet')[['cell_name', 'Cell_ID_DepMap', 'Cell_ID_Cellosaur']]
    clm = clm.drop_duplicates()
    assert clm.cell_name.is_unique, 'cell_line_metadata: one cell_name with several IDs'
    sel = sel.merge(clm, on='cell_name', how='left', validate='m:1')
    assert sel.Cell_ID_Cellosaur.eq(sel.cell_line).all(), 'obs cell_line differs from cell_line_metadata Cellosaurus ID'
    model = pd.read_csv(DEPMAP / 'Model.csv', usecols=['ModelID', 'CCLEName', 'StrippedCellLineName'], dtype=str)
    ids = sorted(set(sel.Cell_ID_DepMap.dropna()) & set(model.ModelID))
    cn = depmap_rows(DEPMAP / 'OmicsCNGene.csv', ids)
    id2line = clm.dropna(subset=['Cell_ID_DepMap']).set_index('Cell_ID_DepMap').cell_name
    assert id2line.index.is_unique
    cn.index = cn.index.map(id2line)

    n = sel.groupby('cell_name').agg(cells=('cell_id', 'size'), depmap_id=('Cell_ID_DepMap', 'first')).reset_index()
    n = n.rename(columns={'cell_name': 'line'})
    n.insert(1, 'pool', POOL)
    n['has_cn'] = n.line.isin(cn.index)
    n['enough'] = n.cells >= MIN_CELLS
    n['chosen'] = n.has_cn & n.enough
    flags = prior_ids(model)
    for k, s in flags.items():
        n[k] = n.depmap_id.isin(s)
    n['primary'] = n.chosen & ~n[list(flags)].any(axis=1)

    rng = np.random.default_rng(SEED)
    take = []
    for line in sorted(n.loc[n.chosen, 'line']):
        c = np.sort(sel.loc[sel.cell_name.eq(line), 'cell_id'].to_numpy())
        take.extend(rng.choice(c, size=min(len(c), MAX_CELLS), replace=False) if len(c) > MAX_CELLS else c)
    take = np.sort(np.array(take, dtype=object))
    n['cells_used'] = n.line.map(sel[sel.cell_id.isin(set(take))].groupby('cell_name').size()).fillna(0).astype(int)
    n.to_csv(out / 'lines.tsv', sep='\t', index=False)
    print(f'lines {len(n)}; chosen {int(n.chosen.sum())}; primary {int(n.primary.sum())}; cells used {len(take)}', flush=True)
    print('primary lines:', ', '.join(sorted(n.loc[n.primary, 'line'])), flush=True)

    gm = pd.read_parquet(src / 'metadata/gene_metadata.parquet')
    assert gm.token_id.is_unique
    sym_all = np.unique(gm.gene_symbol.astype(str).to_numpy())
    tok2col = pd.Series(np.searchsorted(sym_all, gm.gene_symbol.astype(str).to_numpy()), index=gm.token_id.to_numpy())

    cid, tok, val, lens = read_cells(src, set(take))
    assert len(cid) == len(take) and set(cid) == set(take), f'{len(cid)} rows read for {len(take)} cells'
    row = np.repeat(np.arange(len(cid)), lens)
    known = np.isin(tok, tok2col.index.to_numpy())
    unknown = pd.Series(tok[~known]).value_counts()
    v = val[known]
    if not (np.all(v >= 0) and np.all(v == np.round(v))):
        raise ValueError('expressions of gene tokens are not non-negative integers')
    pairs = pd.DataFrame({'r': row[known], 't': tok[known]})
    assert not pairs.duplicated().any(), 'a token appears twice in one cell'
    X = sparse.csr_matrix((v.astype(np.int64), (row[known], tok2col.reindex(tok[known]).to_numpy())),
                          shape=(len(cid), len(sym_all)))
    X.sum_duplicates()
    order = np.argsort(cid)
    cid, X = cid[order], X[order].tocsr()
    assert np.array_equal(cid, take)
    X = X.astype(np.int32)
    qc = dict(tokens_without_gene_metadata={int(k): int(c) for k, c in unknown.items()},
              cells=len(cid), genes=len(sym_all), gene_tokens=len(gm), pass_filter=pass_value)

    meta = sel.set_index('cell_id').loc[cid]
    cells = pd.DataFrame(dict(cell_id=cid, pool=POOL, line=meta.cell_name.to_numpy(), margin=1.0))
    cells['umi'] = np.asarray(X.sum(axis=1)).ravel()
    cells['library'] = cells.umi.astype(float)
    tscp = meta.tscp_count.astype(int).to_numpy()
    qc['cells_umi_equal_obs_tscp_count'] = int((cells.umi.to_numpy() == tscp).sum())
    qc['median_umi_per_cell'] = float(cells.umi.median())
    meta.reset_index()[['cell_id', 'plate', 'sample', 'BARCODE_SUB_LIB_ID', 'cell_name', 'cell_line', 'Cell_ID_DepMap',
                        'pass_filter', 'tscp_count', 'gene_count', 'obs_row']].to_csv(out / 'cells_meta.tsv.gz', sep='\t', index=False)

    names = sorted(n.loc[n.chosen, 'line'])
    perm = np.random.default_rng(SEED + 0).permutation(len(names))
    half = {names[i] for i in perm[:len(names) // 2]}
    rng = np.random.default_rng(SEED)
    for r in (0, 1):
        c = cells.copy()
        c['role'], c['ref_group'], c['ref_split'], c['specimen'], c['design'] = 'other', '', '', '', 'unused'
        c['dataset'] = 'pool' + c.pool
        ref_lines = half if r == 0 else set(names) - half
        ref = c.line.isin(ref_lines)
        q = c.line.isin(set(names) - ref_lines)
        c.loc[ref, 'role'], c.loc[ref, 'ref_group'], c.loc[ref, 'design'] = 'reference', 'pool_ref', 'reference'
        c.loc[ref, 'ref_split'] = np.where(rng.random(ref.sum()) < 0.7, 'fit', 'holdout')
        c.loc[ref, 'specimen'] = f'pool{POOL}_ref'
        c.loc[q, 'role'], c.loc[q, 'design'] = 'query', 'single_line'
        c.loc[q, 'specimen'] = c.loc[q, 'line']
        c['patient'] = c.specimen
        c['truth_line'] = c.line
        c.to_csv(out / f'cells_r{r}.tsv', sep='\t', index=False)
        print(f'rotation {r}:', c.groupby('design').size().to_dict(), flush=True)
    cn.T.to_csv(out / 'depmap_cn.tsv.gz', sep='\t')
    pos = pd.read_csv(POSITIONS, sep='\t')
    pos[pos.gene.isin(sym_all)].to_csv(out / 'positions.tsv', sep='\t', index=False)
    qc['genes_with_position'] = int(pos.gene.isin(sym_all).sum())
    with h5py.File(out / 'counts.h5', 'w') as h:
        h.create_dataset('genes', data=sym_all.astype(h5py.string_dtype()))
        h.create_dataset('cell_ids', data=cells.cell_id.to_numpy().astype(h5py.string_dtype()))
        g = h.create_group('X')
        for k, a in [('data', X.data.astype(np.int32)), ('indices', X.indices.astype(np.int32)), ('indptr', X.indptr)]:
            g.create_dataset(k, data=a, compression='gzip')
    used = n.loc[n.chosen, 'cells_used']
    qc['cells_used_per_line_min_max'] = [int(used.min()), int(used.max())]
    (out / 'qc.json').write_text(json.dumps(qc, indent=1))
    print(json.dumps(qc), flush=True)


if __name__ == '__main__':
    main(*sys.argv[1:])
