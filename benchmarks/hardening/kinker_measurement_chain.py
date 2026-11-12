"""Bounded primary-chain audit of the two locked Kinker development cases."""
import argparse
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from haplodose.provenance import start_run, finish_run


CASES = ('2313287_STOMACH', '42MGBA_CENTRAL_NERVOUS_SYSTEM')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = Path(os.environ['HAPLODOSE_DATA_ROOT'])
    inp = root / 'results/kinker_input_v2_20261101_1551'
    scorer_path = Path(__file__).parents[1] / 'kinker/evaluate.py'
    spec = importlib.util.spec_from_file_location('kinker_chain_scorer', scorer_path)
    scorer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(scorer)
    inputs = dict(counts=inp/'counts.h5', truth=inp/'depmap_cn.tsv.gz',
                  positions=inp/'positions.tsv', scorer=scorer_path,
                  preparation=Path(__file__).parents[1]/'kinker/prepare.py',
                  design=args.task/'DESIGN.md')
    runs = [root/f'results/kinker_field_20261102_0529/haplodose_r{r}' for r in range(2)]
    for r, run in enumerate(runs):
        for name in ('cells.tsv.gz', 'genes.tsv.gz', 'bins.tsv', 'bin_counts.npy', 'profiles.tsv.gz', 'qc.json'):
            inputs[f'r{r}_{name}'] = run/name
    manifest = start_run(args.out, inputs, dict(cases=CASES),
        dict(scope='Locked retrospective development mechanism audit, no performance inference',
             counts='CPM-derived integer reconstruction; original molecules not established',
             DNA='Different assay/culture, relative reference measurement; not cell-level truth'),
        mode='kinker_measurement_chain')
    try:
        X, genes, cells = scorer.load_counts(inputs['counts'])
        scorer.unique_ids(cells, 'input cells'); scorer.unique_ids(genes, 'input genes')
        assert X.shape == (len(cells), len(genes))
        assert np.isfinite(X.data).all() and (X.data >= 0).all()
        assert np.array_equal(X.data, np.round(X.data))
        records, weights, support, identities = [], [], [], []
        for rotation, run in enumerate(runs):
            meta = pd.read_csv(run/'cells.tsv.gz', sep='\t', keep_default_na=False)
            used = pd.read_csv(run/'genes.tsv.gz', sep='\t')
            profiles = pd.read_csv(run/'profiles.tsv.gz', sep='\t')
            bins = pd.read_csv(run/'bins.tsv', sep='\t')
            scorer.unique_ids(meta.cell_id, 'run cells'); scorer.unique_ids(used.gene, 'used genes')
            ci = scorer.checked_indexer(cells, meta.cell_id, 'run cells')
            gi = scorer.checked_indexer(genes, used.gene, 'used genes')
            Y = np.load(run/'bin_counts.npy', mmap_mode='r')
            assert Y.shape == (len(meta), len(bins))
            # Independent sparse count aggregation, without package binning helpers.
            operator = sparse.csr_matrix((np.ones(len(used), dtype=np.int64),
                (np.arange(len(used)), used.bin.to_numpy())), shape=(len(used), len(bins)))
            for case in CASES:
                q = meta.specimen.eq(case) & meta.role.eq('query') & meta.design.eq('single_line')
                if not q.any():
                    continue
                query = meta.loc[q]
                assert query.truth_line.nunique() == query.dataset.nunique() == 1
                dataset = query.dataset.iloc[0]
                ref = meta.dataset.eq(dataset) & meta.role.eq('reference') & meta.ref_split.eq('fit')
                selected = np.flatnonzero(q | ref)
                raw = X[ci[selected]]
                full = np.asarray(raw.sum(axis=1)).ravel()
                reconstructed = np.asarray((raw[:, gi] @ operator).toarray())
                assert np.array_equal(reconstructed, Y[selected])
                assert np.array_equal(full, meta.iloc[selected].umi.to_numpy())
                identities.append(dict(case=case, rotation=rotation, input_rows=ci[selected].tolist(),
                    input_gene_columns=gi.tolist(), run_rows=selected.tolist(),
                    ordered_cells=meta.iloc[selected].cell_id.tolist(), ordered_genes=used.gene.tolist(),
                    gene_bins=used.bin.tolist(), reference_cells=meta.loc[ref, 'cell_id'].tolist(),
                    query_cells=query.cell_id.tolist()))
                local = meta.iloc[selected].copy()
                local['retained_umi'] = reconstructed.sum(axis=1)
                local['full_umi'] = full
                totals = local.loc[local.role.eq('reference')].groupby('line')[['retained_umi', 'full_umi']].sum()
                for line, row in totals.iterrows():
                    weights.append(dict(case=case, rotation=rotation, line=line,
                        full_umi=float(row.full_umi), retained_umi=float(row.retained_umi),
                        full_weight=float(row.full_umi/totals.full_umi.sum()),
                        retained_weight=float(row.retained_umi/totals.retained_umi.sum())))
                n = query.clone.value_counts()
                p = profiles.loc[profiles.specimen.eq(case)]
                assert not p[['clone', 'bin']].duplicated().any()
                missing = sorted(set(n.index)-set(p.clone))
                for col in ('log_ratio', 'log_dose', 'log_ratio_gene'):
                    if col not in p:
                        continue
                    pivot = p.pivot(index='bin', columns='clone', values=col).reindex(index=bins.bin)
                    # Missing clone or missing/nonfinite value never represents a measured zero.
                    values = pivot.reindex(columns=n.index).to_numpy()
                    covered = (np.isfinite(values)*n.to_numpy()[None]).sum(axis=1)
                    for j, covered_cells in enumerate(covered):
                        support.append(dict(case=case, rotation=rotation, quantity=col, bin=int(bins.bin.iloc[j]),
                            query_cells=len(query), represented_cells=int(covered_cells)))
                records.append(dict(case=case, rotation=rotation, dataset=dataset,
                    query_cells=len(query), reference_cells=int(ref.sum()), used_genes=len(used), bins=len(bins),
                    checked_count_entries=int(reconstructed.size),
                    statuses=query.status.value_counts().to_dict(), clone_counts=n.to_dict(),
                    missing_profile_clones=missing,
                    missing_profile_cells=int(n.reindex(missing).sum()),
                    max_reference_weight_difference=float(np.max(np.abs(
                        totals.full_umi/totals.full_umi.sum()-totals.retained_umi/totals.retained_umi.sum())))))
        assert sorted(r['case'] for r in records) == sorted(CASES)
        pd.DataFrame(weights).to_csv(args.out/'reference_weights.tsv', sep='\t', index=False)
        pd.DataFrame(support).to_csv(args.out/'profile_support.tsv', sep='\t', index=False)
        (args.out/'identities.json').write_text(json.dumps(identities, indent=2)+'\n')
        (args.out/'report.json').write_text(json.dumps(dict(cases=records,
            count_scale='Reconstructed from published CPM; integrality does not prove molecular counts',
            interpretation='No method change or DNA-selected reference reweighting'), indent=2)+'\n')
    except BaseException as error:
        finish_run(args.out, manifest, error=error)
        raise
    finish_run(args.out, manifest, backend='CPU NumPy/SciPy, no inference')


if __name__ == '__main__':
    main()
