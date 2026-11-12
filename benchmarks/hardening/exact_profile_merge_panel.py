"""Matched exposed-panel regression of duplicate repair and research weight units."""
import argparse
import importlib.util
import json
from pathlib import Path

from haplodose import clones, pipeline
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    p.add_argument('--condition', choices=['matched', 'state', 'two_clone'], required=True)
    p.add_argument('--start', type=int, required=True)
    p.add_argument('--stop', type=int, required=True)
    a = p.parse_args()
    assert 0 <= a.start < a.stop <= 100
    a.out.mkdir(exist_ok=False)
    root = Path(__file__).parent
    gene = load('panel_gene', root/'gene_thinning_hook.py')
    stage = load('panel_stage', root/'membership_weighted_stages.py')
    veto = load('panel_veto', root/'gene_veto_weights.py')
    scorer = load('panel_score', root/'gene_veto_fresh.py')
    old = a.task.parent/'gene_veto_fresh_20261109'
    rows = []
    for path in sorted(old.glob('chunk_*_v1/results.json')):
        rows.extend(json.loads(path.read_text()))
    cases = sorted((r for r in rows if r['condition']==a.condition and r['arm']=='candidate'
                    and a.start <= r['rep'] < a.stop), key=lambda r: r['rep'])
    assert [r['rep'] for r in cases] == list(range(a.start, a.stop))
    inputs = dict(script=Path(__file__), design=a.task/'PANEL_DESIGN.md', prior_audit=old/'evaluation_audit_v1/report.json')
    for row in cases:
        for name in ('counts.h5', 'model_cells.tsv', 'positions.tsv', 'truth_cells.tsv'):
            inputs[f'{row["rep"]}:{name}'] = Path(row['case'])/'input'/name
    rec = start_run(a.out, inputs, dict(condition=a.condition, start=a.start, stop=a.stop),
                    dict(scope='Previously exposed simulations: regression only; biological n=0'), mode='exact_merge_panel')
    result = []
    try:
        assert json.loads(inputs['prior_audit'].read_text())['all_checks_passed']
        import pandas as pd
        import numpy as np
        for row in cases:
            case = Path(row['case'])
            for units in ('raw', 'effective'):
                out = a.out/f'seed_{row["seed"]}'/units
                out.mkdir(parents=True)
                for name in ('folds', 'veto'):
                    (out/name).mkdir()
                gt, st, vt = [], [], []
                with gene.gene_native(gt, out/'folds', case_seeds=(row['seed'],)), \
                        stage.weighted_stages('both', st), \
                        veto.observe_weights(vt, out/'veto', apply=True, weight_units=units):
                    pipeline.run(case/'input/counts.h5', case/'input/model_cells.tsv', case/'input/positions.tsv', out/'fit',
                        params=clones.CloneParams(split_mode='thin'), chunk=97, counts_kind='raw', genome_build='hg38',
                        coordinate_system='one-based-closed', gene_namespace='symbol')
                # Verify the same gene folds, including historical raw candidate.
                for fold in (out/'folds').glob('*.npz'):
                    previous = np.load(case/'candidate/folds'/fold.name)
                    current = np.load(fold)
                    assert set(previous.files)==set(current.files)
                    assert all(np.array_equal(previous[k], current[k]) for k in previous.files)
                write_json(dict(gene=gt, stage=st, veto=vt), out/'trace.json')
                scores = scorer.score(out/'fit', pd.read_csv(case/'input/truth_cells.tsv', sep='\t'), a.condition)
                scores.update(rep=row['rep'], seed=row['seed'], units=units, out=str(out), prior=row)
                result.append(scores)
                write_json(result, a.out/'results.json')
        finish_run(a.out, rec, backend='CPU isolated corrected wheel; exact historical input/fold pairing')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
