"""Audit fixed-fold stage ablations without changing original scientific criteria."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--task', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(exist_ok=False)
    path = Path(__file__).with_name('offset_objective_audit.py')
    spec = importlib.util.spec_from_file_location('independent_original_criteria', path)
    reference = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(reference)
    old = a.task.parent/'gene_thinning_20261109/gene_v1'
    inputs = dict(script=Path(__file__), scorer=path, design=a.task/'DESIGN.md',
        reference_summary=old/'panel/summary.json', reference_trace=old/'trace.json')
    for arm, folder in [('reference', old)]+[(s, a.task/f'{s}_v1') for s in ('confirmation', 'merge', 'both')]:
        for entry in json.loads((folder/'trace.json').read_text()):
            inputs[f'{arm}:fold:{entry["model"]}'] = Path(entry['path'])
        for path in sorted(folder.rglob('*')):
            if path.is_file() and path.name in ('summary.json', 'stage_trace.json', 'qc.json', 'cells.tsv.gz',
                    'differences.tsv', 'run_manifest.json', 'trace.json'):
                inputs[f'{arm}:{path.relative_to(folder)}'] = path
    for seed in (20261031, 20261107, 20261108, 20261109, 20261110):
        inputs[f'{seed}:truth'] = a.task.parent/f'recovery_20261107/recovery_panel_v1/seed_{seed}/data/cells.tsv'
    rec = start_run(a.out, inputs, {}, dict(scope='Fixed development ablation and original criteria, no calibration'),
        mode='gene_thinning_weighted_audit')
    try:
        before = json.loads((old/'trace.json').read_text())
        results, outcomes = [], []
        for stage in ('confirmation', 'merge', 'both'):
            folder = a.task/f'{stage}_v1'
            summary, decisions = reference.recovery(a.task, stage, 'v1')
            assert summary['stage'] == stage
            outcomes.extend(decisions)
            current = json.loads((folder/'trace.json').read_text())
            assert len(before) == len(current) == 15
            for previous, candidate in zip(before, current):
                for name in ('model', 'case_seed', 'cells', 'genes', 'bins', 'seed', 'share', 'oracle_gene_size', 'veto_calls'):
                    assert previous[name] == candidate[name]
                with np.load(previous['path']) as b, np.load(candidate['path']) as c:
                    assert set(b.files) == set(c.files)
                    for name in b.files:
                        assert np.array_equal(b[name], c[name]), (stage, candidate['model'], name)
            traces = []
            for row in summary['cases']:
                entries = json.loads((folder/f'panel/seed_{row["seed"]}/stage_trace.json').read_text())
                if stage in ('confirmation', 'both'):
                    assert any(e['stage'] == 'confirmation' for e in entries)
                traces.append(dict(seed=row['seed'], confirmation_calls=sum(e['stage'] == 'confirmation' for e in entries),
                    merge_pairs=sum(e['stage'] == 'merge' for e in entries)))
            results.append(dict(stage=stage, passed=sum(d['passed'] for d in decisions), total=5,
                all15_folds_identical=True, trace=traces))
        write_json(dict(all_checks_passed=True, independent_criteria_checked=90, outcomes=outcomes,
            arms=results, folds_compared=45, biological_specimens=0, production_changed=False,
            scientific_criteria_changed=False, calibration_status='Not run; original whole-specimen5%/power gate unchanged'),
            a.out/'report.json')
    except BaseException as error:
        finish_run(a.out, rec, error=error)
        raise
    finish_run(a.out, rec, backend='CPU independent90 criteria and45 exact fold comparisons')


if __name__ == '__main__':
    main()
