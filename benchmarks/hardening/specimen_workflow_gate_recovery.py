"""Apply a frozen complete-workflow gate without changing recovery assertions."""
import argparse
import importlib.util
import json
from pathlib import Path
import shutil

import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    native = args.task/'raw_candidate_v1'
    gate = load('frozen_workflow_gate', Path(__file__).with_name('specimen_workflow_gate.py'))
    scorer = load('unchanged_recovery_scoring', Path(__file__).with_name('recovery_baseline.py'))
    assertions = root/'dev/recovery_20261107/baseline/tests/test_split_mode.py'
    tests = load('original_recovery_assertions', assertions)
    banks_path = args.task/'bank_v1/banks.json'
    previous_path = root/'dev/depth_stratified_gene_noise_20261109/both_v1/report.json'
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md', gate=Path(gate.__file__),
        scorer=Path(scorer.__file__), assertions=assertions, banks=banks_path,
        raw_report=native/'report.json', previous_report=previous_path)
    for seed in scorer.SEEDS:
        for name in ('cells.tsv.gz', 'differences.tsv', 'clones.tsv', 'qc.json'):
            inputs[f'{seed}:{name}'] = native/f'panel/seed_{seed}/fit'/name
        inputs[f'{seed}:truth_view'] = native/f'panel/seed_{seed}/assertion_view/cells.tsv.gz'
    rec = start_run(args.out, inputs, dict(statistic='native final max pair sum abs(z)',
        rank_threshold='1/87 against each of two frozen null banks'),
        dict(scope='Additional RNA distinction support; failed gate means unresolved structure, not certified single genomic clone'),
        mode='specimen_workflow_gate_recovery')
    try:
        previous = json.loads(previous_path.read_text())['cases']
        raw = json.loads(inputs['raw_report'].read_text())['cases']
        assert len(raw) == len(previous) == 5
        for before, now in zip(previous, raw):
            assert before['seed'] == now['seed'] and before['criteria'] == now['criteria']
            assert before['original_assertion_passed'] == now['original_assertion_passed']
        banks = json.loads(banks_path.read_text())['candidate']
        results = []
        for record in raw:
            seed = record['seed']
            fit = native/f'panel/seed_{seed}/fit'
            value = gate.statistic(fit)
            decision = gate.decision(value, [banks['matched'], banks['state']])
            cells = pd.read_csv(inputs[f'{seed}:truth_view'], sep='\t')
            differences = pd.read_csv(fit/'differences.tsv', sep='\t')
            profiles = pd.read_csv(fit/'clones.tsv', sep='\t')
            names = set(profiles.loc[profiles.specimen.eq('S1') & ~profiles.normal_like, 'clone'])
            if len(names) > 1 and not decision['passed']:
                # Keep the tumour-call decision; abstain on its extra CNA structure.
                # This single unresolved RNA population is not a genomic clone claim.
                selected = cells.specimen.eq('S1') & cells.malignant_call
                cells.loc[selected, 'clone'] = 'S1_unresolved'
                cells.loc[selected, 'status'] = 'unresolved_workflow_evidence'
                supported_pair = differences.specimen.eq('S1') & differences.clone_a.isin(names) & differences.clone_b.isin(names)
                differences = differences.loc[~supported_pair]
            view = args.out/f'seed_{seed}'
            view.mkdir()
            cells.to_csv(view/'cells.tsv.gz', sep='\t', index=False)
            differences.to_csv(view/'differences.tsv', sep='\t', index=False)
            shutil.copyfile(fit/'qc.json', view/'qc.json')
            passed = True
            try:
                tests.test_thin_records_mode_and_recovers_clones({'thin': view})
            except AssertionError:
                passed = False
            results.append(dict(seed=seed, raw_max_pair_score=value, gate=decision,
                original_assertion_passed=passed, criteria=scorer.criteria(cells, differences)))
        retained = sum(r['original_assertion_passed'] for r in results)
        write_json(dict(raw_all30_criterion_values_reproduce=True, raw_recovered=4,
            gated_recovered=retained, cases=results, continuation_recovery_gate_passed=retained == 4,
            large_validation_authorized_by_design=retained == 4,
            biological_specimens=0, production_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU frozen-bank gate and unchanged original recovery assertions')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
