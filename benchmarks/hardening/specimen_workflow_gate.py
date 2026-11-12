"""Empirical specimen gate calibrated on complete-workflow null scores."""
import argparse
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def statistic(folder, specimen='S1'):
    profiles = pd.read_csv(folder/'clones.tsv', sep='\t')
    assert pd.api.types.is_bool_dtype(profiles.normal_like)
    profiles = profiles[profiles.specimen.eq(specimen) & ~profiles.normal_like]
    assert profiles.clone.is_unique
    names = set(profiles.clone)
    if len(names) <= 1:
        return 0.
    differences = pd.read_csv(folder/'differences.tsv', sep='\t')
    selected = differences[differences.specimen.eq(specimen) & differences.clone_a.isin(names) & differences.clone_b.isin(names)]
    if not len(selected):
        return 0.
    if not np.isfinite(selected.z).all():
        raise ValueError('Nonfinite reported contrast score')
    selected = selected.assign(score=np.abs(selected.z.to_numpy()))
    return float(selected.groupby(['clone_a', 'clone_b']).score.sum().max())


def decision(value, banks):
    if not np.isfinite(value) or value < 0:
        raise ValueError('Finite nonnegative complete-workflow score required')
    ranks = []
    for bank in banks:
        bank = np.asarray(bank, float)
        if bank.ndim != 1 or len(bank) != 86 or not np.isfinite(bank).all() or np.any(bank < 0):
            raise ValueError('Frozen86-score null bank required')
        ranks.append(Fraction(1+int(np.sum(bank >= value)), len(bank)+1))
    if len(ranks) != 2:
        raise ValueError('Both matched and state null banks required')
    return dict(passed=bool(max(ranks) <= Fraction(1, 87)),
        inclusive_rank_p=[float(p) for p in ranks],
        interpretation='Conditional additional RNA distinction; not genomic validation')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'DESIGN.md')
    ranges = [(0, 22), (22, 44), (44, 65), (65, 86)]
    for condition in ('matched', 'state'):
        for start, stop in ranges:
            inputs[f'{condition}{start}'] = args.task/f'cal_{condition}_{start:02d}_v1/results.json'
    rec = start_run(args.out, inputs, dict(statistic='max final tumour-pair sum absolute regional z',
        null_bank_size=86, target=.05, simultaneous_calibration_tolerance_family=4),
        dict(scope='Complete native workflows; conditional on two fixed simulated null laws, not biological error control'),
        mode='specimen_workflow_gate_bank')
    try:
        # Derive the minimal n, exactly: four calibration families, each tail <= .05/4.
        assert Fraction(19, 20)**86 <= Fraction(1, 80)
        assert Fraction(19, 20)**85 > Fraction(1, 80)
        assert not decision(2., [np.full(86, 2.), np.full(86, 1.)])['passed']
        assert decision(3., [np.full(86, 2.), np.full(86, 1.)])['passed']
        assert not decision(0., [np.zeros(86), np.zeros(86)])['passed']
        results = []
        for name, path in inputs.items():
            if name not in ('script', 'design'):
                results.extend(json.loads(path.read_text()))
        assert len(results) == 344
        banks = {}
        for arm in ('production', 'candidate'):
            banks[arm] = {}
            for condition in ('matched', 'state'):
                rows = sorted([r for r in results if r['arm'] == arm and r['condition'] == condition], key=lambda r: r['replicate'])
                assert [r['replicate'] for r in rows] == list(range(86))
                assert all(r['seed'] == 202611094000+r['replicate'] for r in rows)
                scores = [r['raw_max_pair_score'] for r in rows]
                assert np.isfinite(scores).all() and np.all(np.asarray(scores) >= 0)
                banks[arm][condition] = scores
        write_json(banks, args.out/'banks.json')
        write_json(dict(condition_datasets=172, independent_base_seeds=86, complete_workflows=344,
            simultaneous_calibration_tolerance_family=4,
            probability_of_any_calibration_tail_exceeding_5pct_at_most=float(4*Fraction(19, 20)**86),
            thresholds={arm: {condition: max(scores) for condition, scores in per_condition.items()} for arm, per_condition in banks.items()},
            biological_specimens=0, empirical_validation='not yet run', production_changed=False), args.out/'report.json')
        finish_run(args.out, rec, backend='CPU exact tolerance calibration from complete-workflow null banks')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
