"""Fixed five-comparator shared-reference uncertainty for RNA stacking."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args = p.parse_args();args.out.mkdir(exist_ok=False)
    root = args.task.parents[1]
    old_names = ('pseudobulk','pseudobulk_sm5','pseudobulk_mdl','predictive_smoothed','predictive_partition')
    inputs = dict(script=Path(__file__),design=args.task/'REGRESSION_DESIGN.md',
        identities=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1/predictions_locked.json',
        groups=root/'dev/rna_predictive_smoothing_20261109/characterization_v1/RNA_choices_locked.json',
        old=root/'dev/predictive_partition_20261109/kinker_v1/scores.tsv',
        statistics=Path(__file__).with_name('statistics.py'))
    for start in (0,46,92,138):inputs[f'scores{start}'] = args.task/f'regression_{start:03d}_v3/scores.tsv'
    rec = start_run(args.out,inputs,dict(family_size=5,n_boot=20000,seed=2026111001),
        dict(scope='Exposed regression and fixed reference groups, not prospective or calibrated coverage'),mode='rna_stacked_summary')
    try:
        scores = pd.concat([pd.read_csv(inputs[f'scores{s}'],sep='\t',float_precision='round_trip') for s in (0,46,92,138)],ignore_index=True)
        assert len(scores) == 181*6 and not scores.duplicated(['case','method']).any()
        old = pd.read_csv(inputs['old'],sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(old.reset_index(drop=True),scores[scores.method.isin(old_names)].reset_index(drop=True),check_exact=True)
        cases = [r['case'] for r in json.loads(inputs['identities'].read_text())['cases']]
        matrix = scores.pivot(index='case',columns='method',values='pearson').loc[cases]
        assert np.isfinite(matrix.to_numpy()).all() and matrix.shape == (181,6)
        groups = {r['case']:r['dependence_group'] for r in json.loads(inputs['groups'].read_text())['widths']}
        labels = np.array([groups[c] for c in cases]);assert len(np.unique(labels)) == 8
        spec = importlib.util.spec_from_file_location('stacked_statistics',inputs['statistics'])
        stats = importlib.util.module_from_spec(spec);spec.loader.exec_module(stats)
        records = []
        for name in old_names:
            r = stats.paired_mean_comparison(matrix.stacked_smoothed,matrix[name],clusters=labels,
                seed=2026111001,n_boot=20000,family_size=5)
            r.update(a='stacked_smoothed',b=name,endpoint='pearson',independent_evaluation=False)
            records.append(r)
        records = stats.adjust_family(records)
        scores.to_csv(args.out/'scores.tsv',sep='\t',index=False)
        matrix.to_csv(args.out/'paired_profiles.tsv',sep='\t')
        write_json(dict(cases=181,old_scores_exact=905,methods=matrix.mean().to_dict(),
            comparator_family=5,contrasts=records,reference_groups=8,conditional_reference_bootstrap=True,
            independent_evaluation=False,biological_calibration=False,production_changed=False,
            default_promoted=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU fixed family of five paired means, whole-reference-group bootstrap')
    except BaseException as error:
        finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
