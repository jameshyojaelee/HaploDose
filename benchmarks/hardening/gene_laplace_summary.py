"""Fixed six-baseline conditional uncertainty for the frozen gene-Laplace stack."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    args=p.parse_args()
    args.out.mkdir(exist_ok=False)
    root=args.task.parents[1]
    names=('pseudobulk','pseudobulk_sm5','pseudobulk_mdl','predictive_smoothed','predictive_partition','stacked_smoothed')
    inputs=dict(script=Path(__file__),design=args.task/'REGRESSION_DESIGN.md',
        model=Path(__file__).with_name('gene_laplace_profile.py'),
        scorer=Path(__file__).with_name('statistics.py'),
        identities=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1/predictions_locked.json',
        groups=root/'dev/rna_predictive_smoothing_20261109/characterization_v1/RNA_choices_locked.json',
        old=root/'dev/rna_stacked_smoothing_20261110/summary_v3/scores.tsv')
    for start in (0,46,92,138):
        inputs[f'scores{start}']=args.task/f'regression_{start:03d}_v2/scores.tsv'
    rec=start_run(args.out,inputs,dict(family_size=6,n_boot=20000,seed=2026111003),
        dict(scope='Conditional exposed regression; eight shared groups, no calibrated coverage or independent evaluation'),mode='gene_laplace_summary')
    try:
        scores=pd.concat([pd.read_csv(inputs[f'scores{s}'],sep='\t',float_precision='round_trip')for s in (0,46,92,138)],ignore_index=True)
        assert len(scores)==181*8 and not scores.duplicated(['case','method']).any()
        old=pd.read_csv(inputs['old'],sep='\t',float_precision='round_trip')
        pd.testing.assert_frame_equal(old,scores[scores.method.isin(names)].reset_index(drop=True),check_exact=True)
        cases=[r['case']for r in json.loads(inputs['identities'].read_text())['cases']]
        matrix=scores.pivot(index='case',columns='method',values='pearson').loc[cases]
        assert np.isfinite(matrix[list(names)+['gene_laplace_stack']].to_numpy()).all()
        groups={r['case']:r['dependence_group']for r in json.loads(inputs['groups'].read_text())['widths']}
        labels=np.array([groups[c]for c in cases])
        assert len(np.unique(labels))==8
        spec=importlib.util.spec_from_file_location('gene_laplace_statistics',inputs['scorer'])
        stats=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(stats)
        records=[]
        for name in names:
            record=stats.paired_mean_comparison(matrix.gene_laplace_stack,matrix[name],clusters=labels,
                seed=2026111003,n_boot=20000,family_size=6)
            record.update(a='gene_laplace_stack',b=name,endpoint='pearson',independent_evaluation=False)
            records.append(record)
        records=stats.adjust_family(records)
        scores.to_csv(args.out/'scores.tsv',sep='\t',index=False)
        matrix.to_csv(args.out/'paired_profiles.tsv',sep='\t')
        write_json(dict(cases=181,previous_scores_exact=1086,means=matrix.mean().to_dict(),
            primary='gene_laplace_stack',standalone_is_ablation_only=True,family_size=6,
            contrasts=records,reference_groups=8,conditional_bootstrap=True,
            independent_evaluation=False,biological_calibration=False,default_promoted=False),args.out/'report.json')
        finish_run(args.out,rec,backend='CPU fixed six paired means, whole-reference-group bootstrap')
    except BaseException as error:
        finish_run(args.out,rec,error=error)
        raise


if __name__=='__main__':
    main()
