"""Independent exact calibration ranks and scalar correlation audit of bulk stages."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.provenance import start_run, finish_run, file_record


def threshold_reference(scores, neutral, columns):
    mass = defaultdict(Fraction); contributors = 0
    for j in columns:
        values = [float(v) for v in scores[:,j][neutral[:,j]] if np.isfinite(v)]
        if not values: continue
        contributors += 1
        for value,n in Counter(values).items(): mass[value] += Fraction(n,len(values))
    if not contributors: return np.nan
    cumulative = Fraction()
    for value in sorted(mass):
        cumulative += mass[value]
        if cumulative >= Fraction(contributors,20): return value
    raise AssertionError('Calibration probability mass incomplete')


def ranks(values):
    groups = Counter(map(float,values)); lookup = {}; position = 1
    for value,n in sorted(groups.items()):
        lookup[value] = (2*position+n-1)/2; position += n
    return np.array([lookup[float(value)] for value in values])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True); parser.add_argument('--task',type=Path,required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    stage = args.task/'stage_v1'; panel = args.task.parent/'hardening_20261106_2119/profile_v6'
    reference_path = Path(__file__).with_name('pearson_reference.py')
    spec = importlib.util.spec_from_file_location('independent_pearson',reference_path)
    reference = importlib.util.module_from_spec(spec); spec.loader.exec_module(reference)
    inputs = dict(script=Path(__file__),design=args.task/'STAGE_AUDIT_DESIGN.md',stage=stage/'run_manifest.json',
        reference=reference_path,support=panel/'raw_support_and_truth.npz',split=panel/'split_manifest.tsv')
    run = start_run(args.out,inputs,{},dict(scope='Independent numerical/scoring audit; no new inferential claim'),mode='bulk_posterior_stage_audit')
    try:
        reference.verify()
        manifest = json.loads((stage/'run_manifest.json').read_text()); assert manifest['status']=='complete'
        for record in manifest['output']['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
        with np.load(stage/'window_stages.npz') as archive: arrays = {key:archive[key] for key in archive.files}
        with np.load(panel/'raw_support_and_truth.npz',allow_pickle=True) as archive: frozen = {key:archive[key] for key in archive.files}
        for key,old in (('samples','samples'),('common','common'),('loss','loss'),('neutral','neutral')):
            assert np.array_equal(arrays[key],frozen[old])
        names = arrays['samples'].tolist(); common = arrays['common']
        split = pd.read_csv(panel/'split_manifest.tsv',sep='\t',dtype=str).set_index('sample').loc[names]
        ci = np.flatnonzero(split.split.eq('calibration')); assert len(ci)==24
        score = pd.read_csv(stage/'patient_stage_scores.tsv',sep='\t',float_precision='round_trip').set_index(['sample','method'])
        assert len(score)==288 and score.index.is_unique
        # Two exact fixtures: unequal patient window counts and tied quantiles.
        fixture = np.array([[0.,2.],[1.,2.],[1.,3.]])
        assert threshold_reference(fixture,np.ones_like(fixture,bool),[0,1]) == 0.
        assert threshold_reference(np.ones((2,2)),np.ones((2,2),bool),[0,1]) == 1.
        checks = []; pearson_checked = spearman_checked = 0
        for key,method in (('posterior','hd_posterior'),('segment','hd_dose'),('ratio','hd_lr')):
            values = arrays[key]; threshold = threshold_reference(values,arrays['neutral']&common,ci)
            for j,name in enumerate(names):
                row = score.loc[(name,method)]; mask = common[:,j]
                assert row.callable_windows==int(mask.sum()) and row.loss_threshold==threshold and row.threshold==threshold
                for field,truth in (('recall_loss',arrays['loss']),('fpr_loss',arrays['neutral'])):
                    eligible = truth[:,j]&mask
                    expected = sum(float(x)<threshold for x in values[eligible,j])/int(eligible.sum()) if eligible.any() else np.nan
                    assert (np.isnan(expected) and pd.isna(row[field])) or row[field]==expected
                x,y = values[mask,j],frozen['dna'][mask,j]
                valid = len(x)>=3 and np.ptp(x)>0 and np.ptp(y)>0
                if valid:
                    direct = reference.reference(x,y)
                    assert direct['float32_lower']<=row.pearson<=direct['float32_upper']
                    rank_reference = reference.reference(ranks(x),ranks(y))
                    assert abs(row.rho-rank_reference['value'])<=rank_reference['double_roundoff_budget']
                    pearson_checked += 1; spearman_checked += 1
                else: assert pd.isna(row.pearson) and pd.isna(row.rho)
            checks.append(dict(method=method,exact_calibration_threshold=threshold,patients=len(names)))
        summary = pd.read_csv(stage/'descriptive_summary.tsv',sep='\t',float_precision='round_trip').set_index(['split','method'])
        summaries_checked = 0
        plain = score.reset_index()
        for (role,method),reported in summary.iterrows():
            selected = plain[plain.split.eq(role)&plain.method.eq(method)]
            assert len(selected)==reported.n_samples==24
            for metric,column in (('pearson','mean_pearson'),('recall_loss','mean_loss_recall'),('fpr_loss','mean_loss_fpr')):
                values = selected[metric].dropna().tolist()
                expected = math.fsum(values)/len(values) if values else np.nan
                if not values: assert pd.isna(reported[column]); continue
                budget = len(values)*np.finfo(float).eps*(1+max(map(abs,values)))
                assert abs(reported[column]-expected)<=budget
                summaries_checked += 1
        thresholds = {row['method']:row['exact_calibration_threshold'] for row in checks}
        evaluation = np.flatnonzero(split.split.eq('evaluation')); assert len(evaluation)==24
        post_hit = arrays['posterior'][:,evaluation]<thresholds['hd_posterior']
        segment_hit = arrays['segment'][:,evaluation]<thresholds['hd_dose']
        raw_hit = arrays['ratio'][:,evaluation]<thresholds['hd_lr']
        contingencies = {}
        for label,mask in (('loss',arrays['loss']),('neutral',arrays['neutral'])):
            eligible = mask[:,evaluation]&common[:,evaluation]
            contingencies[label] = dict(both=int((eligible&post_hit&segment_hit).sum()),
                posterior_only=int((eligible&post_hit&~segment_hit).sum()),
                segment_only=int((eligible&~post_hit&segment_hit).sum()),
                neither=int((eligible&~post_hit&~segment_hit).sum()),eligible_windows=int(eligible.sum()))
        loss = arrays['loss'][:,evaluation]&common[:,evaluation]
        contingencies['raw_vs_posterior_loss'] = dict(both=int((loss&raw_hit&post_hit).sum()),
            raw_only=int((loss&raw_hit&~post_hit).sum()),posterior_only=int((loss&~raw_hit&post_hit).sum()),
            neither=int((loss&~raw_hit&~post_hit).sum()),eligible_windows=int(loss.sum()))
        (args.out/'report.json').write_text(json.dumps(dict(exact_patient_loss_checks=288,
            independent_Pearson_checks=pearson_checked,independent_Spearman_checks=spearman_checked,
            exact_calibration_thresholds=checks,descriptive_summary_checks=summaries_checked,
            frozen_callable_and_truth_masks_exact=True,comparison_window_call_contingencies=contingencies,
            window_counts_are_not_independent_biological_specimens=True,new_tests_or_method_selection=0,
            biological_calibration=False,scientific_threshold_changed=False),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error)
        raise


if __name__=='__main__': main()
