"""Post-freeze shared-reference uncertainty for the new Kinker contrast."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=Path('/gpfs/commons/home/jameslee/Haplodose')
    identities=root/'dev/reference_weighted_pseudobulk_20261109/regression_v1/predictions_locked.json'
    previous=root/'dev/rna_predictive_smoothing_20261109/characterization_v1/RNA_choices_locked.json'
    inputs=dict(script=Path(__file__),design=args.task/'UNCERTAINTY_DESIGN.md',identities=identities,
        groups=previous,scores=args.task/'kinker_v1/scores.tsv',
        prior_scores=root/'dev/rna_predictive_smoothing_20261109/regression_v1/scores.tsv')
    run=start_run(args.out,inputs,dict(draws=20000,seed=2026110912,contrast_family=2),
        dict(scope='Post hoc conditional regression characterization, not prospective or calibrated coverage'),
        mode='predictive_partition_uncertainty')
    try:
        cases=[r['case'] for r in json.loads(identities.read_text())['cases']]
        groups={r['case']:r['dependence_group'] for r in json.loads(previous.read_text())['widths']}
        assert len(cases)==len(set(cases))==181 and set(cases)==set(groups)
        labels=np.array([groups[c] for c in cases]);clusters=np.unique(labels);assert len(clusters)==8
        matrix=pd.read_csv(inputs['scores'],sep='\t',float_precision='round_trip').pivot(index='case',columns='method',values='pearson').loc[cases]
        old=pd.read_csv(inputs['prior_scores'],sep='\t',float_precision='round_trip').pivot(index='case',columns='method',values='pearson').loc[cases]
        for name in old.columns:np.testing.assert_array_equal(matrix[name],old[name])
        rng=np.random.default_rng(2026110912)
        draws=rng.integers(len(clusters),size=(20000,len(clusters)))
        size=np.array([sum(labels==c) for c in clusters]);results=[];stored={}
        for baseline in ('predictive_smoothed','pseudobulk_sm5'):
            delta=(matrix.predictive_partition-matrix[baseline]).to_numpy()
            totals=np.array([math.fsum(delta[labels==c]) for c in clusters])
            bootstrap=totals[draws].sum(axis=1)/size[draws].sum(axis=1)
            # Independent scalar arithmetic checks the ratio estimand for the
            # first100 sampled cluster sets without changing its definition.
            for row,value in zip(draws[:100],bootstrap[:100]):
                reference=math.fsum(float(totals[j]) for j in row)/sum(int(size[j]) for j in row)
                allowance=64*np.finfo(float).eps*(1+abs(reference))
                assert abs(value-reference)<=allowance
            ordered=np.sort(bootstrap);endpoints=(.0125,.9875)
            interval=np.quantile(ordered,endpoints).tolist();brackets=[]
            for p in endpoints:
                lo=max(0,int(binom.ppf(.00625,len(ordered),p))-1)
                hi=min(len(ordered)-1,int(binom.ppf(.99375,len(ordered),p)))
                brackets.append([float(ordered[lo]),float(ordered[hi])])
            stored[baseline]=bootstrap
            results.append(dict(baseline=baseline,mean_paired_difference=float(delta.mean()),
                percentile_97p5_interval=interval,Monte_Carlo_quantile_brackets=brackets,
                cluster_effects=[dict(group=int(c),cases=int(size[j]),mean=float(totals[j]/size[j])) for j,c in enumerate(clusters)],
                working_interval_above_zero=bool(interval[0]>0),
                interpretation='Approximate conditional bootstrap; two-contrast Bonferroni family, not a biological coverage guarantee'))
        np.savez_compressed(args.out/'bootstrap.npz',**stored)
        write_json(dict(cases=181,independent_clusters_assumed=8,group_sizes=size.tolist(),contrasts=results,
            family_correction='Bonferroni two contrasts; 97.5% marginal percentile intervals',
            MC_family='four endpoint brackets; two-sided alpha.0125 each',
            independent_evaluation=False,biological_calibration=False,production_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU fixed-score shared-reference cluster resampling')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
