"""Audit fixed support, untouched comparators and paired aggregation/reference effects."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--previous',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    output=a.task/'factorial_v1';manifest=json.loads((output/'run_manifest.json').read_text())
    if manifest['status']!='complete':raise ValueError('Incomplete aggregation/reference ablation')
    for record in manifest['input']['files'].values():
        if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()!=record['sha256']:raise ValueError('Frozen input changed')
    old=pd.read_csv(a.previous/'patient_scores.tsv',sep='\t');new=pd.read_csv(output/'patient_scores.tsv',sep='\t')
    keys=['sample','method']
    if new[keys].duplicated().any() or len(new)!=48*13:raise ValueError('Missing/duplicated comparison units')
    retained=new[new.method.isin(old.method.unique())]
    pd.testing.assert_frame_equal(retained.sort_values(keys).reset_index(drop=True),old.sort_values(keys).reset_index(drop=True),
                                  check_dtype=False,rtol=0,atol=64*np.finfo(float).eps)
    comparisons=pd.read_csv(output/'comparisons.tsv',sep='\t')
    if len(comparisons)!=54 or not comparisons.family_size.eq(54).all():raise ValueError('Comparator family changed')
    evaluation=new[new.split.eq('evaluation')];effects=[]
    for row in comparisons.to_dict('records'):
        x=evaluation[evaluation.method.eq(row['a'])].set_index('sample')[row['endpoint']]
        y=evaluation[evaluation.method.eq(row['b'])].set_index('sample')[row['endpoint']].reindex(x.index)
        values=(x-y).dropna().to_numpy(np.longdouble);estimate=float(values.mean())
        allowance=2*len(values)*np.finfo(float).eps*max(1.,float(np.max(np.abs(values))))
        if len(values)!=row['n'] or abs(estimate-row['mean_diff'])>allowance:raise ValueError('Paired effect mismatch')
        effects.append(dict(endpoint=row['endpoint'],a=row['a'],b=row['b'],n=len(values),effect=estimate))
    methods=[]
    for name,group in evaluation.groupby('method',sort=False):
        fields=['recall_loss','fpr_loss','precision_loss','rho_window','pearson_window','classified_truth_coverage']
        methods.append(dict(method=name,means={c:float(group[c].mean()) for c in fields},
                            evaluable={c:int(group[c].notna().sum()) for c in fields}))
    base=json.loads((a.task/'baseline/source_manifest.json').read_text());unchanged=0
    for relative,digest in base.items():
        if hashlib.sha256((a.task/'baseline'/relative).read_bytes()).hexdigest()!=digest:raise ValueError('Baseline changed')
        if relative.startswith(('python/','tests/','crates/')):
            if hashlib.sha256((a.task.parents[1]/relative).read_bytes()).hexdigest()!=digest:raise ValueError('Production source changed')
            unchanged+=1
    report=dict(input_hashes_verified=len(manifest['input']['files']),baseline_files_preserved=len(base),
        production_test_native_files_unchanged=unchanged,original_comparators_unchanged=10,
        paired_effects=effects,method_summaries=methods,
        scope='Independent effect/count and preservation checks; bootstrap coverage is not validated')
    (a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)


if __name__=='__main__':main()
