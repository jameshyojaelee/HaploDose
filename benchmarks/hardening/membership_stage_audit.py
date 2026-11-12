"""Independent specimen decisions, bound inversion and frozen-source preservation."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import binom


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def records(p):return [json.loads(x) for x in p.read_text().splitlines()]

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
 p.add_argument('--previous',type=Path,required=True);p.add_argument('--panel-root',type=Path,required=True)
 a=p.parse_args();a.out.mkdir(exist_ok=False);report={}
 for label,root in [('initial',a.task/'baseline'),('previous_final',a.previous/'final_source')]:
    manifest=json.loads((root/'source_manifest.json').read_text())
    for file,h in manifest.items():
        if digest(root/file)!=h:raise ValueError('Frozen source changed: '+str(root/file))
    report[label+'_files_preserved']=len(manifest)
 fixed=('tests/test_split_mode.py','tests/test_pipeline.py','python/haplodose/simulate.py',
        'benchmarks/hardening/clones.py')
 for f in fixed:
    if digest(a.task/'baseline'/f)!=digest(a.task/'candidate_v1'/f):
        raise ValueError('Scientific assertion/generator changed')
 report['original_assertions_unchanged']=True
 historical=json.loads((a.previous/'preservation_audit.json').read_text())
 for label in ('selected_historical_results','unchanged_criteria_and_generators',
               'frozen_panel_inputs','builds'):
    for record in historical[label]:
        if digest(Path(record['path']))!=record['sha256']:
            raise ValueError('Previously frozen evidence changed: '+record['path'])
    report['preserved_'+label]=len(historical[label])
 prior_package=Path(json.loads((a.previous/'validation_22850847'/'installed_artifact.json').read_text())['package']['path']).parent
 for record in historical['current_package_matches_tested_source']:
    relative=Path(record['path']).relative_to(a.task.parents[1]/'python/haplodose')
    if digest(prior_package/relative)!=record['sha256']:
        raise ValueError('Prior installed package changed: '+str(relative))
 report['prior_installed_package_files_preserved']=len(historical['current_package_matches_tested_source'])
 old_manifest=json.loads((a.panel_root/'manifest.json').read_text());n=0
 for row in old_manifest['cases']:
    for r in row['input_files']+[row['model_metadata']]:
        if digest(Path(r['path']))!=r['sha256']:raise ValueError('Frozen panel input changed')
        n+=1
 report['panel_input_files_preserved']=n
 comparisons=[]
 for seed in (20261031,20261107,20261108,20261109,20261110):
    old=a.previous/'final_panel_v1'/f'seed_{seed}'/'thin'
    new=a.task/'panel_hard_v1'/f'seed_{seed}'/'fit'
    before,after=np.load(old/'posteriors.npy'),np.load(new/'posteriors.npy')
    if not np.array_equal(before,after,equal_nan=True):raise ValueError('Hard control changed probabilities')
    for name in ('cells.tsv.gz','differences.tsv','segments.tsv'):
        pd.testing.assert_frame_equal(pd.read_csv(old/name,sep='\t'),pd.read_csv(new/name,sep='\t'))
    comparisons.append(dict(seed=seed,posterior_exact=True,scientific_outputs_identical=True))
 report['hard_control']=comparisons
 summary=a.task/'null_summary_v1'/'summary.json'
 if summary.exists():
    saved=json.loads(summary.read_text());checks=[];evaluation=json.loads(summary.with_name('evaluation.json').read_text())
    shards=a.task/'null_shards_v1'
    for stage in ('hard','merge','both'):
        cal={name:[r['score'] for r in records(shards/stage/name/'records.jsonl')] for name in ('cal_nb','cal_state')}
        for row in (r for r in evaluation if r['stage']==stage):
            bank=None if row['method']=='raw' else cal['cal_nb' if row['method']=='NB_calibrated_gate' else 'cal_state']
            prob=None if bank is None else (1+sum(x>=row['score'] for x in bank))/200
            expected=row['raw_groups']>1 and (bank is None or (row['score']>0 and prob<=.05))
            if expected!=row['event'] or prob!=row['rank_p']:raise ValueError('Decision mismatch')
    for row in saved['summary']:
        subset=[r for r in evaluation if all(r[k]==row[k] for k in ('stage','condition','method'))]
        k,n=sum(r['event'] for r in subset),len(subset)
        if (k,n)!=(row['events'],row['n']):raise ValueError('Count mismatch')
        alpha=.05/27
        if row['direction']=='upper':
            value=1. if k==n else brentq(lambda x:binom.cdf(k,n,x)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
        else:value=0. if k==0 else brentq(lambda x:binom.sf(k-1,n,x)-alpha,0.,1.,xtol=np.finfo(float).eps,rtol=4*np.finfo(float).eps)
        difference=abs(value-row['simultaneous_bound'])
        # Independent beta inverse and binomial tail inversion; root solve is
        # requested to machine precision, tolerate accumulated tail arithmetic.
        if difference>64*np.finfo(float).eps:raise ValueError('Bound discrepancy')
        checks.append(dict(stage=row['stage'],condition=row['condition'],method=row['method'],events=k,n=n,
                           independent_bound=value,absolute_discrepancy=difference))
    report['independent_bounds']=checks
 else:report['null_audit']='not run: summary absent'
 (a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)

if __name__=='__main__':main()
