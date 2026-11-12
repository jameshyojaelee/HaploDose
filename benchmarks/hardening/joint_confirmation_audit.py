"""Independent preservation, stage-chain and numerical-result audit."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from haplodose.provenance import file_record


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=a.task.parents[1];report={}
    for name,path in [('baseline',a.task/'baseline'),('prior_final',root/'dev/mixture_confirmation_20261108/final_source')]:
        manifest=json.loads((path/'source_manifest.json').read_text())
        for relative,h in manifest.items():
            if digest(path/relative)!=h:raise ValueError('Frozen source changed')
        report[name+'_preserved']=len(manifest)
    baseline=json.loads((a.task/'baseline/source_manifest.json').read_text());unchanged=[]
    for relative,h in baseline.items():
        if relative.startswith(('python/','tests/','crates/')):
            if digest(root/relative)!=h:raise ValueError('Package/original test changed')
            unchanged.append(relative)
    report['unchanged_package_tests_native_files']=len(unchanged)
    # Follow the bounded historical inventory rather than inventing a complete
    # historical-data preservation claim.
    historical=json.loads((root/'dev/membership_20261107/preservation_audit.json').read_text())
    for name in ('selected_historical_results','unchanged_criteria_and_generators','frozen_panel_inputs','builds'):
        for record in historical[name]:
            if digest(record['path'])!=record['sha256']:raise ValueError('Historical artifact changed')
        report['preserved_'+name]=len(historical[name])
    cases=[]
    for stage in ('initial','matched'):
        summary=json.loads((a.task/f'panel_{stage}_v1/summary.json').read_text())
        for row in summary['cases']:
            seed=row['seed'];directory=a.task/f'panel_{stage}_v1/seed_{seed}'
            control=root/f'dev/confirmation_20261107/panel_merge_v1/seed_{seed}/fit'
            exact=np.array_equal(np.load(directory/'fit/posteriors.npy'),np.load(control/'posteriors.npy'),equal_nan=True)
            criteria=all(v['passed'] for v in row['criteria'].values())
            if criteria!=row['original_assertion_passed']:raise ValueError('Criteria and original assertion disagree')
            calls=[]
            for event in json.loads((directory/'stage_trace.json').read_text()):
                if event['stage']!='exact_confirmation':continue
                for call in event['calls']:
                    for segment in call['exact_segments']:
                        if segment['min_bin_posterior']<.95 or segment['n_test_bins']<5:
                            raise ValueError('Original scientific criterion violated')
                    sizes=call.get('per_stratum_per_group')
                    if sizes is not None and len(call['selected_cells'])!=2*sum(sizes):
                        raise ValueError('Depth-matching count conservation failed')
                    calls.append(dict(source=event['source'],accepted=event['accepted'],context=call['context'],
                                      original_segments=len(call['original_segments']),exact_segments=len(call['exact_segments']),
                                      log_evidence=call.get('log_evidence')))
            cases.append(dict(stage=stage,seed=seed,original_assertion_passed=criteria,
                              control_posteriors_identical=exact,calls=calls))
    report['cases']=cases
    bulk=a.task/'raw_operator_v2'
    if (bulk/'summary.json').exists():
        manifest=json.loads((bulk/'run_manifest.json').read_text())
        if manifest['status']!='complete':raise ValueError('Bulk analysis incomplete')
        for record in manifest['input']['files'].values():
            if digest(record['path'])!=record['sha256']:raise ValueError('Bulk input changed')
        scores=pd.read_csv(bulk/'patient_scores.tsv',sep='\t')
        if scores[['sample','method']].duplicated().any():raise ValueError('Duplicated scoring unit')
        if scores.groupby('method').size().nunique()!=1 or len(scores)!=48*7:raise ValueError('Comparator disappeared')
        comparisons=pd.read_csv(bulk/'comparisons.tsv',sep='\t')
        if len(comparisons)!=24 or not comparisons.family_size.eq(24).all():raise ValueError('Family changed')
        numeric=json.loads((bulk/'pearson_precision.json').read_text())
        report['bulk']=dict(methods=7,patients=48,contrasts=24,preserved_inputs=len(manifest['input']['files']),
            maximum_old_rebuilt_pearson_difference=max(abs(r['historical']-r['reconstructed']) for r in numeric),
            maximum_double_extended_difference=max(r['double_extended_difference'] for r in numeric))
        evaluation=scores[scores.split.eq('evaluation')]
        report['bulk']['method_summaries']=[]
        for method,group in evaluation.groupby('method',sort=False):
            report['bulk']['method_summaries'].append(dict(method=method,patients=len(group),
                means={c:float(group[c].mean()) for c in ('recall_loss','fpr_loss','precision_loss','rho_window','pearson_window','classified_truth_coverage')},
                evaluable={c:int(group[c].notna().sum()) for c in ('recall_loss','fpr_loss','precision_loss','rho_window','pearson_window')}))
        for row in comparisons.to_dict('records'):
            x=evaluation[evaluation.method.eq(row['a'])].set_index('sample')[row['endpoint']]
            y=evaluation[evaluation.method.eq(row['b'])].set_index('sample')[row['endpoint']].reindex(x.index)
            difference=(x-y).dropna().to_numpy(np.longdouble)
            expected=float(difference.sum()/len(difference))
            allowance=2*len(difference)*np.finfo(float).eps*max(1.,float(np.max(np.abs(difference))))
            if len(difference)!=row['n'] or abs(expected-row['mean_diff'])>allowance:
                raise ValueError('Paired effect/count reconstruction failed')
        report['bulk']['paired_effects_independently_recounted']=24
        # Supplemental inventory is explicitly after execution. Saved DNA,
        # masks and original scores already matched before operator comparison.
        paths=[root/'data/external/ensembl_grch37_r75/Homo_sapiens.GRCh37.75.gtf.gz',
               root/'python/haplodose/data/centromeres_hg38.tsv',
               Path('/gpfs/commons/home/jameslee/CRC/data/external/bian_cnv_annotations/cytoBand.txt.gz')]
        paths.extend(Path('/gpfs/commons/home/jameslee/CRC/data/external/tcga_pancan_cn')/f for f in
            ('TCGA_mastercalls.abs_tables_JSedit.fixed.txt','TCGA_mastercalls.abs_segtabs.fixed.txt'))
        report['bulk']['truth_source_inventory_after_execution']=[file_record(path) for path in paths]
    else:raise ValueError('Bulk results absent')
    (a.out/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True)


if __name__=='__main__':main()
