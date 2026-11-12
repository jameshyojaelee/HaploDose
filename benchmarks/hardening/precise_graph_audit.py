"""Independent recovery recount and processor-decision audit for graph precision."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True)
    p.add_argument('--panel-root',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);rows=[];hashes=0
    for allocation in ('default','pooled_exact'):
        root=a.task/('recovery_'+allocation+'_v1');summary=json.loads((root/'panel/summary.json').read_text())
        for directory in (root,root/'panel'):
            manifest=json.loads((directory/'run_manifest.json').read_text())
            if manifest['status']!='complete':raise AssertionError('Incomplete allocation arm')
            for record in manifest['input']['files'].values():
                if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()!=record['sha256']:
                    raise AssertionError('Source/input changed')
                hashes+=1
        if [r['seed'] for r in summary['cases']]!=[20261031,20261107,20261108,20261109,20261110]:
            raise AssertionError('Recovery population changed')
        for r in summary['cases']:
            directory=root/'panel'/f"seed_{r['seed']}"
            observed=pd.read_csv(directory/'fit/cells.tsv.gz',sep='\t')
            truth=pd.read_csv(a.panel_root/f"seed_{r['seed']}/data/cells.tsv",sep='\t')
            cells=observed.merge(truth[['cell_id','truth_clone']],on='cell_id',validate='one_to_one').to_dict('records')
            differences=pd.read_csv(directory/'fit/differences.tsv',sep='\t').to_dict('records')
            normal=[x for x in cells if x['specimen']=='S1' and x['truth_clone']=='normal']
            tumour=[x for x in cells if x['specimen']=='S1' and x['truth_clone']!='normal']
            selected=[x for x in tumour if x['truth_clone'] in ('clone1','clone2') and x['malignant_call']]
            names=sorted({x['truth_clone'] for x in selected});groups=sorted({x['clone'] for x in selected})
            table=Counter((x['truth_clone'],x['clone']) for x in selected)
            diagonal=lambda ns:sum(table[(name,group)] for name,group in zip(ns,groups))
            accuracy=max(diagonal(names),diagonal(names[::-1]))/len(selected) if selected else None
            chrom=sorted({x['chromosome'] for x in differences if x['specimen']=='S1' and
                          not x['clone_a'].endswith('_c0') and not x['clone_b'].endswith('_c0')})
            values=dict(normal_fpr=sum(x['malignant_call'] for x in normal)/len(normal),
                        tumour_tpr=sum(x['malignant_call'] for x in tumour)/len(tumour),s1_two_groups=len(groups),
                        s2_one_group=len({x['clone'] for x in cells if x['specimen']=='S2' and x['malignant_call']}),
                        assignment_accuracy=accuracy,difference_chromosomes=chrom)
            tests=dict(normal_fpr=values['normal_fpr']<.1,tumour_tpr=values['tumour_tpr']>.9,s1_two_groups=len(groups)==2,
                       s2_one_group=values['s2_one_group']==1,assignment_accuracy=accuracy is not None and accuracy>.65,
                       difference_chromosomes=chrom==['chr12'])
            for key,value in values.items():
                if value!=r['criteria'][key]['value'] or tests[key]!=r['criteria'][key]['passed']:
                    raise AssertionError('Independent criterion recount differs')
            if all(tests.values())!=r['original_assertion_passed']:raise AssertionError('Original assertion disagreement')
            rows.append(dict(allocation=allocation,seed=r['seed'],passed=all(tests.values()),values=values,
                             original_assertions_unchanged=True))
    processors=[]
    for replicate in (150,167):
        dirs=[a.task/f'{host}_precision_v1/replicate_{replicate}' for host in ('intel','amd')]
        summaries=[json.loads((d/'summary.json').read_text()) for d in dirs]
        if summaries[0]['input_sha256']!=summaries[1]['input_sha256']:raise AssertionError('Unpaired processor inputs')
        if not all(s['observer_identity'] for s in summaries):raise AssertionError('Instrumentation affected result')
        with np.load(dirs[0]/'result.npz') as x,np.load(dirs[1]/'result.npz') as y:
            np.testing.assert_array_equal(x['R'].argmax(axis=1),y['R'].argmax(axis=1))
            delta=float(np.max(np.abs(x['R']-y['R'])))
        traces=[s['trace'] for s in summaries]
        if len(traces[0])!=len(traces[1]):raise AssertionError('Processor search paths differ')
        graph_count=0;eligible_count=0
        for first,second in zip(*traces):
            if first['stage']!=second['stage']:raise AssertionError('Processor trace identity differs')
            for record in (first,second):
                artifact=record['arrays']
                if hashlib.sha256(Path(artifact['path']).read_bytes()).hexdigest()!=artifact['sha256']:
                    raise AssertionError('Trace arrays changed')
            if first['stage']=='leiden_groups':
                with np.load(first['arrays']['path']) as x,np.load(second['arrays']['path']) as y:
                    for name in ('idx','edges','groups','raw_membership'):
                        np.testing.assert_array_equal(x[name],y[name])
                    if np.any(x['idx']==np.arange(len(x['idx']))[:,None]):raise AssertionError('Self edge retained')
                graph_count+=1
            if first['stage']=='event_fit':
                if first.get('eligible')!=second.get('eligible'):raise AssertionError('Processor eligibility differs')
                for record in (first,second):
                    if record['returned']:
                        with np.load(record['arrays']['path']) as z:
                            sizes=[int(np.sum(z['q']<=.5)),int(np.sum(z['q']>.5))]
                        if sizes!=record['sizes']:raise AssertionError('Incorrect proposal support')
                        eligible=record['gain']>record['budget'] and min(sizes)>=30
                        if eligible!=record['eligible']:raise AssertionError('Incorrect eligibility decision')
                eligible_count+=int(bool(first.get('eligible')))
        processors.append(dict(replicate=replicate,identical_graphs=graph_count,identical_MAP_assignments=True,
            eligible_proposals=eligible_count,maximum_posterior_discrepancy=delta,
            interpretation='Posterior decimals differ; no relaxed tolerance or bitwise probability certificate'))
    base=json.loads((a.task/'baseline/source_manifest.json').read_text());unchanged=0
    for path,digest in base.items():
        if hashlib.sha256((a.task/'baseline'/path).read_bytes()).hexdigest()!=digest:raise AssertionError('Baseline changed')
        if path.startswith(('python/','tests/','crates/')):
            if hashlib.sha256((a.task.parents[1]/path).read_bytes()).hexdigest()!=digest:raise AssertionError('Core/test/native source changed')
            unchanged+=1
    (a.out/'report.json').write_text(json.dumps(dict(cases=rows,input_hashes_verified=hashes,
        processor_comparisons=processors,baseline_files_preserved=len(base),
        core_test_native_unchanged=unchanged,scope='Numerical graph correction; no independent clone/error or power guarantee'),indent=2))
    print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
