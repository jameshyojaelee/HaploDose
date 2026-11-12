"""Independent original-criterion recount for all allocation-screen outputs."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--control',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);rows=[];hashes=0;scale_checks=0
    for allocation,fraction in [('half',.5),('one_third',1/3),('two_thirds',2/3)]:
        root=a.task/(allocation+'_v1');summary=json.loads((root/'summary.json').read_text())
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
            cells=pd.read_csv(directory/'assertion_view/cells.tsv.gz',sep='\t').to_dict('records')
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
            trace=json.loads((directory/'stage_trace.json').read_text())
            for e in trace:
                for call in e.get('calls',[]):
                    if call['confirmation_to_proposal_scale']!=(1-fraction)/fraction:
                        raise AssertionError('Confirmation mean/shape share not transported')
                    scale_checks+=1
            if allocation=='half':
                old=a.control/f"seed_{r['seed']}"/'fit'
                np.testing.assert_array_equal(np.load(directory/'fit/posteriors.npy'),np.load(old/'posteriors.npy'))
                for name in ('cells.tsv.gz','segments.tsv','differences.tsv','clones.tsv'):
                    pd.testing.assert_frame_equal(pd.read_csv(directory/'fit'/name,sep='\t'),
                                                  pd.read_csv(old/name,sep='\t'),check_dtype=False,check_exact=True)
            rows.append(dict(allocation=allocation,seed=r['seed'],passed=all(tests.values()),values=values,
                             original_assertions_unchanged=True))
    base=json.loads((a.task/'baseline/source_manifest.json').read_text());unchanged=0
    for path,digest in base.items():
        if hashlib.sha256((a.task/'baseline'/path).read_bytes()).hexdigest()!=digest:raise AssertionError('Baseline changed')
        if path.startswith(('python/','tests/','crates/')):
            if hashlib.sha256((a.task.parents[1]/path).read_bytes()).hexdigest()!=digest:raise AssertionError('Core/test/native source changed')
            unchanged+=1
    (a.out/'report.json').write_text(json.dumps(dict(cases=rows,input_hashes_verified=hashes,
        confirmation_scale_checks=scale_checks,half_control_exact_cases=5,baseline_files_preserved=len(base),
        core_test_native_unchanged=unchanged,scope='Allocation screen only; no independent clone/error or power guarantee'),indent=2))
    print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
