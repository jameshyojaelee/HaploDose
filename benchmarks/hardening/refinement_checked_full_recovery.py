"""Run all five original recovery cases only after the causal rescue passes."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
from pathlib import Path
import sys

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def load(name):
    spec=importlib.util.spec_from_file_location('refinement_full_'+name,Path(__file__).with_name(name+'.py'))
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


class IdentifiedRecords(list):
    def append(self, record):
        caller=sys._getframe(1)
        assert caller.f_code.co_name=='checked'
        model=caller.f_locals['self']
        super().append(dict(record,source_binned_counts=str(model.rows.Y.filename)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    candidate=args.task/'candidate_v1'
    assert json.loads((candidate/'run_manifest.json').read_text())['status']=='complete'
    result=json.loads((candidate/'report.json').read_text())
    rows=[r for r in result['results']if r['arm']=='pooled_prevalence']
    assert len(rows)==1 and rows[0]['reported_groups']==2 and rows[0]['exact_chr12_recovery']
    old=args.task.parent/'phase8_membership_trace_20261110/trace_v1/replay/results.json'
    previous=[r for r in json.loads(old.read_text())if r['arm']=='pooled_prevalence']
    assert len(previous)==1 and rows[0]['assignment_correct']>previous[0]['assignment_correct']
    assert result['production_parity']['exact']
    assert any(not r['passed']for r in result['preview_records'])
    method=load('refinement_checked_split')
    for name in ('sham_v1','candidate_v1'):
        manifest=json.loads((args.task/name/'run_manifest.json').read_text())
        assert manifest['input']['files']['method']['sha256']==file_record(method.__file__)['sha256']
    assert json.loads((args.task/'sham_v1/report.json').read_text())['sham_parity']['exact']
    args.out.mkdir(exist_ok=False)
    source=load('pooled_prevalence_workflow');original=source.safeguard
    inputs=dict(script=Path(__file__),method=Path(method.__file__),design=args.task/'FULL_DESIGN.md',
                sham=args.task/'sham_v1/report.json',candidate=candidate/'report.json',
                original_results=old,original_full=args.task.parent/'pooled_prevalence_workflow_20261110/candidate_v2/report.json')
    run=start_run(args.out,inputs,dict(cases=5,original_criteria=30),
                  dict(scope='All original exposed recovery criteria; not independent validation or error control'),
                  mode='refinement_checked_full_recovery')
    records=IdentifiedRecords()
    @contextmanager
    def intervene(*a,**kw):
        with original(*a,**kw):
            with method.checked_splits(records,apply=True):yield
    source.safeguard=intervene;argv=sys.argv
    try:
        # The historical runner's own prerequisites/configuration remain tied
        # to its original task. Only the output and declared wrapper change.
        historical=args.task.parent/'pooled_prevalence_workflow_20261110'
        sys.argv=[str(source.__file__),'--task',str(historical),'--out',str(args.out/'workflow')]
        source.main()
        report=json.loads((args.out/'workflow/report.json').read_text())
        assert len(report['cases'])==5
        passed=bool(report['recovered']==5 and all(r['original_assertion_passed']
            and all(v['passed']for v in r['criteria'].values())for r in report['cases']))
        write_json(dict(recovered=report['recovered'],cases=report['cases'],
            all_original_criteria_passed=passed,preview_records=records,
            original_scientific_criteria_changed=False,false_split_control_established=False,
            independent_evaluation=False,biological_specimens=0,promotion=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU full original recovery with unchanged refinement checks')
    except BaseException as error:
        write_json(records,args.out/'preview_records_failed.json');finish_run(args.out,run,error=error);raise
    finally:
        source.safeguard=original;sys.argv=argv


if __name__=='__main__':main()
