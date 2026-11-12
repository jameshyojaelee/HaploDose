"""Audit old/new allele counterexamples, package checks and preserved source."""
import argparse
import hashlib
import json
from pathlib import Path
import re

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    validation=args.task/'validation_23002524'
    inputs=dict(script=Path(__file__),baseline=args.task/'baseline/source_manifest.json',
                candidate=args.task/'candidate_v2/source_manifest.json',artifact=validation/'installed_artifact.json',
                pytest=validation/'pytest.txt',rust=validation/'cargo_test.txt')
    for label in ('old','corrected'):
        inputs[label]=args.task/(label+'_v1')/'report.json'
        inputs[label+'_manifest']=inputs[label].parent/'run_manifest.json'
        manifest=json.loads(inputs[label+'_manifest'].read_text())
        assert manifest['status']=='complete'
        for group in ('input','output'):
            for record in manifest[group]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
    run=start_run(args.out,inputs,{},dict(scope='Computational correction; no biological validation'),
                  mode='allele_cell_model_audit')
    try:
        old=json.loads(inputs['old'].read_text());new=json.loads(inputs['corrected'].read_text())
        assert not old['model_consistent'] and new['model_consistent']
        assert old['checks'][0]==new['checks'][0], 'Legacy None counterexample changed'
        assert not old['checks'][1]['proposal_matches'] and not old['checks'][1]['balance_matches']
        assert new['checks'][1]['proposal_matches'] and new['checks'][1]['balance_matches']
        for name in ('expected_proposal','independent_log_odds','expected_balanced','independent_balance_probability'):
            assert old['checks'][1][name]==new['checks'][1][name]
        baseline=json.loads(inputs['baseline'].read_text());candidate=json.loads(inputs['candidate'].read_text())
        for directory,inventory in ((args.task/'baseline',baseline),(args.task/'candidate_v2',candidate)):
            for name,digest in inventory.items():
                assert hashlib.sha256((directory/name).read_bytes()).hexdigest()==digest
        changed=sorted(name for name in set(baseline)|set(candidate)
                       if name.startswith(('python/','crates/','tests/')) and baseline.get(name)!=candidate.get(name))
        assert changed==['python/haplodose/alleles.py','python/haplodose/clones.py','tests/test_allele_cell_model.py']
        assert '348 passed, 1 skipped' in inputs['pytest'].read_text()
        assert re.search(r'test result: ok\. 4 passed;',inputs['rust'].read_text())
        artifact=json.loads(inputs['artifact'].read_text())
        for key in ('wheel','native'):
            assert file_record(artifact[key]['path'])['sha256']==artifact[key]['sha256']
        (args.out/'report.json').write_text(json.dumps(dict(legacy_probe_exact=True,
            explicit_beta_binomial_counterexamples_corrected=2,python_passed=348,rust_passed=4,gpu_only_skipped=1,
            changed_core_files=changed,baseline_source_files_preserved=len(baseline),
            scientific_thresholds_or_defaults_changed=False,biological_validation=False,
            wheel_sha256=artifact['wheel']['sha256']),indent=2))
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
