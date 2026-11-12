"""Independent binary-partition identity audit; no scientific decisions change."""
import argparse
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.provenance import start_run, finish_run, file_record


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    args=p.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parents[1];prior=root/'dev/rna_population_baseline_20261109'
    model_path=Path(__file__).with_name('rna_population_baseline.py')
    generator_path=root/'dev/recovery_next_20261107/rank_source_v1/benchmarks/hardening/global_rank_generator.py'
    model=load('corrected_label_baseline',model_path);generator=load('frozen_label_generator',generator_path)
    inputs=dict(script=Path(__file__),model=model_path,generator=generator_path,design=args.task/'LABEL_DESIGN.md',
        probe=args.task/'probe_v2/report.json',previous=prior/'summary_v2/report.json')
    rows=[]
    for part in ('calibration','evaluation'):
        directory=prior/(part+'_v2');m=json.loads((directory/'run_manifest.json').read_text());assert m['status']=='complete'
        record=m['output']['files']['records.jsonl'];assert file_record(record['path'])['sha256']==record['sha256']
        inputs[part+':manifest']=directory/'run_manifest.json';inputs[part+':records']=Path(record['path'])
        rows.extend(json.loads(line) for line in Path(record['path']).read_text().splitlines())
    run=start_run(args.out,inputs,dict(cases=998,refits='Exactly8 independently selected numerical probe cases'),
        dict(scope='Label-invariant diagnostic correction; no threshold, method or scientific criterion change'),mode='rna_population_label_audit')
    try:
        fixture_checks=0
        for a,b in itertools.product(itertools.product((0,1),repeat=5),repeat=2):
            literal=sum(x!=y for x,y in zip(a,b));aligned=min(literal,5-literal)
            exhaustive=min(sum(x!=permutation[y] for x,y in zip(a,b)) for permutation in ((0,1),(1,0)))
            pair_disagreements=sum((a[i]==a[j])!=(b[i]==b[j]) for i in range(5) for j in range(i+1,5))
            assert exhaustive==aligned and pair_disagreements==literal*(5-literal)
            fixture_checks+=1
        assert len(rows)==998
        table=[]
        for row in rows:
            d=row['refinement_MAP_changes'];assert 0<=d<=600
            table.append(dict(role=row['role'],replicate=row['replicate'],literal_label_changes=d,
                label_invariant_membership_changes=min(d,600-d),pairwise_disagreements=d*(600-d),
                pure_label_reversal=d==600,score_change=row['refinement_score_change']))
        by_key={(r['role'],r['replicate']):r for r in rows}
        probe=json.loads(inputs['probe'].read_text());checks=[]
        for case in probe['cases']:
            key=(case['role'],case['replicate']);old=by_key[key]
            condition,_,_=model.ROLES[key[0]]
            y,reference,_,_=generator.generate(np.random.default_rng(np.random.SeedSequence(old['seed_key'])),condition,600,240)
            assert generator.array_hash(y)==old['input_sha256']['counts']
            corrected=model.fit_population(y,reference)
            assert file_record(old['arrays']['path'])['sha256']==old['arrays']['sha256']
            with np.load(old['arrays']['path']) as arrays:
                assert np.array_equal(corrected['probability'],arrays['probability'])
                assert np.array_equal(corrected['loading'],arrays['loading'])
            assert corrected['score']==old['score'] and corrected['coarse_score']==old['coarse_score']
            literal=old['refinement_MAP_changes']
            assert corrected['refinement_literal_label_changes']==literal
            assert corrected['refinement_MAP_changes']==min(literal,600-literal)
            assert corrected['refinement_pairwise_disagreements']==literal*(600-literal)
            checks.append(dict(role=key[0],replicate=key[1],scores_and_probabilities_exact=True,
                old_literal_changes=literal,corrected_partition_changes=corrected['refinement_MAP_changes']))
        frame=pd.DataFrame(table);frame.to_csv(args.out/'all_cases.tsv',sep='\t',index=False)
        frame.groupby('role').agg(cases=('replicate','size'),literal_changed=('literal_label_changes',lambda x:int((x>0).sum())),
            membership_changed=('label_invariant_membership_changes',lambda x:int((x>0).sum())),
            pure_reversals=('pure_label_reversal','sum')).to_csv(args.out/'by_role.tsv',sep='\t')
        report=dict(cases=len(frame),exhaustive_binary_fixture_checks=fixture_checks,exact_refits=checks,
            cases_with_literal_label_changes=int((frame.literal_label_changes>0).sum()),
            cases_with_label_invariant_membership_changes=int((frame.label_invariant_membership_changes>0).sum()),
            pure_label_reversal_cases=int(frame.pure_label_reversal.sum()),
            maximum_label_invariant_cell_changes=int(frame.label_invariant_membership_changes.max()),
            numerical_likelihood_sensitivity_resolved=False,scientific_criteria_changed=False,
            old_scores_and_specimen_decisions_changed=False,method_promoted=False)
        (args.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        for record in run['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        finish_run(args.out,run,backend='cpu')
    except Exception as error:finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
