"""Independently check held-out perturbations and unlabelled cell partitions."""
import argparse
import json
from pathlib import Path

import numpy as np

from haplodose.provenance import start_run, finish_run, file_record


def canonical(proposals):
    result=[]
    for proposal in proposals:
        labels=np.asarray(proposal['labels'])
        assert labels.shape==(120,) and set(labels)=={0,1}
        groups=tuple(sorted(tuple(np.flatnonzero(labels==g).tolist()) for g in (0,1)))
        result.append(groups)
    return sorted(result)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    source=args.task/'probe_v1'
    manifest=json.loads((source/'run_manifest.json').read_text())
    assert manifest['status']=='complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
    inputs=dict(script=Path(__file__),manifest=source/'run_manifest.json',
                report=source/'report.json',details=source/'details.json',posterior=source/'posteriors.npz')
    for name in ('odd_balanced','odd_same_orientation','odd_reverse_orientation'):
        inputs[name]=source/(name+'_inputs.npz')
    run=start_run(args.out,inputs,{},dict(scope='Local computational leakage; not global independent confirmation'),
                  mode='allele_confirmation_leakage_audit')
    try:
        names=['odd_balanced','odd_same_orientation','odd_reverse_orientation']
        baseline=None
        for i,name in enumerate(names):
            with np.load(inputs[name]) as data:
                a,d=data['A'],data['D']
            assert a.shape==d.shape==(120,40) and np.array_equal(d,np.full((120,40),20.))
            assert np.array_equal(a[:,::2],np.repeat(np.repeat([10.,12.,16.],40)[:,None],20,axis=1))
            assert np.array_equal(a[:,1::2],np.full((120,20),[10.,19.,1.][i]))
            if baseline is None:
                baseline=a[:,::2].copy()
            assert np.array_equal(a[:,::2],baseline)
        details=json.loads(inputs['details'].read_text());summary=[]
        with np.load(inputs['posterior']) as post:
            for arm in ('current','proposal_only_reference'):
                first=details[names[0]+':'+arm]
                partition=canonical(first['proposals'])
                for name in names:
                    record=details[name+':'+arm]
                    candidate=canonical(record['proposals'])
                    change=not np.array_equal(post[name+':'+arm],post[names[0]+':'+arm])
                    if arm=='proposal_only_reference':
                        assert record==first and not change
                    summary.append(dict(arm=arm,scenario=name,proposals=len(record['proposals']),
                        group_sizes=[sorted([len(g) for g in pair]) for pair in candidate],
                        unlabelled_partitions_identical=candidate==partition,
                        posterior_changed=change,kappa=record['kappa'],nu=record['nu']))
        (args.out/'report.json').write_text(json.dumps(dict(scenarios_verified=3,
            proposal_inputs_identical=True,comparisons=summary,
            local_partition_leakage=any(r['arm']=='current' and not r['unlabelled_partitions_identical'] for r in summary),
            subset_reference_invariant=True,biological_validation=False,
            full_workflow_independence_established=False),indent=2))
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
