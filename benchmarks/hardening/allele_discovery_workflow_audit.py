"""Independent reconstruction of complete-discovery SNP isolation artifacts."""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run, file_record


def fingerprint(directory):
    with np.load(directory/'discovery.npz') as data:
        arrays={name:data[name] for name in data.files}
    digests={}
    for name,value in arrays.items():
        value=np.ascontiguousarray(value)
        digests[name]=hashlib.sha256(str(value.dtype).encode()+json.dumps(value.shape).encode()+value.tobytes()).hexdigest()
    meta=json.loads((directory/'discovery.json').read_text())
    assert file_record(directory/'search_log.json')['sha256']==meta['search_log_sha256']
    result=dict(arrays=digests,metadata=hashlib.sha256(json.dumps(meta,sort_keys=True,allow_nan=False).encode()).hexdigest())
    assert result==json.loads((directory/'fingerprint.json').read_text())
    return arrays,meta,result


def partition(r):
    labels=r.argmax(axis=1)
    return sorted(tuple(np.flatnonzero(labels==k)) for k in np.unique(labels))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md')
    for model in ('binomial','beta'):
        directory=args.task/(model+'_v1')
        inputs[model+':manifest']=directory/'run_manifest.json'
        manifest=json.loads(inputs[model+':manifest'].read_text())
        assert manifest['status']=='complete'
        for group in ('input','output'):
            for record in manifest[group]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        inputs[model+':cases']=directory/'cases.json'
        inputs[model+':report']=directory/'report.json'
    run=start_run(args.out,inputs,{},dict(scope='Computational isolation, not error/power or biological validation'),
                  mode='allele_discovery_workflow_audit')
    try:
        reports=[];contrasts_checked=0;arrays_checked=0
        for model in ('binomial','beta'):
            root=args.task/(model+'_v1')
            cases=json.loads(inputs[model+':cases'].read_text())
            assert len(cases)==8
            baseline={};count_summaries=[]
            for case in cases:
                name,arm=case['scenario'],case['arm']
                directory=root/(name+'_'+arm)
                arrays,meta,value=fingerprint(directory)
                assert value==case['discovery']
                arrays_checked+=len(arrays)
                selected=np.array(meta['selected_snp_indices'])
                assert selected.tolist()==list(range(0,40,2)) if arm=='discovery_subset' else selected.tolist()==list(range(40))
                with np.load(root/(name+'_alleles.npz')) as data:
                    a,d=data['A'],data['D']
                labels=arrays['R'].argmax(axis=1)
                group_sizes=np.bincount(labels,minlength=arrays['R'].shape[1]).tolist()
                assert group_sizes==case['group_sizes']
                eligible=[k for k,n in enumerate(group_sizes) if n>=30]
                reserved=json.loads((directory/'reserved.json').read_text())
                expected_pairs=list(itertools.combinations(eligible,2))
                assert len(reserved['contrasts'])==len(meta['regions'])*len(expected_pairs)
                positions=np.arange(40)*1000+1
                for contrast in reserved['contrasts']:
                    region=contrast['region']
                    assert region in meta['regions']
                    sites=[j for j in range(1,40,2) if region['start']<positions[j]<=region['end']]
                    assert sites==contrast['snps']
                    pair=tuple(contrast['candidate_pair'])
                    assert pair in expected_pairs
                    aa=[];dd=[]
                    for k in pair:
                        cells=np.flatnonzero(labels==k)
                        assert len(cells)==contrast['n_cells'][len(aa)]
                        aa.append([sum(float(a[c,j]) for c in cells) for j in sites])
                        dd.append([sum(float(d[c,j]) for c in cells) for j in sites])
                    assert aa==contrast['alternate_counts'] and dd==contrast['total_counts']
                    assert contrast['status']=='raw_reserved_summary_only'
                    contrasts_checked+=1
                if not reserved['contrasts']:
                    assert reserved['reason']=='no_frozen_region_and_eligible_pair'
                # Independent dense weighted sums check the separate all-SNP
                # descriptive pooling, with a cell-count-scaled summation bound.
                with np.load(directory/'descriptive_all_snp.npz') as descriptive:
                    for key,data in (('Ak',a),('Dk',d)):
                        expected=arrays['R'].astype(float).T@data
                        unit=np.finfo(float).eps/2
                        gamma=120*unit/(1-120*unit)
                        allowance=2*gamma*(1+abs(expected))
                        assert np.all(abs(descriptive[key]-expected)<=allowance)
                key=(arm,name.startswith('absent'))
                if key not in baseline:
                    baseline[key]=(value,arrays['R'])
                unchanged=value==baseline[key][0]
                if arm=='discovery_subset':
                    assert unchanged
                count_summaries.append(dict(scenario=name,arm=arm,group_sizes=group_sizes,
                    eligible_groups=len(eligible),reserved_contrasts=len(reserved['contrasts']),
                    discovery_identical=unchanged,
                    unlabelled_partition_identical=partition(arrays['R'])==partition(baseline[key][1])))
                # Reconstruct arrays AFTER observing the descriptive/reserved
                # outputs; verify agreement with the pre-stage fingerprint.
                assert fingerprint(directory)[2]==value
            repeat=fingerprint(root/'baseline_repeat')[2]
            assert repeat==baseline[('discovery_subset',False)][0]
            reports.append(dict(cell_model=model,cases=count_summaries,baseline_repeat_identical=True))
        write_json(dict(full_discovery_fits_checked=18,case_arrays_checked=arrays_checked,
                        reserved_contrasts_independently_summed=contrasts_checked,models=reports,
                        conditional_binned_isolation=True,public_raw_workflow_validated=False,
                        biological_validation=False,calibrated=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
