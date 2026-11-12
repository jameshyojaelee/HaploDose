"""Independently reconstruct held-out ASE contrasts and permutation decisions."""
import argparse
from decimal import Decimal,localcontext
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def fingerprint(array):
    array=np.ascontiguousarray(array)
    return hashlib.sha256(str(array.dtype).encode()+json.dumps(array.shape).encode()+array.tobytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    inputs=dict(script=Path(__file__),design=args.task/'DESIGN.md')
    for model in ('binomial','beta'):
        inputs[model]=args.task/(model+'_v1/run_manifest.json')
    run=start_run(args.out,inputs,{},dict(scope='Numerical/identity audit, no biological validation'),mode='allele_shared_state_audit')
    try:
        cases_checked=0;permutations_checked=0;maximum_discrepancy=0.;summaries=[]
        with localcontext() as context:
            context.prec=80
            zero_upper=float(1-((Decimal(1)/120).ln()/40).exp())
            all_lower=float(((Decimal(1)/120).ln()/40).exp())
        for model in ('binomial','beta'):
            root=args.task/(model+'_v1');manifest=json.loads(inputs[model].read_text())
            assert manifest['status']=='complete' and manifest['backend']=='cpu'
            for group in ('input','output'):
                for record in manifest[group]['files'].values():
                    assert file_record(record['path'])['sha256']==record['sha256']
            cases=json.loads((root/'cases.json').read_text());report=json.loads((root/'report.json').read_text())
            roles=json.loads((root/'roles.json').read_text());independent=[]
            for s in range(40):
                payload='\0'.join(['haplodose-allele-role-v1','20261108','hg38','chr1',str(s*1000+1)])
                if int(hashlib.blake2b(payload.encode()).hexdigest()[:16],16)<2**63:
                    independent.append(s)
            assert roles['discovery']==independent
            assert sorted(roles['discovery']+roles['reserved'])==list(range(40))
            for case in cases:
                directory=root/f"{case['replicate']:03d}_{case['condition']}"
                test=json.loads((directory/'reserved_test.json').read_text())
                meta=json.loads((directory/'discovery.json').read_text())
                arrays=np.load(directory/'discovery.npz');data=np.load(directory/'audit_inputs.npz')
                assert {name:fingerprint(arrays[name]) for name in arrays.files}==case['frozen_discovery']['arrays']
                assert hashlib.sha256(json.dumps(meta,sort_keys=True,allow_nan=False).encode()).hexdigest()==case['frozen_discovery']['metadata']
                assert file_record(directory/'search_log.json')['sha256']==meta['search_log_sha256']
                assert np.array_equal(arrays['R'].argmax(axis=1),data['labels'])
                assert np.all(data['D']==20) and np.all(data['A']==np.floor(data['A']))
                assert np.all((data['A']>=0)&(data['A']<=data['D']))
                other=args.task/(('beta' if model=='binomial' else 'binomial')+'_v1')/directory.name
                with np.load(other/'audit_inputs.npz') as paired:
                    for name in ('Y','A','D','truth'):
                        assert np.array_equal(data[name],paired[name])
                assert np.array_equal(data['permutations'][0],np.arange(120))
                assert np.array_equal(np.sort(data['permutations'],axis=1),np.broadcast_to(np.arange(120),(1000,120)))
                labels=data['labels'];groups=np.bincount(labels,minlength=arrays['R'].shape[1])
                eligible=[int(k) for k,n in enumerate(groups) if n>=30]
                assert eligible==test['eligible'] and list(map(list,itertools.combinations(eligible,2)))==test['pairs']
                regions=[]
                for region in meta['regions']:
                    sites=[s for s in roles['reserved'] if region['chromosome']=='chr1' and region['start']<s*1000+1<=region['end']]
                    if sites:
                        regions.append(dict(region=region,sites=sites))
                assert regions==test['regions']
                fractions=np.empty((120,len(regions)))
                for c in range(120):
                    for r,region in enumerate(regions):
                        numerator=sum(int(data['A'][c,s]) for s in region['sites'])
                        denominator=sum(int(data['D'][c,s]) for s in region['sites'])
                        fractions[c,r]=numerator/denominator
                assert np.array_equal(fractions,data['fractions'])
                maxima=[]
                for permutation in data['permutations']:
                    largest=0.
                    for left,right in test['pairs']:
                        for r in range(len(regions)):
                            lm=sum(float(fractions[permutation[c],r]) for c in range(120) if labels[c]==left)/int(groups[left])
                            rm=sum(float(fractions[permutation[c],r]) for c in range(120) if labels[c]==right)/int(groups[right])
                            largest=max(largest,abs(lm-rm))
                    maxima.append(largest)
                maxima=np.array(maxima)
                # Scalar summation versus vectorized mean: gamma_(n-1) error
                # for means of fractions in [0,1], twice for a difference.
                epsilon=np.finfo(float).eps/2
                allowance=4*119*epsilon/(1-119*epsilon)
                discrepancy=float(np.max(np.abs(maxima-data['maxima'])))
                assert discrepancy<=allowance
                maximum_discrepancy=max(maximum_discrepancy,discrepancy)
                available=bool(test['pairs'] and regions)
                # Intersections with the rounding envelope cannot be silently
                # classified. Exact all-zero abstentions are explicit ties.
                if available:
                    distances=maxima[1:]-maxima[0]
                    assert np.all(np.abs(distances)>2*allowance), 'Unresolved permutation tie near floating boundary'
                    p=(1+sum(value>=maxima[0] for value in maxima[1:]))/1000
                else:
                    assert np.all(maxima==0);p=1.
                assert p==test['p']==case['p'] and available==case['available']
                assert case['positive']==(p<=.05)
                truth=data['truth'];assert np.bincount(truth).tolist()==[60,60]
                table=[[sum(int(labels[c]==k and truth[c]==t) for c in range(120)) for t in (0,1)] for k in range(len(groups))]
                matched=max((table[left][0]+table[right][1] for left,right in itertools.permutations(range(len(groups)),2)),default=max(table[0]))
                assert matched/120==case['latent_all_cell_matching']
                assert groups.tolist()==case['group_sizes']
                cases_checked+=1;permutations_checked+=len(maxima)
            for summary in report['summaries']:
                selected=[case for case in cases if case['condition']==summary['condition']]
                hits=sum(case['positive'] for case in selected);assert hits in (0,40)
                assert summary['n']==40 and summary['positives']==hits
                assert summary['available']==sum(case['available'] for case in selected)
                lower=all_lower if hits==40 else 0.;upper=zero_upper if hits==0 else 1.
                for name,value in (('simultaneous_lower',lower),('simultaneous_upper',upper)):
                    assert abs(summary[name]-value)<=8*abs(np.spacing(value))
                assert summary['mean_latent_matching']==float(np.mean([case['latent_all_cell_matching'] for case in selected]))
                summaries.append(dict(model=model,**summary))
        write_json(dict(cases_checked=cases_checked,permutation_statistics_checked=permutations_checked,
                        maximum_scalar_vector_discrepancy=maximum_discrepancy,
                        simultaneous_extreme_bounds_independent_decimal=True,summaries=summaries,
                        biological_validation=False,production_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
