"""Independent closed-form audit of the actual frozen bank-tail diagnostics."""
import argparse
from bisect import bisect_left
from fractions import Fraction
import json
import math
from pathlib import Path

from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    directory=a.task/'exact_v1'
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/'AUDIT_DESIGN.md',
                            exact_run=directory/'run_manifest.json'),{},
                  dict(scope='Independent exact arithmetic; does not repair the scientific criterion'),mode='joint_bank_exchangeability_audit')
    try:
        manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
        for kind in ('input','output'):
            for record in manifest[kind]['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        records=json.loads((directory/'score_file_identities.json').read_text());data={}
        for key,record in records.items():
            assert file_record(record['path'])['sha256']==record['sha256']
            arm,role,start=key.split(':');values=[]
            with Path(record['path']).open() as handle:
                for line in handle:
                    value=json.loads(line)
                    values.append((value['replicate'],value['score'],value['raw_candidate_count']))
            assert [r[0] for r in values]==list(range(int(start),int(start)+len(values)))
            data.setdefault((arm,role),[]).extend(values)
        result=json.loads((directory/'report.json').read_text());probabilities=[];checks=[]
        for row in result['diagnostics']:
            arm=row['arm'];state=row['condition']=='correlated_state'
            cal=sorted(r[1] for r in data[(arm,'cal_state' if state else 'cal_nb')])
            evaluation=data[(arm,'eval_state' if state else 'eval_nb')]
            assert len(cal)==199 and len(evaluation)==200
            limit=int(Fraction(str(row['alpha']))*200);h=limit-1
            observed=sum(s>0 and groups>1 and 1+len(cal)-bisect_left(cal,s)<=limit for _,s,groups in evaluation)
            assert observed==row['events']
            pooled=sorted(cal+[r[1] for r in evaluation],reverse=True)
            n=len(pooled);m=observed+h;denominator=math.comb(n,199)
            if observed==0:
                probability=Fraction(1);reference='P[K>=0]=1'
            else:
                # This condition is checked on actual scores. If a tie crosses
                # this boundary, the simple hypergeometric identity is invalid.
                assert pooled[m-1]>0 and pooled[m-1]>pooled[m]
                numerator=sum(math.comb(m,j)*math.comb(n-m,199-j)
                              for j in range(h+1) if 0<=199-j<=n-m)
                probability=Fraction(numerator,denominator)
                reference='Exact hypergeometric subset identity, strict positive score boundary'
            assert probability==Fraction(row['exact_tail'])
            assert sum(map(int,row['distribution'].values()))==denominator==int(row['allocation_denominator'])
            probabilities.append(probability)
            checks.append(dict(arm=arm,condition=row['condition'],alpha=row['alpha'],events=observed,
                               reference=reference,exact_tail=str(probability)))
        assert len(checks)==16
        ordered=sorted(probabilities)
        for probability,row in zip(probabilities,result['diagnostics']):
            adjusted=min(Fraction(1),max((16-rank)*value for rank,value in enumerate(ordered) if value<=probability))
            assert adjusted==Fraction(row['exact_holm_tail'])
        (a.out/'report.json').write_text(json.dumps(dict(checks=checks,score_files_verified=len(records),
                  frozen_events_reconstructed=16,exact_probabilities_verified=16,exact_holm_values_verified=16,
                  new_draws=0,scientific_criterion_changed=False,biological_validation=False),indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
