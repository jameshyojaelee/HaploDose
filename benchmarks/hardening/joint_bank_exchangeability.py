"""Exact finite-bank diagnostic conditional on frozen null-score multisets."""
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from haplodose.provenance import start_run,finish_run,file_record


def choose(n,k):
    return math.comb(n,k) if 0<=k<=n else 0


def allocation_counts(scores,ncal,h):
    """Count all calibration subsets; equal-score blocks enter the rank together.

    A positive evaluation score is accepted iff its number of calibration scores
    greater than or equal to it is <=h. No continuous-score approximation.
    """
    n=len(scores);blocks=Counter(x for x in scores if x>0)
    active={0:1};processed=0;result=Counter()
    for score,b in sorted(blocks.items(),reverse=True):
        following={}
        for previous,ways in active.items():
            for current in range(min(b,ncal-previous)+1):
                multiplicity=ways*choose(b,current);total=previous+current
                if total>h:
                    result[processed-previous]+=multiplicity*choose(n-processed-b,ncal-total)
                else:
                    following[total]=following.get(total,0)+multiplicity
        active=following;processed+=b
    for cal,ways in active.items():
        result[processed-cal]+=ways*choose(n-processed,ncal-cal)
    result=Counter({k:v for k,v in result.items() if v})
    assert sum(result.values())==choose(n,ncal)
    return result


def exhaustive(scores,ncal,h):
    result=Counter()
    for subset in itertools.combinations(range(len(scores)),ncal):
        calibration=[scores[i] for i in subset];selected=set(subset)
        accepted=sum(value>0 and sum(c>=value for c in calibration)<=h
                     for i,value in enumerate(scores) if i not in selected)
        result[accepted]+=1
    return result


def checks():
    cases=0
    for scores in ([0,0,0,0],[-1,0,1,1,2,3],[1,2,3,4,5,6],[0,1,1,1,2,2,3,4]):
        for ncal in range(1,len(scores)):
            for h in range(ncal):
                distribution=allocation_counts(scores,ncal,h)
                assert distribution==exhaustive(scores,ncal,h)
                cases+=1
                if len(set(scores))==len(scores) and min(scores)>0:
                    n=len(scores);denominator=choose(n,ncal)
                    for s in range(1,n-ncal+1):
                        # K>=s iff the first s+h ranks contain <=h calibration
                        # labels: independent hypergeometric counting identity.
                        m=s+h
                        tail=sum(choose(m,j)*choose(n-m,ncal-j) for j in range(h+1))
                        assert sum(v for k,v in distribution.items() if k>=s)==tail
    return cases


def read_scores(task,arm,role):
    scores=[];inputs={}
    for start in (0,50,100,150):
        directory=task/'chunks_v1'/arm/f'{role}_{start}'
        manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
        record=manifest['output']['files']['records.jsonl'];assert file_record(record['path'])['sha256']==record['sha256']
        inputs[f'{arm}:{role}:{start}']=Path(record['path'])
        with Path(record['path']).open() as handle:
            for line in handle:
                row=json.loads(line)
                score=row['score'];assert math.isfinite(score)
                assert score<=0 or row['raw_candidate_count']>1
                scores.append(dict(score=score,groups=row['raw_candidate_count'],replicate=row['replicate']))
    assert [r['replicate'] for r in scores]==list(range(len(scores)))
    return scores,inputs


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--frozen',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/'DESIGN.md',
                            summary=a.frozen/'summary_v1/summary.json',audit=a.frozen/'audit_v1/report.json'),
                  dict(diagnostic_family=16,calibration=199,evaluation=200,alpha=[.05,.005,.01,.025]),
                  dict(scope='Post hoc cause diagnostic; unchanged failed scientific family153 and all frozen decisions'),
                  mode='joint_bank_exchangeability')
    try:
        for name in ('summary_v1','audit_v1'):
            manifest=json.loads((a.frozen/name/'run_manifest.json').read_text())
            assert manifest['status']=='complete'
            for record in manifest['output']['files'].values():
                assert file_record(record['path'])['sha256']==record['sha256']
        fixtures=checks();results=[];files={}
        summary=json.loads((a.frozen/'summary_v1/summary.json').read_text())
        assert summary['family']==153
        for arm in ('joint_proposal','joint_full'):
            for condition,cal_role,eval_role,method in (
                    ('matched_nb','cal_nb','eval_nb','NB_calibrated_gate'),
                    ('correlated_state','cal_state','eval_state','state_calibrated_gate')):
                cal,records=read_scores(a.frozen,arm,cal_role);files.update(records)
                evaluation,records=read_scores(a.frozen,arm,eval_role);files.update(records)
                assert len(cal)==199 and len(evaluation)==200
                for alpha in (.05,.005,.01,.025):
                    limit=Fraction(str(alpha))*200
                    assert limit.denominator==1
                    h=int(limit)-1
                    observed=sum(row['score']>0 and row['groups']>1
                                 and 1+sum(c['score']>=row['score'] for c in cal)<=int(limit)
                                 for row in evaluation)
                    prior=[r for r in summary['summary'] if (r['arm'],r['condition'],r['method'],r['alpha'])==(arm,condition,method,alpha)]
                    assert len(prior)==1 and prior[0]['events']==observed
                    distribution=allocation_counts([r['score'] for r in cal+evaluation],199,h)
                    denominator=choose(399,199)
                    numerator=sum(v for k,v in distribution.items() if k>=observed)
                    probability=Fraction(numerator,denominator)
                    results.append(dict(arm=arm,condition=condition,alpha=alpha,events=observed,
                                        exact_tail=str(probability),tail_float=float(probability),
                                        distribution={str(k):str(v) for k,v in sorted(distribution.items())},
                                        allocation_denominator=str(denominator)))
        assert len(results)==16
        probabilities=[Fraction(r['exact_tail']) for r in results];maximum=Fraction(0)
        for rank,i in enumerate(sorted(range(16),key=lambda i:probabilities[i])):
            maximum=max(maximum,min(Fraction(1),probabilities[i]*(16-rank)))
            results[i].update(exact_holm_tail=str(maximum),holm_tail=float(maximum))
        (a.out/'score_file_identities.json').write_text(json.dumps({key:file_record(path) for key,path in files.items()},indent=2)+'\n')
        (a.out/'report.json').write_text(json.dumps(dict(exhaustive_fixtures=fixtures,diagnostics=results,
                    diagnostic_family=16,original_scientific_family=153,biological_specimens=0,
                    hypothesis='Exchangeable workflow scores under the same stated null law; conditional on pooled scores',
                    limitation='Failure to reject is not proof of exchangeability or biological validity; no correction promoted',
                    new_inference=False,new_draws=0,frozen_decisions_unchanged=True),indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
