"""Exact fixed-cohort mean-score bounds, preserving missing competitor outputs."""
import argparse
import csv
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def bounds(primary, comparator):
    if not primary or len(primary)!=len(comparator):
        raise ValueError('Nonempty aligned populations required')
    lo=[];hi=[]
    for a,b in zip(primary,comparator):
        if not math.isfinite(a) or not -1<=a<=1:
            raise ValueError('Complete finite primary correlation required')
        a=Fraction.from_float(a)
        if b is None:
            lo.append(a-1);hi.append(a+1)
        else:
            if not math.isfinite(b) or not -1<=b<=1:
                raise ValueError('Comparator correlation outside its mathematical range')
            difference=a-Fraction.from_float(b);lo.append(difference);hi.append(difference)
    return sum(lo)/len(lo),sum(hi)/len(hi)


def outward(value, direction):
    rounded=float(value)
    if (direction<0 and Fraction.from_float(rounded)>value or
            direction>0 and Fraction.from_float(rounded)<value):
        rounded=math.nextafter(rounded,-math.inf if direction<0 else math.inf)
    return rounded


def checks():
    primary=[.25,.75,.5]; comparator=[.125,None,None]
    lo,hi=bounds(primary,comparator)
    values=[]
    for a,b in itertools.product((-1.,0.,1.),repeat=2):
        values.append(bounds(primary,[.125,a,b])[0])
    assert lo==min(values) and hi==max(values)
    assert bounds(primary,primary)==(0,0)
    for value in (lo,hi,Fraction(1,3),-Fraction(1,3)):
        assert Fraction.from_float(outward(value,-1))<=value
        assert Fraction.from_float(outward(value,1))>=value
    return dict(exhaustive_endpoint_completion=True,complete_identity=True,outward_rounding=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--task',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    root=args.task.parents[1]
    field=root/'dev/gene_laplace_profile_20261110/field_v1'
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md',
                  scores=field/'scores.tsv',field_report=field/'report.json',
                  field_manifest=field/'run_manifest.json',
                  population_coverage=root/'dev/native_stacked_profile_20261110/field_v2/population_coverage.tsv'),
                  dict(primary='gene_laplace_stack',fixed_cases=181,fixed_contrasts=26),
                  dict(scope='Post hoc fixed-cohort deterministic sensitivity, no new significance or promotion test'),
                  mode='field_missing_bounds')
    try:
        assert json.loads((field/'run_manifest.json').read_text())['status']=='complete'
        previous=json.loads((field/'report.json').read_text())
        assert previous['family_size']==26 and previous['previous_scores_exact']==4344
        with (field/'scores.tsv').open()as stream:
            rows=list(csv.DictReader(stream,delimiter='\t'))
        assert len(rows)==5068
        lookup={}
        for row in rows:
            key=(row['case'],row['method'],row['endpoint']);assert key not in lookup
            value=float(row['score']) if row['score'] else None
            lookup[key]=value if value is not None and math.isfinite(value) else None
        cases=sorted({r['case']for r in rows});methods=sorted({r['method']for r in rows})
        assert len(cases)==181 and len(methods)==14 and 'gene_laplace_stack'in methods
        results=[]
        for endpoint in ('pearson_1mb','spearman_arm'):
            primary=[lookup[c,'gene_laplace_stack',endpoint]for c in cases]
            assert all(v is not None for v in primary)
            for method in methods:
                if method=='gene_laplace_stack':continue
                comparator=[lookup[c,method,endpoint]for c in cases]
                lo,hi=bounds(primary,comparator)
                results.append(dict(endpoint=endpoint,comparator=method,cases=len(cases),
                    known_scores=sum(v is not None for v in comparator),
                    missing_cases=[c for c,v in zip(cases,comparator)if v is None],
                    mean_difference_lower=outward(lo,-1),mean_difference_upper=outward(hi,1),
                    exact_lower=str(lo),exact_upper=str(hi),
                    positive_for_every_defined_completion=bool(lo>0),
                    population_superiority_test=False))
        assert len(results)==26
        write_json(dict(checks=checks(),contrasts=results,existing_family_size=26,
            original_formal_verdicts_unchanged=True,original_coverage_failures_unchanged=True,
            intervals_are='Exact fixed-cohort identification ranges over hypothetical defined missing scores, not confidence intervals',
            biological_replication_added=0,independent_evaluation=False,
            scientific_criteria_changed=False,default_promotion=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU exact rational arithmetic on frozen field scores')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
