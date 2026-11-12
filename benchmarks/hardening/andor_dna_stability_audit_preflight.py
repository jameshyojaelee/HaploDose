"""Hand-calculable checks of the independent DNA count and callable-flow audit."""
import argparse
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    source=Path(__file__).with_name('andor_dna_stability_audit.py')
    spec=importlib.util.spec_from_file_location('audit',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    run=start_run(args.out,dict(script=Path(__file__),reference=source,
                  design=args.task/'AUDIT_DESIGN.md'),{},dict(scope='Integer audit fixtures, biological n0'),
                  mode='andor_dna_stability_audit_preflight')
    try:
        bins=pd.DataFrame(dict(chromosome=['chr1']*5,bin=range(5)))
        counts=args.out/'counts.tsv'
        counts.write_text('002\tchr1\t2\t5\n001\tchr1\t0\t2\n002\tchrX\t0\t7\n')
        matrix,totals=module.independent_counts(counts,['001','002'],bins)
        np.testing.assert_array_equal(matrix,[[2,0,0,0,0],[0,0,5,0,0]])
        assert dict(totals)=={('002','chr1'):5,('001','chr1'):2,('002','chrX'):7}
        ann=pd.DataFrame(dict(chromosome=['chr1']*4,bin=range(4),
                             n_frac=[1,0,0,1],blacklist=[0]*4,mappability=[1]*4))
        old=np.array([[2,3,4,5,6]],dtype=np.int64)
        flow=np.array([[0,1,2,3,4]],dtype=np.int64)
        new=old-flow+np.array([[1,2,3,4,0]])
        result=module.callable_flow(old,new,flow,bins,ann)
        assert {key:int(value[0]) for key,value in result.items()}==dict(
            entering_callable=3,leaving_callable=1,staying_callable=2,staying_outside_callable=4)
        rejected=0
        for bad in (flow+1,flow[:,::-1],np.zeros_like(flow)):
            try:
                module.callable_flow(old,new,bad,bins,ann)
            except AssertionError:
                rejected+=1
            else:
                raise AssertionError('Invalid boundary-flow evidence accepted')
        counts.write_text('001\tchr1\t0\t2\n001\tchr1\t0\t3\n')
        try:
            module.independent_counts(counts,['001','002'],bins)
        except AssertionError:
            pass
        else:
            raise AssertionError('Ambiguous duplicate count accepted')
        write_json(dict(integer_reader_exact=True,literal_IDs_preserved=True,
                        callable_flow_hand_calculated=True,invalid_flows_rejected=rejected,
                        duplicate_rejected=True,biological_specimens=0),args.out/'report.json')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
