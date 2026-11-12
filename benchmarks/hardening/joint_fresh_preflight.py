"""Independent bounds for every possible event count, before new outcomes."""
import argparse
import importlib.util
import json
from pathlib import Path
from haplodose.provenance import start_run,finish_run


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(exist_ok=False)
    path=Path(__file__).with_name('joint_fresh_summary.py');spec=importlib.util.spec_from_file_location('fresh_numerics',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rec=start_run(a.out,dict(script=Path(__file__),statistics=path),dict(n=200,family=153),
        dict(scope='Pure numerical check over all event counts, no inference outcomes'),mode='fresh_bound_preflight')
    try:
        checks=[]
        for power in (False,True):
            previous=-1.
            for k in range(201):
                value=module.bound(k,200,power);lo,hi=module.independent_enclosure(k,200,power)
                if not lo<=value<=hi or value<previous:raise AssertionError(f'Bound check failed at {k},power={power}: {lo},{value},{hi}')
                lower_precision=module.bound(k,200,power,precision=60)
                if value!=lower_precision:raise AssertionError('Decimal precision did not converge after float conversion')
                previous=value;checks.append(dict(k=k,power=power,lo=lo,value=value,hi=hi,
                    scipy_inverse_beta_value=module.scipy_bound(k,200,power),decimal60_value=lower_precision))
        (a.out/'checks.json').write_text(json.dumps(checks,indent=2));finish_run(a.out,rec,backend='cpu')
        print(json.dumps(dict(independent_binomial_bounds=len(checks),status='passed')),flush=True)
    except Exception as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
