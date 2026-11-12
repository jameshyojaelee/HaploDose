"""Check installed numerical correction against its same-processor research hook."""
import argparse
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
import jax
from haplodose import clones
from haplodose.provenance import start_run,finish_run,file_record


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--trace-task',type=Path,required=True);p.add_argument('--host',choices=['intel','amd'],required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    if jax.default_backend()!='cpu':raise RuntimeError('CPU parity requires the actual CPU inference backend')
    source=a.trace_task/f'{a.host}_precision_v1'
    parent=json.loads((source/'run_manifest.json').read_text())
    if parent['status']!='complete':raise ValueError('Incomplete research reference')
    inputs=dict(script=__file__,reference_manifest=source/'run_manifest.json')
    for replicate in (150,167):
        for name in ('input.npz','result.npz','summary.json'):
            inputs[f'{replicate}:{name}']=source/f'replicate_{replicate}'/name
    params=clones.CloneParams(n_programmes=0,spatial_neighbors=0,seed=20261107)
    rec=start_run(a.out,inputs,asdict(params),dict(scope='Installed-wheel versus research-hook identity on two locked cases',
                    host=a.host,backend=str(jax.devices())),mode='precision_wheel_parity')
    try:
        rows=[]
        for replicate in (150,167):
            with np.load(inputs[f'{replicate}:input.npz']) as z:
                Y,base,phi,reset=(z[k].copy() for k in ('counts','reference','phi','reset'))
            result=clones.fit_specimen(clones.Rows(Y,np.arange(len(Y))),Y.sum(axis=1,dtype=np.float64),
                                      base[:,None],reset,phi,params)
            with np.load(inputs[f'{replicate}:result.npz']) as old:
                np.testing.assert_array_equal(result['R'],old['R'])
                np.testing.assert_array_equal(result['doses'],old['doses'])
            np.savez_compressed(a.out/f'{replicate}.npz',R=result['R'],doses=result['doses'])
            rows.append(dict(replicate=replicate,posterior_exact=True,dose_exact=True))
        for record in rec['input']['files'].values():
            if file_record(record['path'])['sha256']!=record['sha256']:raise ValueError('Reference changed')
        (a.out/'summary.json').write_text(json.dumps(dict(cases=rows,scope='Same-processor artifact identity, not biological validation'),indent=2))
        finish_run(a.out,rec,backend='cpu')
    except Exception as error:
        finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
