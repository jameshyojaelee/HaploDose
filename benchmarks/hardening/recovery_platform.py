"""Fixed-node numerical replay and actual transitive-library evidence."""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def load():
    path=Path(__file__).with_name('recovery_kernel_replay.py')
    spec=importlib.util.spec_from_file_location('platform_kernel_replay',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def libraries(out):
    records=[];archive=out/'system_libraries';archive.mkdir(exist_ok=False)
    paths={line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
        if len(line.split())>=6 and line.split()[-1].startswith('/')}
    for text in sorted(paths):
        path=Path(text)
        if not path.is_file() or '.so' not in path.name:continue
        digest=hashlib.sha256(path.read_bytes()).hexdigest()
        system=path.is_relative_to('/lib') or path.is_relative_to('/lib64') or path.is_relative_to('/usr/lib')
        row=dict(path=str(path),sha256=digest,bytes=path.stat().st_size,system=system)
        if system:
            names={path.name}
            result=subprocess.run(['readelf','-d',str(path)],capture_output=True,text=True,check=True)
            match=re.search(r'\(SONAME\).*?\[([^\]]+)\]',result.stdout)
            if match:names.add(match[1])
            for name in names:
                target=archive/name
                if target.exists():assert hashlib.sha256(target.read_bytes()).hexdigest()==digest
                else:shutil.copyfile(path,target)
            row['archive_names']=sorted(names)
        records.append(row)
    write_json(dict(libraries=records,kernel=platform.uname()._asdict(),
        os_release=Path('/etc/os-release').read_text(),executable=os.readlink('/proc/self/exe')),
        out/'platform.json')


def files(left,right):
    result=[]
    for path in sorted(left.iterdir()):
        other=right/path.name
        if path.suffix=='.npy':
            a,b=np.load(path),np.load(other)
            try:np.testing.assert_array_equal(a,b);exact=True
            except AssertionError:exact=False
            result.append(dict(file=path.name,exact=exact,shape_left=list(a.shape),shape_right=list(b.shape)))
        elif path.name.endswith(('.tsv','.tsv.gz')):
            a,b=pd.read_csv(path,sep='\t'),pd.read_csv(other,sep='\t')
            try:pd.testing.assert_frame_equal(a,b,check_exact=True);exact=True
            except AssertionError:exact=False
            result.append(dict(file=path.name,exact=exact,shape_left=list(a.shape),shape_right=list(b.shape)))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    parser.add_argument('--worker',action='store_true')
    parser.add_argument('--right-name',choices=('node010_v1','runtime003_on010_v1'),default='node010_v1')
    args=parser.parse_args()
    if args.worker:
        if os.environ.get('PYTHONHASHSEED')!='0':raise RuntimeError('Fixed hash seed required before interpreter startup')
        module=load();module.worker(args.out,args.task);libraries(args.out);return
    args.out.mkdir(exist_ok=False)
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'DESIGN.md',
        observer=Path(__file__).with_name('recovery_kernel_replay.py')),
        dict(nodes=['ne1dc6-003','ne1dc6-010']),dict(scope='Fixed platform mechanism contrast, no gate replacement'),mode='recovery_platform')
    try:
        left,right=args.task/'node003_v1',args.task/args.right_name;module=load()
        comparison=module.compare(left,right)
        root=args.task.parents[1]
        historical=root/'dev/depth_stratified_gene_noise_20261109/both_v1/panel/seed_20261108/fit'
        recent=root/'dev/depth_conditioned_genes_20261109/control_v1/candidate/panel/seed_20261108/fit'
        final={}
        for name,path in [('node003',left),(args.right_name.removesuffix('_v1'),right)]:
            fit=path/'seed_20261108/fit'
            final[name]=dict(historical=files(historical,fit),recent=files(recent,fit))
        a,b=[json.loads((p/'platform.json').read_text()) for p in (left,right)]
        aa={Path(r['path']).name:r for r in a['libraries']};bb={Path(r['path']).name:r for r in b['libraries']}
        changed=[dict(name=k,left=aa[k],right=bb[k]) for k in sorted(aa.keys()&bb.keys()) if aa[k]['sha256']!=bb[k]['sha256']]
        write_json(dict(numerical_comparison=comparison,final_file_comparisons=final,
            loaded_library_differences=changed,loaded_only_left=sorted(aa.keys()-bb.keys()),loaded_only_right=sorted(bb.keys()-aa.keys()),
            exact_array_assertion='Original np.testing.assert_array_equal; matching NaNs are equal; no tolerance',
            historical_exact_gate_replaced=False,scientific_criteria_changed=False,production_changed=False),args.out/'report.json')
        finish_run(args.out,run,backend='CPU exact traces and transitive system-library comparison')
    except BaseException as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
