"""Stage only exact existing conda RNA-tool dependencies, offline and by copy."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import urlparse

from haplodose.provenance import start_run, finish_run, file_record

SOURCE=Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    mamba=Path(shutil.which('micromamba')).resolve()
    run=start_run(args.out,dict(script=Path(__file__),design=args.task/'TOOLS_DESIGN.md',micromamba=mamba),
                  dict(source_prefix=str(SOURCE),offline=True,copy=True),
                  dict(scope='Isolated tool staging only, not alignment or biological validation'),mode='scone_tool_environment')
    try:
        before={p.name:file_record(p)['sha256'] for p in sorted((SOURCE/'conda-meta').glob('*')) if p.is_file()}
        installed={}
        for path in (SOURCE/'conda-meta').glob('*.json'):
            item=json.loads(path.read_text());name=item['name'];assert name not in installed
            installed[name]=item
        roots={'star':'2.7.11b','samtools':'1.22.1','subread':'2.1.1'}
        assert all(installed[name]['version']==version for name,version in roots.items())
        selected=set();virtual=set();pending=list(roots)
        while pending:
            name=pending.pop()
            if name in selected:continue
            if name.startswith('__'):virtual.add(name);continue
            assert name in installed,f'Missing installed dependency {name}'
            selected.add(name)
            for spec in installed[name].get('depends',[]):
                match=re.match(r'^([A-Za-z0-9_.-]+)',spec);assert match
                pending.append(match.group(1))
        items=[installed[name] for name in sorted(selected)];urls=[]
        for item in items:
            url=item['url'];parts=urlparse(url)
            assert parts.scheme=='https' and parts.hostname in ('conda.anaconda.org','repo.anaconda.com')
            assert not parts.username and not parts.password and not parts.query and not parts.fragment
            assert item['fn']==Path(parts.path).name
            urls.append(url+'#'+item['md5'])
        explicit=args.out/'conda_explicit.txt';explicit.write_text('@EXPLICIT\n'+'\n'.join(urls)+'\n')
        (args.out/'source_packages.json').write_text(json.dumps(items,indent=2)+'\n')
        destination=args.out/'environment'
        command=[str(mamba),'create','--yes','--offline','--copy','--prefix',str(destination),'--file',str(explicit)]
        with (args.out/'create.stdout').open('w') as out,(args.out/'create.stderr').open('w') as error:
            subprocess.run(command,stdout=out,stderr=error,check=True)
        staged={json.loads(p.read_text())['name']:json.loads(p.read_text()) for p in (destination/'conda-meta').glob('*.json')}
        assert set(staged)==selected
        for item in items:
            for key in ('name','version','build','url','md5'):
                assert staged[item['name']][key]==item[key],f'Staged package mismatch: {item["name"]}/{key}'
        after={p.name:file_record(p)['sha256'] for p in sorted((SOURCE/'conda-meta').glob('*')) if p.is_file()}
        assert before==after,'Source conda environment changed'
        inventory={}
        for path in sorted(destination.rglob('*')):
            if path.is_symlink():assert path.resolve().is_relative_to(destination.resolve()),f'Escaping tool symlink: {path}'
            if path.is_file():inventory[str(path.relative_to(destination))]=file_record(path)
        (args.out/'tool_files.json').write_text(json.dumps(inventory,indent=2)+'\n')
        env={**os.environ,'LD_LIBRARY_PATH':str(destination/'lib')};versions={}
        for name,option,expected in [('STAR-plain','--version','2.7.11b'),('samtools','--version','1.22.1'),
                                     ('featureCounts','-v','2.1.1')]:
            executable=destination/'bin'/name
            result=subprocess.run([str(executable),option],env=env,capture_output=True,text=True,timeout=30,check=True)
            text=result.stdout+result.stderr;assert expected in text
            dependencies=subprocess.run(['/usr/bin/ldd',str(executable)],env=env,capture_output=True,text=True,timeout=30,check=True)
            assert 'not found' not in dependencies.stdout+dependencies.stderr
            versions[name]=dict(executable=file_record(executable),argv=[str(executable),option],stdout=result.stdout,
                                stderr=result.stderr,ldd=dependencies.stdout,ldd_stderr=dependencies.stderr)
        (args.out/'report.json').write_text(json.dumps(dict(environment=str(destination),roots=roots,
                exact_conda_packages=len(items),virtual_requirements=sorted(virtual),source_conda_meta=before,
                destination_package_identities_exact=True,source_conda_meta_unchanged=True,
                command=command,versions=versions,alignment_run=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
