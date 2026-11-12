"""Private, checksum-identified BWA index for the verified DNA reference."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    bwa=Path(os.environ['HAPLODOSE_BWA_BIN']).resolve()
    wrapper=Path(__file__).resolve().parents[2]/'scripts/scone_bwa_readonly.sbatch'
    inputs=dict(script=Path(__file__),design=args.task/'DNA_INDEX_DESIGN.md',wrapper=wrapper,bwa=bwa,
                modules=Path(os.environ['HAPLODOSE_BWA_MODULE_RECORD']),
                reference=args.task/'reference_v3/run_manifest.json',preparation=args.task/'prepare_audit_v1/run_manifest.json',
                tools=args.task/'tool_environment_v1/run_manifest.json')
    run=start_run(args.out,inputs,dict(index_algorithm='bwtsw',module='BWA/0.7.18-GCCcore-12.3.0'),
                  dict(scope='DNA reference preparation only, no alignment or truth'),mode='scone_dna_index')
    try:
        for name in ('reference','preparation','tools'):
            manifest=json.loads(inputs[name].read_text());assert manifest['status']=='complete'
            for item in manifest['output']['files'].values():assert file_record(item['path'])['sha256']==item['sha256']
        reference=json.loads((args.task/'reference_v3/report.json').read_text())
        item=reference['reference']['Homo_sapiens.GRCh38.dna.primary_assembly.fa'];source=Path(item['path'])
        private=args.out/'GRCh38_normalized.fa'
        with (source.parent.parent/'.setup_reference.lock').open('rb') as lock:
            fcntl.flock(lock,fcntl.LOCK_SH|fcntl.LOCK_NB)
            assert not private.exists();shutil.copyfile(source,private)
            copied=file_record(private);assert copied['sha256']==item['sha256'] and copied['bytes']==item['bytes']
        env={**os.environ,'LD_LIBRARY_PATH':os.environ['HAPLODOSE_BWA_LD_LIBRARY_PATH']}
        version=subprocess.run([str(bwa)],env=env,capture_output=True,text=True,timeout=30)
        assert version.returncode==1 and re.search(r'Version: 0\.7\.18(?:-|\s)',version.stdout+version.stderr)
        ldd=subprocess.run(['/usr/bin/ldd',str(bwa)],env=env,capture_output=True,text=True,timeout=30,check=True)
        assert 'not found' not in ldd.stdout+ldd.stderr
        libraries={Path(match).resolve() for match in re.findall(r'(/[^\s()]+)',ldd.stdout)}
        runtime=dict(executable=file_record(bwa),version_stdout=version.stdout,version_stderr=version.stderr,
            LD_LIBRARY_PATH=env['LD_LIBRARY_PATH'],ldd=ldd.stdout,libraries=[file_record(path) for path in sorted(libraries)])
        (args.out/'native_runtime.json').write_text(json.dumps(runtime,indent=2)+'\n')
        command=[str(bwa),'index','-a','bwtsw','-p',str(args.out/'GRCh38'),str(private)]
        started=time.monotonic()
        with (args.out/'bwa.stdout').open('w') as out,(args.out/'bwa.stderr').open('w') as error:
            subprocess.run(command,env=env,stdout=out,stderr=error,check=True)
        elapsed=time.monotonic()-started;artifacts={}
        for suffix in ('amb','ann','bwt','pac','sa'):
            path=args.out/('GRCh38.'+suffix);assert path.is_file() and path.stat().st_size>0
            artifacts[suffix]=file_record(path)
        tools=json.loads((args.task/'tool_environment_v1/report.json').read_text());prefix=Path(tools['environment'])
        samtools=prefix/'bin/samtools';expected=tools['versions']['samtools']['executable']
        assert file_record(samtools)['sha256']==expected['sha256']
        faidx=[str(samtools),'faidx',str(private)]
        with (args.out/'faidx.stdout').open('w') as out,(args.out/'faidx.stderr').open('w') as error:
            subprocess.run(faidx,env={**os.environ,'LD_LIBRARY_PATH':str(prefix/'lib')},stdout=out,stderr=error,check=True)
        lengths={line.split('\t')[0]:int(line.split('\t')[1]) for line in Path(str(private)+'.fai').read_text().splitlines()}
        assert lengths=={name:value['length'] for name,value in reference['contigs'].items()}
        assert file_record(private)['sha256']==copied['sha256']
        assert file_record(bwa)['sha256']==runtime['executable']['sha256']
        assert all(file_record(record['path'])['sha256']==record['sha256'] for record in runtime['libraries'])
        (args.out/'report.json').write_text(json.dumps(dict(private_reference=copied,index_prefix=str(args.out/'GRCh38'),
            native_runtime=runtime,index_files=artifacts,commands=[command,faidx],index_seconds=elapsed,
            FAI_contigs_exact=len(lengths),alignment_run=False,biological_null_specimens_verified=0),indent=2)+'\n')
        finish_run(args.out,run,backend='cpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':main()
