"""Execute the unchanged GPU-only wheel regression with an explicit no-skip gate."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

from haplodose.device import select
from haplodose.io import write_json
from haplodose.provenance import start_run,finish_run,file_record


def descendants(parent):
    result=subprocess.run(['ps','-e','-o','pid=,ppid='],capture_output=True,text=True,check=True)
    relations=[tuple(map(int,line.split())) for line in result.stdout.splitlines()]
    found={parent}
    while True:
        expanded=found|{pid for pid,ppid in relations if ppid in found}
        if expanded==found:
            return found
        found=expanded


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    test=Path(__file__).resolve().parents[2]/'tests/test_device.py'
    run=start_run(args.out,dict(script=Path(__file__),unchanged_test=test,
                  installed_artifact=Path(os.environ['HAPLODOSE_VALIDATED_ARTIFACT'])),{},
                  dict(scope='Existing simulated CPU/GPU whole-workflow regression; no biological validation'),mode='wheel_gpu_parity')
    try:
        assert select('gpu')=='gpu'
        import jax
        assert jax.default_backend()=='gpu' and all(device.platform=='gpu' for device in jax.devices())
        write_json(dict(backend=jax.default_backend(),devices=[str(device) for device in jax.devices()]),args.out/'backend.json')
        command=[sys.executable,'-I','-m','pytest','-q','-rA',str(test)+'::test_gpu_matches_cpu',
                 '--import-mode=importlib','-p','no:cacheprovider',
                 '--basetemp',str(args.out/'pytest_tmp'),'--junitxml',str(args.out/'junit.xml')]
        observations=[]
        with (args.out/'pytest.txt').open('w') as stream:
            process=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
            while process.poll() is None:
                owned=descendants(process.pid)
                result=subprocess.run(['nvidia-smi','--query-compute-apps=pid,used_gpu_memory',
                                       '--format=csv,noheader,nounits'],capture_output=True,text=True,check=True)
                total=0
                for row in csv.reader(result.stdout.splitlines()):
                    if len(row)==2 and int(row[0]) in owned:
                        total+=int(row[1])
                observations.append(dict(monotonic_seconds=time.monotonic(),owned_compute_memory_MiB=total))
                time.sleep(5)
            write_json(observations,args.out/'gpu_memory_samples.json')
            if process.returncode:
                raise subprocess.CalledProcessError(process.returncode,command)
        suite=ET.parse(args.out/'junit.xml').getroot()
        cases=list(suite.iter('testcase'))
        assert len(cases)==1
        assert not any(list(case.iter(tag)) for case in cases for tag in ('failure','error','skipped'))
        backends={}
        for name in ('cpu','gpu'):
            matches=list((args.out/'pytest_tmp').glob(f'**/{name}/qc.json'));assert len(matches)==1
            qc=json.loads(matches[0].read_text())
            assert qc['device']==qc['backend']==name
            backends[name]=dict(path=str(matches[0]),device=qc['device'],backend=qc['backend'])
        write_json(dict(existing_test_passed=True,tests_executed=1,skipped=0,backends_inside_inference=backends,
                        maximum_observed_owned_compute_memory_MiB=max((row['owned_compute_memory_MiB'] for row in observations),default=None),
                        memory_sampling_seconds=5,peak_between_samples_unobserved=True,
                        biological_specimens=0,scope='Default simulated CPU/GPU workflow only; no allele/spatial/bulk guarantee',
                        scientific_thresholds='Unchanged tests/test_device.py::test_gpu_matches_cpu'),args.out/'report.json')
        write_json({str(path.relative_to(args.out)):file_record(path) for path in sorted(args.out.rglob('*'))
                    if path.is_file() and path not in (args.out/'run_manifest.json',args.out/'output_files.json')},args.out/'output_files.json')
        finish_run(args.out,run,backend='gpu')
    except Exception as error:
        finish_run(args.out,run,error=error);raise


if __name__=='__main__':
    main()
