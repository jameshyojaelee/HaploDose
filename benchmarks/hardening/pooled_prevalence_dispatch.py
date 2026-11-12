"""Submit the frozen error/power pilot only after actual full recovery passes."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import subprocess

from haplodose.io import write_json
from haplodose.provenance import start_run, finish_run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--task',type=Path,required=True)
    args=parser.parse_args();args.out.mkdir(exist_ok=False)
    task=args.task.resolve();root=task.parents[1]
    workflow=task.parent/'pooled_prevalence_workflow_20261110'
    inputs=dict(script=Path(__file__),design=task/'PILOT_DESIGN.md',
        recovery=workflow/'candidate_v2/report.json',
        recovery_manifest=workflow/'candidate_v2/run_manifest.json',
        control=workflow/'strict_control_v1/report.json',
        numerics=task/'numerics_v1/report.json',
        pilot_source=task/'source_v1/source_manifest.json',
        method_source=workflow/'source_v2/source_manifest.json')
    rec=start_run(args.out,inputs,dict(workers=6,paired_workflows=60),
        dict(scope='Execution gate only; no new inference or generation'),mode='pooled_prevalence_dispatch')
    submitted=[]
    try:
        report=json.loads(inputs['recovery'].read_text())
        assert json.loads(inputs['recovery_manifest'].read_text())['status']=='complete'
        assert report['recovered']==5 and report['full_workflow_including_merge']
        assert len(report['cases'])==5
        assert all(r['original_assertion_passed']and all(c['passed']for c in r['criteria'].values())
                   for r in report['cases'])
        assert json.loads(inputs['control'].read_text())['strict_control']['exact']
        assert json.loads(inputs['numerics'].read_text())['exact_rational_confidence_inversion']
        for key in ('pilot_source','method_source'):
            for relative,digest in json.loads(inputs[key].read_text()).items():
                assert hashlib.sha256((inputs[key].parent/relative).read_bytes()).hexdigest()==digest
        write_json(dict(recovery=5,all_original_criteria=True,merge_included=True,
            exact_control=True,confidence_numerics=True),args.out/'GATES_PASSED.json')
        environment=root/'dev/gene_laplace_profile_20261110/environment'
        runtime=root/'dev/recovery_platform_20261109/runtime003'
        artifact=root/'dev/gene_laplace_profile_20261110/validation_23156464/installed_artifact.json'
        launcher=root/'dev/adjacent_reference_noise_20261110/analysis_runtime.sbatch'
        # A second dispatcher must never duplicate a partially submitted batch.
        with (task/'dispatch.lock').open('a')as lock:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            journal=json.loads((task/'jobs.json').read_text())
            assert not any(r['name'].startswith('pilot_')for r in journal)
            workers=[]
            for index,(condition,start)in enumerate((c,s)for c in('matched','state','two_clone')for s in(0,5)):
                account=('innovation','nslab')[index%2]
                name=f'pilot_{condition}_{start:02d}_v1'
                command=['sbatch','--parsable','--partition=cpu','--constraint=v6&amd',
                    '--nodelist=ne1dc6-005',f'--account={account}',f'--qos={account}',
                    '--cpus-per-task=4','--mem=32G','--time=90:00:00',
                    '--job-name=hd_pool_'+condition+'_'+str(start),
                    f'--output={task}/{name}_%j.log',
                    f'--export=ALL,HAPLODOSE_GLIBC_RUNTIME={runtime},HAPLODOSE_VALIDATED_ARTIFACT={artifact}',
                    str(launcher),str(task),str(environment),str(task/'source_v1'),
                    'benchmarks/hardening/pooled_prevalence_pilot.py',str(task/name),
                    '--task',str(task),'--condition',condition,'--start',str(start),'--stop',str(start+5)]
                job=subprocess.check_output(command,text=True).strip().split(';')[0]
                assert job.isdigit()
                entry=dict(name=name,job_id=job,command=command)
                journal.append(entry);submitted.append(entry);workers.append(job)
                write_json(journal,task/'jobs.json')
                write_json(submitted,args.out/'submitted.json')
                print(json.dumps(entry),flush=True)
            command=['sbatch','--parsable','--partition=cpu','--constraint=v6&amd',
                '--nodelist=ne1dc6-005','--account=nslab','--qos=nslab',
                '--cpus-per-task=1','--mem=4G','--time=01:00:00',
                '--job-name=hd_pool_pilot_summary',f'--output={task}/summary_%j.log',
                '--dependency=afterok:'+':'.join(workers),
                f'--export=ALL,HAPLODOSE_GLIBC_RUNTIME={runtime},HAPLODOSE_VALIDATED_ARTIFACT={artifact}',
                str(launcher),str(task),str(environment),str(task/'source_v1'),
                'benchmarks/hardening/pooled_prevalence_pilot_summary.py',str(task/'summary_v1'),
                '--task',str(task)]
            job=subprocess.check_output(command,text=True).strip().split(';')[0]
            assert job.isdigit()
            entry=dict(name='pilot_summary',job_id=job,command=command)
            journal.append(entry);submitted.append(entry)
            write_json(journal,task/'jobs.json');write_json(submitted,args.out/'submitted.json')
            print(json.dumps(entry),flush=True)
        finish_run(args.out,rec,backend='CPU report-verified Slurm dispatch; no counts read')
    except BaseException as error:
        write_json(submitted,args.out/'submitted.json');finish_run(args.out,rec,error=error);raise


if __name__=='__main__':main()
