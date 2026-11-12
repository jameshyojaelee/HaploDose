"""Stage only six audited missing archives; preserve frozen inputs and outcomes."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.parse import urlparse
from urllib.request import Request,urlopen

from haplodose.provenance import start_run,finish_run,file_record


def verified(directory):
    manifest=json.loads((directory/'run_manifest.json').read_text());assert manifest['status']=='complete'
    for group in ('input','output'):
        for record in manifest[group]['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--mode',choices=('stage','audit'),required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False);root=a.task.parents[1]
    original=root/'data/external/dntrseq_gse144296';audit=a.task/'audit_v1';metadata=a.task/'metadata_v2'
    verified(audit);check=json.loads((audit/'report.json').read_text());assert check['all_checks_passed'] and check['missing_files']==6
    gaps=json.loads((audit/'missing_files.json').read_text());assert len(gaps)==6 and sum(x['expected_bytes'] for x in gaps)==830629012
    inputs=dict(script=Path(__file__),design=a.task/'REPAIR_DESIGN.md',audit=audit/'run_manifest.json',missing=audit/'missing_files.json',
                file_presence=metadata/'file_presence.json',frozen_manifest=original/'manifest.tsv')
    stage=a.task/'repaired_inputs_v1'
    if a.mode=='audit':
        verified(stage);inputs.update(stage_manifest=stage/'run_manifest.json',overlay=stage/'resolved_inputs.tsv',repair_report=stage/'report.json')
    rec=start_run(a.out,inputs,dict(mode=a.mode,only_audited_missing_files=True),
        dict(scope='Separate raw-archive completion,not biological validation or historical result update'),mode='dntr_missing_inputs')
    try:
        if a.mode=='stage':
            repaired=[]
            for gap in gaps:
                relative=Path(gap['file']);assert not relative.is_absolute() and '..' not in relative.parts
                assert not (original/relative).exists(),'Frozen missing file changed concurrently; do not overwrite it'
                url=urlparse(gap['url']);assert url.scheme=='https' and url.hostname=='ftp.sra.ebi.ac.uk' and url.path.startswith('/vol1/fastq/')
                target=a.out/'files'/relative;target.parent.mkdir(parents=True,exist_ok=True);partial=target.with_name(target.name+'.partial')
                md5=hashlib.md5();sha256=hashlib.sha256();size=0
                with urlopen(Request(gap['url'],headers={'User-Agent':'HaploDose-bounded-archive-repair'}),timeout=45) as response:
                    assert urlparse(response.url).hostname=='ftp.sra.ebi.ac.uk'
                    with partial.open('xb') as handle:
                        while True:
                            block=response.read(1024*1024)
                            if not block:break
                            size+=len(block);assert size<=gap['expected_bytes'],'Response exceeds the frozen archive size'
                            md5.update(block);sha256.update(block);handle.write(block)
                assert size==gap['expected_bytes'] and md5.hexdigest()==gap['expected_md5'],gap['file']
                assert not target.exists(),'Unexpected new-output collision; do not overwrite it'
                partial.rename(target)
                repaired.append(dict(run=gap['run'],logical_file=gap['file'],path=str(target),url=gap['url'],
                    bytes=size,md5=md5.hexdigest(),sha256=sha256.hexdigest(),sequence_records_decoded=False))
                (a.out/'repaired_files.json').write_text(json.dumps(repaired,indent=2)+'\n')
                print(json.dumps(dict(repaired=len(repaired),total=6,file=gap['file'],bytes=size)),flush=True)
            by_file={r['logical_file']:r for r in repaired}
            archived={r['file']:r for r in csv.DictReader((original/'manifest.tsv').open(),delimiter='\t')}
            files=json.loads((metadata/'file_presence.json').read_text());rows=[]
            for entry in files:
                name=entry['file']
                if name in by_file:
                    item=by_file[name];path=Path(item['path']);size=item['bytes'];md5=item['md5'];status='new_archive_md5_verified'
                else:
                    item=archived[name];path=original/name;size=int(item['bytes']);md5=item['md5'];status='frozen_archive_size_and_historical_md5_record_only'
                    assert entry['archived_metadata_matches'] and entry['size_matches']
                assert path.is_file() and path.stat().st_size==size
                rows.append(dict(run=entry['run'],logical_file=name,path=str(path),bytes=size,expected_md5=md5,
                    current_content_md5_verified=name in by_file,status=status))
            assert len(rows)==9738 and len({r['logical_file'] for r in rows})==len(rows)
            with (a.out/'resolved_inputs.tsv').open('w') as handle:
                writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
            report=dict(all_staged_files_verified=True,archives_repaired=6,bytes_repaired=sum(r['bytes'] for r in repaired),
                resolved_archive_paths=len(rows),frozen_paths_modified=0,frozen_manifests_modified=0,
                historical_results_changed=False,old_archives_current_content_checksums_recomputed=False,
                sequence_records_decoded=False,rna_or_dna_measurements_analyzed=False,
                biological_null_specimens_added=0,cohort_admitted_for_independent_evaluation=False)
        else:
            repaired=json.loads((stage/'repaired_files.json').read_text());assert len(repaired)==6
            result=subprocess.run(['/usr/bin/md5sum','--',*[r['path'] for r in repaired]],check=True,capture_output=True,text=True)
            (a.out/'independent_md5sum.txt').write_text(result.stdout)
            observed={line[34:]:line[:32] for line in result.stdout.splitlines()}
            expected={str(stage/'files'/r['file']):r for r in gaps};assert set(observed)==set(expected)
            for path,gap in expected.items():
                assert observed[path]==gap['expected_md5'] and Path(path).stat().st_size==gap['expected_bytes']
                assert not (original/gap['file']).exists(),'Original frozen input namespace must remain untouched'
            rows=list(csv.DictReader((stage/'resolved_inputs.tsv').open(),delimiter='\t'));assert len(rows)==9738
            repaired_rows=[r for r in rows if r['current_content_md5_verified']=='True'];assert len(repaired_rows)==6
            assert {r['path'] for r in repaired_rows}==set(expected)
            for row in rows:assert Path(row['path']).is_file() and Path(row['path']).stat().st_size==int(row['bytes'])
            report=dict(all_checks_passed=True,independent_md5sum_verified=6,resolved_archive_paths_checked=len(rows),
                frozen_missing_paths_still_missing=6,frozen_manifest_hash_unchanged=True,
                sequence_records_decoded=False,biological_null_specimens_added=0)
        for record in rec['input']['files'].values():assert file_record(record['path'])['sha256']==record['sha256']
        (a.out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        finish_run(a.out,rec,backend='CPU/network raw archive copy and checksum only' if a.mode=='stage' else 'Independent system md5sum and exact path/size audit')
    except BaseException as error:finish_run(a.out,rec,error=error);raise


if __name__=='__main__':main()
