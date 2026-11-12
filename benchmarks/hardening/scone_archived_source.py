"""Inspect only relevant code members from the author-deposited Zenodo archive."""
import argparse
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlparse
from urllib.request import Request,urlopen
import zipfile

from haplodose.provenance import start_run,finish_run,file_record

WANTED=('Barcodes/scONEseq_barcode.fasta','RawDataProcessing/Separation_bash.sh',
        'RawDataProcessing/1_process_reads.sh','RawDataProcessing/3_DNA_UMI_mapping.sh',
        'RawDataProcessing/4_RNA_UMI_quantification.sh')


def fetch(url,cap):
    assert urlparse(url).hostname=='zenodo.org'
    with urlopen(Request(url,headers={'User-Agent':'HaploDose-primary-source'}),timeout=45) as response:
        assert urlparse(response.url).hostname=='zenodo.org'
        content=response.read(cap+1)
    assert len(content)<=cap
    return content


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--task',type=Path,required=True);p.add_argument('--metadata-task',type=Path,required=True)
    a=p.parse_args();a.out.mkdir(exist_ok=False)
    run=start_run(a.out,dict(script=Path(__file__),design=a.task/'ARCHIVED_SOURCE_DESIGN.md',
                            current_barcode=a.task/'barcodes.json',current_driver=a.task/'driver.json',
                            current_sources=a.metadata_task/'source_records.json'),{},
                  dict(scope='Archived author code only; no external code execution, matrix or clone-label inspection'),mode='scone_archived_source')
    try:
        raw=fetch('https://zenodo.org/api/records/6796059',1024**2)
        (a.out/'zenodo_record.json').write_bytes(raw);record=json.loads(raw)
        files=[r for r in record['files'] if r['key']=='scONE-seq-data-processing-main.zip']
        assert len(files)==1;item=files[0];assert item['size']<=20*1024**2
        content=fetch(item['links']['self'],20*1024**2)
        assert len(content)==item['size']
        assert item['checksum'].startswith('md5:')
        assert hashlib.md5(content).hexdigest()==item['checksum'].split(':',1)[1]
        (a.out/'author_code.zip').write_bytes(content)
        current={r['path']:r['content'].encode() for r in json.loads((a.metadata_task/'source_records.json').read_text())['records']}
        current['Barcodes/scONEseq_barcode.fasta']=json.loads((a.task/'barcodes.json').read_text())['content'].encode()
        current['RawDataProcessing/1_process_reads.sh']=json.loads((a.task/'driver.json').read_text())['content'].encode()
        members=[]
        with zipfile.ZipFile(io.BytesIO(content)) as archive:
            names=archive.namelist();assert len(names)==len(set(names))
            (a.out/'archive_member_names.json').write_text(json.dumps(names,indent=2)+'\n')
            for wanted in WANTED:
                matching=[name for name in names if name==wanted or name.endswith('/'+wanted)]
                assert len(matching)==1
                info=archive.getinfo(matching[0]);assert info.file_size<=100000
                text=archive.read(info)
                destination=a.out/Path(wanted).name
                destination.write_bytes(text)
                members.append(dict(path=wanted,zip_member=matching[0],bytes=len(text),
                                    sha256=hashlib.sha256(text).hexdigest(),
                                    byte_identical_to_inspected_current_source=text==current[wanted]))
        (a.out/'report.json').write_text(json.dumps(dict(archive_metadata_checksum_verified=True,
                    record_id=6796059,archive_sha256=hashlib.sha256(content).hexdigest(),members=members,
                    executed_external_code=False,measurement_members_opened=False,exact_executed_version_verified=False,
                    limitation='A deposited code version is not proof of the commands used for a specific library'),indent=2)+'\n')
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
