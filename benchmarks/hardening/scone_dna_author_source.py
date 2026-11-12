"""Extract bounded DNA-processing source from the already verified author archive."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

from haplodose.provenance import start_run, finish_run, file_record


WANTED = ('DNA_related/README.md', 'DNA_related/process_segnorm.R',
          'DNA_related/Calculate_integerCNV.R', 'DNA_related/mergeLevels_multi_segmentation.R',
          'DNA_related/Transfer_CHISEL_mbaf_to_ginkgo.R', 'DNA_related/extract_baf_from_chisel.R',
          'RawDataProcessing/2_DNA_BWA_mapping.sh', 'Required_packages.txt')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(exist_ok=False)
    parent = args.task.parent/'scone_raw_chain_20261109/archived_source_v1'
    run = start_run(args.out, dict(script=Path(__file__), design=args.task/'DNA_AUTHOR_SOURCE_DESIGN.md',
        prior=parent/'run_manifest.json'), {}, dict(scope='Author DNA source only; no code execution or measurements'),
        mode='scone_dna_author_source')
    try:
        prior = json.loads((parent/'run_manifest.json').read_text()); assert prior['status'] == 'complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        raw = (parent/'author_code.zip').read_bytes()
        report = json.loads((parent/'report.json').read_text())
        assert hashlib.sha256(raw).hexdigest() == report['archive_sha256']
        records = []
        with zipfile.ZipFile(parent/'author_code.zip') as archive:
            names = archive.namelist(); assert len(names) == len(set(names))
            for wanted in WANTED:
                matches = [name for name in names if name == wanted or name.endswith('/'+wanted)]
                assert len(matches) == 1
                member = archive.getinfo(matches[0]); assert member.file_size <= 1024**2
                content = archive.read(matches[0]); content.decode('utf-8')
                target = args.out/'source'/wanted; target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(content)
                assert target.read_bytes() == archive.read(matches[0])
                records.append(dict(path=wanted, zip_member=matches[0], bytes=len(content),
                    lines=len(content.splitlines()), sha256=hashlib.sha256(content).hexdigest()))
        (args.out/'report.json').write_text(json.dumps(dict(archive_sha256=report['archive_sha256'],
            selected_members=records, external_code_executed=False, assay_values_parsed=False,
            historical_execution_verified=False, DNA_CNA_truth_constructed=False), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
