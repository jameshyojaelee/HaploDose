"""Retrieve one primary-linked GSAMD candidate and assess address metadata only."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
from urllib.request import Request, urlopen
import zipfile

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(mode=0o700, exist_ok=False)
    helper = Path(__file__).with_name('onek1k_parent_addresses.py')
    inputs = dict(script=Path(__file__), helper=helper, design=args.task/'CUSTOM_RESOURCE_DESIGN.md',
                  parent=args.task/'parent_addresses_v1/run_manifest.json')
    for path in (args.task/'idat_v1').glob('*.idat.gz'): inputs[path.name] = path
    run = start_run(args.out, inputs, {}, dict(scope='Custom manifest candidate resource metadata only'),
                    mode='onek1k_custom_resource')
    try:
        prior = json.loads(inputs['parent'].read_text()); assert prior['status'] == 'complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        url = 'https://glimdna.org/assets/gsamd-24v2-0_20024620_a1.zip'
        target = args.out/'candidate.zip'
        with urlopen(Request(url, headers={'User-Agent':'HaploDose-resource-audit'}), timeout=180) as response:
            resolved = response.url; content_type = response.headers.get('Content-Type')
            with target.open('xb') as handle:
                total = 0
                while block := response.read(1024*1024):
                    total += len(block); assert total <= 128*1024*1024
                    handle.write(block)
        (args.out/'retrieval.json').write_text(json.dumps(dict(url=url, resolved_url=resolved,
            content_type=content_type, file=file_record(target)), indent=2)+'\n')
        with zipfile.ZipFile(target) as archive:
            listing = archive.infolist()
            assert sum(item.file_size for item in listing) <= 512*1024*1024
            (args.out/'archive_members.json').write_text(json.dumps([dict(name=i.filename,bytes=i.file_size) for i in listing], indent=2)+'\n')
            chosen = [i for i in listing if Path(i.filename).name.lower() == 'gsamd-24v2-0_20024620_a1.csv']
            assert len(chosen) == 1, 'Expected exact candidate CSV name absent or ambiguous'
            path = args.out/Path(chosen[0].filename).name
            path.write_bytes(archive.read(chosen[0]))
        spec = importlib.util.spec_from_file_location('address_reference', helper)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
        assert module.fixtures() == 2
        channels = [module.addresses(gzip.decompress(p.read_bytes())) for p in sorted((args.task/'idat_v1').glob('*.idat.gz'))]
        assert len(channels) == 2 and channels[0] == channels[1]
        observed = channels[0]; names = set(); required = set(); missing = []; total = 0
        with path.open(encoding='utf-8-sig', newline='') as handle:
            records = csv.reader(handle)
            for row in records:
                if row == ['[Assay]']: break
            else: raise AssertionError('Missing assay section')
            keys = next(records)
            with (args.out/'probe_support.tsv').open('x') as output:
                writer = csv.writer(output, delimiter='\t')
                columns = ['IlmnID','Name','AddressA_ID','AddressB_ID','SNP','Chr','MapInfo','GenomeBuild','IlmnStrand','RefStrand']
                assert set(columns) <= set(keys)
                writer.writerow(columns+['addresses_present'])
                for row in records:
                    if row and row[0].startswith('['): break
                    if not row: continue
                    item = dict(zip(keys, row)); assert item['Name'] not in names
                    names.add(item['Name'])
                    values = {int(item[key]) for key in ('AddressA_ID','AddressB_ID') if item[key] not in ('','0')}
                    assert values and all(v > 0 for v in values); required |= values
                    absent = values-observed
                    if absent: missing.append(dict(name=item['Name'],addresses=sorted(absent)))
                    writer.writerow([item[key] for key in columns]+[int(not absent)]); total += 1
        assert total
        (args.out/'report.json').write_text(json.dumps(dict(candidate=file_record(path),
            assay_probes=total, required_addresses=len(required), required_missing=len(required-observed),
            idat_addresses_outside_candidate=len(observed-required), missing_probes=missing,
            full_address_compatibility=not missing, independent_audit='Not run',
            allele_probe_sequence_identity='Not established by address overlap',
            corresponding_BPM_EGT='Not established', intensity_arrays_decoded=False,
            cohort_admitted=False, biological_null_specimens_verified=0), indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__': main()
