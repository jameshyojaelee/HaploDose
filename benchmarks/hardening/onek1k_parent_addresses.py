"""Metadata-only parent GSA resource and locked IDAT address compatibility check."""
import argparse
from collections import Counter
import csv
import gzip
import io
import json
from pathlib import Path
import struct
from urllib.request import Request, urlopen
import zipfile

from haplodose.provenance import start_run, finish_run, file_record


def addresses(raw):
    assert raw[:4] == b'IDAT'
    version, count = struct.unpack_from('<qi', raw, 4)
    assert version == 3 and 0 < count <= 256
    entries = [struct.unpack_from('<Hq', raw, 16+10*i) for i in range(count)]
    offsets = dict(entries)
    assert len(offsets) == count and len(set(offsets.values())) == count
    assert all(16+10*count <= v < len(raw) for v in offsets.values())
    def field(identifier):
        start = offsets[identifier]
        end = min([v for v in offsets.values() if v > start] or [len(raw)])
        return raw[start:end]
    number, = struct.unpack('<i', field(1000))
    values = field(102)
    assert number > 0 and len(values) == number*4
    result = struct.unpack('<'+'I'*number, values)
    assert len(set(result)) == number
    return set(result)


def fixtures():
    # Three independently encoded addresses; intensity bytes are uninterpreted.
    raw = bytearray(120)
    raw[:16] = b'IDAT'+struct.pack('<qi', 3, 3)
    raw[16:46] = struct.pack('<HqHqHq', 1000, 60, 102, 64, 104, 76)
    raw[60:76] = bytes.fromhex('03000000010000000001000000000100')
    raw[76:] = b'\xff'*44
    assert addresses(raw) == {1, 256, 65536}
    bad = bytearray(raw); bad[72:76] = bad[64:68]
    try:
        addresses(bad)
    except AssertionError:
        pass
    else:
        raise AssertionError('Duplicate IDAT probe address accepted')
    return 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(mode=0o700, exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'PARENT_ADDRESS_DESIGN.md',
                  pilot_manifest=args.task/'idat_v1/run_manifest.json')
    for name in ('step_1_idat_to_gtc_sh', 'step_2_gtc_to_vcf_sh', 'README_md'):
        inputs[name] = args.task/('loy_author_'+name+'.json')
    for path in (args.task/'idat_v1').glob('*.idat.gz'):
        inputs[path.name] = path
    run = start_run(args.out, inputs, dict(sample='GSM5901579'),
        dict(scope='Resource/address metadata only; not intensity, genotype or CNA validation'),
        mode='onek1k_parent_addresses')
    try:
        checked = fixtures()
        prior = json.loads(inputs['pilot_manifest'].read_text())
        assert prior['status'] == 'complete'
        for item in prior['output']['files'].values():
            assert file_record(item['path'])['sha256'] == item['sha256']
        archives = {}; members = {}
        base = 'https://webdata.illumina.com/downloads/productfiles/global-screening-array/v2-0/'
        for kind, suffix in [('csv', 'manifest-file-csv'), ('bpm', 'manifest-file-bpm'), ('egt', 'cluster-file')]:
            url = base+'gsa-24-v2-0-A1-'+suffix+'.zip'
            target = args.out/(kind+'.zip')
            with urlopen(Request(url, headers={'User-Agent': 'HaploDose-resource-audit'}), timeout=180) as response:
                with target.open('xb') as handle:
                    size = 0
                    while block := response.read(1024*1024):
                        size += len(block)
                        assert size <= 128*1024*1024, 'Archive cap exceeded'
                        handle.write(block)
                resolved = response.url
            with zipfile.ZipFile(target) as archive:
                listing = archive.infolist()
                assert sum(item.file_size for item in listing) <= 512*1024*1024
                chosen = [item for item in listing if item.filename.lower().endswith('.'+kind)]
                assert len(chosen) == 1
                member = chosen[0]
                path = args.out/Path(member.filename).name
                path.write_bytes(archive.read(member))
            members[kind] = path
            archives[kind] = dict(url=url, resolved_url=resolved, archive=file_record(target),
                member=file_record(path), archive_member=member.filename)
        channels = {}
        for path in sorted((args.task/'idat_v1').glob('*.idat.gz')):
            with gzip.open(path, 'rb') as handle:
                raw = handle.read(33554433)
            assert len(raw) <= 33554432
            channels[path.name] = addresses(raw)
        assert len(channels) == 2
        sets = list(channels.values())
        assert sets[0] == sets[1], 'Green/red address correspondence differs'
        observed = sets[0]
        text = members['csv'].read_text(encoding='utf-8-sig')
        lines = text.splitlines()
        start = lines.index('[Assay]')+1
        end = next((i for i in range(start, len(lines)) if lines[i].startswith('[')), len(lines))
        reader = csv.DictReader(io.StringIO('\n'.join(lines[start:end])))
        required = {'IlmnID', 'Name', 'AddressA_ID', 'AddressB_ID', 'SNP', 'Chr', 'MapInfo', 'GenomeBuild', 'IlmnStrand', 'RefStrand'}
        assert required.issubset(reader.fieldnames)
        rows = list(reader); assert rows
        names = [row['Name'] for row in rows]
        requested = set(); support = Counter(); builds = Counter(); chromosomes = Counter()
        output = args.out/'parent_probe_support.tsv'
        with output.open('x') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['IlmnID', 'Name', 'AddressA_ID', 'AddressB_ID', 'SNP', 'Chr', 'MapInfo', 'GenomeBuild', 'IlmnStrand', 'RefStrand', 'addresses_present'])
            for row in rows:
                values = [int(row[key]) for key in ('AddressA_ID', 'AddressB_ID') if row[key].strip() not in ('', '0')]
                assert values and all(value > 0 for value in values)
                requested.update(values)
                present = all(value in observed for value in values)
                support['complete' if present else 'missing_address'] += 1
                builds[row['GenomeBuild']] += 1; chromosomes[row['Chr']] += 1
                writer.writerow([row[key] for key in ['IlmnID', 'Name', 'AddressA_ID', 'AddressB_ID', 'SNP', 'Chr', 'MapInfo', 'GenomeBuild', 'IlmnStrand', 'RefStrand']]+[int(present)])
        report = dict(fixtures_passed=checked, sample='GSM5901579', archives=archives,
            idat_address_counts={name:len(value) for name,value in channels.items()},
            green_red_addresses_exact=True, parent_rows=len(rows), unique_parent_names=len(set(names)),
            parent_required_addresses=len(requested), required_present=len(requested & observed),
            required_missing=len(requested-observed), idat_addresses_outside_parent=len(observed-requested),
            probe_support=dict(support), genome_builds=dict(builds), chromosomes=dict(chromosomes),
            allele_probe_sequence_correspondence='Not established by address overlap alone',
            bpm_egt_binary_correspondence='Not yet independently parsed',
            intensity_arrays_decoded=False, genotypes_called=False, cohort_admitted=False,
            biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
