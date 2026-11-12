"""Bounded IDAT metadata-only probe of a prespecified array pair."""
import argparse
import gzip
import io
import json
from pathlib import Path
import struct
from urllib.request import Request, urlopen

from haplodose.provenance import start_run, finish_run, file_record


def metadata(raw):
    assert raw[:4] == b'IDAT' and len(raw) >= 16, 'Invalid IDAT signature'
    version, count = struct.unpack_from('<qi', raw, 4)
    assert version == 3 and 0 < count <= 256, 'Unsupported IDAT layout'
    header_end = 16 + 10*count
    assert header_end <= len(raw), 'Truncated directory'
    entries = [struct.unpack_from('<Hq',raw,16+10*i) for i in range(count)]
    assert len({k for k,v in entries}) == count, 'Duplicate directory IDs'
    assert all(header_end <= v < len(raw) for k,v in entries), 'Invalid directory offsets'
    assert len({v for k,v in entries}) == count, 'Overlapping field starts'
    offsets = dict(entries)
    output = dict(idat_version=version, fields_present=sorted(offsets))
    for identifier, name in [(200,'mid_block_count'),(401,'snp_manifest'),(402,'barcode'),
                             (403,'chip_type'),(404,'position'),(405,'beadset'),
                             (407,'description'),(510,'vln'),(1000,'probe_count')]:
        if identifier not in offsets:
            output[name] = None; continue
        start = offsets[identifier]
        end = min([v for v in offsets.values() if v > start] or [len(raw)])
        if identifier in (200,1000):
            assert start+4 <= end, 'Truncated count metadata'
            output[name] = struct.unpack_from('<i',raw,start)[0]
            assert output[name] >= 0
            continue
        length = 0; shift = 0; cursor = start
        while True:
            assert cursor < end and shift <= 28, 'Invalid metadata string prefix'
            byte = raw[cursor]; cursor += 1
            length |= (byte & 127) << shift
            if not byte & 128: break
            shift += 7
        assert length <= 4096 and cursor+length <= end, 'Invalid metadata string length'
        output[name] = raw[cursor:cursor+length].decode('utf-8')
    return output


def fixtures():
    # Independently specified bytes: directory at16, string at64, probe count
    # at80, and an uninterpreted intensity field at96.
    raw = bytearray(128)
    raw[:16] = b'IDAT\x03\x00\x00\x00\x00\x00\x00\x00\x03\x00\x00\x00'
    raw[16:46] = bytes.fromhex('91014000000000000000e803500000000000000068006000000000000000')
    raw[64:72] = b'\x07GSA-v2A'
    raw[80:84] = b'\xd2\x04\x00\x00'
    raw[96:] = bytes([255])*32
    parsed = metadata(bytes(raw))
    assert parsed['snp_manifest']=='GSA-v2A' and parsed['probe_count']==1234
    for offset, value in [(26,b'\x91\x01'),(18,b'\xff'*8)]:
        bad = bytearray(raw);bad[offset:offset+len(value)] = value
        try: metadata(bytes(bad))
        except AssertionError: pass
        else: raise AssertionError('Invalid directory fixture accepted')
    # Two-byte variable length encoding exercises strings beyond127 bytes.
    long = bytearray(400)
    long[:16] = b'IDAT\x03\x00\x00\x00\x00\x00\x00\x00\x01\x00\x00\x00'
    long[16:26] = bytes.fromhex('93014000000000000000')
    long[64:266] = b'\xc8\x01'+b'A'*200
    assert metadata(bytes(long))['chip_type']=='A'*200
    return 4


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--task',type=Path,required=True)
    p.add_argument('--preserved',action='store_true')
    a=p.parse_args();a.out.mkdir(mode=0o700,exist_ok=False)
    inputs=dict(script=Path(__file__),gate=a.task/'IDAT_HEADER_GATE.md',
                reference=a.task/'idat_format_reference.json',pilot=a.task/'metadata_v1/report.json')
    if a.preserved:
        inputs['extension_gate']=a.task/'IDAT_EXTENSION_GATE.md'
        inputs['original_manifest']=a.task/'idat_v1/run_manifest.json'
        prior=json.loads(inputs['original_manifest'].read_text());assert prior['status']=='complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256']==record['sha256']
        for path in (a.task/'idat_v1').glob('*.idat.gz'):inputs[path.name]=path
    run=start_run(a.out,inputs,dict(sample='GSM5901579',compressed_cap=8388608,decoded_cap=33554432),
                  dict(scope='Array resource metadata; intensity/probe arrays not decoded'),mode='onek1k_idat_metadata')
    try:
        checked=fixtures()
        sample=next(s for s in json.loads(inputs['pilot'].read_text())['samples'] if s['accession']=='GSM5901579')
        links=[v for k,values in sample['fields'].items() if k.startswith('supplementary_file') for v in values]
        assert len(links)==2 and all(s.endswith('.idat.gz') for s in links)
        records=[]
        for deposited in links:
            assert deposited.startswith('ftp://ftp.ncbi.nlm.nih.gov/geo/')
            url=deposited.replace('ftp://','https://',1)
            name=url.rsplit('/',1)[1]
            if a.preserved:
                compressed=inputs[name].read_bytes()
            else:
                with urlopen(Request(url,headers={'User-Agent':'HaploDose-resource-metadata'}),timeout=90) as response:
                    compressed=response.read(8388609)
            assert len(compressed)<=8388608, 'Compressed resource cap exceeded'
            if not a.preserved:
                path=a.out/name;path.write_bytes(compressed);path.chmod(0o600)
            with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:raw=stream.read(33554433)
            assert len(raw)<=33554432, 'Decoded resource cap exceeded'
            records.append(dict(deposited_url=deposited,requested_url=url,
                compressed_bytes=len(compressed),decoded_bytes=len(raw),metadata=metadata(raw)))
        (a.out/'report.json').write_text(json.dumps(dict(fixtures_passed=checked,files=records,
            preserved_inputs=a.preserved,
            intensity_arrays_decoded=False,genotypes_called=False,cohort_admitted=False,
            biological_null_specimens_verified=0),indent=2))
        finish_run(a.out,run,backend='cpu')
    except Exception as error:
        finish_run(a.out,run,error=error);raise


if __name__=='__main__':main()
