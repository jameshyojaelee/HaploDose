"""Independently verify every parent-address support row against preserved IDATs."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np

from haplodose.provenance import start_run, finish_run, file_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    args = parser.parse_args(); args.out.mkdir(mode=0o700, exist_ok=False)
    inputs = dict(script=Path(__file__), design=args.task/'PARENT_ADDRESS_DESIGN.md',
                  manifest=args.task/'parent_addresses_v1/run_manifest.json')
    for path in (args.task/'idat_v1').glob('*.idat.gz'):
        inputs[path.name] = path
    run = start_run(args.out, inputs, {}, dict(scope='Independent address evidence; no intensity access'),
                    mode='onek1k_parent_addresses_audit')
    try:
        prior = json.loads(inputs['manifest'].read_text()); assert prior['status'] == 'complete'
        for record in prior['output']['files'].values():
            assert file_record(record['path'])['sha256'] == record['sha256']
        report = json.loads((args.task/'parent_addresses_v1/report.json').read_text())
        assert np.frombuffer(bytes.fromhex('010000000001000000000100'), dtype='<u4').tolist() == [1,256,65536]
        channels = []
        for path in sorted((args.task/'idat_v1').glob('*.idat.gz')):
            raw = gzip.decompress(path.read_bytes())
            assert raw[:4] == b'IDAT' and int.from_bytes(raw[4:12], 'little') == 3
            count = int.from_bytes(raw[12:16], 'little')
            fields = {}
            for offset in range(16, 16+10*count, 10):
                key = int.from_bytes(raw[offset:offset+2], 'little')
                value = int.from_bytes(raw[offset+2:offset+10], 'little')
                assert key not in fields; fields[key] = value
            n = int.from_bytes(raw[fields[1000]:fields[1000]+4], 'little')
            start = fields[102]
            stop = min([v for v in fields.values() if v > start] or [len(raw)])
            assert stop-start == 4*n
            identifiers = np.frombuffer(raw, dtype='<u4', count=n, offset=start)
            assert len(np.unique(identifiers)) == n == report['idat_address_counts'][path.name]
            channels.append(set(identifiers.tolist()))
        assert len(channels) == 2 and channels[0] == channels[1]
        observed = channels[0]
        csv_path = Path(report['archives']['csv']['member']['path'])
        with csv_path.open(encoding='utf-8-sig', newline='') as handle:
            records = csv.reader(handle)
            for row in records:
                if row == ['[Assay]']: break
            else: raise AssertionError('Missing assay section')
            keys = next(records)
            with (args.task/'parent_addresses_v1/parent_probe_support.tsv').open(newline='') as support:
                exported = csv.DictReader(support, delimiter='\t')
                required = set(); missing = []; names = set(); total = 0
                for row in records:
                    if row and row[0].startswith('['): break
                    if not row: continue
                    source = dict(zip(keys, row)); target = next(exported)
                    assert all(source[key] == target[key] for key in target if key != 'addresses_present')
                    assert source['Name'] not in names; names.add(source['Name'])
                    addresses = {int(source[key]) for key in ('AddressA_ID', 'AddressB_ID') if source[key] not in ('', '0')}
                    assert addresses; required |= addresses
                    absent = addresses-observed
                    assert int(target['addresses_present']) == int(not absent)
                    if absent: missing.append(dict(name=source['Name'], addresses=sorted(absent)))
                    total += 1
                assert next(exported, None) is None
        assert total == report['parent_rows'] == report['unique_parent_names']
        assert len(required) == report['parent_required_addresses']
        assert len(required-observed) == report['required_missing']
        assert len(required & observed) == report['required_present']
        assert len(observed-required) == report['idat_addresses_outside_parent']
        assert len(missing) == report['probe_support'].get('missing_address', 0)
        assert total-len(missing) == report['probe_support']['complete']
        result = dict(rows_exact=total, missing_parent_probes=missing,
            full_parent_address_compatibility=not missing,
            intensity_arrays_decoded=False, genotypes_called=False, cohort_admitted=False,
            biological_null_specimens_verified=0)
        (args.out/'report.json').write_text(json.dumps(result, indent=2)+'\n')
        finish_run(args.out, run, backend='cpu')
    except Exception as error:
        finish_run(args.out, run, error=error)
        raise


if __name__ == '__main__':
    main()
