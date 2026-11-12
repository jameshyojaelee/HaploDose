"""Bounded integrity diagnosis of the first failed DNTR archive transfer."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from haplodose.provenance import file_record, finish_run, start_run


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--task', type=Path, required=True)
    parser.add_argument('--all-headers', action='store_true')
    args = parser.parse_args()
    args.out.mkdir(exist_ok=False)
    missing = args.task / 'audit_v1/missing_files.json'
    gap = json.loads(missing.read_text())[0]
    assert gap['run'] == 'SRR12860636'
    partial = args.task / 'repaired_inputs_v1/files' / (gap['file'] + '.partial')
    assert partial.stat().st_size == 802
    rec = start_run(args.out, dict(script=Path(__file__), missing=missing,
        partial=partial, design=args.task/('TRANSFER_DIAGNOSTIC_FOLLOWUP.md'
            if args.all_headers else 'TRANSFER_DIAGNOSTIC.md')),
        dict(max_response_bytes=4096, methods=['HEAD', 'GET'], all_headers=args.all_headers),
        dict(scope='Archive transport diagnosis only; no sequence decoding'),
        mode='dntr_transfer_diagnostic')
    try:
        body = partial.read_bytes()
        report = dict(expected_bytes=gap['expected_bytes'],
            expected_md5=gap['expected_md5'], partial_bytes=len(body),
            partial_md5=hashlib.md5(body).hexdigest(),
            partial_gzip_magic=body[:2] == b'\x1f\x8b',
            partial_text=body.decode('utf-8', errors='replace') if body.lstrip().startswith(b'<') else None,
            requests=[], sequence_records_decoded=False)
        for method in ('HEAD', 'GET'):
            headers = {'User-Agent': 'HaploDose-bounded-archive-repair'}
            if method == 'GET':
                headers['Range'] = 'bytes=0-4095'
            request = Request(gap['url'], headers=headers, method=method)
            try:
                response = urlopen(request, timeout=45)
            except HTTPError as error:
                response = error
            with response:
                data = response.read(4096) if method == 'GET' else b''
                report['requests'].append(dict(method=method, url=response.url,
                    status=response.status, headers={k: v for k, v in response.headers.items()
                        if k.lower() in {'content-type', 'content-length', 'content-range',
                            'location', 'server', 'last-modified', 'etag'}},
                    bytes_inspected=len(data), gzip_magic=data[:2] == b'\x1f\x8b',
                    text=data.decode('utf-8', errors='replace') if data.lstrip().startswith(b'<') else None))
        if args.all_headers:
            report['remaining_mate_headers'] = []
            for other in json.loads(missing.read_text())[1:]:
                with urlopen(Request(other['url'], method='HEAD'), timeout=45) as response:
                    report['remaining_mate_headers'].append(dict(run=other['run'],
                        expected_bytes=other['expected_bytes'], url=response.url,
                        status=response.status, headers={k: v for k, v in response.headers.items()
                            if k.lower() in {'content-type', 'content-length', 'location'}}))
            command = ['/usr/bin/curl', '--silent', '--show-error', '--head',
                '--location', '--max-time', '45', '--max-redirs', '3',
                '--proto', '=https', gap['url']]
            result = subprocess.run(command, capture_output=True, text=True, check=True)
            assert len(result.stdout) < 65536
            report['independent_curl'] = dict(command=command, returncode=result.returncode,
                headers=result.stdout, stderr=result.stderr)
        assert file_record(partial)['sha256'] == rec['input']['files']['partial']['sha256']
        (args.out/'report.json').write_text(json.dumps(report, indent=2)+'\n')
        finish_run(args.out, rec, backend='CPU bounded HTTP headers and archive-byte integrity')
    except BaseException as error:
        finish_run(args.out, rec, error=error)
        raise


if __name__ == '__main__':
    main()
