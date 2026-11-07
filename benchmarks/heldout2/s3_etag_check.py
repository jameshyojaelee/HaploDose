"""Check a file against an S3 ETag. A plain ETag is the file md5; a multipart ETag "<hex>-<n>" is the md5 of the
concatenated per-part md5 digests. The part size is found by trying common sizes (1-64 MiB) that give n parts.
Usage: python -I s3_etag_check.py <file> <etag>   -> prints the md5 of the whole file and 'etag:match' or exits 1."""
import hashlib
import os
import sys

path, etag = sys.argv[1], sys.argv[2].strip('"')
size = os.path.getsize(path)
whole = hashlib.md5()
if "-" not in etag:
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            whole.update(b)
    ok = whole.hexdigest() == etag
else:
    target, n = etag.split("-")
    n = int(n)
    sizes = [m << 20 for m in range(1, 65) if -(-size // (m << 20)) == n]
    digests = {s: [] for s in sizes}
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):   # 1 MiB blocks; every candidate part size is a multiple
            whole.update(b)
            for s in sizes:
                d = digests[s]
                if not d or d[-1][1] == s:
                    d.append([hashlib.md5(), 0])
                d[-1][0].update(b)
                d[-1][1] += len(b)
    ok = any(hashlib.md5(b"".join(h.digest() for h, _ in digests[s])).hexdigest() == target for s in sizes)
print(whole.hexdigest(), "etag:match" if ok else "etag:MISMATCH")
sys.exit(0 if ok else 1)
