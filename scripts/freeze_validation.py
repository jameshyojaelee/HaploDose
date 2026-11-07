"""Capture a coherent, bounded source tree for compute validation, including dirty files."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


def inventory(root):
    paths = [root / p for p in ['AGENTS.md', 'README.md', 'pyproject.toml', '.gitignore']]
    for name in ['python', 'crates', 'tests', 'benchmarks', 'scripts', 'slurm', 'docs']:
        for directory, dirs, files in os.walk(root / name):
            dirs[:] = [d for d in dirs if d not in {'target', '__pycache__', '.git', 'schmid_repo', '.pytest_cache'}]
            paths.extend(Path(directory) / f for f in files if Path(f).suffix not in {
                '.pyc', '.pdf', '.png', '.h5', '.h5ad', '.rds', '.RData', '.bam', '.npy', '.npz'})
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(set(paths)) if p.is_file()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    root, out = args.root.resolve(), args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    before = inventory(root)
    for relative in before:
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / relative, target)
    after = inventory(root)
    if before != after or any(hashlib.sha256((out / p).read_bytes()).hexdigest() != h for p, h in before.items()):
        raise RuntimeError('Source changed during snapshot; snapshot is incomplete, use a new destination')
    (out / 'source_manifest.json').write_text(json.dumps(before, sort_keys=True, indent=2) + '\n')
    for filename, command in [('git_status.txt', ['git', 'status', '--porcelain=v1']),
                              ('git_head.txt', ['git', 'rev-parse', 'HEAD'])]:
        (out / filename).write_bytes(subprocess.check_output(command, cwd=root))
    (out / 'SNAPSHOT_COMPLETE').write_text(hashlib.sha256(json.dumps(before, sort_keys=True).encode()).hexdigest() + '\n')
    print(json.dumps({'snapshot': str(out), 'files': len(before)}))


if __name__ == '__main__':
    main()
