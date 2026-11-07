"""Verify preserved evidence and record the source delta from the initial dirty tree."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--task-root', type=Path, required=True)
    ap.add_argument('--data-root', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parents[2]
    baseline = json.loads((a.task_root / 'baseline_manifest.json').read_text())
    current = json.loads((source / 'source_manifest.json').read_text())
    frozen = json.loads((a.task_root / 'frozen_artifact_checksums.json').read_text())
    original_ok = {p: digest(a.task_root / 'baseline' / p) == v['sha256'] for p, v in baseline.items()}
    source_ok = {p: digest(source / p) == v for p, v in current.items()}
    frozen_ok = {p: digest(a.data_root / p) == v for p, v in frozen.items()}
    validated = json.loads((a.task_root / 'validation_22814288/installed_artifact.json').read_text())
    wheel_ok = digest(Path(validated['wheel']['path'])) == validated['wheel']['sha256']
    native_ok = digest(Path(validated['native']['path'])) == validated['native']['sha256']
    v5 = json.loads((a.task_root / 'candidate_v5/source_manifest.json').read_text())
    package_ok = {p: current.get(p) == h for p, h in v5.items()
                  if p.startswith(('python/haplodose/', 'crates/'))}
    changes, patch = [], []
    for p in sorted(set(baseline) | set(current)):
        old, new = baseline.get(p, {}).get('sha256'), current.get(p)
        if old == new:
            continue
        status = 'added' if old is None else 'removed' if new is None else 'modified'
        changes.append(dict(path=str(a.data_root / p), status=status, before=old, after=new))
        before = (a.task_root / 'baseline' / p).read_text().splitlines(True) if old else []
        after = (source / p).read_text().splitlines(True) if new else []
        patch.extend(difflib.unified_diff(before, after, fromfile='baseline/' + p, tofile='candidate/' + p))
    (a.out / 'changes.json').write_text(json.dumps(changes, indent=2) + '\n')
    (a.out / 'implementation.patch').write_text(''.join(patch))
    checks = dict(initial_dirty_snapshot=original_ok, final_source=source_ok,
                  selected_frozen_artifacts=frozen_ok, package_matches_built_v5=package_ok,
                  wheel=wheel_ok, installed_native=native_ok)
    passed = all(original_ok.values()) and all(source_ok.values()) and all(frozen_ok.values())
    passed = passed and all(package_ok.values()) and wheel_ok and native_ok
    (a.out / 'summary.json').write_text(json.dumps(dict(passed=passed, checks=checks,
        source=str(source), source_manifest_sha256=digest(source / 'source_manifest.json'),
        changed_files=len(changes), scope='Listed source and frozen files only; no full historical data rehash'), indent=2) + '\n')
    if not passed:
        raise RuntimeError('Integrity verification failed; inspect summary.json')


if __name__ == '__main__':
    main()
