#!/usr/bin/env bash
#SBATCH --job-name=haplodose_validate
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --constraint=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
# Build/test a snapshot without touching its source or the shared environment.
# Usage: sbatch --output=TASK/validate_%j.log hardening_validate.sbatch TASK SNAPSHOT ENV [pytest args]
set -euo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Run this launcher on an allocated compute node.' >&2; exit 2; }
TASK=$(realpath "$1")
SNAPSHOT=$(realpath "$2")
ENV=$(realpath "$3")
shift 3
case "$ENV/" in "$TASK/"*) ;; *) echo 'Validation environment must be inside TASK.' >&2; exit 2 ;; esac
[[ -f $SNAPSHOT/pyproject.toml && -x $ENV/bin/python ]] || { echo 'Missing snapshot or isolated Python.' >&2; exit 2; }
OUT=$TASK/validation_$SLURM_JOB_ID
[[ ! -e $OUT ]] || { echo "Refusing to reuse $OUT" >&2; exit 2; }
mkdir "$OUT"
record_exit() {
  status=$?
  printf 'exit\t%s\nend\t%s\n' "$status" "$(date -Is)" >> "$OUT/job.tsv"
  trap - EXIT
  exit "$status"
}
trap record_exit EXIT
printf 'snapshot\t%s\nenvironment\t%s\nhost\t%s\nstart\t%s\n' "$SNAPSHOT" "$ENV" "$(hostname)" "$(date -Is)" > "$OUT/job.tsv"
# Prevent concurrent wheel installs/tests in the same isolated prefix.
exec 9>"$ENV/.hardening-validate.lock"
flock -n 9 || { echo 'Isolated environment is already being validated.' >&2; exit 2; }
PY=$ENV/bin/python
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_PROTECTED_PATHS="$SNAPSHOT:/gpfs/commons/home/jameslee/Haplodose/results:/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export PATH="$ENV/bin:$HOME/.cargo/bin:$PATH" PYO3_PYTHON=$PY CARGO_NET_OFFLINE=true
export CARGO_TARGET_DIR=$OUT/cargo_target
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
"$PY" -I - "$SNAPSHOT" "$OUT/source" <<'PY'
import hashlib, json, shutil, sys
from pathlib import Path
source, destination = map(Path, sys.argv[1:])
inventory = {}
for path in sorted(source.rglob('*')):
    if path.is_file():
        inventory[str(path.relative_to(source))] = hashlib.sha256(path.read_bytes()).hexdigest()
(destination.parent / 'snapshot_sha256.json').write_text(json.dumps(inventory, indent=2) + '\n')
shutil.copytree(source, destination, ignore=shutil.ignore_patterns('*.so', '__pycache__', 'target', '.pytest_cache'))
for relative, digest in inventory.items():
    path = source / relative
    if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
        raise RuntimeError(f'Snapshot changed during copy: {relative}')
    copied = destination / relative
    if copied.is_file() and hashlib.sha256(copied.read_bytes()).hexdigest() != digest:
        raise RuntimeError(f'Copied snapshot differs: {relative}')
PY
mkdir "$OUT/wheels" "$OUT/neutral"
rustc --version > "$OUT/rustc.txt"
(cd "$OUT/source/crates/core" && cargo test --release --offline --locked) > "$OUT/cargo_test.txt" 2>&1
"$PY" -I -m pip wheel --no-index --no-build-isolation --no-deps --wheel-dir "$OUT/wheels" "$OUT/source" > "$OUT/wheel_build.txt" 2>&1
WHEELS=("$OUT"/wheels/haplodose-*.whl)
[[ ${#WHEELS[@]} == 1 && -f ${WHEELS[0]} ]] || { echo 'Expected exactly one wheel.' >&2; exit 2; }
"$PY" -I -m pip install --no-index --no-deps --force-reinstall "${WHEELS[0]}" > "$OUT/wheel_install.txt" 2>&1
cd "$OUT/neutral"
"$PY" -I - "$ENV" "${WHEELS[0]}" > "$OUT/installed_artifact.json" <<'PY'
import hashlib, importlib.metadata, json, sys
from pathlib import Path
import haplodose
import haplodose._core as native
prefix, wheel = map(lambda x: Path(x).resolve(), sys.argv[1:])
artifacts = {}
for name, module in [('package', haplodose), ('native', native)]:
    path = Path(module.__file__).resolve()
    if not path.is_relative_to(prefix):
        raise RuntimeError(f'{name} imported from outside the isolated environment: {path}')
    artifacts[name] = dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
artifacts['wheel'] = dict(path=str(wheel), sha256=hashlib.sha256(wheel.read_bytes()).hexdigest())
artifacts['python'] = sys.executable
artifacts['versions'] = {d.metadata['Name']: d.version for d in importlib.metadata.distributions() if d.metadata['Name']}
print(json.dumps(artifacts, indent=2))
PY
"$PY" -I -m pip freeze > "$OUT/pip_freeze.txt"
if [[ $# == 0 ]]; then set -- "$OUT/source/tests"; fi
# Relative test paths are interpreted against the copied snapshot, never the live checkout.
ARGS=()
for arg in "$@"; do
  case "$arg" in tests|tests/*) ARGS+=("$OUT/source/$arg") ;; *) ARGS+=("$arg") ;; esac
done
/usr/bin/time -v -o "$OUT/pytest_resources.txt" "$PY" -I -m pytest --basetemp "$OUT/pytest_tmp" --import-mode=importlib -p no:cacheprovider -q "${ARGS[@]}" > "$OUT/pytest.txt" 2>&1
