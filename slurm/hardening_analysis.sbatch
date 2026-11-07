#!/usr/bin/env bash
#SBATCH --job-name=haplodose_evidence
#SBATCH --account=innovation
#SBATCH --qos=innovation
#SBATCH --partition=cpu
#SBATCH --constraint=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
# Usage: sbatch [resources] script TASK ENV SNAPSHOT SCRIPT OUTPUT [script args]
set -euo pipefail
[[ -n ${SLURM_JOB_ID:-} ]] || { echo 'Use an allocated compute node.' >&2; exit 2; }
TASK=$(realpath "$1")
ENV=$(realpath "$2")
SOURCE=$(realpath "$3")
SCRIPT=$4
OUTPUT=$5
shift 5
case "$ENV/" in "$TASK/"*) ;; *) echo 'Environment must be isolated under TASK.' >&2; exit 2 ;; esac
[[ -f $SOURCE/SNAPSHOT_COMPLETE && -f $SOURCE/$SCRIPT ]] || { echo 'Incomplete source snapshot.' >&2; exit 2; }
[[ -f ${HAPLODOSE_VALIDATED_ARTIFACT:-} ]] || { echo 'Supply HAPLODOSE_VALIDATED_ARTIFACT from a recorded wheel build/import check.' >&2; exit 2; }
RECORD=$TASK/analysis_$SLURM_JOB_ID
mkdir "$RECORD"
record_exit() {
  status=$?
  printf 'exit\t%s\nend\t%s\n' "$status" "$(date -Is)" >> "$RECORD/job.tsv"
  trap - EXIT
  exit "$status"
}
trap record_exit EXIT
printf 'source\t%s\nscript\t%s\noutput\t%s\nhost\t%s\nstart\t%s\n' "$SOURCE" "$SCRIPT" "$OUTPUT" "$(hostname)" "$(date -Is)" > "$RECORD/job.tsv"
exec 9>"$ENV/.hardening-validate.lock"
flock -sn 9 || { echo 'Environment is being changed by a validation job.' >&2; exit 2; }
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
export XLA_PYTHON_CLIENT_PREALLOCATE=false
export HAPLODOSE_WHEEL_ONLY=1 HAPLODOSE_PYTHON=$ENV/bin/python
export HAPLODOSE_DATA_ROOT=/gpfs/commons/home/jameslee/Haplodose
export HAPLODOSE_PROTECTED_PATHS="$SOURCE:$HAPLODOSE_DATA_ROOT/results:$HAPLODOSE_DATA_ROOT/data:/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets"
export LD_LIBRARY_PATH="$ENV/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export PATH="$ENV/bin:$HOME/.cargo/bin:$PATH"
if [[ ${SLURM_JOB_PARTITION:-cpu} == gpu ]]; then
  export JAX_PLATFORMS=cuda
  nvidia-smi > "$RECORD/gpu_before.txt"
else
  export JAX_PLATFORMS=cpu
fi
cd "$RECORD"
"$ENV/bin/python" -I - "$SOURCE" "$ENV" "$HAPLODOSE_VALIDATED_ARTIFACT" "$SCRIPT" "$OUTPUT" "$@" > "$RECORD/execution.json" <<'PY'
import hashlib, importlib.metadata, json, os, sys
from pathlib import Path
import haplodose
import haplodose._core as native
source, environment = map(Path, sys.argv[1:3])
validated_path = Path(sys.argv[3])
validated = json.loads(validated_path.read_text())
test_report = validated_path.parent / 'pytest.txt'
test_text = test_report.read_text() if test_report.is_file() else ''
test_evidence = dict(path=str(test_report),
    sha256=hashlib.sha256(test_text.encode()).hexdigest() if test_text else None,
    final_line=test_text.strip().splitlines()[-1] if test_text.strip() else 'not available',
    scope='Build/import identity does not imply the test suite passed; inspect the linked report.')
if hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest() != validated['native']['sha256']:
    raise RuntimeError('Native binary differs from the validated wheel')
if hashlib.sha256(Path(validated['wheel']['path']).read_bytes()).hexdigest() != validated['wheel']['sha256']:
    raise RuntimeError('Validated wheel changed')
expected = json.loads((source / 'source_manifest.json').read_text())
for relative, digest in expected.items():
    if hashlib.sha256((source / relative).read_bytes()).hexdigest() != digest:
        raise RuntimeError(f'Changed frozen source: {relative}')
for module in (haplodose, native):
    if not Path(module.__file__).resolve().is_relative_to(environment):
        raise RuntimeError(f'Import escaped isolated environment: {module.__file__}')
# Assert that the wheel matches the frozen Python source. Native binary identity is recorded too.
package = Path(haplodose.__file__).parent
for relative, digest in expected.items():
    if relative.startswith('python/haplodose/') and relative.endswith(('.py', '.tsv')):
        path = package / relative.removeprefix('python/haplodose/')
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError(f'Installed source does not match snapshot: {relative}')
print(json.dumps(dict(argv=sys.argv[4:], source=str(source), python=sys.executable,
    validated_artifact=str(validated_path), wheel=validated['wheel'], test_evidence=test_evidence,
    package=str(package), version=haplodose.__version__,
    native_sha256=hashlib.sha256(Path(native.__file__).read_bytes()).hexdigest(),
    snapshot=(source/'SNAPSHOT_COMPLETE').read_text().strip(),
    versions={d.metadata.get('Name',''): d.version for d in importlib.metadata.distributions()},
    allocation={k:v for k,v in os.environ.items() if k.startswith('SLURM_')}), indent=2))
PY
/usr/bin/time -v -o "$RECORD/resources.txt" "$ENV/bin/python" -I "$SOURCE/$SCRIPT" --out "$OUTPUT" "$@" > "$RECORD/output.txt" 2>&1
if [[ ${SLURM_JOB_PARTITION:-cpu} == gpu ]]; then nvidia-smi > "$RECORD/gpu_after.txt"; fi
