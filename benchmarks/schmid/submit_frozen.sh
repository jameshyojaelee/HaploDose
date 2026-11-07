#!/usr/bin/env bash
# Submit a claim 5 job so that later edits to benchmarks/schmid cannot change it (team-lead rule, 2026-11-04).
# Copies benchmarks/schmid into CODE_OUT/code/benchmarks_schmid at submission, writes the sha256 of every copied file
# to CODE_OUT/code/code_sha256.txt, optionally checks the copy against a frozen manifest, and submits the COPY of the
# sbatch script with CODE_DIR pointing at the copy (scripts resolve their helpers through CODE_DIR).
# Usage: submit_frozen.sh CODE_OUT SBATCH_REL_PATH [sbatch options ...] -- [script args ...]
#   SBATCH_REL_PATH is relative to benchmarks/schmid (e.g. score/run_score_all.sbatch).
#   FROZEN_MANIFEST=<code_sha256.txt>: refuse to submit unless the copy matches it byte for byte.
# Prints the job id.
set -euo pipefail
CODE_OUT=$1 REL=$2; shift 2
OPTS=(); while [ $# -gt 0 ] && [ "$1" != "--" ]; do OPTS+=("$1"); shift; done
[ "${1:-}" = "--" ] && shift
H=/gpfs/commons/home/jameslee/Haplodose
C=$CODE_OUT/code/benchmarks_schmid
test ! -e $CODE_OUT/code
mkdir -p $C
(cd $H/benchmarks/schmid && find . -type f ! -path '*/__pycache__/*' ! -name '*.new' -print0 | cpio -0pdm --quiet $C)
chmod -R a-w $C
(cd $C && find . -type f -print0 | sort -z | xargs -0 sha256sum) > $CODE_OUT/code/code_sha256.txt
if [ -n "${FROZEN_MANIFEST:-}" ]; then
  diff <(awk '{print $1, $2}' $FROZEN_MANIFEST) <(awk '{print $1, $2}' $CODE_OUT/code/code_sha256.txt) > $CODE_OUT/code/frozen_diff.txt \
    || { echo "COPY DIFFERS FROM FROZEN MANIFEST (see $CODE_OUT/code/frozen_diff.txt); not submitted" >&2; exit 4; }
  cp $FROZEN_MANIFEST $CODE_OUT/code/frozen_manifest_used.txt
fi
# CODE_DIR travels in the environment, so an --export=ALL,... in OPTS keeps it
export CODE_DIR=$C
J=$(sbatch --parsable "${OPTS[@]}" $C/$REL "$@")
echo -e "$J\t$REL\t$(date -Is)\t$C" >> $CODE_OUT/code/submitted.tsv
echo $J
