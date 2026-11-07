#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --job-name=masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 4 SPECTRUM: one tool on one patient, then the frozen scorer. Run once (docs/STATUS.md claim 4).
# - haplodose: frozen v-next (results/vnext_freeze_20261103_0744, SHA-256 checked) with the 'cov' flags, no alleles,
#   default seed. haplodose_alleles: the same plus --alleles (descriptive only; not in the paired test).
# - infercnvpy, copykat, scevan: benchmarks/comparators runners, --rule cluster, one run over all the patient's
#   specimens (clones and malignant calls are per specimen in cnv_post).
# - numbat: run_numbat.py per specimen (--reference-scope patient, --genome hg38, the patient's fit reference cells),
#   then merge_runs.py joins the specimens into one layout.
# Every benchmark/comparator script is checked against PLAN/code_sha256.txt first; any mismatch stops the job.
# A failed tool leaves no scorable run; score.py then writes status no_output (paired_test scores it 0.5).
# Usage: sbatch --output=ROOT/logs/<patient>_<tool>_%j.out [--cpus-per-task/--mem] run_tool.sbatch ROOT PATIENT TOOL
#   ROOT holds input/ (prepare_inputs.py), alleles/ (pileup_phase.sbatch); writes runs/<patient>/<tool>/, scores/.
set -euo pipefail
ROOT=$1 P=$2 TOOL=$3
H=/gpfs/commons/home/jameslee/Haplodose
PLAN=$H/dev/spectrum_run_plan_20261103_0750
FREEZE=$H/results/vnext_freeze_20261103_0744
SCORER=$H/results/spectrum_scorer_20261103_amend2_1338   # amendments 2026-11-03 13:30 and 13:38, docs/STATUS.md
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
(cd $H && sha256sum --quiet -c $PLAN/code_sha256.txt) || { echo "eval code hash mismatch" >&2; exit 2; }
(cd $H && sha256sum --quiet -c $SCORER/eval_code_sha256.txt) || { echo "scorer hash mismatch" >&2; exit 2; }
(cd $FREEZE/code && sha256sum --quiet -c ../code_sha256.txt) || { echo "frozen HaploDose hash mismatch" >&2; exit 2; }
IN=$ROOT/input/$P A=$ROOT/alleles/$P OUT=$ROOT/runs/$P/$TOOL DEST=$ROOT/scores/$P/$TOOL
test ! -e "$DEST" || { echo "already scored: $DEST"; exit 0; }
mkdir -p "$OUT"
echo -e "host\t$(hostname)\njob\t$SLURM_JOB_ID\ntool\t$TOOL\npatient\t$P\ncpu\t$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2)\nstart\t$(date -Is)" > "$OUT/job.tsv"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
ARGS=(--counts "$IN/counts.h5" --cells "$IN/cells.tsv" --positions "$IN/positions.tsv")
RUN=
case $TOOL in
  haplodose|haplodose_alleles)
    EXTRA=(); [[ $TOOL == haplodose_alleles ]] && EXTRA=(--alleles "$A/allele_counts.tsv.gz")
    ( export PYTHONPATH=$FREEZE/code/python JAX_PLATFORMS=cpu XLA_PYTHON_CLIENT_PREALLOCATE=false
      export HAPLODOSE_PROTECTED_PATHS=/gpfs/commons/home/jameslee/CRC:/gpfs/commons/groups:/gpfs/commons/datasets
      $HD -s -c "import haplodose; print(haplodose.__file__)" > "$OUT/haplodose_code_path.txt"
      /usr/bin/time -v $HD -s -m haplodose.cli run "${ARGS[@]}" --out "$OUT/run" "${EXTRA[@]}" \
        --split-mode evenodd --split-null-draws 19 --split-null-model cov ) 2> "$OUT/time.txt" \
      || { echo "failed $TOOL $P"; tail -20 "$OUT/time.txt"; }
    [[ -f $OUT/run/qc.json ]] && RUN=$OUT/run ;;
  infercnvpy|copykat|scevan)
    unset PYTHONPATH PYTHONHOME R_LIBS R_LIBS_USER
    if [[ $TOOL == infercnvpy ]]; then CMD=(benchmarks/comparators/run_infercnvpy.py)
    else CMD=(benchmarks/comparators/run_r_field.py --tool $TOOL); fi
    /usr/bin/time -v $CNV -s "${CMD[@]}" "${ARGS[@]}" --out-parent "$OUT" --rule cluster 2> "$OUT/time.txt" \
      || { echo "failed $TOOL $P"; tail -20 "$OUT/time.txt"; }
    RUN=$(ls -d "$OUT/${TOOL}_2"*/ 2>/dev/null | tail -1); RUN=${RUN%/}
    [[ -n $RUN && -f $RUN/qc.json ]] || RUN= ;;
  numbat)
    unset PYTHONPATH PYTHONHOME R_LIBS R_LIBS_USER
    PAIRS=()
    for S in $(awk -F'\t' 'NR>1 && $4=="query" {print $2}' "$IN/cells.tsv" | sort -u); do
      D=$OUT/$S; mkdir -p "$D"
      /usr/bin/time -v $CNV -s benchmarks/comparators/run_numbat.py "${ARGS[@]}" --alleles "$A/${S}_allele_counts.tsv.gz" \
        --specimen "$S" --out-parent "$D" --genome hg38 --reference-scope patient 2> "$D/time.txt" \
        || { echo "numbat failed $S"; tail -20 "$D/time.txt"; }
      R=$(ls -d "$D"/numbat_2*/ 2>/dev/null | tail -1)
      PAIRS+=("$S=${R%/}")
    done
    $HD -s benchmarks/spectrum/merge_runs.py "$OUT/merged" "${PAIRS[@]}"
    [[ -f $OUT/merged/qc.json ]] && RUN=$OUT/merged ;;
  *) echo "unknown tool $TOOL" >&2; exit 2 ;;
esac
echo -e "end\t$(date -Is)\nrun\t${RUN:-none}" >> "$OUT/job.tsv"
# Score (a missing run is scored no_output, as the protocol says).
$HD -s benchmarks/spectrum/score.py "${RUN:-$OUT/none}" "$SCORER/dna_clones/$P/clone_bins.tsv" "$P" hg19 "$DEST"
