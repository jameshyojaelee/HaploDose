#!/usr/bin/env bash
#SBATCH --job-name=kinker_numbat
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020 Numbat comparator for one pool: stage B (phase_pool.R: per-line genotype + Eagle2 phasing from the
# stage-A channel pileups), then stage C: Numbat (benchmarks/comparators/run_numbat.py, genome hg19) on each
# single-line specimen of the pool in both rotations. Reference = the pool's fit reference lines
# (--reference-scope dataset), the same cells the other tools use. Each run's cell table keeps only its own
# specimen as single_line, so evaluate.py --tool numbat scores each specimen once.
# Usage: sbatch -J kinker_numbat_<pool> numbat_pool.sbatch POOL PILEUP_ROOT OUT_ROOT [SPECIMENS]
#   SPECIMENS: optional space-separated subset (smoke test). Channels: the pool's titles under PILEUP_ROOT.
set -euo pipefail
POOL=$1 PU=$2 ROOT=$3 ONLY=${4:-}
IN=results/kinker_input_v2_20261101_1551
P=data/external/numbat_hg19_panel
ENV=/gpfs/commons/home/jameslee/CRC/data/interim/numbat_env
CNV=/gpfs/commons/home/jameslee/micromamba/envs/scrna_cnv/bin/python
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$ROOT/alleles" "$ROOT/cells" "$ROOT/runs"
CHANNELS=$(cd "$PU" && ls -d Pool*_"$POOL" | tr '\n' ' ')
for ch in $CHANNELS; do test -s "$PU/$ch/cellSNP.tag.DP.mtx" || { echo "pileup missing: $ch" >&2; exit 1; }; done
( unset LD_LIBRARY_PATH; export PATH="$ENV/bin:$PATH" R_LIBS_USER="$ENV/lib/R/library" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  /usr/bin/time -v "$ENV/bin/Rscript" --vanilla benchmarks/kinker/alleles/phase_pool.R "$PU" "$ROOT/alleles" \
    "$P/genetic_map_hg19_withX.txt.gz" "$ENV/bin/eagle" "$P/panel" 8 $CHANNELS ) 2> "$ROOT/phase_pool$POOL.log"
for r in 0 1; do
  for S in $(awk -F"\t" -v p="$POOL" 'NR>1 && $11=="single_line" && $2==p {print $10}' "$IN/cells_r$r.tsv" | sort -u); do
    [[ -z $ONLY || " $ONLY " == *" $S "* ]] || continue
    A=$ROOT/alleles/$S/${S}_allele_counts.tsv.gz       # specimen = line for single-line queries
    [[ -s $A ]] || { echo "no alleles for $S" >&2; continue; }
    C=$ROOT/cells/r${r}_$S.tsv
    test -f "$C" || /gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s -c "
import pandas as pd
c = pd.read_csv('$IN/cells_r$r.tsv', sep='\t', dtype={'pool': str})
o = c.design.eq('single_line') & c.specimen.ne('$S')
c.loc[o, ['role', 'ref_group', 'ref_split', 'specimen', 'patient']] = ['other', '', '', '', '']
c.loc[o, 'design'] = 'unused'
c.to_csv('$C', sep='\t', index=False)"
    D=$ROOT/runs/r${r}_pool$POOL/$S
    ls -d "$D"/numbat_* >/dev/null 2>&1 && continue
    mkdir -p "$D"
    /usr/bin/time -v $CNV -s benchmarks/comparators/run_numbat.py --counts "$IN/counts.h5" --cells "$C" \
      --positions "$IN/positions.tsv" --alleles "$A" --specimen "$S" --out-parent "$D" --genome hg19 \
      --reference-scope dataset 2> "$D/time.txt" || { echo "numbat failed $S (r$r)"; tail -20 "$D/time.txt"; }
  done
done
$CNV -s -m pip freeze > "$ROOT/pip_freeze_scrna_cnv.txt"
