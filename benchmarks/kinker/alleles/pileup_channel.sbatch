#!/usr/bin/env bash
#SBATCH --job-name=kinker_pileup
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/kinker/%x_%j.err
# Kinker 2020 Numbat inputs, stage A: one cellsnp-lite pileup per channel BAM over all of that channel's Kinker
# cells (every line in the pool), with the arguments numbat's pileup_and_phase.R uses in 10x mode. Stage B
# (phase_pool.R) splits it per line, since each line is a different person and is genotyped and phased alone.
# A channel BAM is ENA title Pool<channel>_<pool>; Kinker cell_id = <barcode>-<channel>-<pool>; BAM CB = <barcode>-1
# (checked on SRR12557906, Pool6_6). Waits for the BAM's .ok (md5 checked by fetch_bams.sbatch).
# Usage: sbatch -J kinker_pileup_<title> pileup_channel.sbatch RUN TITLE OUT_ROOT   (writes OUT_ROOT/<title>/)
set -euo pipefail
RUN=$1 TITLE=$2 ROOT=$3
BAM=data/external/kinker2020_bams/$RUN.bam
CELLS=results/kinker_input_v2_20261101_1551/cells_r0.tsv   # same cells and ids in both rotations
SNPVCF=${SNPVCF:-data/external/numbat_hg19_panel/snps_AF5e2_chr1to22_hg19.vcf.gz}   # override: smoke tests only
ENV=/gpfs/commons/home/jameslee/CRC/data/interim/numbat_env
OUT=$ROOT/$TITLE
# ALLOW_PARTIAL=1 (smoke tests only): a BAM still downloading, with a complete index; SNPVCF must then cover only
# chromosomes already in the file.
[[ ${ALLOW_PARTIAL:-0} == 1 || ( -e $BAM.ok && -e $BAM.bai.ok ) ]] || { echo "BAM not complete: $BAM" >&2; exit 1; }
[[ ! -e $OUT/cellSNP.tag.DP.mtx ]] || { echo "exists: $OUT"; exit 0; }
mkdir -p "$OUT"
CH=${TITLE#Pool}; CH=${CH%_*}; POOL=${TITLE##*_}
# Barcodes in BAM form and the map back to Kinker cell ids and lines.
/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python -s - "$CELLS" "$CH" "$POOL" "$OUT" <<'EOF'
import sys
import pandas as pd
cells, ch, pool, out = sys.argv[1:]
c = pd.read_csv(cells, sep='\t', dtype=str)
c = c[c.pool.eq(pool) & c.cell_id.str.endswith(f'-{ch}-{pool}') & c.line.notna()]
m = pd.DataFrame({'bam_barcode': c.cell_id.str.split('-').str[0] + '-1', 'cell_id': c.cell_id, 'line': c.line})
assert m.bam_barcode.is_unique and len(m) > 0, (len(m), m.bam_barcode.duplicated().sum())
m.to_csv(f'{out}/cell_map.tsv', sep='\t', index=False)
m.bam_barcode.to_csv(f'{out}/barcodes.tsv', index=False, header=False)
print(f'{len(m)} cells, {m.line.nunique()} lines')
EOF
export PATH="$ENV/bin:$PATH" OMP_NUM_THREADS=1
unset LD_LIBRARY_PATH
cellsnp-lite --version 2>&1 | head -1 > "$OUT/versions.txt" || true
/usr/bin/time -v cellsnp-lite -s "$BAM" -b "$OUT/barcodes.tsv" -O "$OUT" -R "$SNPVCF" -p 8 --minMAF 0 --minCOUNT 2 \
  --UMItag Auto --cellTAG CB 2> "$OUT/cellsnp_time.txt"
test -s "$OUT/cellSNP.tag.DP.mtx"
echo "done $TITLE: $(grep -vc '^#' "$OUT/cellSNP.base.vcf") SNPs"
