#!/usr/bin/env bash
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --job-name=masld-test
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Claim 4 SPECTRUM Numbat allele inputs for one patient: numbat's own pileup_and_phase.R (10x mode) over all of the
# patient's sample BAMs at once (one person: joint genotype, Eagle2 phasing per chromosome), on the cells in the
# patient's input (prepare_inputs.py barcodes/<sample>.tsv). Panel: 1000 Genomes phase 3, GRCh38 (the CellRanger
# BAMs are GRCh38): external/schmid_env/ref_20261102_2215/{1000G_hg38, genetic_map_hg38_withX.txt.gz,
# genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz}. Kinker used the same 1000G phase 3 panel on hg19 because its
# BAMs are hg19. Pileup arguments are pileup_and_phase.R's own (as Kinker stage A: --minMAF 0 --minCOUNT 2,
# UMItag Auto, cellTAG CB).
# Then each sample's allele table is renamed from BAM barcodes to input cell ids (<sample>_<barcode without -1>):
# OUT/<patient>/<sample>_allele_counts.tsv.gz (Numbat, per specimen) and OUT/<patient>/allele_counts.tsv.gz (all
# samples; HaploDose descriptive allele run).
# Usage: sbatch --output=OUT/logs/<patient>_%j.out pileup_phase.sbatch INPUT_ROOT CR_ROOT OUT PATIENT
set -euo pipefail
IN=$1 CR=$2 ROOT=$3 P=$4
REF=external/schmid_env/ref_20261102_2215
ENV=/gpfs/commons/home/jameslee/CRC/data/interim/numbat_env
HD=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
OUT=$ROOT/$P
[[ ! -s $OUT/allele_counts.tsv.gz ]] || { echo "exists: $OUT"; exit 0; }
mkdir -p "$OUT"
SAMPLES=() BAMS=() BCS=()
for f in "$IN/$P"/barcodes/*.tsv; do
  s=$(basename "$f" .tsv)
  test -e "$CR/$s/done.ok" || { echo "CellRanger not done: $s" >&2; exit 1; }
  SAMPLES+=("$s"); BAMS+=("$CR/$s/outs/possorted_genome_bam.bam"); BCS+=("$f")
done
join() { local IFS=,; echo "$*"; }
( unset LD_LIBRARY_PATH PYTHONPATH; export PATH="$ENV/bin:$PATH" R_LIBS_USER="$ENV/lib/R/library" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
  cellsnp-lite --version 2>&1 | head -1 > "$OUT/versions.txt" || true
  /usr/bin/time -v "$ENV/bin/Rscript" --vanilla "$ENV/lib/R/library/numbat/bin/pileup_and_phase.R" \
    --label "$P" --samples "$(join "${SAMPLES[@]}")" --bams "$(join "${BAMS[@]}")" --barcodes "$(join "${BCS[@]}")" \
    --gmap "$REF/genetic_map_hg38_withX.txt.gz" --eagle "$ENV/bin/eagle" \
    --snpvcf "$REF/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz" --paneldir "$REF/1000G_hg38" \
    --outdir "$OUT/numbat_out" --ncores 16 ) 2> "$OUT/pileup_and_phase.log"
$HD -s - "$OUT" "${SAMPLES[@]}" <<'EOF'
import sys
from pathlib import Path
import pandas as pd
out, samples = Path(sys.argv[1]), sys.argv[2:]
parts = []
for s in samples:
    a = pd.read_csv(out / 'numbat_out' / f'{s}_allele_counts.tsv.gz', sep='\t', dtype={'CHROM': str})
    assert len(a) > 0 and a.cell.str.endswith('-1').all(), s   # 0 rows: e.g. a contig-name mismatch in the pileup
    a['cell'] = s + '_' + a.cell.str[:-2]
    a.to_csv(out / f'{s}_allele_counts.tsv.gz', sep='\t', index=False, compression='gzip')
    parts.append(a)
    print(s, a.cell.nunique(), 'cells with het SNP counts', len(a), 'rows')
pd.concat(parts).to_csv(out / 'allele_counts.tsv.gz', sep='\t', index=False, compression='gzip')
EOF
