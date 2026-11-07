#!/usr/bin/env bash
#SBATCH --account=innovation
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# Save the claim 5 external reference cell IDs (STATUS item 5, amendment f). The COLO320 draw is made by R 4.1.0
# in schmid_r41 (the paper's R version) from n alone; the Python step applies it to the atlas obs in row order.
# Usage: sbatch --output=OUT/select_%j.out benchmarks/schmid/reference/run_select.sbatch OUT
set -euo pipefail
OUT=$1
H=/gpfs/commons/home/jameslee/Haplodose
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
PY=/gpfs/commons/home/jameslee/micromamba/envs/haplodose/bin/python
# n is counted from the atlas obs first (Python), then R draws the indices
N=$($PY -I -c "
import sys; sys.argv=['x']
sys.path.insert(0, '$H/benchmarks/schmid/reference')
from select_reference_cells import obs, GUT
g = obs(GUT, ['Diagnosis', 'Region', 'category'])
print(int((g.Diagnosis.eq('Healthy adult') & g.Region.eq('LargeInt') & g.category.eq('Epithelial')).sum()))")
echo "n_eligible $N"
$MM run -n schmid_r41 Rscript -e "cat(R.version.string, RNGkind(), '\n'); set.seed(1); writeLines(as.character(sample(1:$N, 5000)), '$OUT/colo320_r_indices.txt')" > $OUT/r_draw.log 2>&1
$PY -I $H/benchmarks/schmid/reference/select_reference_cells.py $OUT $OUT/colo320_r_indices.txt
echo SELECT_DONE
