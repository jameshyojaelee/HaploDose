#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# COLO320 scorer fixture (2026-11-05, after the CNVkit genomic_format fix in score_one.R). No HaploDose output is read.
# A copyKat-derived fixture (make_fixture.R) stands in for both HaploDose arms on the real COLO320 gate. The map adds the
# BCC06post gate anchor with HaploDose paths "none", as in the scoring shards. Expect: run_score_all exits 0; COLO320
# scored in both arms with the CNVkit truth; HaploDose rows == copyKat rows in every flagged table; HaploDose gate_miss NA;
# copyKat gate_miss FALSE; CaSpER gate_miss TRUE (auc_losses gap 0.052 in the gate table).
# Usage: run_fixture_colo320.sbatch OUT
set -euo pipefail
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=COLO320 G=results/schmid_gate3_COLO320_20261104_0641 GA=results/schmid_gate3_BCC06post_20261104_0641
H=/gpfs/commons/home/jameslee/Haplodose
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
test ! -e $OUT/fixture; mkdir -p $OUT/fixture
$MM run -n schmid_r41 Rscript $S/score/test/make_fixture.R $H/$G/pipeline $DS $OUT/fixture > $OUT/make_fixture.log 2>&1
for f in count_matrix.txt sample_annotation.txt ref_groups.txt; do md5sum $(readlink -f $H/$G/pipeline/data/input_$DS/$f) | sed "s#  .*#  /x/$f#"; done > $OUT/fixture/input_md5.txt
F=${OUT#$H/}/fixture
printf "run_id\tdataset\trole\tgate_dir\thd_vnext\thd_sens1424\n%s\t%s\tprimary\t%s\t%s\t%s\nBCC06post_gateanchor\tBCC06post\tprimary\t%s\tnone\tnone\n" \
  $DS $DS $G $F $F $GA > $OUT/map.tsv
rc=0; bash $S/score/run_score_all.sbatch $OUT/score $OUT/map.tsv > $OUT/case.out 2>&1 || rc=$?
$MM run -n schmid_r41 Rscript $S/score/test/check_fixture_colo320.R $OUT/score $rc > $OUT/check_fixture.log 2>&1 || true
cat $OUT/check_fixture.log
