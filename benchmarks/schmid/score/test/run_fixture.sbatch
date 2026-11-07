#!/usr/bin/env bash
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
# End-to-end fixture tests of run_score_all.sbatch on one passed gate (DS), with a copyKat-derived fixture in place of
# both HaploDose arms (no HaploDose output is read). Cases:
#   A: runs DS (real gate, primary) + DS_miss (gate table copy with copykat pearson gap set to 0.07) +
#      DS_nogate (gate dir without a table) + DS_badinput (HaploDose input md5 differs from the gate). Expect: job completes; HaploDose rows == copyKat rows; DS_miss scored with
#      copykat gate_miss TRUE in gate_check.tsv and in every flagged table; DS_nogate and DS_badinput skipped and
#      reported with their reasons, the other runs scored.                                                [(2), (3), integrity]
#      Also DS_nonumbat (gate table copy with both Numbat rows' ours and difference set to NA, as for a failed Numbat):
#      expect numbat and numbat_cnv gate_miss TRUE with max_abs_gap NA (not -Inf), flagged in every table.  [MCF7 rule]
#   B: only DS_miss as primary. Expect: stop before scoring (GATE_CHECK_STOP, exit 3), nothing scored.       [(1)]
#   C (if DS2 GATE_DIR2 given; DS2 = MCF7): DS real (primary) + DS2 line_overlap with a copyKat fixture built from
#      GATE_DIR2. Expect: exit 0; DS2 scored without the Numbat input (no Numbat rows in any table); HaploDose == copyKat;
#      Numbat outputs gate_miss TRUE with gap NA in gate_methods.tsv.                          [STATUS 2026-11-04 ~16:45]
# Usage: run_fixture.sbatch OUT DS GATE_DIR [DS2 GATE_DIR2]
set -euo pipefail
S=${CODE_DIR:-/gpfs/commons/home/jameslee/Haplodose/benchmarks/schmid}
OUT=$1 DS=$2 G=$3 DS2=${4:-} G2=${5:-}
H=/gpfs/commons/home/jameslee/Haplodose
MM=/gpfs/commons/home/jameslee/.local/bin/micromamba
export MAMBA_ROOT_PREFIX=/gpfs/commons/home/jameslee/micromamba R_LIBS_USER=/nonexistent R_PROFILE_USER=/dev/null
test ! -e $OUT/fixture; mkdir -p $OUT/fixture
$MM run -n schmid_r41 Rscript $S/score/test/make_fixture.R $H/$G/pipeline $DS $OUT/fixture > $OUT/make_fixture.log 2>&1
for f in count_matrix.txt sample_annotation.txt ref_groups.txt; do md5sum $(readlink -f $H/$G/pipeline/data/input_$DS/$f) | sed "s#  .*#  /x/$f#"; done > $OUT/fixture/input_md5.txt
R=${OUT#$H/}; F=$R/fixture
mkdir -p $OUT/gate_miss $OUT/gate_none $OUT/fixture_bad $OUT/gate_nonumbat
ln -s $OUT/fixture/haplodose $OUT/fixture_bad/haplodose
awk 'NR==1{$1="00000000000000000000000000000000"} {print $1"  "$2}' $OUT/fixture/input_md5.txt > $OUT/fixture_bad/input_md5.txt
ln -s $H/$G/pipeline $OUT/gate_miss/pipeline; ln -s $H/$G/pipeline $OUT/gate_none/pipeline
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $2=="copykat" && $3=="pearson"{$4=$5+0.07; $6=0.07} {print}' $H/$G/gate_comparison.tsv > $OUT/gate_miss/gate_comparison.tsv
test "$(awk -F'\t' '$6==0.07' $OUT/gate_miss/gate_comparison.tsv | wc -l)" = 1
ln -s $H/$G/pipeline $OUT/gate_nonumbat/pipeline
awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && ($2=="numbat" || $2=="numbat_cnv"){$4=""; $6=""} {print}' $H/$G/gate_comparison.tsv > $OUT/gate_nonumbat/gate_comparison.tsv
test "$(awk -F'\t' 'NR>1 && $6==""' $OUT/gate_nonumbat/gate_comparison.tsv | wc -l)" = 6
hdr="run_id\tdataset\trole\tgate_dir\thd_vnext\thd_sens1424"
printf "$hdr\n%s\t%s\tprimary\t%s\t%s\t%s\n%s_miss\t%s\tprimary\t%s\t%s\t%s\n%s_nogate\t%s\tprimary\t%s\t%s\t%s\n%s_badinput\t%s\tprimary\t%s\t%s\t%s\n" \
  $DS $DS $G $F $F  $DS $DS $R/gate_miss $F $F  $DS $DS $R/gate_none $F $F  $DS $DS $G $R/fixture_bad $R/fixture_bad > $OUT/map_A.tsv
printf "%s_nonumbat\t%s\tprimary\t%s\t%s\t%s\n" $DS $DS $R/gate_nonumbat $F $F >> $OUT/map_A.tsv
printf "$hdr\n%s_miss\t%s\tprimary\t%s\t%s\t%s\n" $DS $DS $R/gate_miss $F $F > $OUT/map_B.tsv
rcA=0; bash $S/score/run_score_all.sbatch $OUT/score_A $OUT/map_A.tsv > $OUT/case_A.out 2>&1 || rcA=$?
rcB=0; bash $S/score/run_score_all.sbatch $OUT/score_B $OUT/map_B.tsv > $OUT/case_B.out 2>&1 || rcB=$?
rcC=NA
if [ -n "$DS2" ]; then
  mkdir -p $OUT/fixture_C
  $MM run -n schmid_r41 Rscript $S/score/test/make_fixture.R $H/$G2/pipeline $DS2 $OUT/fixture_C > $OUT/make_fixture_C.log 2>&1
  for f in count_matrix.txt sample_annotation.txt ref_groups.txt; do md5sum $(readlink -f $H/$G2/pipeline/data/input_$DS2/$f) | sed "s#  .*#  /x/$f#"; done > $OUT/fixture_C/input_md5.txt
  printf "$hdr\n%s\t%s\tprimary\t%s\t%s\t%s\n%s\t%s\tline_overlap\t%s\t%s\t%s\n" $DS $DS $G $F $F  $DS2 $DS2 $G2 $R/fixture_C $R/fixture_C > $OUT/map_C.tsv
  rcC=0; bash $S/score/run_score_all.sbatch $OUT/score_C $OUT/map_C.tsv > $OUT/case_C.out 2>&1 || rcC=$?
fi
echo -e "case_A_exit\t$rcA\ncase_B_exit\t$rcB\ncase_C_exit\t$rcC" > $OUT/case_exits.tsv
$MM run -n schmid_r41 Rscript $S/score/test/check_fixture.R $OUT $DS $rcA $rcB $rcC ${DS2:-none} > $OUT/check_fixture.log 2>&1 || true
cat $OUT/check_fixture.log
