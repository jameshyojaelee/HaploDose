#!/usr/bin/env bash
#SBATCH --job-name=h2_fetch_schmid
#SBATCH --account=innovation
#SBATCH --partition=io
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=72:00:00
#SBATCH --chdir=/gpfs/commons/home/jameslee/Haplodose
#SBATCH --output=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.out
#SBATCH --error=/gpfs/commons/home/jameslee/Haplodose/logs/heldout2/%x_%j.err
# Claim 5 (Schmid et al. 2025 benchmark) raw inputs. Accessions from the authors' repo at commit
# 319ea34ddfc46bdd4c0cb1ed534f4ece0f5c44ca (benchmarks/heldout2/schmid_repo/preprocessing_datasets/download_*.sh):
#   COLO320-HSR  SRR12900668 (GSE160148, 10x Multiome GEX, cellranger --chemistry=ARC-v1)
#   BCC06        SRR8315756  (GSE123813 su006_pre_Tumor_RNA);  BCC06post SRR8315759 (su006_post_Tumor_RNA)
#   MCF7         SRR10018149 (MCF7-AA-t0)
#   MM199        SRR10574405-12 (GSE141299; ENA serves one fastq per run, so SRA-Toolkit fasterq-dump
#                --split-files --include-technical gives R1/R2/I1 as in the authors' fastq-dump --split-files)
#   A375, HCT116 DNTR-seq SRP245376 (Zachariadis 2020; all 4869 runs: RNA and scDNA plates, both lines)
# Plus GEO processed files used by the authors' prepare_*.R (BCC count matrix + annotation, MM199 sciCNV matrices,
# A375 reference GSE151091, HCT116 reference GSE95435). Cell Ranger is not run here.
# Total ENA 152 GB (benchmarks/heldout2/schmid_runs.tsv) + MM199 ~11 GB compressed + GEO ~0.2 GB.
set -euo pipefail
D=data/external/schmid2025_raw
bash benchmarks/heldout2/ena_fetch.sh benchmarks/heldout2/schmid_runs.tsv $D 4 || echo "ENA part incomplete"
G=https://ftp.ncbi.nlm.nih.gov/geo/series
for u in $G/GSE123nnn/GSE123813/suppl/GSE123813_bcc_scRNA_counts.txt.gz \
         $G/GSE123nnn/GSE123813/suppl/GSE123813_bcc_all_metadata.txt.gz \
         $G/GSE141nnn/GSE141299/suppl/GSE141299_RAW.tar \
         $G/GSE151nnn/GSE151091/suppl/GSE151091_raw_matrix.csv.gz \
         $G/GSE151nnn/GSE151091/suppl/GSE151091_Metadata.csv.gz \
         $G/GSE95nnn/GSE95435/suppl/GSE95435_P150057_full_gene_count_table.txt.gz; do
  bash benchmarks/heldout2/geo_fetch.sh $u $D/geo
done
# MM199 via SRA (prefetch verifies the .sra md5; fasterq-dump splits technical reads)
module load SRA-Toolkit/3.2.1
M=$D/fastq/MM199; mkdir -p $M/sra $M/tmp
for r in $(cat benchmarks/heldout2/schmid_mm199_sra_runs.txt); do
  [[ -e $M/$r.ok ]] && continue
  prefetch --max-size 200G --output-directory $M/sra $r
  sra=$(ls $M/sra/$r/$r.sra* | head -1)
  vdb-validate $sra
  fasterq-dump --split-files --include-technical --threads 4 --temp $M/tmp --outdir $M $sra
  for x in $M/${r}_*.fastq; do gzip -f "$x"; done
  for x in $M/${r}_*.fastq.gz; do
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "${x#$D/}" "sra:$r" "$(stat -c %s $x)" "$(md5sum $x | cut -d' ' -f1)" \
      "local_after_vdb-validate" "$(date -Iseconds)"
  done > $M/$r.ok
done
cat $D/geo/*.ok $M/*.ok >> $D/manifest.tsv
# GSE144296 view (claim 4 false-split control) = DNTRseq_ subset of this manifest (added 2026-11-02 after job 22612843 started)
{ head -1 $D/manifest.tsv; grep -P "^fastq/DNTRseq_" $D/manifest.tsv; } > data/external/dntrseq_gse144296/manifest.tsv
echo done
