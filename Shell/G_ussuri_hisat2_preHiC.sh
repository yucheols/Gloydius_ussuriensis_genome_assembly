#!/bin/bash
#SBATCH --job-name=hisat2_preHiC
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=180:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
export PATH=/home/yshin/mendel-nas1/miniconda3/bin:$PATH
conda activate hisat2

# path to draft
draft_path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito

# path to RNAseq reads
rna_path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads

# output directory
outpath=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/hisat2_preHiC
mkdir -p ${outpath}

# build a genome index
hisat2-build ${draft_path}/Gloydius_ussuriensis_AMNH_21010_noMito.fa ${outpath}/Gloydius_ussuriensis_AMNH_21010_noMito

# loop over the RNA reads per tissue type and map to the reference genome using HiSat2
for type in ${rna_path}/*; do
  if [ -d "$type" ]; then
    echo "Mapping RNA-seq reads from the following tissue type: $type"
    hisat2 -p ${SLURM_CPUS_PER_TASK} -x ${outpath}/Gloydius_ussuriensis_AMNH_21010_noMito \
    -1 ${rna_path}/AMNH_21010_${type}_R1_paired.fastq.gz -2 ${rna_path}/AMNH_21010_${type}_R2_paired.fastq.gz | \
    samtools sort -o ${outpath}/${type}.bam
  fi
done

# index the bam file
samtools index ${outpath}/rnaseq.bam

# print when done
echo "[DONE] HISAT2 mapping complete."
echo "[INFO] BAMs in: ${outpath}/bam"
echo "[INFO] Logs in: ${outpath}/logs"