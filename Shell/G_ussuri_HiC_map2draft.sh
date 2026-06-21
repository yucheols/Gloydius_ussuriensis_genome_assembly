#!/bin/bash
#SBATCH --job-name=hic_map
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# make sure the job will stop if any step fails
set -euo pipefail

# set directories
indir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined"
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding"

# set variables
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
R1="${indir}/Gloydius_ussuriensis_HiC_R1_trimmed.fastq.gz"
R2="${indir}/Gloydius_ussuriensis_HiC_R2_trimmed.fastq.gz"
BAM_PREFIX="${outdir}/Gloydius_ussuriensis_HiC_to_draft"

# run bwa mem
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} "${REF}" "${R1}" "${R2}" | \
  samtools view -@ 16 -bS - | \
  samtools sort -@ 16 -o "${BAM_PREFIX}.sorted.bam"

samtools index "${BAM_PREFIX}.sorted.bam"

# name-sort for fixmate
samtools sort -@ 16 -n \
  -o "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.sorted.bam"

# add mate information
samtools fixmate -@ 16 -m \
  "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

# coordinate-sort again
samtools sort -@ 16 \
  -o "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

# mark duplicates
samtools markdup -@ 16 \
  "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.markdup.bam"

samtools index "${BAM_PREFIX}.markdup.bam"

# produce basic mapping summary
samtools flagstat "${BAM_PREFIX}.markdup.bam" > "${BAM_PREFIX}.markdup.flagstat.txt"