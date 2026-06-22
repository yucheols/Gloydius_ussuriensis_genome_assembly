#!/bin/bash
#SBATCH --job-name=hic_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# make sure the job will stop if any step fails
set -euo pipefail

# set directories
indir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined"
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding"

# set variables
REF="${outdir}/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
R1="${indir}/Gloydius_ussuriensis_HiC_R1_trimmed.fastq.gz"
R2="${indir}/Gloydius_ussuriensis_HiC_R2_trimmed.fastq.gz"
BAM_PREFIX="${outdir}/Gloydius_ussuriensis_HiC_to_draft"

# run bwa mem
echo "Starting BWA mapping: $(date)"
bwa mem -5SP -t "${SLURM_CPUS_PER_TASK}" "$REF" "$R1" "$R2" | \
  samtools view -@ 16 -bS - | \
  samtools sort -@ 16 -o "${BAM_PREFIX}.sorted.bam" -

echo "Finished BWA mapping and coordinate sort: $(date)"
ls -lh "${BAM_PREFIX}.sorted.bam"

# check sorted BAM
echo "Checking sorted BAM: $(date)"
samtools quickcheck -v "${BAM_PREFIX}.sorted.bam"

# index sorted BAM with CSI index for large chromosome-scale scaffolds
echo "Indexing sorted BAM with CSI index: $(date)"
samtools index -@ 16 -c "${BAM_PREFIX}.sorted.bam"

# name-sort for fixmate
echo "Starting name sort: $(date)"
samtools sort -@ 16 -n \
  -o "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.sorted.bam"

echo "Finished name sort: $(date)"
ls -lh "${BAM_PREFIX}.namesort.bam"

# add mate information
echo "Starting fixmate: $(date)"
samtools fixmate -@ 16 -m \
  "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

echo "Finished fixmate: $(date)"
ls -lh "${BAM_PREFIX}.fixmate.bam"

# coordinate-sort again
echo "Starting coordinate sort after fixmate: $(date)"
samtools sort -@ 16 \
  -o "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

echo "Finished coordinate sort after fixmate: $(date)"
ls -lh "${BAM_PREFIX}.fixmate.coordsort.bam"

# mark duplicates
echo "Starting markdup: $(date)"
samtools markdup -@ 16 \
  "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.markdup.bam"

echo "Finished markdup: $(date)"
ls -lh "${BAM_PREFIX}.markdup.bam"

# index final duplicate-marked BAM with CSI index
echo "Indexing final markdup BAM with CSI index: $(date)"
samtools index -@ 16 -c "${BAM_PREFIX}.markdup.bam"

# produce basic mapping summary
echo "Running flagstat: $(date)"
samtools flagstat "${BAM_PREFIX}.markdup.bam" > "${BAM_PREFIX}.markdup.flagstat.txt"

echo "Final files:"
ls -lh "${BAM_PREFIX}.markdup.bam"*
ls -lh "${BAM_PREFIX}.markdup.flagstat.txt"

echo "Done: $(date)"