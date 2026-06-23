#!/bin/bash
#SBATCH --job-name=hic_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=72:00:00
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
tmpdir="${outdir}/tmp_hic_map_${SLURM_JOB_ID}"

# make directories for temp files
sort1_tmp="${tmpdir}/sort1"
namesort_tmp="${tmpdir}/namesort"
coordsort_tmp="${tmpdir}/coordsort"

mkdir -p "${sort1_tmp}" "${namesort_tmp}" "${coordsort_tmp}"

# make a directory for bwa logs
mkdir -p "${outdir}/logs"

# set variables
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
R1="${indir}/Gloydius_ussuriensis_HiC_R1.fastq.gz"
R2="${indir}/Gloydius_ussuriensis_HiC_R2.fastq.gz"
BAM_PREFIX="${outdir}/Gloydius_ussuriensis_HiC_to_draft"

# basic checks
echo "Starting BWA mapping: $(date)"
echo "REF = ${REF}"
echo "R1  = ${R1}"
echo "R2  = ${R2}"
echo "TMP = ${tmpdir}"

# quota/filesystem reporting
quota -s || true
df -h "${outdir}" || true

# flag if inputs are missing
[[ -s "${REF}" ]] || { echo "ERROR: missing reference ${REF}"; exit 1; }
[[ -s "${R1}" ]] || { echo "ERROR: missing R1 ${R1}"; exit 1; }
[[ -s "${R2}" ]] || { echo "ERROR: missing R2 ${R2}"; exit 1; }

# run bwa mem
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} "${REF}" "${R1}" "${R2}" 2> "${outdir}/logs/bwa_mem.log" | \
  samtools view -@ 16 -bS - | \
  samtools sort -@ 16 -m 4G -T "${sort1_tmp}/sort1" -o "${BAM_PREFIX}.sorted.bam" -

samtools index "${BAM_PREFIX}.sorted.bam"

# name-sort for fixmate
echo "Starting name-sort: $(date)"
ls -ld "${namesort_tmp}"
df -h "${namesort_tmp}" || true
quota -s || true

samtools sort -@ 16 -m 4G -T "${namesort_tmp}/namesort" -n \
  -o "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.sorted.bam"

# add mate information
samtools fixmate -@ 16 -m \
  "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

# coordinate-sort again
samtools sort -@ 16 -m 4G -T "${coordsort_tmp}/coordsort" \
  -o "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

# mark duplicates
samtools markdup -@ 16 \
  "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.markdup.bam"

samtools index "${BAM_PREFIX}.markdup.bam"

# produce basic mapping summary
samtools flagstat "${BAM_PREFIX}.markdup.bam" > "${BAM_PREFIX}.markdup.flagstat.txt"

# check final BAM integrity
samtools quickcheck -v "${BAM_PREFIX}.markdup.bam"

# clean temp files only if everything succeeded
rm -rf "${tmpdir}"

echo "Finished Hi-C mapping: $(date)"
ls -lh "${BAM_PREFIX}.markdup.bam" "${BAM_PREFIX}.markdup.bam.bai" "${BAM_PREFIX}.markdup.flagstat.txt"