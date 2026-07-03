#!/bin/bash
#SBATCH --job-name=sex_chr_depth
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate scaffolding

set -euo pipefail

# project dir
project="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/sex_chr_coverage"

# chromosome level assembly
ref="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"

# low coverage whole genome sequencing reads
reads="/home/yshin/mendel-nas1/ussuri_popgen/Genome_seq/04.Trimmed"

# metadata and window files
metadata="${project}/sample_metadata.tsv"
windows="${project}/windows_100kb.bed"

# slurm threads
threads="${SLURM_CPUS_PER_TASK}"

# chooses which sample an array task should process
sample=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n+1 {print $1}' "$metadata")

r1="${reads}/${sample}_R1.trimmed.fq.gz"
r2="${reads}/${sample}_R2.trimmed.fq.gz"

echo "Sample: $sample"
echo "R1: $r1"
echo "R2: $r2"

# spit out errors if read files are missing
if [[ ! -s "$r1" ]]; then
    echo "ERROR: missing R1: $r1"
    exit 1
fi

if [[ ! -s "$r2" ]]; then
    echo "ERROR: missing R2: $r2"
    exit 1
fi

# create output & other directories
mkdir -p "${project}/bam" "${project}/depth" "${project}/tmp"

tmpdir="${project}/tmp/${sample}"
mkdir -p "$tmpdir"

if command -v bwa-mem2 >/dev/null 2>&1; then
    aligner="bwa-mem2 mem"
else
    aligner="bwa mem"
fi

echo "Using aligner: $aligner"

# Name-sort for duplicate marking
$aligner -t "$threads" \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    "$ref" "$r1" "$r2" \
    | samtools sort -n -@ "$threads" -T "${tmpdir}/${sample}.name" \
      -o "${tmpdir}/${sample}.name.bam" -

# Fix mate information
samtools fixmate -m \
    "${tmpdir}/${sample}.name.bam" \
    "${tmpdir}/${sample}.fixmate.bam"

# Coordinate sort
samtools sort -@ "$threads" -T "${tmpdir}/${sample}.coord" \
    -o "${tmpdir}/${sample}.coord.bam" \
    "${tmpdir}/${sample}.fixmate.bam"

# Remove duplicates
samtools markdup -r -@ "$threads" \
    "${tmpdir}/${sample}.coord.bam" \
    "${project}/bam/${sample}.markdup.bam"

samtools index -@ "$threads" "${project}/bam/${sample}.markdup.bam"

# Windowed depth.
# MAPQ 30 keeps mainly confidently mapped reads.
mosdepth -t "$threads" \
    --mapq 30 \
    --by "$windows" \
    "${project}/depth/${sample}" \
    "${project}/bam/${sample}.markdup.bam"

# Basic mapping statistics
samtools flagstat -@ "$threads" \
    "${project}/bam/${sample}.markdup.bam" \
    > "${project}/bam/${sample}.flagstat.txt"

rm -rf "$tmpdir"

echo "Finished sample: $sample"