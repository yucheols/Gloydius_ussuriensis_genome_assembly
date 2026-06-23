#!/bin/bash
#SBATCH --job-name=get_venom_data
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=25:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

set -euo pipefail

### commands start here ###
# path to SRA toolkit
export PATH=$PWD/sratoolkit.3.4.1-alma_linux64/bin:$PATH

# set directories
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland"
sradir="${basedir}/ncbi_seq"
fastqdir="${basedir}/fastq"
tmpdir="${basedir}/tmp"

mkdir -p "$sradir" "$fastqdir" "$tmpdir"

# download SRA files
prefetch SRR35908235 -O "$sradir"
prefetch SRR35908238 -O "$sradir"

# convert .sra files to FASTQ files
fasterq-dump "$sradir/SRR35908235/SRR35908235.sra" \
  --split-files \
  --threads 12 \
  --temp "$tmpdir" \
  -O "$fastqdir"

fasterq-dump "$sradir/SRR35908238/SRR35908238.sra" \
  --split-files \
  --threads 12 \
  --temp "$tmpdir" \
  -O "$fastqdir"

# gzip FASTQ files
gzip "$fastqdir"/*.fastq