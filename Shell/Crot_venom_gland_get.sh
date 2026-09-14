#!/bin/bash
#SBATCH --job-name=Crot_get_venom_data
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
# obtain venom gland RNA-seq data for Crotalus adamanteus from NCBI SRA database
# the specimen code for these RNA-seq data is DRR0105, which matches the genome individual as per Hogan et al. 2024 PNAS paper

# activate conda env
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate sra_tools

# set directories
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/crot_venom_gland"
sradir="${basedir}/ncbi_seq"
fastqdir="${basedir}/fastq"
tmpdir="${basedir}/tmp"

mkdir -p $basedir "$sradir" "$fastqdir" "$tmpdir"

# download SRA files
prefetch SRR12915693 -O "$sradir"   # C. adamanteus DRR0105 right venom gland
prefetch SRR12915694 -O "$sradir"   # C. adamanteus DRR0105 left venom gland

# convert .sra files to FASTQ files
fasterq-dump "$sradir/SRR12915693/SRR12915693.sra" \
  --split-files \
  --threads 12 \
  --temp "$tmpdir" \
  -O "$fastqdir"

fasterq-dump "$sradir/SRR12915694/SRR12915694.sra" \
  --split-files \
  --threads 12 \
  --temp "$tmpdir" \
  -O "$fastqdir"

# gzip FASTQ files
gzip "$fastqdir"/*.fastq