#!/bin/bash
#SBATCH --job-name=venom_gland_qc
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

set -euo pipefail

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# set directory
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland"
cd "$basedir"

# make output dir
mkdir -p qc_fastqc

# run fastqc
fastqc fastq/*.fastq.gz \
  -o qc_fastqc \
  -t 12

# activate multiqc conda env 
conda activate multiqc

# run multiqc
multiqc qc_fastqc \
  -o qc_fastqc