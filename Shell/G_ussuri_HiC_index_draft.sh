#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# set directory
indir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/draft

# run bwa
bwa index ${indir}/Gloydius_ussuriensis_AMNH_21010_noMito.fa
samtools faidx ${indir}/Gloydius_ussuriensis_AMNH_21010_noMito.fa