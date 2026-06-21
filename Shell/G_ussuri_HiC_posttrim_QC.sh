#!/bin/bash
#SBATCH --job-name=fastqc_hic
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# cd into working directory
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined
mkdir -p posttrim_qc


fastqc \
  Gloydius_ussuriensis_HiC_R1_trimmed.fastq.gz \
  Gloydius_ussuriensis_HiC_R2_trimmed.fastq.gz \
  -o posttrim_qc \
  -t 8

# activate multiqc conda env 
conda activate multiqc

# run multiqc
multiqc -o posttrim_qc --filename "HiC_posttrim_QC" posttrim_qc