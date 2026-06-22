#!/bin/bash
#SBATCH --job-name=fastp_HiC
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# path to HiC reads
path_to_seq=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined

# run fastp
fastp \
  -i ${path_to_seq}/Gloydius_ussuriensis_HiC_R1.fastq.gz \
  -I ${path_to_seq}/Gloydius_ussuriensis_HiC_R2.fastq.gz \
  -o ${path_to_seq}/Gloydius_ussuriensis_HiC_R1_trimmed.fastq.gz \
  -O ${path_to_seq}/Gloydius_ussuriensis_HiC_R2_trimmed.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 30 \
  --thread 16 \
  --html Gloydius_ussuriensis_HiC_fastp.html \
  --json Gloydius_ussuriensis_HiC_fastp.json
