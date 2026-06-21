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

# set directory
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
R1="Gloydius_ussuriensis_HiC_R1_trimmed.fastq.gz"
R2="Gloydius_ussuriensis_HiC_R2_trimmed.fastq.gz"
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding"

# run bwa mem
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} ${REF} ${R1} ${R2} | \
  samtools view -@ 16 -bS - | \
  samtools sort -@ 16 -o ${outdir}/Gloydius_ussuriensis_HiC_to_draft.sorted.bam

samtools index ${outdir}/Gloydius_ussuriensis_HiC_to_draft.sorted.bam