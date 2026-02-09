#!/bin/bash
#SBATCH --job-name=blobMapping_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate genome_assembly

# designate paths to assembly and fastq
path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# mapping
# -@ 8 means "use 8 threads" 
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${path}/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa \
  ${path}/AMNH_21010_HiFi.fastq.gz | samtools view -@ 8 -b \
  | samtools sort -@ 8 -o ussuri_aln_sorted.bam

# index BAM 
samtools index ussuri_aln_sorted.bam 

