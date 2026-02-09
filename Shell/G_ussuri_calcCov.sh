#!/bin/bash
#SBATCH --job-name=calcCov_ussuri
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=compute
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate samtools

# cd into the working directory
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# calculate genome wide coverage 
samtools coverage mapping/ussuri_aln_sorted.bam \
| awk 'NR>1 {sum += $3*$7; len += $3} END {print "Mean depth =", sum/len "x"}'

# calculate how much of the genome (%) is covered >= 30x
samtools depth -a mapping/ussuri_aln_sorted.bam \
| awk '{cov=$3; total++; if(cov>=30) c30++} END {print "≥30× =", 100*c30/total "%"}'