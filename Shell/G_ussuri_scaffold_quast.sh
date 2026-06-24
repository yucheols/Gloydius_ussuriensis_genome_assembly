#!/bin/bash
#SBATCH --job-name=quast_scaffold
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the quast conda environment
source ~/.bash_profile
conda activate quast

# path to assembly
path_to_asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/yahs_out

# output directory
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/yahs_out/final_QC/quast
mkdir -p ${out_dir}

# run quast
quast.py -t ${SLURM_CPUS_PER_TASK} ${path_to_asm}/Gloydius_ussuriensis_AMNH_21010_yahs_scaffolds_final.fa -o ${out_dir} 