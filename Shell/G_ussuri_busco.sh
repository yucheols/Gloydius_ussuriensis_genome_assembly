#!/bin/sh
#SBATCH --job-name=yshin_G_ussuri_busco
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100gb
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/busco/out/busco-%j-%x.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/busco/err/busco-%j-%x.err

# conda init

source ~/.bash_profile
conda activate mytools

cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/busco/out

G_ussuriensis_assembly