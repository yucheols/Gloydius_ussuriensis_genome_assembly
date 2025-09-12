#!/bin/sh
#SBATCH --job-name=yshin_hifiasm_G_ussuri
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/out/yshin_hifiasm_G_ussuri_%j-%x.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/err/yshin_hifiasm_G_ussuri_%j-%x.err

#conda init

source ~/.bash_profile
conda activate mytools

hifiasm -o /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/out/Gloydius_ussuriensis_v1.asm -t /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/24GUHW001.hifireads.fastq.gz