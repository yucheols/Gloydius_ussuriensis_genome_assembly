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

G_ussuriensis_assembly="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/out/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa"
busco -m genome -i $G_ussuriensis_assembly -o G_ussuriensis_BUSCO -l /home/yshin/mendel-nas1/snake_genome_ass/busco/sauropsida_odb10 -f --metaeuk --offline --download_path /home/yshin/mendel-nas1/snake_genome_ass/busco