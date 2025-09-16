#!/bin/bash
#SBATCH --job-name=yshin_G_ussuri_kmer
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=60gb
#SBATCH --time=50:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/jellyfish/out/jellyfish-%j-%x.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/jellyfish/err/jellyfish-%j-%x.err

source ~/.bash_profile
conda activate mytools

# create a temp directory
TMPDIR=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/tmp_$SLURM_JOB_ID
mkdir -p $TMPDIR

# decompress in tempdir
zcat /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/24GUHW001.hifireads.fastq.gz > $TMPDIR/24GUHW001.hifireads.fastq

# run jellyfish
jellyfish count -m 21 -s 1G -o /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/jellyfish/out/Gloydius_ussuriensis_kmer.jf $TMPDIR/24GUHW001.hifireads.fastq
jellyfish histo -t 30 /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/jellyfish/out/Gloydius_ussuriensis_kmer.jf > /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/jellyfish/out/Gloydius_ussuriensis_kmer.histo