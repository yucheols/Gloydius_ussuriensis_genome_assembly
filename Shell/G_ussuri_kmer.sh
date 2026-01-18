#!/bin/bash
#SBATCH --job-name=kmer_ussuri
#SBATCH --nodes=1
#SBATCH --mem=650G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=48
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%j_%x.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%j_%x.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate mytools

# create a temp directory
TMPDIR=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/tmp_$SLURM_JOB_ID
mkdir -p $TMPDIR

# decompress in tempdir
zcat /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz > $TMPDIR/AMNH_21010_HiFi.fastq

# specify output directory
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/

# run jellyfish
jellyfish count -m 21 -s 1G -o ${outdir}/Gloydius_ussuriensis_kmer.jf $TMPDIR/AMNH_21010_HiFi.fastq
jellyfish histo -t ${SLURM_CPUS_PER_TASK} ${outdir}/Gloydius_ussuriensis_kmer.jf > ${outdir}/Gloydius_ussuriensis_kmer.histo