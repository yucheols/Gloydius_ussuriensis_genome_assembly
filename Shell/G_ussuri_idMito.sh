#!/bin/bash
#SBATCH --job-name=idMito_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# load blast as module
module load NCBI/blast-2.10.1+

# path to assembly and mito
mito_ref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/NC_026553.1.fa
asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/blast_db/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa

# output directory
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# run blast to id mito
blastn -query ${mito_ref} -db ${asm} -outfmt '6 sseqid pident length bitscore evalue' \
  | awk '$2 >= 1000' \
  | sort -k2,2nr \
  | cut -f1 \
  | uniq \
  > ${outdir}/mito_contigs.txt
