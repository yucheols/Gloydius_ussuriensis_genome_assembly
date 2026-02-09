#!/bin/bash
#SBATCH --job-name=blobView_ussuri
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
conda activate blobtools

# set paths to different inputs
asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa
aln_sorted_bam
megablast
nodes=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/data/nodes.dmp
names=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/data/names.dmp

# output dir
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/blob_out

# run blobtools
blobtools create -i ${asm} -b ${aln_sorted_bam} -t ${megablast} --nodes ${nodes} --names ${names} \
  -o ${outdir}/ussuri_blob_results

# view results and plot
#blobtools view -i result.blobDB.json
#blobtools plot -i result.blobDB.json