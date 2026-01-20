#!/bin/bash
#SBATCH --job-name=mito_ussuri
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda environment
source ~/.bash_profile
conda activate mito_assembly

##############
#  1. setup  #
##############

# make the run fail loudly if something is broken
set -euo pipefail

# set some paths as variables
path_to_hifi=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ
path_to_mitoref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_out

# set mapping quality score
# set this to 60 to drop numts and other crap that we do not want
mq=60


################
#  2. mapping  #
################

# print this message when starting
echo "start mapping reads to reference..."

# run minimap2
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${path_to_mitoref}/NC_026553.1.fa \
  ${path_to_hifi}/AMNH_21010_HiFi.fastq.gz | \
  samtools view -@ 8 -b -F 4 | \
  samtools sort -@ 8 -T ${out_dir}/tmp.mito -o ${out_dir}/mito_map.bam
samtools index ${out_dir}/mito_map.bam

#############################
#  3. extract mapped reads  #
#############################

# print this message when starting
echo "extract mapped reads for mapping quality score >= ${mq}"

# use samtools
samtools view -h ${out_dir}/mito_map.bam | \
  awk -v mq=${mq} 'BEGIN{OFS="\t"} $1 ~ /^@/ || $5 >= mq' | \
  samtools view -b -F 4 -o ${out_dir}/mito_map.MQ.bam

samtools index ${out_dir}/mito_map.MQ.bam

# Convert to FASTQ for assembly
samtools fastq ${out_dir}/mito_map.MQ.bam | gzip > ${out_dir}/mito_reads.MQ.fastq.gz

# print out this message
echo "mito reads written as FASTQ: mito_reads.MQ.fastq.gz"


#######################
#  4. assemble reads  #
#######################

# print this message when starting
echo "assemble reads with flye..."

# run flye
flye --pacbio-hifi ${out_dir}/mito_reads.MQ.fastq.gz \
  --out-dir ${out_dir}/flye_mito -t ${SLURM_CPUS_PER_TASK} --genome-size 20k

mito_asm=${out_dir}/flye_mito/assembly.fasta

# print out this message
echo "mitochondrial genome assembly complete! assembly written to: ${mito_asm}"


###############################################################
#  5. map mito reads back to mito assembly and compute depth  #
###############################################################

# print this message when starting
echo "remap mito reads back to assembled contigs..."

# run minimap2
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${mito_asm} ${out_dir}/mito_reads.MQ.fastq.gz | \
  samtools sort -@ 8 -o ${out_dir}/mito_realign.bam
samtools index ${out_dir}/mito_realign.bam

# compute depth
echo "compute depth..."

# run samtools depth
samtools depth ${out_dir}/mito_realign.bam > ${out_dir}/mito_depth.tsv
awk '{sum+=$3; n++} END{print "mean_depth\t"sum/n"\npositions\t"n}' ${out_dir}/mito_depth.tsv > ${out_dir}/depth_summary.txt


###########################################
#  6. validate by comparing to reference  #
###########################################

# print this message when starting
echo "validate assembly by comparing to reference..."

# run minimap2
minimap2 -x asm5 ${path_to_mitoref}/NC_026553.1.fa ${mito_asm} > ${out_dir}/mito_vs_ref.paf


###################################
# print this when all is complete #
###################################
echo "mitochondrial genome assembly pipeline complete!"