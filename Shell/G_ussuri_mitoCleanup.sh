#!/bin/bash
#SBATCH --job-name=mitoCleanup_ussuri
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

# make the run fail loudly if something is broken
set -euo pipefail

# set paths as variables
wkdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_cleanup
path_to_hifi=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ
path_to_mitoref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref
path_to_assembly=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_cleanup

# cd into working directory
cd ${wkdir}


####################################
#  1. identify likely mito contig  #
####################################

# print this message when starting
echo "ID likely mito contig..."

# list contigs by size
seqkit fx2tab -n -l ${path_to_assembly}/assembly.fasta | sort -k2,2n > contig_len.txt
cat contig_len.txt

# pick contig that is within the exptected size range of snake mitogenome
mito_contig=$(awk '$2>=16000 && $2<=19000 {print $1}' contig_len.txt)

if [[ -z $mito_contig ]]; then
  echo "ERROR! No contig in expected mitogenome size range"
  exit 1
fi

echo "selected likely mito contig: ${mito_contig}, moving on..."

# extract the selected contig
seqkit grep -n -p ${mito_contig} ${path_to_assembly}/assembly.fasta > mito_extract.fasta
seqkit stats mito_extract.fasta

# validate against the reference
minimap2 -x asm5 ${path_to_mitoref}/NC_026553.1.fa mito_extract.fasta > validate.paf
head validate.paf


###########################
#  2. remap & reassemble  #
###########################

# print this message when starting
echo "remap and reassemble..."

# mkdir for flye output
mkdir -p flye_mito_2

# remap raw hifi reads to extracted mito contig
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi mito_extract.fasta ${path_to_hifi}/AMNH_21010_HiFi.fastq.gz | \
  samtools view -b -F 4 -q 60 -o mito_remap.bam

samtools fastq mito_remap.bam | gzip > mito_remap.fastq.gz

# reassemble
flye --pacbio-hifi mito_remap.fastq.gz \
  --out-dir flye_mito_2 \
  --genome-size 20k \
  -t ${SLURM_CPUS_PER_TASK}

# print this when flye run is complete
echo "reassembly stats"
seqkit stats flye_mito_2/assembly.fasta

### print this message when all runs are complete
echo "remap and reassembly completed"

