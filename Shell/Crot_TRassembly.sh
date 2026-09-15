#!/bin/bash
#SBATCH --job-name=Crot_TRassembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

##### activate conda environment
# TRassembly.py is a part of ToxCodAn-Genome
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate ToxcodanGenome
set -euo pipefail

##### force Python UTF-8 output
export LANG=C
export LC_ALL=C
export PYTHONIOENCODING=UTF-8
export PYTHONUTF8=1


##### set directories
dir_TRassembly="/home/yshin/mendel-nas1/ToxCodAn-Genome/bin"
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/crot_venom_gland"
indir="${basedir}/trimmed_fastq"
pooldir="${basedir}/pooled_trimmed_fastq"
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/crot_toxin_gene_annotation"
assembly_out="${outdir}/SRR12915693_SRR12915694_TRassembly"

mkdir -p \
    "${pooldir}" \
    "${outdir}"


##### input genome
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny/Crotalus_adamanteus/Crotalus_adamanteus.genome.fa"

##### input venom-gland RNA-seq reads
SRR12915693_R1="${indir}/SRR12915693_1.trimmed.fastq.gz"
SRR12915693_R2="${indir}/SRR12915693_2.trimmed.fastq.gz"

SRR12915694_R1="${indir}/SRR12915694_1.trimmed.fastq.gz"
SRR12915694_R2="${indir}/SRR12915694_2.trimmed.fastq.gz"


##### pool left and right venom-gland RNA-seq libraries
#
# SRR12915693 and SRR12915694 are the two venom glands
# sampled from the same individual.
#
# _1 = sequencing mate R1
# _2 = sequencing mate R2

pooled_R1="${pooldir}/Crot_venom_gland_pooled_R1.fastq.gz"
pooled_R2="${pooldir}/Crot_venom_gland_pooled_R2.fastq.gz"

echo "Pooling R1 reads..."

zcat \
    "${SRR12915693_R1}" \
    "${SRR12915694_R1}" \
    | gzip -c \
    > "${pooled_R1}"

echo "Pooling R2 reads..."

zcat \
    "${SRR12915693_R2}" \
    "${SRR12915694_R2}" \
    | gzip -c \
    > "${pooled_R2}"


##### verify pooled files
ls -lh "${pooled_R1}" "${pooled_R2}"


##### run TRassembly.py
#
# default method = both
#
# this performs:
#   genome-guided assembly
#   de novo assembly

cd "${dir_TRassembly}"

python TRassembly.py \
    -g "${genome}" \
    -r "${pooled_R1},${pooled_R2}" \
    -o "${assembly_out}" \
    -m both \
    -c "${SLURM_CPUS_PER_TASK}" \
    -M 250G


##### completion message
echo
echo "TRassembly completed successfully."
echo "The output files are located in:"
echo "${assembly_out}"
echo

ls -lh "${assembly_out}"