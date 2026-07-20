#!/bin/bash
#SBATCH --job-name=ToxCodAn-Genome
#SBATCH --nodes=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env for mapping
source ~/.bash_profile
conda activate ToxcodanGenome

# force UTF-8 output in non-interactive SLURM jobs
# locale available on all standard Linux nodes
export LANG=C
export LC_ALL=C

# force Python itself to use UTF-8
export PYTHONIOENCODING=UTF-8
export PYTHONUTF8=1

set -euo pipefail

# set toxcodan-genome path
dir_toxcodan_genome="/home/yshin/mendel-nas1/ToxCodAn-Genome/bin"

# allow ToxCodAn-Genome to find helper scripts
export PATH="${dir_toxcodan_genome}:${PATH}"

# set input paths
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
db_dir="/home/yshin/mendel-nas1/ToxCodAn-Genome/Databases/Viperidae_db_appended.fasta"

# output dir for all toxin gene annotation steps
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added"
mkdir -p ${outdir}

# output dir for ToxCodAn-Genome
mkdir -p ${outdir}/SRR35908235_ToxCodAn-Genome

# run ToxCodAn-Genome
python ${dir_toxcodan_genome}/toxcodan-genome.py \
    -g ${genome} \
    -d ${db_dir} \
    -C ${outdir}/G_ussuriensis_VG_toxins.toxcodan.fasta \
    -o ${outdir}/SRR35908235_ToxCodAn-Genome \
    -c ${SLURM_CPUS_PER_TASK}

echo "ToxCodAn-Genome completed successfully."
echo "The output files are located in ${outdir}/SRR35908235_ToxCodAn-Genome"
ls ${outdir}/SRR35908235_ToxCodAn-Genome