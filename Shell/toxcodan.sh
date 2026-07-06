#!/bin/bash
#SBATCH --job-name=ToxCodAn-Genome
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env for mapping
source ~/.bash_profile
conda activate ToxcodanGenome

# force UTF-8 output in non-interactive SLURM jobs
export LANG=C.UTF-8
export LC_ALL=C.UTF-8
export PYTHONIOENCODING=UTF-8
export PYTHONUTF8=1

set -euo pipefail

# toxcodan-genome path
dir_toxcodan="/home/yshin/mendel-nas1/ToxCodAn-Genome/bin"
cd ${dir_toxcodan}

# set input paths
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
venom_read_1="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq/SRR35908235_1.trimmed.fastq.gz"
venom_read_2="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq/SRR35908235_2.trimmed.fastq.gz"
db_dir="/home/yshin/mendel-nas1/ToxCodAn-Genome/Databases/Viperidae_db.fasta"

# output dir
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation"
mkdir -p ${outdir}

# run ToxCodAn-Genome
python toxcodan-genome.py \
    -g ${genome} \
    -d ${db_dir} \
    -r ${venom_read_1},${venom_read_2} \
    -o ${outdir} \
    -c ${SLURM_CPUS_PER_TASK}
