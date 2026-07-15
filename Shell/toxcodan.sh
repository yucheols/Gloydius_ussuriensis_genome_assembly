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

# Trinity memory for ToxCodAn/Trinity
#export TRINITY_MEM=250G

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
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added"
mkdir -p ${outdir}

# Patch TRassembly.py default Trinity memory to increase memory allocation for Trinity assembly step. 
# This is necessary because the default 8G is insufficient for large genomes, and ToxCodAn does not allow passing -M to TRassembly.py directly.
#if grep -q "8G" "${dir_toxcodan}/TRassembly.py"; then
#    cp "${dir_toxcodan}/TRassembly.py" "${dir_toxcodan}/TRassembly.py.bak_before_memory_patch"
#    sed -i "s/default *= *['\"]8G['\"]/default='${TRINITY_MEM}'/g" "${dir_toxcodan}/TRassembly.py"
#    echo "Patched TRassembly.py Trinity memory default to ${TRINITY_MEM}"
#else
#    echo "TRassembly.py does not contain default 8G; checking current memory setting:"
#fi

#grep -n "max_memory\|default=.*G\|-M" "${dir_toxcodan}/TRassembly.py" || true

# run ToxCodAn-Genome
python toxcodan-genome.py \
    -g ${genome} \
    -d ${db_dir} \
    -r ${venom_read_1},${venom_read_2} \
    -o ${outdir} \
    -c ${SLURM_CPUS_PER_TASK}