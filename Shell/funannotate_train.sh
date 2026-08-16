#!/bin/bash
#SBATCH --job-name=funannotate_train
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=120:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate environment and set up analyses
source ~/.bash_profile
conda activate funannotate
set -euo pipefail

# avoid system library conflicts
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# use a short temporary directory and clean it up on exit
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT

# set paths for input and output
# softmakesked genome assembly
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# path to paired-end RNA-seq reads from organ tissues and venom gland
tissue_reads="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads"
venom_reads="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq"

# make output dir for all funannotate results
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate/"
mkdir -p "$outdir"

# Train to align RNA-seq data, run Trinity, and then run PASA
funannotate train -i ${genome} -o ${outdir} \
    --left ${tissue_reads}/AMNH_21010_Ht_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lr_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Ky_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Me_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lg_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Skin_R1_paired.fastq.gz ${venom_reads}/SRR35908235_1.trimmed.fastq.gz \
    --right ${tissue_reads}/AMNH_21010_Ht_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lr_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Ky_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Me_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lg_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Skin_R2_paired.fastq.gz ${venom_reads}/SRR35908235_2.trimmed.fastq.gz \
    --species "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --stranded RF \
    --max_intronlen 100000 \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --memory 250G \
    --no_trimmomatic
