#!/bin/bash
#SBATCH --job-name=funannotate_predict
#SBATCH --nodes=1
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# funannotate predict
#
# GlimmerHMM is disabled because the previous run produced
# 431,192 GlimmerHMM predictions, causing an excessively large
# EVM input and contributing to an OOM failure.
# ============================================================

# activate environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# avoid system library conflicts
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="${GENEMARK_PATH}:${PATH}"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# temporary directory
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT


### set paths for input and output
# input genome
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# output directory
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"

# run funannotate predict
echo
echo "============================================================"
echo "Starting funannotate predict"
echo "GlimmerHMM weight = 0 (disabled)"
echo "============================================================"
echo

funannotate predict \
    -i "${genome}" \
    -o "${outdir}" \
    --species "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --busco_db tetrapoda \
    --organism other \
    --busco_seed_species Taeniopygia_guttata \
    --max_intronlen 100000 \
    --repeats2evm \
    --weights glimmerhmm:0 \
    --cpus "${SLURM_CPUS_PER_TASK}"

# print when finished
echo
echo "============================================================"
echo "funannotate predict completed successfully"
echo "End time: $(date)"
echo "============================================================"