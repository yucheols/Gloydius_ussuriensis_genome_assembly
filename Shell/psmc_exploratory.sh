#!/bin/bash
#SBATCH --job-name=psmc_exploratory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC/slurm_logs/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC/slurm_logs/slurm-%x_%j.err

# ================================================================
# PSMC demographic inference
#
# individual: AMNH_21010
#
# input:
#   Repeat-masked diploid consensus
#   Autosomes only
#   Illumina depth filter = 4-24x
#
# initial parameterization:
#   -N25
#   -t15
#   -r5
#   -p "4+25*2+4+6"
#
# NOTE:
# This is the initial exploratory fit. The interval pattern will
# be evaluated after the run before bootstrapping.
# ================================================================


# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate psmc

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"
INDIR="${PSMCROOT}/03_psmc"
INPUT="${INDIR}/AMNH_21010.autosomes.dp4-24.repeatmasked.psmcfa"
OUTPUT="${INDIR}/AMNH_21010.autosomes.dp4-24.repeatmasked.psmc"


# ------------------------------------------------------------
# check input
# ------------------------------------------------------------
if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: missing or empty PSMCFA:"
    echo "${INPUT}"
    exit 1
fi


# ------------------------------------------------------------
# run PSMC
# ------------------------------------------------------------
echo "============================================================"
echo "Running PSMC"
echo "============================================================"
echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
echo

psmc \
    -N25 \
    -t15 \
    -r5 \
    -p "4+25*2+4+6" \
    -o "${OUTPUT}" \
    "${INPUT}"


# ------------------------------------------------------------
# validate output
# ------------------------------------------------------------
if [[ ! -s "${OUTPUT}" ]]; then
    echo "ERROR: PSMC output is empty"
    exit 1
fi

echo
echo "============================================================"
echo "PSMC completed"
echo "============================================================"

ls -lh "${OUTPUT}"

echo
echo "Final iteration:"
grep '^RD' "${OUTPUT}" | tail -1

echo
echo "Final theta/rho:"
grep '^TR' "${OUTPUT}" | tail -1

echo
echo "Final inferred recombination summary:"
grep 'n_recomb' "${OUTPUT}" | tail -1

echo "============================================================"