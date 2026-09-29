#!/bin/bash
#SBATCH --job-name=psmc_primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ================================================================
# PSMC demographic inference
#
# individual: AMNH_21010
#
# input:
#   repeat-masked diploid consensus
#   autosomes only
#   Illumina depth filter = 4-24x
#
# parameterization:
#   -N100
#   -t15
#   -r5
#   -p "4+25*2+4+6"
#
# NOTE: this is aclean restart from the default PSMC 
# initialization with up to 100 iterations
#
# the previous PSMC output, if present, 
# is removed before the new run
#
# evaluate convergence before bootstrapping.
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
# remove previous run
# ------------------------------------------------------------
echo "============================================================"
echo "Removing previous PSMC output"
echo "============================================================"

if [[ -e "${OUTPUT}" ]]; then
    echo "Removing:"
    echo "${OUTPUT}"
    rm -f "${OUTPUT}"
else
    echo "No previous output found."
fi

echo


# ------------------------------------------------------------
# run PSMC
# ------------------------------------------------------------
echo "============================================================"
echo "Running PSMC"
echo "============================================================"
echo "Input:          ${INPUT}"
echo "Output:         ${OUTPUT}"
echo "Iterations:     100"
echo "Initial t:      15"
echo "Initial rho/theta: 5"
echo "Pattern:        4+25*2+4+6"
echo "============================================================"
echo

psmc \
    -N100 \
    -t15 \
    -r5 \
    -p "4+25*2+4+6" \
    -o "${OUTPUT}" \
    "${INPUT}"


# ------------------------------------------------------------
# report final results
# ------------------------------------------------------------
echo
echo "============================================================"
echo "PSMC completed"
echo "============================================================"

ls -lh "${OUTPUT}"

echo
echo "Final iteration:"
grep '^RD' "${OUTPUT}" | tail -1

echo
echo "Final log likelihood:"
grep '^LK' "${OUTPUT}" | tail -1

echo
echo "Final QD:"
grep '^QD' "${OUTPUT}" | tail -1

echo
echo "Final theta/rho:"
grep '^TR' "${OUTPUT}" | tail -1

echo
echo "Final inferred recombination summary:"
grep 'n_recomb' "${OUTPUT}" | tail -1

echo
echo "Final parameter set:"
grep '^PA' "${OUTPUT}" | tail -1

echo
echo "============================================================"
echo "PSMC run completed successfully"
echo "============================================================"