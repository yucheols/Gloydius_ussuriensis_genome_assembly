#!/bin/bash
#SBATCH --job-name=psmc_bootstrap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=96:00:00
#SBATCH --partition=compute
#SBATCH --array=1-100%20
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%A_%a.err

# ================================================================
# PSMC bootstrap
#
# individual:
#   AMNH_21010
#
# bootstrap replicates:
#   100
#
# input:
#   split PSMCFA generated with splitfa
#
# parameterization:
#   -N500
#   -t15
#   -r5
#   -p "4+25*2+4+6"
#
# -b:
#   resample split genomic segments with replacement
#
# each SLURM array task performs one independent bootstrap replicate.
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
BOOTDIR="${PSMCROOT}/04_bootstrap"
INPUT="${BOOTDIR}/AMNH_21010.autosomes.dp4-24.repeatmasked.split.psmcfa"
OUTDIR="${BOOTDIR}/replicates"
QCDIR="${BOOTDIR}/qc"

mkdir -p "${OUTDIR}" "${QCDIR}"


# ------------------------------------------------------------
# replicate ID
# ------------------------------------------------------------
REP=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")
OUTPUT="${OUTDIR}/bootstrap_${REP}.psmc"
QC="${QCDIR}/bootstrap_${REP}.qc.tsv"


# ------------------------------------------------------------
# check software/input
# ------------------------------------------------------------
if ! command -v psmc >/dev/null 2>&1; then
    echo "ERROR: psmc not found in PATH"
    exit 1
fi

if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: split PSMCFA missing or empty:"
    echo "${INPUT}"
    exit 1
fi


# ------------------------------------------------------------
# remove previous replicate if present
# ------------------------------------------------------------
if [[ -e "${OUTPUT}" ]]; then
    echo "Removing previous bootstrap replicate:"
    echo "${OUTPUT}"
    rm -f "${OUTPUT}"
fi

rm -f "${QC}"


# ------------------------------------------------------------
# report
# ------------------------------------------------------------
echo "============================================================"
echo "PSMC bootstrap replicate ${REP}"
echo "============================================================"
echo "Input:      ${INPUT}"
echo "Output:     ${OUTPUT}"
echo "Iterations: 500"
echo "Pattern:    4+25*2+4+6"
echo "============================================================"
echo


# ------------------------------------------------------------
# run bootstrap PSMC
# ------------------------------------------------------------
psmc \
    -N500 \
    -t15 \
    -r5 \
    -b \
    -p "4+25*2+4+6" \
    -o "${OUTPUT}" \
    "${INPUT}"


# ------------------------------------------------------------
# validate output
# ------------------------------------------------------------
if [[ ! -s "${OUTPUT}" ]]; then
    echo "ERROR: bootstrap output is missing or empty"
    exit 1
fi

FINAL_RD=$(grep '^RD' "${OUTPUT}" | tail -1 | awk '{print $2}')

if [[ "${FINAL_RD}" -ne 500 ]]; then
    echo "ERROR: bootstrap replicate did not reach iteration 500"
    echo "Final RD: ${FINAL_RD}"
    exit 1
fi


# ------------------------------------------------------------
# convergence/QC summary
# ------------------------------------------------------------
FINAL_LK=$(grep '^LK' "${OUTPUT}" | tail -1 | awk '{print $2}')

LK_DELTA=$(
    grep '^LK' "${OUTPUT}" \
    | tail -2 \
    | awk '
        NR==1 {prev=$2}
        NR==2 {printf "%.9f", $2-prev}
    '
)

FINAL_THETA=$(grep '^TR' "${OUTPUT}" | tail -1 | awk '{print $2}')
FINAL_RHO=$(grep '^TR' "${OUTPUT}" | tail -1 | awk '{print $3}')
FINAL_MT=$(grep '^MT' "${OUTPUT}" | tail -1 | awk '{print $2}')

FINAL_NRECOMB=$(
    grep 'n_recomb' "${OUTPUT}" \
    | tail -1 \
    | sed -E 's/.*n_recomb: ([0-9.eE+-]+).*/\1/'
)

{
    printf "replicate\tfinal_iteration\tfinal_LK\tLK_delta\ttheta0\trho0\tmax_t\tn_recomb\n"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${REP}" \
        "${FINAL_RD}" \
        "${FINAL_LK}" \
        "${LK_DELTA}" \
        "${FINAL_THETA}" \
        "${FINAL_RHO}" \
        "${FINAL_MT}" \
        "${FINAL_NRECOMB}"
} > "${QC}"


# ------------------------------------------------------------
# final report
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Bootstrap ${REP} completed"
echo "============================================================"

cat "${QC}"

echo
ls -lh "${OUTPUT}"

echo "============================================================"