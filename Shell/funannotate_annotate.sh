#!/bin/bash
#SBATCH --job-name=funannotate_annotate
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ------------------------------------------------------------
# activate conda environment
# ------------------------------------------------------------

source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# ------------------------------------------------------------
# environment
# ------------------------------------------------------------

# avoid system libstdc++ conflict
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Funannotate databases
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"

# EggNOG database
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"


# ------------------------------------------------------------
# temporary directory
# ------------------------------------------------------------

export TMPDIR="/tmp/yshin_funannotate_${SLURM_JOB_ID}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT


# ------------------------------------------------------------
# set paths
# ------------------------------------------------------------

FUN_DIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
IPR_XML="${FUN_DIR}/annotate_misc/iprscan.xml"
PROTEINS="${FUN_DIR}/update_results/Gloydius_ussuriensis_AMNH_21010.proteins.fa"


# ------------------------------------------------------------
# run Funannotate annotate
# ------------------------------------------------------------

cd "${FUN_DIR}"

echo "============================================================"
echo "Starting funannotate annotate"
echo "Date: $(date)"
echo "============================================================"
echo

funannotate annotate \
    -i "${FUN_DIR}" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --busco_db tetrapoda \
    --database "${FUNANNOTATE_DB}" \
    --tmpdir "${TMPDIR}"


# ------------------------------------------------------------
# verify annotation output
# ------------------------------------------------------------

ANNOTATE_DIR="${FUN_DIR}/annotate_results"

echo
echo "============================================================"
echo "Checking Funannotate annotation output"
echo "============================================================"
echo

if [[ ! -d "${ANNOTATE_DIR}" ]]; then
    echo "ERROR: annotate_results directory was not created:"
    echo "${ANNOTATE_DIR}"
    exit 1
fi

echo "annotate_results created successfully:"
ls -lh "${ANNOTATE_DIR}"

echo
echo "============================================================"
echo "Funannotate annotate completed"
echo "Finished: $(date)"
echo "============================================================"