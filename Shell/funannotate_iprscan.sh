#!/bin/bash
#SBATCH --job-name=funannotate_iprscan
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=720:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# ------------------------------------------------------------
# set paths
# ------------------------------------------------------------

FUN_DIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
IPR_DIR="/home/yshin/mendel-nas1/interproscan/interproscan-5.78-109.0"
IPRSCAN="${IPR_DIR}/interproscan.sh"


# ------------------------------------------------------------
# run InterProScan through funannotate
# ------------------------------------------------------------

cd "${FUN_DIR}"

echo "Starting InterProScan: $(date)"

funannotate iprscan \
    -i "${FUN_DIR}" \
    -m local \
    --iprscan_path "${IPRSCAN}" \
    -c 4


# ------------------------------------------------------------
# verify output
# ------------------------------------------------------------

IPR_XML="${FUN_DIR}/annotate_misc/iprscan.xml"

echo
echo "============================================================"

if [[ -s "${IPR_XML}" ]]; then
    echo "InterProScan completed successfully."
    echo "Output:"
    ls -lh "${IPR_XML}"
else
    echo "ERROR: Expected InterProScan XML was not produced:"
    echo "${IPR_XML}"
    exit 1
fi


# ------------------------------------------------------------
# print when finished
# ------------------------------------------------------------

echo "Finished: $(date)"
echo "============================================================"