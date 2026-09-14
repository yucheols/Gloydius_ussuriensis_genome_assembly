#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# GENESPACE macrosynteny
# 11 snake genomes
# ============================================================

# activate conda environment
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate genespace

set -euo pipefail

# set path
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"
cd "${GS}"

# print helpful information
echo "============================================================"
echo "GENESPACE genome macrosynteny"
echo "============================================================"
echo "Date:   $(date)"
echo "Node:   $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo


# ------------------------------------------------------------
# software
# ------------------------------------------------------------

echo "===== SOFTWARE ====="

echo "Rscript:"
command -v Rscript

echo "OrthoFinder:"
command -v orthofinder

echo "DIAMOND:"
command -v diamond

echo "MCScanX_h:"
command -v MCScanX_h

echo

Rscript -e '
cat(
    "GENESPACE ",
    as.character(packageVersion("GENESPACE")),
    "\n"
)
'


# ------------------------------------------------------------
# input summary
# ------------------------------------------------------------

echo
echo "===== INPUT FILES ====="

for sp in \
    Argyrophis_diardii \
    Xenopeltis_unicolor \
    Candoia_aspera \
    Elaphe_schrenckii \
    Naja_naja \
    Cerastes_gasperettii \
    Vipera_berus \
    Bothrops_insularis \
    Crotalus_adamanteus \
    Gloydius_shedaoensis \
    Gloydius_ussuriensis
do

    bed="${GS}/bed/${sp}.bed"
    pep="${GS}/peptide/${sp}.fa"

    if [[ ! -s "${bed}" ]]; then
        echo "ERROR: missing BED: ${bed}"
        exit 1
    fi

    if [[ ! -s "${pep}" ]]; then
        echo "ERROR: missing peptide FASTA: ${pep}"
        exit 1
    fi

    nbed=$(wc -l < "${bed}")
    npep=$(grep -c '^>' "${pep}")

    echo -e "${sp}\tBED=${nbed}\tproteins=${npep}"

done


# ------------------------------------------------------------
# run
# ------------------------------------------------------------

echo "============================================================"
echo "STARTING GENESPACE"
echo "============================================================"

Rscript \
    "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/R/genespace.R"


echo
echo "============================================================"
echo "GENESPACE FINISHED"
echo "Date: $(date)"
echo "============================================================"