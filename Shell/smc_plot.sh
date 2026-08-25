#!/bin/bash
#SBATCH --job-name=smc_plot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.err

# ================================================================
# plot SMC++ demographic history
#
# population: Mainland G. ussuriensis
#
# mutation rate used during estimate:
#   1.25e-8 mutations/site/generation
#
# generation time used for plotting:
#   3 years/generation
# ================================================================


# ------------------------------------------------------------
# load Apptainer
# ------------------------------------------------------------
module load Apptainer/apptainer-1.2.5


# ------------------------------------------------------------
# strict bash
# ------------------------------------------------------------
set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR_LINK="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
WORKDIR=$(readlink -f "$WORKDIR_LINK")

SIF="${WORKDIR}/00_env/smcpp_latest.sif"
MODEL="${WORKDIR}/07_models/primary/model.final.json"
OUTDIR="${WORKDIR}/08_plots/primary"

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# parameters
# ------------------------------------------------------------
GENERATION_TIME="3"


# ------------------------------------------------------------
# container paths
# ------------------------------------------------------------
MODEL_IN="/work/07_models/primary/model.final.json"
PLOT_IN="/work/08_plots/primary/G_ussuriensis_mainland_SMCpp.pdf"


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ plot"
echo "============================================================"
echo "Population:          Mainland G. ussuriensis"
echo "Generation time:     ${GENERATION_TIME} years"
echo "Model:               $MODEL"
echo "Output directory:    $OUTDIR"
echo "Start:               $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# plot demographic history
#
# -g 3:
#   convert generations to years using 3 years/generation
#
# -c:
#   also produce CSV containing plotted x/y coordinates
# ------------------------------------------------------------
apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:/work" \
    "$SIF" \
    smc++ plot \
    -g "$GENERATION_TIME" \
    -c \
    "$PLOT_IN" \
    "$MODEL_IN"


# ------------------------------------------------------------
# check outputs
# ------------------------------------------------------------
PLOT="${OUTDIR}/G_ussuriensis_mainland_SMCpp.pdf"
CSV="${OUTDIR}/G_ussuriensis_mainland_SMCpp.csv"


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ plotting completed"
echo "Plot:     $PLOT"
echo "CSV:      $CSV"
echo "Finish:   $(date)"
echo "============================================================"