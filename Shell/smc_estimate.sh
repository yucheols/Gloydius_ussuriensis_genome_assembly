#!/bin/bash
#SBATCH --job-name=smc_estimate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.err

# ================================================================
# SMC++ demographic inference
#
# population: Mainland G. ussuriensis
#
# mutation rate:
#   1.25e-8 mutations/site/generation
#
# generation time:
#   3 years/generation
#
# NOTE:
# generation time is NOT used during smc++ estimate.
# It will be supplied later during smc++ plot.
#
# input:
#   17 autosomal .smc.gz files
#
# output:
#   fitted SMC++ demographic model
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

# resolve symlinked NAS path
WORKDIR=$(readlink -f "$WORKDIR_LINK")

SIF="${WORKDIR}/00_env/smcpp_latest.sif"
SMCDIR="${WORKDIR}/06_smc/primary"
OUTDIR="${WORKDIR}/07_models/primary"

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# parameters
# ------------------------------------------------------------
MU="1.25e-8"


# ------------------------------------------------------------
# collect SMC++ input files
# ------------------------------------------------------------
SMCFILES=("${SMCDIR}"/*.smc.gz)
NSMC=${#SMCFILES[@]}

if [[ "$NSMC" -ne 17 ]]; then
    echo "ERROR: Expected 17 SMC files but found $NSMC"
    exit 1
fi


# ------------------------------------------------------------
# convert host paths to container paths
# ------------------------------------------------------------
SMC_IN=()

for f in "${SMCFILES[@]}"; do

    basename_f=$(basename "$f")

    SMC_IN+=("/work/06_smc/primary/${basename_f}")

done

OUT_IN="/work/07_models/primary"


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ estimate"
echo "============================================================"
echo "Population:             Mainland G. ussuriensis"
echo "Mutation rate:          $MU / site / generation"
echo "Generation time:        3 years"
echo "SMC files:              $NSMC"
echo "Resolved workdir:       $WORKDIR"
echo "Output directory:       $OUTDIR"
echo "Container:              $SIF"
echo "Host:                   $(hostname)"
echo "CPUs:                   ${SLURM_CPUS_PER_TASK}"
echo "Start:                  $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# list input files
# ------------------------------------------------------------
printf '%s\n' "${SMCFILES[@]}"


# ------------------------------------------------------------
# verify SMC++ container
# ------------------------------------------------------------
echo
echo "Checking SMC++ container..."

apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:/work" \
    "$SIF" \
    smc++ --help \
    >/dev/null


# ------------------------------------------------------------
# estimate demographic history
# ------------------------------------------------------------
echo
echo "Starting SMC++ estimate..."

apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:/work" \
    "$SIF" \
    smc++ estimate \
    -o "$OUT_IN" \
    "$MU" \
    "${SMC_IN[@]}"


# ------------------------------------------------------------
# verify final model
# ------------------------------------------------------------
MODEL="${OUTDIR}/model.final.json"

if [[ ! -s "$MODEL" ]]; then
    echo "ERROR: model.final.json was not produced"
    exit 1
fi


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ estimate completed"
echo "Mutation rate:          $MU / site / generation"
echo "Model:                  $MODEL"
echo "Model size:             $(du -h "$MODEL" | cut -f1)"
echo "Finish:                 $(date)"
echo "============================================================"