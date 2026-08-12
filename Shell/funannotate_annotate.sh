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


# ============================================================
# 1. activate Funannotate environment
# ============================================================

source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# avoid system libstdc++ conflict with SignalP/PIL
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"


# ============================================================
# 2. temporary directory
# ============================================================

export TMPDIR="/tmp/yshin_funannotate_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT


# ============================================================
# 3. existing Funannotate predict output
# ============================================================

outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate/G_ussuriensis_funannotate"

if [ ! -d "$outdir" ]; then
    echo "ERROR: Funannotate output directory does not exist:"
    echo "$outdir"
    exit 1
fi

if [ ! -d "$outdir/predict_results" ]; then
    echo "ERROR: predict_results directory not found:"
    echo "$outdir/predict_results"
    exit 1
fi

echo
echo "Using completed Funannotate predict output:"
echo "$outdir"
echo


# ============================================================
# 4. verify SignalP installation and model
# ============================================================

echo "Checking SignalP..."

command -v signalp6

SIGNALP_DIR="${CONDA_PREFIX}/lib/python3.8/site-packages/signalp"
SIGNALP_MODEL="${SIGNALP_DIR}/model_weights/distilled_model_signalp6.pt"

if [ ! -s "$SIGNALP_MODEL" ]; then
    echo
    echo "ERROR: SignalP fast model is missing:"
    echo "$SIGNALP_MODEL"
    exit 1
fi

echo
echo "SignalP model found:"
ls -lh "$SIGNALP_MODEL"
echo


# ============================================================
# 5. test SignalP imports
# ============================================================

python -c 'from PIL import Image; print("PIL OK")'
python -c 'import matplotlib; print("matplotlib OK")'
python -c 'import signalp; print("SignalP import OK")'

echo


# ============================================================
# 6. remove only incomplete SignalP output from failed run
# ============================================================

signalp_old="${outdir}/annotate_misc/signalp"

if [ -d "$signalp_old" ]; then
    echo "Removing incomplete SignalP output from previous failed run:"
    echo "$signalp_old"

    rm -rf "$signalp_old"
fi

echo


# ============================================================
# 7. resume Funannotate annotate
# ============================================================

echo "Running funannotate annotate..."
echo

funannotate annotate \
    -i "$outdir" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --cpus 32 \
    --busco_db tetrapoda \
    --database "$FUNANNOTATE_DB" \
    --tmpdir "$TMPDIR"

echo
echo "Funannotate annotate finished."
echo


# ============================================================
# 8. check final results
# ============================================================

results="${outdir}/annotate_results"

if [ ! -d "$results" ]; then
    echo "ERROR: annotate_results directory was not created."
    exit 1
fi

echo "Final annotation results:"
echo "$results"
echo

ls -lh "$results"

echo
echo "Done."