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


# activate conda anv
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# avoid system libstdc++ conflict with SignalP / PIL
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# create temporary directory 
export TMPDIR="/tmp/yshin_funannotate_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT


# existing funannotate output
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"

# verify InterProScan output
iprxml="${outdir}/annotate_misc/iprscan.xml"

if [ ! -s "$iprxml" ]; then
    echo "ERROR: InterProScan XML is missing or empty:"
    echo "$iprxml"
    exit 1
fi

echo
echo "InterProScan XML found:"
ls -lh "$iprxml"
echo

# verify SignalP installation/model
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

# run funannotate annotate
funannotate annotate \
    -i "$outdir" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --cpus 32 \
    --busco_db tetrapoda \
    --database "$FUNANNOTATE_DB" \
    --iprscan "$iprxml" \
    --tmpdir "$TMPDIR"