#!/bin/bash
#SBATCH --job-name funannotate_fix
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate environment and set up analyses
source ~/.bash_profile
conda activate funannotate
set -euo pipefail

# avoid system library conflicts
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# use a short temporary directory and clean it up on exit
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT

# output directory for funannotate results
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate/"

# run funannotate fix
funannotate fix -i ${outdir}/update_results/Gloydius_ussuriensis.gbk -t ${outdir}/update_results/Gloydius_ussuriensis.tbl