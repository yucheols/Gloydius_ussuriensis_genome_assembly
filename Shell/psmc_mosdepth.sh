#!/bin/bash
#SBATCH --job-name=psmc_mosdepth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate psmc

set -euo pipefail

# setup paths
PSMCROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC"

BAM="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam/AMNH_21010.markdup.bam"
AUTOSOMES="${PSMCROOT}/00_setup/autosomes.bed"
OUTDIR="${PSMCROOT}/01_depth_qc"

# create output directory and cd into it
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# run mosdepth to calculate depth of coverage
# -n          do not write gigantic per-base output
# -t 8        use 8 BAM decompression threads
# -Q 30       require MAPQ ≥ 30
# -F 3844     exclude unmapped, secondary, QC-fail, duplicate, and supplementary alignments
# --by BED    calculate statistics specifically for our 17 autosomes

mosdepth \
    -n \
    -t 8 \
    -Q 30 \
    -F 3844 \
    --by "${AUTOSOMES}" \
    AMNH_21010.autosomes \
    "${BAM}"