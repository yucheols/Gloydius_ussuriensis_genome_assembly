#!/bin/bash
#SBATCH --job-name=smc_mosdepth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-8%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

set -euo pipefail

# --------------------------------------------------
# activate conda
# --------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

# --------------------------------------------------
# directories
# --------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"

SAMPLEFILE="${WORKDIR}/mainland_samples.txt"
BAMDIR="${WORKDIR}/02_bams"
AUTOSOMES="${WORKDIR}/01_reference/autosomes.bed"
OUTDIR="${WORKDIR}/03_qc/mosdepth"

mkdir -p "$OUTDIR"

# --------------------------------------------------
# get sample from array index
# --------------------------------------------------
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLEFILE")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BAM="${BAMDIR}/${sample}.markdup.bam"

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found: ${BAM}.bai"
    exit 1
fi

echo "Sample: $sample"
echo "BAM: $BAM"
echo "Autosome BED: $AUTOSOMES"
echo "Starting: $(date)"

# --------------------------------------------------
# autosomal coverage
# --------------------------------------------------
mosdepth \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --mapq 30 \
    --by "$AUTOSOMES" \
    --no-per-base \
    "${OUTDIR}/${sample}" \
    "$BAM"

echo "Finished: $(date)"