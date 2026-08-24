#!/bin/bash
#SBATCH --job-name=smc_callability
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-8%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ============================================================
# generate depth-based callability masks for SMC++
#
# coverage classes:
#   0x       = NO_COVERAGE
#   1-4x     = LOW_COVERAGE
#   5-25x    = CALLABLE
#   >=26x    = HIGH_COVERAGE
#
# one sample is processed per SLURM array task.
# ============================================================

# ------------------------------------------------------------
# activate conda environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
SAMPLEFILE="${WORKDIR}/mainland_samples.txt"
BAMDIR="${WORKDIR}/02_bams"
OUTDIR="${WORKDIR}/03_qc/mosdepth_callability"

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# get sample corresponding to SLURM array task
# ------------------------------------------------------------
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLEFILE")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


# ------------------------------------------------------------
# BAM
# ------------------------------------------------------------
BAM="${BAMDIR}/${sample}.markdup.bam"

if [[ ! -e "$BAM" ]]; then
    echo "ERROR: BAM not found:"
    echo "$BAM"
    exit 1
fi

if [[ ! -e "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found:"
    echo "${BAM}.bai"
    exit 1
fi


# ------------------------------------------------------------
# print job information
# ------------------------------------------------------------

echo "============================================================"
echo "SMC++ callability mask"
echo "============================================================"
echo "Sample:           $sample"
echo "BAM:              $BAM"
echo "Output directory: $OUTDIR"
echo "Array task:       ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "Host:             $(hostname)"
echo "Start time:       $(date)"
echo "============================================================"


# ------------------------------------------------------------
# define labels for mosdepth quantization
#
# --quantize 0:1:5:26:
#
# bins:
#   Q0 = 0x
#   Q1 = 1-4x
#   Q2 = 5-25x
#   Q3 = >=26x
# ------------------------------------------------------------
export MOSDEPTH_Q0="NO_COVERAGE"
export MOSDEPTH_Q1="LOW_COVERAGE"
export MOSDEPTH_Q2="CALLABLE"
export MOSDEPTH_Q3="HIGH_COVERAGE"


# ------------------------------------------------------------
# run mosdepth
#
# --mapq 30:
#   ignore reads with mapping quality <30
#
# --flag 1796:
#   exclude unmapped, secondary, QC-failed, and duplicate reads
#
# --no-per-base:
#   do not produce huge per-base BED files
#
# --quantize:
#   merge adjacent bases belonging to the same coverage class
# ------------------------------------------------------------
mosdepth \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --mapq 30 \
    --flag 1796 \
    --no-per-base \
    --quantize 0:1:5:26: \
    "${OUTDIR}/${sample}.callability" \
    "$BAM"


# ------------------------------------------------------------
# check expected output
# ------------------------------------------------------------
QUANT="${OUTDIR}/${sample}.callability.quantized.bed.gz"

if [[ ! -s "$QUANT" ]]; then
    echo "ERROR: Expected quantized BED was not produced:"
    echo "$QUANT"
    exit 1
fi


# ------------------------------------------------------------
# report class counts
# ------------------------------------------------------------
echo
echo "Coverage-class intervals:"
zcat "$QUANT" |
awk '
{
    n[$4]++
    bp[$4] += $3-$2
}
END {
    for (x in n) {
        printf "%-15s intervals=%d bp=%d\n", x, n[x], bp[x]
    }
}'


echo
echo "============================================================"
echo "Completed sample: $sample"
echo "Finish time:      $(date)"
echo "============================================================"