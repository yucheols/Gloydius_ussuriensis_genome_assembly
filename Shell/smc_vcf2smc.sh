#!/bin/bash
#SBATCH --job-name=smc_vcf2smc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ================================================================
# convert filtered autosomal VCF to SMC++ input
#
# one autosome per SLURM array task.
#
# population: Mainland G. ussuriensis
#
# samples: 8 mainland individuals
# 
# distinguished individual (-d): AMNH_21172
#    this sample has excellent coverage and the lowest final missingness
#
# because the VCF is unphased, both distinguished lineages
# are taken from the same individual:
#
#   -d AMNH_21172 AMNH_21172
#
# the depth/callability exclusion mask is supplied with --mask.
# ================================================================

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
SIF="${WORKDIR}/00_env/smcpp_latest.sif"
VCF="${WORKDIR}/04_vcf/final/G_ussuriensis_mainland.autosomes.filtered.vcf.gz"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
MASK="${WORKDIR}/05_masks/smcpp_exclude.primary.autosomes.sorted.bed"
MASKDIR="${WORKDIR}/05_masks/by_chrom"
OUTDIR="${WORKDIR}/06_smc/primary"

mkdir -p "$MASKDIR"
mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# samples
# ------------------------------------------------------------
POP="Mainland"
SAMPLES="AMNH_21010,AMNH_21128,AMNH_21147,AMNH_21161,AMNH_21162,AMNH_21164,AMNH_21172,AMNH_21185"
DIST="AMNH_21172"


# ------------------------------------------------------------
# load Apptainer
# ------------------------------------------------------------
module load Apptainer/apptainer-1.2.5


# ------------------------------------------------------------
# get chromosome
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")


# ------------------------------------------------------------
# chromosome-specific mask
#
# SMC++ processes one chromosome at a time, so extract only
# the mask intervals belonging to the current chromosome.
# ------------------------------------------------------------
CHRMASK="${MASKDIR}/${CHR}.smcpp_exclude.bed"

awk -v chr="$CHR" '
BEGIN {OFS="\t"}
$1==chr {
    print $1,$2,$3
}
' "$MASK" > "$CHRMASK"


# ------------------------------------------------------------
# output
# ------------------------------------------------------------
OUT="${OUTDIR}/${CHR}.${DIST}.smc.gz"


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ vcf2smc"
echo "============================================================"
echo "Chromosome:             $CHR"
echo "Population:             $POP"
echo "Distinguished sample:   $DIST"
echo "VCF:                    $VCF"
echo "Mask:                   $CHRMASK"
echo "Container:              $SIF"
echo "Output:                 $OUT"
echo "Host:                   $(hostname)"
echo "Array task:             ${SLURM_ARRAY_TASK_ID}"
echo "Start:                  $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# report chromosome mask size
# ------------------------------------------------------------
MASKED_BP=$(
    awk '
    {
        bp += $3-$2
    }
    END {
        print bp+0
    }' "$CHRMASK"
)

echo "Masked bp on ${CHR}: $MASKED_BP"
echo


# ------------------------------------------------------------
# verify container
# ------------------------------------------------------------
echo "Checking SMC++ container..."

apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$SIF" \
    smc++ --help \
    >/dev/null


# ------------------------------------------------------------
# convert VCF to SMC++ format
#
# -d DIST DIST:
#   for unphased data, both distinguished lineages are drawn
#   from the same diploid individual.
#
# --mask:
#   positions in this BED are marked missing rather than being
#   interpreted as long homozygous sequence.
# ------------------------------------------------------------
apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$SIF" \
    smc++ vcf2smc \
    -d "$DIST" "$DIST" \
    --mask "$CHRMASK" \
    "$VCF" \
    "$OUT" \
    "$CHR" \
    "${POP}:${SAMPLES}"


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Completed chromosome:   $CHR"
echo "Distinguished sample:   $DIST"
echo "Output:                 $OUT"
echo "Output size:            $(du -h "$OUT" | cut -f1)"
echo "Finish:                 $(date)"
echo "============================================================"