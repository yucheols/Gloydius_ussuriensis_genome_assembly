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


# ------------------------------------------------------------
# load Apptainer
# ------------------------------------------------------------
module load Apptainer/apptainer-1.2.5


# ------------------------------------------------------------
# activate conda environment for bgzip/tabix
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools


# ------------------------------------------------------------
# strict bash
# ------------------------------------------------------------
set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR_LINK="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"

# resolve symlinked NAS path to its real filesystem location
WORKDIR=$(readlink -f "$WORKDIR_LINK")

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
CHRMASK_GZ="${CHRMASK}.gz"

awk -v chr="$CHR" '
BEGIN {OFS="\t"}
$1==chr {
    print $1,$2,$3
}
' "$MASK" > "$CHRMASK"

bgzip -f "$CHRMASK"
tabix -f -p bed "$CHRMASK_GZ"


# ------------------------------------------------------------
# host output
# ------------------------------------------------------------
OUT="${OUTDIR}/${CHR}.${DIST}.smc.gz"


# ------------------------------------------------------------
# paths as seen INSIDE the container
#
# the real project directory is mounted as /work
# ------------------------------------------------------------
VCF_IN="/work/04_vcf/final/G_ussuriensis_mainland.autosomes.filtered.vcf.gz"
MASK_IN="/work/05_masks/by_chrom/${CHR}.smcpp_exclude.bed.gz"
OUT_IN="/work/06_smc/primary/${CHR}.${DIST}.smc.gz"


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
echo "Resolved workdir:       $WORKDIR"
echo "Host VCF:               $VCF"
echo "Container VCF:          $VCF_IN"
echo "Host mask:              $CHRMASK_GZ"
echo "Container mask:         $MASK_IN"
echo "Container:              $SIF"
echo "Host output:            $OUT"
echo "Container output:       $OUT_IN"
echo "Host:                   $(hostname)"
echo "Array task:             ${SLURM_ARRAY_TASK_ID}"
echo "Start:                  $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# report chromosome mask size
# ------------------------------------------------------------
MASKED_BP=$(
    zcat "$CHRMASK_GZ" |
    awk '
    {
        bp += $3-$2
    }
    END {
        print bp+0
    }'
)

echo "Masked bp on ${CHR}: $MASKED_BP"
echo


# ------------------------------------------------------------
# verify that the container can see the input files
# ------------------------------------------------------------
echo "Checking input files inside container..."

apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:/work" \
    "$SIF" \
    ls -lh \
    "$VCF_IN" \
    "${VCF_IN}.tbi" \
    "$MASK_IN" \
    "${MASK_IN}.tbi"


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
    --bind "${WORKDIR}:/work" \
    "$SIF" \
    smc++ vcf2smc \
    -d "$DIST" "$DIST" \
    --mask "$MASK_IN" \
    "$VCF_IN" \
    "$OUT_IN" \
    "$CHR" \
    "${POP}:${SAMPLES}"


# ------------------------------------------------------------
# verify output
# ------------------------------------------------------------
if [[ ! -s "$OUT" ]]; then
    echo "ERROR: SMC++ output was not created:"
    echo "$OUT"
    exit 1
fi

gzip -t "$OUT"


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