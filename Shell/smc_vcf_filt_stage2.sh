#!/bin/bash
#SBATCH --job-name=smc_vcf_filt_2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ================================================================
# SMC++ VCF filtering - Stage 2
#
# starting data:
#   stage-1 filtered autosomal VCFs
#
# stage-1 already applied:
#   shared callable genome
#   SNPs only
#   biallelic sites only
#   QUAL >= 20
#   genotype DP 5-25
#   genotype GQ >= 20
#   failed genotypes set to ./.
#
# stage-2 filter:
#   retain sites with no more than 1 missing sample
#
# therefore:
#   0 missing samples -> KEEP
#   1 missing sample  -> KEEP
#   >=2 missing       -> REMOVE
#
# this requires at least 7 of 8 samples called at each SNP.
# ================================================================


# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
INDIR="${WORKDIR}/04_vcf/filtered_stage1"
STAGE1STATDIR="${WORKDIR}/04_vcf/filtered_stage1_stats"
OUTDIR="${WORKDIR}/04_vcf/filtered_stage2"
STATDIR="${WORKDIR}/04_vcf/filtered_stage2_stats"

mkdir -p "$OUTDIR"
mkdir -p "$STATDIR"


# ------------------------------------------------------------
# get chromosome
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")


# ------------------------------------------------------------
# files
# ------------------------------------------------------------
IN="${INDIR}/${CHR}.filtered_stage1.vcf.gz"
OUT="${OUTDIR}/${CHR}.filtered_stage2.vcf.gz"
STAGE1STATS="${STAGE1STATDIR}/${CHR}.filtered_stage1.stats.txt"
STATS="${STATDIR}/${CHR}.filtered_stage2.stats.txt"

# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ VCF filtering - Stage 2"
echo "============================================================"
echo "Chromosome:       $CHR"
echo "Input:            $IN"
echo "Output:           $OUT"
echo "Missingness rule: N_MISSING <= 1"
echo "Minimum called:   7 of 8 samples"
echo "Array task:       ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "Host:             $(hostname)"
echo "Start:            $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# stage-2 filtering
#
# N_MISSING is calculated by bcftools from the GT fields.
#
# keep sites where:
#
#     N_MISSING <= 1
#
# this retains sites with:
#     8/8 called
#     7/8 called
#
# and removes sites with:
#     <=6/8 called
# ------------------------------------------------------------
bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -i 'N_MISSING<=1' \
    -Oz \
    -o "$OUT" \
    "$IN"


# ------------------------------------------------------------
# index output
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$OUT"


# ------------------------------------------------------------
# generate statistics
# ------------------------------------------------------------
bcftools stats \
    -s - \
    "$OUT" \
    > "$STATS"


# ------------------------------------------------------------
# validate number of samples
# ------------------------------------------------------------
NSAMPLES=$(bcftools query -l "$OUT" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: Stage-2 VCF contains $NSAMPLES samples"
    exit 1
fi


# ------------------------------------------------------------
# get Stage-1 record count
# ------------------------------------------------------------
if [[ -s "$STAGE1STATS" ]]; then

    STAGE1_N=$(
        awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "$STAGE1STATS"
    )

else

    STAGE1_N="NA"

fi


# ------------------------------------------------------------
# get Stage-2 record count
# ------------------------------------------------------------
STAGE2_N=$(
    awk -F'\t' \
    '$1=="SN" && $3=="number of records:" {print $4}' \
    "$STATS"
)


# ------------------------------------------------------------
# confirm no sites with >=2 missing samples remain
# ------------------------------------------------------------
BAD_MISSING=$(
    bcftools view \
        -H \
        -i 'N_MISSING>1' \
        "$OUT" |
    wc -l
)

if [[ "$BAD_MISSING" -ne 0 ]]; then
    echo "ERROR: $BAD_MISSING sites with >1 missing sample remain"
    exit 1
fi


# ------------------------------------------------------------
# confirm output still contains only SNPs
# ------------------------------------------------------------
NON_SNP=$(
    bcftools view \
        -H \
        -V snps \
        "$OUT" |
    wc -l
)

if [[ "$NON_SNP" -ne 0 ]]; then
    echo "ERROR: $NON_SNP non-SNP records remain"
    exit 1
fi


# ------------------------------------------------------------
# confirm output remains biallelic
# ------------------------------------------------------------
MULTI=$(
    bcftools view \
        -H \
        -m3 \
        "$OUT" |
    wc -l
)

if [[ "$MULTI" -ne 0 ]]; then
    echo "ERROR: $MULTI multiallelic records remain"
    exit 1
fi


# ------------------------------------------------------------
# calculate number and percent removed
# ------------------------------------------------------------
if [[ "$STAGE1_N" != "NA" && "$STAGE1_N" -gt 0 ]]; then

    REMOVED=$((STAGE1_N - STAGE2_N))

    PCT_RETAINED=$(
        awk \
            -v kept="$STAGE2_N" \
            -v total="$STAGE1_N" \
            'BEGIN {printf "%.2f",100*kept/total}'
    )

    PCT_REMOVED=$(
        awk \
            -v removed="$REMOVED" \
            -v total="$STAGE1_N" \
            'BEGIN {printf "%.2f",100*removed/total}'
    )

else

    REMOVED="NA"
    PCT_RETAINED="NA"
    PCT_REMOVED="NA"

fi


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Completed chromosome:     $CHR"
echo "Stage-1 records:          $STAGE1_N"
echo "Stage-2 records:          $STAGE2_N"
echo "Removed:                  $REMOVED"
echo "Percent retained:         ${PCT_RETAINED}%"
echo "Percent removed:          ${PCT_REMOVED}%"
echo "Samples:                  $NSAMPLES"
echo "Sites with >1 missing:    $BAD_MISSING"
echo "Non-SNP records:          $NON_SNP"
echo "Multiallelic records:     $MULTI"
echo "Output:                   $OUT"
echo "Finish:                   $(date)"
echo "============================================================"