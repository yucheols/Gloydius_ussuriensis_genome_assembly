#!/bin/bash
#SBATCH --job-name=smc_vcf_filt_1
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
# SMC++ VCF filtering - Stage 1
#
# 1) starting data:
#   raw joint-called autosomal VCFs
#
# 2) site filters:
#   shared callable genome only
#   SNPs only
#   biallelic sites only
#   QUAL >= 20
#
# 3) genotype filters:
#   DP 5-25
#   GQ >= 20
#
# genotypes failing DP/GQ are changed to ./.
#
# this script DOES NOT remove sites based on missingness just yet.
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
INDIR="${WORKDIR}/04_vcf/raw"
OUTDIR="${WORKDIR}/04_vcf/filtered_stage1"
STATDIR="${WORKDIR}/04_vcf/filtered_stage1_stats"
CALLABLE="${WORKDIR}/05_masks/shared_callable_depth.autosomes.merged.bed"
RAWSTATDIR="${WORKDIR}/04_vcf/raw_stats"

mkdir -p "$OUTDIR"
mkdir -p "$STATDIR"


# ------------------------------------------------------------
# get chromosome
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")


# ------------------------------------------------------------
# files
# ------------------------------------------------------------
RAW="${INDIR}/${CHR}.raw.vcf.gz"
OUT="${OUTDIR}/${CHR}.filtered_stage1.vcf.gz"
STATS="${STATDIR}/${CHR}.filtered_stage1.stats.txt"
RAWSTATS="${RAWSTATDIR}/${CHR}.raw.stats.txt"


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------

echo
echo "============================================================"
echo "SMC++ VCF filtering - Stage 1"
echo "============================================================"
echo "Chromosome:     $CHR"
echo "Raw VCF:        $RAW"
echo "Callable BED:   $CALLABLE"
echo "Output:         $OUT"
echo "Array task:     ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:           ${SLURM_CPUS_PER_TASK}"
echo "Host:           $(hostname)"
echo "Start:          $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# FILTERING
#
# 1) Step 1:
#   Remove INFO/MQ because its raw header definition causes
#   the bcftools MQ sanity warning and we do not use this tag.
#
# 2) Step 2:
#   Retain positions within shared callable regions.
#
# 3)  Step 3:
#   Retain SNPs only.
#
# 4) Step 4:
#   Retain exactly two alleles (REF + one ALT).
#
# 5) Step 5:
#   Require site QUAL >=20.
#
# 6) Step 6:
#   Change an individual genotype to missing (./.) if:
#
#       DP <5
#       DP >25
#       DP missing
#       GQ <20
#       GQ missing
#
# IMPORTANT:
# Single "|" operators are intentional. For FORMAT fields,
# this evaluates the conditions at the individual-sample level.
# ------------------------------------------------------------
bcftools annotate \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -x INFO/MQ \
    -Ou \
    "$RAW" \
| \
bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -T "$CALLABLE" \
    --targets-overlap 0 \
    -v snps \
    -m2 \
    -M2 \
    -i 'QUAL>=20' \
    -Ou \
| \
bcftools +setGT \
    -Ou \
    -- \
    -t q \
    -n . \
    -i 'FMT/DP="." | FMT/DP<5 | FMT/DP>25 | FMT/GQ="." | FMT/GQ<20' \
| \
bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "$OUT"


# ------------------------------------------------------------
# index filtered VCF
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$OUT"


# ------------------------------------------------------------
# statistics
#
# -s - tells bcftools stats to generate per-sample statistics
# for all samples.
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
    echo "ERROR: Filtered VCF contains $NSAMPLES samples"
    exit 1
fi


# ------------------------------------------------------------
# obtain raw and filtered record counts from stats
# ------------------------------------------------------------
if [[ -s "$RAWSTATS" ]]; then

    RAW_N=$(
        awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "$RAWSTATS"
    )

else

    RAW_N="NA"

fi


FILTERED_N=$(
    awk -F'\t' \
    '$1=="SN" && $3=="number of records:" {print $4}' \
    "$STATS"
)


# ------------------------------------------------------------
# confirm output contains only biallelic SNPs
# ------------------------------------------------------------
NON_SNP=$(
    bcftools view \
        -H \
        -V snps \
        "$OUT" |
    wc -l
)

MULTI=$(
    bcftools view \
        -H \
        -m3 \
        "$OUT" |
    wc -l
)

if [[ "$NON_SNP" -ne 0 ]]; then
    echo "ERROR: $NON_SNP non-SNP records remain"
    exit 1
fi

if [[ "$MULTI" -ne 0 ]]; then
    echo "ERROR: $MULTI multiallelic records remain"
    exit 1
fi


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Completed chromosome:   $CHR"
echo "Raw records:            $RAW_N"
echo "Stage-1 records:        $FILTERED_N"
echo "Samples:                $NSAMPLES"
echo "Non-SNP records:        $NON_SNP"
echo "Multiallelic records:   $MULTI"
echo "Output:                 $OUT"
echo "Finish:                 $(date)"
echo "============================================================"