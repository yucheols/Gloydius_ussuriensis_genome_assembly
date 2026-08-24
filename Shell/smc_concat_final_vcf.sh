#!/bin/bash
#SBATCH --job-name=smc_vcf_concat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.err

# ================================================================
# Concatenate Stage-2 autosomal VCFs and generate final QC stats
# ================================================================

source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"

CHRFILE="${WORKDIR}/01_reference/autosomes.txt"

INDIR="${WORKDIR}/04_vcf/filtered_stage2"

OUTDIR="${WORKDIR}/04_vcf/final"
QCDIR="${WORKDIR}/04_vcf/final_qc"

VCFLIST="${OUTDIR}/stage2_vcf_list.txt"

FINALVCF="${OUTDIR}/G_ussuriensis_mainland.autosomes.filtered.vcf.gz"

STATS="${QCDIR}/G_ussuriensis_mainland.autosomes.filtered.stats.txt"

mkdir -p "$OUTDIR"
mkdir -p "$QCDIR"


# ------------------------------------------------------------
# create VCF list in chromosome order
# ------------------------------------------------------------
> "$VCFLIST"

while read -r chr; do

    VCF="${INDIR}/${chr}.filtered_stage2.vcf.gz"

    if [[ ! -s "$VCF" ]]; then
        echo "ERROR: Missing VCF:"
        echo "$VCF"
        exit 1
    fi

    if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
        echo "ERROR: Missing VCF index:"
        echo "$VCF"
        exit 1
    fi

    echo "$VCF" >> "$VCFLIST"

done < "$CHRFILE"


# ------------------------------------------------------------
# verify 17 chromosome VCFs
# ------------------------------------------------------------
NVCF=$(wc -l < "$VCFLIST")

if [[ "$NVCF" -ne 17 ]]; then
    echo "ERROR: Expected 17 VCFs but found $NVCF"
    exit 1
fi


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo "============================================================"
echo "SMC++ final autosomal VCF"
echo "============================================================"
echo "Number of chromosomes: $NVCF"
echo "Output:                $FINALVCF"
echo "Start:                 $(date)"
echo "============================================================"


# ------------------------------------------------------------
# concatenate chromosomes
# ------------------------------------------------------------
bcftools concat \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f "$VCFLIST" \
    -Oz \
    -o "$FINALVCF"


# ------------------------------------------------------------
# index
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$FINALVCF"


# ------------------------------------------------------------
# final statistics including per-sample statistics
# ------------------------------------------------------------
bcftools stats \
    -s - \
    "$FINALVCF" \
    > "$STATS"


# ------------------------------------------------------------
# validation
# ------------------------------------------------------------
NSAMPLES=$(bcftools query -l "$FINALVCF" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: Expected 8 samples but found $NSAMPLES"
    exit 1
fi


NRECORDS=$(
    awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "$STATS"
)


NON_SNP=$(
    bcftools view \
        -H \
        -V snps \
        "$FINALVCF" |
    wc -l
)


MULTI=$(
    bcftools view \
        -H \
        -m3 \
        "$FINALVCF" |
    wc -l
)


BAD_MISSING=$(
    bcftools view \
        -H \
        -i 'N_MISSING>1' \
        "$FINALVCF" |
    wc -l
)


# ------------------------------------------------------------
# final report
# ------------------------------------------------------------
echo
echo "============================================================"
echo "FINAL VCF VALIDATION"
echo "============================================================"
echo "Records:                  $NRECORDS"
echo "Samples:                  $NSAMPLES"
echo "Non-SNP records:          $NON_SNP"
echo "Multiallelic records:     $MULTI"
echo "Sites with >1 missing:    $BAD_MISSING"
echo "Final VCF:                $FINALVCF"
echo "Stats:                    $STATS"
echo "Finish:                   $(date)"
echo "============================================================"