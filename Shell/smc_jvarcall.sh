#!/bin/bash
#SBATCH --job-name=smc_jvarcall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ============================================================
# Joint variant calling for SMC++
#
# Eight mainland G. ussuriensis individuals are jointly called
# One autosome is processed per SLURM array task
#
# This is RAW variant calling
# Final SNP/DP/QUAL/missingness filtering is done afterward
# ============================================================

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
REF="${WORKDIR}/01_reference/G_ussuriensis.softmasked.fa"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
BAMLIST="${WORKDIR}/02_bams/mainland_bams.txt"
OUTDIR="${WORKDIR}/04_vcf/raw"
STATDIR="${WORKDIR}/04_vcf/raw_stats"

mkdir -p "$OUTDIR"
mkdir -p "$STATDIR"


# ------------------------------------------------------------
# check input files
# ------------------------------------------------------------
if [[ ! -s "$REF" ]]; then
    echo "ERROR: Reference not found:"
    echo "$REF"
    exit 1
fi

if [[ ! -s "${REF}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found:"
    echo "${REF}.fai"
    exit 1
fi

if [[ ! -s "$CHRFILE" ]]; then
    echo "ERROR: Autosome list not found:"
    echo "$CHRFILE"
    exit 1
fi

if [[ ! -s "$BAMLIST" ]]; then
    echo "ERROR: BAM list not found:"
    echo "$BAMLIST"
    exit 1
fi

# ------------------------------------------------------------
# verify number of BAMs
# ------------------------------------------------------------
NBAMS=$(grep -cve '^[[:space:]]*$' "$BAMLIST")

if [[ "$NBAMS" -ne 8 ]]; then
    echo "ERROR: Expected 8 BAMs but found $NBAMS"
    echo "BAM list: $BAMLIST"
    exit 1
fi

echo "Found $NBAMS BAMs in BAM list."


# ------------------------------------------------------------
# get chromosome for this array task
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")

if [[ -z "$CHR" ]]; then
    echo "ERROR: No chromosome found for array task:"
    echo "${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


# ------------------------------------------------------------
# output files
# ------------------------------------------------------------
VCF="${OUTDIR}/${CHR}.raw.vcf.gz"
STATS="${STATDIR}/${CHR}.raw.stats.txt"


# ------------------------------------------------------------
# print job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Joint variant calling"
echo "============================================================"
echo "Chromosome:       $CHR"
echo "Number of BAMs:   $NBAMS"
echo "Reference:        $REF"
echo "BAM list:         $BAMLIST"
echo "Output:           $VCF"
echo "Array task:       ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "Start time:       $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# check BAMs before calling
# ------------------------------------------------------------
while read -r bam; do

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: BAM not found:"
        echo "$bam"
        exit 1
    fi

    samtools quickcheck "$bam"

done < "$BAMLIST"


# ------------------------------------------------------------
# joint variant calling
#
# mpileup:
#
#   -f       reference genome
#   -r       current chromosome
#   -b       BAM list containing all 8 individuals
#   -q 30    require mapping quality >=30
#   -Q 20    require base quality >=20
#   -d 1000  allow up to 1000 reads/sample/site
#
#   --ff:
#     exclude unmapped, secondary, QC-failed,
#     and duplicate-marked reads
#
#   -a:
#     retain per-individual depth (DP)
#     and allele depth (AD)
#
# call:
#
#   -m       multiallelic calling model
#   -v       output variant sites only
#   -f GQ    include genotype quality
#
# no final SNP or depth filtering is performed here.
# ------------------------------------------------------------

bcftools mpileup \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f "$REF" \
    -r "$CHR" \
    -b "$BAMLIST" \
    -q 30 \
    -Q 20 \
    -d 1000 \
    --ff UNMAP,SECONDARY,QCFAIL,DUP \
    -a FORMAT/DP,FORMAT/AD \
    -Ou \
| \
bcftools call \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -m \
    -v \
    -f GQ \
    -Oz \
    -o "$VCF"


# ------------------------------------------------------------
# index VCF
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$VCF"


# ------------------------------------------------------------
# generate raw VCF statistics
# ------------------------------------------------------------
bcftools stats "$VCF" > "$STATS"


# ------------------------------------------------------------
# quick output validation
# ------------------------------------------------------------
echo
echo "Samples in VCF:"

bcftools query -l "$VCF"

NSAMPLES=$(bcftools query -l "$VCF" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: VCF contains $NSAMPLES samples instead of 8."
    exit 1
fi


# ------------------------------------------------------------
# count raw variants
# ------------------------------------------------------------
NVAR=$(bcftools view -H "$VCF" | wc -l)

echo
echo "Raw variant records: $NVAR"


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Completed chromosome: $CHR"
echo "Raw variants:         $NVAR"
echo "Output:               $VCF"
echo "Finish time:          $(date)"
echo "============================================================"