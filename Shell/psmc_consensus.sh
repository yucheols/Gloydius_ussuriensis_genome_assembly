#!/bin/bash
#SBATCH --job-name=psmc_consensus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================================
# generate diploid consensus FASTQ for PSMC
# individual: AMNH_21010
#
# input: Illumina WGS mapped to the AMNH_21010 chromosome-level assembly
# regions: 17 autosomes only
# filtering:
#   minimum mapping quality = 30
#   minimum base quality    = 20
#   minimum depth           = 4x
#   maximum depth           = 24x
#
# output: diploid consensus FASTQ for downstream fq2psmcfa
# ============================================================================


# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate psmc

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"

REF="${PROJECT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
BAM="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam/AMNH_21010.markdup.bam"
AUTOSOMES="${PSMCROOT}/00_setup/autosomes.bed"
OUTDIR="${PSMCROOT}/02_consensus"

mkdir -p "${OUTDIR}"


# ------------------------------------------------------------
# parameters
# ------------------------------------------------------------
SAMPLE="AMNH_21010"

MIN_MQ=30
MIN_BQ=20

MIN_DP=4
MAX_DP=24

CONSENSUS_MQ=30

# small region used to verify that current bcftools output
# is compatible with vcfutils.pl vcf2fq
TEST_REGION="G_ussuri_chr18:1-1000000"


# ------------------------------------------------------------
# output names
# ------------------------------------------------------------
PREFIX="${OUTDIR}/${SAMPLE}.autosomes.dp${MIN_DP}-${MAX_DP}"
TEST_VCF="${OUTDIR}/${SAMPLE}.preflight.vcf"
TEST_FQ="${OUTDIR}/${SAMPLE}.preflight.fq"
DIPLOID_FQ="${PREFIX}.diploid.fq.gz"


# ------------------------------------------------------------
# check inputs
# ------------------------------------------------------------

echo "============================================================"
echo "PSMC diploid consensus generation"
echo "============================================================"
echo "Sample:        ${SAMPLE}"
echo "Reference:     ${REF}"
echo "BAM:           ${BAM}"
echo "Autosomes:     ${AUTOSOMES}"
echo "MAPQ cutoff:   ${MIN_MQ}"
echo "BaseQ cutoff:  ${MIN_BQ}"
echo "Depth range:   ${MIN_DP}-${MAX_DP}x"
echo "============================================================"
echo

for FILE in "${REF}" "${REF}.fai" "${BAM}" "${BAM}.bai" "${AUTOSOMES}"; do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: missing or empty input: ${FILE}" >&2
        exit 1
    fi
done

for PROGRAM in bcftools samtools vcfutils.pl fq2psmcfa; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} not found in PATH" >&2
        exit 1
    fi
done

samtools quickcheck -v "${BAM}"


# ------------------------------------------------------------
# record software versions
# ------------------------------------------------------------
{
    echo "bcftools"
    bcftools --version | head -2
    echo
    echo "samtools"
    samtools --version | head -2
    echo
    echo "vcfutils.pl"
    command -v vcfutils.pl
    echo
    echo "fq2psmcfa"
    command -v fq2psmcfa
} > "${OUTDIR}/software_versions.txt"


# ------------------------------------------------------------
# preflight test
#
# vcfutils.pl vcf2fq was written around the original bcftools
# consensus-caller annotations. test a small autosomal region
# before processing the complete 1.38-Gb autosomal genome.
# ------------------------------------------------------------
echo
echo "Running preflight test on ${TEST_REGION} ..."
echo

bcftools mpileup \
    -f "${REF}" \
    -r "${TEST_REGION}" \
    -q "${MIN_MQ}" \
    -Q "${MIN_BQ}" \
    --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
    -Ou \
    "${BAM}" \
| bcftools call \
    -c \
    -Ov \
    -o "${TEST_VCF}"


# ------------------------------------------------------------
# verify single sample
# ------------------------------------------------------------
N_SAMPLES=$(bcftools query -l "${TEST_VCF}" | wc -l)

if [[ "${N_SAMPLES}" -ne 1 ]]; then
    echo "ERROR: expected one sample in VCF but found ${N_SAMPLES}" >&2
    bcftools query -l "${TEST_VCF}" >&2
    exit 1
fi

echo "VCF sample:"
bcftools query -l "${TEST_VCF}"


# ------------------------------------------------------------
# verify fields required by vcf2fq
# ------------------------------------------------------------
for TAG in DP MQ FQ; do
    if ! grep -q "##INFO=<ID=${TAG}," "${TEST_VCF}"; then
        echo "ERROR: INFO/${TAG} is absent from bcftools output." >&2
        echo "Current bcftools output is not directly compatible with the" >&2
        echo "legacy vcfutils.pl vcf2fq workflow." >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# verify vcf2fq conversion
# ------------------------------------------------------------
vcfutils.pl vcf2fq \
    -d "${MIN_DP}" \
    -D "${MAX_DP}" \
    -Q "${CONSENSUS_MQ}" \
    "${TEST_VCF}" \
    > "${TEST_FQ}"

if [[ ! -s "${TEST_FQ}" ]]; then
    echo "ERROR: preflight diploid FASTQ is empty" >&2
    exit 1
fi

echo
echo "Preflight test passed."
echo


# ------------------------------------------------------------
# generate full autosomal diploid consensus
#
# IMPORTANT:
# Do not use bcftools call -v here.
# vcfutils.pl vcf2fq requires an all-sites VCF so that both
# homozygous and heterozygous regions are represented.
# ------------------------------------------------------------
echo "Generating full autosomal diploid consensus ..."
echo

bcftools mpileup \
    -f "${REF}" \
    -R "${AUTOSOMES}" \
    -q "${MIN_MQ}" \
    -Q "${MIN_BQ}" \
    --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
    -Ou \
    "${BAM}" \
| bcftools call \
    -c \
    -Ov \
| vcfutils.pl vcf2fq \
    -d "${MIN_DP}" \
    -D "${MAX_DP}" \
    -Q "${CONSENSUS_MQ}" \
| gzip -c \
    > "${DIPLOID_FQ}"


# ------------------------------------------------------------
# validate output
# ------------------------------------------------------------
gzip -t "${DIPLOID_FQ}"

if [[ ! -s "${DIPLOID_FQ}" ]]; then
    echo "ERROR: diploid consensus FASTQ is empty" >&2
    exit 1
fi


# ------------------------------------------------------------
# report
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Consensus generation completed"
echo "============================================================"
echo "Output:"
echo "${DIPLOID_FQ}"
echo
ls -lh "${DIPLOID_FQ}"
echo "============================================================"