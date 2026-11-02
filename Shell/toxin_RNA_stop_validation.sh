#!/bin/bash
#SBATCH --job-name=toxin_RNA_stop_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# directories
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

# input files
REF="${ROOT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
RNA_BAM="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.sorted.bam"
BED="${WORKDIR}/toxin_internal_stops.bed"

# output
OUT="${WORKDIR}/toxin_internal_stops_SRR35908235.mpileup"
TMP="${OUT}.tmp.${SLURM_JOB_ID:-$$}"

cd "$WORKDIR"

trap 'rm -f "$TMP"' EXIT

# ------------------------------------------------------------
# RNA pileup using indexed genomic regions
#
# BED: 0-based, half-open
# mpileup -r: 1-based, inclusive
#
# -a: include zero-depth positions in requested regions
# -q 20: minimum mapping quality
# -Q 20: minimum base quality
# -d 100000: maximum depth
# ------------------------------------------------------------

: > "$TMP"

count=0

while IFS=$'\t' read -r chrom start end gene score strand; do

    region="${chrom}:$((start + 1))-${end}"

    echo "Processing ${gene}: ${region}"

    samtools mpileup \
        -f "$REF" \
        -r "$region" \
        -q 20 \
        -Q 20 \
        -d 100000 \
        -a \
        "$RNA_BAM" \
        >> "$TMP"

    count=$((count + 1))

done < "$BED"


# validate output
lines=$(wc -l < "$TMP")

echo
echo "===== RNA STOP VALIDATION ====="
echo "Stop codons processed: $count"
echo "Pileup positions:       $lines"

if [[ "$count" -ne 14 || "$lines" -ne 42 ]]; then
    echo "ERROR: Expected 14 codons and 42 positions." >&2
    exit 1
fi

mv -f "$TMP" "$OUT"

echo "Output: $OUT"
echo "Status: PASS"