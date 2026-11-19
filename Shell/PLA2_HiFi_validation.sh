#!/bin/bash
#SBATCH --job-name=PLA2_HiFi_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# set paths
ROOT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'

REF="${ROOT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
HIFI="${ROOT}/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz"

OUTDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/all120_QC/PLA2_HiFi_validation"

mkdir -p "$OUTDIR"

REGION="${OUTDIR}/PLA2_chr16_region.fasta"
PAF="${OUTDIR}/PLA2_HiFi_alignments.paf"
SUMMARY="${OUTDIR}/PLA2_HiFi_support_summary.tsv"
BRIDGES="${OUTDIR}/PLA2_HiFi_bridging_reads.tsv"

# validate required tools
for tool in samtools minimap2 python; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Missing $tool" >&2
        exit 1
    }
done

[[ -s "$REF" && -s "$HIFI" ]] || {
    echo 'ERROR: Missing reference or HiFi FASTQ.' >&2
    exit 1
}

# ------------------------------------------------------------
# STEP 1: extract the genomic interval
# ------------------------------------------------------------

echo '===== STEP 1: EXTRACT REFERENCE REGION ====='

# grab the chr16 region spanning PLA2-5, PLA2-6, and flanking sequence
samtools faidx "$REF" 'G_ussuri_chr16:289000-310000' > "$REGION"

echo "Reference: $REGION"

# ------------------------------------------------------------
# STEP 2: align original HiFi reads
# ------------------------------------------------------------

echo
echo '===== STEP 2: TARGETED HIFI ALIGNMENT ====='

minimap2 \
    -x map-hifi \
    -c \
    --secondary=no \
    -t "${SLURM_CPUS_PER_TASK:-8}" \
    "$REGION" \
    "$HIFI" \
    > "$PAF" \
    2> "${OUTDIR}/minimap2.log"

echo "PAF file: $PAF"

# ------------------------------------------------------------
# STEP 3: evaluate spanning HiFi reads
# ------------------------------------------------------------

echo
echo '===== STEP 3: TEST TANDEM DUPLICATION SUPPORT ====='

python - "$PAF" "$SUMMARY" "$BRIDGES" <<'PY'
import sys
import csv
import re
from collections import defaultdict

paf_file, summary_file, bridges_file = sys.argv[1:]

# Original genomic coordinates (1-based inclusive)
region_start = 289000

pla2_5 = (293721, 295345)
pla2_6 = (302420, 304043)

# Require alignment 250 bp beyond both genes.
flank = 250

bridge_left = pla2_5[0] - flank
bridge_right = pla2_6[1] + flank

# Minimum alignment identity.
min_identity = 0.97

# Maximum single deletion allowed within a bridging alignment.
max_deletion = 100

aligned = set()
span_5 = set()
span_6 = set()
bridge_any = set()
bridge_high_quality = {}

def parse_cigar(text):
    operations = re.findall(r'(\d+)([MIDNSHP=X])', text)

    if ''.join(n + op for n, op in operations) != text:
        return None

    return [(int(n), op) for n, op in operations]

with open(paf_file) as f:

    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip('\n').split('\t')

        if len(fields) < 12:
            continue

        qname = fields[0]

        try:
            qlength = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            alignment_length = int(fields[10])
            mapq = int(fields[11])
        except ValueError:
            continue

        if alignment_length == 0:
            continue

        tags = {}

        for tag in fields[12:]:
            parts = tag.split(':', 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]

        # Exclude non-primary PAF alignments.
        if tags.get('tp', 'P') != 'P':
            continue

        aligned.add(qname)

        # Convert target-local PAF coordinates to genome coordinates.
        # PAF target start is 0-based; target end is exclusive.
        genome_start = region_start + tstart
        genome_end = region_start + tend - 1

        identity = matches / alignment_length

        if (
            genome_start <= pla2_5[0]
            and genome_end >= pla2_5[1]
        ):
            span_5.add(qname)

        if (
            genome_start <= pla2_6[0]
            and genome_end >= pla2_6[1]
        ):
            span_6.add(qname)

        if not (
            genome_start <= bridge_left
            and genome_end >= bridge_right
        ):
            continue

        bridge_any.add(qname)

        cigar_text = tags.get('cg')

        if not cigar_text:
            continue

        operations = parse_cigar(cigar_text)

        if operations is None:
            continue

        # Reject alignments with large reference gaps.
        large_gap = any(
            op == 'N' or (op == 'D' and n > max_deletion)
            for n, op in operations
        )

        if large_gap or identity < min_identity:
            continue

        record = {
            'read': qname,
            'read_length': qlength,
            'query_start': qstart,
            'query_end': qend,
            'genome_start': genome_start,
            'genome_end': genome_end,
            'identity': identity,
            'mapq': mapq,
            'cigar': cigar_text
        }

        previous = bridge_high_quality.get(qname)

        if previous is None or identity > previous['identity']:
            bridge_high_quality[qname] = record

# ------------------------------------------------------------
# Save summary
# ------------------------------------------------------------

statistics = [
    ('Distinct_HiFi_reads_aligned_to_target', len(aligned)),
    ('Reads_spanning_PLA2_5', len(span_5)),
    ('Reads_spanning_PLA2_6', len(span_6)),
    ('Reads_spanning_both_plus_250bp', len(bridge_any)),
    ('High_identity_continuous_bridging_reads', len(bridge_high_quality))
]

with open(summary_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerow(['Metric', 'Count'])
    writer.writerows(statistics)

with open(bridges_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Read_ID',
        'Read_length',
        'Query_start_0based',
        'Query_end_exclusive',
        'Alignment_start',
        'Alignment_end',
        'Alignment_identity_percent',
        'MAPQ_local_reference',
        'CIGAR'
    ])

    for read in sorted(bridge_high_quality):
        r = bridge_high_quality[read]

        writer.writerow([
            r['read'],
            r['read_length'],
            r['query_start'],
            r['query_end'],
            r['genome_start'],
            r['genome_end'],
            round(100*r['identity'], 3),
            r['mapq'],
            r['cigar']
        ])

print('\n===== PLA2 TANDEM DUPLICATION HIFI SUPPORT =====')

for name, value in statistics:
    print('{}\t{}'.format(name, value))

print('\n===== BRIDGING READS (FIRST 20) =====')
print('Read_ID\tAlignment_start\tAlignment_end\tIdentity_%')

for read in sorted(bridge_high_quality)[:20]:
    r = bridge_high_quality[read]

    print('{}\t{}\t{}\t{:.3f}'.format(
        read,
        r['genome_start'],
        r['genome_end'],
        100*r['identity']
    ))

print('\nSummary:', summary_file)
print('Bridging reads:', bridges_file)
PY

echo
echo '===== COMPLETED ====='
echo 'Status: PASS'