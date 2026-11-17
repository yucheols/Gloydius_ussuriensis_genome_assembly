#!/bin/bash
#SBATCH --job-name=toxin_gene_model_review
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# ------------------------------------------------------------
# directories and inputs
# ------------------------------------------------------------
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
OUTDIR="${WORKDIR}/gene_model_review"

REF="${ROOT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
GTF="${WORKDIR}/toxin_annotation.gtf"
INVENTORY="${WORKDIR}/toxin_transcript_evidence_paths.tsv"

TRANSCRIPTS="${OUTDIR}/regional_trinity_transcripts.fasta"
METADATA="${OUTDIR}/regional_trinity_metadata.tsv"

BAM="${OUTDIR}/regional_trinity_genome_alignment.bam"

SUMMARY="${OUTDIR}/gene_model_review_summary.tsv"
DETAILS="${OUTDIR}/gene_model_review_transcripts.tsv"
BED="${OUTDIR}/gene_model_review_transcripts.bed"
MODEL_GTF="${OUTDIR}/original_toxin_models.gtf"

mkdir -p "$OUTDIR"

for tool in python minimap2 samtools; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Missing tool: $tool" >&2
        exit 1
    }
done

for file in "$REF" "$GTF" "$INVENTORY"; do
    [[ -s "$file" ]] || {
        echo "ERROR: Missing input: $file" >&2
        exit 1
    }
done

# ------------------------------------------------------------
# STEP 1: collect regional Trinity transcripts
# ------------------------------------------------------------

echo "===== STEP 1: COLLECT REGIONAL TRANSCRIPTS ====="

python - "$GTF" "$INVENTORY" "$TRANSCRIPTS" "$METADATA" "$MODEL_GTF" <<'PY'
import sys
import csv
import re
import hashlib
from collections import defaultdict

gtf_file, inventory_file, output_fasta, output_metadata, output_gtf = sys.argv[1:]

genes = [
    'SVSP-16',
    'Ficolin-4',
    'PLA2-4',
    'CTL-10',
    'HYAL-4',
    'CTL-16'
]

# Save original gene models
found = set()

with open(gtf_file) as source, open(output_gtf, 'w') as output:
    for line in source:
        match = re.search(r'gene_id "([^"]+)"', line)
        if match and match.group(1) in genes:
            output.write(line)
            found.add(match.group(1))

missing = set(genes) - found

if missing:
    raise RuntimeError(
        'Missing toxin models: ' + ','.join(sorted(missing))
    )

# Read the inventory produced during Step 4
paths = defaultdict(set)

with open(inventory_file) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        if row['Gene'] in genes:
            if row['Evidence_type'] == 'Trinity_region':
                paths[row['Gene']].add(row['Path'])

def read_fasta(path):
    name = None
    parts = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts)

                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    if name is not None:
        yield name, ''.join(parts)

# Deduplicate identical sequences across regional FASTA files
sequences = {}
sources = defaultdict(set)

for gene in genes:

    for path in sorted(paths[gene]):

        import os

        if not os.path.isfile(path):
            raise RuntimeError(
                'Missing Trinity FASTA: ' + path
            )

        for original_id, sequence in read_fasta(path):

            sequence = sequence.upper().replace('U', 'T')

            if not sequence:
                continue

            digest = hashlib.sha256(
                sequence.encode('ascii')
            ).hexdigest()

            sequences[digest] = sequence
            sources[digest].add((gene, original_id, path))

with open(output_fasta, 'w') as fasta, \
     open(output_metadata, 'w', newline='') as meta:

    writer = csv.writer(meta, delimiter='\t')

    writer.writerow([
        'Transcript_ID',
        'Length_bp',
        'Source_genes',
        'Original_IDs',
        'Source_paths'
    ])

    for index, digest in enumerate(sorted(sequences), 1):

        transcript_id = 'TX{:06d}'.format(index)
        sequence = sequences[digest]

        annotations = sorted(sources[digest])

        source_genes = sorted(set(x[0] for x in annotations))
        original_ids = sorted(set(x[1] for x in annotations))
        source_paths = sorted(set(x[2] for x in annotations))

        fasta.write(
            '>{}\n{}\n'.format(transcript_id, sequence)
        )

        writer.writerow([
            transcript_id,
            len(sequence),
            ','.join(source_genes),
            ','.join(original_ids),
            ';'.join(source_paths)
        ])

print('Trinity transcripts collected:', len(sequences))

for gene in genes:
    count = sum(
        gene in set(x[0] for x in annotations)
        for annotations in sources.values()
    )

    print('{}\t{}'.format(gene, count))

if not sequences:
    raise RuntimeError('No regional Trinity transcripts found')
PY


# ------------------------------------------------------------
# STEP 2: align assembled transcripts to the ref genome
#
# this helps to identify paralogous placements
# do not force transcript strand with -uf
# ------------------------------------------------------------

echo
echo "===== STEP 2: SPLICE-AWARE TRANSCRIPT ALIGNMENT ====="

minimap2 \
    -ax splice:hq \
    -N 10 \
    -t "${SLURM_CPUS_PER_TASK:-8}" \
    "$REF" \
    "$TRANSCRIPTS" \
    2> "${OUTDIR}/minimap2.log" | \
samtools sort \
    -@ 4 \
    -m 1G \
    -o "$BAM" -

samtools index "$BAM"

echo "Alignment complete."


# ------------------------------------------------------------
# STEP 3: compare aligned transcripts with toxin gene models
# ------------------------------------------------------------

echo
echo "===== STEP 3: GENE-MODEL COMPARISON ====="

python - "$GTF" "$METADATA" "$BAM" "$SUMMARY" "$DETAILS" "$BED" <<'PY'
import sys
import csv
import re
import subprocess
from collections import defaultdict

gtf_file, metadata_file, bam_file, summary_file, details_file, bed_file = sys.argv[1:]

genes = [
    'SVSP-16',
    'Ficolin-4',
    'PLA2-4',
    'CTL-10',
    'HYAL-4',
    'CTL-16'
]

models = defaultdict(list)
metadata = {}

with open(metadata_file) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        metadata[row['Transcript_ID']] = row

with open(gtf_file) as f:
    for line in f:
        if line.startswith('#'):
            continue

        fields = line.rstrip('\n').split('\t')

        if len(fields) < 9 or fields[2] != 'CDS':
            continue

        match = re.search(
            r'gene_id "([^"]+)"',
            fields[8]
        )

        if not match:
            continue

        gene = match.group(1)

        if gene in genes:
            models[gene].append({
                'chrom': fields[0],
                'start': int(fields[3]),
                'end': int(fields[4]),
                'strand': fields[6]
            })

def merge(intervals):
    result = []

    for start, end in sorted(intervals):
        if result and start <= result[-1][1] + 1:
            result[-1][1] = max(result[-1][1], end)
        else:
            result.append([start, end])

    return result

def overlap(a, b):
    return sum(
        max(0, min(x2, y2) - max(x1, y1) + 1)
        for x1, x2 in a
        for y1, y2 in b
    )

def cigar_blocks(pos, cigar):
    operations = re.findall(r'(\d+)([MIDNSHP=X])', cigar)

    if ''.join(a+b for a,b in operations) != cigar:
        return [], 0

    refpos = pos
    block_start = refpos
    aligned_query_bases = 0
    blocks = []

    for number, operation in operations:
        number = int(number)

        if operation in ('M', '=', 'X', 'I'):
            aligned_query_bases += number

        if operation in ('M', 'D', '=', 'X'):
            refpos += number

        elif operation == 'N':
            if refpos > block_start:
                blocks.append((block_start, refpos - 1))

            refpos += number
            block_start = refpos

    if refpos > block_start:
        blocks.append((block_start, refpos - 1))

    return blocks, aligned_query_bases

def inferred_strand(flag, tags):
    ts = tags.get('ts')

    if ts not in ('+', '-'):
        return '.'

    alignment_strand = '-' if flag & 16 else '+'

    if ts == '+':
        return alignment_strand

    return '+' if alignment_strand == '-' else '-'

def junction_matches(expected, observed, tolerance):
    return sum(
        any(
            abs(a-c) <= tolerance and abs(b-d) <= tolerance
            for c,d in observed
        )
        for a,b in expected
    )

model_info = {}

for gene in genes:
    records = models[gene]

    if not records:
        raise RuntimeError('Missing CDS model: ' + gene)

    cds = merge([
        (r['start'], r['end']) for r in records
    ])

    model_info[gene] = {
        'chrom': records[0]['chrom'],
        'strand': records[0]['strand'],
        'cds': cds,
        'cds_bp': sum(b-a+1 for a,b in cds),
        'junctions': [
            (cds[i][1], cds[i+1][0])
            for i in range(len(cds)-1)
        ]
    }

candidates = defaultdict(list)

process = subprocess.Popen(
    ['samtools', 'view', bam_file],
    stdout=subprocess.PIPE,
    universal_newlines=True
)

with open(bed_file, 'w') as bed:

    for line in process.stdout:

        f = line.rstrip('\n').split('\t')

        if len(f) < 11:
            continue

        qname = f[0]
        flag = int(f[1])

        # Exclude unmapped and supplementary alignments.
        # Retain secondary alignments to assess paralogs.
        if flag & (4 | 2048):
            continue

        if qname not in metadata:
            continue

        chrom = f[2]
        pos = int(f[3])
        mapq = int(f[4])
        cigar = f[5]

        blocks, aligned_query_bases = cigar_blocks(
            pos, cigar
        )

        if not blocks:
            continue

        tags = {}

        for field in f[11:]:
            parts = field.split(':', 2)
            if len(parts) == 3:
                tags[parts[0]] = parts[2]

        strand = inferred_strand(flag, tags)

        observed_junctions = [
            (blocks[i][1], blocks[i+1][0])
            for i in range(len(blocks)-1)
        ]

        for gene in genes:

            model = model_info[gene]

            if chrom != model['chrom']:
                continue

            cds_overlap = overlap(
                blocks, model['cds']
            )

            if cds_overlap == 0:
                continue

            exact = junction_matches(
                model['junctions'],
                observed_junctions,
                0
            )

            near = junction_matches(
                model['junctions'],
                observed_junctions,
                2
            )

            cds_coverage = (
                100.0 * cds_overlap /
                model['cds_bp']
            )

            query_coverage = (
                100.0 * aligned_query_bases /
                int(metadata[qname]['Length_bp'])
            )

            primary = not bool(flag & 256)

            if strand == '.':
                strand_status = 'UNRESOLVED'
            elif strand == model['strand']:
                strand_status = 'MATCH'
            else:
                strand_status = 'OPPOSITE'

            row = {
                'gene': gene,
                'transcript': qname,
                'original': metadata[qname]['Original_IDs'],
                'source': metadata[qname]['Source_genes'],
                'mapq': mapq,
                'primary': int(primary),
                'strand': strand,
                'strand_status': strand_status,
                'cds_coverage': cds_coverage,
                'query_coverage': query_coverage,
                'exact': exact,
                'near': near,
                'segments': len(blocks),
                'chrom': chrom,
                'start': blocks[0][0],
                'end': blocks[-1][1]
            }

            candidates[gene].append(row)

            # BED12 transcript blocks for IGV inspection
            chrom_start = blocks[0][0] - 1
            chrom_end = blocks[-1][1]

            block_sizes = [
                end-start+1 for start,end in blocks
            ]

            block_starts = [
                start-1-chrom_start
                for start,end in blocks
            ]

            bed.write('\t'.join(map(str, [
                chrom,
                chrom_start,
                chrom_end,
                gene + '__' + qname,
                mapq,
                strand,
                chrom_start,
                chrom_end,
                0,
                len(blocks),
                ','.join(map(str, block_sizes)),
                ','.join(map(str, block_starts))
            ])) + '\n')

returncode = process.wait()

if returncode != 0:
    raise RuntimeError('samtools view failed')

# Save detailed candidate alignments
with open(details_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene', 'Transcript_ID', 'Original_Trinity_IDs',
        'Source_genes', 'Primary_alignment',
        'MAPQ', 'Inferred_transcript_strand',
        'Strand_status', 'CDS_coverage_percent',
        'Transcript_aligned_percent',
        'Exact_junctions', 'Junctions_within_2bp',
        'Aligned_exon_blocks', 'Chromosome',
        'Alignment_start', 'Alignment_end'
    ])

    for gene in genes:
        for x in candidates[gene]:
            writer.writerow([
                gene, x['transcript'], x['original'],
                x['source'], x['primary'], x['mapq'],
                x['strand'], x['strand_status'],
                round(x['cds_coverage'], 2),
                round(x['query_coverage'], 2),
                x['exact'], x['near'],
                x['segments'], x['chrom'],
                x['start'], x['end']
            ])

# Summarize the most promising candidates
summary = []

for gene in genes:

    model = model_info[gene]
    rows = candidates[gene]

    # Prefer same-strand, primary, well-mapped transcripts.
    # Unknown-strand alignments remain visible for review.
    ranked = sorted(
        rows,
        key=lambda x: (
            x['strand_status'] == 'MATCH',
            x['primary'],
            x['mapq'],
            x['exact'],
            x['near'],
            x['cds_coverage']
        ),
        reverse=True
    )

    best = ranked[0] if ranked else None

    summary.append([
        gene,
        model['chrom'],
        model['strand'],
        model['cds_bp'],
        len(model['junctions']),
        len(set(x['transcript'] for x in rows)),
        len(rows),
        sum(x['strand_status'] == 'MATCH' for x in rows),
        round(max(
            [x['cds_coverage'] for x in rows] or [0]
        ), 2),
        best['transcript'] if best else '.',
        best['strand_status'] if best else '.',
        best['primary'] if best else '.',
        best['mapq'] if best else '.',
        round(best['cds_coverage'], 2) if best else 0,
        best['exact'] if best else 0,
        best['near'] if best else 0,
        'CANDIDATES_FOR_REVIEW' if best else 'NONE_DETECTED'
    ])

with open(summary_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene', 'Chromosome', 'Strand', 'CDS_bp',
        'Expected_junctions', 'Distinct_transcripts',
        'Candidate_alignments', 'Matching_strand_alignments',
        'Maximum_CDS_coverage_percent',
        'Best_transcript', 'Best_strand_status',
        'Best_is_primary', 'Best_MAPQ',
        'Best_CDS_coverage_percent',
        'Best_exact_junctions',
        'Best_junctions_within_2bp',
        'Review_status'
    ])

    writer.writerows(summary)

print('\n===== GENE-MODEL REVIEW SUMMARY =====')
print('Gene\tTranscripts\tBest_CDS_%\tExact\tWithin_2bp\tMAPQ\tStrand')

for row in summary:
    print('{}\t{}\t{}\t{}/{}\t{}/{}\t{}\t{}'.format(
        row[0], row[5], row[13],
        row[14], row[4],
        row[15], row[4],
        row[12], row[10]
    ))

print('\nSummary:', summary_file)
print('Details:', details_file)
print('IGV BED12:', bed_file)
PY

echo
echo "===== COMPLETED ====="
echo "Status: PASS"