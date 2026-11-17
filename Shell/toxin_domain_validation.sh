#!/bin/bash
#SBATCH --job-name=toxin_domain_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
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
# paths
# ------------------------------------------------------------
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
OUTDIR="${WORKDIR}/protein_validation/pfam_targeted"
PFAMDIR="${OUTDIR}/profiles"

PEP="${WORKDIR}/toxin_annotation_pep.fasta"
INPUT="${OUTDIR}/toxin_ORF_fragments.fasta"
METADATA="${OUTDIR}/toxin_ORF_fragments.tsv"

HMM="${OUTDIR}/toxin_families.hmm"
DOMTBLOUT="${OUTDIR}/toxin_domains.domtblout"
SUMMARY="${OUTDIR}/toxin_domain_validation_summary.tsv"
DETAILS="${OUTDIR}/toxin_domain_validation_details.tsv"

mkdir -p "$PFAMDIR"

# verify tools
for tool in python hmmscan hmmpress curl gzip; do
    command -v "$tool" >/dev/null || {
        echo "ERROR: Missing command: $tool" >&2
        exit 1
    }
done

# ------------------------------------------------------------
# STEP 1: obtain five Pfam HMM profiles
#
# reuse cached profiles if present
# only five small models are downloaded
# ------------------------------------------------------------

echo "===== STEP 1: PFAM PROFILES ====="

PFAM_IDS=(
    PF00089
    PF00147
    PF00068
    PF00059
    PF01630
)

for pfam in "${PFAM_IDS[@]}"; do

    profile="${PFAMDIR}/${pfam}.hmm"

    if [[ -s "$profile" ]] && grep -q '^HMMER3/' "$profile"; then
        echo "Cached: $pfam"
        continue
    fi

    echo "Downloading: $pfam"

    url="https://www.ebi.ac.uk/interpro/api/entry/pfam/${pfam}?annotation=hmm"
    tmp="${profile}.download"

    rm -f "$tmp" "${profile}.tmp"

    curl \
        --fail \
        --location \
        --silent \
        --show-error \
        --retry 2 \
        --connect-timeout 10 \
        --max-time 60 \
        "$url" \
        -o "$tmp" || {
            echo "ERROR: Could not download $pfam" >&2
            echo "If compute-node internet is blocked, run this script with --download-only on the login node." >&2
            exit 1
        }

    if gzip -t "$tmp" 2>/dev/null; then
        gzip -dc "$tmp" > "${profile}.tmp"
    else
        cp "$tmp" "${profile}.tmp"
    fi

    if ! grep -q '^HMMER3/' "${profile}.tmp"; then
        echo "ERROR: Downloaded file is not a valid HMM: $pfam" >&2
        exit 1
    fi

    mv "${profile}.tmp" "$profile"
    rm -f "$tmp"

done

# optional: download profiles from login node without running analyses
if [[ "${1:-}" == '--download-only' ]]; then
    echo "All five Pfam profiles cached."
    exit 0
fi

# combine and index HMM profiles
: > "$HMM"

for pfam in "${PFAM_IDS[@]}"; do
    cat "${PFAMDIR}/${pfam}.hmm" >> "$HMM"
done

hmmpress -f "$HMM" > "${OUTDIR}/hmmpress.log"


# ------------------------------------------------------------
# STEP 2: prepare all 10 protein models
#
# split sequences at stop codons
# keep original amino-acid coordinates
# no artificial concatenation across premature stops
# ------------------------------------------------------------

echo
echo "===== STEP 2: PREPARE PROTEIN FRAGMENTS ====="

python - "$PEP" "$INPUT" "$METADATA" <<'PY'
import sys
import re
import csv

source, output_fasta, output_metadata = sys.argv[1:]

genes = [
    'SVSP-1', 'SVSP-16', 'Ficolin-4',
    'PLA2-4', 'CTL-10', 'HYAL-4',
    'CTL-15', 'CTL-16', 'CTL-17', 'CTL-18'
]

records = []
header = None
parts = []

with open(source) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                records.append((header, ''.join(parts).upper()))
            header = line[1:]
            parts = []
        else:
            parts.append(line)

if header is not None:
    records.append((header, ''.join(parts).upper()))

with open(output_fasta, 'w') as fasta, \
     open(output_metadata, 'w', newline='') as meta:

    writer = csv.writer(meta, delimiter='\t')

    writer.writerow([
        'Fragment_ID',
        'Gene',
        'Protein_start',
        'Protein_end',
        'Fragment_length',
        'Internal_stops',
        'First_stop_position'
    ])

    for gene_index, gene in enumerate(genes, 1):

        pattern = re.compile(
            r'(?<![A-Za-z0-9])' +
            re.escape(gene) +
            r'(?![A-Za-z0-9])'
        )

        matches = [
            seq for name, seq in records
            if pattern.search(name)
        ]

        if len(matches) != 1:
            raise RuntimeError(
                '{}: expected one protein, found {}'.format(
                    gene, len(matches)
                )
            )

        original = matches[0]

        # Remove only the annotated final stop.
        seq = (
            original[:-1]
            if original.endswith('*')
            else original
        )

        stops = seq.count('*')
        first_stop = (
            seq.index('*') + 1
            if stops else None
        )

        fragments = seq.split('*')
        cursor = 1
        fragment_count = 0

        for fragment_index, fragment in enumerate(fragments, 1):

            start = cursor
            end = start + len(fragment) - 1

            if fragment:
                fragment_id = 'T{:02d}F{:02d}'.format(
                    gene_index, fragment_index
                )

                fasta.write(
                    '>{}\n{}\n'.format(fragment_id, fragment)
                )

                writer.writerow([
                    fragment_id,
                    gene,
                    start,
                    end,
                    len(fragment),
                    stops,
                    first_stop if first_stop else 'NA'
                ])

                fragment_count += 1

            cursor = end + 2

        print(
            '{}\tstops={}\tfragments={}'.format(
                gene, stops, fragment_count
            )
        )

print('Prepared all 10 toxin gene models.')
PY


# ------------------------------------------------------------
# STEP 3: HMMER domain analysis
# ------------------------------------------------------------

echo
echo "===== STEP 3: HMMER DOMAIN SEARCH ====="

hmmscan \
    --cpu "${SLURM_CPUS_PER_TASK:-4}" \
    --noali \
    -E 10 \
    --domE 10 \
    --incE 10 \
    --incdomE 10 \
    --domtblout "$DOMTBLOUT" \
    "$HMM" \
    "$INPUT" \
    > "${OUTDIR}/hmmscan.log"

[[ -s "$DOMTBLOUT" ]]


# ------------------------------------------------------------
# STEP 4: integrate protein-domain results
# ------------------------------------------------------------

echo
echo "===== STEP 4: DOMAIN VALIDATION SUMMARY ====="

python - "$METADATA" "$DOMTBLOUT" "$SUMMARY" "$DETAILS" <<'PY'
import sys
import csv
from collections import defaultdict

metadata_file, domtbl_file, summary_file, details_file = sys.argv[1:]

genes = [
    'SVSP-1', 'SVSP-16', 'Ficolin-4',
    'PLA2-4', 'CTL-10', 'HYAL-4',
    'CTL-15', 'CTL-16', 'CTL-17', 'CTL-18'
]

expected = {
    'SVSP-1': 'PF00089',
    'SVSP-16': 'PF00089',
    'Ficolin-4': 'PF00147',
    'PLA2-4': 'PF00068',
    'CTL-10': 'PF00059',
    'HYAL-4': 'PF01630',
    'CTL-15': 'PF00059',
    'CTL-16': 'PF00059',
    'CTL-17': 'PF00059',
    'CTL-18': 'PF00059'
}

metadata = {}
gene_metadata = {}

with open(metadata_file) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        metadata[row['Fragment_ID']] = row
        gene_metadata[row['Gene']] = row

hits = defaultdict(list)
details = []

with open(domtbl_file) as f:
    for line in f:

        if line.startswith('#') or not line.strip():
            continue

        fields = line.rstrip('\n').split(None, 22)

        if len(fields) < 22:
            raise RuntimeError('Malformed HMMER domain table row')

        # HMMER hmmscan domtblout format
        pfam_name = fields[0]
        pfam_id = fields[1].split('.')[0]
        hmm_length = int(fields[2])

        fragment_id = fields[3]

        if fragment_id not in metadata:
            raise RuntimeError(
                'Unknown fragment: ' + fragment_id
            )

        m = metadata[fragment_id]
        gene = m['Gene']

        evalue = float(fields[12])

        hmm_start = int(fields[15])
        hmm_end = int(fields[16])
        aa_start = int(fields[17])
        aa_end = int(fields[18])

        model_coverage = (
            100.0 *
            (hmm_end - hmm_start + 1) /
            hmm_length
        )

        # Map fragment coordinates to original protein.
        protein_start = (
            int(m['Protein_start']) + aa_start - 1
        )
        protein_end = (
            int(m['Protein_start']) + aa_end - 1
        )

        stop = m['First_stop_position']

        if stop == 'NA':
            stop_relation = 'no_internal_stop'
        elif protein_end < int(stop):
            stop_relation = 'before_first_stop'
        else:
            stop_relation = 'after_first_stop'

        hit = {
            'gene': gene,
            'pfam': pfam_id,
            'name': pfam_name,
            'evalue': evalue,
            'coverage': model_coverage,
            'start': protein_start,
            'end': protein_end,
            'fragment': fragment_id,
            'relation': stop_relation
        }

        hits[gene].append(hit)

        details.append([
            gene,
            fragment_id,
            pfam_id,
            pfam_name,
            '{:.3g}'.format(evalue),
            round(model_coverage, 2),
            protein_start,
            protein_end,
            stop_relation,
            int(pfam_id == expected[gene]),
            int(evalue <= 1e-3)
        ])

# Full domain-hit table
with open(details_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene',
        'Fragment_ID',
        'Pfam_ID',
        'Pfam_name',
        'Domain_iEvalue',
        'HMM_coverage_percent',
        'Protein_start',
        'Protein_end',
        'Stop_relation',
        'Expected_family',
        'Evalue_le_1e-3'
    ])

    writer.writerows(details)

# Gene-level summary
summary = []

for gene in genes:

    gene_hits = [
        hit for hit in hits[gene]
        if hit['pfam'] == expected[gene]
    ]

    gene_hits.sort(key=lambda x: x['evalue'])

    best = gene_hits[0] if gene_hits else None

    metadata_row = gene_metadata[gene]

    if best:
        supported = (
            'Yes' if best['evalue'] <= 1e-3 else 'No'
        )
    else:
        supported = 'No'

    summary.append([
        gene,
        expected[gene],
        metadata_row['Internal_stops'],
        metadata_row['First_stop_position'],
        best['name'] if best else '.',
        '{:.3g}'.format(best['evalue']) if best else '.',
        round(best['coverage'], 2) if best else 0,
        best['start'] if best else '.',
        best['end'] if best else '.',
        best['relation'] if best else '.',
        supported
    ])

with open(summary_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene',
        'Expected_Pfam',
        'Internal_stops',
        'First_stop_position',
        'Best_domain',
        'Best_domain_iEvalue',
        'HMM_coverage_percent',
        'Protein_start',
        'Protein_end',
        'Stop_relation',
        'Expected_domain_supported'
    ])

    writer.writerows(summary)

print(
    'Gene\tExpected_domain\tEvalue\t'
    'HMM_coverage_%\tStop_relation\tSupported'
)

for row in summary:
    print('{}\t{}\t{}\t{}\t{}\t{}'.format(
        row[0], row[4], row[5], row[6], row[9], row[10]
    ))

print('\nSummary:', summary_file)
print('Details:', details_file)
PY

echo
echo "===== COMPLETED ====="
echo "Status: PASS"