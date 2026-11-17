#!/bin/bash
#SBATCH --job-name=toxin_residue_validation
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

# fix shared-library compatibility
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# set paths
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
OUTDIR="${WORKDIR}/protein_validation/catalytic_residues"
REFDIR="${OUTDIR}/references"

PEP="${WORKDIR}/toxin_annotation_pep.fasta"
PFAM="${WORKDIR}/protein_validation/pfam_targeted/toxin_domain_validation_details.tsv"

SUMMARY="${OUTDIR}/toxin_catalytic_summary.tsv"
DETAILS="${OUTDIR}/toxin_catalytic_residue_details.tsv"

mkdir -p "$REFDIR"

# ------------------------------------------------------------
# STEP 1: validate environment
# ------------------------------------------------------------

echo "===== STEP 1: ENVIRONMENT ====="

[[ -s "$PEP" ]] || {
    echo "ERROR: Missing protein FASTA: $PEP" >&2
    exit 1
}

python -c 'from Bio import pairwise2; import json; print("Biopython: PASS")'

command -v curl >/dev/null || {
    echo "ERROR: curl not found" >&2
    exit 1
}

# ------------------------------------------------------------
# STEP 2: obtain reference sequences and annotations
#
# cached JSON files are reused automatically.
# ------------------------------------------------------------

echo
echo "===== STEP 2: REFERENCE PROTEINS ====="

for accession in P04971 P00623 J3S820; do

    ref="${REFDIR}/${accession}.json"

    if [[ -s "$ref" ]] && python - "$ref" "$accession" <<'PY'
import sys, json
data = json.load(open(sys.argv[1]))
assert data['primaryAccession'] == sys.argv[2]
assert len(data['sequence']['value']) > 100
PY
    then
        echo "Cached: $accession"
        continue
    fi

    echo "Downloading: $accession"

    tmp="${ref}.tmp"
    rm -f "$tmp"

    curl -fLsS \
        --retry 2 \
        --connect-timeout 10 \
        --max-time 45 \
        "https://rest.uniprot.org/uniprotkb/${accession}.json" \
        -o "$tmp" || {
            echo "ERROR: Reference download failed: $accession" >&2
            echo "Run this script with --download-only on the login node." >&2
            exit 1
        }

    python - "$tmp" "$accession" <<'PY'
import sys, json
data = json.load(open(sys.argv[1]))
assert data['primaryAccession'] == sys.argv[2]
assert len(data['sequence']['value']) > 100
print('Validated:', sys.argv[2])
PY

    mv "$tmp" "$ref"
done

if [[ "${1:-}" == '--download-only' ]]; then
    echo "All reference records are cached."
    exit 0
fi

# ------------------------------------------------------------
# STEP 3: reference-guided residue comparison
# ------------------------------------------------------------

echo
echo "===== STEP 3: CATALYTIC RESIDUE ANALYSIS ====="

python - "$PEP" "$REFDIR" "$PFAM" "$SUMMARY" "$DETAILS" <<'PY'
import sys
import os
import re
import csv
import json
import warnings
from collections import defaultdict

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import pairwise2

pep_file, refdir, pfam_file, summary_file, details_file = sys.argv[1:]

genes = [
    'SVSP-1', 'SVSP-16', 'Ficolin-4',
    'PLA2-4', 'CTL-10', 'HYAL-4',
    'CTL-15', 'CTL-16', 'CTL-17', 'CTL-18'
]

enzyme_references = {
    'SVSP-1': 'P04971',
    'SVSP-16': 'P04971',
    'PLA2-4': 'P00623',
    'HYAL-4': 'J3S820'
}

# ------------------------------------------------------------
# Read protein FASTA
# ------------------------------------------------------------

records = []
header = None
parts = []

with open(pep_file) as f:
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

proteins = {}

for gene in genes:
    pattern = re.compile(
        r'(?<![A-Za-z0-9])' +
        re.escape(gene) +
        r'(?![A-Za-z0-9])'
    )

    matches = [
        s for h, s in records if pattern.search(h)
    ]

    if len(matches) != 1:
        raise RuntimeError(
            '{}: expected one protein, found {}'.format(
                gene, len(matches)
            )
        )

    original = matches[0]

    # Remove one annotated terminal stop only.
    proteins[gene] = (
        original[:-1] if original.endswith('*') else original
    )

# ------------------------------------------------------------
# Read previously validated Pfam domain positions
# ------------------------------------------------------------

domain_hits = defaultdict(list)

if os.path.isfile(pfam_file):
    with open(pfam_file) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            if row['Gene'] in genes:
                if float(row['Domain_iEvalue']) <= 1e-3:
                    domain_hits[row['Gene']].append(row)

# ------------------------------------------------------------
# Read reference records
# ------------------------------------------------------------

references = {}

for accession in sorted(set(enzyme_references.values())):
    with open(os.path.join(refdir, accession + '.json')) as f:
        record = json.load(f)

    references[accession] = record

def feature_position(feature, endpoint):
    value = feature.get('location', {}).get(endpoint, {}).get('value')
    return int(value) if value is not None else None

def mature_start(record):
    features = record.get('features', [])

    chains = [
        feature_position(f, 'start')
        for f in features
        if f.get('type') == 'Chain'
        and feature_position(f, 'start') is not None
    ]

    if chains:
        return min(chains)

    signals = [
        feature_position(f, 'end')
        for f in features
        if f.get('type') in ('Signal', 'Signal peptide')
        and feature_position(f, 'end') is not None
    ]

    if signals:
        return max(signals) + 1

    return None

def annotated_sites(record):
    seq = record['sequence']['value'].upper()
    sites = []
    seen = set()

    for feature in record.get('features', []):
        if feature.get('type') not in ('Active site', 'Binding site'):
            continue

        start = feature_position(feature, 'start')
        end = feature_position(feature, 'end')

        if start is None or start != end:
            continue

        if not 1 <= start <= len(seq):
            continue

        if start in seen:
            continue

        seen.add(start)

        sites.append({
            'label': feature['type'].replace(' ', '_') + '_' + str(start),
            'position': start,
            'expected': seq[start - 1],
            'source': 'UniProt_feature'
        })

    return sites

def canonical_sites(record, numbered_sites, max_shift=5):
    seq = record['sequence']['value'].upper()
    start = mature_start(record)

    if start is None:
        return []

    candidates = []

    for shift in range(-max_shift, max_shift + 1):
        positions = [
            start + number - 1 + shift
            for label, number, aa in numbered_sites
        ]

        if not all(1 <= p <= len(seq) for p in positions):
            continue

        if all(
            seq[p - 1] == aa
            for p, (label, number, aa)
            in zip(positions, numbered_sites)
        ):
            candidates.append(positions)

    # Never force uncertain reference numbering.
    if len(candidates) != 1:
        return []

    sites = []

    for (label, number, aa), position in zip(
        numbered_sites, candidates[0]
    ):
        sites.append({
            'label': label,
            'position': position,
            'expected': aa,
            'source': 'validated_reference_numbering'
        })

    return sites

def reference_sites(accession):
    record = references[accession]
    seq = record['sequence']['value'].upper()

    if accession == 'P04971':
        # Verified batroxobin precursor positions (MEROPS).
        positions = [
            ('Catalytic_His', 65, 'H'),
            ('Catalytic_Asp', 110, 'D'),
            ('Catalytic_Ser', 202, 'S')
        ]

        for label, pos, aa in positions:
            if seq[pos - 1] != aa:
                raise RuntimeError(
                    'Batroxobin reference residue mismatch: ' + label
                )

        return [
            {
                'label': label,
                'position': pos,
                'expected': aa,
                'source': 'MEROPS_P04971'
            }
            for label, pos, aa in positions
        ]

    if accession == 'P00623':
        sites = canonical_sites(record, [
            ('Catalytic_His48', 48, 'H'),
            ('Calcium_Asp49', 49, 'D'),
            ('Catalytic_network_Tyr52', 52, 'Y'),
            ('Catalytic_network_Asp99', 99, 'D')
        ])

        return sites if sites else annotated_sites(record)

    if accession == 'J3S820':
        # Snake hyaluronidase numbering is validated against
        # the downloaded reference before mapping.
        sites = canonical_sites(record, [
            ('Substrate_Asp110', 110, 'D'),
            ('Catalytic_Glu112', 112, 'E'),
            ('Substrate_Tyr183', 183, 'Y'),
            ('Substrate_Tyr230', 230, 'Y'),
            ('Substrate_Trp315', 315, 'W')
        ])

        return sites if sites else annotated_sites(record)

    return []

def align_and_map(reference, query):
    # Keep stop positions as X so original protein
    # coordinates remain unchanged.
    query_for_alignment = query.replace('*', 'X')

    aln = pairwise2.align.globalms(
        reference,
        query_for_alignment,
        2, -1, -8, -0.5,
        one_alignment_only=True
    )[0]

    aligned_ref, aligned_query = aln[0], aln[1]

    mapping = {}
    refpos = 0
    querypos = 0
    matches = 0
    compared = 0

    for a, b in zip(aligned_ref, aligned_query):
        if a != '-':
            refpos += 1
        if b != '-':
            querypos += 1

        if a != '-':
            mapping[refpos] = querypos if b != '-' else None

        if a != '-' and b not in ('-', 'X'):
            compared += 1
            if a == b:
                matches += 1

    identity = (
        100.0 * matches / compared
        if compared else 0.0
    )

    return mapping, identity

# ------------------------------------------------------------
# Evaluate all genes
# ------------------------------------------------------------

summary_rows = []
detail_rows = []

for gene in genes:
    seq = proteins[gene]

    first_stop = seq.find('*')
    first_stop = first_stop + 1 if first_stop >= 0 else None

    internal_stops = seq.count('*')
    amino_acids = sum(aa != '*' for aa in seq)

    best_domain = None
    if domain_hits[gene]:
        best_domain = min(
            domain_hits[gene],
            key=lambda x: float(x['Domain_iEvalue'])
        )

    domain_cys = 'NA'

    if best_domain:
        a = int(best_domain['Protein_start'])
        b = int(best_domain['Protein_end'])
        domain_cys = seq[a-1:b].count('C')

    if gene not in enzyme_references:
        summary_rows.append([
            gene, amino_acids, internal_stops,
            first_stop if first_stop else 'NA',
            domain_cys, 'NA', 'NA',
            0, 0, 0, 0, 0,
            'STRUCTURAL_DOMAIN_ONLY'
        ])
        continue

    accession = enzyme_references[gene]
    reference = references[accession]['sequence']['value'].upper()

    mapping, identity = align_and_map(reference, seq)
    sites = reference_sites(accession)

    counts = defaultdict(int)

    for site in sites:
        pos = site['position']
        query_pos = mapping.get(pos)

        if query_pos is None:
            observed = '-'
            status = 'UNMAPPED'
        else:
            observed = seq[query_pos - 1]

            if observed == '*':
                status = 'AT_STOP'
            elif first_stop and query_pos > first_stop:
                status = 'DOWNSTREAM_OF_FIRST_STOP'
            elif observed == site['expected']:
                status = 'CONSERVED_BEFORE_STOP'
            else:
                status = 'ALTERED_BEFORE_STOP'

        counts[status] += 1

        detail_rows.append([
            gene,
            accession,
            site['label'],
            site['source'],
            pos,
            site['expected'],
            query_pos if query_pos is not None else 'NA',
            observed,
            first_stop if first_stop else 'NA',
            status,
            round(identity, 2)
        ])

    if not sites:
        interpretation = 'REFERENCE_SITE_MAPPING_UNRESOLVED'
    elif counts['DOWNSTREAM_OF_FIRST_STOP'] or counts['AT_STOP']:
        interpretation = 'SITES_DOWNSTREAM_OF_STOP'
    elif counts['ALTERED_BEFORE_STOP']:
        interpretation = 'CATALYTIC_SITE_ALTERED'
    elif counts['UNMAPPED']:
        interpretation = 'INCOMPLETE_SITE_MAPPING'
    elif counts['CONSERVED_BEFORE_STOP'] == len(sites):
        interpretation = 'ALL_TESTED_SITES_CONSERVED'
    else:
        interpretation = 'REVIEW'

    summary_rows.append([
        gene,
        amino_acids,
        internal_stops,
        first_stop if first_stop else 'NA',
        domain_cys,
        accession,
        round(identity, 2),
        len(sites),
        counts['CONSERVED_BEFORE_STOP'],
        counts['ALTERED_BEFORE_STOP'],
        counts['DOWNSTREAM_OF_FIRST_STOP'] + counts['AT_STOP'],
        counts['UNMAPPED'],
        interpretation
    ])

# ------------------------------------------------------------
# Write output files
# ------------------------------------------------------------

with open(summary_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene',
        'Residues_aa',
        'Internal_stops',
        'First_stop_position',
        'Pfam_domain_Cys',
        'Reference_accession',
        'Alignment_identity_percent',
        'Sites_checked',
        'Conserved_before_stop',
        'Altered_before_stop',
        'Sites_at_or_after_stop',
        'Unmapped_sites',
        'Site_assessment'
    ])

    writer.writerows(summary_rows)

with open(details_file, 'w', newline='') as f:
    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene',
        'Reference_accession',
        'Site',
        'Site_source',
        'Reference_position',
        'Expected_residue',
        'ToxCodAn_position',
        'Observed_residue',
        'First_stop_position',
        'Status',
        'Alignment_identity_percent'
    ])

    writer.writerows(detail_rows)

print('\n===== CATALYTIC RESIDUE SUMMARY =====')
print('Gene\tSites\tConserved\tAltered\tAfter_stop\tUnmapped\tAssessment')

for row in summary_rows:
    print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
        row[0], row[7], row[8],
        row[9], row[10], row[11], row[12]
    ))

print('\nSummary:', summary_file)
print('Details:', details_file)
PY

echo
echo "===== COMPLETED ====="
echo "Status: PASS"