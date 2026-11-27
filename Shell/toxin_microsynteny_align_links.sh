#!/usr/bin/env bash

# activate conda env
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

set -euo pipefail

# ============================================================
# select informative local alignments for microsynteny
# ============================================================

# 1. directories
ROOT="${ROOT:-$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo}"
SYNTENY="$ROOT/synteny/toxin_microsynteny"

GENES="$SYNTENY/05_annotations/comparative_all_gene_models_in_regions.tsv"
ALIGNMENTS="$SYNTENY/06_local_alignments/filtered_blocks_for_plot.tsv"
OUTPUT="$SYNTENY/06_local_alignments/selected_links_for_plot.tsv"

# 2. verify inputs
command -v python >/dev/null || {
    echo "ERROR: python not found" >&2
    exit 1
}

for FILE in "$GENES" "$ALIGNMENTS"; do
    [[ -s "$FILE" ]] || {
        echo "ERROR: missing $FILE" >&2
        exit 1
    }
done

# 3. select gene-supported alignment links
python - "$GENES" "$ALIGNMENTS" "$OUTPUT" <<'PY'
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

genes_file, alignments_file, output_file = map(Path, sys.argv[1:])

# Maximum toxin links per comparison, for visualization only
MAX_TOXIN_LINKS = {'SVMP': 25, 'PLA2': 12}

ANCHOR_ROLES = {
    'conserved_anchor_candidate',
    'provisional_VWA5B1_anchor_group'
}

TOXIN_ROLES = {
    'putative_venom_toxin',
    'SVMP_like_candidate',
    'PLA2_like_candidate'
}


def read_tsv(path, expected):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = expected - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(
                f'ERROR: {path.name} missing columns: {sorted(missing)}'
            )
        return list(reader), reader.fieldnames


gene_rows, _ = read_tsv(genes_file, {
    'species', 'family', 'chr', 'start', 'end', 'gene_id',
    'plot_role', 'plot_include', 'anchor_symbol'
})

blocks, block_fields = read_tsv(alignments_file, {
    'family', 'reference', 'query', 'reference_chr', 'query_chr',
    'reference_genomic_start', 'reference_genomic_end',
    'query_genomic_start', 'query_genomic_end', 'bitscore'
})

# Index candidate anchors and toxin genes
index = defaultdict(list)

for g in gene_rows:
    if g['plot_include'] != '1':
        continue
    if g['plot_role'] not in ANCHOR_ROLES | TOXIN_ROLES:
        continue

    g['start'] = int(g['start'])
    g['end'] = int(g['end'])

    index[(g['species'], g['family'], g['chr'])].append(g)


def midpoint_gene(species, family, chrom, begin, end):
    """Assign an alignment to one gene containing its midpoint."""
    midpoint = (int(begin) + int(end)) / 2

    matches = [
        g for g in index[(species, family, chrom)]
        if g['start'] <= midpoint <= g['end']
    ]

    return matches[0] if len(matches) == 1 else None


# Retain the strongest HSP for each candidate gene pair
best = {}

for block in blocks:
    family = block['family']

    ref_gene = midpoint_gene(
        block['reference'], family, block['reference_chr'],
        block['reference_genomic_start'],
        block['reference_genomic_end']
    )

    query_gene = midpoint_gene(
        block['query'], family, block['query_chr'],
        block['query_genomic_start'],
        block['query_genomic_end']
    )

    if ref_gene is None or query_gene is None:
        continue

    ref_role = ref_gene['plot_role']
    query_role = query_gene['plot_role']

    # Conserved anchor candidates
    if ref_role in ANCHOR_ROLES and query_role in ANCHOR_ROLES:

        shared = (
            set(ref_gene['anchor_symbol'].split(';')) &
            set(query_gene['anchor_symbol'].split(';'))
        )
        shared.discard('')

        if not shared:
            continue

        link_type = 'same_named_anchor_candidate'
        symbol = sorted(shared)[0]

        identity = (
            'anchor', family, block['query'], symbol
        )

    # Candidate toxin-family homologs
    elif ref_role in TOXIN_ROLES and query_role in TOXIN_ROLES:

        link_type = 'toxin_family_homology_candidate'
        symbol = ''

        identity = (
            'toxin', family, block['query'],
            ref_gene['gene_id'], query_gene['gene_id']
        )

    else:
        continue

    candidate = dict(
        block,
        link_type=link_type,
        reference_gene=ref_gene['gene_id'],
        query_gene=query_gene['gene_id'],
        anchor_symbol=symbol
    )

    if (
        identity not in best
        or float(candidate['bitscore']) >
        float(best[identity]['bitscore'])
    ):
        best[identity] = candidate


# Reduce toxin-link density without removing anchors
groups = defaultdict(list)

for record in best.values():
    key = (
        record['family'],
        record['query'],
        record['link_type']
    )
    groups[key].append(record)

selected = []

for (family, query, link_type), matches in sorted(groups.items()):

    matches.sort(
        key=lambda r: float(r['bitscore']),
        reverse=True
    )

    if link_type == 'toxin_family_homology_candidate':
        matches = matches[:MAX_TOXIN_LINKS[family]]

    selected.extend(matches)


# Export selected alignment links
outfields = block_fields + [
    'link_type',
    'reference_gene',
    'query_gene',
    'anchor_symbol'
]

with output_file.open('w', newline='') as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=outfields,
        delimiter='\t',
        lineterminator='\n'
    )
    writer.writeheader()
    writer.writerows(selected)


# 4. summarize results
print('===== STEP 9: SELECTED MICROSYNTENY LINKS =====')

print(f'Input filtered HSPs: {len(blocks):,}')
print(f'Selected visualization links: {len(selected):,}')

counts = Counter(
    (r['family'], r['query'], r['link_type'])
    for r in selected
)

for family in ('SVMP', 'PLA2'):
    for species in (
        'Cerastes_gasperettii',
        'Crotalus_adamanteus',
        'Naja_naja'
    ):

        anchors = counts[
            (family, species, 'same_named_anchor_candidate')
        ]

        toxins = counts[
            (family, species, 'toxin_family_homology_candidate')
        ]

        print(
            f'{family:<5} {species:<23} '
            f'anchor_links={anchors:>2} '
            f'toxin_links={toxins:>2}'
        )

        if anchors == 0:
            print(
                '  NOTE: no same-named anchor HSP passed '
                'the gene-midpoint rule'
            )

print(f'Output: {output_file}')
print(
    'NOTE: links are selected local HSPs, '
    'not validated orthologs or chained syntenic blocks.'
)
print('===== COMPLETE =====')
PY