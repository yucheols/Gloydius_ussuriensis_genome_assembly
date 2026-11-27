#!/usr/bin/env bash

# activate conda env
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

set -euo pipefail

# ============================================================
# local genomic alignments of SVMP and PLA2 regions
# ============================================================

# 1. directories
ROOT="${ROOT:-$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo}"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
MANIFEST="$SYNTENY/01_reference_inventory/anchor_bounded_regions.tsv"
OUTDIR="$SYNTENY/06_local_alignments"
mkdir -p "$OUTDIR/blastn"

# 2. verify genomic intervals and prepare comparisons
for CMD in python blastn; do
    command -v "$CMD" >/dev/null || { echo "ERROR: $CMD is required" >&2; exit 1; }
done
[[ -s "$MANIFEST" ]] || { echo "ERROR: missing $MANIFEST" >&2; exit 1; }

python - "$MANIFEST" "$OUTDIR" <<'PY'
import csv
import sys
from pathlib import Path

manifest, outdir = map(Path, sys.argv[1:])
with manifest.open(newline='') as handle:
    data = list(csv.DictReader(handle, delimiter='\t'))
required = {'species', 'family', 'chromosome', 'start', 'end', 'fasta'}
if not data or not required.issubset(data[0]):
    raise SystemExit('ERROR: invalid Step 5 manifest columns')
regions = {}
for r in data:
    key = (r['family'], r['species'])
    if key in regions:
        raise SystemExit(f'ERROR: duplicate interval: {key}')
    path = Path(r['fasta'])
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f'ERROR: missing FASTA: {path}')
    sequence_id, length = None, 0
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if sequence_id is not None:
                    raise SystemExit(f'ERROR: multiple sequences in {path}')
                sequence_id = line[1:].split()[0]
            else:
                length += len(line.strip())
    expected = int(r['end']) - int(r['start']) + 1
    if sequence_id is None or length != expected:
        raise SystemExit(f'ERROR: unexpected FASTA length for {key}: {length} != {expected}')
    r['sequence_id'], r['length_bp'] = sequence_id, length
    regions[key] = r

species = ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja')
pairs = []
for family in ('SVMP', 'PLA2'):
    ref = regions.get((family, 'Gloydius_ussuriensis'))
    if ref is None:
        raise SystemExit(f'ERROR: no Gloydius {family} interval')
    for sp in species:
        query = regions.get((family, sp))
        if query is None:
            raise SystemExit(f'ERROR: missing {sp} {family} interval')
        pairs.append(dict(family=family, reference='Gloydius_ussuriensis', query=sp,
                          reference_fasta=ref['fasta'], query_fasta=query['fasta']))

with (outdir / 'alignment_pairs.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(pairs[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(pairs)
with (outdir / 'interval_metadata.tsv').open('w', newline='') as handle:
    fields = ['species', 'family', 'chromosome', 'start', 'end', 'fasta', 'sequence_id', 'length_bp']
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader(); writer.writerows(data)
print(f'Validated {len(regions)} FASTAs; prepared {len(pairs)} Gloydius-centered comparisons')
PY

# 3. run BLASTN local alignments
while IFS=$'\t' read -r FAMILY REF QUERY REF_FASTA QUERY_FASTA; do
    [[ "$FAMILY" == family ]] && continue
    NAME="${FAMILY}_${REF}__${QUERY}"
    echo "===== $NAME ====="
    blastn -task blastn -query "$QUERY_FASTA" -subject "$REF_FASTA" \
        -dust yes -word_size 11 -evalue 1e-10 \
        -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send pident length nident evalue bitscore' \
        -out "$OUTDIR/blastn/${NAME}.tsv"
    echo "HSPs: $(wc -l < "$OUTDIR/blastn/${NAME}.tsv")"
done < "$OUTDIR/alignment_pairs.tsv"

# 4. export alignment evidence and provisional plotting blocks
python - "$OUTDIR" <<'PY'
import csv
import sys
from collections import Counter
from pathlib import Path

outdir = Path(sys.argv[1])
def load(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
meta = {(r['family'], r['species']): r for r in load(outdir / 'interval_metadata.tsv')}
pairs = load(outdir / 'alignment_pairs.tsv')
fields = ['family', 'reference', 'query', 'reference_chr', 'query_chr',
          'start', 'end', 'start2', 'end2', 'strand', 'aligned_bp', 'matching_bp',
          'percent_identity', 'evalue', 'bitscore', 'reference_genomic_start',
          'reference_genomic_end', 'query_genomic_start', 'query_genomic_end']
all_rows, filtered = [], []
for pair in pairs:
    family, ref, qry = pair['family'], pair['reference'], pair['query']
    rmeta, qmeta = meta[(family, ref)], meta[(family, qry)]
    path = outdir / 'blastn' / f'{family}_{ref}__{qry}.tsv'
    counts = Counter()
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 13:
                raise ValueError(f'Malformed BLASTN line in {path}: {line[:80]}')
            qid, sid, qlen, slen, qs, qe, ss, se, identity, alen, matches, evalue, score = f
            qs, qe, ss, se, alen = map(int, (qs, qe, ss, se, alen))
            identity = float(identity)
            if qid != qmeta['sequence_id'] or sid != rmeta['sequence_id']:
                raise ValueError(f'Unexpected FASTA IDs in {path}: {qid} / {sid}')
            if int(qlen) != int(qmeta['length_bp']) or int(slen) != int(rmeta['length_bp']):
                raise ValueError(f'BLASTN length mismatch in {path}')
            a, b = sorted((ss, se)); c, d = sorted((qs, qe))
            record = dict(family=family, reference=ref, query=qry,
                          reference_chr=rmeta['chromosome'], query_chr=qmeta['chromosome'],
                          start=a, end=b, start2=c, end2=d,
                          strand='+' if (qe-qs)*(se-ss) >= 0 else '-',
                          aligned_bp=alen, matching_bp=int(matches),
                          percent_identity=identity, evalue=evalue, bitscore=score,
                          reference_genomic_start=int(rmeta['start'])+a-1,
                          reference_genomic_end=int(rmeta['start'])+b-1,
                          query_genomic_start=int(qmeta['start'])+c-1,
                          query_genomic_end=int(qmeta['start'])+d-1)
            all_rows.append(record); counts['all'] += 1
            minimum = 500 if family == 'SVMP' else 150
            if alen >= minimum and identity >= 75.0:
                filtered.append(record); counts['filtered'] += 1
    print(f'{family} {ref} vs {qry}: {counts["all"]} HSPs, {counts["filtered"]} filtered')

for name, rows in (('all_local_alignment_blocks.tsv', all_rows),
                   ('filtered_blocks_for_plot.tsv', filtered)):
    with (outdir / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print(f'Output: {outdir / name} ({len(rows)} rows)')
print('NOTE: these are local HSPs, not verified orthology or nonredundant synteny blocks.')
print('===== COMPLETE =====')
PY