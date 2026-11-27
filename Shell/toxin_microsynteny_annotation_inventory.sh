# ============================================================
# inventory existing SVMP and PLA2 gene annotations
# ============================================================

# 0. activate conda env
conda activate funannotate

# 1. directories
ROOT="${ROOT:-$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo}"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"
OUTDIR="$SYNTENY/05_annotations"
MANIFEST="$SYNTENY/01_reference_inventory/anchor_bounded_regions.tsv"
ANCHORS="$SYNTENY/01_reference_inventory/selected_anchor_coordinates.tsv"
EVIDENCE="$SYNTENY/02_blast/flanking_homology_evidence"
TOXIN="$ROOT/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/toxin_annotation.gtf"
GUSS_GFF="$ROOT/annotation/funannotate/update_results/Gloydius_ussuriensis_AMNH_21010.gff3"

mkdir -p "$OUTDIR"
command -v python >/dev/null || { echo 'ERROR: python is required' >&2; exit 1; }

# 2. comparative annotation files
CGAS="$ASSEMBLIES/Cerastes_gasperettii/Cerastes_gasperettii.annotation.gff3"
CADAM="$ASSEMBLIES/Crotalus_adamanteus/Crotalus_adamanteus.annotation.gff3"
NAJA="$ASSEMBLIES/Naja_naja/Naja_naja.annotation.gff3"

# 3. verify input files
echo '===== VERIFY INPUT FILES ====='
for FILE in "$MANIFEST" "$ANCHORS" "$TOXIN" "$GUSS_GFF" "$CGAS" "$CADAM" "$NAJA" \
    "$EVIDENCE/anchor_reference_blastp_candidates.tsv" \
    "$EVIDENCE/flanking_blastp_reciprocal_gff_evidence.tsv"; do
    [[ -s "$FILE" ]] || { echo "ERROR: missing or empty input: $FILE" >&2; exit 1; }
    echo "OK: $FILE"
done

# 4. extract existing annotations
python - "$MANIFEST" "$ANCHORS" "$OUTDIR" "$TOXIN" "$EVIDENCE" \
    "$GUSS_GFF" "$CGAS" "$CADAM" "$NAJA" <<'PY'
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

manifest_path, anchors_path, outdir, toxin_path, evidence = map(Path, sys.argv[1:6])
gff_paths = dict(zip(
    ('Gloydius_ussuriensis', 'Cerastes_gasperettii',
     'Crotalus_adamanteus', 'Naja_naja'), map(Path, sys.argv[6:10])
))
comp_species = tuple(s for s in gff_paths if s != 'Gloydius_ussuriensis')
families = ('SVMP', 'PLA2')

# ------------------------------------------------------------
# 1. read genomic interval manifest and selected anchors
# ------------------------------------------------------------
def read_tsv(path, required=()):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        return list(reader)


def write_tsv(name, rows, columns):
    dest = outdir / name
    temp = outdir / f'.{name}.{os.getpid()}.tmp'
    with temp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, dest)
    print(f'OUTPUT: {dest} ({len(rows)} rows)')


regions = {}
for r in read_tsv(manifest_path, ('species', 'family', 'chromosome', 'start', 'end', 'fasta')):
    key = (r['species'], r['family'])
    if key in regions or r['species'] not in gff_paths or r['family'] not in families:
        raise ValueError(f'Duplicate/unexpected Step 5 region: {key}')
    r['start'], r['end'] = int(r['start']), int(r['end'])
    if not (1 <= r['start'] < r['end']):
        raise ValueError(f'Invalid Step 5 interval: {r}')
    if not Path(r['fasta']).is_file() or Path(r['fasta']).stat().st_size == 0:
        raise FileNotFoundError(f'Missing Step 5 FASTA: {r["fasta"]}')
    regions[key] = r
if len(regions) != 8 or any((s, f) not in regions for s in gff_paths for f in families):
    raise ValueError('Expected exactly eight Step 5 genomic intervals')

selected_anchors = defaultdict(list)
anchor_by_gene = defaultdict(list)
for a in read_tsv(anchors_path, ('species', 'family', 'anchor_symbol', 'chromosome',
                                 'start', 'end', 'gene_ids', 'evidence_status')):
    key = (a['species'], a['family'])
    if key not in regions:
        raise ValueError(f'Selected anchor has no region: {key}')
    region = regions[key]
    a['start'], a['end'] = int(a['start']), int(a['end'])
    if (a['chromosome'] != region['chromosome'] or
            not (region['start'] <= a['start'] <= a['end'] <= region['end'])):
        raise ValueError(f'Anchor outside Step 5 interval: {key} {a["anchor_symbol"]}')
    selected_anchors[key].append(a)
    for gid in filter(None, a['gene_ids'].split(';')):
        if not gid.startswith('unannotated_'):
            anchor_by_gene[(a['species'], a['family'], gid)].append(a)
if sum(map(len, selected_anchors.values())) != 32:
    raise ValueError('Expected 32 Step 5 anchor records')

# ------------------------------------------------------------
# 2. annotation parsing and conservative evidence labels
# ------------------------------------------------------------
def parse_gff(text):
    result = {}
    for part in text.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            result[key.strip()] = unquote(value.strip())
    return result


def parse_gtf(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def overlap(a, b, c, d):
    return a <= d and b >= c


def split_ids(value):
    return [v for v in value.split(',') if v]


patterns = {
    'SVMP': re.compile(r'\bSVMP\b|venom.{0,60}metalloprote|metalloprote.{0,60}venom|'
                       r'snake venom metalloproteinase', re.I),
    'PLA2': re.compile(r'\bPLA\s*[-_]?\s*2\b|\bPLA2[A-Z0-9_-]*\b|'
                       r'phospholipase\s*A\s*2|phospholipase\s*A2', re.I)
}
venom_marker = re.compile(r'\bvenom\b|\btoxin\b|\bSVMP\b', re.I)


def classify(text, family, anchor_symbols, gene_id):
    keyword = bool(patterns[family].search(text))
    if anchor_symbols:
        return ('family_keyword' if keyword else 'no_keyword',
                'selected_flanking_anchor_candidate')
    if gene_id in ('FUN_018074', 'FUN_018075') and family == 'PLA2':
        return ('family_keyword' if keyword else 'no_keyword',
                'unresolved_neighboring_gene_models')
    if keyword and venom_marker.search(text):
        return 'family_keyword', 'venom_toxin_candidate_from_annotation'
    if keyword:
        return 'family_keyword', 'family_homolog_candidate_function_unresolved'
    return 'no_keyword', 'other_annotated_gene'


def functional_hint(attributes):
    for key in ('product', 'gene', 'gene_name', 'description', 'Note', 'Name'):
        for attrs in attributes:
            if attrs.get(key):
                return attrs[key][:250]
    return 'uncharacterized'

# ------------------------------------------------------------
# 3. read comparative and Gloydius GFF3 annotations
# ------------------------------------------------------------
def build_models(species, family, region, path):
    chromosome, left, right = region['chromosome'], region['start'], region['end']
    genes, transcripts, cds_features = {}, {}, []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chromosome or f[2] not in ('gene', 'pseudogene', 'mRNA', 'transcript', 'CDS'):
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            if not overlap(start, end, left, right):
                continue
            record = dict(start=start, end=end, strand=f[6], attrs=parse_gff(f[8]))
            if f[2] in ('gene', 'pseudogene'):
                gid = record['attrs'].get('ID') or f'unnamed_gene_{chromosome}_{start}_{end}'
                genes[gid] = record
            elif f[2] in ('mRNA', 'transcript'):
                tid = record['attrs'].get('ID') or f'unnamed_transcript_{chromosome}_{start}_{end}'
                transcripts[tid] = record
            else:
                cds_features.append(record)

    grouped = defaultdict(lambda: {'gene': None, 'transcripts': {}, 'cds': []})
    for gid, gene in genes.items():
        grouped[gid]['gene'] = gene
    tx_gene = {}
    for tid, tx in transcripts.items():
        parents = split_ids(tx['attrs'].get('Parent', '') or tx['attrs'].get('gene_id', ''))
        gid = next((p for p in parents if p in genes), parents[0] if parents else f'no_gene_parent_{tid}')
        tx_gene[tid] = gid
        grouped[gid]['transcripts'][tid] = tx
    cds_by_tx = defaultdict(list)
    for feature in cds_features:
        parents = split_ids(feature['attrs'].get('Parent', ''))
        for parent in parents or [f'no_CDS_parent_{chromosome}_{feature["start"]}_{feature["end"]}']:
            if parent in tx_gene:
                gid = tx_gene[parent]
                cds_by_tx[parent].append(feature)
            elif parent in genes:
                gid = parent
            else:
                gid = f'no_gene_for_CDS_parent_{parent}'
            grouped[gid]['cds'].append(feature)

    gene_rows, model_rows = [], []
    for gid, bundle in grouped.items():
        gene = bundle['gene']
        txs = bundle['transcripts']
        cdss = bundle['cds']
        records = ([gene] if gene else []) + list(txs.values()) + cdss
        if not records:
            continue
        start = gene['start'] if gene else min(r['start'] for r in records)
        end = gene['end'] if gene else max(r['end'] for r in records)
        strand = gene['strand'] if gene else records[0]['strand']
        attrs = [r['attrs'] for r in records]
        combined = ' '.join(' '.join(map(str, d.values())) for d in attrs)
        linked_anchors = anchor_by_gene.get((species, family, gid), [])
        symbols = ';'.join(sorted(set(a['anchor_symbol'] for a in linked_anchors)))
        keyword, category = classify(combined, family, symbols, gid)
        gene_status = 'gene_model' if gene else 'gene_feature_missing_in_GFF_interval'
        boundary = 'partial_interval_overlap' if start < left or end > right else 'entire_gene_span_inside_interval'
        distinct_cds = {(c['start'], c['end'], c['strand']) for c in cdss}
        review_flag = ''
        if species == 'Gloydius_ussuriensis' and family == 'PLA2':
            if gid in {'FUN_018053', 'FUN_018054', 'FUN_018055'}:
                review_flag = 'VWA5B1_three_models_one_provisional_locus_group'
            elif gid in {'FUN_018074', 'FUN_018075'}:
                review_flag = 'neighboring_pair_distinct_models_manual_review'
        row = dict(species=species, family=family, chr=chromosome, start=start, end=end,
                   strand=strand, gene_id=gid, transcript_count=len(txs),
                   transcript_ids=';'.join(sorted(txs)), cds_segments=len(distinct_cds),
                   functional_hint=functional_hint(attrs), classification=keyword,
                   evidence_class=category, selected_anchor_symbols=symbols,
                   anchor_assignment='selected_step5_gene_id' if symbols else '',
                   gene_model_status=gene_status, interval_boundary_status=boundary,
                   review_flag=review_flag, model_source='GFF3')
        gene_rows.append(row)
        if txs:
            for tid, tx in sorted(txs.items()):
                subset = cds_by_tx.get(tid, [])
                model_rows.append(dict(**{
                    'species': species, 'family': family, 'chr': chromosome,
                    'start': tx['start'], 'end': tx['end'], 'strand': tx['strand'],
                    'gene_id': gid, 'model_id': tid, 'model_source': 'transcript',
                    'cds_segments': len({(c['start'], c['end'], c['strand']) for c in subset}),
                    'functional_hint': row['functional_hint'], 'classification': keyword,
                    'evidence_class': category, 'selected_anchor_symbols': symbols,
                    'interval_boundary_status': boundary
                }))
        else:
            model_rows.append(dict(
                species=species, family=family, chr=chromosome, start=start, end=end,
                strand=strand, gene_id=gid, model_id=gid,
                model_source='gene_only' if gene and not cdss else 'CDS_parent_only',
                cds_segments=len(distinct_cds), functional_hint=row['functional_hint'],
                classification=keyword, evidence_class=category,
                selected_anchor_symbols=symbols, interval_boundary_status=boundary
            ))
    return gene_rows, model_rows


genes_all, models_all = [], []
for species, gff in gff_paths.items():
    print(f'\n===== {species} =====')
    for family in families:
        gene_rows, model_rows = build_models(species, family, regions[(species, family)], gff)
        genes_all.extend(gene_rows)
        models_all.extend(model_rows)
        print(f'{family}: {len(gene_rows)} gene records, {len(model_rows)} transcript/model records')

# step 4 reference-anchored BLASTp candidates and reciprocal hits are supporting evidence; neither by itself establishes orthology or toxin function
anchor_hits = defaultdict(list)
for r in read_tsv(evidence / 'anchor_reference_blastp_candidates.tsv',
                  ('family', 'species', 'gene_id', 'subject_protein_id', 'bitscore')):
    if r['gene_id']:
        anchor_hits[(r['species'], r['family'], r['gene_id'])].append(r)

reciprocal_hits = defaultdict(list)
for r in read_tsv(evidence / 'flanking_blastp_reciprocal_gff_evidence.tsv',
                  ('family', 'reference_species', 'gff_gene_id', 'subject_protein_id', 'reciprocal_rank')):
    if r['gff_gene_id']:
        reciprocal_hits[(r['reference_species'], r['family'], r['gff_gene_id'])].append(r)

for row in genes_all:
    key = (row['species'], row['family'], row['gene_id'])
    ah, rh = anchor_hits[key], reciprocal_hits[key]
    row['Step4_anchor_hit_count'] = len(ah)
    row['Step4_best_anchor_bitscore'] = max((float(r['bitscore']) for r in ah if r['bitscore']), default='')
    row['Step4_anchor_query_ids'] = ';'.join(sorted(set(r.get('anchor_query', '') for r in ah if r.get('anchor_query'))))
    row['Step4_reciprocal_hit_count'] = len(rh)
    ranks = [int(r['reciprocal_rank']) for r in rh if r['reciprocal_rank'].isdigit()]
    row['Step4_best_reciprocal_rank'] = min(ranks) if ranks else ''
    row['Step4_reciprocal_query_gene_ids'] = ';'.join(sorted(set(r.get('query_gene_id', '') for r in rh if r.get('query_gene_id'))))

# ------------------------------------------------------------
# 4. recover Gloydius ToxCodAn loci and record curation status
# ------------------------------------------------------------
raw_loci = {}
with toxin_path.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        gid = parse_gtf(f[8]).get('gene_id', '')
        family = next((x for x in families if gid.startswith(x + '-')), None)
        if family is None:
            continue
        key = (family, gid, f[0])
        item = raw_loci.setdefault(key, {'segments': set(), 'strands': set()})
        item['segments'].add((int(f[3]), int(f[4])))
        item['strands'].add(f[6])

excluded = {'SVMP-6': 'redundant_overlapping_model_excluded',
            'PLA2-1': 'PLA2I_like_excluded_from_enzymatic_PLA2_count'}
toxin_loci = []
for (family, gid, chrom), item in sorted(raw_loci.items()):
    region = regions[('Gloydius_ussuriensis', family)]
    start = min(a for a, b in item['segments'])
    end = max(b for a, b in item['segments'])
    if chrom != region['chromosome'] or not overlap(start, end, region['start'], region['end']):
        continue
    fully_contained = region['start'] <= start and end <= region['end']
    if not fully_contained:
        raise ValueError(f'ToxCodAn locus not fully included in Step 5 region: {gid}')
    status = excluded.get(gid, 'core_status_not_established_from_GTF')
    toxin_loci.append(dict(
        species='Gloydius_ussuriensis', family=family, chr=chrom,
        start=start, end=end, strand=';'.join(sorted(item['strands'])),
        gene_id=gid, cds_segments=len(item['segments']),
        curation_status=status, included_in_provisional_count=(gid not in excluded),
        full_locus_in_interval='yes', source_annotation='ToxCodAn_GTF'
    ))

for row in genes_all:
    if row['species'] != 'Gloydius_ussuriensis':
        row['overlapping_ToxCodAn_IDs'] = ''
        continue
    row['overlapping_ToxCodAn_IDs'] = ';'.join(
        sorted(t['gene_id'] for t in toxin_loci
               if t['family'] == row['family'] and t['chr'] == row['chr']
               and overlap(row['start'], row['end'], t['start'], t['end']))
    )

# ------------------------------------------------------------
# 5. export annotation inventories and evidence tables
# ------------------------------------------------------------
genes_all.sort(key=lambda r: (r['species'], r['family'], r['chr'], r['start'], r['gene_id']))
models_all.sort(key=lambda r: (r['species'], r['family'], r['chr'], r['start'], r['model_id']))
toxin_loci.sort(key=lambda r: (r['family'], r['chr'], r['start'], r['gene_id']))

model_cols = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id',
              'model_id', 'model_source', 'cds_segments', 'functional_hint',
              'classification', 'evidence_class', 'selected_anchor_symbols',
              'interval_boundary_status']
gene_cols = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id',
             'transcript_count', 'transcript_ids', 'cds_segments', 'functional_hint',
             'classification', 'evidence_class', 'selected_anchor_symbols',
             'anchor_assignment', 'gene_model_status', 'interval_boundary_status',
             'review_flag', 'model_source', 'Step4_anchor_hit_count',
             'Step4_best_anchor_bitscore', 'Step4_anchor_query_ids',
             'Step4_reciprocal_hit_count', 'Step4_best_reciprocal_rank',
             'Step4_reciprocal_query_gene_ids', 'overlapping_ToxCodAn_IDs']
toxin_cols = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id',
              'cds_segments', 'curation_status', 'included_in_provisional_count',
              'full_locus_in_interval', 'source_annotation']

comparative_models = [r for r in models_all if r['species'] in comp_species]
comparative_genes = [r for r in genes_all if r['species'] in comp_species]

def get_gene_col(r, col):
    return r.get(col, '')

# populate the original model-level output with original column names preserved
write_tsv('comparative_GFF_models_in_regions.tsv', comparative_models, model_cols)
write_tsv('comparative_toxin_keyword_candidates.tsv',
          [r for r in comparative_models if r['classification'] == 'family_keyword'], model_cols)
write_tsv('Guss_ToxCodAn_loci_in_regions.tsv', toxin_loci, toxin_cols)
write_tsv('comparative_GFF_genes_in_regions.tsv', comparative_genes, gene_cols)
write_tsv('Guss_Funannotate_genes_in_regions.tsv',
          [r for r in genes_all if r['species'] == 'Gloydius_ussuriensis'], gene_cols)

# integrated table contains GFF gene records, Gloydius ToxCodAn loci, and selected genomic anchors without a corresponding annotated GFF gene
integrated = []
for g in genes_all:
    integrated.append(dict(g, record_type='GFF_gene_model',
                           source_record_id=g['gene_id'],
                           source_curation_status='',
                           provisional_locus_group=(
                               'Guss_VWA5B1_candidate' if g['review_flag'].startswith('VWA5B1_') else ''
                           )))
for t in toxin_loci:
    integrated.append(dict(t, record_type='ToxCodAn_toxin_locus',
                           source_record_id=t['gene_id'],
                           evidence_class='venom_toxin_model_provisional_curation',
                           selected_anchor_symbols='', anchor_assignment='',
                           overlapping_ToxCodAn_IDs=t['gene_id'],
                           provisional_locus_group='',
                           source_curation_status=t['curation_status']))
for key, anchors in selected_anchors.items():
    for a in anchors:
        matches = [g for g in genes_all if g['species'] == key[0]
                   and g['family'] == key[1]
                   and a['anchor_symbol'] in g['selected_anchor_symbols'].split(';')]
        if matches:
            continue
        integrated.append(dict(
            species=a['species'], family=a['family'], chr=a['chromosome'],
            start=a['start'], end=a['end'], strand=a.get('strand', ''),
            gene_id='', record_type='unannotated_selected_anchor_locus',
            source_record_id=f'{a["species"]}_{a["family"]}_{a["anchor_symbol"]}',
            evidence_class='provisional_genomic_anchor_candidate',
            selected_anchor_symbols=a['anchor_symbol'],
            anchor_assignment='genomic_coordinate_only',
            source_curation_status='not_confirmed_ortholog',
            provisional_locus_group='', review_flag='no_selected_GFF_gene_model'
        ))

integrated.sort(key=lambda r: (r['species'], r['family'], r['chr'], int(r['start']),
                               r['record_type'], r['source_record_id']))
integrated_cols = ['species', 'family', 'chr', 'start', 'end', 'strand',
                   'record_type', 'source_record_id', 'gene_id', 'transcript_count',
                   'transcript_ids', 'cds_segments', 'functional_hint',
                   'classification', 'evidence_class', 'selected_anchor_symbols',
                   'anchor_assignment', 'gene_model_status',
                   'interval_boundary_status', 'provisional_locus_group',
                   'review_flag', 'Step4_anchor_hit_count', 'Step4_best_anchor_bitscore',
                   'Step4_anchor_query_ids', 'Step4_reciprocal_hit_count',
                   'Step4_best_reciprocal_rank', 'Step4_reciprocal_query_gene_ids',
                   'overlapping_ToxCodAn_IDs', 'source_curation_status']
write_tsv('comparative_gene_evidence_inventory.tsv', integrated, integrated_cols)

summary = []
for species in gff_paths:
    for family in families:
        gs = [g for g in genes_all if g['species'] == species and g['family'] == family]
        ts = [t for t in toxin_loci if species == 'Gloydius_ussuriensis' and t['family'] == family]
        summary.append(dict(
            species=species, family=family,
            total_GFF_gene_records=len(gs),
            annotation_keyword_gene_candidates=sum(g['classification'] == 'family_keyword' for g in gs),
            annotation_venom_candidate_genes=sum(g['evidence_class'] == 'venom_toxin_candidate_from_annotation' for g in gs),
            selected_anchor_gene_models=sum(bool(g['selected_anchor_symbols']) for g in gs),
            raw_ToxCodAn_locus_IDs=len(ts) if species == 'Gloydius_ussuriensis' else '',
            provisionally_included_ToxCodAn_IDs=sum(t['included_in_provisional_count'] for t in ts) if species == 'Gloydius_ussuriensis' else '',
            interpretation='requires_manual_toxin_family_confirmation_and_copy_number_review'
        ))
write_tsv('comparative_toxin_candidate_summary.tsv', summary,
          ['species', 'family', 'total_GFF_gene_records',
           'annotation_keyword_gene_candidates', 'annotation_venom_candidate_genes',
           'selected_anchor_gene_models', 'raw_ToxCodAn_locus_IDs',
           'provisionally_included_ToxCodAn_IDs', 'interpretation'])

flags = [
    dict(species='Gloydius_ussuriensis', family='PLA2',
         locus='FUN_018053;FUN_018054;FUN_018055',
         issue='three_separate_GFF_models_with_complementary_VWA5B1_protein_alignments',
         action='retain_original_IDs;evaluate_as_one_provisional_locus;review_splice_support'),
    dict(species='Gloydius_ussuriensis', family='PLA2',
         locus='FUN_018074;FUN_018075',
         issue='neighboring_gene_models_with_different_cross_species_best_matches',
         action='retain_as_distinct_provisional_models;review_exons_and_reference_annotations'),
    dict(species='Gloydius_ussuriensis', family='SVMP',
         locus='G_ussuri_chr9:5383702-5413622',
         issue='ADAM28_like_TBLASTN_candidate_not_validated_ortholog',
         action='review_domains_exon_structure_and_local_synteny'),
    dict(species='Cerastes_gasperettii', family='PLA2',
         locus='venom_PLA2_cluster',
         issue='figure_three_candidates_versus_two_published_venom_PLA2_genes',
         action='compare_gene_IDs_coordinates_CDS_models_and_published_toxin_catalog'),
    dict(species='Gloydius_ussuriensis', family='SVMP;PLA2',
         locus='ToxCodAn_core_set',
         issue='GTF_only_cannot_establish_full_manually_curated_core_set',
         action='reconcile_provisional_IDs_with_curated_locus_list_before_reporting_copies')
]
write_tsv('Step6_annotation_review_flags.tsv', flags,
          ['species', 'family', 'locus', 'issue', 'action'])

# ------------------------------------------------------------
# 6. print summary
# ------------------------------------------------------------
print('\n===== COMPARATIVE TOXIN ANNOTATION SUMMARY =====')
for r in summary:
    print(f"{r['species']:<23} {r['family']:<5} "
          f"genes={r['total_GFF_gene_records']:<4} "
          f"keyword_genes={r['annotation_keyword_gene_candidates']:<3} "
          f"venom_annotation_candidates={r['annotation_venom_candidate_genes']}")
print('\n===== GLOYDIUS TOXCODAN LOCUS INVENTORY =====')
for family in families:
    ts = [t for t in toxin_loci if t['family'] == family]
    n = sum(t['included_in_provisional_count'] for t in ts)
    print(f'{family}: {len(ts)} raw in-region IDs; {n} provisionally included; '
          f'{len(ts)-n} explicitly excluded')
print('\nNOTE: keyword matches, protein homology and provisional loci are NOT validated toxin copies or ortholog calls.')
print('NOTE: comparative toxin-copy completeness and the Cerastes PLA2 discrepancy require manual review.')
print('===== COMPLETE =====')