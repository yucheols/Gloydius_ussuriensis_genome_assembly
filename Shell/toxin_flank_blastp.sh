#!/bin/bash
#SBATCH --job-name=toxin_flank_blastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# SVMP and PLA2 flanking-gene protein homology
# ============================================================

# activate conda env
source ~/.bash_profile
conda activate funannotate

set -euo pipefail
RESUME_BLAST="${RESUME_BLAST:-1}"

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"
ANNOTATION="$ROOT/annotation/funannotate/update_results"
PROTEINS="$ANNOTATION/Gloydius_ussuriensis_AMNH_21010.proteins.fa"
FLANKING="$SYNTENY/05_annotations/Guss_SVMP_PLA2_flanking_genes.tsv"
BLASTDIR="$SYNTENY/02_blast"
GUSS_GENOME="$ROOT/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
GUSS_GFF="$ANNOTATION/Gloydius_ussuriensis_AMNH_21010.gff3"
EVIDENCE="$BLASTDIR/flanking_homology_evidence"

mkdir -p "$BLASTDIR/SVMP" "$BLASTDIR/PLA2" "$EVIDENCE"

# 2. verify inputs and software
printf '\n===== VERIFY INPUTS =====\n'
for FILE in "$PROTEINS" "$FLANKING" "$GUSS_GENOME" "$GUSS_GFF"; do
    if [[ ! -s "$FILE" ]]; then echo "ERROR: missing input: $FILE" >&2; exit 1; fi
    echo "OK: $FILE"
done

for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
    for FILE in "$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa" "$ASSEMBLIES/$SPECIES/${SPECIES}.annotation.gff3"; do
        if [[ ! -s "$FILE" ]]; then echo "ERROR: missing input: $FILE" >&2; exit 1; fi
    done
done

for TOOL in python blastp makeblastdb tblastn samtools sha256sum; do command -v "$TOOL" >/dev/null || exit 1; done
if [[ ! -s "${GUSS_GENOME}.fai" ]]; then samtools faidx "$GUSS_GENOME"; fi

# Reuse alignment files only when query, database source and parameters match.
# Legacy outputs without a matching signature are refreshed once.
blast_signature() {
    local query="$1" source="$2" parameters="$3"
    { sha256sum "$query" "$source"; printf '%s\n' "$parameters"; } | sha256sum | awk '{print $1}'
}

run_cached_blastp() {
    local query="$1" db="$2" source="$3" output="$4" targets="$5" format="$6"
    local signature prior tmp
    signature=$(blast_signature "$query" "$source" "blastp|evalue=1e-5|max_target_seqs=$targets|outfmt=$format")
    prior="$(cat "${output}.signature" 2>/dev/null || true)"
    if [[ "$RESUME_BLAST" == 1 && -f "$output" && "$prior" == "$signature" ]]; then
        echo "Reusing verified BLASTp: $output"
    else
        tmp="${output}.tmp.${SLURM_JOB_ID:-$$}"
        blastp -query "$query" -db "$db" -evalue 1e-5 -max_target_seqs "$targets" -num_threads 2 -outfmt "$format" -out "$tmp"
        mv -f "$tmp" "$output"
        printf '%s\n' "$signature" > "${output}.signature"
    fi
}

run_cached_tblastn() {
    local query="$1" subject="$2" output="$3" format="$4"
    local signature prior tmp
    signature=$(blast_signature "$query" "$subject" "tblastn|evalue=1e-5|seg=yes|outfmt=$format")
    prior="$(cat "${output}.signature" 2>/dev/null || true)"
    if [[ "$RESUME_BLAST" == 1 && -f "$output" && "$prior" == "$signature" ]]; then
        echo "Reusing verified TBLASTN: $output"
    else
        tmp="${output}.tmp.${SLURM_JOB_ID:-$$}"
        tblastn -query "$query" -subject "$subject" -seg yes -evalue 1e-5 -num_threads 2 -outfmt "$format" -out "$tmp"
        mv -f "$tmp" "$output"
        printf '%s\n' "$signature" > "${output}.signature"
    fi
}


# 3. extract Gloydius flanking proteins 
printf '\n===== EXTRACT FLANKING PROTEINS =====\n'
python - "$PROTEINS" "$FLANKING" "$BLASTDIR" <<'PY'
import csv, re, sys
from pathlib import Path
from collections import defaultdict
proteins_file, flanking_file, outdir = map(Path, sys.argv[1:])
proteins = defaultdict(list)
header, seq = None, []
def save():
    if header is None: return
    parts = header.split()
    transcript = parts[0]
    gene = parts[1] if len(parts) > 1 and re.fullmatch(r'FUN_\d+', parts[1]) else re.sub(r'-T\d+$', '', transcript)
    proteins[gene].append((transcript, ''.join(seq)))
with proteins_file.open() as handle:
    for line in handle:
        line = line.strip()
        if line.startswith('>'):
            save()
            header, seq = line[1:], []
        elif line: seq.append(line)
    save()
with flanking_file.open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
for family in ('SVMP', 'PLA2'):
    selected = []
    for side in ('left_flank', 'right_flank'):
        subset = sorted((r for r in rows if r['family'] == family and r['position'] == side), key=lambda r: int(r['distance_bp']))
        selected.extend(subset[:8])
    manifest = []
    fasta_out = outdir / family / f'Guss_{family}_flanking_proteins.faa'
    manifest_out = outdir / family / f'Guss_{family}_query_manifest.tsv'
    with fasta_out.open('w') as fasta:
        for row in selected:
            gene = row['gene_id']
            candidates = proteins.get(gene, [])
            if not candidates:
                print(f'WARNING: No protein found for {gene}')
                continue
            transcript, sequence = next((c for c in candidates if c[0] == f'{gene}-T1'), max(candidates, key=lambda x: len(x[1])))
            if not sequence: continue
            fasta.write(f'>{transcript} {gene}\n')
            for i in range(0, len(sequence), 80): fasta.write(sequence[i:i+80] + '\n')
            manifest.append(dict(gene_id=gene, transcript_id=transcript, position=row['position'], distance_bp=row['distance_bp'], length_aa=len(sequence)))
    if not manifest: raise RuntimeError(f'No flanking proteins extracted for {family}')
    with manifest_out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['gene_id', 'transcript_id', 'position', 'distance_bp', 'length_aa'], delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest)
    print(f'{family}: {len(manifest)} proteins; {fasta_out}')
# Additional PLA2 models for pairwise fragmentation/paralogy diagnostics.
discrepancy = ['FUN_018053', 'FUN_018054', 'FUN_018055', 'FUN_018074', 'FUN_018075', 'FUN_018058']
path = outdir / 'flanking_homology_evidence' / 'PLA2_discrepancy_queries.faa'
manifest = []
with path.open('w') as fasta:
    for gene in discrepancy:
        choices = proteins.get(gene, [])
        if not choices:
            print(f'WARNING: PLA2 discrepancy protein unavailable: {gene}')
            continue
        transcript, sequence = next((x for x in choices if x[0] == f'{gene}-T1'), max(choices, key=lambda x: len(x[1])))
        if not sequence: continue
        fasta.write(f'>{transcript} {gene}\n')
        for i in range(0, len(sequence), 80): fasta.write(sequence[i:i+80] + '\n')
        manifest.append(dict(gene_id=gene, transcript_id=transcript, length_aa=len(sequence)))
with (path.parent / 'PLA2_discrepancy_query_manifest.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=['gene_id', 'transcript_id', 'length_aa'], delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(manifest)
if len(manifest) < 6:
    print(f'WARNING: {6-len(manifest)} of 6 PLA2 discrepancy models lack translated proteins')
print(f'PLA2 model-reconciliation queries: {len(manifest)}; {path}')
PY

# 4. prepare comparative BLAST databases and reciprocal database
printf '\n===== PREPARE PROTEIN BLAST DATABASES =====\n'
for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
    PROTEOME="$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa"
    DB="$BLASTDIR/${SPECIES}_proteins_db"
    if [[ ( ! -s "${DB}.pin" && ! -s "${DB}.pdb" ) || ( -f "${DB}.pin" && "$PROTEOME" -nt "${DB}.pin" ) || ( -f "${DB}.pdb" && "$PROTEOME" -nt "${DB}.pdb" ) ]]; then makeblastdb -in "$PROTEOME" -dbtype prot -out "$DB"; fi
done
GUSS_DB="$EVIDENCE/Gloydius_proteins_db"
if [[ ( ! -s "${GUSS_DB}.pin" && ! -s "${GUSS_DB}.pdb" ) || ( -f "${GUSS_DB}.pin" && "$PROTEINS" -nt "${GUSS_DB}.pin" ) || ( -f "${GUSS_DB}.pdb" && "$PROTEINS" -nt "${GUSS_DB}.pdb" ) ]]; then makeblastdb -in "$PROTEINS" -dbtype prot -out "$GUSS_DB"; fi

# 5. forward BLASTp 
printf '\n===== FORWARD BLASTP =====\n'
for FAMILY in SVMP PLA2; do
    QUERY="$BLASTDIR/$FAMILY/Guss_${FAMILY}_flanking_proteins.faa"
    for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        OUTPUT="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.blastp.tsv"
        run_cached_blastp "$QUERY" "$BLASTDIR/${SPECIES}_proteins_db" "$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa" "$OUTPUT" 20 '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore'
        echo "$FAMILY vs $SPECIES: $(wc -l < "$OUTPUT") HSPs"
    done
done

# 6. extract top forward subject sequences for reciprocal BLASTp
printf '\n===== PREPARE RECIPROCAL QUERIES =====\n'
python - "$BLASTDIR" "$ASSEMBLIES" "$EVIDENCE" <<'PY'
import csv, sys
from pathlib import Path
from collections import defaultdict
blastdir, assemblies, evidence = map(Path, sys.argv[1:])
def fasta_records(path):
    header, lines = None, []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None: yield header.split()[0], ''.join(lines)
                header, lines = line[1:], []
            elif line.strip(): lines.append(line.strip())
        if header is not None: yield header.split()[0], ''.join(lines)
for family in ('SVMP', 'PLA2'):
    for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        path = blastdir / family / f'Guss_{family}_vs_{species}.blastp.tsv'
        hits = defaultdict(list)
        with path.open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) >= 9: hits[row[0]].append((row[1], float(row[8])))
        subjects = set()
        for query, matches in hits.items():
            for subject, score in sorted(matches, key=lambda r: -r[1])[:5]: subjects.add(subject)
        ref = assemblies / species / f'{species}.protein.faa'
        records = {header: seq for header, seq in fasta_records(ref) if header in subjects}
        missing = subjects - records.keys()
        if missing: print(f'WARNING: {species} {family}: {len(missing)} BLAST subjects not recovered from FASTA')
        out = evidence / f'{family}_{species}_reciprocal_queries.faa'
        with out.open('w') as handle:
            for header, seq in sorted(records.items()):
                handle.write(f'>{header}\n')
                for i in range(0, len(seq), 80): handle.write(seq[i:i+80] + '\n')
        print(f'{family} {species}: {len(records)} reciprocal protein queries')
PY

# 7. reciprocal BLASTp against the Gloydius proteome
printf '\n===== RECIPROCAL BLASTP =====\n'
for FAMILY in SVMP PLA2; do
    for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        QUERY="$EVIDENCE/${FAMILY}_${SPECIES}_reciprocal_queries.faa"
        OUTPUT="$EVIDENCE/${FAMILY}_${SPECIES}_reciprocal.blastp.tsv"
        if [[ -s "$QUERY" ]]; then
            run_cached_blastp "$QUERY" "$GUSS_DB" "$PROTEINS" "$OUTPUT" 5 '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore'
        else
            : > "$OUTPUT"
        fi
    done
done

# 8. build protein reference panel from annotated Cerastes and Crotalus proteins
# make ADAM selection deliberately conservative == ambiguous source IDs are reported
printf '\n===== EXTRACT TARGETED ANCHOR REFERENCE PROTEINS =====\n'
python - "$ASSEMBLIES" "$EVIDENCE" <<'PY'
import csv, re, sys
from pathlib import Path
from urllib.parse import unquote
assemblies, evidence = map(Path, sys.argv[1:])
def fasta_records(path):
    header, lines = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if header is not None: yield header, ''.join(lines)
                header, lines = line[1:], []
            elif line: lines.append(line)
        if header is not None: yield header, ''.join(lines)
cer = list(fasta_records(assemblies / 'Cerastes_gasperettii' / 'Cerastes_gasperettii.protein.faa'))
cer_by_id = {h.split()[0]: (h, s) for h, s in cer}
patterns = {'STC1': r'(?i)STC[-_ ]?1|STANNIOCALCIN[-_ ]?1', 'NEFM': r'(?i)(?<![A-Z0-9])NEFM(?![A-Z0-9])|NEUROFILAMENT[^;]*MEDIUM', 'NEFL': r'(?i)(?<![A-Z0-9])NEFL(?![A-Z0-9])|NEUROFILAMENT[^;]*LIGHT', 'ADAM': r'(?i)(?<![A-Z0-9])(?:RNA-|GENE-)?ADAM(?![0-9A-Z])', 'OTUD3': r'(?i)OTUD[-_ ]?3', 'MUL1': r'(?i)(?<![A-Z0-9])MUL[-_ ]?1(?![A-Z0-9])', 'VWA5B1': r'(?i)VWA5B1', 'PLA2G2F': r'(?i)PLA2G2F|PLA2[-_ ]?G2F'}
exact_ids = {'STC1': 'crotalus-ruber-cadam-rna-Stc1-gene14556-R0-R0', 'NEFM': 'rna-NEFM', 'NEFL': 'rna-NEFL', 'OTUD3': 'crotalus-ruber-Otud3-R0', 'MUL1': 'crotalus-ruber-Mul1-R0', 'VWA5B1': 'crotalus-ruber-cadam-rna-Vwa5b1-gene17777-R0-R0', 'PLA2G2F': 'rna-Pla2g2F'}
cer_gff = assemblies / 'Cerastes_gasperettii' / 'Cerastes_gasperettii.annotation.gff3'
cer_gene_coords = {}
with cer_gff.open() as handle:
    for line in handle:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene': continue
        a = {k: unquote(v) for part in f[8].split(';') if '=' in part for k, v in [part.split('=', 1)]}
        if a.get('ID'): cer_gene_coords[a['ID']] = (f[0], f[3], f[4], f[6])
def header_value(header, field):
    m = re.search(r'\b' + re.escape(field) + r'=([^\s]+)', header)
    return m.group(1) if m else ''
def audit_row(symbol, species, status, header='', sequence='', alternatives=0):
    gene_id, seq_id = header_value(header, 'gene'), header_value(header, 'seq_id')
    chrom, start, end, strand = cer_gene_coords.get(gene_id, ('', '', '', '')) if species == 'Cerastes_gasperettii' else ('', '', '', '')
    return dict(anchor_symbol=symbol, species=species, selection_status=status, source_header=header, length_aa=len(sequence) if sequence else '', source_protein_id=header.split()[0] if header else '', source_gene_id=gene_id, source_seq_id=seq_id, comparative_GFF_chr=chrom, comparative_GFF_start=start, comparative_GFF_end=end, comparative_GFF_strand=strand, alternate_header_matches=alternatives, genomic_context_status='gff_gene_mapped' if chrom else 'GFF_link_unresolved_review')
found, audit = [], []
for symbol, pattern in patterns.items():
    matches = [(h, s) for h, s in cer if re.search(pattern, h)]
    desired = exact_ids.get(symbol)
    selected = cer_by_id.get(desired) if desired else None
    if selected:
        h, seq = selected
        found.append(('Cerastes', symbol, h, seq))
        status = 'exact_protein_ID_selected_multiple_header_candidates' if len(matches) > 1 else 'exact_protein_ID_selected'
        audit.append(audit_row(symbol, 'Cerastes_gasperettii', status, h, seq, max(0, len(matches) - 1)))
        if symbol == 'MUL1' and header_value(h, 'seq_id') not in ('CHR17', ''):
            print(f'WARNING: selected Cerastes MUL1 sequence not on expected CHR17: {h}')
        print(f'Cerastes {symbol}: selected exact protein {desired}; {len(matches)} keyword matches')
    else:
        audit.append(audit_row(symbol, 'Cerastes_gasperettii', f'exact_ID_missing_{len(matches)}_alternatives_review', ';'.join(h[:100] for h, _ in matches[:5]), '', len(matches)))
        print(f'WARNING: Cerastes {symbol}: exact protein ID {desired or "not_specified"} absent; {len(matches)} candidates; anchor not auto-selected')
# Resolve Crotalus ADAM28 using the gene ID in its GFF and protein_id on CDS.
sp = 'Crotalus_adamanteus'
crot = list(fasta_records(assemblies / sp / f'{sp}.protein.faa'))
crot_by_id = {h.split()[0]: (h, s) for h, s in crot}
gff = assemblies / sp / f'{sp}.annotation.gff3'
features, parents, aliases = {}, {}, set()
with gff.open() as handle:
    for line in handle:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = {k: unquote(v) for part in f[8].split(';') if '=' in part for k, v in [part.split('=', 1)]}
        fid = a.get('ID', '')
        if fid:
            features[fid] = a
            parents[fid] = a.get('Parent', '').split(',')
        if fid == 'gene-NXF25_014549': aliases.add(fid)
def belongs_to_adam(fid):
    seen = set()
    todo = [fid]
    while todo:
        node = todo.pop()
        if node in seen: continue
        if node == 'gene-NXF25_014549': return True
        seen.add(node)
        todo.extend(parents.get(node, []))
    return False
ids = set()
for fid, a in features.items():
    if belongs_to_adam(fid):
        for key in ('protein_id', 'Name', 'ID'):
            if a.get(key): ids.update(a[key].split(','))
matches = [crot_by_id[k] for k in sorted(ids) if k in crot_by_id]
if len(matches) != 1:
    matches = [(h, s) for h, s in crot if re.search(r'(?i)(?<![A-Z0-9])ADAM[-_ ]?28(?![0-9A-Z])|NXF25_014549', h)]
if len(matches) == 1:
    h, s = matches[0]
    found.append(('Crotalus', 'ADAM28', h, s))
    audit.append(audit_row('ADAM28', sp, 'unique_protein_match', h, s))
else:
    audit.append(audit_row('ADAM28', sp, f'{len(matches)}_protein_matches_review', ';'.join(h[:100] for h, _ in matches[:5]), '', len(matches)))
    print('WARNING: Crotalus ADAM28 protein not uniquely recoverable; Cerastes ADAM may still provide the search')

# Keep protein panels family-specific to avoid irrelevant cross-family matches.
for family in ('SVMP', 'PLA2'):
    members = [('Cerastes', 'STC1'), ('Cerastes', 'NEFM'), ('Cerastes', 'NEFL'), ('Cerastes', 'ADAM'), ('Crotalus', 'ADAM28')] if family == 'SVMP' else [('Cerastes', 'OTUD3'), ('Cerastes', 'MUL1'), ('Cerastes', 'VWA5B1'), ('Cerastes', 'PLA2G2F')]
    selected = [(sp, symbol, h, seq) for sp, symbol, h, seq in found if (sp, symbol) in members]
    with (evidence / f'{family}_anchor_reference_proteins.faa').open('w') as handle:
        for species, symbol, header, seq in selected:
            handle.write(f'>{species}_{symbol}\n')
            for i in range(0, len(seq), 80): handle.write(seq[i:i+80] + '\n')
    if not selected: print(f'WARNING: no reference proteins for {family}')
audit_cols = ['anchor_symbol', 'species', 'selection_status', 'source_header', 'length_aa', 'source_protein_id', 'source_gene_id', 'source_seq_id', 'comparative_GFF_chr', 'comparative_GFF_start', 'comparative_GFF_end', 'comparative_GFF_strand', 'alternate_header_matches', 'genomic_context_status']
with (evidence / 'anchor_reference_protein_audit.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=audit_cols, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(audit)
with (evidence / 'PLA2_anchor_reference_audit.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=audit_cols, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows([r for r in audit if r['anchor_symbol'] in ('OTUD3', 'MUL1', 'VWA5B1', 'PLA2G2F')])
if not any(species == 'Cerastes' and symbol == 'ADAM' for species, symbol, _, _ in found) and not any(symbol == 'ADAM28' for _, symbol, _, _ in found):
    raise RuntimeError('No uniquely recoverable ADAM reference protein; see anchor_reference_protein_audit.tsv and resolve the reference ID')
print('Targeted proteins:', ', '.join(f'{s}_{g}' for s, g, _, _ in found))
PY

# 8b. all-species reference-anchor BLASTp, including Gloydius.
# This detects candidate anchors even when Gloydius query proteins were missing.
printf '\n===== REFERENCE ANCHOR BLASTP IN ALL FOUR SPECIES =====\n'
for FAMILY in SVMP PLA2; do
    QUERY="$EVIDENCE/${FAMILY}_anchor_reference_proteins.faa"
    for SPECIES in Gloydius_ussuriensis Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        OUTPUT="$EVIDENCE/${FAMILY}_reference_anchors_vs_${SPECIES}.blastp.tsv"
        if [[ ! -s "$QUERY" ]]; then
            : > "$OUTPUT"
            continue
        fi
        if [[ "$SPECIES" == 'Gloydius_ussuriensis' ]]; then DB="$GUSS_DB"; else DB="$BLASTDIR/${SPECIES}_proteins_db"; fi
        if [[ "$SPECIES" == 'Gloydius_ussuriensis' ]]; then SOURCE="$PROTEINS"; else SOURCE="$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa"; fi
        run_cached_blastp "$QUERY" "$DB" "$SOURCE" "$OUTPUT" 10 '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore'
        echo "$FAMILY $SPECIES: $(wc -l < "$OUTPUT") anchor BLASTp HSPs"
    done
done

# 9. derive Gloydius search intervals from flanking coordinates obtained in step 2
printf '\n===== EXTRACT GLOYDIUS GENOMIC SEARCH INTERVALS =====\n'
python - "$FLANKING" "${GUSS_GENOME}.fai" "$EVIDENCE/Guss_target_intervals.tsv" <<'PY'
import csv, sys
from collections import defaultdict
flanking, fai, output = sys.argv[1:]
lengths = {}
with open(fai) as handle:
    for line in handle:
        parts = line.split('\t')
        lengths[parts[0]] = int(parts[1])
with open(flanking) as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
regions = []
for fam in ('SVMP', 'PLA2'):
    subset = [r for r in rows if r['family'] == fam and r['chr'] in lengths]
    if not subset: raise RuntimeError(f'Missing genomic coordinates in flanking inventory for {fam}')
    chrs = set(r['chr'] for r in subset)
    if len(chrs) != 1: raise RuntimeError(f'{fam} spans multiple chromosomes in flanking inventory: {chrs}')
    chrom = next(iter(chrs))
    start = max(1, min(int(r['start']) for r in subset) - 20000)
    end = min(lengths[chrom], max(int(r['end']) for r in subset) + 20000)
    regions.append(dict(family=fam, chromosome=chrom, start=start, end=end, search_region=f'{chrom}:{start}-{end}'))
with open(output, 'w', newline='') as handle:
    w = csv.DictWriter(handle, fieldnames=['family', 'chromosome', 'start', 'end', 'search_region'], delimiter='\t', lineterminator='\n')
    w.writeheader()
    w.writerows(regions)
for r in regions: print(r['family'], r['search_region'])
PY

while IFS=$'\t' read -r FAMILY CHR START END REGION; do
    if [[ "$FAMILY" == 'family' ]]; then continue; fi
    REGION="${REGION%$'\r'}"
    GENOME_REGION="$EVIDENCE/Guss_${FAMILY}_search_region.fa"
    samtools faidx "$GUSS_GENOME" "$REGION" > "$GENOME_REGION"
    OUTPUT="$EVIDENCE/Guss_${FAMILY}_anchors.tblastn.tsv"
    if [[ ! -s "$EVIDENCE/${FAMILY}_anchor_reference_proteins.faa" ]]; then : > "$OUTPUT"; continue; fi
    run_cached_tblastn "$EVIDENCE/${FAMILY}_anchor_reference_proteins.faa" "$GENOME_REGION" "$OUTPUT" '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore'
    echo "$FAMILY: $(wc -l < "$OUTPUT") TBLASTN HSPs"
done < "$EVIDENCE/Guss_target_intervals.tsv"

printf '\n===== PLA2 DETAILED MODEL-RECONCILIATION BLASTP =====\n'
QUERY="$EVIDENCE/PLA2_discrepancy_queries.faa"
if [[ -s "$QUERY" ]]; then
    for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        OUTPUT="$EVIDENCE/PLA2_discrepancy_vs_${SPECIES}.blastp.tsv"
        run_cached_blastp "$QUERY" "$BLASTDIR/${SPECIES}_proteins_db" "$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa" "$OUTPUT" 30 '6 qseqid sseqid pident length qlen slen qstart qend sstart send qcovs evalue bitscore'
        echo "PLA2 models vs $SPECIES: $(wc -l < "$OUTPUT") alignment HSPs"
    done
else
    echo "WARNING: no PLA2 discrepancy proteins available"
    for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do : > "$EVIDENCE/PLA2_discrepancy_vs_${SPECIES}.blastp.tsv"; done
fi

# 10. consolidate BLASTp reciprocal evidence, GFF coordinates and TBLASTN hits
printf '\n===== GENERATE REVIEW TABLES =====\n'
python - "$BLASTDIR" "$ASSEMBLIES" "$FLANKING" "$EVIDENCE" "$GUSS_GFF" <<'PY'
import csv, re, sys
from pathlib import Path
from collections import defaultdict
from urllib.parse import unquote
blastdir, assemblies, flanking, evidence, guss_gff = [Path(x) for x in sys.argv[1:]]
def read_tsv(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter='\t'))
def write_tsv(path, rows, cols):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

# Subject protein accession -> annotated GFF locus; records may only map to a CDS.
def build_gff_index(path):
    features, parent, alias = {}, {}, defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9: continue
            a = {k: unquote(v) for part in f[8].split(';') if '=' in part for k, v in [part.split('=', 1)]}
            fid = a.get('ID', '') or f'_record_{len(features)}'
            features[fid] = (f, a)
            parent[fid] = a.get('Parent', '').split(',') if a.get('Parent') else []
            for k in ('ID', 'Name', 'protein_id', 'transcript_id', 'gene_id'):
                for v in a.get(k, '').split(','):
                    if v: alias[v].append(fid)
    return features, parent, alias

def map_gene(subject, index):
    features, parent, alias = index
    matches = alias.get(subject, [])
    if not matches: return ('', '', '', '', '', 'unmapped_GFF_subject')
    resolved = []
    for node in matches:
        seen = set()
        todo = [node]
        while todo:
            cur = todo.pop()
            if cur in seen or cur not in features: continue
            seen.add(cur)
            f, a = features[cur]
            if f[2] == 'gene': resolved.append((cur, f))
            todo.extend(parent.get(cur, []))
    if resolved:
        gid, f = sorted(set((g, tuple(f)) for g, f in resolved), key=lambda x: x[0])[0]
        return (gid, f[0], f[3], f[4], f[6], 'gene_model_mapped')
    f, a = features[matches[0]]
    return (matches[0], f[0], f[3], f[4], f[6], 'feature_only_review')

index = {species: build_gff_index(assemblies / species / f'{species}.annotation.gff3') for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja')}
guss_index = build_gff_index(guss_gff)
guss_genes = []
for fid, (f, a) in guss_index[0].items():
    if f[2] == 'gene':
        guss_genes.append(dict(gene_id=a.get('ID', fid), chr=f[0], start=int(f[3]), end=int(f[4]), strand=f[6], name=a.get('Name', ''), product=a.get('product', '')))

def overlapping_genes(chrom, start, end):
    return sorted((g for g in guss_genes if g['chr'] == chrom and g['start'] <= end and g['end'] >= start), key=lambda g: (g['start'], g['gene_id']))

def nearest_gene(chrom, start, end):
    candidates = (g for g in guss_genes if g['chr'] == chrom)
    return min(((max(0, start - g['end'], g['start'] - end), g) for g in candidates), key=lambda item: item[0], default=(None, None))

flanks = {r['gene_id']: r for r in read_tsv(flanking)}
known = {'FUN_014513':'STC1', 'FUN_014526':'NEFM', 'FUN_014527':'NEFL', 'FUN_018057':'MUL1', 'FUN_018058':'PLA2G2F', 'FUN_018073':'OTUD3'}

forward = []
for family in ('SVMP', 'PLA2'):
    for species in index:
        reciprocal = defaultdict(list)
        rpath = evidence / f'{family}_{species}_reciprocal.blastp.tsv'
        with rpath.open() as handle:
            for r in csv.reader(handle, delimiter='\t'):
                if len(r) >= 9: reciprocal[r[0]].append(r)
        fpath = blastdir / family / f'Guss_{family}_vs_{species}.blastp.tsv'
        with fpath.open() as handle:
            rows = [r for r in csv.reader(handle, delimiter='\t') if len(r) >= 9]
        rows.sort(key=lambda r: (r[0], -float(r[8])))
        ranks = defaultdict(int)
        for r in rows:
            ranks[r[0]] += 1
            if ranks[r[0]] > 20: continue
            query, subject = r[:2]
            source_gene = re.sub(r'-T\d+$', '', query)
            rev_hits = sorted(reciprocal.get(subject, []), key=lambda x: -float(x[8]))
            reciprocal_best = rev_hits[0][1] if rev_hits else ''
            reciprocal_rank = next((i+1 for i, x in enumerate(rev_hits) if re.sub(r'-T\d+$', '', x[1]) == source_gene), '')
            geneid, chrom, start, end, strand, mapping = map_gene(subject, index[species])
            forward.append(dict(family=family, query_gene_id=source_gene, anchor_hint=known.get(source_gene, ''), query_transcript=query, reference_species=species, subject_protein_id=subject, forward_rank=ranks[query], percent_identity=r[2], aligned_aa=r[3], query_length_aa=r[4], subject_length_aa=r[5], query_coverage_pct=r[6], evalue=r[7], bitscore=r[8], reciprocal_best=reciprocal_best, reciprocal_rank=reciprocal_rank, gff_gene_id=geneid, chr=chrom, start=start, end=end, strand=strand, gff_mapping_status=mapping, validation_status='candidate_requires_review'))

cols = ['family', 'query_gene_id', 'anchor_hint', 'query_transcript', 'reference_species', 'subject_protein_id', 'forward_rank', 'percent_identity', 'aligned_aa', 'query_length_aa', 'subject_length_aa', 'query_coverage_pct', 'evalue', 'bitscore', 'reciprocal_best', 'reciprocal_rank', 'gff_gene_id', 'chr', 'start', 'end', 'strand', 'gff_mapping_status', 'validation_status']
write_tsv(evidence / 'flanking_blastp_reciprocal_gff_evidence.tsv', forward, cols)

# Reference-anchored BLASTp across all four species. These are candidates, not called orthologs.
reference_candidates = []
for family in ('SVMP', 'PLA2'):
    for species in ('Gloydius_ussuriensis', 'Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        path = evidence / f'{family}_reference_anchors_vs_{species}.blastp.tsv'
        hits = defaultdict(list)
        with path.open() as handle:
            for r in csv.reader(handle, delimiter='\t'):
                if len(r) >= 9: hits[r[0]].append(r)
        for query, matches in sorted(hits.items()):
            for rank, r in enumerate(sorted(matches, key=lambda x: -float(x[8]))[:10], 1):
                subject = r[1]
                if species == 'Gloydius_ussuriensis':
                    gid = re.sub(r'-T\d+$', '', subject)
                    full_gid, chrom, start, end, strand, mapping = map_gene(subject, guss_index)
                    if mapping == 'unmapped_GFF_subject':
                        full_gid, chrom, start, end, strand, mapping = map_gene(gid, guss_index)
                    if mapping != 'unmapped_GFF_subject':
                        gid = full_gid
                        mapping = 'full_GFF_' + mapping
                    else:
                        g = flanks.get(gid, {})
                        chrom, start, end, strand = (g.get('chr', ''), g.get('start', ''), g.get('end', ''), g.get('strand', ''))
                        mapping = 'Step2_gene_model' if g else 'unmapped_full_GFF_and_Step2'
                else:
                    gid, chrom, start, end, strand, mapping = map_gene(subject, index[species])
                reference_candidates.append(dict(family=family, anchor_query=query, species=species, subject_protein_id=subject, hit_rank=rank, identity_pct=r[2], aligned_aa=r[3], query_length_aa=r[4], subject_length_aa=r[5], query_coverage_pct=r[6], evalue=r[7], bitscore=r[8], gene_id=gid, chr=chrom, start=start, end=end, strand=strand, mapping_status=mapping, validation_status='candidate_requires_review'))

write_tsv(evidence / 'anchor_reference_blastp_candidates.tsv', reference_candidates, ['family', 'anchor_query', 'species', 'subject_protein_id', 'hit_rank', 'identity_pct', 'aligned_aa', 'query_length_aa', 'subject_length_aa', 'query_coverage_pct', 'evalue', 'bitscore', 'gene_id', 'chr', 'start', 'end', 'strand', 'mapping_status', 'validation_status'])

# Genomic TBLASTN source sequences are extracted in reference orientation; convert local to 1-based genome coordinates.
regions = {r['family']: r for r in read_tsv(evidence / 'Guss_target_intervals.tsv')}
all_genome_hits = []
for family, region in regions.items():
    offset = int(region['start']) - 1
    path = evidence / f'Guss_{family}_anchors.tblastn.tsv'
    with path.open() as handle:
        for r in csv.reader(handle, delimiter='\t'):
            if len(r) != 11: continue
            gs, ge = sorted((int(r[7]) + offset, int(r[8]) + offset))
            overlapping = sorted(g['gene_id'] for g in flanks.values() if g['family'] == family and g['chr'] == region['chromosome'] and int(g['start']) <= ge and int(g['end']) >= gs)
            full_overlapping = overlapping_genes(region['chromosome'], gs, ge)
            distance, nearest = nearest_gene(region['chromosome'], gs, ge)
            all_genome_hits.append(dict(family=family, query_anchor=r[0], chr=region['chromosome'], start=gs, end=ge, strand='+' if int(r[7]) <= int(r[8]) else '-', query_start_aa=int(r[5]), query_end_aa=int(r[6]), query_length_aa=int(r[4]), aligned_aa=int(r[3]), percent_identity=r[2], evalue=r[9], bitscore=r[10], overlapping_FUN_gene_ids=';'.join(overlapping), overlapping_full_GFF_gene_ids=';'.join(g['gene_id'] for g in full_overlapping), nearest_full_GFF_gene_id=nearest['gene_id'] if nearest else '', nearest_full_GFF_distance_bp=distance if distance is not None else '', validation_status='tblastn_HSP_not_gene_assignment'))

all_genome_hits.sort(key=lambda x: (x['family'], x['query_anchor'], -float(x['bitscore'])))
write_tsv(evidence / 'Guss_anchor_tblastn_genomic_hits.tsv', all_genome_hits, ['family', 'query_anchor', 'chr', 'start', 'end', 'strand', 'query_start_aa', 'query_end_aa', 'query_length_aa', 'aligned_aa', 'percent_identity', 'evalue', 'bitscore', 'overlapping_FUN_gene_ids', 'overlapping_full_GFF_gene_ids', 'nearest_full_GFF_gene_id', 'nearest_full_GFF_distance_bp', 'validation_status'])

def unique_query_coverage(hits):
    intervals = sorted((min(int(h['query_start_aa']), int(h['query_end_aa'])), max(int(h['query_start_aa']), int(h['query_end_aa']))) for h in hits)
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return sum(end - start + 1 for start, end in merged)

def summarize_cluster(hits):
    hits = sorted(hits, key=lambda h: (int(h['start']), int(h['end'])))
    first = hits[0]
    start, end = min(int(h['start']) for h in hits), max(int(h['end']) for h in hits)
    positions = hits if first['strand'] == '+' else list(reversed(hits))
    qstarts = [int(h['query_start_aa']) for h in positions]
    collinear = all(b >= a - 20 for a, b in zip(qstarts, qstarts[1:]))
    covered_aa = unique_query_coverage(hits)
    query_len = int(first['query_length_aa'])
    models = overlapping_genes(first['chr'], start, end)
    distance, nearest = nearest_gene(first['chr'], start, end)
    return dict(family=first['family'], query_anchor=first['query_anchor'], chr=first['chr'], start=start, end=end, strand=first['strand'], span_bp=end-start+1, hsp_count=len(hits), unique_query_aa=covered_aa, query_length_aa=query_len, unique_query_coverage_pct=round(100*covered_aa/query_len, 2) if query_len else '', collinear_query_order='yes' if collinear else 'no', highest_hsp_identity_pct=max(float(h['percent_identity']) for h in hits), best_hsp_bitscore=max(float(h['bitscore']) for h in hits), sum_hsp_bitscore=round(sum(float(h['bitscore']) for h in hits), 1), overlapping_full_GFF_gene_ids=';'.join(g['gene_id'] for g in models), full_GFF_model_status='overlaps_annotated_gene' if models else 'no_overlapping_annotated_gene', nearest_full_GFF_gene_id=nearest['gene_id'] if nearest else '', nearest_full_GFF_distance_bp=distance if distance is not None else '', validation_status='candidate_HSP_cluster_requires_manual_review')

clusters = []
by_anchor = defaultdict(list)
for h in all_genome_hits:
    by_anchor[(h['family'], h['query_anchor'], h['chr'], h['strand'])].append(h)
for key, hits in by_anchor.items():
    ordered = sorted(hits, key=lambda h: (int(h['start']), int(h['end'])))
    group = []
    current_end = -1
    for h in ordered:
        if group and int(h['start']) - current_end > 7500:
            clusters.append(summarize_cluster(group))
            group = []
        group.append(h)
        current_end = max(current_end, int(h['end']))
    if group: clusters.append(summarize_cluster(group))
clusters.sort(key=lambda c: (c['family'], c['query_anchor'], c['chr'], c['start']))
cluster_cols = ['family', 'query_anchor', 'chr', 'start', 'end', 'strand', 'span_bp', 'hsp_count', 'unique_query_aa', 'query_length_aa', 'unique_query_coverage_pct', 'collinear_query_order', 'highest_hsp_identity_pct', 'best_hsp_bitscore', 'sum_hsp_bitscore', 'overlapping_full_GFF_gene_ids', 'full_GFF_model_status', 'nearest_full_GFF_gene_id', 'nearest_full_GFF_distance_bp', 'validation_status']
write_tsv(evidence / 'Guss_anchor_tblastn_candidate_loci.tsv', clusters, cluster_cols)
adam_clusters = [c for c in clusters if c['family'] == 'SVMP' and c['query_anchor'].split('_', 1)[-1] in ('ADAM', 'ADAM28')]
write_tsv(evidence / 'Guss_SVMP_ADAM28_candidate_loci.tsv', adam_clusters, cluster_cols)

# Anchor validation checklist: no automatic orthology claims based on one BLAST HSP.
expected = {'SVMP': ['STC1', 'NEFM', 'NEFL', 'ADAM'], 'PLA2': ['OTUD3', 'MUL1', 'VWA5B1', 'PLA2G2F']}
reference_queries = {}
for family in ('SVMP', 'PLA2'):
    with (evidence / f'{family}_anchor_reference_proteins.faa').open() as handle:
        reference_queries[family] = set(line[1:].strip().split()[0] for line in handle if line.startswith('>'))

def matches_anchor(query, symbol):
    label = query.split('_', 1)[-1]
    return label == symbol or (symbol == 'ADAM' and label == 'ADAM28')

review = []
for family, symbols in expected.items():
    for symbol in symbols:
        candidates = sorted(gid for gid, sym in known.items() if sym == symbol and gid in flanks and flanks[gid]['family'] == family)
        queries = sorted(q for q in reference_queries[family] if matches_anchor(q, symbol))
        hits = [r for r in all_genome_hits if r['family'] == family and r['query_anchor'] in queries]
        loci = [c for c in clusters if c['family'] == family and c['query_anchor'] in queries]
        ranked_loci = sorted(loci, key=lambda c: (c['collinear_query_order'] == 'yes', int(c['unique_query_aa']), int(c['hsp_count']), float(c['sum_hsp_bitscore'])), reverse=True)
        best_locus = ranked_loci[0] if ranked_loci else {}
        loc = f"{best_locus['chr']}:{best_locus['start']}-{best_locus['end']}" if best_locus else ''
        for species in ('Gloydius_ussuriensis', 'Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
            p = [r for r in reference_candidates if r['species'] == species and r['family'] == family and r['anchor_query'] in queries]
            p.sort(key=lambda r: -float(r['bitscore']))
            best = p[0] if p else {}
            review.append(dict(family=family, species=species, anchor_symbol=symbol, Gloydius_model_hint=';'.join(candidates) if species == 'Gloydius_ussuriensis' else '', top_protein_id=best.get('subject_protein_id', ''), top_gene_id=best.get('gene_id', ''), top_blastp_identity=best.get('identity_pct', ''), top_blastp_query_coverage=best.get('query_coverage_pct', ''), reference_tblastn_queries=';'.join(queries) if species == 'Gloydius_ussuriensis' else '', tblastn_hsp_count=len(hits) if species == 'Gloydius_ussuriensis' else '', best_tblastn_bitscore=max((float(r['bitscore']) for r in hits), default=0) if species == 'Gloydius_ussuriensis' else '', top_tblastn_candidate_locus=loc if species == 'Gloydius_ussuriensis' else '', top_locus_collinear_query_order=best_locus.get('collinear_query_order', '') if species == 'Gloydius_ussuriensis' else '', top_locus_unique_query_coverage_pct=best_locus.get('unique_query_coverage_pct', '') if species == 'Gloydius_ussuriensis' else '', top_locus_full_GFF_gene_ids=best_locus.get('overlapping_full_GFF_gene_ids', '') if species == 'Gloydius_ussuriensis' else '', top_locus_GFF_status=best_locus.get('full_GFF_model_status', '') if species == 'Gloydius_ussuriensis' else '', validation_status='manual_homology_and_synteny_review_required'))

write_tsv(evidence / 'flanking_anchor_review_checklist.tsv', review, ['family', 'species', 'anchor_symbol', 'Gloydius_model_hint', 'top_protein_id', 'top_gene_id', 'top_blastp_identity', 'top_blastp_query_coverage', 'reference_tblastn_queries', 'tblastn_hsp_count', 'best_tblastn_bitscore', 'top_tblastn_candidate_locus', 'top_locus_collinear_query_order', 'top_locus_unique_query_coverage_pct', 'top_locus_full_GFF_gene_ids', 'top_locus_GFF_status', 'validation_status'])
print('Forward BLASTp + reciprocal + GFF evidence rows:', len(forward))
print('Targeted genomic TBLASTN HSPs:', len(all_genome_hits))
print('Reference anchored BLASTp candidate HSPs:', len(reference_candidates))
print('Gloydius candidate TBLASTN loci:', len(clusters))
print('ADAM/ADAM28 candidate loci:', len(adam_clusters))
print('Checklist entries:', len(review))

# Detailed subject-coordinate comparisons for potentially split or tandem PLA2 models.
pla2_alignments = []
guss_by_id = {g['gene_id']: g for g in guss_genes}
for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
    path = evidence / f'PLA2_discrepancy_vs_{species}.blastp.tsv'
    with path.open() as handle:
        for r in csv.reader(handle, delimiter='\t'):
            if len(r) != 13: continue
            qid, subject = r[0], r[1]
            gene = re.sub(r'-T\d+$', '', qid)
            g = guss_by_id.get(gene)
            sgid, schrom, ss, se, sstrand, smap = map_gene(subject, index[species])
            pla2_alignments.append(dict(species=species, Gloydius_gene_id=gene, Gloydius_transcript=qid, Gloydius_chr=g['chr'] if g else '', Gloydius_start=g['start'] if g else '', Gloydius_end=g['end'] if g else '', Gloydius_strand=g['strand'] if g else '', subject_protein_id=subject, subject_gene_id=sgid, subject_chr=schrom, subject_genomic_start=ss, subject_genomic_end=se, subject_mapping_status=smap, identity_pct=r[2], aligned_aa=r[3], query_length_aa=r[4], subject_length_aa=r[5], query_start_aa=r[6], query_end_aa=r[7], subject_start_aa=r[8], subject_end_aa=r[9], query_coverage_pct=r[10], evalue=r[11], bitscore=r[12]))
pla2_alignments.sort(key=lambda x: (x['species'], x['Gloydius_gene_id'], -float(x['bitscore'])))
aln_cols = ['species', 'Gloydius_gene_id', 'Gloydius_transcript', 'Gloydius_chr', 'Gloydius_start', 'Gloydius_end', 'Gloydius_strand', 'subject_protein_id', 'subject_gene_id', 'subject_chr', 'subject_genomic_start', 'subject_genomic_end', 'subject_mapping_status', 'identity_pct', 'aligned_aa', 'query_length_aa', 'subject_length_aa', 'query_start_aa', 'query_end_aa', 'subject_start_aa', 'subject_end_aa', 'query_coverage_pct', 'evalue', 'bitscore']
write_tsv(evidence / 'PLA2_model_discrepancy_alignments.tsv', pla2_alignments, aln_cols)

reconciliation = []
pairs = [('VWA5B1_like', 'FUN_018053', 'FUN_018054'), ('VWA5B1_like', 'FUN_018053', 'FUN_018055'), ('VWA5B1_like', 'FUN_018054', 'FUN_018055'), ('downstream_pair', 'FUN_018074', 'FUN_018075')]
for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
    rows = [r for r in pla2_alignments if r['species'] == species]
    best_per_gene_subject = {}
    best_per_gene = {}
    for r in rows:
        key = (r['Gloydius_gene_id'], r['subject_protein_id'])
        if key not in best_per_gene_subject or float(r['bitscore']) > float(best_per_gene_subject[key]['bitscore']): best_per_gene_subject[key] = r
        if r['Gloydius_gene_id'] not in best_per_gene or float(r['bitscore']) > float(best_per_gene[r['Gloydius_gene_id']]['bitscore']): best_per_gene[r['Gloydius_gene_id']] = r
    for label, a, b in pairs:
        aa = [r for (g, subject), r in best_per_gene_subject.items() if g == a]
        bb = [r for (g, subject), r in best_per_gene_subject.items() if g == b]
        common = set(r['subject_protein_id'] for r in aa) & set(r['subject_protein_id'] for r in bb)
        if common:
            shared = max(common, key=lambda p: float(best_per_gene_subject[(a, p)]['bitscore']) + float(best_per_gene_subject[(b, p)]['bitscore']))
            ra, rb = best_per_gene_subject[(a, shared)], best_per_gene_subject[(b, shared)]
            a1, a2 = sorted((int(ra['subject_start_aa']), int(ra['subject_end_aa'])))
            b1, b2 = sorted((int(rb['subject_start_aa']), int(rb['subject_end_aa'])))
            overlap = max(0, min(a2, b2) - max(a1, b1) + 1)
            shorter = min(a2-a1+1, b2-b1+1)
            overlap_pct = round(100 * overlap / shorter, 2) if shorter else ''
            subject_span_union = (a2-a1+1) + (b2-b1+1) - overlap
            evidence_flag = 'shared_subject_overlapping_segments' if overlap else 'shared_subject_complementary_segments'
        else:
            ra, rb = best_per_gene.get(a, {}), best_per_gene.get(b, {})
            shared, overlap, overlap_pct, subject_span_union = '', '', '', ''
            evidence_flag = 'no_shared_subject_in_top_30' if ra and rb else 'missing_one_or_both_query_proteins_or_hits'
        ga, gb = guss_by_id.get(a), guss_by_id.get(b)
        if ga and gb and ga['chr'] == gb['chr']:
            genomic_gap = max(0, ga['start'] - gb['end'] - 1, gb['start'] - ga['end'] - 1)
            genomic_overlap = max(0, min(ga['end'], gb['end']) - max(ga['start'], gb['start']) + 1)
        else:
            genomic_gap, genomic_overlap = '', ''
        if evidence_flag == 'shared_subject_complementary_segments':
            interpretation = 'possible_fragmented_annotation_or_divergent_paralogs_review_exons'
        elif evidence_flag == 'shared_subject_overlapping_segments':
            interpretation = 'shared_protein_region_possible_independent_paralogs_or_redundant_models'
        else:
            interpretation = 'unresolved_check_reference_and_gene_models'
        reconciliation.append(dict(group=label, species=species, Gloydius_gene_a=a, Gloydius_gene_b=b, Gloydius_genomic_gap_bp=genomic_gap, Gloydius_gene_overlap_bp=genomic_overlap, gene_a_best_subject=best_per_gene.get(a, {}).get('subject_protein_id', ''), gene_b_best_subject=best_per_gene.get(b, {}).get('subject_protein_id', ''), shared_subject=shared, gene_a_subject_aa=f"{ra['subject_start_aa']}-{ra['subject_end_aa']}" if ra else '', gene_b_subject_aa=f"{rb['subject_start_aa']}-{rb['subject_end_aa']}" if rb else '', shared_subject_aa_overlap=overlap, overlap_pct_of_shorter_subject_segment=overlap_pct, subject_aa_union_span=subject_span_union, evidence_flag=evidence_flag, interpretation=interpretation, validation_status='manual_model_review_required'))
reconcile_cols = ['group', 'species', 'Gloydius_gene_a', 'Gloydius_gene_b', 'Gloydius_genomic_gap_bp', 'Gloydius_gene_overlap_bp', 'gene_a_best_subject', 'gene_b_best_subject', 'shared_subject', 'gene_a_subject_aa', 'gene_b_subject_aa', 'shared_subject_aa_overlap', 'overlap_pct_of_shorter_subject_segment', 'subject_aa_union_span', 'evidence_flag', 'interpretation', 'validation_status']
write_tsv(evidence / 'PLA2_model_reconciliation.tsv', reconciliation, reconcile_cols)

# One row per species x PLA2 anchor for transfer to Steps 6 and 7.
pla2_presence = []
for symbol in ('OTUD3', 'MUL1', 'VWA5B1', 'PLA2G2F'):
    queries = sorted(q for q in reference_queries['PLA2'] if matches_anchor(q, symbol))
    for species in ('Gloydius_ussuriensis', 'Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        hits = sorted((r for r in reference_candidates if r['family'] == 'PLA2' and r['species'] == species and r['anchor_query'] in queries), key=lambda r: -float(r['bitscore']))
        global_top = hits[0] if hits else {}
        target = regions['PLA2']
        in_target = [r for r in hits if species == 'Gloydius_ussuriensis' and r['chr'] == target['chromosome'] and r['start'] and r['end'] and int(r['start']) <= int(target['end']) and int(r['end']) >= int(target['start'])]
        top = in_target[0] if in_target else global_top
        target_status = 'inside_search_interval' if in_target else ('outside_or_unmapped_search_interval' if species == 'Gloydius_ussuriensis' and hits else '')
        loci = sorted((c for c in clusters if c['family'] == 'PLA2' and c['query_anchor'] in queries), key=lambda c: (c['collinear_query_order']=='yes', float(c['unique_query_coverage_pct']), float(c['sum_hsp_bitscore'])), reverse=True)
        locus = loci[0] if loci else {}
        if not queries: status = 'reference_not_recovered'
        elif species == 'Gloydius_ussuriensis' and hits and not in_target: status = 'protein_hit_outside_or_unmapped_target_region_review'
        elif top.get('gene_id') and top.get('chr'): status = 'protein_hit_with_mapped_gene_model'
        elif hits: status = 'protein_hit_mapping_unresolved'
        elif species == 'Gloydius_ussuriensis' and locus: status = 'genomic_tblastn_candidate_only'
        else: status = 'no_detected_protein_candidate'
        pla2_presence.append(dict(family='PLA2', anchor_symbol=symbol, species=species, reference_queries=';'.join(queries), evidence_status=status, top_subject_protein=top.get('subject_protein_id', ''), top_subject_gene_id=top.get('gene_id', ''), top_subject_chr=top.get('chr', ''), top_subject_start=top.get('start', ''), top_subject_end=top.get('end', ''), protein_identity_pct=top.get('identity_pct', ''), protein_query_coverage_pct=top.get('query_coverage_pct', ''), protein_hit_count=len(hits), top_protein_target_region_status=target_status, global_best_hit_protein=global_top.get('subject_protein_id', ''), Gloydius_top_tblastn_locus=f"{locus['chr']}:{locus['start']}-{locus['end']}" if species=='Gloydius_ussuriensis' and locus else '', Gloydius_locus_query_coverage_pct=locus.get('unique_query_coverage_pct', '') if species=='Gloydius_ussuriensis' else '', Gloydius_locus_GFF_models=locus.get('overlapping_full_GFF_gene_ids', '') if species=='Gloydius_ussuriensis' else '', validation_status='synteny_and_orthology_review_required'))
pla2_presence_cols = ['family', 'anchor_symbol', 'species', 'reference_queries', 'evidence_status', 'top_subject_protein', 'top_subject_gene_id', 'top_subject_chr', 'top_subject_start', 'top_subject_end', 'protein_identity_pct', 'protein_query_coverage_pct', 'protein_hit_count', 'top_protein_target_region_status', 'global_best_hit_protein', 'Gloydius_top_tblastn_locus', 'Gloydius_locus_query_coverage_pct', 'Gloydius_locus_GFF_models', 'validation_status']
write_tsv(evidence / 'PLA2_anchor_presence_matrix.tsv', pla2_presence, pla2_presence_cols)
print('PLA2 detailed alignment HSPs:', len(pla2_alignments))
print('PLA2 pairwise model-review comparisons:', len(reconciliation))
print('PLA2 anchor presence rows:', len(pla2_presence))

print('NOTE: no ortholog identities have been automatically validated.')
PY

# 11. best-hit report
printf '\n===== BEST FORWARD BLASTP HITS =====\n'
python - "$BLASTDIR" <<'PY'
import csv, sys
from pathlib import Path
from collections import defaultdict
blastdir = Path(sys.argv[1])
for family in ('SVMP', 'PLA2'):
    print(f'\n===== {family} =====')
    for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        path = blastdir / family / f'Guss_{family}_vs_{species}.blastp.tsv'
        hits = defaultdict(list)
        with path.open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) >= 9: hits[row[0]].append(row)
        print(f'\n--- {species} ---')
        if not hits: print('No hits passing threshold.')
        for query in sorted(hits):
            best = max(hits[query], key=lambda r: float(r[8]))
            print(f'{query:<20} {best[1]:<28} ID={float(best[2]):5.1f}% Qcov={float(best[6]):5.1f}% E={best[7]}')
print('\n===== COMPLETE =====')
PY

echo
printf 'Additional Step 4 evidence: %s\n' "$EVIDENCE"
printf 'Step 6/7 input for review: %s\n' "$EVIDENCE/flanking_blastp_reciprocal_gff_evidence.tsv"
printf 'Genomic candidates: %s\n' "$EVIDENCE/Guss_anchor_tblastn_genomic_hits.tsv"
printf 'Anchor checklist: %s\n' "$EVIDENCE/flanking_anchor_review_checklist.tsv"
printf 'ADAM28 candidate loci: %s\n' "$EVIDENCE/Guss_SVMP_ADAM28_candidate_loci.tsv"
printf 'PLA2 anchor selection: %s\n' "$EVIDENCE/PLA2_anchor_reference_audit.tsv"
printf 'PLA2 model alignments: %s\n' "$EVIDENCE/PLA2_model_discrepancy_alignments.tsv"
printf 'PLA2 model reconciliation: %s\n' "$EVIDENCE/PLA2_model_reconciliation.tsv"
