#!/bin/bash
#SBATCH --job-name=toxin_all_QC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
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

# fix shared-library compatibility
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

# ------------------------------------------------------------
# directories and inputs
# ------------------------------------------------------------

ROOT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
OUTDIR="${WORKDIR}/all120_QC"

GTF="${WORKDIR}/toxin_annotation.gtf"
PEP="${WORKDIR}/toxin_annotation_pep.fasta"
CDS="${WORKDIR}/toxin_annotation_cds.fasta"

RNA_BAM="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.sorted.bam"
PFAM_HMM="${WORKDIR}/protein_validation/pfam_targeted/toxin_families.hmm"

MODEL_JSON="${OUTDIR}/model_metadata.json"
CDS_BED="${OUTDIR}/all_CDS_union.bed"
GENE_BED="${OUTDIR}/all_gene_spans.bed"
SP_FASTA="${OUTDIR}/signalp_input.fasta"
HMM_FASTA="${OUTDIR}/hmm_fragments.fasta"
DEPTH_TSV="${OUTDIR}/venom_CDS_depth.tsv"
HMM_TBL="${OUTDIR}/targeted_pfam.domtblout"
SIGNAL_DIR="${OUTDIR}/signalp6_${SLURM_JOB_ID:-manual}"

SUMMARY="${OUTDIR}/toxin_all120_QC.tsv"
FLAGGED="${OUTDIR}/toxin_all120_QC_flagged.tsv"
JUNCTIONS="${OUTDIR}/toxin_all120_RNA_junctions.tsv"
DOMAINS="${OUTDIR}/toxin_all120_targeted_domains.tsv"

mkdir -p "$OUTDIR"

# validate inputs and programs
for file in "$GTF" "$PEP" "$CDS" "$RNA_BAM"; do
    [[ -s "$file" ]] || {
        echo "ERROR: Missing input: $file" >&2
        exit 1
    }
done

for program in python samtools signalp6; do
    command -v "$program" >/dev/null || {
        echo "ERROR: $program unavailable" >&2
        exit 1
    }
done

python -c 'from Bio.Seq import Seq; from PIL import Image; import matplotlib; import signalp; print("Python and SignalP imports: PASS")'


# ------------------------------------------------------------
# STEP 1: GTF / FASTA integrity and model inventory
# ------------------------------------------------------------

echo '===== STEP 1: BUILD 120-GENE MODEL INVENTORY ====='

python - "$GTF" "$PEP" "$CDS" "$MODEL_JSON" "$CDS_BED" "$GENE_BED" "$SP_FASTA" "$HMM_FASTA" <<'PY'
import sys, re, csv, json, hashlib
from collections import defaultdict
from Bio.Seq import Seq

gtf, pepfile, cdsfile, jsonfile, cdsbed, genebed, spfile, hmmfile = sys.argv[1:]

def fasta_records(path):
    records, name, parts = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, ''.join(parts).upper()))
                name, parts = line[1:], []
            else:
                parts.append(line)
    if name is not None:
        records.append((name, ''.join(parts).upper()))
    return records

def attribute(text, key):
    m = re.search(r'(?:^|;)\s*' + re.escape(key) + r'\s+"([^"]+)"', text)
    return m.group(1) if m else None

def merge(intervals):
    result = []
    for a, b in sorted(intervals):
        if result and a <= result[-1][1] + 1:
            result[-1][1] = max(result[-1][1], b)
        else:
            result.append([a, b])
    return result

def choose_record(records, gene):
    pattern = re.compile(r'(?<![A-Za-z0-9])' + re.escape(gene) + r'(?![A-Za-z0-9])')
    hits = [(head, seq) for head, seq in records if pattern.search(head)]
    return (hits[0][1] if len(hits) == 1 else None), len(hits)

models = defaultdict(lambda: defaultdict(list))

with open(gtf) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue
        cols = line.rstrip('\n').split('\t')
        if len(cols) < 9 or cols[2] != 'CDS':
            continue
        gene = attribute(cols[8], 'gene_id')
        if not gene:
            raise RuntimeError('CDS without gene_id: ' + line[:120])
        transcript = attribute(cols[8], 'transcript_id') or gene
        models[gene][transcript].append(dict(chrom=cols[0], start=int(cols[3]), end=int(cols[4]), strand=cols[6], phase=cols[7]))

pep_records, cds_records = fasta_records(pepfile), fasta_records(cdsfile)
output, cds_intervals, gene_spans = [], defaultdict(list), []
sp_count, hmm_count = 0, 0

with open(spfile, 'w') as sp, open(hmmfile, 'w') as hf:
    for index, gene in enumerate(sorted(models), 1):
        isoforms = models[gene]

        # Flag multiple isoforms and screen the longest CDS model.
        def score(item):
            return sum(b['end'] - b['start'] + 1 for b in item[1])

        tx_id, features = sorted(isoforms.items(), key=lambda x: (-score(x), x[0]))[0]

        chroms = set(x['chrom'] for x in features)
        strands = set(x['strand'] for x in features)

        if len(chroms) != 1 or len(strands) != 1:
            raise RuntimeError('Mixed chromosome/strand in ' + gene)

        chrom, strand = features[0]['chrom'], features[0]['strand']
        ordered = sorted(features, key=lambda x: x['start'], reverse=(strand == '-'))
        segments = merge((x['start'], x['end']) for x in features)
        cds_length = sum(b - a + 1 for a, b in segments)
        gtf_len = sum(x['end'] - x['start'] + 1 for x in features)
        junc = [[segments[j][1], segments[j+1][0]] for j in range(len(segments) - 1)]

        phases = [x['phase'] for x in ordered]
        phase_status = 'NOT_AVAILABLE'

        if all(x in ('0', '1', '2') for x in phases):
            p0 = int(phases[0])
            consumed = ordered[0]['end'] - ordered[0]['start'] + 1 - p0
            phase_status = 'PASS'

            for segment in ordered[1:]:
                expected = (-consumed) % 3
                if int(segment['phase']) != expected:
                    phase_status = 'FAIL'
                consumed += segment['end'] - segment['start'] + 1

        protein, protein_hits = choose_record(pep_records, gene)
        dna, cds_hits = choose_record(cds_records, gene)

        translated_status = 'NOT_ASSESSED'
        dna_frame = 'NA'

        if dna is not None:
            dna = dna.replace('U', 'T')
            dna_frame = len(dna) % 3

            if protein is not None:
                translated = str(Seq(dna[:(len(dna) // 3) * 3]).translate())
                translated_status = 'MATCH' if translated.rstrip('*') == protein.rstrip('*') else 'MISMATCH'

        if protein is not None:
            without_final = protein[:-1] if protein.endswith('*') else protein

            # Distinguish genuine internal stops from tandem terminal stops.
            internal_stops = without_final.rstrip('*').count('*')
            terminal_stop_run = len(protein) - len(protein.rstrip('*'))
            first_internal_stop = without_final.find('*') + 1 if internal_stops else None

            if internal_stops:
                stop_class = 'INTERNAL_STOP'
            elif terminal_stop_run >= 2:
                stop_class = 'TANDEM_TERMINAL_STOP'
            elif terminal_stop_run == 1:
                stop_class = 'SINGLE_TERMINAL_STOP'
            else:
                stop_class = 'NO_TERMINAL_STOP_IN_FASTA'

            length_aa = sum(x != '*' for x in protein)
            prefix = protein.split('*')[0]

            if len(prefix) >= 20:
                sp_id = 'Q{:05d}'.format(index)
                sp.write('>{}\n{}\n'.format(sp_id, prefix))
                sp_count += 1
            else:
                sp_id = None

            for number, part in enumerate(without_final.split('*'), 1):
                if len(part) >= 12:
                    hf.write('>Q{:05d}F{:03d}\n{}\n'.format(index, number, part))
                    hmm_count += 1

        else:
            internal_stops, terminal_stop_run, first_internal_stop = None, None, None
            stop_class, length_aa, sp_id = 'UNMATCHED_PROTEIN', None, None

        model = dict(
            gene=gene,
            chrom=chrom,
            strand=strand,
            transcript=tx_id,
            isoforms=len(isoforms),
            segments=segments,
            cds_bp=cds_length,
            gtf_cds_feature_bp=gtf_len,
            junctions=junc,
            phases=phases,
            phase_status=phase_status,
            protein=protein,
            protein_records=protein_hits,
            protein_length_aa=length_aa,
            internal_stops=internal_stops,
            terminal_stop_run=terminal_stop_run,
            first_internal_stop=first_internal_stop,
            stop_class=stop_class,
            cds_records=cds_hits,
            cds_fasta_bp=len(dna) if dna else None,
            cds_frame_mod3=dna_frame,
            translation_status=translated_status,
            signalp_id=sp_id,
            protein_sha256=hashlib.sha256(protein.rstrip('*').encode()).hexdigest() if protein else None,
            query_id='Q{:05d}'.format(index)
        )

        output.append(model)
        cds_intervals[chrom].extend(segments)
        gene_spans.append((chrom, min(a for a, b in segments) - 1, max(b for a, b in segments)))

with open(cdsbed, 'w') as f:
    for chrom in sorted(cds_intervals):
        for start, end in merge(cds_intervals[chrom]):
            f.write('{}\t{}\t{}\n'.format(chrom, start - 1, end))

gene_span_union = defaultdict(list)

for chrom, start0, end in gene_spans:
    gene_span_union[chrom].append((start0 + 1, end))

with open(genebed, 'w') as f:
    for chrom in sorted(gene_span_union):
        for a, b in merge(gene_span_union[chrom]):
            f.write('{}\t{}\t{}\n'.format(chrom, a - 1, b))

with open(jsonfile, 'w') as f:
    json.dump(output, f, indent=2)

print('GTF genes:', len(output), '| protein FASTA records:', len(pep_records), '| CDS FASTA records:', len(cds_records))
print('SignalP inputs:', sp_count, '| HMMER protein fragments:', hmm_count)

if len(output) != 120:
    print('WARNING: GTF contained {} distinct genes, not 120.'.format(len(output)))
PY


# ------------------------------------------------------------
# STEP 2: venom RNA-seq coverage
# ------------------------------------------------------------

echo '===== STEP 2: RNA-SEQ COVERAGE ====='

samtools depth \
    -b "$CDS_BED" \
    -q 20 \
    -Q 20 \
    -G 0xF04 \
    "$RNA_BAM" > "$DEPTH_TSV"


# ------------------------------------------------------------
# STEP 3: SignalP 6
# ------------------------------------------------------------

echo '===== STEP 3: SIGNALP 6 ====='

if grep -q '^>' "$SP_FASTA"; then

    signalp6 \
        --fastafile "$SP_FASTA" \
        --organism eukarya \
        --output_dir "$SIGNAL_DIR" \
        --format none \
        --mode fast \
        --write_procs 1 \
        --torch_num_threads "${SLURM_CPUS_PER_TASK:-4}"

    [[ -s "${SIGNAL_DIR}/prediction_results.txt" ]] || {
        echo 'ERROR: SignalP result missing' >&2
        exit 1
    }

else
    echo 'WARNING: No eligible SignalP input proteins.'
    mkdir -p "$SIGNAL_DIR"
    : > "${SIGNAL_DIR}/prediction_results.txt"
fi


# ------------------------------------------------------------
# STEP 4: Targeted Pfam analysis
#
# Five families:
# SVSP, Ficolin, PLA2, CTL, HYAL
# ------------------------------------------------------------

echo '===== STEP 4: TARGETED PFAM ====='

if [[ -s "$PFAM_HMM" ]]; then

    command -v hmmscan >/dev/null || {
        echo 'ERROR: hmmscan not available' >&2
        exit 1
    }

    if [[ ! -s "${PFAM_HMM}.h3p" ]]; then
        hmmpress -f "$PFAM_HMM" > "${OUTDIR}/hmmpress.log"
    fi

    hmmscan \
        --cpu "${SLURM_CPUS_PER_TASK:-4}" \
        --noali \
        --domtblout "$HMM_TBL" \
        "$PFAM_HMM" \
        "$HMM_FASTA" > "${OUTDIR}/hmmscan.log"

else
    echo 'WARNING: Five-family Pfam HMM not found; domain columns will be NOT_ASSESSED.'
    : > "$HMM_TBL"
fi


# ------------------------------------------------------------
# STEP 5: integrate all QC results
# ------------------------------------------------------------

echo '===== STEP 5: INTEGRATED QC ====='

python - "$MODEL_JSON" "$DEPTH_TSV" "$RNA_BAM" "$GENE_BED" "${SIGNAL_DIR}/prediction_results.txt" "$HMM_TBL" "$SUMMARY" "$FLAGGED" "$JUNCTIONS" "$DOMAINS" <<'PY'
import sys, re, os, csv, json, subprocess
from collections import defaultdict, Counter

jsonfile, depthfile, bam, gene_bed, spfile, domfile, summary_file, flagged_file, junction_file, domain_file = sys.argv[1:]

models = json.load(open(jsonfile))

# ------------------------------------------------------------
# RNA base coverage
# ------------------------------------------------------------

intervals = defaultdict(list)

for m in models:
    for a, b in m['segments']:
        intervals[m['chrom']].append((a, b, m['gene']))

for chrom in intervals:
    intervals[chrom].sort()

counts = defaultdict(lambda: [0, 0])
state = {}

with open(depthfile) as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')

        if len(fields) < 3:
            continue

        chrom, pos, dep = fields[0], int(fields[1]), int(fields[2])

        if chrom not in state:
            state[chrom] = [0, []]

        pointer, active = state[chrom]
        source = intervals.get(chrom, [])

        while pointer < len(source) and source[pointer][0] <= pos:
            active.append(source[pointer])
            pointer += 1

        active = [x for x in active if x[1] >= pos]

        for a, b, gene in active:
            counts[gene][0] += 1
            counts[gene][1] += dep

        state[chrom] = [pointer, active]

# ------------------------------------------------------------
# Annotated splice-junction support
# ------------------------------------------------------------

junction_map = defaultdict(list)

for m in models:
    for left, right in m['junctions']:
        junction_map[(m['chrom'], left, right)].append(m['gene'])

junction_counts = Counter()

proc = subprocess.Popen(
    ['samtools', 'view', '-L', gene_bed, '-q', '20', '-F', '0xF04', bam],
    stdout=subprocess.PIPE,
    universal_newlines=True
)

for line in proc.stdout:
    f = line.rstrip('\n').split('\t')

    if len(f) < 11 or f[5] == '*':
        continue

    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', f[5])]

    if ''.join(str(n)+op for n, op in ops) != f[5]:
        continue

    refpos = int(f[3])

    for i, (n, op) in enumerate(ops):

        if op in ('M', 'D', '=', 'X'):
            refpos += n

        elif op == 'N':
            left, right = refpos - 1, refpos + n

            if (
                i > 0 and i + 1 < len(ops)
                and ops[i-1][1] in ('M', '=', 'X')
                and ops[i+1][1] in ('M', '=', 'X')
                and ops[i-1][0] >= 8
                and ops[i+1][0] >= 8
            ):
                for gene in junction_map.get((f[2], left, right), []):
                    junction_counts[(gene, left, right)] += 1

            refpos += n

if proc.wait() != 0:
    raise RuntimeError('samtools view failed')

# ------------------------------------------------------------
# SignalP 6 results
# ------------------------------------------------------------

signal = {}

with open(spfile) as f:
    for line in f:
        if not line.strip() or line.startswith('#'):
            continue

        fields = line.rstrip('\n').split('\t')

        if len(fields) < 4:
            raise RuntimeError('Malformed SignalP line: ' + line[:120])

        if fields[0] in signal:
            raise RuntimeError('Duplicate SignalP result: ' + fields[0])

        signal[fields[0]] = (
            fields[1],
            fields[3],
            fields[4] if len(fields) > 4 else '.'
        )

expected_sp = set(
    m['signalp_id'] for m in models if m['signalp_id']
)

if set(signal) != expected_sp:
    raise RuntimeError(
        'SignalP ID mismatch. Missing: {} Extra: {}'.format(
            sorted(expected_sp-set(signal)),
            sorted(set(signal)-expected_sp)
        )
    )

# ------------------------------------------------------------
# Pfam domain matches
# ------------------------------------------------------------

expected_pfam = {
    'SVSP': 'PF00089',
    'Ficolin': 'PF00147',
    'PLA2': 'PF00068',
    'CTL': 'PF00059',
    'HYAL': 'PF01630'
}

by_query = {m['query_id']: m for m in models}
domains = defaultdict(list)

with open(domfile) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue

        parts = line.split(None, 22)

        if len(parts) < 22:
            raise RuntimeError('Malformed hmmscan domtblout line')

        fragment = parts[3]
        matched = re.fullmatch(r'(Q\d{5})F(\d{3})', fragment)

        if not matched:
            continue

        model = by_query.get(matched.group(1))

        if not model:
            continue

        if model['protein'] and model['protein'].endswith('*'):
            fragments = model['protein'][:-1].split('*')
        else:
            fragments = (model['protein'] or '').split('*')

        fragment_number = int(matched.group(2))

        if not 1 <= fragment_number <= len(fragments):
            continue

        offset = sum(len(s)+1 for s in fragments[:fragment_number-1])
        a, b = offset + int(parts[17]), offset + int(parts[18])

        accession = parts[1].split('.')[0]
        name = parts[0]
        ie = float(parts[12])

        hmm_coverage = (
            100.0 *
            (int(parts[16])-int(parts[15])+1) /
            int(parts[2])
        )

        first_stop = model['first_internal_stop']

        if first_stop is None or b < first_stop:
            stop_side = 'BEFORE_FIRST_INTERNAL_STOP'
        else:
            stop_side = 'AFTER_FIRST_INTERNAL_STOP'

        domains[model['gene']].append((
            accession, name, ie, hmm_coverage, a, b, stop_side
        ))

with open(domain_file, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')

    w.writerow([
        'Gene', 'Pfam', 'Name', 'Domain_iEvalue',
        'HMM_coverage_pct', 'Protein_start',
        'Protein_end', 'Relative_to_first_internal_stop'
    ])

    for gene in sorted(domains):
        for hit in sorted(domains[gene], key=lambda h: h[2]):
            w.writerow([
                gene, hit[0], hit[1],
                '{:.3g}'.format(hit[2]),
                round(hit[3], 2),
                hit[4], hit[5], hit[6]
            ])

# ------------------------------------------------------------
# Identical protein predictions
# ------------------------------------------------------------

duplicate_groups = defaultdict(list)

for m in models:
    if m['protein_sha256']:
        duplicate_groups[m['protein_sha256']].append(m['gene'])

# ------------------------------------------------------------
# Coding-sequence overlap between gene models
# ------------------------------------------------------------

overlap_info = defaultdict(list)

for i, a in enumerate(models):
    for b in models[i+1:]:
        if a['chrom'] != b['chrom']:
            continue

        bp = sum(
            max(0, min(x2,y2)-max(x1,y1)+1)
            for x1,x2 in a['segments']
            for y1,y2 in b['segments']
        )

        if bp:
            fraction = 100.0 * bp / min(a['cds_bp'], b['cds_bp'])

            overlap_info[a['gene']].append((
                b['gene'], bp, fraction,
                a['strand'] == b['strand']
            ))

            overlap_info[b['gene']].append((
                a['gene'], bp, fraction,
                a['strand'] == b['strand']
            ))

# ------------------------------------------------------------
# Write junction-level evidence
# ------------------------------------------------------------

with open(junction_file, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')

    w.writerow([
        'Gene', 'Chromosome',
        'Left_exon_end', 'Right_exon_start',
        'Supporting_alignments_MAPQ20_anchor8'
    ])

    for m in models:
        for left, right in m['junctions']:
            w.writerow([
                m['gene'], m['chrom'],
                left, right,
                junction_counts[(m['gene'],left,right)]
            ])

# ------------------------------------------------------------
# Build all-gene summary
# ------------------------------------------------------------

headers = [
    'Gene',
    'Family',
    'Chromosome',
    'Strand',
    'Isoforms_in_GTF',
    'CDS_segments',
    'CDS_bp',
    'GTF_CDS_feature_bp',
    'Protein_records_matched',
    'Protein_residues',
    'Stop_class',
    'Internal_stops',
    'Terminal_stop_run',
    'First_internal_stop',
    'CDS_FASTA_records_matched',
    'CDS_FASTA_bp',
    'CDS_frame_mod3',
    'Translation_vs_protein',
    'GTF_phase_check',
    'Venom_RNA_CDS_coverage_pct',
    'Venom_RNA_mean_CDS_depth',
    'RNA_junctions_supported_ge3',
    'RNA_junctions_total',
    'RNA_junction_supporting_alignments',
    'SignalP_prediction',
    'SignalP_SP_probability',
    'SignalP_cleavage',
    'Expected_targeted_Pfam',
    'Pfam_status',
    'Best_Pfam_domain_iEvalue',
    'Best_Pfam_HMM_coverage_pct',
    'Pfam_stop_relation',
    'Identical_protein_other_genes',
    'Same_strand_CDS_overlap_genes',
    'Structural_review_flags',
    'Venom_RNA_note'
]

rows = []

for m in models:
    gene = m['gene']
    family = gene.rsplit('-', 1)[0] if '-' in gene else gene

    pfam = expected_pfam.get(family)
    pfam_status, best_ie, best_cov, best_side = (
        'NOT_ASSESSED', '.', '.', '.'
    )

    if pfam is not None:
        if os.path.getsize(domfile) == 0:
            pfam_status = 'NOT_ASSESSED'
        else:
            hits = [
                x for x in domains[gene]
                if x[0] == pfam
            ]

            hits.sort(key=lambda x: x[2])
            good = [x for x in hits if x[2] <= 1e-3]

            pfam_status = (
                'SUPPORTED' if good
                else 'NO_SIGNIFICANT_HIT'
            )

            if hits:
                best_ie = '{:.3g}'.format(hits[0][2])
                best_cov = round(hits[0][3], 2)
                best_side = hits[0][6]

    c, s = counts[gene]

    cov = round(100.0*c/m['cds_bp'], 2) if m['cds_bp'] else 0
    mean = round(s/m['cds_bp'], 2) if m['cds_bp'] else 0

    jvals = [
        junction_counts[(gene,a,b)]
        for a,b in m['junctions']
    ]

    supported = sum(n >= 3 for n in jvals)

    sp = signal.get(
        m['signalp_id'],
        ('NOT_TESTED', '.', '.')
    )

    twins = [
        g for g in duplicate_groups.get(m['protein_sha256'], [])
        if g != gene
    ] if m['protein_sha256'] else []

    same_overlap = [
        g for g,bp,frac,same in overlap_info[gene]
        if same and bp >= 30
    ]

    significant_overlap = any(
        same and frac >= 50
        for g,bp,frac,same in overlap_info[gene]
    )

    # Structural flags: independent of RNA expression.
    flags = []

    if m['isoforms'] > 1:
        flags.append('MULTIPLE_GTF_ISOFORMS')

    if m['protein_records'] != 1:
        flags.append('PROTEIN_MAPPING_UNCERTAIN')

    if m['cds_records'] != 1:
        flags.append('CDS_FASTA_MAPPING_UNCERTAIN')

    if m['internal_stops'] is not None and m['internal_stops'] > 0:
        flags.append('INTERNAL_STOP')

    if m['cds_frame_mod3'] != 'NA' and m['cds_frame_mod3'] != 0:
        flags.append('CDS_LENGTH_NOT_DIVISIBLE_BY_3')

    if (
        m['cds_fasta_bp'] is not None
        and m['gtf_cds_feature_bp']-m['cds_fasta_bp'] not in (0,3)
    ):
        flags.append('GTF_CDS_LENGTH_DISCREPANCY')

    if m['translation_status'] == 'MISMATCH':
        flags.append('CDS_PROTEIN_TRANSLATION_MISMATCH')

    if m['phase_status'] == 'FAIL':
        flags.append('GTF_PHASE_INCONSISTENT')

    if twins:
        flags.append('IDENTICAL_PREDICTED_PROTEIN')

    if significant_overlap:
        flags.append('LARGE_SAME_STRAND_CDS_OVERLAP')

    if pfam_status == 'NO_SIGNIFICANT_HIT':
        flags.append('EXPECTED_TARGETED_PFAM_NOT_DETECTED')

    # Venom RNA notes remain separate from structural flags.
    rna_notes = []

    if cov < 20:
        rna_notes.append('LOW_CDS_COVERAGE_THIS_VENOM_SAMPLE')

    if jvals and not supported:
        rna_notes.append('NO_ANNOTATED_JUNCTION_WITH_3_READS')

    rows.append([
        gene,
        family,
        m['chrom'],
        m['strand'],
        m['isoforms'],
        len(m['segments']),
        m['cds_bp'],
        m['gtf_cds_feature_bp'],
        m['protein_records'],
        m['protein_length_aa'] if m['protein_length_aa'] is not None else '.',
        m['stop_class'],
        m['internal_stops'] if m['internal_stops'] is not None else '.',
        m['terminal_stop_run'] if m['terminal_stop_run'] is not None else '.',
        m['first_internal_stop'] if m['first_internal_stop'] is not None else '.',
        m['cds_records'],
        m['cds_fasta_bp'] if m['cds_fasta_bp'] is not None else '.',
        m['cds_frame_mod3'],
        m['translation_status'],
        m['phase_status'],
        cov,
        mean,
        supported,
        len(jvals),
        sum(jvals),
        sp[0],
        sp[1],
        sp[2],
        pfam or 'NOT_IN_FIVE_PROFILE_PANEL',
        pfam_status,
        best_ie,
        best_cov,
        best_side,
        ','.join(sorted(twins)) or '.',
        ','.join(sorted(set(same_overlap))) or '.',
        ';'.join(flags) or 'PASS',
        ';'.join(rna_notes) or '.'
    ])

# ------------------------------------------------------------
# Write QC tables
# ------------------------------------------------------------

with open(summary_file, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(headers)
    w.writerows(rows)

with open(flagged_file, 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(headers)
    w.writerows(r for r in rows if r[-2] != 'PASS')

# ------------------------------------------------------------
# Print concise job summary
# ------------------------------------------------------------

print('\n===== ALL-TOXIN AUTOMATED QC =====')

print(
    'Genes:', len(rows),
    '| Structural review flags:',
    sum(r[-2] != 'PASS' for r in rows)
)

print(
    'SignalP SP predictions:',
    sum(r[24] == 'SP' for r in rows)
)

print(
    'Proteins with genuine internal stops:',
    sum(
        isinstance(m['internal_stops'], int)
        and m['internal_stops'] > 0
        for m in models
    )
)

print(
    'Tandem terminal stop sequences:',
    sum(
        m['stop_class'] == 'TANDEM_TERMINAL_STOP'
        for m in models
    )
)

print('Gene\tStructural_review_flags\tVenom_RNA_note')

for r in rows:
    if r[-2] != 'PASS':
        print('{}\t{}\t{}'.format(
            r[0], r[-2], r[-1]
        ))

print('\nMain table:', summary_file)
print('Structural flags:', flagged_file)
print('Splice junctions:', junction_file)
print('Targeted Pfam hits:', domain_file)
PY

echo
echo '===== ALL-TOXIN QC COMPLETED ====='
echo "Main QC: $SUMMARY"
echo "Flagged: $FLAGGED"