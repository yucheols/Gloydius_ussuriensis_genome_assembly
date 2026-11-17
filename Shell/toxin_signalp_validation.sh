#!/bin/bash
#SBATCH --job-name=toxin_signalp_validation
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

# activate conda env
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# ------------------------------------------------------------
# STEP 0: fix shared-library compatibility
# ------------------------------------------------------------

# prioritize Conda's C++ runtime over the older system libraries
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

# restrict CPU thread usage
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"

echo "===== ENVIRONMENT VALIDATION ====="
echo "Conda environment: $CONDA_DEFAULT_ENV"
echo "Conda prefix:      $CONDA_PREFIX"
echo "Python:            $(command -v python)"
echo "SignalP 6:         $(command -v signalp6)"

if [[ ! -s "${CONDA_PREFIX}/lib/libstdc++.so.6" ]]; then
    echo "ERROR: Conda libstdc++.so.6 not found." >&2
    exit 1
fi

# verify that the original import error is resolved
if ! python -c 'from PIL import Image; import matplotlib; import signalp; print("SignalP dependencies: PASS")'; then
    echo "ERROR: SignalP dependencies could not be loaded." >&2
    echo "Check the Conda C++ runtime and GLIBCXX_3.4.29." >&2
    exit 1
fi


# ------------------------------------------------------------
# directories
# ------------------------------------------------------------
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
OUTDIR="${WORKDIR}/protein_validation"

# input
PEP="${WORKDIR}/toxin_annotation_pep.fasta"

# intermediate and final outputs
INPUT="${OUTDIR}/toxin_signalp6_input.fasta"
METADATA="${OUTDIR}/toxin_signalp6_metadata.tsv"
SIGNAL_OUT="${OUTDIR}/signalp6_${SLURM_JOB_ID:-manual}"
SUMMARY="${OUTDIR}/toxin_signalp6_summary.tsv"
TMP_SUMMARY="${SUMMARY}.tmp.${SLURM_JOB_ID:-$$}"

mkdir -p "$OUTDIR"
cd "$WORKDIR"

trap 'rm -f "$TMP_SUMMARY"' EXIT

if [[ ! -s "$PEP" ]]; then
    echo "ERROR: Protein FASTA not found: $PEP" >&2
    exit 1
fi


# -------------------------------------------------------------------------------------
# STEP 1: prepare all 10 protein sequences
#
# SignalP requires an amino acid sequence
# for disrupted proteins, use only the sequence receding the first internal stop codon
#
# original protein sequences are not modified
# -------------------------------------------------------------------------------------

echo
echo "===== STEP 1: PREPARE PROTEINS ====="

python - "$PEP" "$INPUT" "$METADATA" <<'PY'
import sys
import re
import csv

source, input_file, metadata_file = sys.argv[1:]

genes = [
    'SVSP-1',
    'SVSP-16',
    'Ficolin-4',
    'PLA2-4',
    'CTL-10',
    'HYAL-4',
    'CTL-15',
    'CTL-16',
    'CTL-17',
    'CTL-18'
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

with open(input_file, 'w') as fasta, \
     open(metadata_file, 'w', newline='') as meta:

    writer = csv.writer(meta, delimiter='\t')

    writer.writerow([
        'SignalP_ID',
        'Gene',
        'Residues_aa',
        'Internal_stops',
        'First_stop_position',
        'SignalP_input_length'
    ])

    for index, gene in enumerate(genes, 1):

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
                '{}: expected one sequence, found {}'.format(
                    gene, len(matches)
                )
            )

        seq = matches[0][:-1] if matches[0].endswith('*') else matches[0]

        internal_stops = seq.count('*')
        first_stop = seq.find('*')

        if first_stop >= 0:
            signal_seq = seq[:first_stop]
        else:
            signal_seq = seq

        if not signal_seq:
            raise RuntimeError(
                'Empty SignalP sequence: ' + gene
            )

        # Use simple identifiers for SignalP compatibility.
        signal_id = 'TOXIN{:02d}'.format(index)

        fasta.write(
            '>{}\n{}\n'.format(signal_id, signal_seq)
        )

        writer.writerow([
            signal_id,
            gene,
            sum(aa != '*' for aa in seq),
            internal_stops,
            first_stop + 1 if first_stop >= 0 else 'NA',
            len(signal_seq)
        ])

        print(
            '{}\tinput_length={}\tinternal_stops={}'.format(
                gene,
                len(signal_seq),
                internal_stops
            )
        )

print('Prepared 10 protein sequences.')
PY


# ------------------------------------------------------------------------
# STEP 2: run SignalP 6
#
# Eukaryotic proteins; fast prediction mode; no individual sequence plots
# ------------------------------------------------------------------------

echo
echo "===== STEP 2: SIGNALP 6 ====="

signalp6 \
    --fastafile "$INPUT" \
    --organism eukarya \
    --output_dir "$SIGNAL_OUT" \
    --format none \
    --mode fast \
    --write_procs 1 \
    --torch_num_threads "${SLURM_CPUS_PER_TASK:-4}"

RESULTS="${SIGNAL_OUT}/prediction_results.txt"

if [[ ! -s "$RESULTS" ]]; then
    echo "ERROR: SignalP prediction_results.txt missing." >&2
    exit 1
fi

echo "SignalP 6 predictions completed."


# ------------------------------------------------------------------------
# STEP 3: generate final summary
# ------------------------------------------------------------------------

echo
echo "===== STEP 3: GENERATE SUMMARY ====="

python - "$METADATA" "$RESULTS" "$TMP_SUMMARY" <<'PY'
import sys
import csv

metadata_file, signal_file, summary_file = sys.argv[1:]

metadata = {}
gene_order = []

with open(metadata_file) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        signal_id = row['SignalP_ID']
        metadata[signal_id] = row
        gene_order.append(signal_id)

if len(metadata) != 10:
    raise RuntimeError(
        'Expected metadata for 10 proteins.'
    )

predictions = {}

with open(signal_file) as f:
    for line in f:

        if not line.strip() or line.startswith('#'):
            continue

        fields = line.rstrip('\n').split('\t')

        # Fallback for whitespace-separated output.
        if len(fields) < 4:
            fields = line.strip().split(None, 4)

        if len(fields) < 4:
            raise RuntimeError(
                'Unexpected SignalP output row: ' + line.strip()
            )

        signal_id = fields[0].strip()

        if signal_id not in metadata:
            raise RuntimeError(
                'Unexpected SignalP identifier: ' + signal_id
            )

        if signal_id in predictions:
            raise RuntimeError(
                'Duplicate prediction: ' + signal_id
            )

        prediction = fields[1].strip()
        other_probability = float(fields[2])
        sp_probability = float(fields[3])

        cleavage = (
            fields[4].strip()
            if len(fields) > 4 and fields[4].strip()
            else '.'
        )

        predictions[signal_id] = {
            'prediction': prediction,
            'other': other_probability,
            'sp': sp_probability,
            'cleavage': cleavage
        }

if set(predictions) != set(metadata):
    missing = set(metadata) - set(predictions)
    raise RuntimeError(
        'Missing SignalP predictions: ' +
        ','.join(sorted(missing))
    )

rows = []

for signal_id in gene_order:

    m = metadata[signal_id]
    p = predictions[signal_id]

    rows.append([
        m['Gene'],
        m['Residues_aa'],
        m['Internal_stops'],
        m['First_stop_position'],
        m['SignalP_input_length'],
        p['prediction'],
        '{:.6f}'.format(p['other']),
        '{:.6f}'.format(p['sp']),
        p['cleavage']
    ])

with open(summary_file, 'w', newline='') as f:

    writer = csv.writer(f, delimiter='\t')

    writer.writerow([
        'Gene',
        'Residues_aa',
        'Internal_stops',
        'First_stop_position',
        'SignalP_input_length',
        'SignalP_prediction',
        'OTHER_probability',
        'SP_probability',
        'Cleavage_site'
    ])

    writer.writerows(rows)

print(
    'Gene\tInternal_stops\tSignalP_prediction\t'
    'SP_probability\tCleavage_site'
)

for row in rows:
    print('{}\t{}\t{}\t{}\t{}'.format(
        row[0],
        row[2],
        row[5],
        row[7],
        row[8]
    ))

print('\nValidated 10 protein predictions.')
PY


# ------------------------------------------------------------------------------------- 
# STEP 4: validate and publish final output
# -------------------------------------------------------------------------------------

echo
echo "===== STEP 4: OUTPUT VALIDATION ====="

lines=$(wc -l < "$TMP_SUMMARY")

if [[ "$lines" -ne 11 ]]; then
    echo "ERROR: Expected 10 proteins plus header." >&2
    exit 1
fi

mv -f "$TMP_SUMMARY" "$SUMMARY"

echo
echo "===== COMPLETED ====="
echo "Proteins analyzed: 10"
echo "SignalP output:   $SIGNAL_OUT"
echo "Final summary:    $SUMMARY"
echo "Status: PASS"