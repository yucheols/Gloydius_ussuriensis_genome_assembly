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
# activate environment
# ============================================================
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

# ============================================================
# run blastp to identify SVMP and PLA2 flanking gene homologs
# ============================================================

# activate conda environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"
ANNOTATION="$ROOT/annotation/funannotate/update_results"
PROTEINS="$ANNOTATION/Gloydius_ussuriensis_AMNH_21010.proteins.fa"
FLANKING="$SYNTENY/05_annotations/Guss_SVMP_PLA2_flanking_genes.tsv"
BLASTDIR="$SYNTENY/02_blast"

mkdir -p "$BLASTDIR/SVMP" "$BLASTDIR/PLA2"

# 2. verify input files
echo "===== VERIFY INPUTS ====="
for FILE in "$PROTEINS" "$FLANKING"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing input: $FILE"
        exit 1
    fi
    echo "OK: $FILE"
done

for TOOL in python blastp makeblastdb
do
    command -v "$TOOL"
done

# 3. extract flanking gene proteins
echo
echo "===== EXTRACT FLANKING PROTEINS ====="

python - "$PROTEINS" "$FLANKING" "$BLASTDIR" <<'PY'

import sys
import csv
import re
from pathlib import Path
from collections import defaultdict

proteins_file = Path(sys.argv[1])
flanking_file = Path(sys.argv[2])
outdir = Path(sys.argv[3])

# Read protein FASTA
proteins = defaultdict(list)

header = None
seq = []

def save():
    if header is None:
        return

    parts = header.split()
    transcript = parts[0]

    gene = (
        parts[1]
        if len(parts) > 1
        else re.sub(r"-T\d+$", "", transcript)
    )

    proteins[gene].append(
        (transcript, "".join(seq))
    )

with proteins_file.open() as handle:
    for line in handle:
        line = line.strip()

        if line.startswith(">"):
            save()
            header = line[1:]
            seq = []
        elif line:
            seq.append(line)

    save()

# Read flanking-gene inventory
with flanking_file.open() as handle:
    rows = list(
        csv.DictReader(handle, delimiter="\t")
    )

for family in ("SVMP", "PLA2"):

    selected = []

    for side in ("left_flank", "right_flank"):

        subset = [
            row for row in rows
            if row["family"] == family
            and row["position"] == side
        ]

        subset.sort(
            key=lambda x: int(x["distance_bp"])
        )

        selected.extend(subset[:8])

    fasta_out = (
        outdir / family /
        f"Guss_{family}_flanking_proteins.faa"
    )

    manifest_out = (
        outdir / family /
        f"Guss_{family}_query_manifest.tsv"
    )

    manifest = []

    with fasta_out.open("w") as fasta:

        for row in selected:

            gene = row["gene_id"]

            candidates = proteins.get(gene, [])

            if not candidates:
                print(
                    f"WARNING: No protein found for {gene}"
                )
                continue

            # Prefer T1 when available.
            t1 = [
                x for x in candidates
                if x[0] == f"{gene}-T1"
            ]

            if t1:
                transcript, sequence = t1[0]
            else:
                transcript, sequence = max(
                    candidates,
                    key=lambda x: len(x[1])
                )

            if not sequence:
                print(
                    f"WARNING: Empty protein for {gene}"
                )
                continue

            fasta.write(f">{transcript} {gene}\n")

            for i in range(0, len(sequence), 80):
                fasta.write(sequence[i:i+80] + "\n")

            manifest.append({
                "gene_id": gene,
                "transcript_id": transcript,
                "position": row["position"],
                "distance_bp": row["distance_bp"],
                "length_aa": len(sequence)
            })

    with manifest_out.open("w", newline="") as handle:

        fields = [
            "gene_id",
            "transcript_id",
            "position",
            "distance_bp",
            "length_aa"
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(manifest)

    if not manifest:
        raise ValueError(
            f"No proteins recovered for {family}"
        )

    print(f"\n===== {family} =====")

    for row in manifest:
        print(
            f"{row['transcript_id']:<20}"
            f"{row['length_aa']:>5} aa  "
            f"{row['position']}"
        )

    print(f"\nTotal proteins: {len(manifest)}")
    print(f"FASTA: {fasta_out}")
    print(f"Manifest: {manifest_out}")

PY

# 4. prepare comparative protein BLAST databases
echo
echo "===== PREPARE BLAST DATABASES ====="

for SPECIES in \
    Cerastes_gasperettii \
    Crotalus_adamanteus \
    Naja_naja
do
    PROTEOME="$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa"

    DB="$BLASTDIR/${SPECIES}_proteins_db"

    if [[ ! -s "$PROTEOME" ]]; then
        echo "ERROR: Missing proteome: $PROTEOME"
        exit 1
    fi

    echo
    echo "Species: $SPECIES"

    if [[ ! -s "${DB}.pin" && ! -s "${DB}.pdb" ]]; then

        makeblastdb \
            -in "$PROTEOME" \
            -dbtype prot \
            -out "$DB"

    else
        echo "Existing BLAST database found."
    fi
done

# 5. run comparative BLASTP
echo
echo "===== RUN BLASTP ====="

for FAMILY in SVMP PLA2
do
    QUERY="$BLASTDIR/$FAMILY/Guss_${FAMILY}_flanking_proteins.faa"

    for SPECIES in \
        Cerastes_gasperettii \
        Crotalus_adamanteus \
        Naja_naja
    do
        DB="$BLASTDIR/${SPECIES}_proteins_db"

        OUTPUT="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.blastp.tsv"

        echo
        echo "$FAMILY vs $SPECIES"

        blastp \
            -query "$QUERY" \
            -db "$DB" \
            -evalue 1e-5 \
            -max_target_seqs 20 \
            -num_threads 2 \
            -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' \
            -out "$OUTPUT"

        echo "Output: $OUTPUT"
    done
done

# 6. summarize best BLASTP hits
echo
echo "===== BEST BLASTP HITS ====="

python - "$BLASTDIR" <<'PY'
import sys
import csv
from pathlib import Path
from collections import defaultdict

blastdir = Path(sys.argv[1])

species_list = (
    "Cerastes_gasperettii",
    "Crotalus_adamanteus",
    "Naja_naja"
)

for family in ("SVMP", "PLA2"):

    print(f"\n===== {family} =====")

    for species in species_list:

        path = (
            blastdir / family /
            f"Guss_{family}_vs_{species}.blastp.tsv"
        )

        hits = defaultdict(list)

        with path.open() as handle:
            reader = csv.reader(handle, delimiter="\t")

            for row in reader:
                if len(row) < 9:
                    continue

                hits[row[0]].append(row)

        print(f"\n--- {species} ---")

        if not hits:
            print("No hits passing the search threshold.")
            continue

        for query in sorted(hits):

            best = max(
                hits[query],
                key=lambda x: float(x[8])
            )

            print(
                f"{query:<20} "
                f"{best[1]:<28} "
                f"ID={float(best[2]):5.1f}% "
                f"Qcov={float(best[6]):5.1f}% "
                f"E={best[7]}"
            )

print("\n===== COMPLETE =====")
PY