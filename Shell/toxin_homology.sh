#!/bin/bash
#SBATCH --job-name=toxin_homology
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
# comparative SVMP and PLA2 BLASTN and TBLASTN searches
# ============================================================

# 1. environment
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

set -euo pipefail

# 2. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
TOXIN="$ROOT/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
PHYLO="$TOXIN/phylogeny"
CDS="$TOXIN/toxin_annotation_cds.fasta"
INTERVALS="$SYNTENY/03_genomic_intervals"
BLASTDIR="$SYNTENY/02_blast"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# 3. verify and extract query sequences
echo "===== PREPARE QUERY SEQUENCES ====="

python - "$PHYLO" "$CDS" "$BLASTDIR" <<'PY'

import sys
from pathlib import Path

phylo = Path(sys.argv[1])
cds_file = Path(sys.argv[2])
outdir = Path(sys.argv[3])

def read_fasta(path):
    seqs = {}
    name = None
    parts = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if name is not None:
                    if name in seqs:
                        raise ValueError(f"Duplicate ID: {name}")
                    seqs[name] = "".join(parts)

                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)

    if name is not None:
        if name in seqs:
            raise ValueError(f"Duplicate ID: {name}")
        seqs[name] = "".join(parts)

    return seqs

all_cds = read_fasta(cds_file)

for family in ("SVMP", "PLA2"):

    proteins = read_fasta(
        phylo / f"Guss_{family}_core.faa"
    )

    missing = [
        gene for gene in proteins
        if gene not in all_cds
    ]

    if missing:
        raise ValueError(
            f"Missing {family} CDS: {missing}"
        )

    dest = (
        outdir / family /
        f"Guss_{family}_core_cds.fna"
    )

    with dest.open("w") as handle:

        for gene in proteins:
            seq = all_cds[gene]

            if not seq:
                raise ValueError(f"Empty CDS: {gene}")

            handle.write(f">{gene}\n")

            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

    print(
        f"{family}: {len(proteins)} proteins, "
        f"{len(proteins)} matching CDS sequences"
    )

    print(f"CDS output: {dest}")

PY

# 4. build local genomic BLAST databases
echo
echo "===== BUILD GENOMIC BLAST DATABASES ====="

for FAMILY in SVMP PLA2
do
    mkdir -p "$BLASTDIR/$FAMILY/db"

    for SPECIES in \
        Gloydius_ussuriensis \
        Cerastes_gasperettii \
        Crotalus_adamanteus \
        Naja_naja
    do
        GENOME="$INTERVALS/$FAMILY/${SPECIES}_${FAMILY}_anchor_bounded.fa"

        DB="$BLASTDIR/$FAMILY/db/${SPECIES}"

        if [[ ! -s "$GENOME" ]]; then
            echo "ERROR: Missing genomic interval: $GENOME"
            exit 1
        fi

        makeblastdb \
            -in "$GENOME" \
            -dbtype nucl \
            -out "$DB" \
            > /dev/null

        echo "OK: $FAMILY $SPECIES"
    done
done

# 5. run TBLASTN and BLASTN
echo
echo "===== RUN HOMOLOGY SEARCHES ====="

for FAMILY in SVMP PLA2
do
    PROTEINS="$PHYLO/Guss_${FAMILY}_core.faa"
    NUCLEOTIDES="$BLASTDIR/$FAMILY/Guss_${FAMILY}_core_cds.fna"

    for SPECIES in \
        Gloydius_ussuriensis \
        Cerastes_gasperettii \
        Crotalus_adamanteus \
        Naja_naja
    do
        DB="$BLASTDIR/$FAMILY/db/${SPECIES}"

        OUT_TBLASTN="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.tblastn.tsv"

        OUT_BLASTN="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.blastn.tsv"

        echo
        echo "===== $FAMILY vs $SPECIES ====="

        tblastn \
            -query "$PROTEINS" \
            -db "$DB" \
            -evalue 1e-5 \
            -num_threads "$THREADS" \
            -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send sframe evalue bitscore' \
            -out "$OUT_TBLASTN"

        blastn \
            -task blastn \
            -query "$NUCLEOTIDES" \
            -db "$DB" \
            -evalue 1e-5 \
            -num_threads "$THREADS" \
            -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send sstrand evalue bitscore' \
            -out "$OUT_BLASTN"

        echo "TBLASTN HSPs: $(wc -l < "$OUT_TBLASTN")"
        echo "BLASTN HSPs:  $(wc -l < "$OUT_BLASTN")"
    done
done

# 6. verify output files
echo
echo "===== SEARCH OUTPUT SUMMARY ====="

for FAMILY in SVMP PLA2
do
    for SPECIES in \
        Gloydius_ussuriensis \
        Cerastes_gasperettii \
        Crotalus_adamanteus \
        Naja_naja
    do
        for METHOD in tblastn blastn
        do
            FILE="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.${METHOD}.tsv"

            echo "$FAMILY $SPECIES $METHOD: $(wc -l < "$FILE") HSPs"
        done
    done
done

echo
echo "===== COMPLETE ====="