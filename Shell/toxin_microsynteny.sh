#!/bin/bash
#SBATCH --job-name=toxin_microsynteny
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================================
# local genomic alignments
#
# gene families: SVMP and PLA2
# species: G. ussuriensis, Cerastes gasperettii, Crotalus adamanteus, N. naja
# Method: minimap2 asm20
# Comparisons: 12 pairwise alignments
# ============================================================================

# activate conda environment
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate


# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
MANIFEST="$SYNTENY/01_reference_inventory/anchor_bounded_regions.tsv"
INTERVALS="$SYNTENY/03_genomic_intervals"
OUTDIR="$SYNTENY/06_local_alignments"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$OUTDIR/paf"

# 2. validate sequences and prepare pairwise comparisons
echo "===== VALIDATE ANCHOR-BOUNDED FASTAS ====="

python - "$MANIFEST" "$INTERVALS" "$OUTDIR" <<'PY'

import csv
import sys

from pathlib import Path
from itertools import combinations

manifest = Path(sys.argv[1])
intervals = Path(sys.argv[2])
outdir = Path(sys.argv[3])

species = [
    "Gloydius_ussuriensis",
    "Cerastes_gasperettii",
    "Crotalus_adamanteus",
    "Naja_naja"
]

families = ["SVMP", "PLA2"]

def fasta_info(path):

    name = None
    length = 0

    with path.open() as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if name is not None:
                    raise ValueError(
                        f"Multiple FASTA sequences: {path}"
                    )

                name = line[1:].split()[0]

            else:

                if name is None:
                    raise ValueError(
                        f"Missing FASTA header: {path}"
                    )

                length += len(line)

    if name is None:
        raise ValueError(f"Empty FASTA: {path}")

    return name, length


with manifest.open() as handle:

    indexed = {
        (r["species"], r["family"]): r
        for r in csv.DictReader(
            handle, delimiter="\t"
        )
    }

metadata = []
lookups = {}

for family in families:

    for sp in species:

        key = (sp, family)

        if key not in indexed:
            raise ValueError(
                f"Missing manifest row: {key}"
            )

        r = indexed[key]

        fasta = (
            intervals /
            family /
            f"{sp}_{family}_anchor_bounded.fa"
        )

        if not fasta.is_file():
            raise FileNotFoundError(fasta)

        header, bp = fasta_info(fasta)

        expected = (
            int(r["end"]) -
            int(r["start"]) + 1
        )

        if bp != expected:
            raise ValueError(
                f"Length mismatch: {sp} {family}: "
                f"FASTA={bp}, manifest={expected}"
            )

        item = {
            "species": sp,
            "family": family,
            "chr": r["chromosome"],
            "genomic_start_1based": r["start"],
            "genomic_end_1based": r["end"],
            "sequence_id": header,
            "length_bp": bp,
            "fasta": str(fasta)
        }

        metadata.append(item)
        lookups[key] = item

        print(
            f"{family:5s} "
            f"{sp:25s} "
            f"{bp:>10,} bp"
        )

# Export interval information

with (
    outdir / "interval_metadata.tsv"
).open("w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=list(metadata[0]),
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(metadata)

# Export alignment pair manifest

with (
    outdir / "alignment_pairs.tsv"
).open("w", newline="") as handle:

    writer = csv.writer(
        handle, delimiter="\t"
    )

    writer.writerow([
        "family",
        "reference",
        "query",
        "reference_fasta",
        "query_fasta"
    ])

    for family in families:

        for a, b in combinations(species, 2):

            writer.writerow([
                family,
                a,
                b,
                lookups[(a, family)]["fasta"],
                lookups[(b, family)]["fasta"]
            ])

PY


# 3. run minimap2 genomic alignments
echo
echo "===== RUN MINIMAP2 (ASM20) ====="

while IFS=$'\t' read -r FAMILY REF QUERY REF_FASTA QUERY_FASTA
do

    if [[ "$FAMILY" == "family" ]]; then
        continue
    fi

    PAF="$OUTDIR/paf/${FAMILY}_${REF}__${QUERY}.paf"

    echo
    echo "$FAMILY: $REF vs $QUERY"

    minimap2 \
        -x asm20 \
        -c \
        --cs=short \
        -N 50 \
        -t "$THREADS" \
        "$REF_FASTA" \
        "$QUERY_FASTA" \
        > "$PAF"

    echo "PAF records: $(wc -l < "$PAF")"

done < "$OUTDIR/alignment_pairs.tsv"

# 4. summarize genomic alignment blocks
echo
echo "===== SUMMARIZE ALIGNMENT BLOCKS ====="

python - "$OUTDIR" <<'PY'

import csv
import sys

from pathlib import Path
from collections import Counter

outdir = Path(sys.argv[1])

# ------------------------------------------------------------
# read interval metadata
# ------------------------------------------------------------

with (
    outdir / "interval_metadata.tsv"
).open() as handle:

    metadata = {
        (r["species"], r["family"]): r
        for r in csv.DictReader(
            handle, delimiter="\t"
        )
    }

with (
    outdir / "alignment_pairs.tsv"
).open() as handle:

    pairs = list(
        csv.DictReader(
            handle, delimiter="\t"
        )
    )

# ------------------------------------------------------------
# initialize output tables
# ------------------------------------------------------------

fields = [
    "family",
    "reference",
    "query",
    "reference_chr",
    "query_chr",
    "seq_id",
    "seq_id2",
    "start",
    "end",
    "start2",
    "end2",
    "strand",
    "alignment_type",
    "mapq",
    "aligned_bp",
    "matching_bp",
    "match_fraction",
    "reference_genomic_start",
    "reference_genomic_end",
    "query_genomic_start",
    "query_genomic_end"
]

all_blocks = []
plot_blocks = []

# ------------------------------------------------------------
# process PAF files
# ------------------------------------------------------------

for pair in pairs:

    family = pair["family"]
    reference = pair["reference"]
    query = pair["query"]

    ref = metadata[(reference, family)]
    qry = metadata[(query, family)]

    paf = (
        outdir /
        "paf" /
        f"{family}_{reference}__{query}.paf"
    )

    counts = Counter()

    with paf.open() as handle:

        for line in handle:

            if not line.strip():
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) < 12:
                raise ValueError(
                    f"Invalid PAF record: {paf}"
                )

            (
                qname, qlen, qs, qe,
                strand,
                tname, tlen, ts, te,
                nmatch, alen, mapq
            ) = f[:12]

            (
                qlen, qs, qe,
                tlen, ts, te,
                nmatch, alen, mapq
            ) = map(
                int,
                (
                    qlen, qs, qe,
                    tlen, ts, te,
                    nmatch, alen, mapq
                )
            )

            # Verify sequence identity and length

            if (
                qname != qry["sequence_id"]
                or tname != ref["sequence_id"]
            ):
                raise ValueError(
                    f"Sequence mismatch: {paf}"
                )

            if (
                qlen != int(qry["length_bp"])
                or tlen != int(ref["length_bp"])
            ):
                raise ValueError(
                    f"Sequence length mismatch: {paf}"
                )

            # Retrieve PAF alignment type

            tags = {
                x[:2]: x[5:]
                for x in f[12:]
                if len(x) > 5
                and x[2] == ":"
                and x[4] == ":"
            }

            alignment_type = tags.get(
                "tp", "unknown"
            )

            match_fraction = (
                nmatch / alen
                if alen else 0
            )

            # Convert to 1-based inclusive coordinates

            record = {
                "family": family,
                "reference": reference,
                "query": query,
                "reference_chr": ref["chr"],
                "query_chr": qry["chr"],
                "seq_id": f"{reference}_{family}",
                "seq_id2": f"{query}_{family}",
                "start": ts + 1,
                "end": te,
                "start2": qs + 1,
                "end2": qe,
                "strand": strand,
                "alignment_type": alignment_type,
                "mapq": mapq,
                "aligned_bp": alen,
                "matching_bp": nmatch,
                "match_fraction": f"{match_fraction:.4f}",
                "reference_genomic_start":
                    int(ref["genomic_start_1based"]) + ts,
                "reference_genomic_end":
                    int(ref["genomic_start_1based"]) + te - 1,
                "query_genomic_start":
                    int(qry["genomic_start_1based"]) + qs,
                "query_genomic_end":
                    int(qry["genomic_start_1based"]) + qe - 1
            }

            all_blocks.append(record)
            counts["all"] += 1

            if alignment_type == "P":
                counts["primary"] += 1

            # Provisional visualization thresholds

            min_length = (
                500 if family == "SVMP"
                else 150
            )

            if (
                alignment_type == "P"
                and alen >= min_length
                and match_fraction >= 0.70
            ):

                plot_blocks.append(record)
                counts["plot"] += 1

    print(
        f"{family} | "
        f"{reference} vs {query} | "
        f"All={counts['all']} | "
        f"Primary={counts['primary']} | "
        f"Plot={counts['plot']}"
    )

# ------------------------------------------------------------
# export results
# ------------------------------------------------------------

def write_table(path, rows):

    with path.open("w", newline="") as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(rows)


write_table(
    outdir / "all_local_alignment_blocks.tsv",
    all_blocks
)

write_table(
    outdir / "primary_blocks_for_plot.tsv",
    plot_blocks
)

print("\n===== ALIGNMENT SUMMARY =====")

print(
    "Total alignment records:",
    len(all_blocks)
)

print(
    "Primary blocks passing provisional thresholds:",
    len(plot_blocks)
)

print("\n===== OUTPUT FILES =====")

print(
    outdir / "all_local_alignment_blocks.tsv"
)

print(
    outdir / "primary_blocks_for_plot.tsv"
)

print("\n===== COMPLETE =====")
PY