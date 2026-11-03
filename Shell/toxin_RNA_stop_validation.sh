#!/bin/bash
#SBATCH --job-name=toxin_RNA_stop_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# directories
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

# input files
REF="${ROOT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
RNA_BAM="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.sorted.bam"
BED="${WORKDIR}/toxin_internal_stops.bed"

# final outputs
OUT_CODONS="${WORKDIR}/toxin_RNA_stop_codons_SRR35908235.final.tsv"
OUT_BASES="${WORKDIR}/toxin_RNA_stop_bases_SRR35908235.final.tsv"

cd "$WORKDIR"

# validate inputs
for file in \
    "$REF" \
    "${REF}.fai" \
    "$RNA_BAM" \
    "${RNA_BAM}.bai" \
    "$BED"; do

    if [[ ! -s "$file" ]]; then
        echo "ERROR: Missing or empty input: $file" >&2
        exit 1
    fi
done

samtools quickcheck -v "$RNA_BAM"

# temporary outputs
TMP_DIR=$(mktemp -d "${WORKDIR}/.RNA_stop_final.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

echo "===== RNA STOP-CODON VALIDATION ====="

python3 - \
    "$REF" \
    "$RNA_BAM" \
    "$BED" \
    "${TMP_DIR}/codons.tsv" \
    "${TMP_DIR}/bases.tsv" <<'PY'

import sys
import csv
import re
import subprocess
from bisect import bisect_left
from collections import Counter, defaultdict

REF, BAM, BED, OUT_CODONS, OUT_BASES = sys.argv[1:]

MIN_BQ = 20
STOPS = {"TAA", "TAG", "TGA"}
DNA = set("ACGT")

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")

def revcomp(seq):
    return seq.translate(
        str.maketrans("ACGT", "TGCA")
    )[::-1]

def select_observation(store, qname, allele, quality):
    """
    Keep one observation per RNA fragment (QNAME).

    If both mates cover a position, retain the
    higher-quality observation. Equal-quality
    disagreements are marked ambiguous.
    """
    old = store.get(qname)

    if old is None or quality > old[1]:
        store[qname] = (allele, quality)

    elif quality == old[1] and allele != old[0]:
        store[qname] = ("N", quality)


# ------------------------------------------------------------
# read BED and reference stop codons
# ------------------------------------------------------------

targets = []

with open(BED) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        if len(f) != 6:
            raise ValueError("Expected six-column BED")

        chrom, start, end, name, score, strand = f

        start = int(start)
        end = int(end)

        if end - start != 3 or strand not in ("+", "-"):
            raise ValueError(f"Invalid stop interval: {name}")

        region = f"{chrom}:{start+1}-{end}"

        result = subprocess.check_output(
            ["samtools", "faidx", REF, region],
            text=True
        )

        genomic_codon = "".join(
            result.splitlines()[1:]
        ).upper()

        if len(genomic_codon) != 3:
            raise ValueError(f"Reference retrieval failed: {name}")

        coding_codon = (
            revcomp(genomic_codon)
            if strand == "-"
            else genomic_codon
        )

        if coding_codon not in STOPS:
            raise ValueError(
                f"{name}: expected stop codon, "
                f"observed {coding_codon}"
            )

        targets.append({
            "name": name,
            "chrom": chrom,
            "start": start + 1,
            "end": end,
            "strand": strand,
            "genomic_ref": genomic_codon,
            "coding_ref": coding_codon
        })

if len(targets) != 14:
    raise ValueError(
        f"Expected 14 stop codons, found {len(targets)}"
    )

# each target has three nucleotide positions.
position_index = {}
positions_by_chrom = defaultdict(list)

for i, target in enumerate(targets):
    for j in range(3):
        pos = target["start"] + j
        key = (target["chrom"], pos)

        if key in position_index:
            raise ValueError(f"Overlapping BED position: {key}")

        position_index[key] = (i, j)
        positions_by_chrom[target["chrom"]].append(pos)

for chrom in positions_by_chrom:
    positions_by_chrom[chrom].sort()


# ---------------------------------------------------------------------------------------------
# accumulate high-quality RNA nucleotide observations
# we are not using mpileup; therefore reference skips (N) are never counted as actual RNA bases
# ---------------------------------------------------------------------------------------------

base_calls = [
    [dict(), dict(), dict()]
    for _ in targets
]

codon_calls = [dict() for _ in targets]
skip_fragments = [set() for _ in targets]

cmd = [
    "samtools", "view",
    "-M",
    "-L", BED,
    "-F", "0xF04",
    "-q", "20",
    BAM
]

print("Reading indexed RNA alignments...", flush=True)

proc = subprocess.Popen(
    cmd,
    stdout=subprocess.PIPE,
    text=True
)

alignments_processed = 0

for line in proc.stdout:

    f = line.rstrip("\n").split("\t")

    if len(f) < 11:
        continue

    qname = f[0]
    chrom = f[2]
    refpos = int(f[3])
    cigar = f[5]
    sequence = f[9].upper()
    qualities = f[10]

    if chrom not in positions_by_chrom:
        continue

    if sequence == "*" or qualities == "*":
        continue

    alignments_processed += 1

    querypos = 0
    positions = positions_by_chrom[chrom]

    # calls obtained from this individual alignment
    read_calls = defaultdict(dict)

    for length_str, op in CIGAR_RE.findall(cigar):

        length = int(length_str)

        if op in ("M", "=", "X"):

            block_end = refpos + length

            k = bisect_left(positions, refpos)

            while k < len(positions) and positions[k] < block_end:

                genomic_pos = positions[k]
                query_index = querypos + genomic_pos - refpos

                if query_index < len(sequence):
                    base = sequence[query_index]
                    bq = ord(qualities[query_index]) - 33

                    if base in DNA and bq >= MIN_BQ:
                        idx, j = position_index[
                            (chrom, genomic_pos)
                        ]
                        read_calls[idx][j] = (base, bq)

                k += 1

            refpos += length
            querypos += length

        elif op == "N":

            block_end = refpos + length

            k = bisect_left(positions, refpos)

            while k < len(positions) and positions[k] < block_end:
                idx, j = position_index[
                    (chrom, positions[k])
                ]
                skip_fragments[idx].add(qname)
                k += 1

            refpos += length

        elif op == "D":
            refpos += length

        elif op in ("I", "S"):
            querypos += length

        # H and P consume neither reference nor query.

    for idx, calls in read_calls.items():

        # Individual nucleotide observations
        for j, (base, bq) in calls.items():
            select_observation(
                base_calls[idx][j],
                qname,
                base,
                bq
            )

        # complete three-base codon on the same read
        if len(calls) == 3:

            genomic_codon = "".join(
                calls[j][0] for j in range(3)
            )

            coding_codon = (
                revcomp(genomic_codon)
                if targets[idx]["strand"] == "-"
                else genomic_codon
            )

            total_quality = sum(
                calls[j][1] for j in range(3)
            )

            select_observation(
                codon_calls[idx],
                qname,
                coding_codon,
                total_quality
            )

if proc.wait() != 0:
    raise RuntimeError("samtools view failed")

print(
    f"Alignments examined: {alignments_processed:,}",
    flush=True
)


# ------------------------------------------------------------
# write final codon-level report
# ------------------------------------------------------------

with open(OUT_CODONS, "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Stop",
        "Chromosome",
        "Strand",
        "Genomic_positions",
        "Reference_coding_codon",
        "RNA_fragments_full_codon",
        "Reference_stop_fragments",
        "Other_stop_fragments",
        "Nonstop_fragments",
        "Reference_stop_percent",
        "Alternative_codons",
        "Reference_skip_fragments"
    ])

    for i, target in enumerate(targets):

        counts = Counter(
            codon
            for codon, quality in codon_calls[i].values()
            if len(codon) == 3
            and set(codon).issubset(DNA)
        )

        ref_stop = target["coding_ref"]

        total = sum(counts.values())
        ref_count = counts[ref_stop]

        other_stop = sum(
            n for codon, n in counts.items()
            if codon in STOPS and codon != ref_stop
        )

        nonstop = sum(
            n for codon, n in counts.items()
            if codon not in STOPS
        )

        alt_codons = ",".join(
            f"{codon}:{n}"
            for codon, n in counts.most_common()
            if codon != ref_stop
        ) or "."

        percent = (
            f"{100 * ref_count / total:.2f}"
            if total else "NA"
        )

        writer.writerow([
            target["name"],
            target["chrom"],
            target["strand"],
            f'{target["start"]}-{target["end"]}',
            ref_stop,
            total,
            ref_count,
            other_stop,
            nonstop,
            percent,
            alt_codons,
            len(skip_fragments[i])
        ])


# ------------------------------------------------------------
# write final nucleotide-level report
# ------------------------------------------------------------

with open(OUT_BASES, "w", newline="") as out:

    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Stop",
        "Chromosome",
        "Position",
        "Reference_base",
        "RNA_base_depth",
        "A",
        "C",
        "G",
        "T",
        "Reference_count",
        "Alternative_count",
        "Reference_percent"
    ])

    for i, target in enumerate(targets):

        for j in range(3):

            counts = Counter(
                base
                for base, quality in base_calls[i][j].values()
                if base in DNA
            )

            ref_base = target["genomic_ref"][j]
            total = sum(counts.values())

            ref_count = counts[ref_base]
            alt_count = total - ref_count

            percent = (
                f"{100 * ref_count / total:.2f}"
                if total else "NA"
            )

            writer.writerow([
                target["name"],
                target["chrom"],
                target["start"] + j,
                ref_base,
                total,
                counts["A"],
                counts["C"],
                counts["G"],
                counts["T"],
                ref_count,
                alt_count,
                percent
            ])

print("Stop codons analyzed:", len(targets))
print("Nucleotide positions:", len(targets) * 3)
print("Codon-level report:", OUT_CODONS)
print("Nucleotide-level report:", OUT_BASES)
PY

# create outputs only after successful analysis
[[ $(wc -l < "${TMP_DIR}/codons.tsv") -eq 15 ]]
[[ $(wc -l < "${TMP_DIR}/bases.tsv") -eq 43 ]]

mv -f "${TMP_DIR}/codons.tsv" "$OUT_CODONS"
mv -f "${TMP_DIR}/bases.tsv" "$OUT_BASES"

echo
echo "===== COMPLETED ====="
echo "Codon-level results: $OUT_CODONS"
echo "Base-level results:  $OUT_BASES"
echo "Status: PASS"