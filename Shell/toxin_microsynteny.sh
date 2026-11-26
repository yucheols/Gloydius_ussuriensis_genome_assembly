
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
# method: blastn local nucleotide alignments
# comparisons: 12 pairwise alignments
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

mkdir -p "$OUTDIR/blastn"

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
        handle,
        delimiter="\t",
        lineterminator="\n"
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

# 3. run blastn genomic alignments
mkdir -p "$OUTDIR/logs"

echo
echo "===== RUN BLASTN GENOMIC ALIGNMENTS ====="

# Retain all BLASTN HSPs; no best-hit or top-hit truncation.
# Query = second species; subject = first species in manifest.
# -subject reproduces the successful diagnostic BLASTN approach.

while IFS=$'\t' read -r FAMILY REF QUERY REF_FASTA QUERY_FASTA
do

    if [[ "$FAMILY" == "family" ]]; then
        continue
    fi

    BASE="${FAMILY}_${REF}__${QUERY}"

    RESULT="$OUTDIR/blastn/${BASE}.tsv"
    LOG="$OUTDIR/logs/${BASE}_blastn.log"

    echo
    echo "===== $FAMILY: $REF vs $QUERY ====="

    blastn \
        -task blastn \
        -query "$QUERY_FASTA" \
        -subject "$REF_FASTA" \
        -dust no \
        -word_size 11 \
        -evalue 1e-5 \
        -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send pident length nident evalue bitscore' \
        -out "$RESULT" \
        2> "$LOG"

    # Stop immediately if BLASTN fails.
    STATUS=$?

    if (( STATUS != 0 )); then
        echo "ERROR: BLASTN failed for $BASE (exit code $STATUS)" >&2
        cat "$LOG" >&2
        exit "$STATUS"
    fi

    HITS=$(wc -l < "$RESULT")
    echo "BLASTN HSPs: $HITS"

done < "$OUTDIR/alignment_pairs.tsv"

# 4. summarize genomic alignment blocks
echo
echo "===== SUMMARIZE ALIGNMENT BLOCKS ====="

python - "$OUTDIR" <<'PYCODE'

import csv
import sys

from pathlib import Path
from collections import Counter

outdir = Path(sys.argv[1])

# ------------------------------------------------------------
# read interval metadata
# ------------------------------------------------------------

with (outdir / "interval_metadata.tsv").open() as handle:

    metadata = {
        (r["species"], r["family"]): r
        for r in csv.DictReader(handle, delimiter="\t")
    }

with (outdir / "alignment_pairs.tsv").open() as handle:
    pairs = list(csv.DictReader(handle, delimiter="\t"))

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
    "aligned_bp",
    "matching_bp",
    "percent_identity",
    "match_fraction",
    "evalue",
    "bitscore",
    "reference_genomic_start",
    "reference_genomic_end",
    "query_genomic_start",
    "query_genomic_end"
]

all_blocks = []
plot_blocks = []

# ------------------------------------------------------------
# process BLASTN files
# ------------------------------------------------------------

for pair in pairs:

    family = pair["family"]
    reference = pair["reference"]
    query = pair["query"]

    ref = metadata[(reference, family)]
    qry = metadata[(query, family)]

    result = (
        outdir / "blastn" /
        f"{family}_{reference}__{query}.tsv"
    )

    counts = Counter()

    with result.open() as handle:

        for line in handle:

            if not line.strip():
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) != 13:
                raise ValueError(
                    f"Invalid BLASTN record in {result}: {line[:100]}"
                )

            (
                qname, sname, qlen, slen,
                qstart, qend, sstart, send,
                pident, alen, nident, evalue, bitscore
            ) = f

            qlen, slen = int(qlen), int(slen)
            qstart, qend = int(qstart), int(qend)
            sstart, send = int(sstart), int(send)
            alen, nident = int(alen), int(nident)
            pident, bitscore = float(pident), float(bitscore)

            if (qname != qry["sequence_id"] or
                    sname != ref["sequence_id"]):
                raise ValueError(
                    f"Sequence ID mismatch in {result}: {qname}, {sname}"
                )

            if (qlen != int(qry["length_bp"]) or
                    slen != int(ref["length_bp"])):
                raise ValueError(
                    f"FASTA length mismatch in {result}"
                )

            # BLASTN coordinates are 1-based inclusive;
            # descending coordinates indicate reverse alignment.
            ref_start, ref_end = sorted((sstart, send))
            qry_start, qry_end = sorted((qstart, qend))

            strand = (
                "+" if (qend - qstart) * (send - sstart) >= 0
                else "-"
            )

            record = {
                "family": family,
                "reference": reference,
                "query": query,
                "reference_chr": ref["chr"],
                "query_chr": qry["chr"],
                "seq_id": f"{reference}_{family}",
                "seq_id2": f"{query}_{family}",
                "start": ref_start,
                "end": ref_end,
                "start2": qry_start,
                "end2": qry_end,
                "strand": strand,
                "aligned_bp": alen,
                "matching_bp": nident,
                "percent_identity": f"{pident:.3f}",
                "match_fraction": f"{pident / 100:.4f}",
                "evalue": evalue,
                "bitscore": f"{bitscore:.1f}",
                "reference_genomic_start":
                    int(ref["genomic_start_1based"]) + ref_start - 1,
                "reference_genomic_end":
                    int(ref["genomic_start_1based"]) + ref_end - 1,
                "query_genomic_start":
                    int(qry["genomic_start_1based"]) + qry_start - 1,
                "query_genomic_end":
                    int(qry["genomic_start_1based"]) + qry_end - 1
            }

            all_blocks.append(record)
            counts["all"] += 1

            # Provisional visualization thresholds
            min_length = 500 if family == "SVMP" else 150

            if alen >= min_length and pident >= 70.0:
                plot_blocks.append(record)
                counts["plot"] += 1

    print(
        f"{family} | {reference} vs {query} | "
        f"All BLASTN HSPs={counts['all']} | "
        f"Filtered BLASTN HSPs={counts['plot']}"
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
    outdir / "filtered_blocks_for_plot.tsv",
    plot_blocks
)

print("\n===== ALIGNMENT SUMMARY =====")
print("Total BLASTN HSPs:", len(all_blocks))
print(
    "Filtered BLASTN HSPs passing provisional thresholds:",
    len(plot_blocks)
)

print("\n===== OUTPUT FILES =====")
print(outdir / "all_local_alignment_blocks.tsv")
print(outdir / "filtered_blocks_for_plot.tsv")

print("\n===== COMPLETE =====")

PYCODE

# 5. identify nonredundant and collinear conserved blocks
echo
echo "===== IDENTIFY CONSERVED COLLINEAR BLOCKS ====="

python - "$OUTDIR" <<'PYBLOCKS'

import csv
import sys

from collections import Counter, defaultdict
from pathlib import Path

outdir = Path(sys.argv[1])

# Provisional criteria. These affect only the additional block tables;
# the unfiltered BLASTN HSPs and Step 4 tables remain unchanged.
FAMILY_SETTINGS = {
    "SVMP": dict(
        seed_bp=1000,
        seed_identity=75.0,
        max_gap=25000,
        min_block_bp=2500,
        single_bp=4000
    ),
    "PLA2": dict(
        seed_bp=300,
        seed_identity=75.0,
        max_gap=5000,
        min_block_bp=600,
        single_bp=1500
    ),
}

BIN_BP = 2000
REDUNDANT_OVERLAP = 0.85
BLOCK_OVERLAP = 0.50


def read_tsv(path):

    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def bins(start, end, width=BIN_BP):

    return range(
        (start - 1) // width,
        (end - 1) // width + 1
    )


def overlap_fraction_shorter(a, b, start_field, end_field):

    a0, a1 = a[start_field], a[end_field]
    b0, b1 = b[start_field], b[end_field]

    overlap = max(
        0,
        min(a1, b1) - max(a0, b0) + 1
    )

    return overlap / min(
        a1 - a0 + 1,
        b1 - b0 + 1
    )


def matching_candidates(item, ref_bins, qry_bins, width=BIN_BP):

    # Existing hits whose intervals intersect the same bins
    # on BOTH genomic axes.

    r = set()
    q = set()

    for k in bins(item["start"], item["end"], width):
        r.update(ref_bins.get(k, ()))

    if not r:
        return ()

    for k in bins(item["start2"], item["end2"], width):
        q.update(qry_bins.get(k, ()))

    return r.intersection(q)


def nonredundant_hsps(hsps):

    # Remove near-identical alternative HSPs at the same
    # pair of loci, preserving distinct paralogous loci.

    ref_bins = defaultdict(set)
    qry_bins = defaultdict(set)

    kept = []

    ordered = sorted(
        hsps,
        key=lambda h: (
            -h["bitscore_num"],
            -h["aligned_bp_num"],
            h["start"],
            h["start2"]
        )
    )

    for h in ordered:

        redundant = False

        for idx in matching_candidates(
            h, ref_bins, qry_bins
        ):

            k = kept[idx]

            if h["strand"] != k["strand"]:
                continue

            length_ratio = (
                min(
                    h["aligned_bp_num"],
                    k["aligned_bp_num"]
                ) /
                max(
                    h["aligned_bp_num"],
                    k["aligned_bp_num"]
                )
            )

            if length_ratio < 0.65:
                continue

            if (
                overlap_fraction_shorter(
                    h, k, "start", "end"
                ) >= REDUNDANT_OVERLAP
                and
                overlap_fraction_shorter(
                    h, k, "start2", "end2"
                ) >= REDUNDANT_OVERLAP
            ):
                redundant = True
                break

        if redundant:
            continue

        idx = len(kept)
        kept.append(h)

        for k in bins(h["start"], h["end"]):
            ref_bins[k].add(idx)

        for k in bins(h["start2"], h["end2"]):
            qry_bins[k].add(idx)

    return sorted(
        kept,
        key=lambda h: (
            h["start"],
            h["start2"],
            h["strand"]
        )
    )


def oriented_query(h):

    # Both axes increase together for a collinear chain.
    # On the minus strand, reverse the query direction.

    if h["strand"] == "+":
        return h["start2"], h["end2"]

    return -h["end2"], -h["start2"]


def make_chains(hsps, settings):

    # Greedy, strand-consistent local chaining.
    # Chains are conserved-block candidates,
    # not proven orthologous segments.

    candidates = []

    for strand in ("+", "-"):

        subset = sorted(
            (
                h for h in hsps
                if h["strand"] == strand
            ),
            key=lambda h: (
                h["start"],
                oriented_query(h)[0]
            )
        )

        chains = []
        active = []

        max_gap = settings["max_gap"]

        for h in subset:

            q0, q1 = oriented_query(h)

            active = [
                c for c in active
                if (
                    c["last_ref_end"] +
                    max_gap + 50
                ) >= h["start"]
            ]

            best = None
            best_cost = None

            for c in active:

                gap_ref = (
                    h["start"] -
                    c["last_ref_end"] - 1
                )

                gap_qry = (
                    q0 -
                    c["last_qend"] - 1
                )

                if not (
                    -50 <= gap_ref <= max_gap
                    and
                    -50 <= gap_qry <= max_gap
                ):
                    continue

                if abs(gap_ref - gap_qry) > max_gap / 2:
                    continue

                cost = (
                    abs(gap_ref - gap_qry)
                    +
                    0.1 * (
                        max(0, gap_ref) +
                        max(0, gap_qry)
                    )
                )

                if (
                    best_cost is None
                    or cost < best_cost
                ):
                    best = c
                    best_cost = cost

            if best is None:

                best = {
                    "members": [],
                    "start": h["start"],
                    "end": h["end"],
                    "start2": h["start2"],
                    "end2": h["end2"],
                    "last_ref_end": h["end"],
                    "last_qend": q1,
                    "strand": strand,
                    "aligned_sum": 0,
                    "matching_sum": 0,
                    "bitscore_sum": 0.0,
                    "longest_hsp": 0
                }

                chains.append(best)
                active.append(best)

            else:

                best["end"] = max(
                    best["end"],
                    h["end"]
                )

                best["start2"] = min(
                    best["start2"],
                    h["start2"]
                )

                best["end2"] = max(
                    best["end2"],
                    h["end2"]
                )

                best["last_ref_end"] = max(
                    best["last_ref_end"],
                    h["end"]
                )

                best["last_qend"] = max(
                    best["last_qend"],
                    q1
                )

            best["members"].append(h)
            best["aligned_sum"] += h["aligned_bp_num"]
            best["matching_sum"] += h["matching_bp_num"]
            best["bitscore_sum"] += h["bitscore_num"]

            best["longest_hsp"] = max(
                best["longest_hsp"],
                h["aligned_bp_num"]
            )

        for c in chains:

            if (
                c["aligned_sum"] >= settings["min_block_bp"]
                and (
                    len(c["members"]) >= 2
                    or
                    c["longest_hsp"] >= settings["single_bp"]
                )
            ):
                candidates.append(c)

    return candidates


def distinct_chains(chains):

    # Remove nested duplicate chains mapping to the
    # same pair of genomic neighborhoods.

    keep = []

    ref_bins = defaultdict(set)
    qry_bins = defaultdict(set)

    ordered = sorted(
        chains,
        key=lambda c: -c["bitscore_sum"]
    )

    for c in ordered:

        is_duplicate = False

        for idx in matching_candidates(
            c, ref_bins, qry_bins
        ):

            other = keep[idx]

            if (
                c["strand"] == other["strand"]
                and
                overlap_fraction_shorter(
                    c, other, "start", "end"
                ) >= 0.85
                and
                overlap_fraction_shorter(
                    c, other, "start2", "end2"
                ) >= 0.85
            ):
                is_duplicate = True
                break

        if is_duplicate:
            continue

        idx = len(keep)
        keep.append(c)

        for k in bins(c["start"], c["end"]):
            ref_bins[k].add(idx)

        for k in bins(c["start2"], c["end2"]):
            qry_bins[k].add(idx)

    return sorted(
        keep,
        key=lambda c: (
            c["start"],
            c["start2"],
            c["strand"]
        )
    )


def flag_overlapping_chains(chains):

    # Flag potentially ambiguous blocks, preserving them
    # for manual examination of toxin paralogs.

    ref_bins = defaultdict(set)
    qry_bins = defaultdict(set)

    for i, c in enumerate(chains):

        c["flags"] = set()

        for k in bins(c["start"], c["end"], 10000):
            ref_bins[k].add(i)

        for k in bins(c["start2"], c["end2"], 10000):
            qry_bins[k].add(i)

    for i, c in enumerate(chains):

        candidate_idxs = set()

        for k in bins(c["start"], c["end"], 10000):
            candidate_idxs.update(ref_bins[k])

        for k in bins(c["start2"], c["end2"], 10000):
            candidate_idxs.update(qry_bins[k])

        for j in candidate_idxs:

            if j <= i:
                continue

            other = chains[j]

            r = overlap_fraction_shorter(
                c, other, "start", "end"
            )

            q = overlap_fraction_shorter(
                c, other, "start2", "end2"
            )

            if (
                r >= BLOCK_OVERLAP
                and q < BLOCK_OVERLAP
            ):

                c["flags"].add("reference_multimapping")
                other["flags"].add("reference_multimapping")

            elif (
                q >= BLOCK_OVERLAP
                and r < BLOCK_OVERLAP
            ):

                c["flags"].add("query_multimapping")
                other["flags"].add("query_multimapping")

            elif (
                r >= BLOCK_OVERLAP
                and q >= BLOCK_OVERLAP
            ):

                c["flags"].add(
                    "competing_overlapping_blocks"
                )
                other["flags"].add(
                    "competing_overlapping_blocks"
                )


pairs = read_tsv(outdir / "alignment_pairs.tsv")
raw_filtered = read_tsv(
    outdir / "filtered_blocks_for_plot.tsv"
)

by_pair = defaultdict(list)

for idx, rec in enumerate(raw_filtered, 1):

    rec = dict(rec)
    rec["hsp_id"] = f"HSP{idx:07d}"

    for name in (
        "start",
        "end",
        "start2",
        "end2"
    ):
        rec[name] = int(rec[name])

    rec["aligned_bp_num"] = int(rec["aligned_bp"])
    rec["matching_bp_num"] = int(rec["matching_bp"])
    rec["bitscore_num"] = float(rec["bitscore"])
    rec["identity_num"] = float(rec["percent_identity"])

    by_pair[
        (
            rec["family"],
            rec["reference"],
            rec["query"]
        )
    ].append(rec)

original_fields = (
    list(raw_filtered[0].keys())
    if raw_filtered else []
)

original_fields = [
    k for k in original_fields
    if k not in (
        "hsp_id",
        "aligned_bp_num",
        "matching_bp_num",
        "bitscore_num",
        "identity_num"
    )
]

hsp_fields = ["hsp_id"] + original_fields

block_fields = [
    "block_id",
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
    "reference_genomic_start",
    "reference_genomic_end",
    "query_genomic_start",
    "query_genomic_end",
    "hsp_count",
    "sum_aligned_bp",
    "longest_hsp_bp",
    "mean_identity_pct",
    "sum_bitscore",
    "reference_span_bp",
    "query_span_bp",
    "mapping_status"
]

nonredundant = []
blocks = []
members = []
summaries = []

for pair in pairs:

    family = pair["family"]
    ref = pair["reference"]
    query = pair["query"]

    settings = FAMILY_SETTINGS[family]

    source = by_pair.get(
        (family, ref, query), []
    )

    kept = nonredundant_hsps(source)

    for h in kept:
        nonredundant.append(
            {k: h[k] for k in hsp_fields}
        )

    seeds = [
        h for h in kept
        if (
            h["aligned_bp_num"] >= settings["seed_bp"]
            and
            h["identity_num"] >= settings["seed_identity"]
        )
    ]

    chains = distinct_chains(
        make_chains(seeds, settings)
    )

    flag_overlapping_chains(chains)

    status_counts = Counter()

    for i, c in enumerate(chains, 1):

        status = (
            ";".join(sorted(c["flags"]))
            if c["flags"]
            else "unique_candidate"
        )

        block_id = (
            f"{family}_{ref}__{query}_B{i:05d}"
        )

        h0 = c["members"][0]

        ref_offset = (
            int(h0["reference_genomic_start"]) -
            h0["start"]
        )

        qry_offset = (
            int(h0["query_genomic_start"]) -
            h0["start2"]
        )

        block = {
            "block_id": block_id,
            "family": family,
            "reference": ref,
            "query": query,
            "reference_chr": h0["reference_chr"],
            "query_chr": h0["query_chr"],
            "seq_id": h0["seq_id"],
            "seq_id2": h0["seq_id2"],
            "start": c["start"],
            "end": c["end"],
            "start2": c["start2"],
            "end2": c["end2"],
            "strand": c["strand"],
            "reference_genomic_start":
                ref_offset + c["start"],
            "reference_genomic_end":
                ref_offset + c["end"],
            "query_genomic_start":
                qry_offset + c["start2"],
            "query_genomic_end":
                qry_offset + c["end2"],
            "hsp_count": len(c["members"]),
            "sum_aligned_bp": c["aligned_sum"],
            "longest_hsp_bp": c["longest_hsp"],
            "mean_identity_pct":
                f"{100 * c['matching_sum'] / c['aligned_sum']:.2f}",
            "sum_bitscore":
                f"{c['bitscore_sum']:.1f}",
            "reference_span_bp":
                c["end"] - c["start"] + 1,
            "query_span_bp":
                c["end2"] - c["start2"] + 1,
            "mapping_status": status
        }

        blocks.append(block)

        if status == "unique_candidate":
            status_counts["unique"] += 1
        else:
            status_counts["review"] += 1

        for h in c["members"]:
            members.append({
                "block_id": block_id,
                "hsp_id": h["hsp_id"]
            })

    summary = {
        "family": family,
        "reference": ref,
        "query": query,
        "filtered_hsps": len(source),
        "redundant_hsps_removed":
            len(source) - len(kept),
        "nonredundant_hsps": len(kept),
        "seed_hsps": len(seeds),
        "collinear_blocks": len(chains),
        "unique_candidate_blocks":
            status_counts["unique"],
        "blocks_for_review":
            status_counts["review"]
    }

    summaries.append(summary)

    print(
        f"{family:5s} "
        f"{ref:25s} vs {query:25s} | "
        f"Filtered={len(source):6d} "
        f"Dedup={len(kept):6d} "
        f"Seeds={len(seeds):6d} "
        f"Blocks={len(chains):5d} "
        f"Unique={status_counts['unique']:4d} "
        f"Review={status_counts['review']:4d}"
    )

write_tsv(
    outdir / "nonredundant_alignment_hsps.tsv",
    nonredundant,
    hsp_fields
)

write_tsv(
    outdir / "conserved_collinear_blocks.tsv",
    blocks,
    block_fields
)

write_tsv(
    outdir / "conserved_blocks_unique_candidates.tsv",
    [
        b for b in blocks
        if b["mapping_status"] == "unique_candidate"
    ],
    block_fields
)

write_tsv(
    outdir / "conserved_blocks_multimapping_review.tsv",
    [
        b for b in blocks
        if b["mapping_status"] != "unique_candidate"
    ],
    block_fields
)

write_tsv(
    outdir / "conserved_block_members.tsv",
    members,
    ["block_id", "hsp_id"]
)

write_tsv(
    outdir / "conserved_block_summary.tsv",
    summaries,
    list(summaries[0]) if summaries else []
)

print("\n===== CONSERVED BLOCK SUMMARY =====")
print(
    "Nonredundant filtered HSPs:",
    len(nonredundant)
)
print(
    "Collinear block candidates:",
    len(blocks)
)
print(
    "Unique-mapping candidates:",
    sum(
        b["mapping_status"] == "unique_candidate"
        for b in blocks
    )
)
print(
    "Potentially duplicated/ambiguous blocks:",
    sum(
        b["mapping_status"] != "unique_candidate"
        for b in blocks
    )
)

print("\n===== COMPLETE =====")

PYBLOCKS