#!/usr/bin/env python3

from itertools import product

from Bio import SeqIO
from collections import defaultdict
from pathlib import Path
import csv
import re
import sys


# ============================================================
# set paths
# ============================================================

CONSREF = Path(
    "/home/yshin/mendel-nas1/snake_genome_ass/"
    "G_ussuriensis_Chromo/mitogenome_curation/conspecific_ref"
)

ALIGNMENT = Path("/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitogenome_curation/mitogenomes_merged.mafft.fasta")

GB_FILES = {
    "AMNH_21010_Ref": Path(
        "/home/yshin/mendel-nas1/ussuri_popgen/"
        "mitogenome/AMNH_21010_Ref.gb"
    ),
    "OR680782": CONSREF / "OR680782.gb",
    "KP262412": CONSREF / "KP262412.gb",
    "NC_026553": CONSREF / "NC_026553.gb",
}

# output directory for boundary-audit results
OUTDIR = Path("/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitogenome_curation/boundary_mismatch")
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_LONG = OUTDIR / "feature_boundaries_long.tsv"
OUT_WIDE = OUTDIR / "feature_boundaries_wide.tsv"
OUT_OVERLAPS = OUTDIR / "adjacent_feature_gaps_overlaps.tsv"


# ============================================================
# define helpers
# ============================================================

def strip_version(name):
    """
    OR680782.1 -> OR680782
    NC_026553.1 -> NC_026553
    """
    return re.sub(r"\.\d+$", "", name)


def clean_name(x):
    if x is None:
        return ""
    return str(x).strip()


def normalize_gene_name(name):
    """
    Normalize common mitochondrial gene-name variants.
    """
    if not name:
        return ""

    x = name.upper().replace("-", "").replace("_", "").replace(" ", "")

    aliases = {
        "COI": "COX1",
        "CO1": "COX1",
        "COXI": "COX1",

        "COII": "COX2",
        "CO2": "COX2",
        "COXII": "COX2",

        "COIII": "COX3",
        "CO3": "COX3",
        "COXIII": "COX3",

        "CYTB": "CYTB",
        "COB": "CYTB",

        "NAD1": "ND1",
        "NAD2": "ND2",
        "NAD3": "ND3",
        "NAD4": "ND4",
        "NAD4L": "ND4L",
        "NAD5": "ND5",
        "NAD6": "ND6",

        "ATPASE6": "ATP6",
        "ATPASE8": "ATP8",

        "RRNS": "12S",
        "12S": "12S",
        "12SRRNA": "12S",

        "RRNL": "16S",
        "16S": "16S",
        "16SRRNA": "16S",
    }

    return aliases.get(x, name)


def feature_base_label(feature):
    """
    Produce a biologically useful label before duplicate-numbering.
    """

    ftype = feature.type

    gene = clean_name(
        feature.qualifiers.get("gene", [""])[0]
    )

    product = clean_name(
        feature.qualifiers.get("product", [""])[0]
    )

    # ------------------------
    # protein-coding genes
    # ------------------------
    if ftype == "CDS":
        label = gene if gene else product
        return normalize_gene_name(label)

    # ------------------------
    # rRNAs
    # ------------------------
    if ftype == "rRNA":
        candidate = gene if gene else product

        x = candidate.lower()

        if (
            "rrns" in x
            or "12s" in x
            or "small subunit" in x
        ):
            return "12S"

        if (
            "rrnl" in x
            or "16s" in x
            or "large subunit" in x
        ):
            return "16S"

        return normalize_gene_name(candidate)

    # ------------------------
    # tRNAs
    # ------------------------
    if ftype == "tRNA":

        if gene:
            label = gene.replace(" ", "_")
        elif product:
            label = product.replace(" ", "_")
        else:
            return "tRNA"

        label = re.sub(r"([_-]?)[12]$", "", label)

        return label

    # ------------------------
    # control regions
    # ------------------------
    if ftype == "D-loop":
        return "D-loop"

    if ftype == "rep_origin":
        return "rep_origin"

    return ftype


def get_qualifier(feature, key):
    values = feature.qualifiers.get(key, [])
    if not values:
        return ""
    return ",".join(str(x) for x in values)


def map_sequence_to_alignment(aligned_seq):
    """
    Return:
        ref_position (1-based ungapped)
            ->
        alignment_column (1-based)

    Example:
       A C - T
       1 2   3

    mapping:
       1 -> 1
       2 -> 2
       3 -> 4
    """

    mapping = {}

    ungapped_pos = 0

    for aln_col, base in enumerate(aligned_seq, start=1):

        if base not in "-.":
            ungapped_pos += 1
            mapping[ungapped_pos] = aln_col

    return mapping


# ============================================================
# input validation
# ============================================================

if not ALIGNMENT.exists():
    sys.exit(
        f"ERROR: alignment not found:\n{ALIGNMENT}"
    )

for name, gb in GB_FILES.items():
    if not gb.exists():
        sys.exit(
            f"ERROR: GenBank file missing for {name}:\n{gb}"
        )


# ============================================================
# read mafft alignment
# ============================================================

alignment_records = list(
    SeqIO.parse(ALIGNMENT, "fasta")
)

if not alignment_records:
    sys.exit("ERROR: alignment contains no sequences.")


alignment_length = len(alignment_records[0].seq)

for rec in alignment_records:
    if len(rec.seq) != alignment_length:
        sys.exit(
            "ERROR: alignment sequences have different lengths."
        )


# build aliases such that OR680782.1 matches OR680782
alignment_by_name = {}

for rec in alignment_records:
    raw = rec.id
    base = strip_version(raw)

    alignment_by_name[raw] = rec
    alignment_by_name[base] = rec


print(f"Alignment: {ALIGNMENT}")
print(f"Alignment length: {alignment_length}")
print(f"Alignment sequences: {len(alignment_records)}")
print()


# ============================================================
# read GenBank records and validate against MAFFT rows
# ============================================================

gb_records = {}
position_maps = {}

for accession, gb_path in GB_FILES.items():

    gb_record = SeqIO.read(gb_path, "genbank")
    gb_records[accession] = gb_record

    # Locate corresponding sequence in alignment
    aln_record = None

    for key in [
        accession,
        strip_version(accession),
        gb_record.id,
        strip_version(gb_record.id),
    ]:
        if key in alignment_by_name:
            aln_record = alignment_by_name[key]
            break

    if aln_record is None:
        available = sorted(
            set(rec.id for rec in alignment_records)
        )

        sys.exit(
            f"\nERROR: Could not find {accession} in alignment.\n"
            f"Available IDs:\n" +
            "\n".join(available)
        )

    ungapped_alignment = (
        str(aln_record.seq)
        .replace("-", "")
        .replace(".", "")
        .upper()
    )

    gb_sequence = str(gb_record.seq).upper()

    print(
        f"{accession}: "
        f"GB={len(gb_sequence)} bp; "
        f"alignment ungapped={len(ungapped_alignment)} bp"
    )

    if ungapped_alignment != gb_sequence:

        # Find first mismatch if lengths are equal
        msg = (
            f"\nERROR: Alignment row for {accession} "
            f"does not exactly match its GenBank sequence."
        )

        if len(ungapped_alignment) == len(gb_sequence):
            for i, (a, b) in enumerate(
                zip(ungapped_alignment, gb_sequence),
                start=1
            ):
                if a != b:
                    msg += (
                        f"\nFirst mismatch at nucleotide {i}: "
                        f"alignment={a}, GenBank={b}"
                    )
                    break
        else:
            msg += (
                f"\nLengths differ: "
                f"{len(ungapped_alignment)} vs "
                f"{len(gb_sequence)}"
            )

        sys.exit(msg)

    print("  sequence match: PASS")

    position_maps[accession] = map_sequence_to_alignment(
        str(aln_record.seq)
    )

print()
print("All four MAFFT rows exactly match their GenBank sequences.")
print()


# ============================================================
# extract functional annotation features
# ============================================================

KEEP_TYPES = {
    "CDS",
    "rRNA",
    "tRNA",
    "D-loop",
    "rep_origin",
}


all_features = []

for accession, record in gb_records.items():

    accession_features = []

    for feature in record.features:

        if feature.type not in KEEP_TYPES:
            continue

        raw_start = int(feature.location.start) + 1
        raw_end = int(feature.location.end)

        strand_value = feature.location.strand

        if strand_value == 1:
            strand = "+"
        elif strand_value == -1:
            strand = "-"
        else:
            strand = "."

        # Biological 5' and 3' boundaries
        if strand == "-":
            bio_start = raw_end
            bio_end = raw_start
        else:
            bio_start = raw_start
            bio_end = raw_end

        mapping = position_maps[accession]

        try:
            alignment_raw_start = mapping[raw_start]
            alignment_raw_end = mapping[raw_end]

            alignment_5prime = mapping[bio_start]
            alignment_3prime = mapping[bio_end]

        except KeyError as e:
            sys.exit(
                f"ERROR mapping {accession} "
                f"{feature.type} position {e} "
                f"onto alignment."
            )

        gene = get_qualifier(feature, "gene")
        product = get_qualifier(feature, "product")

        item = {
            "accession": accession,
            "feature_type": feature.type,
            "base_label": feature_base_label(feature),

            "gene": gene,
            "product": product,

            "raw_start": raw_start,
            "raw_end": raw_end,
            "strand": strand,
            "feature_length": len(feature.location),

            "bio_5prime": bio_start,
            "bio_3prime": bio_end,

            "alignment_start": alignment_raw_start,
            "alignment_end": alignment_raw_end,

            "alignment_5prime": alignment_5prime,
            "alignment_3prime": alignment_3prime,

            "codon_start": get_qualifier(
                feature, "codon_start"
            ),

            "transl_table": get_qualifier(
                feature, "transl_table"
            ),
        }

        accession_features.append(item)

    # --------------------------------------------------------
    # number duplicate labels by genomic order
    #
    # e.g.
    # tRNA-Leu_1
    # tRNA-Leu_2
    # D-loop_1
    # D-loop_2
    # --------------------------------------------------------

    grouped = defaultdict(list)

    for item in accession_features:
        grouped[item["base_label"]].append(item)

    for label, entries in grouped.items():

        entries.sort(key=lambda x: x["raw_start"])

        if len(entries) == 1:
            entries[0]["feature_label"] = label

        else:
            for i, item in enumerate(entries, start=1):
                item["feature_label"] = f"{label}_{i}"

    all_features.extend(accession_features)


# ============================================================
# write long-format feature table
# ============================================================

long_columns = [
    "accession",
    "feature_label",
    "feature_type",
    "gene",
    "product",
    "raw_start",
    "raw_end",
    "strand",
    "feature_length",
    "bio_5prime",
    "bio_3prime",
    "alignment_start",
    "alignment_end",
    "alignment_5prime",
    "alignment_3prime",
    "codon_start",
    "transl_table",
]


all_features.sort(
    key=lambda x: (
        x["accession"],
        x["raw_start"],
        x["raw_end"],
    )
)


with open(OUT_LONG, "w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=long_columns,
        delimiter="\t",
    )

    writer.writeheader()

    for item in all_features:
        writer.writerow(
            {k: item.get(k, "") for k in long_columns}
        )


# ============================================================
# write wide comparison table
# ============================================================

ACCESSION_ORDER = [
    "AMNH_21010_Ref",
    "OR680782",
    "KP262412",
    "NC_026553",
]


# feature_label -> accession -> item
feature_matrix = defaultdict(dict)

for item in all_features:
    feature_matrix[item["feature_label"]][
        item["accession"]
    ] = item


wide_columns = [
    "feature_label",
    "feature_type",
]

for acc in ACCESSION_ORDER:
    wide_columns.extend([
        f"{acc}_raw_start",
        f"{acc}_raw_end",
        f"{acc}_strand",
        f"{acc}_length",
        f"{acc}_aln_start",
        f"{acc}_aln_end",
        f"{acc}_aln_5prime",
        f"{acc}_aln_3prime",
    ])

wide_columns.extend([
    "all_alignment_starts_same",
    "all_alignment_ends_same",
])


# sort features according to AMNH genomic order where possible
def feature_sort_key(label):

    amnh = feature_matrix[label].get("AMNH_21010_Ref")

    if amnh:
        return (0, amnh["raw_start"])

    starts = [
        x["raw_start"]
        for x in feature_matrix[label].values()
    ]

    return (1, min(starts) if starts else 10**9)


with open(OUT_WIDE, "w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=wide_columns,
        delimiter="\t",
    )

    writer.writeheader()

    for label in sorted(
        feature_matrix.keys(),
        key=feature_sort_key,
    ):

        row = {
            "feature_label": label,
            "feature_type": "",
        }

        starts = []
        ends = []

        for acc in ACCESSION_ORDER:

            item = feature_matrix[label].get(acc)

            if item is None:
                continue

            row["feature_type"] = item["feature_type"]

            row[f"{acc}_raw_start"] = item["raw_start"]
            row[f"{acc}_raw_end"] = item["raw_end"]
            row[f"{acc}_strand"] = item["strand"]
            row[f"{acc}_length"] = item["feature_length"]

            row[f"{acc}_aln_start"] = item[
                "alignment_start"
            ]
            row[f"{acc}_aln_end"] = item[
                "alignment_end"
            ]

            row[f"{acc}_aln_5prime"] = item[
                "alignment_5prime"
            ]
            row[f"{acc}_aln_3prime"] = item[
                "alignment_3prime"
            ]

            starts.append(item["alignment_start"])
            ends.append(item["alignment_end"])

        row["all_alignment_starts_same"] = (
            "YES"
            if len(starts) == len(ACCESSION_ORDER)
            and len(set(starts)) == 1
            else "NO"
        )

        row["all_alignment_ends_same"] = (
            "YES"
            if len(ends) == len(ACCESSION_ORDER)
            and len(set(ends)) == 1
            else "NO"
        )

        writer.writerow(row)


# ============================================================
# report adjacent gaps / overlaps
# ============================================================

overlap_columns = [
    "accession",
    "feature_1",
    "feature_1_type",
    "feature_1_start",
    "feature_1_end",

    "feature_2",
    "feature_2_type",
    "feature_2_start",
    "feature_2_end",

    "gap_bp",
    "overlap_bp",
    "relationship",
]


with open(OUT_OVERLAPS, "w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=overlap_columns,
        delimiter="\t",
    )

    writer.writeheader()

    for accession in ACCESSION_ORDER:

        feats = [
            x for x in all_features
            if x["accession"] == accession
        ]

        feats.sort(
            key=lambda x: (
                x["raw_start"],
                x["raw_end"],
            )
        )

        for f1, f2 in zip(feats[:-1], feats[1:]):

            # Number of unannotated bases between features
            gap = f2["raw_start"] - f1["raw_end"] - 1

            if gap > 0:
                overlap = 0
                relationship = "gap"

            elif gap == 0:
                overlap = 0
                relationship = "adjacent"

            else:
                overlap = -gap
                gap = 0
                relationship = "overlap"

            writer.writerow({
                "accession": accession,

                "feature_1": f1["feature_label"],
                "feature_1_type": f1["feature_type"],
                "feature_1_start": f1["raw_start"],
                "feature_1_end": f1["raw_end"],

                "feature_2": f2["feature_label"],
                "feature_2_type": f2["feature_type"],
                "feature_2_start": f2["raw_start"],
                "feature_2_end": f2["raw_end"],

                "gap_bp": gap,
                "overlap_bp": overlap,
                "relationship": relationship,
            })


# ============================================================
# done
# ============================================================

print("Wrote:")
print(f"  {OUT_LONG}")
print(f"  {OUT_WIDE}")
print(f"  {OUT_OVERLAPS}")