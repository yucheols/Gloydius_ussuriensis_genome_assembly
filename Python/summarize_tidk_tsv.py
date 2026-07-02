#!/usr/bin/env python3

import csv
import os
import sys
from collections import defaultdict

# ------------------------------------------------------------
# Summarize terminal telomeric repeat signal from TIDK output
# for the final curated Gloydius ussuriensis chromosome assembly.
#
# Usage:
#   python summarize_tidk_terminal_curated.py TIDK_WINDOWS_TSV
#
# Example:
#   python summarize_tidk_terminal_curated.py \
#     final_QC/tidk/Gloydius_curated_TTAGGG_telomeric_repeat_windows.tsv \
#     > final_QC/tidk/Gloydius_curated_terminal_telomere_summary.tsv
# ------------------------------------------------------------

# Required input: TIDK windows TSV generated from the curated assembly
if len(sys.argv) != 2:
    sys.stderr.write(
        "Usage:\n"
        "  python summarize_tidk_terminal_curated.py TIDK_WINDOWS_TSV\n\n"
        "Example:\n"
        "  python summarize_tidk_terminal_curated.py "
        "final_QC/tidk/Gloydius_curated_TTAGGG_telomeric_repeat_windows.tsv "
        "> final_QC/tidk/Gloydius_curated_terminal_telomere_summary.tsv\n"
    )
    sys.exit(1)

tidk_tsv = sys.argv[1]

# Final curated assembly index
curated_dir = "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated"
fai = os.path.join(
    curated_dir,
    "Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa.fai"
)

# Settings
terminal_window = 50000   # summarize first and last 50 kb
threshold = 100           # counts >= this are treated as telomere-positive

# ------------------------------------------------------------
# Check input files
# ------------------------------------------------------------

if not os.path.exists(tidk_tsv):
    sys.stderr.write(f"ERROR: TIDK TSV not found:\n{tidk_tsv}\n")
    sys.exit(1)

if not os.path.exists(fai):
    sys.stderr.write(f"ERROR: FASTA index not found:\n{fai}\n")
    sys.stderr.write("Run: samtools faidx Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa\n")
    sys.exit(1)

# ------------------------------------------------------------
# Read chromosome lengths from curated assembly .fai
# ------------------------------------------------------------

lengths = {}

with open(fai) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 2:
            continue
        lengths[fields[0]] = int(fields[1])

# ------------------------------------------------------------
# Store summed terminal repeat counts
# ------------------------------------------------------------

left_counts = defaultdict(int)
right_counts = defaultdict(int)

seen_ids = set()

with open(tidk_tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")

    required_cols = {"id", "window", "forward_repeat_number", "reverse_repeat_number"}
    missing_cols = required_cols - set(reader.fieldnames or [])

    if missing_cols:
        sys.stderr.write("ERROR: TIDK TSV is missing required columns:\n")
        sys.stderr.write(", ".join(sorted(missing_cols)) + "\n")
        sys.stderr.write("Found columns:\n")
        sys.stderr.write(", ".join(reader.fieldnames or []) + "\n")
        sys.exit(1)

    for row in reader:
        chrom = row["id"]
        seen_ids.add(chrom)

        if chrom not in lengths:
            continue

        window_end = int(row["window"])
        count = int(row["forward_repeat_number"]) + int(row["reverse_repeat_number"])
        chrom_len = lengths[chrom]

        # First 50 kb
        if window_end <= terminal_window:
            left_counts[chrom] += count

        # Last 50 kb
        if window_end > chrom_len - terminal_window:
            right_counts[chrom] += count

# ------------------------------------------------------------
# Warn if TIDK names do not match curated FASTA names
# ------------------------------------------------------------

matched = set(lengths).intersection(seen_ids)

if len(matched) == 0:
    sys.stderr.write(
        "WARNING: No TIDK sequence IDs match the curated FASTA .fai names.\n"
        "This usually means the TIDK file was generated from the old 69-scaffold assembly.\n"
        "Rerun TIDK on the curated FASTA:\n"
        f"{os.path.join(curated_dir, 'Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa')}\n"
    )

# ------------------------------------------------------------
# Print terminal telomere summary
# ------------------------------------------------------------

print(
    "chromosome\tlength_bp\tsize_Mb\t"
    "left_50kb_count\tright_50kb_count\t"
    "left_positive\tright_positive\tinterpretation"
)

for chrom, length in sorted(lengths.items(), key=lambda x: x[1], reverse=True):
    left = left_counts[chrom]
    right = right_counts[chrom]

    left_pos = left >= threshold
    right_pos = right >= threshold

    if left_pos and right_pos:
        interpretation = "both ends"
    elif left_pos and not right_pos:
        interpretation = "left end only"
    elif not left_pos and right_pos:
        interpretation = "right end only"
    else:
        interpretation = "no strong terminal signal"

    print(
        f"{chrom}\t{length}\t{length / 1e6:.2f}\t"
        f"{left}\t{right}\t"
        f"{left_pos}\t{right_pos}\t{interpretation}"
    )