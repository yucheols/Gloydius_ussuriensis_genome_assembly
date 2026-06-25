#!/usr/bin/env python3

import csv
from collections import defaultdict

# inputs
tidk_tsv = "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/find_telomere/Gloydius_ussuriensis_tidk_telomeric_repeat_windows.tsv"
fai = "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa.fai"

# settings
terminal_window = 50000   # summarize first and last 50 kb
threshold = 100           # counts >= this are treated as telomere-positive

# read scaffold lengths from .fai
lengths = {}
with open(fai) as f:
    for line in f:
        fields = line.rstrip("\n").split("\t")
        lengths[fields[0]] = int(fields[1])

# store summed terminal repeat counts
left_counts = defaultdict(int)
right_counts = defaultdict(int)

with open(tidk_tsv) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        scaffold = row["id"]
        window_end = int(row["window"])

        if scaffold not in lengths:
            continue

        count = int(row["forward_repeat_number"]) + int(row["reverse_repeat_number"])
        scaffold_len = lengths[scaffold]

        # first 50 kb
        if window_end <= terminal_window:
            left_counts[scaffold] += count

        # last 50 kb
        if window_end > scaffold_len - terminal_window:
            right_counts[scaffold] += count

# print table sorted by scaffold size
print("scaffold\tlength_bp\tsize_Mb\tleft_50kb_count\tright_50kb_count\tinterpretation")

for scaffold, length in sorted(lengths.items(), key=lambda x: x[1], reverse=True):
    left = left_counts[scaffold]
    right = right_counts[scaffold]

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
        f"{scaffold}\t{length}\t{length/1e6:.1f}\t"
        f"{left}\t{right}\t{interpretation}"
    )