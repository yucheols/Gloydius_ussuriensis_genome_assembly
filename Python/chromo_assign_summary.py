#!/usr/bin/env python3

import sys
from collections import defaultdict

paf = sys.argv[1]

# qname -> qlen
query_lengths = {}

# qname -> ref -> list of query intervals
intervals = defaultdict(lambda: defaultdict(list))

def merge_intervals(ints):
    if not ints:
        return []
    ints = sorted(ints)
    merged = [ints[0]]
    for start, end in ints[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged

def interval_total(ints):
    return sum(end - start for start, end in ints)

with open(paf) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        qname = fields[0]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])

        rname = fields[5]
        aln_len = int(fields[10])
        mapq = int(fields[11])

        # filter small/noisy alignments
        if aln_len < 10000:
            continue
        if mapq < 5:
            continue

        query_lengths[qname] = qlen
        intervals[qname][rname].append((qstart, qend))

print(
    "query_scaffold\tquery_length\tbest_reference_chromosome\t"
    "best_unique_query_bp\ttotal_unique_query_bp\t"
    "best_fraction_of_aligned\tbest_fraction_of_scaffold"
)

for qname in sorted(query_lengths, key=lambda x: query_lengths[x], reverse=True):
    ref_intervals = intervals[qname]

    if not ref_intervals:
        continue

    ref_unique = {}
    all_intervals = []

    for rname, ints in ref_intervals.items():
        merged = merge_intervals(ints)
        bp = interval_total(merged)
        ref_unique[rname] = bp
        all_intervals.extend(ints)

    total_unique = interval_total(merge_intervals(all_intervals))
    best_ref, best_bp = max(ref_unique.items(), key=lambda x: x[1])

    best_fraction_of_aligned = best_bp / total_unique if total_unique else 0
    best_fraction_of_scaffold = best_bp / query_lengths[qname] if query_lengths[qname] else 0

    print(
        f"{qname}\t{query_lengths[qname]}\t{best_ref}\t"
        f"{best_bp}\t{total_unique}\t"
        f"{best_fraction_of_aligned:.4f}\t{best_fraction_of_scaffold:.4f}"
    )