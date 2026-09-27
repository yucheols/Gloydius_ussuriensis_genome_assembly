#!/usr/bin/env python3

import sys

fasta = sys.argv[1]
keep_file = sys.argv[2]

with open(keep_file) as f:
    keep = {x.strip() for x in f if x.strip()}

chrom = None
pos = 0
run_start = None

def close_run(chrom, start, end):
    if chrom is not None and start is not None and end > start:
        print(f"{chrom}\t{start}\t{end}")

with open(fasta) as f:
    for line in f:
        line = line.rstrip("\n")

        if line.startswith(">"):
            if chrom in keep and run_start is not None:
                close_run(chrom, run_start, pos)

            chrom = line[1:].split()[0]
            pos = 0
            run_start = None
            continue

        if chrom not in keep:
            pos += len(line)
            continue

        for base in line:
            is_lower = base.islower()

            if is_lower and run_start is None:
                run_start = pos
            elif not is_lower and run_start is not None:
                close_run(chrom, run_start, pos)
                run_start = None

            pos += 1

    if chrom in keep and run_start is not None:
        close_run(chrom, run_start, pos)