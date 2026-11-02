#!/usr/bin/env python3

import csv
import re
import subprocess
from collections import Counter, defaultdict

ROOT = (
    "/home/yshin/mendel-nas1/snake_genome_ass/"
    "G_ussuriensis_Chromo/annotation/"
    "toxin_gene_annotation_RNAdata_added/"
)

GTF = ROOT + "SRR35908235_ToxCodAn-Genome/toxin_annotation.gtf"
BAM = ROOT + "venom_RNA_mapping/venom_genome.sorted.bam"

TARGETS = [
    "SVSP-1", "SVSP-16", "Ficolin-4", "PLA2-4",
    "CTL-10", "HYAL-4", "CTL-15", "CTL-16",
    "CTL-17", "CTL-18"
]

MATCH_OPS = {"M", "=", "X"}
REF_OPS = {"M", "=", "X", "D", "N"}
MIN_ANCHOR = 8

models = defaultdict(list)

with open(GTF) as handle:
    for line in handle:
        if line.startswith("#") or not line.strip():
            continue

        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS":
            continue

        match = re.search(r'gene_id "([^"]+)"', f[8])
        if not match or match.group(1) not in TARGETS:
            continue

        models[match.group(1)].append(
            (f[0], int(f[3]), int(f[4]), f[6])
        )

expected_rows = []
alternative_rows = []

for gene in TARGETS:

    exons = models[gene]
    if not exons:
        raise RuntimeError(f"Missing GTF records: {gene}")

    chroms = {x[0] for x in exons}
    strands = {x[3] for x in exons}

    if len(chroms) != 1 or len(strands) != 1:
        raise RuntimeError(f"Inconsistent GTF: {gene}")

    chrom = exons[0][0]
    strand = exons[0][3]

    ordered = sorted(
        exons,
        key=lambda x: x[1],
        reverse=(strand == "-")
    )

    expected = {}

    for i in range(len(ordered) - 1):
        pair = sorted(ordered[i:i+2], key=lambda x: x[1])
        left, right = pair

        key = (left[2], right[1])
        expected[key] = f"CDS{i+1}-CDS{i+2}"

    region_start = min(x[1] for x in exons)
    region_end = max(x[2] for x in exons)

    region = f"{chrom}:{region_start}-{region_end}"

    all_junctions = Counter()
    high_junctions = Counter()

    cmd = [
        "samtools", "view",
        "-F", "0xF04",
        BAM, region
    ]

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        text=True
    )

    for line in proc.stdout:

        f = line.rstrip("\n").split("\t")

        if len(f) < 6 or f[5] == "*":
            continue

        mapq = int(f[4])
        refpos = int(f[3])

        cigar = f[5]
        ops = [
            (int(n), op)
            for n, op in re.findall(
                r'(\d+)([MIDNSHP=X])', cigar
            )
        ]

        if not ops:
            continue

        for i, (length, op) in enumerate(ops):

            if op == "N":

                left = refpos - 1
                right = refpos + length

                # Require aligned anchors on both sides.
                anchored = (
                    i > 0
                    and i + 1 < len(ops)
                    and ops[i-1][1] in MATCH_OPS
                    and ops[i+1][1] in MATCH_OPS
                    and ops[i-1][0] >= MIN_ANCHOR
                    and ops[i+1][0] >= MIN_ANCHOR
                )

                if (
                    anchored
                    and left >= region_start
                    and right <= region_end
                ):
                    key = (left, right)

                    all_junctions[key] += 1

                    if mapq >= 20:
                        high_junctions[key] += 1

            if op in REF_OPS:
                refpos += length

    if proc.wait() != 0:
        raise RuntimeError(f"samtools failed: {gene}")

    for (left, right), pair_name in expected.items():

        expected_rows.append([
            gene, chrom, strand, pair_name,
            left, right, right-left-1,
            all_junctions[(left, right)],
            high_junctions[(left, right)]
        ])

    for (left, right), count in high_junctions.items():

        if (left, right) in expected or count < 3:
            continue

        alternative_rows.append([
            gene, chrom, strand,
            left, right, right-left-1,
            all_junctions[(left, right)],
            count
        ])

    supported = sum(
        high_junctions[j] >= 3 for j in expected
    )

    print(
        f"{gene}: {supported}/{len(expected)} "
        "predicted junctions with >=3 MAPQ20 alignments"
    )

with open(
    "toxin_RNA_expected_junctions.tsv", "w", newline=""
) as out:

    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "Gene", "Chromosome", "Strand", "CDS_connection",
        "Left_exon_end", "Right_exon_start", "Intron_bp",
        "Primary_support", "MAPQ20_support"
    ])
    writer.writerows(expected_rows)

with open(
    "toxin_RNA_alternative_junctions.tsv", "w", newline=""
) as out:

    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "Gene", "Chromosome", "Strand",
        "Left_position", "Right_position", "Intron_bp",
        "Primary_support", "MAPQ20_support"
    ])
    writer.writerows(alternative_rows)

print("\nCompleted splice-junction analysis.")
print("Expected junctions:", len(expected_rows))
print("Alternative junctions:", len(alternative_rows))