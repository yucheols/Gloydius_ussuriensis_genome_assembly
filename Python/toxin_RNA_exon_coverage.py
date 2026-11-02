#!/usr/bin/env python3

import csv
import re
import subprocess
from collections import defaultdict

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------

GTF = "toxin_annotation.gtf"

BAM = (
    "/home/yshin/mendel-nas1/snake_genome_ass/"
    "G_ussuriensis_Chromo/annotation/"
    "toxin_gene_annotation_RNAdata_added/"
    "venom_RNA_mapping/venom_genome.sorted.bam"
)

TARGETS = [
    "SVSP-1",
    "SVSP-16",
    "Ficolin-4",
    "PLA2-4",
    "CTL-10",
    "HYAL-4",
    "CTL-15",
    "CTL-16",
    "CTL-17",
    "CTL-18"
]

EXON_OUT = "toxin_RNA_validation_CDS_exons.tsv"
GENE_OUT = "toxin_RNA_validation_CDS_summary.tsv"

# ------------------------------------------------------------
# Extract CDS coordinates from ToxCodAn GTF
# ------------------------------------------------------------

models = defaultdict(list)

with open(GTF) as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9 or fields[2] != "CDS":
            continue

        match = re.search(r'gene_id "([^"]+)"', fields[8])

        if not match:
            continue

        gene = match.group(1)

        if gene not in TARGETS:
            continue

        models[gene].append({
            "chrom": fields[0],
            "start": int(fields[3]),
            "end": int(fields[4]),
            "strand": fields[6]
        })

missing = [g for g in TARGETS if g not in models]

if missing:
    raise RuntimeError(f"Missing target genes: {missing}")

# ------------------------------------------------------------
# Calculate per-CDS RNA coverage
# ------------------------------------------------------------

exon_rows = []
gene_rows = []

for gene in TARGETS:

    exons = models[gene]

    chroms = {e["chrom"] for e in exons}
    strands = {e["strand"] for e in exons}

    if len(chroms) != 1 or len(strands) != 1:
        raise RuntimeError(f"Inconsistent GTF records: {gene}")

    chrom = next(iter(chroms))
    strand = next(iter(strands))

    # Genomic order for depth calculation
    genomic_exons = sorted(exons, key=lambda e: e["start"])

    # Transcript order for reporting exon numbers
    transcript_exons = sorted(
        exons,
        key=lambda e: e["start"],
        reverse=(strand == "-")
    )

    exon_rank = {
        (e["start"], e["end"]): i
        for i, e in enumerate(transcript_exons, start=1)
    }

    region_start = min(e["start"] for e in exons)
    region_end = max(e["end"] for e in exons)

    region = f"{chrom}:{region_start}-{region_end}"

    # -a: include zero-coverage positions
    # -q: minimum base quality
    # -Q: minimum mapping quality
    # -G: exclude unmapped, secondary, QC-failed,
    #     duplicate, and supplementary alignments

    cmd = [
        "samtools", "depth",
        "-a",
        "-q", "20",
        "-Q", "20",
        "-G", "0xF04",
        "-r", region,
        BAM
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    stats = [
        {"covered": 0, "depth_sum": 0}
        for _ in genomic_exons
    ]

    idx = 0
    observed_positions = 0

    for line in result.stdout.splitlines():

        fields = line.split("\t")

        if len(fields) < 3:
            continue

        pos = int(fields[1])
        depth = int(fields[2])
        observed_positions += 1

        while (
            idx < len(genomic_exons)
            and pos > genomic_exons[idx]["end"]
        ):
            idx += 1

        if idx >= len(genomic_exons):
            continue

        exon = genomic_exons[idx]

        if exon["start"] <= pos <= exon["end"]:
            stats[idx]["depth_sum"] += depth

            if depth > 0:
                stats[idx]["covered"] += 1

    expected_positions = region_end - region_start + 1

    if observed_positions != expected_positions:
        raise RuntimeError(
            f"{gene}: expected {expected_positions} depth rows, "
            f"received {observed_positions}"
        )

    total_length = 0
    total_covered = 0
    total_depth_sum = 0

    for exon, stat in zip(genomic_exons, stats):

        length = exon["end"] - exon["start"] + 1
        covered = stat["covered"]
        depth_sum = stat["depth_sum"]

        total_length += length
        total_covered += covered
        total_depth_sum += depth_sum

        exon_rows.append([
            gene,
            chrom,
            exon_rank[(exon["start"], exon["end"])],
            exon["start"],
            exon["end"],
            strand,
            length,
            covered,
            round(100 * covered / length, 2),
            round(depth_sum / length, 2)
        ])

    gene_rows.append([
        gene,
        chrom,
        strand,
        len(exons),
        total_length,
        total_covered,
        round(100 * total_covered / total_length, 2),
        round(total_depth_sum / total_length, 2)
    ])

    print(f"Processed: {gene}")

# ------------------------------------------------------------
# Write results
# ------------------------------------------------------------

with open(EXON_OUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Gene",
        "Chromosome",
        "CDS_rank",
        "Start",
        "End",
        "Strand",
        "CDS_length",
        "Covered_bases",
        "Coverage_percent",
        "Mean_depth"
    ])

    writer.writerows(exon_rows)

with open(GENE_OUT, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")

    writer.writerow([
        "Gene",
        "Chromosome",
        "Strand",
        "CDS_segments",
        "Total_CDS_bp",
        "Covered_bases",
        "Coverage_percent",
        "Mean_depth"
    ])

    writer.writerows(gene_rows)

print()
print("Completed RNA CDS coverage analysis.")
print(f"Exon-level results: {EXON_OUT}")
print(f"Gene-level results: {GENE_OUT}")