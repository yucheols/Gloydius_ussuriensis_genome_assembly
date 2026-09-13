#!/usr/bin/env python3

import os
import re
import glob
from collections import Counter


BASE = "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"


def parse_gff3_attributes(s):
    out = {}

    for item in s.strip().split(";"):
        if not item:
            continue

        if "=" in item:
            key, value = item.split("=", 1)
            out[key] = value

    return out


def parse_gtf_attributes(s):
    out = {}

    for match in re.finditer(r'(\S+)\s+"([^"]*)"', s):
        out[match.group(1)] = match.group(2)

    return out


def read_fasta_headers(path):
    primary = set()
    gene_tags = set()

    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()

            if not header:
                continue

            first = header.split()[0]
            primary.add(first)

            m = re.search(r'(?:^|\s)gene=([^\s]+)', header)

            if m:
                gene_tags.add(m.group(1))

    return primary, gene_tags


def inspect_annotation(path):

    ext = os.path.splitext(path)[1].lower()

    gene_ids = set()
    transcript_ids = set()
    cds_ids = set()

    protein_ids = set()
    orig_protein_ids = set()

    mrna_parents = set()
    cds_parents = set()

    gene_attributes = set()

    feature_counts = Counter()
    source_counts = Counter()

    with open(path) as fh:

        for line in fh:

            if not line.strip() or line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) != 9:
                continue

            source = f[1]
            feature = f[2]

            feature_counts[feature] += 1
            source_counts[source] += 1

            if ext == ".gtf":
                a = parse_gtf_attributes(f[8])
            else:
                a = parse_gff3_attributes(f[8])

            # ------------------------------------------
            # gene identifiers
            # ------------------------------------------

            if feature == "gene":

                if "ID" in a:
                    gene_ids.add(a["ID"])

                if "gene_id" in a:
                    gene_ids.add(a["gene_id"])

            # ------------------------------------------
            # transcript identifiers
            # ------------------------------------------

            if feature in ("mRNA", "transcript"):

                if "ID" in a:
                    transcript_ids.add(a["ID"])

                if "transcript_id" in a:
                    transcript_ids.add(a["transcript_id"])

                if "Parent" in a:
                    for x in a["Parent"].split(","):
                        mrna_parents.add(x)

                if "gene_id" in a:
                    mrna_parents.add(a["gene_id"])

            # GTF often stores transcript_id on CDS directly
            if "transcript_id" in a:
                transcript_ids.add(a["transcript_id"])

            if "gene_id" in a:
                gene_ids.add(a["gene_id"])

            # ------------------------------------------
            # CDS identifiers
            # ------------------------------------------

            if feature == "CDS":

                if "ID" in a:
                    cds_ids.add(a["ID"])

                if "Parent" in a:
                    for x in a["Parent"].split(","):
                        cds_parents.add(x)

                if "transcript_id" in a:
                    cds_parents.add(a["transcript_id"])

            # ------------------------------------------
            # protein mappings
            # ------------------------------------------

            if "protein_id" in a:
                protein_ids.add(a["protein_id"])

            if "orig_protein_id" in a:
                orig_protein_ids.add(a["orig_protein_id"])

            if "gene" in a:
                gene_attributes.add(a["gene"])

    return {
        "genes": gene_ids,
        "transcripts": transcript_ids,
        "cds": cds_ids,
        "protein_ids": protein_ids,
        "orig_protein_ids": orig_protein_ids,
        "mrna_parents": mrna_parents,
        "cds_parents": cds_parents,
        "gene_attributes": gene_attributes,
        "feature_counts": feature_counts,
        "source_counts": source_counts,
    }


for directory in sorted(glob.glob(os.path.join(BASE, "*"))):

    if not os.path.isdir(directory):
        continue

    species = os.path.basename(directory)

    anns = (
        glob.glob(os.path.join(directory, "*.gff3")) +
        glob.glob(os.path.join(directory, "*.gtf"))
    )

    fastas = glob.glob(os.path.join(directory, "*.faa"))

    if len(anns) != 1 or len(fastas) != 1:
        print()
        print("=" * 72)
        print(species)
        print("=" * 72)
        print("WARNING: expected exactly one annotation and one protein FASTA")
        print("annotations:", anns)
        print("proteins:", fastas)
        continue

    annotation = anns[0]
    fasta = fastas[0]

    A = inspect_annotation(annotation)
    fasta_ids, fasta_gene_tags = read_fasta_headers(fasta)

    def overlap(x):
        return len(fasta_ids & x)

    print()
    print("=" * 72)
    print(species)
    print("=" * 72)

    print(f"Annotation : {os.path.basename(annotation)}")
    print(f"Proteins   : {os.path.basename(fasta)}")

    print()
    print("ANNOTATION COUNTS")
    print(f"  gene IDs                    : {len(A['genes']):,}")
    print(f"  transcript IDs              : {len(A['transcripts']):,}")
    print(f"  CDS IDs                     : {len(A['cds']):,}")
    print(f"  protein_id values           : {len(A['protein_ids']):,}")
    print(f"  orig_protein_id values      : {len(A['orig_protein_ids']):,}")
    print(f"  CDS Parent IDs              : {len(A['cds_parents']):,}")

    print()
    print("FASTA COUNTS")
    print(f"  protein FASTA IDs           : {len(fasta_ids):,}")
    print(f"  FASTA gene= tags            : {len(fasta_gene_tags):,}")

    print()
    print("FASTA PRIMARY-ID OVERLAPS")
    print(f"  with annotation gene IDs    : {overlap(A['genes']):,}")
    print(f"  with transcript IDs         : {overlap(A['transcripts']):,}")
    print(f"  with CDS IDs                : {overlap(A['cds']):,}")
    print(f"  with protein_id             : {overlap(A['protein_ids']):,}")
    print(f"  with orig_protein_id        : {overlap(A['orig_protein_ids']):,}")
    print(f"  with CDS Parent IDs         : {overlap(A['cds_parents']):,}")

    if fasta_gene_tags:
        print()
        print("FASTA gene= TAG OVERLAPS")
        print(
            f"  gene= vs annotation genes   : "
            f"{len(fasta_gene_tags & A['genes']):,}"
        )
        print(
            f"  gene= vs gene attributes    : "
            f"{len(fasta_gene_tags & A['gene_attributes']):,}"
        )

    print()
    print("FEATURE COUNTS")
    for key, value in A["feature_counts"].most_common():
        print(f"  {key:25s} {value:,}")

    print()
    print("ANNOTATION SOURCES")
    for key, value in A["source_counts"].most_common():
        print(f"  {key:25s} {value:,}")