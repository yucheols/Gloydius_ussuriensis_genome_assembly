#!/usr/bin/env python3

import sys
from collections import defaultdict


if len(sys.argv) != 6:
    sys.exit(
        "Usage:\n"
        "  make_genespace_representative.py "
        "<annotation.gff3> <proteins.fa> <out.bed> <out.fa> <mapping.tsv>"
    )


gff_file = sys.argv[1]
protein_file = sys.argv[2]
bed_out = sys.argv[3]
fasta_out = sys.argv[4]
map_out = sys.argv[5]


# ------------------------------------------------------------
# helper: parse GFF3 attributes
# ------------------------------------------------------------

def parse_attributes(attr_string):
    attrs = {}

    for field in attr_string.strip().split(";"):
        if not field:
            continue

        if "=" in field:
            key, value = field.split("=", 1)
            attrs[key] = value

    return attrs


# ------------------------------------------------------------
# read GFF3
#
# store:
#   gene coordinates
#   transcript -> parent gene relationship
# ------------------------------------------------------------

gene_coords = {}
transcript_to_gene = {}

with open(gff_file) as handle:

    for line in handle:

        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        chrom = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]

        attrs = parse_attributes(fields[8])

        if feature == "gene":

            gene_id = attrs.get("ID")

            if gene_id:
                gene_coords[gene_id] = (
                    chrom,
                    start,
                    end,
                    strand
                )

        elif feature in ("mRNA", "transcript"):

            transcript_id = attrs.get("ID")
            parent = attrs.get("Parent")

            if transcript_id and parent:

                # GFF3 can theoretically contain multiple parents.
                # Funannotate transcripts here have one parent gene.
                parent_gene = parent.split(",")[0]

                transcript_to_gene[transcript_id] = parent_gene


# ------------------------------------------------------------
# read protein FASTA
# ------------------------------------------------------------

proteins = {}

current_id = None
current_seq = []

with open(protein_file) as handle:

    for line in handle:

        line = line.rstrip("\n")

        if line.startswith(">"):

            if current_id is not None:
                proteins[current_id] = "".join(current_seq)

            current_id = line[1:].split()[0]
            current_seq = []

        else:
            current_seq.append(line.strip())

    if current_id is not None:
        proteins[current_id] = "".join(current_seq)


# ------------------------------------------------------------
# associate proteins with genes
# ------------------------------------------------------------

gene_proteins = defaultdict(list)

unmapped_proteins = []

for transcript_id, sequence in proteins.items():

    if transcript_id not in transcript_to_gene:
        unmapped_proteins.append(transcript_id)
        continue

    gene_id = transcript_to_gene[transcript_id]

    gene_proteins[gene_id].append(
        (
            transcript_id,
            sequence
        )
    )


# ------------------------------------------------------------
# select longest protein per gene
#
# ties are resolved deterministically by transcript ID
# ------------------------------------------------------------

representatives = {}

for gene_id, candidates in gene_proteins.items():

    candidates = sorted(
        candidates,
        key=lambda x: (-len(x[1]), x[0])
    )

    representatives[gene_id] = candidates[0]


# ------------------------------------------------------------
# retain only representative genes that also have coordinates
# ------------------------------------------------------------

valid_genes = [
    gene_id
    for gene_id in representatives
    if gene_id in gene_coords
]


# sort by chromosome, then genomic position
valid_genes.sort(
    key=lambda gene_id: (
        gene_coords[gene_id][0],
        gene_coords[gene_id][1],
        gene_coords[gene_id][2],
        gene_id
    )
)


# ------------------------------------------------------------
# write BED
#
# standard BED:
#   chromosome
#   0-based start
#   end
#   gene ID
#
# peptide FASTA IDs are renamed to gene IDs so BED column 4
# matches FASTA IDs exactly.
# ------------------------------------------------------------

with open(bed_out, "w") as bed:

    for gene_id in valid_genes:

        chrom, start, end, strand = gene_coords[gene_id]

        bed_start = start - 1

        bed.write(
            f"{chrom}\t"
            f"{bed_start}\t"
            f"{end}\t"
            f"{gene_id}\n"
        )


# ------------------------------------------------------------
# write representative protein FASTA
# ------------------------------------------------------------

with open(fasta_out, "w") as fasta:

    for gene_id in valid_genes:

        transcript_id, sequence = representatives[gene_id]

        fasta.write(f">{gene_id}\n")

        for i in range(0, len(sequence), 60):
            fasta.write(sequence[i:i+60] + "\n")


# ------------------------------------------------------------
# write mapping table
# ------------------------------------------------------------

with open(map_out, "w") as out:

    out.write(
        "gene_id\t"
        "representative_transcript\t"
        "protein_length_aa\t"
        "chromosome\t"
        "start\t"
        "end\t"
        "strand\t"
        "n_isoforms\n"
    )

    for gene_id in valid_genes:

        transcript_id, sequence = representatives[gene_id]
        chrom, start, end, strand = gene_coords[gene_id]

        out.write(
            f"{gene_id}\t"
            f"{transcript_id}\t"
            f"{len(sequence)}\t"
            f"{chrom}\t"
            f"{start}\t"
            f"{end}\t"
            f"{strand}\t"
            f"{len(gene_proteins[gene_id])}\n"
        )


# ------------------------------------------------------------
# summary
# ------------------------------------------------------------

genes_without_protein = set(gene_coords) - set(representatives)

representatives_without_coords = (
    set(representatives) - set(gene_coords)
)


print("============================================================")
print("GENESPACE representative-protein preparation")
print("============================================================")
print(f"GFF genes:                         {len(gene_coords):,}")
print(f"GFF transcripts:                   {len(transcript_to_gene):,}")
print(f"Input proteins:                    {len(proteins):,}")
print(f"Genes represented by proteins:     {len(representatives):,}")
print(f"GENESPACE loci written:            {len(valid_genes):,}")
print(f"Genes without protein:              {len(genes_without_protein):,}")
print(f"Proteins without transcript map:    {len(unmapped_proteins):,}")
print(f"Representative genes without coord: {len(representatives_without_coords):,}")
print("============================================================")

if unmapped_proteins:
    print("\nFirst unmapped proteins:")
    for x in unmapped_proteins[:10]:
        print(x)