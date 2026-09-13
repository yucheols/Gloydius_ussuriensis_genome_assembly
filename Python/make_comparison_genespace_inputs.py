#!/usr/bin/env python3

import os
import re
import glob
from collections import defaultdict


BASE = "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"
GS   = "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

BED_DIR = os.path.join(GS, "bed")
PEP_DIR = os.path.join(GS, "peptide")
MAP_DIR = os.path.join(GS, "mapping")

os.makedirs(BED_DIR, exist_ok=True)
os.makedirs(PEP_DIR, exist_ok=True)
os.makedirs(MAP_DIR, exist_ok=True)


# ------------------------------------------------------------
# helper functions
# ------------------------------------------------------------

def parse_gff3_attributes(s):
    out = {}

    for item in s.strip().split(";"):
        if not item:
            continue

        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v

    return out


def parse_gtf_attributes(s):
    return {
        m.group(1): m.group(2)
        for m in re.finditer(r'(\S+)\s+"([^"]*)"', s)
    }


def natural_key(s):
    return [
        int(x) if x.isdigit() else x
        for x in re.split(r"(\d+)", s)
    ]


def read_fasta(path):
    records = {}

    current = None
    seq = []

    with open(path) as fh:

        for line in fh:
            line = line.rstrip("\n")

            if line.startswith(">"):

                if current is not None:
                    records[current] = "".join(seq)

                current = line[1:].split()[0]
                seq = []

            else:
                seq.append(line.strip())

        if current is not None:
            records[current] = "".join(seq)

    return records


# ------------------------------------------------------------
# parse annotation
# ------------------------------------------------------------

def parse_annotation(path):

    is_gtf = path.endswith(".gtf")

    gene_coords = {}
    transcript_coords = {}

    transcript_to_gene = {}

    # protein ID from CDS -> transcript ID(s)
    protein_to_transcript = defaultdict(set)

    with open(path) as fh:

        for line in fh:

            if not line.strip() or line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) != 9:
                continue

            chrom = f[0]
            feature = f[2]
            start = int(f[3])
            end = int(f[4])
            strand = f[6]

            if is_gtf:
                a = parse_gtf_attributes(f[8])
            else:
                a = parse_gff3_attributes(f[8])

            # ------------------------------------------------
            # GENE
            # ------------------------------------------------

            if feature == "gene":

                gene_id = a.get("gene_id") if is_gtf else a.get("ID")

                if gene_id:
                    gene_coords[gene_id] = (
                        chrom,
                        start,
                        end,
                        strand
                    )

            # ------------------------------------------------
            # TRANSCRIPT
            # ------------------------------------------------

            if feature in ("mRNA", "transcript"):

                if is_gtf:
                    tid = a.get("transcript_id")
                    gid = a.get("gene_id")

                else:
                    tid = a.get("ID")

                    parent = a.get("Parent")
                    gid = parent.split(",")[0] if parent else None

                if tid:

                    transcript_coords[tid] = (
                        chrom,
                        start,
                        end,
                        strand
                    )

                    if gid:
                        transcript_to_gene[tid] = gid

            # ------------------------------------------------
            # CDS
            #
            # Needed both for NCBI protein_id mapping and for
            # GTF records such as ToxCodAn where transcript
            # relationships may be encoded directly on CDS.
            # ------------------------------------------------

            if feature == "CDS":

                if is_gtf:

                    tid = a.get("transcript_id")
                    gid = a.get("gene_id")
                    protein_id = a.get("protein_id")

                    if tid and gid:
                        transcript_to_gene[tid] = gid

                else:

                    parent = a.get("Parent")
                    tid = parent.split(",")[0] if parent else None
                    protein_id = a.get("protein_id")

                if protein_id and tid:
                    protein_to_transcript[protein_id].add(tid)

    return (
        gene_coords,
        transcript_coords,
        transcript_to_gene,
        protein_to_transcript
    )


# ------------------------------------------------------------
# process one species
# ------------------------------------------------------------

def process_species(species, annotation, protein_fasta):

    proteins = read_fasta(protein_fasta)

    (
        gene_coords,
        transcript_coords,
        transcript_to_gene,
        protein_to_transcript
    ) = parse_annotation(annotation)

    # locus -> list of protein candidates
    candidates = defaultdict(list)

    unmapped = []
    ambiguous = []

    direct_transcript = 0
    via_protein_id = 0
    transcript_as_locus = 0

    for protein_id, sequence in proteins.items():

        transcript_id = None
        method = None

        # ------------------------------------------
        # Pattern 1:
        # FASTA ID is transcript/mRNA ID
        # ------------------------------------------

        if (
            protein_id in transcript_coords
            or protein_id in transcript_to_gene
        ):

            transcript_id = protein_id
            method = "fasta_id_is_transcript"

            direct_transcript += 1

        # ------------------------------------------
        # Pattern 2:
        # FASTA ID is NCBI protein_id
        # ------------------------------------------

        elif protein_id in protein_to_transcript:

            tids = sorted(protein_to_transcript[protein_id])

            if len(tids) != 1:
                ambiguous.append(
                    (protein_id, ",".join(tids))
                )
                continue

            transcript_id = tids[0]
            method = "fasta_id_is_protein_id"

            via_protein_id += 1

        else:
            unmapped.append(protein_id)
            continue

        # ------------------------------------------
        # resolve biological gene/locus
        # ------------------------------------------

        gene_id = transcript_to_gene.get(transcript_id)

        if gene_id and gene_id in gene_coords:

            locus_id = gene_id
            coords = gene_coords[gene_id]

        elif transcript_id in transcript_coords:

            # Important for Crotalus viridis:
            # its annotation has mRNAs but no gene features.
            locus_id = transcript_id
            coords = transcript_coords[transcript_id]

            transcript_as_locus += 1

        else:
            unmapped.append(protein_id)
            continue

        candidates[locus_id].append(
            (
                protein_id,
                transcript_id,
                sequence,
                method,
                coords
            )
        )

    # --------------------------------------------------------
    # longest protein per locus
    # --------------------------------------------------------

    representatives = {}

    for locus_id, values in candidates.items():

        values = sorted(
            values,
            key=lambda x: (
                -len(x[2]),
                x[0],
                x[1]
            )
        )

        representatives[locus_id] = values[0]

    # --------------------------------------------------------
    # sort genomic loci
    # --------------------------------------------------------

    loci = list(representatives)

    loci.sort(
        key=lambda locus: (
            natural_key(representatives[locus][4][0]),
            representatives[locus][4][1],
            representatives[locus][4][2],
            locus
        )
    )

    bed_path = os.path.join(
        BED_DIR,
        f"{species}.bed"
    )

    pep_path = os.path.join(
        PEP_DIR,
        f"{species}.fa"
    )

    map_path = os.path.join(
        MAP_DIR,
        f"{species}.longest_isoform.tsv"
    )

    # --------------------------------------------------------
    # BED
    # --------------------------------------------------------

    with open(bed_path, "w") as out:

        for locus in loci:

            protein_id, transcript_id, sequence, method, coords = \
                representatives[locus]

            chrom, start, end, strand = coords

            # GFF/GTF = 1-based
            # BED     = 0-based start
            out.write(
                f"{chrom}\t"
                f"{start - 1}\t"
                f"{end}\t"
                f"{locus}\n"
            )

    # --------------------------------------------------------
    # peptide FASTA
    #
    # Rename representative protein to locus ID so that
    # peptide FASTA ID exactly matches BED column 4.
    # --------------------------------------------------------

    with open(pep_path, "w") as out:

        for locus in loci:

            protein_id, transcript_id, sequence, method, coords = \
                representatives[locus]

            out.write(f">{locus}\n")

            for i in range(0, len(sequence), 60):
                out.write(sequence[i:i+60] + "\n")

    # --------------------------------------------------------
    # mapping table
    # --------------------------------------------------------

    with open(map_path, "w") as out:

        out.write(
            "locus_id\t"
            "representative_protein_id\t"
            "representative_transcript_id\t"
            "protein_length_aa\t"
            "chromosome\t"
            "start\t"
            "end\t"
            "strand\t"
            "n_protein_isoforms\t"
            "mapping_method\n"
        )

        for locus in loci:

            protein_id, transcript_id, sequence, method, coords = \
                representatives[locus]

            chrom, start, end, strand = coords

            out.write(
                f"{locus}\t"
                f"{protein_id}\t"
                f"{transcript_id}\t"
                f"{len(sequence)}\t"
                f"{chrom}\t"
                f"{start}\t"
                f"{end}\t"
                f"{strand}\t"
                f"{len(candidates[locus])}\t"
                f"{method}\n"
            )

    # --------------------------------------------------------
    # summary
    # --------------------------------------------------------

    genes_with_protein = len(
        set(candidates) & set(gene_coords)
    )

    explicit_genes_without_protein = (
        len(set(gene_coords) - set(candidates))
    )

    print()
    print("=" * 72)
    print(species)
    print("=" * 72)

    print(f"Input proteins:                     {len(proteins):,}")
    print(f"Annotation gene features:           {len(gene_coords):,}")
    print(f"Annotation transcript features:     {len(transcript_coords):,}")
    print(f"Mapped directly by transcript ID:   {direct_transcript:,}")
    print(f"Mapped through CDS protein_id:      {via_protein_id:,}")
    print(f"Representative loci written:        {len(loci):,}")
    print(f"Explicit genes represented:         {genes_with_protein:,}")
    print(f"Explicit genes without protein:     {explicit_genes_without_protein:,}")
    print(f"Transcript-as-locus assignments:    {transcript_as_locus:,}")
    print(f"Unmapped proteins:                  {len(unmapped):,}")
    print(f"Ambiguous protein mappings:         {len(ambiguous):,}")

    if unmapped:
        print("\nFirst unmapped proteins:")
        for x in unmapped[:10]:
            print(" ", x)

    if ambiguous:
        print("\nFirst ambiguous mappings:")
        for x in ambiguous[:10]:
            print(" ", x[0], "->", x[1])

    return len(loci), len(unmapped), len(ambiguous)


# ------------------------------------------------------------
# run all comparison species
# ------------------------------------------------------------

summary = []

for directory in sorted(glob.glob(os.path.join(BASE, "*"))):

    if not os.path.isdir(directory):
        continue

    species = os.path.basename(directory)

    anns = (
        glob.glob(os.path.join(directory, "*.gff3"))
        +
        glob.glob(os.path.join(directory, "*.gtf"))
    )

    fastas = glob.glob(
        os.path.join(directory, "*.protein.faa")
    )

    if len(anns) != 1 or len(fastas) != 1:

        print()
        print(f"SKIPPING {species}")
        print("Annotation candidates:", anns)
        print("Protein candidates:", fastas)

        continue

    loci, unmapped, ambiguous = process_species(
        species,
        anns[0],
        fastas[0]
    )

    summary.append(
        (
            species,
            loci,
            unmapped,
            ambiguous
        )
    )


print()
print()
print("=" * 72)
print("FINAL SUMMARY")
print("=" * 72)
print(
    "species\t"
    "representative_loci\t"
    "unmapped_proteins\t"
    "ambiguous_mappings"
)

for row in summary:
    print("\t".join(map(str, row)))