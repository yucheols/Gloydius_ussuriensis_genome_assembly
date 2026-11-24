## 10) Toxin gene evolution
In this section, we will use predicted venom toxin genes annotated on our reference genome and perform several different analyses. We will focus on two toxin families that are known to be abundant in the venom protein of G. ussuriensis: snake venom metalloproteases (SVMPs) and phospholipase A2 (PLA2). These are also known to be prominent components of the crotaline venom toxin repertoire more generally.

### Protein phylogenies of toxin-associated gene families
We will start by analyzing the evolutionary histories of SVMPs and PLA2s. In our core toxin annotation (containing 110 curated genomic loci), we have 17 SVMP proteins with genomic loci predicted on chromosome 9, and 5 PLA2 proteins with genomic loci predicted on chromosome 16. 

We will start by validating these predictions. The script below will inspect the toxin CDS FASTA files, SVMP/PLA2 sequence counting, protein length and sequence integrity checks:
```sh
# make sure to activate the funannotate conda env to use python

# ============================================================
# SVMP and PLA2 protein sequence QC
# ============================================================

# 1. set directories and input files
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

PEP="$ROOT/toxin_annotation_pep.fasta"
CDS="$ROOT/toxin_annotation_cds.fasta"
GTF="$ROOT/toxin_annotation.gtf"

OUTDIR="$ROOT/phylogeny"
mkdir -p "$OUTDIR"

# 2. verify input files
echo "===== INPUT FILES ====="

for file in "$PEP" "$CDS" "$GTF"
do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing file: $file"
        exit 1
    fi
    ls -lh "$file"
done


# 3. inspect protein FASTA headers
echo
echo "===== FIRST 15 FASTA HEADERS ====="
grep '^>' "$PEP" | head -n 15 || true

echo
echo "===== SVMP HEADERS ====="
grep '^>SVMP-' "$PEP" || true

echo
echo "===== PLA2 HEADERS ====="
grep '^>PLA2-' "$PEP" || true


# 4. protein sequence QC and extraction
python - "$PEP" "$OUTDIR" <<'PY'

import sys
import csv
from pathlib import Path
from collections import Counter

pep_file = Path(sys.argv[1])
outdir = Path(sys.argv[2])

families = ("SVMP", "PLA2")

# ------------------------------------------------------------
# Read protein FASTA
# ------------------------------------------------------------

sequences = {}
headers = {}

name = None
parts = []

def save_record():
    if name is not None:
        if name in sequences:
            raise ValueError(f"Duplicate FASTA ID: {name}")
        sequences[name] = "".join(parts).upper()

with pep_file.open() as handle:
    for line in handle:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            save_record()
            name = line[1:].split()[0]
            headers[name] = line[1:]
            parts = []
        else:
            if name is None:
                raise ValueError("Invalid FASTA format")
            parts.append(line)

    save_record()

print("\n===== TOTAL PROTEIN RECORDS =====")
print(len(sequences))

# ------------------------------------------------------------
# Define preliminary sequence QC criteria
# ------------------------------------------------------------

valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

qc_rows = []

# ------------------------------------------------------------
# Assess and export SVMP and PLA2 proteins
# ------------------------------------------------------------

for family in families:

    selected = {
        gene: seq
        for gene, seq in sequences.items()
        if gene.startswith(family + "-")
    }

    print(f"\n===== {family} PROTEIN QC =====")
    print(f"{'Gene ID':<25} {'Length (aa)':>12}  Status")
    print("-" * 65)

    output_fasta = outdir / f"Guss_{family}_all.faa"

    with output_fasta.open("w") as handle:

        for gene, seq in sorted(selected.items()):

            protein = seq.rstrip("*")

            flags = []

            if not protein:
                flags.append("EMPTY_SEQUENCE")

            if "*" in protein:
                flags.append("INTERNAL_STOP")

            if any(aa not in valid_aa for aa in protein):
                flags.append("NONSTANDARD_AA")

            status = ",".join(flags) if flags else "PASS"

            print(
                f"{gene:<25} "
                f"{len(protein):>12}  "
                f"{status}"
            )

            # Preserve all sequences, including flagged models.
            # Remove only terminal stop symbols for protein FASTA.
            handle.write(f">{gene}\n")

            for i in range(0, len(protein), 80):
                handle.write(protein[i:i+80] + "\n")

            qc_rows.append({
                "gene_id": gene,
                "family": family,
                "length_aa": len(protein),
                "qc_status": status
            })

    print(f"\nTotal {family} proteins: {len(selected)}")
    print(f"FASTA: {output_fasta}")

# ------------------------------------------------------------
# Export QC table
# ------------------------------------------------------------

qc_file = outdir / "Guss_SVMP_PLA2_protein_QC.tsv"

with qc_file.open("w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "gene_id",
            "family",
            "length_aa",
            "qc_status"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(qc_rows)

# ------------------------------------------------------------
# Summarize QC
# ------------------------------------------------------------

print("\n===== QC SUMMARY =====")

for status, count in Counter(
    row["qc_status"] for row in qc_rows
).items():
    print(f"{status}: {count}")

print("\n===== OUTPUT FILES =====")

print(outdir / "Guss_SVMP_all.faa")
print(outdir / "Guss_PLA2_all.faa")
print(qc_file)

print("\npreliminary protein QC completed.")

PY
```
The output will show that there are a total of 120 protein records, 20 of which are SVMPs and 7 are PLA2s. 

During our curation step to generate a "core set" of annotated toxins, we excluded SVMP-6 because it was a redundant annotation with respect to SVMP-5, which had a longer amino acid sequence (487 aa vs. 430 aa of SVMP-6) and had stronger RNA transcript evidence. We further excluded SVMP-16 and SVMP-17 because these were either overlapping or uncertain gene models.  

Regarding PLA2s, we excluded PLA2-1 because a BLASTp run suggested that this is a PLA2 inhibitor-like protein rather than an enzymatic PLA2. In addtion, we excluded  


### Toxin gene microsynteny with *Crotalus adamanteus* assembly