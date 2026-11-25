## 10) Toxin gene evolution
In this section, we will use predicted venom toxin genes annotated on our reference genome and perform several different analyses. We will focus on two toxin families that are known to be abundant in the venom protein of G. ussuriensis: snake venom metalloproteases (SVMPs) and phospholipase A2 (PLA2). These are also known to be prominent components of the crotaline venom toxin repertoire more generally.

### SVMP and PLA2 phylogenies
#### Step 1: Validation of predicted proteins
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

Regarding PLA2s, we excluded PLA2-1 because a BLASTp run suggested that this is a PLA2 inhibitor-like protein rather than an enzymatic PLA2. In addtion, we excluded PLA2-4 because it contained a premature stop codon in its predicted CDS. 

Thus, this gives us 17 predicted SVMPs and 5 PLA2s to work with. One thing that is notable from our basic predicted protein QC is that the amino acid length is quite variable for SVMPs, ranging from 241aa to 641aa, while PLA2 lengths are more constrained around 118aa-154aa.

Another important thing to keep in mind before moving forward: `PASS` indicates that the amino acid sequence contains no obvious character-level problems. It does not independently establish protein completeness, correct domain architecture, catalytic activity, or venom function.

#### Step 2: Generate curated protein FASTA files
Run the script below to generate curated protein FASTA files that can be used for phylogenetic analyses. The script will generate the two working FASTAs, one for SVMPs and another for PLA2s, validate the expected counts, and identify identical protein sequences. However, it will not collapse identical proteins.
```sh
# ============================================================
# prepare core SVMP and PLA2 protein datasets
# ============================================================
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

OUTDIR="$ROOT/phylogeny"
export OUTDIR

python - <<'PY'
import os
import csv
from pathlib import Path
from collections import defaultdict

outdir = Path(os.environ["OUTDIR"])

# ----------------------------------------------
# manual curation decisions
# ----------------------------------------------

excluded = {
    "SVMP-6": "Redundant with SVMP-5",
    "PLA2-1": "PLA2 inhibitor-like candidate"
}

noncore = {
    "SVMP-16": "Uncertain gene model",
    "SVMP-17": "Uncertain gene model",
    "PLA2-4": "Internal stop codon"
}

expected = {
    "SVMP": 17,
    "PLA2": 5
}

valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

# ----------------------------------------------
# read FASTA
# ----------------------------------------------

def read_fasta(path):
    records = {}
    name = None
    seq = []

    def save():
        if name is not None:
            if name in records:
                raise ValueError(f"Duplicate ID: {name}")
            records[name] = "".join(seq).upper()

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                name = line[1:].split()[0]
                seq = []
            else:
                if name is None:
                    raise ValueError("Invalid FASTA")
                seq.append(line)
        save()

    return records

def write_fasta(path, records):
    with path.open("w") as handle:
        for name, seq in records.items():
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

# ----------------------------------------------
# process protein families
# ----------------------------------------------

manifest = []

for family in ("SVMP", "PLA2"):

    source = outdir / f"Guss_{family}_all.faa"

    if not source.is_file():
        raise FileNotFoundError(source)

    records = read_fasta(source)

    core = {}
    others = {}

    for gene, original_seq in sorted(records.items()):

        seq = original_seq.rstrip("*")

        if gene in excluded:
            status = "EXCLUDED"
            reason = excluded[gene]
        elif gene in noncore:
            status = "NONCORE"
            reason = noncore[gene]
        else:
            status = "CORE"
            reason = "Retained in curated core set"

        if status == "CORE":
            invalid = set(seq) - valid_aa
            if not seq or invalid:
                raise ValueError(
                    f"Core sequence QC failed: {gene}, "
                    f"invalid residues: {sorted(invalid)}"
                )
            core[gene] = seq
        else:
            others[gene] = seq

        manifest.append({
            "gene_id": gene,
            "family": family,
            "length_aa": len(seq),
            "status": status,
            "reason": reason
        })

    if len(core) != expected[family]:
        raise ValueError(
            f"{family}: expected {expected[family]} "
            f"core proteins, found {len(core)}"
        )

    core_file = outdir / f"Guss_{family}_core.faa"
    review_file = outdir / f"Guss_{family}_excluded_review.faa"

    write_fasta(core_file, core)
    write_fasta(review_file, others)

    print(f"\n===== {family} =====")
    print(f"Original predictions: {len(records)}")
    print(f"Core proteins: {len(core)}")
    print(f"Excluded/noncore: {len(others)}")
    print(f"Core FASTA: {core_file}")

    # identify identical protein sequences
    groups = defaultdict(list)

    for gene, seq in core.items():
        groups[seq].append(gene)

    identical = [
        genes for genes in groups.values()
        if len(genes) > 1
    ]

    print("\nIdentical protein groups:")

    if identical:
        for genes in identical:
            print("  " + ", ".join(genes))
    else:
        print("  None")

# ----------------------------------------------
# export curation manifest
# ----------------------------------------------

manifest_file = outdir / "Guss_SVMP_PLA2_curation.tsv"

with manifest_file.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "gene_id",
            "family",
            "length_aa",
            "status",
            "reason"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(manifest)

print("\n===== COMPLETE =====")
print("Curation manifest:", manifest_file)
PY

# ----------------------------------------------
# check available domain annotation software
# ----------------------------------------------

echo
echo "===== DOMAIN ANNOTATION SOFTWARE ====="

for tool in interproscan.sh hmmscan hmmsearch
do
    if command -v "$tool" >/dev/null 2>&1; then
        command -v "$tool"
    else
        echo "$tool: not found in current PATH"
    fi
done
```

The output of this script will show that there are no identical protein pairs among the 17 retained SVMPs, and one identical protein pairs among the 5 retained PLA2s. 

#### Step 3: Protein domain identification using Pfam and HMMER
Here, we will use HMMER and Pfam to annotate the protein domains in our core SVMP and PLA2 sets. Pfam contains models representing conserved protein families and domains. What HMMER basically does is comparing the predicted aa sequences against these Pfam models to identify statistically supported domain matches.  
```sh

# ============================================================
# Pfam domain annotation for SVMP and PLA2
# ============================================================

conda activate funannotate

export FUNANNOTATE_DB="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"

# ------------------------------------------------------------
# 1. directories
# ------------------------------------------------------------

ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

OUTDIR="$ROOT/phylogeny"
mkdir -p "$OUTDIR"

# ------------------------------------------------------------
# 2. locate Pfam database
# ------------------------------------------------------------

echo "===== PFAM DATABASE ====="

echo "FUNANNOTATE_DB: ${FUNANNOTATE_DB:-NOT_SET}"

PFAM="${PFAM_HMM:-}"

if [[ -z "$PFAM" && -n "${FUNANNOTATE_DB:-}" ]]; then
    PFAM="$FUNANNOTATE_DB/Pfam-A.hmm"
fi

if [[ -z "$PFAM" || ! -f "$PFAM" ]]; then
    echo "ERROR: Pfam-A.hmm not found."
    echo "Set PFAM_HMM to the full database path."
    exit 1
fi

echo "Database: $PFAM"

# hmmscan requires a pressed HMM database
for ext in h3f h3i h3m h3p
do
    if [[ ! -f "${PFAM}.${ext}" ]]; then
        echo "ERROR: Missing HMM index: ${PFAM}.${ext}"
        echo "Database needs hmmpress indexing."
        exit 1
    fi
done

# ------------------------------------------------------------
# 3. run hmmscan
# ------------------------------------------------------------

for FAMILY in SVMP PLA2
do
    echo
    echo "===== ${FAMILY} DOMAIN ANNOTATION ====="

    INPUT="$OUTDIR/Guss_${FAMILY}_core.faa"
    DOMTBL="$OUTDIR/Guss_${FAMILY}_Pfam.domtblout"
    REPORT="$OUTDIR/Guss_${FAMILY}_Pfam_report.txt"

    if [[ ! -s "$INPUT" ]]; then
        echo "ERROR: Missing protein FASTA: $INPUT"
        exit 1
    fi

    hmmscan \
        --cpu 2 \
        --cut_ga \
        --noali \
        --domtblout "$DOMTBL" \
        "$PFAM" \
        "$INPUT" \
        > "$REPORT"

    echo "Completed: $FAMILY"
done

# ------------------------------------------------------------
# 4. parse domain annotations
# ------------------------------------------------------------

python - "$OUTDIR" <<'PY'

import sys
import csv
from pathlib import Path
from collections import defaultdict

outdir = Path(sys.argv[1])

families = ("SVMP", "PLA2")

domain_names = {
    "PF01421": "Reprolysin",
    "PF00200": "Disintegrin",
    "PF08516": "ADAM_CR",
    "PF00068": "Phospholip_A2_1"
}

all_hits = []
summary = []

for family in families:

    fasta = outdir / f"Guss_{family}_core.faa"
    domtbl = outdir / f"Guss_{family}_Pfam.domtblout"

    # Read protein IDs and lengths
    proteins = {}
    name = None
    seq = []

    def save():
        if name is not None:
            proteins[name] = len("".join(seq))

    with fasta.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                save()
                name = line[1:].split()[0]
                seq = []
            elif line:
                seq.append(line)
        save()

    hits = defaultdict(list)

    # Parse HMMER domtblout
    with domtbl.open() as handle:
        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.split(maxsplit=22)

            accession = fields[1].split(".")[0]
            gene = fields[3]

            row = {
                "gene_id": gene,
                "family": family,
                "pfam_id": accession,
                "domain": fields[0],
                "ali_start": int(fields[17]),
                "ali_end": int(fields[18]),
                "i_evalue": fields[12],
                "domain_score": fields[13]
            }

            all_hits.append(row)
            hits[gene].append(row)

    print(f"\n===== {family} DOMAIN SUMMARY =====")

    for gene in sorted(proteins):

        gene_hits = hits.get(gene, [])

        descriptions = [
            f"{h['pfam_id']}:{h['domain']}"
            for h in gene_hits
        ]

        domain_string = (
            "; ".join(descriptions)
            if descriptions else "NO_PFAM_HIT"
        )

        print(
            f"{gene:<12} "
            f"{proteins[gene]:>4} aa  "
            f"{domain_string}"
        )

        summary.append({
            "gene_id": gene,
            "family": family,
            "length_aa": proteins[gene],
            "n_domains": len(gene_hits),
            "domains": domain_string
        })

# Export detailed hits
detail_file = outdir / "Guss_SVMP_PLA2_Pfam_domains.tsv"

with detail_file.open("w", newline="") as handle:
    fields = [
        "gene_id", "family", "pfam_id",
        "domain", "ali_start", "ali_end",
        "i_evalue", "domain_score"
    ]
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(all_hits)

# Export per-protein summary
summary_file = outdir / "Guss_SVMP_PLA2_domain_summary.tsv"

with summary_file.open("w", newline="") as handle:
    fields = [
        "gene_id", "family", "length_aa",
        "n_domains", "domains"
    ]
    writer = csv.DictWriter(
        handle, fieldnames=fields, delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(summary)

print("\n===== OUTPUT FILES =====")
print(detail_file)
print(summary_file)

print("\nDomain annotation completed.")

PY
```
The output of this script will look something like:
```sh
===== SVMP DOMAIN SUMMARY =====
SVMP-1        624 aa  PF01421:Reprolysin; PF08516:ADAM_CR; PF13688:Reprolysin_5; PF00200:Disintegrin; PF01562:Pep_M12B_propep; PF13582:Reprolysin_3; PF13574:Reprolysin_2; PF13583:Reprolysin_4
SVMP-10       482 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF01562:Pep_M12B_propep; PF13582:Reprolysin_3; PF00200:Disintegrin; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-11       611 aa  PF01421:Reprolysin; PF00200:Disintegrin; PF13688:Reprolysin_5; PF08516:ADAM_CR; PF01562:Pep_M12B_propep; PF13582:Reprolysin_3; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-12       614 aa  PF01421:Reprolysin; PF00200:Disintegrin; PF13688:Reprolysin_5; PF08516:ADAM_CR; PF01562:Pep_M12B_propep; PF13582:Reprolysin_3; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-13       487 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF00200:Disintegrin; PF13582:Reprolysin_3; PF01562:Pep_M12B_propep; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-14       487 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF01562:Pep_M12B_propep; PF13582:Reprolysin_3; PF00200:Disintegrin; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-15       247 aa  PF01562:Pep_M12B_propep; PF01421:Reprolysin
SVMP-18       603 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF00200:Disintegrin; PF01562:Pep_M12B_propep; PF08516:ADAM_CR; PF13582:Reprolysin_3; PF13574:Reprolysin_2; PF13583:Reprolysin_4
SVMP-19       612 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF13582:Reprolysin_3; PF08516:ADAM_CR; PF00200:Disintegrin; PF01562:Pep_M12B_propep; PF13574:Reprolysin_2; PF13583:Reprolysin_4
SVMP-2        618 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF00200:Disintegrin; PF08516:ADAM_CR; PF01562:Pep_M12B_propep; PF13582:Reprolysin_3; PF13574:Reprolysin_2; PF13583:Reprolysin_4
SVMP-20       612 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF00200:Disintegrin; PF13582:Reprolysin_3; PF01562:Pep_M12B_propep; PF08516:ADAM_CR; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-3        532 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF00200:Disintegrin; PF13582:Reprolysin_3; PF01562:Pep_M12B_propep; PF13574:Reprolysin_2; PF13583:Reprolysin_4
SVMP-4        241 aa  PF01562:Pep_M12B_propep
SVMP-5        487 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF01562:Pep_M12B_propep; PF00200:Disintegrin; PF13582:Reprolysin_3; PF13583:Reprolysin_4
SVMP-7        476 aa  PF01421:Reprolysin; PF01562:Pep_M12B_propep; PF13688:Reprolysin_5; PF00200:Disintegrin; PF13582:Reprolysin_3; PF13583:Reprolysin_4
SVMP-8        478 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF01562:Pep_M12B_propep; PF00200:Disintegrin; PF13582:Reprolysin_3; PF13583:Reprolysin_4; PF13574:Reprolysin_2
SVMP-9        482 aa  PF01421:Reprolysin; PF13688:Reprolysin_5; PF13582:Reprolysin_3; PF01562:Pep_M12B_propep; PF00200:Disintegrin; PF13583:Reprolysin_4; PF13574:Reprolysin_2

===== PLA2 DOMAIN SUMMARY =====
PLA2-2        118 aa  PF00068:Phospholip_A2_1
PLA2-3        140 aa  PF00068:Phospholip_A2_1
PLA2-5        138 aa  PF00068:Phospholip_A2_1
PLA2-6        138 aa  PF00068:Phospholip_A2_1
PLA2-7        138 aa  PF00068:Phospholip_A2_1
```
Regarding SVMPs, SVMP-3, 5, 7, 8, 9, 10, 13, 14 are similar to SVMP P-II class because they possess a propeptide domain, reprolysin/metalloprotease domain, and a disintegrin domain. These SVMPs are 603aa-624aa in length. On the other hand, SVMP-1, 2, 11, 12, 18, 19, 20 can be classified as SVMP P-III class becasue they possess a propeptide domain, reprolysin/metalloprotease domain, disintegrin/DIS-like domain, and a cys-rich domain (ADAM-CR). These SVMPs range in length from 476aa to 532aa. See the below figure from Casewell et al. (2011) for the linear domain structure of SVMPs.
![alt text](/docs/etc/svmp_domains.png) 

SVMP-4 (241aa) and SVMP-15 (247aa) require special attention. First, SVMP-4 only has a propeptide domain and has no significant matches to the canonical metalloprotease domain. SVMP-15 has both domainss but its short length raises the question of whether this metalloprotease domain is complete or truncated.

### Toxin gene microsynteny with *Crotalus adamanteus* assembly