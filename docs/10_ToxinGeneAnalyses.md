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

Regarding PLA2s, all five PLA2s in our curated toxin annotation set had significant Pfam matches to PLA2 family, and have lengths ranging from 118aa to 140aa. Although PLA2-2 is shorter than the rest, this in and of itself is not problematic.

Before extracting protein sequences for phylogenetic analyses, we will run one additional diagnostic step using the HMMER result we already obtained. The below script will check the following:
 - Presence/absence of metalloprotease domain in SVMP-4 and SVMP-15
 - Whether the P-II- and P-III-like proteins have the expected domain order
 - Whether PLA2-2 contains a near-complete PF00068 domain despite its shorter length
```sh
# ============================================================
# inspect SVMP and PLA2 domain coordinates and completeness
# ============================================================

ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

OUTDIR="$ROOT/phylogeny"

python - "$OUTDIR" <<'PY'
import sys
import csv
from pathlib import Path

outdir = Path(sys.argv[1])

target_domains = {
    "PF01421",
    "PF00200",
    "PF08516",
    "PF01562",
    "PF00068"
}

rows = []

for family in ("SVMP", "PLA2"):

    domtbl = outdir / f"Guss_{family}_Pfam.domtblout"

    with domtbl.open() as handle:
        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            f = line.split(maxsplit=22)
            pfam = f[1].split(".")[0]

            if pfam not in target_domains:
                continue

            hmm_length = int(f[2])
            hmm_start = int(f[15])
            hmm_end = int(f[16])

            coverage = (
                100 * (hmm_end - hmm_start + 1)
                / hmm_length
            )

            rows.append({
                "family": family,
                "gene_id": f[3],
                "pfam_id": pfam,
                "domain": f[0],
                "aa_start": int(f[17]),
                "aa_end": int(f[18]),
                "hmm_coverage_pct": round(coverage, 1),
                "i_evalue": f[12]
            })

rows.sort(
    key=lambda r: (
        r["family"],
        r["gene_id"],
        r["aa_start"]
    )
)

output = outdir / "Guss_SVMP_PLA2_domain_coordinates.tsv"

with output.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=list(rows[0]),
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(rows)

for family in ("SVMP", "PLA2"):
    print(f"\n===== {family} DOMAIN COORDINATES =====")

    for r in rows:
        if r["family"] != family:
            continue

        print(
            f"{r['gene_id']:<10} "
            f"{r['pfam_id']:<9} "
            f"{r['aa_start']:>4}-"
            f"{r['aa_end']:<4} "
            f"HMM coverage: "
            f"{r['hmm_coverage_pct']:>5.1f}% "
            f"E-value: {r['i_evalue']}"
        )

print("\nOutput:", output)
PY
```
The output will show that SVMP-15 has only 18.5% HMM coverage for the Pfam domain `PF01421`, which is a reprolysin domain. Regarding SVMP-4, there is 100% match to the Pfam domain `PF01562` (reprolysin propeptide). However, there is no reprolysin/metalloprotease domain.

For all five PLA2s, the HMM coverage is > 90% for the Pfam domain `PF00068`, which corresponds to PLA2. PLA2-1 is noteworthy compared to other PLA2s because its domain begins at amino acid 1, whereas the other proteins have approximately 17 residues before the detected domain. This suggests that PLA2-2 may lack the N-terminal signal peptide or that its annotation begins at the mature protein sequence. Alternatively, it could be an incomplete N-terminal prediction. However, its strong Pfam match (E-value 4.4 × 10⁻³²) and 93.6% HMM coverage support its identification as a PLA₂-family protein. We will thus retain all five PLA₂ sequences for our primary analysis, while flagging PLA2-2 for inspection.

#### Step 4: Prepare primary protein sequences
Now, we will prepare primary protein sequences for phylogenetic analyses. For SVMPs, we will extract the metalloprotease domain (`PF01421`) that is conserved across all SVMP classes. Our SVMPs frequently include `PF01421`, `PF13574`, `PF13582`, `PF13583`, and `PF13688` domains. This does not mean that these are five different domains. Instead, these are different Pfam profiles that recognize related metalloproteinase sequences. For example, a single SVMP protein might match multiple Reprolysin HMM profiles. 

For PLA2s, we will initialy retain the full predicted proteins, because strictly extracting `PF00068` would discard some C-terminal residues that may contain useful phylogenetic information. We can account for the apparently missing N-terminal region of PLA2-2 during alignment.

Run the script below:
```sh
# ============================================================
# prepare SVMP and PLA2 sequences for phylogenetic analyses
# ============================================================

# ------------------------------------------------------------
# 1. directories
# ------------------------------------------------------------

ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

OUTDIR="$ROOT/phylogeny"

# ------------------------------------------------------------
# 2. extract phylogenetic sequences
# ------------------------------------------------------------

python - "$OUTDIR" <<'PY'

import sys
import csv
from pathlib import Path

outdir = Path(sys.argv[1])

def read_fasta(path):
    sequences = {}
    name = None
    parts = []

    def save():
        if name is not None:
            if name in sequences:
                raise ValueError(f"Duplicate FASTA ID: {name}")
            sequences[name] = "".join(parts)

    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        save()

    return sequences

def write_fasta(path, sequences):
    with path.open("w") as handle:
        for gene, seq in sequences.items():
            handle.write(f">{gene}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

# ------------------------------------------------------------
# 3. extract complete SVMP metalloproteinase domains
# ------------------------------------------------------------

svmp = read_fasta(outdir / "Guss_SVMP_core.faa")

domtbl = outdir / "Guss_SVMP_Pfam.domtblout"

selected = {}
excluded = {}

with domtbl.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith("#"):
            continue

        f = line.split(maxsplit=22)

        pfam = f[1].split(".")[0]
        gene = f[3]

        if pfam != "PF01421":
            continue

        if gene not in svmp:
            raise ValueError(f"Unknown SVMP: {gene}")

        hmm_length = int(f[2])
        hmm_start = int(f[15])
        hmm_end = int(f[16])

        aa_start = int(f[17])
        aa_end = int(f[18])

        coverage = (
            hmm_end - hmm_start + 1
        ) / hmm_length

        if coverage >= 0.90:

            if gene in selected:
                raise ValueError(
                    f"Multiple qualifying domains: {gene}"
                )

            selected[gene] = svmp[gene][
                aa_start - 1:aa_end
            ]

for gene in svmp:
    if gene not in selected:
        excluded[gene] = svmp[gene]

if len(selected) != 15:
    raise ValueError(
        f"Expected 15 SVMP domains, found {len(selected)}"
    )

svmp_file = outdir / "Guss_SVMP_PF01421_primary.faa"

write_fasta(svmp_file, selected)

# ------------------------------------------------------------
# 4. prepare PLA2 proteins
# ------------------------------------------------------------

pla2 = read_fasta(outdir / "Guss_PLA2_core.faa")

if len(pla2) != 5:
    raise ValueError(
        f"Expected 5 PLA2 proteins, found {len(pla2)}"
    )

pla2_file = outdir / "Guss_PLA2_full_primary.faa"

write_fasta(pla2_file, pla2)

# ------------------------------------------------------------
# 5. print summary
# ------------------------------------------------------------

print("\n===== SVMP PRIMARY DATASET =====")

for gene, seq in sorted(selected.items()):
    print(f"{gene:<12} {len(seq):>4} aa")

print("\nExcluded from primary SVMP alignment:")

for gene in sorted(excluded):
    print(gene)

print("\n===== PLA2 PRIMARY DATASET =====")

for gene, seq in sorted(pla2.items()):
    print(f"{gene:<12} {len(seq):>4} aa")

print("\n===== OUTPUT FILES =====")
print(svmp_file)
print(pla2_file)

print("\nsequence preparation for phylogenetic analyses completed.")

PY
```

#### Step 5: Collect homologous sequences from other species

### Toxin gene microsynteny with other venomous snakes
Here, we will explore the organization of SVMPs and PLA2s in different venom snake lineages. We will use our *G. ussuriensis* assembly and chromosome-level genomes of *Crotalus adamanteus*, *Cerastes gasperettii*, and *Naja naja*. For this analysis, we will closely follow the methods used by Mochales-Riaño et al. (2025). As outlined in that study, we will use BLASTn to identify the genomic locations of homologous SVMPs and PLA2s and their nontoxin paralogs/conserved flanking genes. 

These genomes have their respective toxin gene anootation, and we already downloaded these assemblies for chromosomal macrosynteny analysis (see `11_Synteny.md`. I actually ran syteny analysis first, but documentation-wise I tought it would be better placed after these analyses). These assemblies are placed under the following directory:
```sh
/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny
```  
We will create the following directory structure:
```sh
synteny/
├── assemblies_synteny/          # existing reference genomes
│
└── toxin_microsynteny/
    ├── 01_reference_inventory/
    ├── 02_blast/
    │   ├── SVMP/
    │   └── PLA2/
    ├── 03_genomic_intervals/
    │   ├── SVMP/
    │   └── PLA2/
    ├── 04_alignments/
    │   ├── SVMP/
    │   └── PLA2/
    ├── 05_annotations/
    └── 06_figures/
```

Use the script below to create this directory stucture:
```sh

# ============================================================
# create SVMP and PLA2 comparative microsynteny directories
# ============================================================

# 1. directories

ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"

SYNTENY="$ROOT/synteny/toxin_microsynteny"

# 2. create directories

mkdir -p \
    "$SYNTENY/01_reference_inventory" \
    "$SYNTENY/02_blast/SVMP" \
    "$SYNTENY/02_blast/PLA2" \
    "$SYNTENY/03_genomic_intervals/SVMP" \
    "$SYNTENY/03_genomic_intervals/PLA2" \
    "$SYNTENY/04_alignments/SVMP" \
    "$SYNTENY/04_alignments/PLA2" \
    "$SYNTENY/05_annotations" \
    "$SYNTENY/06_figures"

# 3. verify directory structure

echo "===== MICROSYNTENY DIRECTORY STRUCTURE ====="

find "$SYNTENY" -type d | sort

echo
echo "===== COMPLETE ====="
```

#### Step 1: Locate the files needed
First, let's identify the actual FASTA files, assembly versions, and available GFF3/GTF annotations. Below is a read-only script that will do the job. This will not create any new folders or files.
```sh
# ============================================================
# inspect reference genomes and annotations
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"

# 2. inspect reference species
for SPECIES in \
    Cerastes_gasperettii \
    Crotalus_adamanteus \
    Naja_naja
do
    echo
    echo "============================================================"
    echo "$SPECIES"
    echo "============================================================"

    DIR="$ASSEMBLIES/$SPECIES"

    if [[ ! -d "$DIR" ]]; then
        echo "ERROR: Species directory missing."
        continue
    fi

    echo
    echo "===== DIRECTORY CONTENTS ====="

    ls -lah "$DIR"

    echo
    echo "===== GENOME AND ANNOTATION FILES ====="

    find "$DIR" -maxdepth 4 -type f \
        \( -iname '*.fa' \
        -o -iname '*.fasta' \
        -o -iname '*.fna' \
        -o -iname '*.fa.gz' \
        -o -iname '*.fasta.gz' \
        -o -iname '*.fna.gz' \
        -o -iname '*.gff' \
        -o -iname '*.gff3' \
        -o -iname '*.gtf' \
        -o -iname '*.gff.gz' \
        -o -iname '*.gff3.gz' \
        -o -iname '*.gtf.gz' \
        -o -iname '*.faa' \
        -o -iname '*.faa.gz' \
        -o -iname '*.gbff' \
        -o -iname '*.gbff.gz' \) \
        -print | sort

done

# 3. inspect existing assembly inventories
echo
echo "===== ASSEMBLY MANIFEST ====="

head -n 12 "$ASSEMBLIES/assembly_manifest.tsv"

echo
echo "===== ANNOTATION INVENTORY ====="

head -n 12 "$ASSEMBLIES/annotation_protein_inventory.tsv"


# 4. inspect G. ussuriensis toxin annotation
echo
echo "===== G. USSURIENSIS TOXIN ANNOTATION ====="

TOXIN="$ROOT/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

ls -lh \
    "$TOXIN/toxin_annotation.gtf" \
    "$TOXIN/toxin_annotation_cds.fasta" \
    "$TOXIN/toxin_annotation_pep.fasta"

echo
echo "===== COMPLETE ====="
```
The output of this script will show that our three comparison taxa all have the genome FASTA, gneome annotation GFF3, and protein FASTA files. It will also show that we have (as expected) *G. ussuriensis* ToxCodAn GTF, CDS FASTA, and protein FASTA.

#### Step 2: Extract genomic coordinates of *G. ussuriensis* SVMPs and PLA2s
This step will establish the focal genomic regions against which we will compare the other three species. Unlike the phylogenetic analysis, microsynteny should include our complete locus inventory, including SVMP-4, SVMP-15, and PLA2-4 where their genomic positions are supported.
```sh
# ============================================================
# inspect G. ussuriensis SVMP and PLA2 gene coordinates
# ============================================================

# activate funannotate conda env to access python
conda activate funannotate

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
TOXIN="$ROOT/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
GTF="$TOXIN/toxin_annotation.gtf"

# 2. inspect GTF structure
echo "===== GTF STRUCTURE ====="

grep -v '^#' "$GTF" | head -n 5

echo
echo "===== FEATURE TYPES ====="

awk -F '\t' '
    $0 !~ /^#/ && NF >= 9 {
        count[$3]++
    }
    END {
        for (feature in count)
            print feature, count[feature]
    }
' "$GTF" | sort

# 3. extract SVMP and PLA2 gene coordinates
echo
echo "===== TOXIN GENE COORDINATES ====="

python - "$GTF" <<'PY'

import sys
import re
from collections import defaultdict

gtf = sys.argv[1]

loci = defaultdict(
    lambda: {
        "chr": set(),
        "start": float("inf"),
        "end": 0,
        "features": set()
    }
)

with open(gtf) as handle:

    for line in handle:

        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) < 9:
            continue

        chrom = fields[0]
        feature = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        attributes = fields[8]

        match = re.search(
            r'(?:SVMP|PLA2)-\d+',
            attributes
        )

        if not match:
            continue

        gene = match.group(0)

        loci[gene]["chr"].add(chrom)
        loci[gene]["start"] = min(
            loci[gene]["start"], start
        )
        loci[gene]["end"] = max(
            loci[gene]["end"], end
        )
        loci[gene]["features"].add(feature)

def sort_key(gene):
    family, number = gene.split("-")
    return (family, int(number))

print(
    f"{'Gene':<12}"
    f"{'Chromosome':<24}"
    f"{'Start':>14}"
    f"{'End':>14}"
    f"{'Span (bp)':>14}"
)

for gene in sorted(loci, key=sort_key):

    row = loci[gene]

    span = row["end"] - row["start"] + 1

    print(
        f"{gene:<12}"
        f"{','.join(sorted(row['chr'])):<24}"
        f"{row['start']:>14,}"
        f"{row['end']:>14,}"
        f"{span:>14,}"
    )

print()
print("Total loci:", len(loci))

print(
    "SVMP loci:",
    sum(g.startswith("SVMP-") for g in loci)
)

print(
    "PLA2 loci:",
    sum(g.startswith("PLA2-") for g in loci)
)

if not loci:
    print(
        "WARNING: No matching gene IDs found. "
        "Inspect the GTF attributes shown above."
    )

PY

# 4. locate G. ussuriensis reference assembly
echo
echo "===== G. USSURIENSIS ASSEMBLY CANDIDATES ====="

find "$ROOT" -maxdepth 6 \
    \( -type f -o -type l \) \
    \( -name '*.fai' \
       -o -name '*genome*.fasta' \
       -o -name '*genome*.fa' \
       -o -name '*assembly*.fasta' \
       -o -name '*assembly*.fa' \) \
    -not -path '*/assemblies_synteny/*' \
    -print 2>/dev/null | head -n 40

echo
echo "===== COMPLETE ====="
```
The coordinate extraction output will tell us that within our G. ussuriensis assembly, 20 SVMPs spanning 747.2 kb are clustered on chromosome 9 (genomic intervals 5,375,758 bp - 6,122,976 bp) and 6 PLA2s spanning 34.5 kb (genomic intervals 278,191 bp - 312,734bp) are clustered on chromosome 16. A seventh PLA₂-labeled model, PLA2-1, is located separately on chromosome 10. However, recall that we reclassified this loci as PLA2 inhibitor-like protein rather than enzymatic PLA2 (see `09_ToxinAnnotationEval.md`).

Also, recall some important considerations from our toxin curation step. The 20 SVMP predictions should not be interpreted as 20 independently supported functional genes. In our previous curation, SVMP-5 and SVMP-6 represented overlapping alternative candidates, with SVMP-5 retained. SVMP-16 and SVMP-17 also have overlapping genomic coordinates and we flagged these as being uncertain gene models. Likewise, the PLA₂ array contains PLA2-4, which had a problematic gene model, and PLA2-5 and PLA2-6, which are distinct tandem loci despite encoding identical proteins.

#### Step 3: Verify the reference assembly and prepare genomic intervals
For this analysis, we will use unmasked genome (i.e. input FASTA we gave to Earl Grey). Recall that this genome is located in:
```sh
/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa
```

We will use the previously identified toxin clusters and initially include 250 kb of flanking sequence on each side. So for SVMPs and PLA2s, this becomes:
 - SVMPs: Chromosome 9. Genomic interval including the flanking sequences = 5,125,758–6,372,976 bp. Total length = 1,247,219 bp
 - PLA2s: Chromosome 16. Genomic interval including the flanking sequences = 28,191–562,734 bp. Total length = 534,544 bp

Use the script below will extract these genomic intervals from our *G. ussuriensis* reference genome.
```sh
# ============================================================
# extract G. ussuriensis SVMP and PLA2 genomic intervals
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
GENOME="$ROOT/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
OUTDIR="$ROOT/synteny/toxin_microsynteny/03_genomic_intervals/Gloydius_ussuriensis"

mkdir -p "$OUTDIR"

# 2. verify reference genome
echo "===== REFERENCE GENOME ====="

for FILE in "$GENOME" "${GENOME}.fai"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing $FILE"
        return 1 2>/dev/null || exit 1
    fi
done

command -v samtools

# 3. extract genomic intervals
echo
echo "===== EXTRACT SVMP REGION ====="

samtools faidx \
    "$GENOME" \
    "G_ussuri_chr9:5125758-6372976" \
    > "$OUTDIR/Guss_SVMP_region.fa"

echo
echo "===== EXTRACT PLA2 REGION ====="

samtools faidx \
    "$GENOME" \
    "G_ussuri_chr16:28191-562734" \
    > "$OUTDIR/Guss_PLA2_region.fa"

# 4. verify extracted sequences
echo
echo "===== EXTRACTED SEQUENCES ====="

for FILE in \
    "$OUTDIR/Guss_SVMP_region.fa" \
    "$OUTDIR/Guss_PLA2_region.fa"
do
    echo
    echo "$FILE"

    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing or empty output"
        continue
    fi

    head -n 1 "$FILE"

    awk '
        !/^>/ {
            total += length($0)
        }
        END {
            print "Sequence length:", total, "bp"
        }
    ' "$FILE"
done

# 5. locate genome gene annotations
echo
echo "===== G. USSURIENSIS GENE ANNOTATIONS ====="

find "$ROOT/annotation/funannotate" -maxdepth 5 -type f \
    \( -name '*.gff3' -o -name '*.gff' -o -name '*.gtf' \) \
    -print 2>/dev/null | head -n 30

echo
echo "===== COMPLETE ====="
```
This script will put the genomic region .fa file under the directory `/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny/03_genomic_intervals/Gloydius_ussuriensis`. The filenames will appear as `Guss_SVMP_region.fa` and `Guss_PLA2_region.fa`. Our directory structure has a separate folder for each gene family. We will thus move these files to their respective directory:
```sh
# ============================================================
# organize SVMP and PLA2 genomic interval files
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
SOURCE="$SYNTENY/03_genomic_intervals/Gloydius_ussuriensis"

# 2. organize genomic intervals
for FAMILY in SVMP PLA2
do
    SRC="$SOURCE/Guss_${FAMILY}_region.fa"
    DEST="$SYNTENY/03_genomic_intervals/$FAMILY/Guss_${FAMILY}_region.fa"

    mkdir -p "$(dirname "$DEST")"

    if [[ -s "$SRC" && ! -e "$DEST" ]]; then
        mv "$SRC" "$DEST"
    elif [[ -s "$DEST" ]]; then
        echo "Destination already exists: $DEST"
    else
        echo "WARNING: Missing genomic interval for $FAMILY"
    fi
done

# 3. verify files
echo "===== GENOMIC INTERVALS ====="

find "$SYNTENY/03_genomic_intervals" \
    -type f -name '*.fa' -print | sort
```

Before moving on to BLASTn, let's identify the flanking genes around the toxin clusters to establish that we're comparing homologous chromosomal neighborhoods, not merely regions containing similar toxin genes. For this, we will use our Funannotate annotation to identify conserved genes surrounding the *G. ussuriensis* SVMP and PLA₂ arrays.

We should first confirm that the GFF3 contains gene features on the expected chromosomes.
```sh
# ============================================================
# inspect Funannotate genes in SVMP and PLA2 regions
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
GFF="$ROOT/annotation/funannotate/update_results/Gloydius_ussuriensis_AMNH_21010.gff3"

# 2. inspect annotation
echo "===== FUNANNOTATE ANNOTATION ====="

ls -lh "$GFF"

echo
echo "===== EXAMPLE GENE RECORDS ====="

awk -F '\t' '
    $0 !~ /^#/ && $3 == "gene" {
        print
        if (++n == 5) exit
    }
' "$GFF"

# 3. inspect genes around toxin clusters
python - "$GFF" <<'PY'

import sys
from collections import Counter

gff = sys.argv[1]

regions = {
    "SVMP": ("G_ussuri_chr9", 5125758, 6372976),
    "PLA2": ("G_ussuri_chr16", 28191, 562734)
}

for family, (chrom, start, end) in regions.items():

    genes = []
    features = Counter()

    with open(gff) as handle:
        for line in handle:

            if not line.strip() or line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) < 9 or f[0] != chrom:
                continue

            a = int(f[3])
            b = int(f[4])

            if a > end or b < start:
                continue

            features[f[2]] += 1

            if f[2] == "gene":
                genes.append((a, b, f[6], f[8]))

    genes.sort()

    print(f"\n===== {family} REGION =====")

    print(f"Chromosome: {chrom}")
    print(f"Interval: {start:,}-{end:,}")
    print(f"Annotated genes: {len(genes)}")
    print(f"Feature types: {dict(features)}")

    print("\nFirst 5 genes:")

    for gene in genes[:5]:
        print(*gene, sep="\t")

    print("\nLast 5 genes:")

    for gene in genes[-5:]:
        print(*gene, sep="\t")

print("\n===== COMPLETE =====")
PY
```

The results should confirm that our Funannotate GFF3 contains gene models surrounding both toxin clusters. Now, we will move on to identifying what these gene models are. The script below will handle this task and create `05_annotations/Guss_SVMP_PLA2_flanking_genes.tsv` as an output.
```sh
# ============================================================
# identify SVMP and PLA2 flanking genes
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
GFF="$ROOT/annotation/funannotate/update_results/Gloydius_ussuriensis_AMNH_21010.gff3"
OUTDIR="$ROOT/synteny/toxin_microsynteny/05_annotations"

mkdir -p "$OUTDIR"

# 2. extract Funannotate gene annotations
python - "$GFF" "$OUTDIR" <<'PY'

import sys
import csv
from pathlib import Path
from urllib.parse import unquote
from collections import defaultdict

gff = Path(sys.argv[1])
outdir = Path(sys.argv[2])

regions = {
    "SVMP": {
        "chr": "G_ussuri_chr9",
        "region_start": 5125758,
        "region_end": 6372976,
        "cluster_start": 5375758,
        "cluster_end": 6122976
    },
    "PLA2": {
        "chr": "G_ussuri_chr16",
        "region_start": 28191,
        "region_end": 562734,
        "cluster_start": 278191,
        "cluster_end": 312734
    }
}

def parse_attributes(text):
    result = {}
    for part in text.strip().split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            result[key] = unquote(value)
    return result

genes = {}
transcripts = defaultdict(list)

with gff.open() as handle:
    for line in handle:

        if not line.strip() or line.startswith("#"):
            continue

        f = line.rstrip("\n").split("\t")

        if len(f) < 9:
            continue

        if f[2] not in ("gene", "mRNA"):
            continue

        attr = parse_attributes(f[8])

        if f[2] == "gene":

            gene_id = attr.get("ID")

            if not gene_id:
                continue

            genes[gene_id] = {
                "chr": f[0],
                "start": int(f[3]),
                "end": int(f[4]),
                "strand": f[6],
                "gene_name": attr.get("Name", ""),
                "product": attr.get("product", ""),
                "note": attr.get("Note", "")
            }

        elif f[2] == "mRNA":

            for parent in attr.get("Parent", "").split(","):
                if parent:
                    transcripts[parent].append(attr)

# Recover functional descriptions from transcripts
for gene_id, gene in genes.items():

    for tx in transcripts.get(gene_id, []):

        if not gene["product"]:
            gene["product"] = tx.get("product", "")

        if not gene["note"]:
            gene["note"] = tx.get("Note", "")

# 3. classify genes relative to toxin clusters

rows = []

for family, r in regions.items():

    for gene_id, gene in genes.items():

        if gene["chr"] != r["chr"]:
            continue

        start = gene["start"]
        end = gene["end"]

        if (
            end < r["region_start"]
            or start > r["region_end"]
        ):
            continue

        if end < r["cluster_start"]:
            position = "left_flank"
            distance = r["cluster_start"] - end

        elif start > r["cluster_end"]:
            position = "right_flank"
            distance = start - r["cluster_end"]

        else:
            position = "cluster_overlap"
            distance = 0

        rows.append({
            "family": family,
            "gene_id": gene_id,
            "chr": gene["chr"],
            "start": start,
            "end": end,
            "strand": gene["strand"],
            "position": position,
            "distance_bp": distance,
            "gene_name": gene["gene_name"],
            "product": gene["product"],
            "note": gene["note"]
        })

rows.sort(
    key=lambda x: (
        x["family"],
        x["start"]
    )
)

# 4. export annotation table

output = outdir / "Guss_SVMP_PLA2_flanking_genes.tsv"

with output.open("w", newline="") as handle:

    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "family",
            "gene_id",
            "chr",
            "start",
            "end",
            "strand",
            "position",
            "distance_bp",
            "gene_name",
            "product",
            "note"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

# 5. print closest flanking genes

for family in ("SVMP", "PLA2"):

    print(f"\n===== {family} REGION =====")

    subset = [
        r for r in rows
        if r["family"] == family
    ]

    for side in ("left_flank", "right_flank"):

        flank = sorted(
            [r for r in subset if r["position"] == side],
            key=lambda x: x["distance_bp"]
        )

        print(f"\n===== {side.upper()} =====")

        if not flank:
            print("No genes identified in working interval.")

        for r in flank[:8]:

            description = (
                r["product"]
                or r["gene_name"]
                or r["note"]
                or "UNANNOTATED"
            )

            print(
                f"{r['gene_id']:<14}"
                f"{r['start']:>10,}-"
                f"{r['end']:<10,} "
                f"{r['strand']} "
                f"{r['distance_bp']:>9,} bp  "
                f"{description[:90]}"
            )

    print(
        "\nFunannotate genes overlapping cluster:",
        sum(r["position"] == "cluster_overlap" for r in subset)
    )

print("\n===== OUTPUT FILE =====")
print(output)

PY
```
The output should indicate that we have successfully identified candidate flanking genes around both *G. ussuriensis* toxin clusters, but their functional identities remain unresolved (flanking genes marked as either `hypothetical protein` or `UNANNOTATED`). Therefore, the next step is to identify their homologs in *Cerastes gasperettii*, *Crotalus adamanteus*, and *Naja naja*. Because we already have annotated proteomes for these three comparative species, we will protein homology searches to establish conserved flanking genes.
```sh
# ============================================================
# identify protein sequences for SVMP and PLA2 flanking genes
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
ANNOTATION="$ROOT/annotation/funannotate"

# 2. locate protein FASTA files
echo "===== FUNANNOTATE PROTEIN FILES ====="

find "$ANNOTATION" -maxdepth 4 -type f \
    \( -iname '*proteins*.fa' \
       -o -iname '*proteins*.fasta' \
       -o -iname '*proteins*.faa' \) \
    -print | sort

# 3. inspect protein FASTA headers

echo
echo "===== PROTEIN FASTA HEADERS ====="

for FILE in \
    "$ANNOTATION/update_results/"*proteins* \
    "$ANNOTATION/annotate_results/"*proteins*
do
    if [[ ! -f "$FILE" ]]; then
        continue
    fi

    echo
    echo "FILE: $FILE"

    grep -m 5 '^>' "$FILE" || true
done

# 4. search candidate flanking genes
echo
echo "===== FLANKING GENE PROTEIN IDS ====="

for FILE in \
    "$ANNOTATION/update_results/"*proteins* \
    "$ANNOTATION/annotate_results/"*proteins*
do
    if [[ ! -f "$FILE" ]]; then
        continue
    fi

    echo
    echo "FILE: $FILE"

    grep '^>' "$FILE" |
        grep -E \
        'FUN_01451[123]|FUN_01452[67]|FUN_01805[78]|FUN_01807[23]' \
        || true
done

# 5. check homology search software
echo
echo "===== HOMOLOGY SEARCH SOFTWARE ====="

for TOOL in blastp makeblastdb blastn seqkit gffread
do
    if command -v "$TOOL" >/dev/null 2>&1; then
        echo "$TOOL: $(command -v "$TOOL")"
    else
        echo "$TOOL: NOT_FOUND"
    fi
done

echo
echo "===== COMPLETE ====="
```
The output of this script will verify that we have the Funannotate protein sequences needed for our flanking-gene analysis, and that all the necessary BLAST tools are installed. We will use the protein sequences from `annotation/funannotate/update_results/Gloydius_ussuriensis_AMNH_21010.proteins.fa`.

#### Step 4: BLASTp of flanking genes
Our immediate objective now is to identify candidate homologs of the *G. ussuriensis* flanking genes in the three comparison taxa.

The approach will be:
 1. Extract proteins corresponding to the closest eight annotated genes on each side of each toxin cluster, where available.
 2. Build BLAST databases from the existing comparative proteomes.
 3. Run BLASTp separately for SVMP and PLA₂ flanking-gene sets.
 4. Record the top 20 matches per query for subsequent orthology assessment.

 The script below will do the job:
 ```sh
 #!/bin/bash
#SBATCH --job-name=toxin_flank_blastp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================
# activate environment
# ============================================================
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

# ============================================================
# run blastp to identify SVMP and PLA2 flanking gene homologs
# ============================================================

# activate conda environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"
ANNOTATION="$ROOT/annotation/funannotate/update_results"
PROTEINS="$ANNOTATION/Gloydius_ussuriensis_AMNH_21010.proteins.fa"
FLANKING="$SYNTENY/05_annotations/Guss_SVMP_PLA2_flanking_genes.tsv"
BLASTDIR="$SYNTENY/02_blast"

mkdir -p "$BLASTDIR/SVMP" "$BLASTDIR/PLA2"

# 2. verify input files
echo "===== VERIFY INPUTS ====="
for FILE in "$PROTEINS" "$FLANKING"
do
    if [[ ! -s "$FILE" ]]; then
        echo "ERROR: Missing input: $FILE"
        exit 1
    fi
    echo "OK: $FILE"
done

for TOOL in python blastp makeblastdb
do
    command -v "$TOOL"
done

# 3. extract flanking gene proteins
echo
echo "===== EXTRACT FLANKING PROTEINS ====="

python - "$PROTEINS" "$FLANKING" "$BLASTDIR" <<'PY'

import sys
import csv
import re
from pathlib import Path
from collections import defaultdict

proteins_file = Path(sys.argv[1])
flanking_file = Path(sys.argv[2])
outdir = Path(sys.argv[3])

# Read protein FASTA
proteins = defaultdict(list)

header = None
seq = []

def save():
    if header is None:
        return

    parts = header.split()
    transcript = parts[0]

    gene = (
        parts[1]
        if len(parts) > 1
        else re.sub(r"-T\d+$", "", transcript)
    )

    proteins[gene].append(
        (transcript, "".join(seq))
    )

with proteins_file.open() as handle:
    for line in handle:
        line = line.strip()

        if line.startswith(">"):
            save()
            header = line[1:]
            seq = []
        elif line:
            seq.append(line)

    save()

# Read flanking-gene inventory
with flanking_file.open() as handle:
    rows = list(
        csv.DictReader(handle, delimiter="\t")
    )

for family in ("SVMP", "PLA2"):

    selected = []

    for side in ("left_flank", "right_flank"):

        subset = [
            row for row in rows
            if row["family"] == family
            and row["position"] == side
        ]

        subset.sort(
            key=lambda x: int(x["distance_bp"])
        )

        selected.extend(subset[:8])

    fasta_out = (
        outdir / family /
        f"Guss_{family}_flanking_proteins.faa"
    )

    manifest_out = (
        outdir / family /
        f"Guss_{family}_query_manifest.tsv"
    )

    manifest = []

    with fasta_out.open("w") as fasta:

        for row in selected:

            gene = row["gene_id"]

            candidates = proteins.get(gene, [])

            if not candidates:
                print(
                    f"WARNING: No protein found for {gene}"
                )
                continue

            # Prefer T1 when available.
            t1 = [
                x for x in candidates
                if x[0] == f"{gene}-T1"
            ]

            if t1:
                transcript, sequence = t1[0]
            else:
                transcript, sequence = max(
                    candidates,
                    key=lambda x: len(x[1])
                )

            if not sequence:
                print(
                    f"WARNING: Empty protein for {gene}"
                )
                continue

            fasta.write(f">{transcript} {gene}\n")

            for i in range(0, len(sequence), 80):
                fasta.write(sequence[i:i+80] + "\n")

            manifest.append({
                "gene_id": gene,
                "transcript_id": transcript,
                "position": row["position"],
                "distance_bp": row["distance_bp"],
                "length_aa": len(sequence)
            })

    with manifest_out.open("w", newline="") as handle:

        fields = [
            "gene_id",
            "transcript_id",
            "position",
            "distance_bp",
            "length_aa"
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(manifest)

    if not manifest:
        raise ValueError(
            f"No proteins recovered for {family}"
        )

    print(f"\n===== {family} =====")

    for row in manifest:
        print(
            f"{row['transcript_id']:<20}"
            f"{row['length_aa']:>5} aa  "
            f"{row['position']}"
        )

    print(f"\nTotal proteins: {len(manifest)}")
    print(f"FASTA: {fasta_out}")
    print(f"Manifest: {manifest_out}")

PY

# 4. prepare comparative protein BLAST databases
echo
echo "===== PREPARE BLAST DATABASES ====="

for SPECIES in \
    Cerastes_gasperettii \
    Crotalus_adamanteus \
    Naja_naja
do
    PROTEOME="$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa"

    DB="$BLASTDIR/${SPECIES}_proteins_db"

    if [[ ! -s "$PROTEOME" ]]; then
        echo "ERROR: Missing proteome: $PROTEOME"
        exit 1
    fi

    echo
    echo "Species: $SPECIES"

    if [[ ! -s "${DB}.pin" && ! -s "${DB}.pdb" ]]; then

        makeblastdb \
            -in "$PROTEOME" \
            -dbtype prot \
            -out "$DB"

    else
        echo "Existing BLAST database found."
    fi
done

# 5. run comparative BLASTP
echo
echo "===== RUN BLASTP ====="

for FAMILY in SVMP PLA2
do
    QUERY="$BLASTDIR/$FAMILY/Guss_${FAMILY}_flanking_proteins.faa"

    for SPECIES in \
        Cerastes_gasperettii \
        Crotalus_adamanteus \
        Naja_naja
    do
        DB="$BLASTDIR/${SPECIES}_proteins_db"

        OUTPUT="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.blastp.tsv"

        echo
        echo "$FAMILY vs $SPECIES"

        blastp \
            -query "$QUERY" \
            -db "$DB" \
            -evalue 1e-5 \
            -max_target_seqs 20 \
            -num_threads 2 \
            -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' \
            -out "$OUTPUT"

        echo "Output: $OUTPUT"
    done
done

# 6. summarize best BLASTP hits
echo
echo "===== BEST BLASTP HITS ====="

python - "$BLASTDIR" <<'PY'
import sys
import csv
from pathlib import Path
from collections import defaultdict

blastdir = Path(sys.argv[1])

species_list = (
    "Cerastes_gasperettii",
    "Crotalus_adamanteus",
    "Naja_naja"
)

for family in ("SVMP", "PLA2"):

    print(f"\n===== {family} =====")

    for species in species_list:

        path = (
            blastdir / family /
            f"Guss_{family}_vs_{species}.blastp.tsv"
        )

        hits = defaultdict(list)

        with path.open() as handle:
            reader = csv.reader(handle, delimiter="\t")

            for row in reader:
                if len(row) < 9:
                    continue

                hits[row[0]].append(row)

        print(f"\n--- {species} ---")

        if not hits:
            print("No hits passing the search threshold.")
            continue

        for query in sorted(hits):

            best = max(
                hits[query],
                key=lambda x: float(x[8])
            )

            print(
                f"{query:<20} "
                f"{best[1]:<28} "
                f"ID={float(best[2]):5.1f}% "
                f"Qcov={float(best[6]):5.1f}% "
                f"E={best[7]}"
            )

print("\n===== COMPLETE =====")
PY
 ```