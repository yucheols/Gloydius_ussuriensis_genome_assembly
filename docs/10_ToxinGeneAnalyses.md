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
# SVMP and PLA2 flanking-gene protein homology
# ============================================================

# activate conda env
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
GUSS_GENOME="$ROOT/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
EVIDENCE="$BLASTDIR/flanking_homology_evidence"

mkdir -p "$BLASTDIR/SVMP" "$BLASTDIR/PLA2" "$EVIDENCE"

# 2. verify inputs and software
printf '\n===== VERIFY INPUTS =====\n'
for FILE in "$PROTEINS" "$FLANKING" "$GUSS_GENOME"; do
    if [[ ! -s "$FILE" ]]; then echo "ERROR: missing input: $FILE" >&2; exit 1; fi
    echo "OK: $FILE"
done

for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
    for FILE in "$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa" "$ASSEMBLIES/$SPECIES/${SPECIES}.annotation.gff3"; do
        if [[ ! -s "$FILE" ]]; then echo "ERROR: missing input: $FILE" >&2; exit 1; fi
    done
done

for TOOL in python blastp makeblastdb tblastn samtools; do command -v "$TOOL" >/dev/null || exit 1; done
if [[ ! -s "${GUSS_GENOME}.fai" ]]; then samtools faidx "$GUSS_GENOME"; fi

# 3. extract Gloydius flanking proteins 
printf '\n===== EXTRACT FLANKING PROTEINS =====\n'
python - "$PROTEINS" "$FLANKING" "$BLASTDIR" <<'PY'
import csv, re, sys
from pathlib import Path
from collections import defaultdict
proteins_file, flanking_file, outdir = map(Path, sys.argv[1:])
proteins = defaultdict(list)
header, seq = None, []
def save():
    if header is None: return
    parts = header.split()
    transcript = parts[0]
    gene = parts[1] if len(parts) > 1 and re.fullmatch(r'FUN_\d+', parts[1]) else re.sub(r'-T\d+$', '', transcript)
    proteins[gene].append((transcript, ''.join(seq)))
with proteins_file.open() as handle:
    for line in handle:
        line = line.strip()
        if line.startswith('>'):
            save()
            header, seq = line[1:], []
        elif line: seq.append(line)
    save()
with flanking_file.open() as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
for family in ('SVMP', 'PLA2'):
    selected = []
    for side in ('left_flank', 'right_flank'):
        subset = sorted((r for r in rows if r['family'] == family and r['position'] == side), key=lambda r: int(r['distance_bp']))
        selected.extend(subset[:8])
    manifest = []
    fasta_out = outdir / family / f'Guss_{family}_flanking_proteins.faa'
    manifest_out = outdir / family / f'Guss_{family}_query_manifest.tsv'
    with fasta_out.open('w') as fasta:
        for row in selected:
            gene = row['gene_id']
            candidates = proteins.get(gene, [])
            if not candidates:
                print(f'WARNING: No protein found for {gene}')
                continue
            transcript, sequence = next((c for c in candidates if c[0] == f'{gene}-T1'), max(candidates, key=lambda x: len(x[1])))
            if not sequence: continue
            fasta.write(f'>{transcript} {gene}\n')
            for i in range(0, len(sequence), 80): fasta.write(sequence[i:i+80] + '\n')
            manifest.append(dict(gene_id=gene, transcript_id=transcript, position=row['position'], distance_bp=row['distance_bp'], length_aa=len(sequence)))
    if not manifest: raise RuntimeError(f'No flanking proteins extracted for {family}')
    with manifest_out.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['gene_id', 'transcript_id', 'position', 'distance_bp', 'length_aa'], delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest)
    print(f'{family}: {len(manifest)} proteins; {fasta_out}')
PY

# 4. prepare comparative BLAST databases and reciprocal database
printf '\n===== PREPARE PROTEIN BLAST DATABASES =====\n'
for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
    PROTEOME="$ASSEMBLIES/$SPECIES/${SPECIES}.protein.faa"
    DB="$BLASTDIR/${SPECIES}_proteins_db"
    if [[ ! -s "${DB}.pin" && ! -s "${DB}.pdb" ]]; then makeblastdb -in "$PROTEOME" -dbtype prot -out "$DB"; fi
done
GUSS_DB="$EVIDENCE/Gloydius_proteins_db"
if [[ ! -s "${GUSS_DB}.pin" && ! -s "${GUSS_DB}.pdb" ]]; then makeblastdb -in "$PROTEINS" -dbtype prot -out "$GUSS_DB"; fi

# 5. forward BLASTp 
printf '\n===== FORWARD BLASTP =====\n'
for FAMILY in SVMP PLA2; do
    QUERY="$BLASTDIR/$FAMILY/Guss_${FAMILY}_flanking_proteins.faa"
    for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        OUTPUT="$BLASTDIR/$FAMILY/Guss_${FAMILY}_vs_${SPECIES}.blastp.tsv"
        blastp -query "$QUERY" -db "$BLASTDIR/${SPECIES}_proteins_db" -evalue 1e-5 -max_target_seqs 20 -num_threads 2 -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' -out "$OUTPUT"
        echo "$FAMILY vs $SPECIES: $(wc -l < "$OUTPUT") HSPs"
    done
done

# 6. extract top forward subject sequences for reciprocal BLASTp
printf '\n===== PREPARE RECIPROCAL QUERIES =====\n'
python - "$BLASTDIR" "$ASSEMBLIES" "$EVIDENCE" <<'PY'
import csv, sys
from pathlib import Path
from collections import defaultdict
blastdir, assemblies, evidence = map(Path, sys.argv[1:])
def fasta_records(path):
    header, lines = None, []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None: yield header.split()[0], ''.join(lines)
                header, lines = line[1:], []
            elif line.strip(): lines.append(line.strip())
        if header is not None: yield header.split()[0], ''.join(lines)
for family in ('SVMP', 'PLA2'):
    for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        path = blastdir / family / f'Guss_{family}_vs_{species}.blastp.tsv'
        hits = defaultdict(list)
        with path.open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) >= 9: hits[row[0]].append((row[1], float(row[8])))
        subjects = set()
        for query, matches in hits.items():
            for subject, score in sorted(matches, key=lambda r: -r[1])[:5]: subjects.add(subject)
        ref = assemblies / species / f'{species}.protein.faa'
        records = {header: seq for header, seq in fasta_records(ref) if header in subjects}
        missing = subjects - records.keys()
        if missing: print(f'WARNING: {species} {family}: {len(missing)} BLAST subjects not recovered from FASTA')
        out = evidence / f'{family}_{species}_reciprocal_queries.faa'
        with out.open('w') as handle:
            for header, seq in sorted(records.items()):
                handle.write(f'>{header}\n')
                for i in range(0, len(seq), 80): handle.write(seq[i:i+80] + '\n')
        print(f'{family} {species}: {len(records)} reciprocal protein queries')
PY

# 7. reciprocal BLASTp against the Gloydius proteome
printf '\n===== RECIPROCAL BLASTP =====\n'
for FAMILY in SVMP PLA2; do
    for SPECIES in Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        QUERY="$EVIDENCE/${FAMILY}_${SPECIES}_reciprocal_queries.faa"
        OUTPUT="$EVIDENCE/${FAMILY}_${SPECIES}_reciprocal.blastp.tsv"
        if [[ -s "$QUERY" ]]; then
            blastp -query "$QUERY" -db "$GUSS_DB" -evalue 1e-5 -max_target_seqs 5 -num_threads 2 -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' -out "$OUTPUT"
        else
            : > "$OUTPUT"
        fi
    done
done

# 8. build protein reference panel from annotated Cerastes and Crotalus proteins
# make ADAM selection deliberately conservative == ambiguous source IDs are reported
printf '\n===== EXTRACT TARGETED ANCHOR REFERENCE PROTEINS =====\n'
python - "$ASSEMBLIES" "$EVIDENCE" <<'PY'
import csv, re, sys
from pathlib import Path
from urllib.parse import unquote
assemblies, evidence = map(Path, sys.argv[1:])
def fasta_records(path):
    header, lines = None, []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                if header is not None: yield header, ''.join(lines)
                header, lines = line[1:], []
            elif line: lines.append(line)
        if header is not None: yield header, ''.join(lines)
cer = list(fasta_records(assemblies / 'Cerastes_gasperettii' / 'Cerastes_gasperettii.protein.faa'))
patterns = {'STC1': r'(?i)STC[-_ ]?1|STANNIOCALCIN[-_ ]?1', 'NEFM': r'(?i)(?<![A-Z0-9])NEFM(?![A-Z0-9])|NEUROFILAMENT[^;]*MEDIUM', 'NEFL': r'(?i)(?<![A-Z0-9])NEFL(?![A-Z0-9])|NEUROFILAMENT[^;]*LIGHT', 'ADAM': r'(?i)(?<![A-Z0-9])(?:RNA-|GENE-)?ADAM(?![0-9A-Z])', 'OTUD3': r'(?i)OTUD[-_ ]?3', 'MUL1': r'(?i)(?<![A-Z0-9])MUL[-_ ]?1(?![A-Z0-9])', 'VWA5B1': r'(?i)VWA5B1', 'PLA2G2F': r'(?i)PLA2G2F|PLA2[-_ ]?G2F'}
found, audit = [], []
for symbol, pattern in patterns.items():
    matches = [(h, s) for h, s in cer if re.search(pattern, h)]
    if len(matches) == 1:
        h, s = matches[0]
        found.append(('Cerastes', symbol, h, s))
        audit.append((symbol, 'Cerastes_gasperettii', 'unique_header_match', h, len(s)))
    else:
        audit.append((symbol, 'Cerastes_gasperettii', f'{len(matches)}_header_matches_review', ';'.join(h[:100] for h, _ in matches[:5]), ''))
        print(f'WARNING: Cerastes {symbol}: {len(matches)} candidate FASTA records; not automatically choosing one')

# Resolve Crotalus ADAM28 using the gene ID in its GFF and protein_id on CDS.
sp = 'Crotalus_adamanteus'
crot = list(fasta_records(assemblies / sp / f'{sp}.protein.faa'))
crot_by_id = {h.split()[0]: (h, s) for h, s in crot}
gff = assemblies / sp / f'{sp}.annotation.gff3'
features, parents, aliases = {}, {}, set()
with gff.open() as handle:
    for line in handle:
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9: continue
        a = {k: unquote(v) for part in f[8].split(';') if '=' in part for k, v in [part.split('=', 1)]}
        fid = a.get('ID', '')
        if fid:
            features[fid] = a
            parents[fid] = a.get('Parent', '').split(',')
        if fid == 'gene-NXF25_014549': aliases.add(fid)
def belongs_to_adam(fid):
    seen = set()
    todo = [fid]
    while todo:
        node = todo.pop()
        if node in seen: continue
        if node == 'gene-NXF25_014549': return True
        seen.add(node)
        todo.extend(parents.get(node, []))
    return False
ids = set()
for fid, a in features.items():
    if belongs_to_adam(fid):
        for key in ('protein_id', 'Name', 'ID'):
            if a.get(key): ids.update(a[key].split(','))
matches = [crot_by_id[k] for k in sorted(ids) if k in crot_by_id]
if len(matches) != 1:
    matches = [(h, s) for h, s in crot if re.search(r'(?i)(?<![A-Z0-9])ADAM[-_ ]?28(?![A-Z0-9])|NXF25_014549', h)]
if len(matches) == 1:
    h, s = matches[0]
    found.append(('Crotalus', 'ADAM28', h, s))
    audit.append(('ADAM28', sp, 'unique_protein_match', h, len(s)))
else:
    audit.append(('ADAM28', sp, f'{len(matches)}_protein_matches_review', ';'.join(h[:100] for h, _ in matches[:5]), ''))
    print('WARNING: Crotalus ADAM28 protein not uniquely recoverable; Cerastes ADAM may still provide the search')

# Keep protein panels family-specific to avoid irrelevant cross-family matches.
for family in ('SVMP', 'PLA2'):
    members = [('Cerastes', 'STC1'), ('Cerastes', 'NEFM'), ('Cerastes', 'NEFL'), ('Cerastes', 'ADAM'), ('Crotalus', 'ADAM28')] if family == 'SVMP' else [('Cerastes', 'OTUD3'), ('Cerastes', 'MUL1'), ('Cerastes', 'VWA5B1'), ('Cerastes', 'PLA2G2F')]
    selected = [(sp, symbol, h, seq) for sp, symbol, h, seq in found if (sp, symbol) in members]
    with (evidence / f'{family}_anchor_reference_proteins.faa').open('w') as handle:
        for species, symbol, header, seq in selected:
            handle.write(f'>{species}_{symbol}\n')
            for i in range(0, len(seq), 80): handle.write(seq[i:i+80] + '\n')
    if not selected: print(f'WARNING: no reference proteins for {family}')
with (evidence / 'anchor_reference_protein_audit.tsv').open('w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t')
    writer.writerow(['anchor_symbol', 'species', 'selection_status', 'source_header', 'length_aa'])
    writer.writerows(audit)
if not any(species == 'Cerastes' and symbol == 'ADAM' for species, symbol, _, _ in found) and not any(symbol == 'ADAM28' for _, symbol, _, _ in found):
    raise RuntimeError('No uniquely recoverable ADAM reference protein; see anchor_reference_protein_audit.tsv and resolve the reference ID')
print('Targeted proteins:', ', '.join(f'{s}_{g}' for s, g, _, _ in found))
PY

# 8b. all-species reference-anchor BLASTp, including Gloydius.
# This detects candidate anchors even when Gloydius query proteins were missing.
printf '\n===== REFERENCE ANCHOR BLASTP IN ALL FOUR SPECIES =====\n'
for FAMILY in SVMP PLA2; do
    QUERY="$EVIDENCE/${FAMILY}_anchor_reference_proteins.faa"
    for SPECIES in Gloydius_ussuriensis Cerastes_gasperettii Crotalus_adamanteus Naja_naja; do
        OUTPUT="$EVIDENCE/${FAMILY}_reference_anchors_vs_${SPECIES}.blastp.tsv"
        if [[ ! -s "$QUERY" ]]; then
            : > "$OUTPUT"
            continue
        fi
        if [[ "$SPECIES" == 'Gloydius_ussuriensis' ]]; then DB="$GUSS_DB"; else DB="$BLASTDIR/${SPECIES}_proteins_db"; fi
        blastp -query "$QUERY" -db "$DB" -evalue 1e-5 -max_target_seqs 10 -num_threads 2 -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' -out "$OUTPUT"
        echo "$FAMILY $SPECIES: $(wc -l < "$OUTPUT") anchor BLASTp HSPs"
    done
done

# 9. derive Gloydius search intervals from flanking coordinates obtained in step 2
printf '\n===== EXTRACT GLOYDIUS GENOMIC SEARCH INTERVALS =====\n'
python - "$FLANKING" "${GUSS_GENOME}.fai" "$EVIDENCE/Guss_target_intervals.tsv" <<'PY'
import csv, sys
from collections import defaultdict
flanking, fai, output = sys.argv[1:]
lengths = {}
with open(fai) as handle:
    for line in handle:
        parts = line.split('\t')
        lengths[parts[0]] = int(parts[1])
with open(flanking) as handle: rows = list(csv.DictReader(handle, delimiter='\t'))
regions = []
for fam in ('SVMP', 'PLA2'):
    subset = [r for r in rows if r['family'] == fam and r['chr'] in lengths]
    if not subset: raise RuntimeError(f'Missing genomic coordinates in flanking inventory for {fam}')
    chrs = set(r['chr'] for r in subset)
    if len(chrs) != 1: raise RuntimeError(f'{fam} spans multiple chromosomes in flanking inventory: {chrs}')
    chrom = next(iter(chrs))
    start = max(1, min(int(r['start']) for r in subset) - 20000)
    end = min(lengths[chrom], max(int(r['end']) for r in subset) + 20000)
    regions.append(dict(family=fam, chromosome=chrom, start=start, end=end, search_region=f'{chrom}:{start}-{end}'))
with open(output, 'w', newline='') as handle:
    w = csv.DictWriter(handle, fieldnames=['family', 'chromosome', 'start', 'end', 'search_region'], delimiter='\t')
    w.writeheader()
    w.writerows(regions)
for r in regions: print(r['family'], r['search_region'])
PY

while IFS=$'\t' read -r FAMILY CHR START END REGION; do
    if [[ "$FAMILY" == 'family' ]]; then continue; fi
    GENOME_REGION="$EVIDENCE/Guss_${FAMILY}_search_region.fa"
    samtools faidx "$GUSS_GENOME" "$REGION" > "$GENOME_REGION"
    OUTPUT="$EVIDENCE/Guss_${FAMILY}_anchors.tblastn.tsv"
    if [[ ! -s "$EVIDENCE/${FAMILY}_anchor_reference_proteins.faa" ]]; then : > "$OUTPUT"; continue; fi
    tblastn -query "$EVIDENCE/${FAMILY}_anchor_reference_proteins.faa" -subject "$GENOME_REGION" -seg yes -evalue 1e-5 -num_threads 2 -outfmt '6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore' -out "$OUTPUT"
    echo "$FAMILY: $(wc -l < "$OUTPUT") TBLASTN HSPs"
done < "$EVIDENCE/Guss_target_intervals.tsv"

# 10. consolidate BLASTp reciprocal evidence, GFF coordinates and TBLASTN hits
printf '\n===== GENERATE REVIEW TABLES =====\n'
python - "$BLASTDIR" "$ASSEMBLIES" "$FLANKING" "$EVIDENCE" <<'PY'
import csv, re, sys
from pathlib import Path
from collections import defaultdict
from urllib.parse import unquote
blastdir, assemblies, flanking, evidence = [Path(x) for x in sys.argv[1:]]
def read_tsv(path):
    with path.open() as handle: return list(csv.DictReader(handle, delimiter='\t'))
def write_tsv(path, rows, cols):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)

# Subject protein accession -> annotated GFF locus; records may only map to a CDS.
def build_gff_index(path):
    features, parent, alias = {}, {}, defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith('#'): continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9: continue
            a = {k: unquote(v) for part in f[8].split(';') if '=' in part for k, v in [part.split('=', 1)]}
            fid = a.get('ID', '') or f'_record_{len(features)}'
            features[fid] = (f, a)
            parent[fid] = a.get('Parent', '').split(',') if a.get('Parent') else []
            for k in ('ID', 'Name', 'protein_id', 'transcript_id', 'gene_id'):
                for v in a.get(k, '').split(','):
                    if v: alias[v].append(fid)
    return features, parent, alias

def map_gene(subject, index):
    features, parent, alias = index
    matches = alias.get(subject, [])
    if not matches: return ('', '', '', '', '', 'unmapped_GFF_subject')
    resolved = []
    for node in matches:
        seen = set()
        todo = [node]
        while todo:
            cur = todo.pop()
            if cur in seen or cur not in features: continue
            seen.add(cur)
            f, a = features[cur]
            if f[2] == 'gene': resolved.append((cur, f))
            todo.extend(parent.get(cur, []))
    if resolved:
        gid, f = sorted(set((g, tuple(f)) for g, f in resolved), key=lambda x: x[0])[0]
        return (gid, f[0], f[3], f[4], f[6], 'gene_model_mapped')
    f, a = features[matches[0]]
    return (matches[0], f[0], f[3], f[4], f[6], 'feature_only_review')

index = {species: build_gff_index(assemblies / species / f'{species}.annotation.gff3') for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja')}
flanks = {r['gene_id']: r for r in read_tsv(flanking)}
known = {'FUN_014513':'STC1', 'FUN_014526':'NEFM', 'FUN_014527':'NEFL', 'FUN_018057':'MUL1', 'FUN_018058':'PLA2G2F', 'FUN_018073':'OTUD3'}

forward = []
for family in ('SVMP', 'PLA2'):
    for species in index:
        reciprocal = defaultdict(list)
        rpath = evidence / f'{family}_{species}_reciprocal.blastp.tsv'
        with rpath.open() as handle:
            for r in csv.reader(handle, delimiter='\t'):
                if len(r) >= 9: reciprocal[r[0]].append(r)
        fpath = blastdir / family / f'Guss_{family}_vs_{species}.blastp.tsv'
        with fpath.open() as handle:
            rows = [r for r in csv.reader(handle, delimiter='\t') if len(r) >= 9]
        rows.sort(key=lambda r: (r[0], -float(r[8])))
        ranks = defaultdict(int)
        for r in rows:
            ranks[r[0]] += 1
            if ranks[r[0]] > 20: continue
            query, subject = r[:2]
            source_gene = re.sub(r'-T\d+$', '', query)
            rev_hits = sorted(reciprocal.get(subject, []), key=lambda x: -float(x[8]))
            reciprocal_best = rev_hits[0][1] if rev_hits else ''
            reciprocal_rank = next((i+1 for i, x in enumerate(rev_hits) if re.sub(r'-T\d+$', '', x[1]) == source_gene), '')
            geneid, chrom, start, end, strand, mapping = map_gene(subject, index[species])
            forward.append(dict(family=family, query_gene_id=source_gene, anchor_hint=known.get(source_gene, ''), query_transcript=query, reference_species=species, subject_protein_id=subject, forward_rank=ranks[query], percent_identity=r[2], aligned_aa=r[3], query_length_aa=r[4], subject_length_aa=r[5], query_coverage_pct=r[6], evalue=r[7], bitscore=r[8], reciprocal_best=reciprocal_best, reciprocal_rank=reciprocal_rank, gff_gene_id=geneid, chr=chrom, start=start, end=end, strand=strand, gff_mapping_status=mapping, validation_status='candidate_requires_review'))

cols = ['family', 'query_gene_id', 'anchor_hint', 'query_transcript', 'reference_species', 'subject_protein_id', 'forward_rank', 'percent_identity', 'aligned_aa', 'query_length_aa', 'subject_length_aa', 'query_coverage_pct', 'evalue', 'bitscore', 'reciprocal_best', 'reciprocal_rank', 'gff_gene_id', 'chr', 'start', 'end', 'strand', 'gff_mapping_status', 'validation_status']
write_tsv(evidence / 'flanking_blastp_reciprocal_gff_evidence.tsv', forward, cols)

# Reference-anchored BLASTp across all four species. These are candidates, not called orthologs.
reference_candidates = []
for family in ('SVMP', 'PLA2'):
    for species in ('Gloydius_ussuriensis', 'Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        path = evidence / f'{family}_reference_anchors_vs_{species}.blastp.tsv'
        hits = defaultdict(list)
        with path.open() as handle:
            for r in csv.reader(handle, delimiter='\t'):
                if len(r) >= 9: hits[r[0]].append(r)
        for query, matches in sorted(hits.items()):
            for rank, r in enumerate(sorted(matches, key=lambda x: -float(x[8]))[:10], 1):
                subject = r[1]
                if species == 'Gloydius_ussuriensis':
                    gid = re.sub(r'-T\d+$', '', subject)
                    g = flanks.get(gid, {})
                    chrom, start, end, strand = (g.get('chr', ''), g.get('start', ''), g.get('end', ''), g.get('strand', ''))
                    mapping = 'Step2_gene_model' if g else 'outside_Step2_inventory'
                else:
                    gid, chrom, start, end, strand, mapping = map_gene(subject, index[species])
                reference_candidates.append(dict(family=family, anchor_query=query, species=species, subject_protein_id=subject, hit_rank=rank, identity_pct=r[2], aligned_aa=r[3], query_length_aa=r[4], subject_length_aa=r[5], query_coverage_pct=r[6], evalue=r[7], bitscore=r[8], gene_id=gid, chr=chrom, start=start, end=end, strand=strand, mapping_status=mapping, validation_status='candidate_requires_review'))

write_tsv(evidence / 'anchor_reference_blastp_candidates.tsv', reference_candidates, ['family', 'anchor_query', 'species', 'subject_protein_id', 'hit_rank', 'identity_pct', 'aligned_aa', 'query_length_aa', 'subject_length_aa', 'query_coverage_pct', 'evalue', 'bitscore', 'gene_id', 'chr', 'start', 'end', 'strand', 'mapping_status', 'validation_status'])

# Genomic TBLASTN source sequences are extracted in reference orientation; convert local to 1-based genome coordinates.
regions = {r['family']: r for r in read_tsv(evidence / 'Guss_target_intervals.tsv')}
all_genome_hits = []
for family, region in regions.items():
    offset = int(region['start']) - 1
    path = evidence / f'Guss_{family}_anchors.tblastn.tsv'
    with path.open() as handle:
        for r in csv.reader(handle, delimiter='\t'):
            if len(r) != 11: continue
            gs, ge = sorted((int(r[7]) + offset, int(r[8]) + offset))
            overlapping = sorted(g['gene_id'] for g in flanks.values() if g['family'] == family and g['chr'] == region['chromosome'] and int(g['start']) <= ge and int(g['end']) >= gs)
            all_genome_hits.append(dict(family=family, query_anchor=r[0], chr=region['chromosome'], start=gs, end=ge, strand='+' if int(r[7]) <= int(r[8]) else '-', query_start_aa=r[5], query_end_aa=r[6], query_length_aa=r[4], aligned_aa=r[3], percent_identity=r[2], evalue=r[9], bitscore=r[10], overlapping_FUN_gene_ids=';'.join(overlapping), validation_status='tblastn_HSP_not_gene_assignment'))

all_genome_hits.sort(key=lambda x: (x['family'], x['query_anchor'], -float(x['bitscore'])))
write_tsv(evidence / 'Guss_anchor_tblastn_genomic_hits.tsv', all_genome_hits, ['family', 'query_anchor', 'chr', 'start', 'end', 'strand', 'query_start_aa', 'query_end_aa', 'query_length_aa', 'aligned_aa', 'percent_identity', 'evalue', 'bitscore', 'overlapping_FUN_gene_ids', 'validation_status'])

# Anchor validation checklist: no automatic orthology claims based on one BLAST HSP.
expected = {'SVMP': ['STC1', 'NEFM', 'NEFL', 'ADAM'], 'PLA2': ['OTUD3', 'MUL1', 'VWA5B1', 'PLA2G2F']}
review = []
for family, symbols in expected.items():
    for symbol in symbols:
        candidates = sorted(gid for gid, sym in known.items() if sym == symbol and gid in flanks and flanks[gid]['family'] == family)
        queries = sorted(set(r['query_anchor'] for r in all_genome_hits if r['family'] == family and r['query_anchor'].endswith('_' + symbol)))
        hits = [r for r in all_genome_hits if r['family'] == family and r['query_anchor'] in queries]
        for species in ('Gloydius_ussuriensis', 'Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
            p = [r for r in reference_candidates if r['species'] == species and r['family'] == family and r['anchor_query'].endswith('_' + symbol)]
            p.sort(key=lambda r: -float(r['bitscore']))
            best = p[0] if p else {}
            review.append(dict(family=family, species=species, anchor_symbol=symbol, Gloydius_model_hint=';'.join(candidates) if species == 'Gloydius_ussuriensis' else '', top_protein_id=best.get('subject_protein_id', ''), top_gene_id=best.get('gene_id', ''), top_blastp_identity=best.get('identity_pct', ''), top_blastp_query_coverage=best.get('query_coverage_pct', ''), reference_tblastn_queries=';'.join(queries) if species == 'Gloydius_ussuriensis' else '', tblastn_hsp_count=len(hits) if species == 'Gloydius_ussuriensis' else '', best_tblastn_bitscore=max((float(r['bitscore']) for r in hits), default=0) if species == 'Gloydius_ussuriensis' else '', validation_status='manual_homology_and_synteny_review_required'))

write_tsv(evidence / 'flanking_anchor_review_checklist.tsv', review, ['family', 'species', 'anchor_symbol', 'Gloydius_model_hint', 'top_protein_id', 'top_gene_id', 'top_blastp_identity', 'top_blastp_query_coverage', 'reference_tblastn_queries', 'tblastn_hsp_count', 'best_tblastn_bitscore', 'validation_status'])
print('Forward BLASTp + reciprocal + GFF evidence rows:', len(forward))
print('Targeted genomic TBLASTN HSPs:', len(all_genome_hits))
print('Reference anchored BLASTp candidate HSPs:', len(reference_candidates))
print('Checklist entries:', len(review))
print('NOTE: no ortholog identities have been automatically validated.')
PY

# 11. best-hit report
printf '\n===== BEST FORWARD BLASTP HITS =====\n'
python - "$BLASTDIR" <<'PY'
import csv, sys
from pathlib import Path
from collections import defaultdict
blastdir = Path(sys.argv[1])
for family in ('SVMP', 'PLA2'):
    print(f'\n===== {family} =====')
    for species in ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja'):
        path = blastdir / family / f'Guss_{family}_vs_{species}.blastp.tsv'
        hits = defaultdict(list)
        with path.open() as handle:
            for row in csv.reader(handle, delimiter='\t'):
                if len(row) >= 9: hits[row[0]].append(row)
        print(f'\n--- {species} ---')
        if not hits: print('No hits passing threshold.')
        for query in sorted(hits):
            best = max(hits[query], key=lambda r: float(r[8]))
            print(f'{query:<20} {best[1]:<28} ID={float(best[2]):5.1f}% Qcov={float(best[6]):5.1f}% E={best[7]}')
print('\n===== COMPLETE =====')
PY

echo
printf 'Additional Step 4 evidence: %s\n' "$EVIDENCE"
printf 'Step 6/7 input for review: %s\n' "$EVIDENCE/flanking_blastp_reciprocal_gff_evidence.tsv"
printf 'Genomic candidates: %s\n' "$EVIDENCE/Guss_anchor_tblastn_genomic_hits.tsv"
printf 'Anchor checklist: %s\n' "$EVIDENCE/flanking_anchor_review_checklist.tsv"
```

The most important outputs from this script are:
  - `flanking_blastp_reciprocal_gff_evidence.tsv`: Forward BLASTp, reciprocal BLASTp, and comparative genomic coordinates; needed for homolog identification in step 6.
  - `anchor_reference_blastp_candidates.tsv`: Candidate anchor proteins identified across all four species; needed for anchor comparison in steps 6 and 7.
  - `Guss_SVMP_ADAM28_candidate_loci.tsv`: Potential ADAM28-like loci, genomic coordinates, coverage, and Funannotate gene-model overlap; needed for missing ADAM investigation in step 6. I want to clarify that we check for this in step 6, because when I ran through this workflow for the first time the microsynteny plot didn't show ADAM for *G. ussuriensis*. We are circling back here to revise the workflow. 
  - `PLA2_model_reconciliation.tsv`: Comparison of potentially fragmented or duplicated Gloydius PLA₂-region gene models; needed for model reconciliation in step 6.
  - `PLA2_anchor_presence_matrix.tsv`: Candidate PLA₂ anchors and their coordinates across the four species; needed for microsynteny annotation in step 6.
  - `flanking_anchor_review_checklist.tsv`: Summary of SVMP and PLA₂ anchor evidence and unresolved assignments	needed for final review in steps 6 and 7.

Inspect the output once the script finishes running:
```sh
# ============================================================
# inspect SVMP and PLA2 homology evidence
# ============================================================

conda activate funannotate

ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
EVIDENCE="$ROOT/synteny/toxin_microsynteny/02_blast/flanking_homology_evidence"

python - "$EVIDENCE" <<'PY'

import csv
import sys
from pathlib import Path

root = Path(sys.argv[1])

def read_file(name):
    path = root / name
    if not path.exists():
        print(f"FILE MISSING: {path}")
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

# ============================================================
# 1. SVMP ADAM28 candidate loci
# ============================================================

print("\n===== 1. ADAM28 CANDIDATE LOCI =====")

rows = read_file("Guss_SVMP_ADAM28_candidate_loci.tsv")

rows = sorted(
    rows,
    key=lambda r: float(r["unique_query_coverage_pct"] or 0),
    reverse=True
)

for r in rows[:15]:
    print(
        f"\n{r['chr']}:{r['start']}-{r['end']} "
        f"({r['strand']})"
        f"\n  Reference: {r['query_anchor']}"
        f"\n  HSPs: {r['hsp_count']}"
        f"\n  Query coverage: {r['unique_query_coverage_pct']}%"
        f"\n  Collinear: {r['collinear_query_order']}"
        f"\n  GFF status: {r['full_GFF_model_status']}"
        f"\n  Overlapping genes: {r['overlapping_full_GFF_gene_ids']}"
        f"\n  Nearest gene: {r['nearest_full_GFF_gene_id']}"
        f"\n  Distance: {r['nearest_full_GFF_distance_bp']} bp"
    )

# ============================================================
# 2. PLA2 reference anchor recovery
# ============================================================

print("\n===== 2. PLA2 REFERENCE ANCHORS =====")

rows = read_file("PLA2_anchor_reference_audit.tsv")

for r in rows:
    print(
        f"\n{r['anchor_symbol']}"
        f"\n  Protein: {r['source_protein_id']}"
        f"\n  Selection: {r['selection_status']}"
        f"\n  Chromosome: {r['source_seq_id']}"
        f"\n  GFF chromosome: {r['comparative_GFF_chr']}"
        f"\n  Length: {r['length_aa']} aa"
        f"\n  Context: {r['genomic_context_status']}"
    )

# ============================================================
# 3. PLA2 anchor presence across four species
# ============================================================

print("\n===== 3. PLA2 ANCHOR PRESENCE =====")

rows = read_file("PLA2_anchor_presence_matrix.tsv")

for r in rows:
    print(
        f"{r['anchor_symbol']:<9} "
        f"{r['species']:<23} "
        f"{r['top_subject_protein']:<25} "
        f"Chr={r['top_subject_chr']:<16} "
        f"Status={r['evidence_status']}"
    )

# ============================================================
# 4. PLA2 gene-model reconciliation
# ============================================================

print("\n===== 4. PLA2 MODEL RECONCILIATION =====")

rows = read_file("PLA2_model_reconciliation.tsv")

for r in rows:
    print(
        f"\n{r['species']}: "
        f"{r['Gloydius_gene_a']} vs {r['Gloydius_gene_b']}"
        f"\n  Genomic gap: {r['Gloydius_genomic_gap_bp']} bp"
        f"\n  Genomic overlap: {r['Gloydius_gene_overlap_bp']} bp"
        f"\n  Shared subject: {r['shared_subject']}"
        f"\n  Subject AA overlap: {r['shared_subject_aa_overlap']}"
        f"\n  Overlap of shorter segment: "
        f"{r['overlap_pct_of_shorter_subject_segment']}%"
        f"\n  Evidence: {r['evidence_flag']}"
    )

# ============================================================
# 5. Overall anchor review
# ============================================================

print("\n===== 5. ANCHOR REVIEW CHECKLIST =====")

rows = read_file("flanking_anchor_review_checklist.tsv")

for r in rows:
    print(
        f"{r['family']:<5} "
        f"{r['anchor_symbol']:<9} "
        f"{r['species']:<23} "
        f"Top={r['top_protein_id']:<25} "
        f"Gene={r['top_gene_id']}"
    )

print("\n===== COMPLETE =====")

PY
```

The final output line of this script should look like:
```sh
===== 5. ANCHOR REVIEW CHECKLIST =====
SVMP  STC1      Gloydius_ussuriensis    Top=FUN_014513-T1             Gene=FUN_014513
SVMP  STC1      Cerastes_gasperettii    Top=crotalus-ruber-cadam-rna-Stc1-gene14556-R0-R0 Gene=gene-24098
SVMP  STC1      Crotalus_adamanteus     Top=KAK9395204.1              Gene=gene-NXF25_014550
SVMP  STC1      Naja_naja               Top=KAG8141316.1              Gene=gene-E2320_006948
SVMP  NEFM      Gloydius_ussuriensis    Top=FUN_014526-T1             Gene=FUN_014526
SVMP  NEFM      Cerastes_gasperettii    Top=rna-NEFM                  Gene=NEFM
SVMP  NEFM      Crotalus_adamanteus     Top=KAK9395179.1              Gene=gene-NXF25_014525
SVMP  NEFM      Naja_naja               Top=KAG8141328.1              Gene=gene-E2320_006958
SVMP  NEFL      Gloydius_ussuriensis    Top=FUN_014527-T1             Gene=FUN_014527
SVMP  NEFL      Cerastes_gasperettii    Top=rna-NEFL                  Gene=NEFL
SVMP  NEFL      Crotalus_adamanteus     Top=KAK9395178.1              Gene=gene-NXF25_014524
SVMP  NEFL      Naja_naja               Top=KAG8141329.1              Gene=gene-E2320_006959
SVMP  ADAM      Gloydius_ussuriensis    Top=FUN_023349-T1             Gene=FUN_023349
SVMP  ADAM      Cerastes_gasperettii    Top=rna-ADAM                  Gene=ADAM
SVMP  ADAM      Crotalus_adamanteus     Top=KAK9395203.1              Gene=gene-NXF25_014549
SVMP  ADAM      Naja_naja               Top=KAG8141317.1              Gene=gene-E2320_006950
PLA2  OTUD3     Gloydius_ussuriensis    Top=FUN_018073-T1             Gene=FUN_018073
PLA2  OTUD3     Cerastes_gasperettii    Top=crotalus-ruber-Otud3-R0   Gene=gene-22973
PLA2  OTUD3     Crotalus_adamanteus     Top=KAK9392176.1              Gene=gene-NXF25_017763
PLA2  OTUD3     Naja_naja               Top=KAG8142454.1              Gene=gene-E2320_006360
PLA2  MUL1      Gloydius_ussuriensis    Top=FUN_018057-T1             Gene=FUN_018057
PLA2  MUL1      Cerastes_gasperettii    Top=crotalus-ruber-Mul1-R0    Gene=gene-22976
PLA2  MUL1      Crotalus_adamanteus     Top=KAK9392180.1              Gene=gene-NXF25_017767
PLA2  MUL1      Naja_naja               Top=KAG8142459.1              Gene=gene-E2320_006356
PLA2  VWA5B1    Gloydius_ussuriensis    Top=FUN_018055-T1             Gene=FUN_018055
PLA2  VWA5B1    Cerastes_gasperettii    Top=crotalus-ruber-cadam-rna-Vwa5b1-gene17777-R0-R0 Gene=gene-23185
PLA2  VWA5B1    Crotalus_adamanteus     Top=KAK9392182.1              Gene=gene-NXF25_017769
PLA2  VWA5B1    Naja_naja               Top=KAG8142461.1              Gene=gene-E2320_006373
PLA2  PLA2G2F   Gloydius_ussuriensis    Top=FUN_018058-T1             Gene=FUN_018058
PLA2  PLA2G2F   Cerastes_gasperettii    Top=rna-Pla2g2F               Gene=Pla2g2F
PLA2  PLA2G2F   Crotalus_adamanteus     Top=KAK9392178.1              Gene=gene-NXF25_017765
PLA2  PLA2G2F   Naja_naja               Top=KAG8142458.1              Gene=gene-E2320_006357
```
The results suggest that we have recovered the PLA₂ flanking anchors at the protein-annotation level, identified a strong unannotated ADAM28-like genomic candidate, and obtained evidence that three Gloydius VWA5B1-associated gene models may represent a fragmented annotation.

#### Step 5: Extract anchor-bounded SVMP and PLA₂ genomic intervals
The script below will extract genomic intervals spanning the conserved anchors, with an additional 20 kb on each side. This will produce eight FASTA sequences: one SVMP and one PLA₂ region for each species.
```sh
# ============================================================
# extract anchor-bounded SVMP and PLA2 genomic intervals
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"
OUTDIR="$SYNTENY/03_genomic_intervals"
INVENTORY="$SYNTENY/01_reference_inventory"
EVIDENCE="$SYNTENY/02_blast/flanking_homology_evidence"
GUSS_GFF="$ROOT/annotation/funannotate/update_results/Gloydius_ussuriensis_AMNH_21010.gff3"
TOXIN="$ROOT/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/toxin_annotation.gtf"

mkdir -p "$OUTDIR/SVMP" "$OUTDIR/PLA2" "$INVENTORY"

for TOOL in python samtools; do
    command -v "$TOOL" >/dev/null || { echo "ERROR: $TOOL not found" >&2; exit 1; }
done

# 2. reference genomes
GUSS="$ROOT/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
CGAS="$ASSEMBLIES/Cerastes_gasperettii/Cerastes_gasperettii.genome.fa"
CADAM="$ASSEMBLIES/Crotalus_adamanteus/Crotalus_adamanteus.genome.fa"
NAJA="$ASSEMBLIES/Naja_naja/Naja_naja.genome.fa"

for FILE in "$GUSS" "$CGAS" "$CADAM" "$NAJA" "$GUSS_GFF" "$TOXIN" \
    "$EVIDENCE/anchor_reference_blastp_candidates.tsv" \
    "$EVIDENCE/flanking_anchor_review_checklist.tsv" \
    "$EVIDENCE/Guss_SVMP_ADAM28_candidate_loci.tsv"; do

    [[ -s "$FILE" ]] || { echo "ERROR: missing or empty input: $FILE" >&2; exit 1; }

done

for GENOME in "$GUSS" "$CGAS" "$CADAM" "$NAJA"; do
    [[ -s "${GENOME}.fai" ]] || samtools faidx "$GENOME"
done

# 3. extract anchor-bounded regions
echo '===== EXTRACT ANCHOR-BOUNDED REGIONS ====='

MANIFEST="$INVENTORY/anchor_bounded_regions.tsv"
ANCHORS="$INVENTORY/selected_anchor_coordinates.tsv"
TOXIN_CHECK="$INVENTORY/Guss_ToxCodAn_region_containment.tsv"
PLAN="$INVENTORY/anchor_interval_extraction_plan.tsv"

# use Step 4 GFF-mapped anchors and TBLASTN candidate loci as coordinate evidence
# expand Gloydius intervals to include retained ToxCodAn toxin loci on that chromosome
# all intervals receive 20 kb padding; ambiguous anchors remain explicitly provisional
python - "$EVIDENCE" "$GUSS_GFF" "$TOXIN" "$MANIFEST" "$ANCHORS" "$TOXIN_CHECK" "$PLAN" \
    "$GUSS" "$CGAS" "$CADAM" "$NAJA" "$OUTDIR" <<'PY'

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

evidence, guss_gff, toxin, manifest, anchor_file, toxin_file, plan = map(Path, sys.argv[1:8])

species_order = (
    'Gloydius_ussuriensis',
    'Cerastes_gasperettii',
    'Crotalus_adamanteus',
    'Naja_naja'
)

fastas = dict(zip(species_order, map(Path, sys.argv[8:12])))
outdir = Path(sys.argv[12])

family_order = ('SVMP', 'PLA2')

expected = {
    'SVMP': ('STC1', 'NEFM', 'NEFL', 'ADAM'),
    'PLA2': ('OTUD3', 'PLA2G2F', 'MUL1', 'VWA5B1')
}

anchor_padding = 20000
excluded_loci = {'SVMP-6', 'PLA2-1'}


def read_tsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write_tsv(path, rows, fields):
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter='\t',
            lineterminator='\n',
            extrasaction='ignore'
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_gff_attributes(text):
    return {
        k: unquote(v)
        for part in text.split(';')
        if '=' in part
        for k, v in [part.split('=', 1)]
    }


def gff_gene_coordinates(path, wanted):
    found = {}

    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue

            f = line.rstrip('\n').split('\t')

            if len(f) < 9 or f[2] != 'gene':
                continue

            gid = parse_gff_attributes(f[8]).get('ID', '')

            if gid in wanted:
                found[gid] = (
                    f[0],
                    int(f[3]),
                    int(f[4]),
                    f[6]
                )

    return found


# read genome chromosome lengths
lengths = {}

for species, fasta in fastas.items():
    lengths[species] = {}

    with Path(str(fasta) + '.fai').open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')

            if len(fields) >= 2:
                lengths[species][fields[0]] = int(fields[1])


# read step 4 anchor review and BLASTp candidate tables
reviews = {
    (r['family'], r['species'], r['anchor_symbol']): r
    for r in read_tsv(evidence / 'flanking_anchor_review_checklist.tsv')
}

candidates = defaultdict(list)

for r in read_tsv(evidence / 'anchor_reference_blastp_candidates.tsv'):
    candidates[
        (r['family'], r['species'], r['subject_protein_id'])
    ].append(r)


# recover the entire candidate Gloydius VWA5B1 locus
vwa_ids = ('FUN_018053', 'FUN_018054', 'FUN_018055')

vwa_models = gff_gene_coordinates(
    guss_gff,
    set(vwa_ids)
)

if set(vwa_models) != set(vwa_ids):
    raise RuntimeError(
        'Missing Gloydius VWA5B1 component models in GFF: ' +
        str(set(vwa_ids) - set(vwa_models))
    )

if (
    len({vwa_models[g][0] for g in vwa_ids}) != 1
    or len({vwa_models[g][3] for g in vwa_ids}) != 1
):
    raise RuntimeError(
        'VWA5B1 component models disagree on chromosome or strand'
    )


# build candidate anchor coordinates using Step 4 evidence
anchors = []

for species in species_order:

    for family in family_order:

        preliminary = []

        for symbol in expected[family]:

            review = reviews.get((family, species, symbol))

            if review is None:
                raise RuntimeError(
                    f'Missing Step 4 anchor checklist entry: '
                    f'{species} {family} {symbol}'
                )

            if (
                species == 'Gloydius_ussuriensis'
                and family == 'SVMP'
                and symbol == 'ADAM'
            ):
                continue

            if (
                species == 'Gloydius_ussuriensis'
                and family == 'PLA2'
                and symbol == 'VWA5B1'
            ):

                chrom = vwa_models[vwa_ids[0]][0]

                preliminary.append(dict(
                    species=species,
                    family=family,
                    anchor_symbol=symbol,
                    chromosome=chrom,
                    start=min(vwa_models[g][1] for g in vwa_ids),
                    end=max(vwa_models[g][2] for g in vwa_ids),
                    strand=vwa_models[vwa_ids[0]][3],
                    gene_ids=';'.join(vwa_ids),
                    protein_ids=';'.join(g + '-T1' for g in vwa_ids),
                    evidence_status='provisional_split_VWA5B1_three_models_not_merged'
                ))

                continue

            protein = review['top_protein_id']

            if not protein:
                raise RuntimeError(
                    f'No Step 4 candidate protein: '
                    f'{species} {family} {symbol}'
                )

            matches = [
                r
                for r in candidates.get((family, species, protein), [])
                if r['chr'] and r['start'] and r['end']
            ]

            if review['top_gene_id']:
                matches = [
                    r for r in matches
                    if r['gene_id'] == review['top_gene_id']
                ]

            if not matches:
                raise RuntimeError(
                    f'No mapped GFF candidate: '
                    f'{species} {family} {symbol} ({protein})'
                )

            matches.sort(
                key=lambda r: float(r['bitscore']),
                reverse=True
            )

            r = matches[0]

            preliminary.append(dict(
                species=species,
                family=family,
                anchor_symbol=symbol,
                chromosome=r['chr'],
                start=int(r['start']),
                end=int(r['end']),
                strand=r['strand'],
                gene_ids=r['gene_id'],
                protein_ids=protein,
                evidence_status='protein_homology_and_GFF_candidate_requires_review'
            ))

        # require the non-ADAM flanking anchors to share a chromosome
        counts = Counter(
            a['chromosome']
            for a in preliminary
            if a['anchor_symbol'] != 'ADAM'
        )

        if len(counts) != 1:
            raise RuntimeError(
                f'Anchor chromosome disagreement: '
                f'{species} {family}: {dict(counts)}'
            )

        chrom = next(iter(counts))

        if chrom not in lengths[species]:
            raise RuntimeError(
                f'{species} chromosome absent from FASTA index: {chrom}'
            )

        # in Gloydius, use the collinear genomic ADAM28-like candidate
        # do not assign ADAM28 based only on its top SVMP-like protein hit
        if (
            species == 'Gloydius_ussuriensis'
            and family == 'SVMP'
        ):

            adam_candidates = [
                r
                for r in read_tsv(
                    evidence / 'Guss_SVMP_ADAM28_candidate_loci.tsv'
                )
                if r['chr'] == chrom
                and r['collinear_query_order'] == 'yes'
                and float(r['unique_query_coverage_pct'] or 0) >= 50
            ]

            if not adam_candidates:
                raise RuntimeError(
                    f'No collinear ADAM28-like genomic candidate on {chrom}'
                )

            adam_candidates.sort(
                key=lambda r: (
                    float(r['unique_query_coverage_pct']),
                    float(r['sum_hsp_bitscore'])
                ),
                reverse=True
            )

            r = adam_candidates[0]

            preliminary.append(dict(
                species=species,
                family=family,
                anchor_symbol='ADAM',
                chromosome=chrom,
                start=int(r['start']),
                end=int(r['end']),
                strand=r['strand'],
                gene_ids='unannotated_ADAM28_like_candidate',
                protein_ids=r['query_anchor'],
                evidence_status='provisional_TBLASTN_candidate_not_confirmed_ortholog'
            ))

        if any(a['chromosome'] != chrom for a in preliminary):
            raise RuntimeError(
                f'Potential nonsyntenic anchor: {species} {family}: ' +
                str([
                    (a['anchor_symbol'], a['chromosome'])
                    for a in preliminary
                ])
            )

        if len(preliminary) != len(expected[family]):
            raise RuntimeError(
                f'Incomplete anchor set for {species} {family}: '
                f'{len(preliminary)}/{len(expected[family])}'
            )

        anchors.extend(preliminary)


# extract Gloydius ToxCodAn CDS locus extents
gtf_loci = {}

with toxin.open() as handle:

    for line in handle:

        if not line.strip() or line.startswith('#'):
            continue

        f = line.rstrip('\n').split('\t')

        if len(f) < 9 or f[2] != 'CDS':
            continue

        match = re.search(
            r'\bgene_id\s+"([^"]+)"',
            f[8]
        )

        if not match:
            continue

        gene = match.group(1)

        family = next(
            (
                fam
                for fam in family_order
                if gene.startswith(fam + '-')
            ),
            None
        )

        if family is None:
            continue

        key = (family, gene, f[0])
        start, end = int(f[3]), int(f[4])

        if key not in gtf_loci:
            gtf_loci[key] = [start, end]
        else:
            gtf_loci[key][0] = min(
                gtf_loci[key][0],
                start
            )

            gtf_loci[key][1] = max(
                gtf_loci[key][1],
                end
            )

if not gtf_loci:
    raise RuntimeError(
        'No SVMP or PLA2 CDS loci found in the ToxCodAn GTF'
    )


# establish genomic intervals using anchors and retained toxin loci
manifest_rows = []
plan_rows = []

for species in species_order:

    for family in family_order:

        selected = [
            a for a in anchors
            if a['species'] == species
            and a['family'] == family
        ]

        chrom = selected[0]['chromosome']

        start = min(a['start'] for a in selected)
        end = max(a['end'] for a in selected)

        toxin_count = 0

        if species == 'Gloydius_ussuriensis':

            retained = [
                (g, lo, hi)
                for (fam, g, c), (lo, hi) in gtf_loci.items()
                if fam == family
                and c == chrom
                and g not in excluded_loci
            ]

            for gene, lo, hi in retained:
                start = min(start, lo)
                end = max(end, hi)

            toxin_count = len(retained)

            if not retained:
                raise RuntimeError(
                    f'No retained {family} ToxCodAn loci found '
                    f'on Gloydius anchor chromosome {chrom}'
                )

        start = max(1, start - anchor_padding)
        end = min(lengths[species][chrom], end + anchor_padding)

        if start >= end:
            raise RuntimeError(
                f'Invalid interval for {species} {family}: '
                f'{chrom}:{start}-{end}'
            )

        if any(
            not (start <= a['start'] <= a['end'] <= end)
            for a in selected
        ):
            raise RuntimeError(
                f'Anchor containment failed for {species} {family}'
            )

        dest = (
            outdir / family /
            f'{species}_{family}_anchor_bounded.fa'
        )

        status = (
            'all_selected_anchors_and_retained_toxin_loci_contained'
            if species == 'Gloydius_ussuriensis'
            else 'all_selected_anchors_contained;comparative_toxins_not_yet_independently_validated'
        )

        row = dict(
            species=species,
            family=family,
            chromosome=chrom,
            start=start,
            end=end,
            length_bp=end-start+1,
            fasta=str(dest),
            anchor_count=len(selected),
            anchor_symbols=';'.join(
                a['anchor_symbol']
                for a in selected
            ),
            anchor_coordinate_order=';'.join(
                a['anchor_symbol']
                for a in sorted(
                    selected,
                    key=lambda a: a['start']
                )
            ),
            Gloydius_retained_toxin_loci=(
                toxin_count
                if species == 'Gloydius_ussuriensis'
                else ''
            ),
            interval_status=status
        )

        manifest_rows.append(row)

        plan_rows.append(dict(
            species=species,
            family=family,
            chromosome=chrom,
            start=start,
            end=end,
            fasta=str(dest),
            expected_bp=end-start+1
        ))


# verify Gloydius ToxCodAn locus containment
region_lookup = {
    (r['species'], r['family']): r
    for r in manifest_rows
}

checks = []

for (family, gene, chrom), (start, end) in sorted(
    gtf_loci.items()
):

    region = region_lookup[
        ('Gloydius_ussuriensis', family)
    ]

    if gene in excluded_loci:
        status = 'excluded_from_core_toxin_set'

    elif chrom != region['chromosome']:
        status = 'other_chromosome_not_this_target_cluster'

    elif (
        region['start'] <= start
        and end <= region['end']
    ):
        status = 'retained_GTF_locus_contained_pending_curation'

    else:
        raise RuntimeError(
            f'Retained {gene} lies outside {family} interval: '
            f'{chrom}:{start}-{end}'
        )

    checks.append(dict(
        family=family,
        gene_id=gene,
        chromosome=chrom,
        start=start,
        end=end,
        region_chr=region['chromosome'],
        region_start=region['start'],
        region_end=region['end'],
        status=status
    ))


# export interval manifests and evidence tables
write_tsv(
    anchor_file,
    anchors,
    [
        'species', 'family', 'anchor_symbol',
        'chromosome', 'start', 'end', 'strand',
        'gene_ids', 'protein_ids', 'evidence_status'
    ]
)

write_tsv(
    toxin_file,
    checks,
    [
        'family', 'gene_id', 'chromosome',
        'start', 'end', 'region_chr',
        'region_start', 'region_end', 'status'
    ]
)

write_tsv(
    manifest,
    manifest_rows,
    [
        'species', 'family', 'chromosome',
        'start', 'end', 'length_bp', 'fasta',
        'anchor_count', 'anchor_symbols',
        'anchor_coordinate_order',
        'Gloydius_retained_toxin_loci',
        'interval_status'
    ]
)

write_tsv(
    plan,
    plan_rows,
    [
        'species', 'family', 'chromosome',
        'start', 'end', 'fasta', 'expected_bp'
    ]
)

for r in manifest_rows:
    print(
        f"{r['species']} {r['family']}: "
        f"{r['chromosome']}:{r['start']}-{r['end']} "
        f"({r['length_bp']:,} bp)"
    )

print('Anchor evidence:', anchor_file)
print('ToxCodAn locus check:', toxin_file)

PY

while IFS=$'\t' read -r SPECIES FAMILY CHR START END OUTPUT EXPECTED; do

    [[ "$SPECIES" == 'species' ]] && continue

    case "$SPECIES" in
        Gloydius_ussuriensis) GENOME="$GUSS" ;;
        Cerastes_gasperettii) GENOME="$CGAS" ;;
        Crotalus_adamanteus) GENOME="$CADAM" ;;
        Naja_naja) GENOME="$NAJA" ;;
        *) echo "ERROR: unknown species: $SPECIES" >&2; exit 1 ;;
    esac

    echo "Extracting $SPECIES $FAMILY: $CHR:$START-$END"

    TMP="${OUTPUT}.tmp.$$"

    samtools faidx "$GENOME" "$CHR:$START-$END" > "$TMP"

    LENGTH=$(awk '!/^>/ {n += length($0)} END {print n+0}' "$TMP")

    if (( LENGTH != EXPECTED )); then

        rm -f -- "$TMP"

        echo "ERROR: length mismatch: $SPECIES $FAMILY ($LENGTH != $EXPECTED)" >&2

        exit 1

    fi

    mv -f -- "$TMP" "$OUTPUT"

    rm -f -- "${OUTPUT}.fai" "${OUTPUT}.gzi"

done < "$PLAN"

# 4. summarize extracted intervals
echo
echo '===== GENOMIC INTERVAL SUMMARY ====='

if command -v column >/dev/null; then
    column -t -s $'\t' "$MANIFEST"
else
    cat "$MANIFEST"
fi

echo
echo '===== OUTPUTS ====='

printf '%s\n' "$MANIFEST" "$ANCHORS" "$TOXIN_CHECK"

echo '===== COMPLETE ====='
```

Inspect the results after running the script:
```sh
# ============================================================
# inspect genomic interval extraction results
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
INVENTORY="$ROOT/synteny/toxin_microsynteny/01_reference_inventory"

# 2. inspect interval manifest
echo
echo "===== ANCHOR-BOUNDED GENOMIC INTERVALS ====="

column -t -s $'\t' "$INVENTORY/anchor_bounded_regions.tsv"

# 3. inspect selected anchors
echo
echo "===== SELECTED ANCHOR COORDINATES ====="

column -t -s $'\t' "$INVENTORY/selected_anchor_coordinates.tsv"

# 4. inspect Gloydius toxin containment
echo
echo "===== GLOYDIUS TOXCODAN LOCUS CONTAINMENT ====="

column -t -s $'\t' "$INVENTORY/Guss_ToxCodAn_region_containment.tsv"

# 5. verify FASTAs and record counts
echo
echo "===== STEP 5 VALIDATION ====="

python - "$INVENTORY" <<'PY'
import csv
import sys
from pathlib import Path
from collections import Counter

root = Path(sys.argv[1])

def load(filename):
    with (root / filename).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

regions = load("anchor_bounded_regions.tsv")
anchors = load("selected_anchor_coordinates.tsv")
toxins = load("Guss_ToxCodAn_region_containment.tsv")

print(f"Extracted regions: {len(regions)} / 8")
print(f"Selected anchors: {len(anchors)} / 32")

print("\n--- FASTA verification ---")

for r in regions:
    path = Path(r["fasta"])
    status = "OK" if path.is_file() and path.stat().st_size > 0 else "MISSING"

    print(
        f"{r['species']:<23} "
        f"{r['family']:<5} "
        f"{int(r['length_bp']):>10,} bp "
        f"{status}"
    )

print("\n--- ToxCodAn containment status ---")

for status, count in Counter(r["status"] for r in toxins).items():
    print(f"{status}: {count}")

print("\n===== COMPLETE =====")
PY
```

#### Step 6: Chromosome-coordinate inventory of the existing SVMP and PLA2 annotations in *Cerastes*, *Crotalus*, and *Naja* 
Now, let's establish an evidence-based inventory of the genes already annotated within our eight SVMP and PLA2 genomic intervals, before searching for additional unannotated toxin loci.

The script below will do the following:
  1) Inventories annotated genes and transcripts from all four species, allowing examination of gene organization without confusing transcript isoforms with independent genes.
  2) Cross-references the conserved flanking anchors from Step 4 and Step 5, including the provisional ADAM28-like locus and the candidate fragmented VWA5B1 annotation.
  3) Identifies preliminary toxin candidates from GFF3 annotations and BLAST evidence, while distinguishing these from nonvenom homologs and conserved anchors.
  4) Inventories Gloydius ToxCodAn loci, preserving their original IDs and explicitly identifying SVMP-6 and PLA2-1 as exclusions from the relevant toxin counts. Other loci remain provisional until reconciled against the complete manually curated set.
  5) Flags unresolved annotations, including FUN_018074–075, the VWA5B1 models, and the Cerastes venom PLA2 copy-number discrepancy.

All outputs are written to:
```sh
$ROOT/synteny/toxin_microsynteny/05_annotations/
```
The output files will be: 
  - `comparative_GFF_models_in_regions.tsv`: Transcript-level inventory for the three comparative species
  - `comparative_toxin_keyword_candidates.tsv`: Preliminary keyword-matched models
  - `Guss_ToxCodAn_loci_in_regions.tsv`: Gloydius ToxCodAn locus inventory
  - `comparative_GFF_genes_in_regions.tsv`: Distinct GFF gene records in the comparative species
  - `Guss_Funannotate_genes_in_regions.tsv`: Gloydius Funannotate gene records
  - `comparative_gene_evidence_inventory.tsv`: Integrated evidence from gene models, selected anchors, and ToxCodAn
  - `comparative_toxin_candidate_summary.tsv`: Per-species, per-family summary
  - `Step6_annotation_review_flags.tsv`: Outstanding manual validation tasks

```sh
# ============================================================
# inventory existing SVMP and PLA2 gene annotations
# ============================================================

# 0. activate conda env
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

# 1. directories
ROOT="${ROOT:-$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo}"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
ASSEMBLIES="$ROOT/synteny/assemblies_synteny"
OUTDIR="$SYNTENY/05_annotations"
MANIFEST="$SYNTENY/01_reference_inventory/anchor_bounded_regions.tsv"
ANCHORS="$SYNTENY/01_reference_inventory/selected_anchor_coordinates.tsv"
EVIDENCE="$SYNTENY/02_blast/flanking_homology_evidence"
TOXIN="$ROOT/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/toxin_annotation.gtf"
GUSS_GFF="$ROOT/annotation/funannotate/update_results/Gloydius_ussuriensis_AMNH_21010.gff3"

mkdir -p "$OUTDIR"
command -v python >/dev/null || { echo 'ERROR: python is required' >&2; exit 1; }

# 2. comparative annotation files
CGAS="$ASSEMBLIES/Cerastes_gasperettii/Cerastes_gasperettii.annotation.gff3"
CADAM="$ASSEMBLIES/Crotalus_adamanteus/Crotalus_adamanteus.annotation.gff3"
NAJA="$ASSEMBLIES/Naja_naja/Naja_naja.annotation.gff3"

# 3. verify input files
echo '===== VERIFY INPUT FILES ====='
for FILE in "$MANIFEST" "$ANCHORS" "$TOXIN" "$GUSS_GFF" "$CGAS" "$CADAM" "$NAJA" \
    "$EVIDENCE/anchor_reference_blastp_candidates.tsv" \
    "$EVIDENCE/flanking_blastp_reciprocal_gff_evidence.tsv"; do
    [[ -s "$FILE" ]] || { echo "ERROR: missing or empty input: $FILE" >&2; exit 1; }
    echo "OK: $FILE"
done

# 4. extract existing annotations
python - "$MANIFEST" "$ANCHORS" "$OUTDIR" "$TOXIN" "$EVIDENCE" \
    "$GUSS_GFF" "$CGAS" "$CADAM" "$NAJA" <<'PY'
import csv
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from urllib.parse import unquote

manifest_path, anchors_path, outdir, toxin_path, evidence = map(Path, sys.argv[1:6])
gff_paths = dict(zip(
    ('Gloydius_ussuriensis', 'Cerastes_gasperettii',
     'Crotalus_adamanteus', 'Naja_naja'), map(Path, sys.argv[6:10])
))
comp_species = tuple(s for s in gff_paths if s != 'Gloydius_ussuriensis')
families = ('SVMP', 'PLA2')

# ------------------------------------------------------------
# 1. read genomic interval manifest and selected anchors
# ------------------------------------------------------------
def read_tsv(path, required=()):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        return list(reader)


def write_tsv(name, rows, columns):
    dest = outdir / name
    temp = outdir / f'.{name}.{os.getpid()}.tmp'
    with temp.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t',
                                lineterminator='\n', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temp, dest)
    print(f'OUTPUT: {dest} ({len(rows)} rows)')


regions = {}
for r in read_tsv(manifest_path, ('species', 'family', 'chromosome', 'start', 'end', 'fasta')):
    key = (r['species'], r['family'])
    if key in regions or r['species'] not in gff_paths or r['family'] not in families:
        raise ValueError(f'Duplicate/unexpected Step 5 region: {key}')
    r['start'], r['end'] = int(r['start']), int(r['end'])
    if not (1 <= r['start'] < r['end']):
        raise ValueError(f'Invalid Step 5 interval: {r}')
    if not Path(r['fasta']).is_file() or Path(r['fasta']).stat().st_size == 0:
        raise FileNotFoundError(f'Missing Step 5 FASTA: {r["fasta"]}')
    regions[key] = r
if len(regions) != 8 or any((s, f) not in regions for s in gff_paths for f in families):
    raise ValueError('Expected exactly eight Step 5 genomic intervals')

selected_anchors = defaultdict(list)
anchor_by_gene = defaultdict(list)
for a in read_tsv(anchors_path, ('species', 'family', 'anchor_symbol', 'chromosome',
                                 'start', 'end', 'gene_ids', 'evidence_status')):
    key = (a['species'], a['family'])
    if key not in regions:
        raise ValueError(f'Selected anchor has no region: {key}')
    region = regions[key]
    a['start'], a['end'] = int(a['start']), int(a['end'])
    if (a['chromosome'] != region['chromosome'] or
            not (region['start'] <= a['start'] <= a['end'] <= region['end'])):
        raise ValueError(f'Anchor outside Step 5 interval: {key} {a["anchor_symbol"]}')
    selected_anchors[key].append(a)
    for gid in filter(None, a['gene_ids'].split(';')):
        if not gid.startswith('unannotated_'):
            anchor_by_gene[(a['species'], a['family'], gid)].append(a)
if sum(map(len, selected_anchors.values())) != 32:
    raise ValueError('Expected 32 Step 5 anchor records')

# ------------------------------------------------------------
# 2. annotation parsing and conservative evidence labels
# ------------------------------------------------------------
def parse_gff(text):
    result = {}
    for part in text.split(';'):
        if '=' in part:
            key, value = part.split('=', 1)
            result[key.strip()] = unquote(value.strip())
    return result


def parse_gtf(text):
    return dict(re.findall(r'(\S+)\s+"([^"]*)"', text))


def overlap(a, b, c, d):
    return a <= d and b >= c


def split_ids(value):
    return [v for v in value.split(',') if v]


patterns = {
    'SVMP': re.compile(r'\bSVMP\b|venom.{0,60}metalloprote|metalloprote.{0,60}venom|'
                       r'snake venom metalloproteinase', re.I),
    'PLA2': re.compile(r'\bPLA\s*[-_]?\s*2\b|\bPLA2[A-Z0-9_-]*\b|'
                       r'phospholipase\s*A\s*2|phospholipase\s*A2', re.I)
}
venom_marker = re.compile(r'\bvenom\b|\btoxin\b|\bSVMP\b', re.I)


def classify(text, family, anchor_symbols, gene_id):
    keyword = bool(patterns[family].search(text))
    if anchor_symbols:
        return ('family_keyword' if keyword else 'no_keyword',
                'selected_flanking_anchor_candidate')
    if gene_id in ('FUN_018074', 'FUN_018075') and family == 'PLA2':
        return ('family_keyword' if keyword else 'no_keyword',
                'unresolved_neighboring_gene_models')
    if keyword and venom_marker.search(text):
        return 'family_keyword', 'venom_toxin_candidate_from_annotation'
    if keyword:
        return 'family_keyword', 'family_homolog_candidate_function_unresolved'
    return 'no_keyword', 'other_annotated_gene'


def functional_hint(attributes):
    for key in ('product', 'gene', 'gene_name', 'description', 'Note', 'Name'):
        for attrs in attributes:
            if attrs.get(key):
                return attrs[key][:250]
    return 'uncharacterized'

# ------------------------------------------------------------
# 3. read comparative and Gloydius GFF3 annotations
# ------------------------------------------------------------
def build_models(species, family, region, path):
    chromosome, left, right = region['chromosome'], region['start'], region['end']
    genes, transcripts, cds_features = {}, {}, []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chromosome or f[2] not in ('gene', 'pseudogene', 'mRNA', 'transcript', 'CDS'):
                continue
            try:
                start, end = int(f[3]), int(f[4])
            except ValueError:
                continue
            if not overlap(start, end, left, right):
                continue
            record = dict(start=start, end=end, strand=f[6], attrs=parse_gff(f[8]))
            if f[2] in ('gene', 'pseudogene'):
                gid = record['attrs'].get('ID') or f'unnamed_gene_{chromosome}_{start}_{end}'
                genes[gid] = record
            elif f[2] in ('mRNA', 'transcript'):
                tid = record['attrs'].get('ID') or f'unnamed_transcript_{chromosome}_{start}_{end}'
                transcripts[tid] = record
            else:
                cds_features.append(record)

    grouped = defaultdict(lambda: {'gene': None, 'transcripts': {}, 'cds': []})
    for gid, gene in genes.items():
        grouped[gid]['gene'] = gene
    tx_gene = {}
    for tid, tx in transcripts.items():
        parents = split_ids(tx['attrs'].get('Parent', '') or tx['attrs'].get('gene_id', ''))
        gid = next((p for p in parents if p in genes), parents[0] if parents else f'no_gene_parent_{tid}')
        tx_gene[tid] = gid
        grouped[gid]['transcripts'][tid] = tx
    cds_by_tx = defaultdict(list)
    for feature in cds_features:
        parents = split_ids(feature['attrs'].get('Parent', ''))
        for parent in parents or [f'no_CDS_parent_{chromosome}_{feature["start"]}_{feature["end"]}']:
            if parent in tx_gene:
                gid = tx_gene[parent]
                cds_by_tx[parent].append(feature)
            elif parent in genes:
                gid = parent
            else:
                gid = f'no_gene_for_CDS_parent_{parent}'
            grouped[gid]['cds'].append(feature)

    gene_rows, model_rows = [], []
    for gid, bundle in grouped.items():
        gene = bundle['gene']
        txs = bundle['transcripts']
        cdss = bundle['cds']
        records = ([gene] if gene else []) + list(txs.values()) + cdss
        if not records:
            continue
        start = gene['start'] if gene else min(r['start'] for r in records)
        end = gene['end'] if gene else max(r['end'] for r in records)
        strand = gene['strand'] if gene else records[0]['strand']
        attrs = [r['attrs'] for r in records]
        combined = ' '.join(' '.join(map(str, d.values())) for d in attrs)
        linked_anchors = anchor_by_gene.get((species, family, gid), [])
        symbols = ';'.join(sorted(set(a['anchor_symbol'] for a in linked_anchors)))
        keyword, category = classify(combined, family, symbols, gid)
        gene_status = 'gene_model' if gene else 'gene_feature_missing_in_GFF_interval'
        boundary = 'partial_interval_overlap' if start < left or end > right else 'entire_gene_span_inside_interval'
        distinct_cds = {(c['start'], c['end'], c['strand']) for c in cdss}
        review_flag = ''
        if species == 'Gloydius_ussuriensis' and family == 'PLA2':
            if gid in {'FUN_018053', 'FUN_018054', 'FUN_018055'}:
                review_flag = 'VWA5B1_three_models_one_provisional_locus_group'
            elif gid in {'FUN_018074', 'FUN_018075'}:
                review_flag = 'neighboring_pair_distinct_models_manual_review'
        row = dict(species=species, family=family, chr=chromosome, start=start, end=end,
                   strand=strand, gene_id=gid, transcript_count=len(txs),
                   transcript_ids=';'.join(sorted(txs)), cds_segments=len(distinct_cds),
                   functional_hint=functional_hint(attrs), classification=keyword,
                   evidence_class=category, selected_anchor_symbols=symbols,
                   anchor_assignment='selected_step5_gene_id' if symbols else '',
                   gene_model_status=gene_status, interval_boundary_status=boundary,
                   review_flag=review_flag, model_source='GFF3')
        gene_rows.append(row)
        if txs:
            for tid, tx in sorted(txs.items()):
                subset = cds_by_tx.get(tid, [])
                model_rows.append(dict(**{
                    'species': species, 'family': family, 'chr': chromosome,
                    'start': tx['start'], 'end': tx['end'], 'strand': tx['strand'],
                    'gene_id': gid, 'model_id': tid, 'model_source': 'transcript',
                    'cds_segments': len({(c['start'], c['end'], c['strand']) for c in subset}),
                    'functional_hint': row['functional_hint'], 'classification': keyword,
                    'evidence_class': category, 'selected_anchor_symbols': symbols,
                    'interval_boundary_status': boundary
                }))
        else:
            model_rows.append(dict(
                species=species, family=family, chr=chromosome, start=start, end=end,
                strand=strand, gene_id=gid, model_id=gid,
                model_source='gene_only' if gene and not cdss else 'CDS_parent_only',
                cds_segments=len(distinct_cds), functional_hint=row['functional_hint'],
                classification=keyword, evidence_class=category,
                selected_anchor_symbols=symbols, interval_boundary_status=boundary
            ))
    return gene_rows, model_rows


genes_all, models_all = [], []
for species, gff in gff_paths.items():
    print(f'\n===== {species} =====')
    for family in families:
        gene_rows, model_rows = build_models(species, family, regions[(species, family)], gff)
        genes_all.extend(gene_rows)
        models_all.extend(model_rows)
        print(f'{family}: {len(gene_rows)} gene records, {len(model_rows)} transcript/model records')

# step 4 reference-anchored BLASTp candidates and reciprocal hits are supporting evidence; neither by itself establishes orthology or toxin function
anchor_hits = defaultdict(list)
for r in read_tsv(evidence / 'anchor_reference_blastp_candidates.tsv',
                  ('family', 'species', 'gene_id', 'subject_protein_id', 'bitscore')):
    if r['gene_id']:
        anchor_hits[(r['species'], r['family'], r['gene_id'])].append(r)

reciprocal_hits = defaultdict(list)
for r in read_tsv(evidence / 'flanking_blastp_reciprocal_gff_evidence.tsv',
                  ('family', 'reference_species', 'gff_gene_id', 'subject_protein_id', 'reciprocal_rank')):
    if r['gff_gene_id']:
        reciprocal_hits[(r['reference_species'], r['family'], r['gff_gene_id'])].append(r)

for row in genes_all:
    key = (row['species'], row['family'], row['gene_id'])
    ah, rh = anchor_hits[key], reciprocal_hits[key]
    row['Step4_anchor_hit_count'] = len(ah)
    row['Step4_best_anchor_bitscore'] = max((float(r['bitscore']) for r in ah if r['bitscore']), default='')
    row['Step4_anchor_query_ids'] = ';'.join(sorted(set(r.get('anchor_query', '') for r in ah if r.get('anchor_query'))))
    row['Step4_reciprocal_hit_count'] = len(rh)
    ranks = [int(r['reciprocal_rank']) for r in rh if r['reciprocal_rank'].isdigit()]
    row['Step4_best_reciprocal_rank'] = min(ranks) if ranks else ''
    row['Step4_reciprocal_query_gene_ids'] = ';'.join(sorted(set(r.get('query_gene_id', '') for r in rh if r.get('query_gene_id'))))

# ------------------------------------------------------------
# 4. recover Gloydius ToxCodAn loci and record curation status
# ------------------------------------------------------------
raw_loci = {}
with toxin_path.open() as handle:
    for line in handle:
        if not line.strip() or line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        gid = parse_gtf(f[8]).get('gene_id', '')
        family = next((x for x in families if gid.startswith(x + '-')), None)
        if family is None:
            continue
        key = (family, gid, f[0])
        item = raw_loci.setdefault(key, {'segments': set(), 'strands': set()})
        item['segments'].add((int(f[3]), int(f[4])))
        item['strands'].add(f[6])

excluded = {'SVMP-6': 'redundant_overlapping_model_excluded',
            'PLA2-1': 'PLA2I_like_excluded_from_enzymatic_PLA2_count'}
toxin_loci = []
for (family, gid, chrom), item in sorted(raw_loci.items()):
    region = regions[('Gloydius_ussuriensis', family)]
    start = min(a for a, b in item['segments'])
    end = max(b for a, b in item['segments'])
    if chrom != region['chromosome'] or not overlap(start, end, region['start'], region['end']):
        continue
    fully_contained = region['start'] <= start and end <= region['end']
    if not fully_contained:
        raise ValueError(f'ToxCodAn locus not fully included in Step 5 region: {gid}')
    status = excluded.get(gid, 'core_status_not_established_from_GTF')
    toxin_loci.append(dict(
        species='Gloydius_ussuriensis', family=family, chr=chrom,
        start=start, end=end, strand=';'.join(sorted(item['strands'])),
        gene_id=gid, cds_segments=len(item['segments']),
        curation_status=status, included_in_provisional_count=(gid not in excluded),
        full_locus_in_interval='yes', source_annotation='ToxCodAn_GTF'
    ))

for row in genes_all:
    if row['species'] != 'Gloydius_ussuriensis':
        row['overlapping_ToxCodAn_IDs'] = ''
        continue
    row['overlapping_ToxCodAn_IDs'] = ';'.join(
        sorted(t['gene_id'] for t in toxin_loci
               if t['family'] == row['family'] and t['chr'] == row['chr']
               and overlap(row['start'], row['end'], t['start'], t['end']))
    )

# ------------------------------------------------------------
# 5. export annotation inventories and evidence tables
# ------------------------------------------------------------
genes_all.sort(key=lambda r: (r['species'], r['family'], r['chr'], r['start'], r['gene_id']))
models_all.sort(key=lambda r: (r['species'], r['family'], r['chr'], r['start'], r['model_id']))
toxin_loci.sort(key=lambda r: (r['family'], r['chr'], r['start'], r['gene_id']))

model_cols = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id',
              'model_id', 'model_source', 'cds_segments', 'functional_hint',
              'classification', 'evidence_class', 'selected_anchor_symbols',
              'interval_boundary_status']
gene_cols = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id',
             'transcript_count', 'transcript_ids', 'cds_segments', 'functional_hint',
             'classification', 'evidence_class', 'selected_anchor_symbols',
             'anchor_assignment', 'gene_model_status', 'interval_boundary_status',
             'review_flag', 'model_source', 'Step4_anchor_hit_count',
             'Step4_best_anchor_bitscore', 'Step4_anchor_query_ids',
             'Step4_reciprocal_hit_count', 'Step4_best_reciprocal_rank',
             'Step4_reciprocal_query_gene_ids', 'overlapping_ToxCodAn_IDs']
toxin_cols = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id',
              'cds_segments', 'curation_status', 'included_in_provisional_count',
              'full_locus_in_interval', 'source_annotation']

comparative_models = [r for r in models_all if r['species'] in comp_species]
comparative_genes = [r for r in genes_all if r['species'] in comp_species]

def get_gene_col(r, col):
    return r.get(col, '')

# populate the original model-level output with original column names preserved
write_tsv('comparative_GFF_models_in_regions.tsv', comparative_models, model_cols)
write_tsv('comparative_toxin_keyword_candidates.tsv',
          [r for r in comparative_models if r['classification'] == 'family_keyword'], model_cols)
write_tsv('Guss_ToxCodAn_loci_in_regions.tsv', toxin_loci, toxin_cols)
write_tsv('comparative_GFF_genes_in_regions.tsv', comparative_genes, gene_cols)
write_tsv('Guss_Funannotate_genes_in_regions.tsv',
          [r for r in genes_all if r['species'] == 'Gloydius_ussuriensis'], gene_cols)

# integrated table contains GFF gene records, Gloydius ToxCodAn loci, and selected genomic anchors without a corresponding annotated GFF gene
integrated = []
for g in genes_all:
    integrated.append(dict(g, record_type='GFF_gene_model',
                           source_record_id=g['gene_id'],
                           source_curation_status='',
                           provisional_locus_group=(
                               'Guss_VWA5B1_candidate' if g['review_flag'].startswith('VWA5B1_') else ''
                           )))
for t in toxin_loci:
    integrated.append(dict(t, record_type='ToxCodAn_toxin_locus',
                           source_record_id=t['gene_id'],
                           evidence_class='venom_toxin_model_provisional_curation',
                           selected_anchor_symbols='', anchor_assignment='',
                           overlapping_ToxCodAn_IDs=t['gene_id'],
                           provisional_locus_group='',
                           source_curation_status=t['curation_status']))
for key, anchors in selected_anchors.items():
    for a in anchors:
        matches = [g for g in genes_all if g['species'] == key[0]
                   and g['family'] == key[1]
                   and a['anchor_symbol'] in g['selected_anchor_symbols'].split(';')]
        if matches:
            continue
        integrated.append(dict(
            species=a['species'], family=a['family'], chr=a['chromosome'],
            start=a['start'], end=a['end'], strand=a.get('strand', ''),
            gene_id='', record_type='unannotated_selected_anchor_locus',
            source_record_id=f'{a["species"]}_{a["family"]}_{a["anchor_symbol"]}',
            evidence_class='provisional_genomic_anchor_candidate',
            selected_anchor_symbols=a['anchor_symbol'],
            anchor_assignment='genomic_coordinate_only',
            source_curation_status='not_confirmed_ortholog',
            provisional_locus_group='', review_flag='no_selected_GFF_gene_model'
        ))

integrated.sort(key=lambda r: (r['species'], r['family'], r['chr'], int(r['start']),
                               r['record_type'], r['source_record_id']))
integrated_cols = ['species', 'family', 'chr', 'start', 'end', 'strand',
                   'record_type', 'source_record_id', 'gene_id', 'transcript_count',
                   'transcript_ids', 'cds_segments', 'functional_hint',
                   'classification', 'evidence_class', 'selected_anchor_symbols',
                   'anchor_assignment', 'gene_model_status',
                   'interval_boundary_status', 'provisional_locus_group',
                   'review_flag', 'Step4_anchor_hit_count', 'Step4_best_anchor_bitscore',
                   'Step4_anchor_query_ids', 'Step4_reciprocal_hit_count',
                   'Step4_best_reciprocal_rank', 'Step4_reciprocal_query_gene_ids',
                   'overlapping_ToxCodAn_IDs', 'source_curation_status']
write_tsv('comparative_gene_evidence_inventory.tsv', integrated, integrated_cols)

summary = []
for species in gff_paths:
    for family in families:
        gs = [g for g in genes_all if g['species'] == species and g['family'] == family]
        ts = [t for t in toxin_loci if species == 'Gloydius_ussuriensis' and t['family'] == family]
        summary.append(dict(
            species=species, family=family,
            total_GFF_gene_records=len(gs),
            annotation_keyword_gene_candidates=sum(g['classification'] == 'family_keyword' for g in gs),
            annotation_venom_candidate_genes=sum(g['evidence_class'] == 'venom_toxin_candidate_from_annotation' for g in gs),
            selected_anchor_gene_models=sum(bool(g['selected_anchor_symbols']) for g in gs),
            raw_ToxCodAn_locus_IDs=len(ts) if species == 'Gloydius_ussuriensis' else '',
            provisionally_included_ToxCodAn_IDs=sum(t['included_in_provisional_count'] for t in ts) if species == 'Gloydius_ussuriensis' else '',
            interpretation='requires_manual_toxin_family_confirmation_and_copy_number_review'
        ))
write_tsv('comparative_toxin_candidate_summary.tsv', summary,
          ['species', 'family', 'total_GFF_gene_records',
           'annotation_keyword_gene_candidates', 'annotation_venom_candidate_genes',
           'selected_anchor_gene_models', 'raw_ToxCodAn_locus_IDs',
           'provisionally_included_ToxCodAn_IDs', 'interpretation'])

flags = [
    dict(species='Gloydius_ussuriensis', family='PLA2',
         locus='FUN_018053;FUN_018054;FUN_018055',
         issue='three_separate_GFF_models_with_complementary_VWA5B1_protein_alignments',
         action='retain_original_IDs;evaluate_as_one_provisional_locus;review_splice_support'),
    dict(species='Gloydius_ussuriensis', family='PLA2',
         locus='FUN_018074;FUN_018075',
         issue='neighboring_gene_models_with_different_cross_species_best_matches',
         action='retain_as_distinct_provisional_models;review_exons_and_reference_annotations'),
    dict(species='Gloydius_ussuriensis', family='SVMP',
         locus='G_ussuri_chr9:5383702-5413622',
         issue='ADAM28_like_TBLASTN_candidate_not_validated_ortholog',
         action='review_domains_exon_structure_and_local_synteny'),
    dict(species='Cerastes_gasperettii', family='PLA2',
         locus='venom_PLA2_cluster',
         issue='figure_three_candidates_versus_two_published_venom_PLA2_genes',
         action='compare_gene_IDs_coordinates_CDS_models_and_published_toxin_catalog'),
    dict(species='Gloydius_ussuriensis', family='SVMP;PLA2',
         locus='ToxCodAn_core_set',
         issue='GTF_only_cannot_establish_full_manually_curated_core_set',
         action='reconcile_provisional_IDs_with_curated_locus_list_before_reporting_copies')
]
write_tsv('Step6_annotation_review_flags.tsv', flags,
          ['species', 'family', 'locus', 'issue', 'action'])

# ------------------------------------------------------------
# 6. print summary
# ------------------------------------------------------------
print('\n===== COMPARATIVE TOXIN ANNOTATION SUMMARY =====')
for r in summary:
    print(f"{r['species']:<23} {r['family']:<5} "
          f"genes={r['total_GFF_gene_records']:<4} "
          f"keyword_genes={r['annotation_keyword_gene_candidates']:<3} "
          f"venom_annotation_candidates={r['annotation_venom_candidate_genes']}")
print('\n===== GLOYDIUS TOXCODAN LOCUS INVENTORY =====')
for family in families:
    ts = [t for t in toxin_loci if t['family'] == family]
    n = sum(t['included_in_provisional_count'] for t in ts)
    print(f'{family}: {len(ts)} raw in-region IDs; {n} provisionally included; '
          f'{len(ts)-n} explicitly excluded')
print('\nNOTE: keyword matches, protein homology and provisional loci are NOT validated toxin copies or ortholog calls.')
print('NOTE: comparative toxin-copy completeness and the Cerastes PLA2 discrepancy require manual review.')
print('===== COMPLETE =====')
```

Run this as a shell script:
```sh
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Shell/"

chmod +x "${ROOT}/toxin_microsynteny_annotation_inventory.sh"
bash "${ROOT}/toxin_microsynteny_annotation_inventory.sh"
```
The output 

After the script finishes running, run the diagnostic script below to examine individual gene IDs before proceeding to Step 7. The primary goal here is to inspect the Step 6 output files and identify potential annotation inconsistencies before we use those results for comparative microsynteny.
```sh
# ============================================================
# review toxin candidates and annotation discrepancies
# ============================================================

# 1. directories
ROOT="$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
OUTDIR="$SYNTENY/05_annotations"
INVENTORY="$SYNTENY/01_reference_inventory"

# 2. inspect Gloydius PLA2 ToxCodAn loci
echo
echo "===== GLOYDIUS PLA2 TOXCODAN LOCI ====="

awk -F '\t' 'NR == 1 || $2 == "PLA2" {print}' \
    "$OUTDIR/Guss_ToxCodAn_loci_in_regions.tsv"

# 3. locate PLA2-1 in the Step 5 containment report
echo
echo "===== PLA2-1 LOCATION AND STATUS ====="

awk -F '\t' 'NR == 1 || $2 == "PLA2-1" {print}' \
    "$INVENTORY/Guss_ToxCodAn_region_containment.tsv"

# 4. inspect Cerastes PLA2 keyword candidates
echo
echo "===== CERASTES PLA2 ANNOTATION CANDIDATES ====="

awk -F '\t' 'NR == 1 || ($1 == "Cerastes_gasperettii" && $2 == "PLA2" && $12 == "family_keyword") {print}' \
    "$OUTDIR/comparative_GFF_genes_in_regions.tsv"

# 5. inspect Crotalus PLA2 keyword candidates
echo
echo "===== CROTALUS PLA2 ANNOTATION CANDIDATES ====="

awk -F '\t' 'NR == 1 || ($1 == "Crotalus_adamanteus" && $2 == "PLA2" && $12 == "family_keyword") {print}' \
    "$OUTDIR/comparative_GFF_genes_in_regions.tsv"

# 6. inspect unresolved annotation flags
echo
echo "===== ANNOTATION REVIEW FLAGS ====="

column -t -s $'\t' "$OUTDIR/Step6_annotation_review_flags.tsv"

echo
echo "===== COMPLETE ====="
```


#### Step 7: Generate one consistent gene-coordinate catalog for microsynteny plotting
The script below will perform data-consolidation and preliminary gene-classification based on outputs generated in previous steps and create one output: `comparative_all_gene_models_in_regions.tsv`
```sh
# ============================================================
# inventory all annotated genes in homologous regions
# ============================================================

# activate conda env to use python
conda activate funannotate

python - <<'PY'
# 1. directories
import csv
import os
from collections import Counter
from pathlib import Path

ROOT = Path(os.environ.get('ROOT', str(Path.home() / 'mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo')))
OUTDIR = ROOT / 'synteny/toxin_microsynteny/05_annotations'
INPUT = OUTDIR / 'comparative_gene_evidence_inventory.tsv'
OUTPUT = OUTDIR / 'comparative_all_gene_models_in_regions.tsv'

# 2. verify input files
if not INPUT.is_file() or INPUT.stat().st_size == 0:
    raise SystemExit(f'ERROR: Missing Step 6 inventory: {INPUT}')

with INPUT.open(newline='') as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    required = {'species', 'family', 'chr', 'start', 'end', 'strand', 'record_type',
                'source_record_id', 'gene_id', 'functional_hint', 'evidence_class',
                'selected_anchor_symbols', 'overlapping_ToxCodAn_IDs',
                'source_curation_status', 'review_flag', 'provisional_locus_group',
                'classification'}
    if required - set(reader.fieldnames or []):
        raise SystemExit(f'ERROR: Missing columns: {sorted(required - set(reader.fieldnames or []))}')
    records = list(reader)
if not records:
    raise SystemExit('ERROR: Step 6 inventory is empty')

# 3. build minimal comparative gene catalog
# Prior Naja homology labels are provisional; no new BLAST is performed.
naja_svmp = {'006949', '006952', '006954', '006955', '006956', '006957'}
naja_short_svmp = {'006951', '006953'}
naja_pla2 = {'006365'}
rows, vwa = [], []

for r in records:
    species, family, gene = r['species'], r['family'], r['gene_id']
    if species not in {'Gloydius_ussuriensis', 'Cerastes_gasperettii',
                       'Crotalus_adamanteus', 'Naja_naja'} or family not in {'SVMP', 'PLA2'}:
        raise SystemExit(f'ERROR: Unexpected species/family: {species} {family}')
    start, end = int(r['start']), int(r['end'])
    if start < 1 or end < start:
        raise SystemExit(f'ERROR: Invalid gene interval: {species} {gene} {start}-{end}')
    rtype = r['record_type']
    anchor = r['selected_anchor_symbols']
    hint = r['functional_hint'].lower()
    status = r['source_curation_status']
    note = r['review_flag'] or ''
    role, include = 'other_annotated_gene', 1

    if rtype == 'ToxCodAn_toxin_locus':
        role = 'putative_venom_toxin'
        if 'excluded' in status or 'redundant' in status or 'pla2i_like' in status.lower():
            role, include = 'excluded_toxin_model', 0
        else:
            note = 'Provisional ToxCodAn locus; core membership not verified here'
    elif rtype == 'unannotated_selected_anchor_locus':
        role = 'unannotated_anchor_candidate'
    elif anchor:
        role = 'conserved_anchor_candidate'
        if species == 'Crotalus_adamanteus' and family == 'PLA2' and 'PLA2G2F' in anchor and 'venom' in hint:
            role = 'ambiguous_anchor_or_venom_PLA2'
            note = 'Conflicting positional PLA2G2F assignment and VENOM-PLA2-2 description'
    elif species == 'Cerastes_gasperettii' and family == 'PLA2' and gene in {'CG-PLA2-gC', 'CG-PLA2-gD'}:
        role = 'putative_venom_toxin'
        note = 'Two explicitly named Cerastes venom PLA2 candidates'
    elif species == 'Naja_naja':
        suffix = gene.rsplit('_', 1)[-1]
        if family == 'SVMP' and suffix in naja_svmp:
            role, note = 'SVMP_like_candidate', 'Provisional from earlier Naja BLASTp analysis'
        elif family == 'SVMP' and suffix in naja_short_svmp:
            role, note = 'short_SVMP_like_model', 'Short prior BLASTp alignments; not a confirmed copy'
        elif family == 'PLA2' and suffix in naja_pla2:
            role, note = 'PLA2_like_candidate', 'Provisional from earlier Naja BLASTp analysis'
        elif family == 'PLA2' and suffix == '006364':
            role, note = 'uncertain_PLA2_family_model', 'Short/partial earlier protein matches'
        elif r['evidence_class'] == 'venom_toxin_candidate_from_annotation':
            role = 'putative_venom_toxin'
    elif r['evidence_class'] == 'venom_toxin_candidate_from_annotation':
        role = 'putative_venom_toxin'
    elif r['classification'] == 'family_keyword':
        role = 'unresolved_family_homolog'

    if (species == 'Gloydius_ussuriensis' and rtype == 'GFF_gene_model'
            and r['overlapping_ToxCodAn_IDs'] and not anchor):
        include = 0
        note = (note + '; ' if note else '') + 'Overlaps separate ToxCodAn annotation; hidden to avoid duplicate plotting'

    if species == 'Gloydius_ussuriensis' and family == 'PLA2' and gene in {'FUN_018053', 'FUN_018054', 'FUN_018055'} and rtype == 'GFF_gene_model':
        vwa.append(r)
        include = 0
        note = 'Original Funannotate component of provisional VWA5B1 locus group'

    rows.append(dict(species=species, family=family, chr=r['chr'], start=start,
                     end=end, strand=r['strand'], gene_id=gene or r['source_record_id'],
                     source=rtype, functional_hint=r['functional_hint'],
                     anchor_symbol=anchor, plot_role=role, plot_include=include,
                     annotation_status=status, review_notes=note))

if len(vwa) != 3:
    raise SystemExit(f'ERROR: Expected three Gloydius VWA5B1 component models; found {len(vwa)}')
if len({r['chr'] for r in vwa}) != 1 or len({r['strand'] for r in vwa}) != 1:
    raise SystemExit('ERROR: VWA5B1 component models disagree on chromosome or strand')

rows.append(dict(species='Gloydius_ussuriensis', family='PLA2', chr=vwa[0]['chr'],
                 start=min(int(r['start']) for r in vwa),
                 end=max(int(r['end']) for r in vwa), strand=vwa[0]['strand'],
                 gene_id='FUN_018053;FUN_018054;FUN_018055',
                 source='provisional_group_of_three_GFF_models',
                 functional_hint='VWA5B1 candidate: complementary protein alignments',
                 anchor_symbol='VWA5B1', plot_role='provisional_VWA5B1_anchor_group',
                 plot_include=1, annotation_status='not_a_merged_gene_model',
                 review_notes='Display group as one candidate region; original 3 models retained above'))

order = {'Gloydius_ussuriensis': 0, 'Cerastes_gasperettii': 1,
         'Crotalus_adamanteus': 2, 'Naja_naja': 3}

rows.sort(key=lambda r: (r['family'], order[r['species']], r['chr'], r['start'], r['end']))

fields = ['species', 'family', 'chr', 'start', 'end', 'strand', 'gene_id', 'source',
          'functional_hint', 'anchor_symbol', 'plot_role', 'plot_include',
          'annotation_status', 'review_notes']

with OUTPUT.open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)

# 4. summarize all species
print('===== MINIMAL COMPARATIVE MICROSYNTENY GENE CATALOG =====')
for family in ('SVMP', 'PLA2'):
    for species in order:
        subset = [r for r in rows if r['family'] == family and r['species'] == species]
        print(f'{species:<23} {family:<5} records={len(subset):>3} displayed={sum(r["plot_include"] for r in subset):>3}')

print('Role counts (displayed records):')
for role, count in sorted(Counter(r['plot_role'] for r in rows if r['plot_include']).items()):
    print(f'  {role}: {count}')

print(f'Output: {OUTPUT}')
print('NOTE: Display counts are not validated toxin copy numbers.')
print('===== COMPLETE =====')
PY
```

#### Step 8: Identify conserved genomic blocks using BLASTn
Now, let's identify conserved nucleotide sequences between the SVMP and PLA2 regions of the four species.

This is distinct from the gene inventories we have already completed:
- Steps 4–7: Identify genes, candidate homologs, and their genomic positions.
- Step 8: Identify local sequence conservation that can support connecting blocks in the microsynteny figures.

In subsequent plotting, we will combine the gene coordinates with the nucleotide-alignment evidence. In the script below, we will use *Gloydius ussuriensis* as the reference, compare each of the other three species separately for SVMP and PLA2, and retain only well-supported alignment blocks for visualization. We do not need whole-genome alignment, additional toxin annotation, or other extensive analyses at this stage.

Run the script below:
```sh
#!/usr/bin/env bash

# activate conda env
source "$HOME/mendel-nas1/miniconda3/etc/profile.d/conda.sh"
conda activate funannotate

set -euo pipefail

# ============================================================
# local genomic alignments of SVMP and PLA2 regions
# ============================================================

# 1. directories
ROOT="${ROOT:-$HOME/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo}"
SYNTENY="$ROOT/synteny/toxin_microsynteny"
MANIFEST="$SYNTENY/01_reference_inventory/anchor_bounded_regions.tsv"
OUTDIR="$SYNTENY/06_local_alignments"
mkdir -p "$OUTDIR/blastn"

# 2. verify genomic intervals and prepare comparisons
for CMD in python blastn; do
    command -v "$CMD" >/dev/null || { echo "ERROR: $CMD is required" >&2; exit 1; }
done
[[ -s "$MANIFEST" ]] || { echo "ERROR: missing $MANIFEST" >&2; exit 1; }

python - "$MANIFEST" "$OUTDIR" <<'PY'
import csv
import sys
from pathlib import Path

manifest, outdir = map(Path, sys.argv[1:])
with manifest.open(newline='') as handle:
    data = list(csv.DictReader(handle, delimiter='\t'))
required = {'species', 'family', 'chromosome', 'start', 'end', 'fasta'}
if not data or not required.issubset(data[0]):
    raise SystemExit('ERROR: invalid Step 5 manifest columns')
regions = {}
for r in data:
    key = (r['family'], r['species'])
    if key in regions:
        raise SystemExit(f'ERROR: duplicate interval: {key}')
    path = Path(r['fasta'])
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f'ERROR: missing FASTA: {path}')
    sequence_id, length = None, 0
    with path.open() as handle:
        for line in handle:
            if line.startswith('>'):
                if sequence_id is not None:
                    raise SystemExit(f'ERROR: multiple sequences in {path}')
                sequence_id = line[1:].split()[0]
            else:
                length += len(line.strip())
    expected = int(r['end']) - int(r['start']) + 1
    if sequence_id is None or length != expected:
        raise SystemExit(f'ERROR: unexpected FASTA length for {key}: {length} != {expected}')
    r['sequence_id'], r['length_bp'] = sequence_id, length
    regions[key] = r

species = ('Cerastes_gasperettii', 'Crotalus_adamanteus', 'Naja_naja')
pairs = []
for family in ('SVMP', 'PLA2'):
    ref = regions.get((family, 'Gloydius_ussuriensis'))
    if ref is None:
        raise SystemExit(f'ERROR: no Gloydius {family} interval')
    for sp in species:
        query = regions.get((family, sp))
        if query is None:
            raise SystemExit(f'ERROR: missing {sp} {family} interval')
        pairs.append(dict(family=family, reference='Gloydius_ussuriensis', query=sp,
                          reference_fasta=ref['fasta'], query_fasta=query['fasta']))

with (outdir / 'alignment_pairs.tsv').open('w', newline='') as handle:
    writer = csv.DictWriter(handle, fieldnames=list(pairs[0]), delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(pairs)
with (outdir / 'interval_metadata.tsv').open('w', newline='') as handle:
    fields = ['species', 'family', 'chromosome', 'start', 'end', 'fasta', 'sequence_id', 'length_bp']
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='ignore')
    writer.writeheader(); writer.writerows(data)
print(f'Validated {len(regions)} FASTAs; prepared {len(pairs)} Gloydius-centered comparisons')
PY

# 3. run BLASTN local alignments
while IFS=$'\t' read -r FAMILY REF QUERY REF_FASTA QUERY_FASTA; do
    [[ "$FAMILY" == family ]] && continue
    NAME="${FAMILY}_${REF}__${QUERY}"
    echo "===== $NAME ====="
    blastn -task blastn -query "$QUERY_FASTA" -subject "$REF_FASTA" \
        -dust yes -word_size 11 -evalue 1e-10 \
        -outfmt '6 qseqid sseqid qlen slen qstart qend sstart send pident length nident evalue bitscore' \
        -out "$OUTDIR/blastn/${NAME}.tsv"
    echo "HSPs: $(wc -l < "$OUTDIR/blastn/${NAME}.tsv")"
done < "$OUTDIR/alignment_pairs.tsv"

# 4. export alignment evidence and provisional plotting blocks
python - "$OUTDIR" <<'PY'
import csv
import sys
from collections import Counter
from pathlib import Path

outdir = Path(sys.argv[1])
def load(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))
meta = {(r['family'], r['species']): r for r in load(outdir / 'interval_metadata.tsv')}
pairs = load(outdir / 'alignment_pairs.tsv')
fields = ['family', 'reference', 'query', 'reference_chr', 'query_chr',
          'start', 'end', 'start2', 'end2', 'strand', 'aligned_bp', 'matching_bp',
          'percent_identity', 'evalue', 'bitscore', 'reference_genomic_start',
          'reference_genomic_end', 'query_genomic_start', 'query_genomic_end']
all_rows, filtered = [], []
for pair in pairs:
    family, ref, qry = pair['family'], pair['reference'], pair['query']
    rmeta, qmeta = meta[(family, ref)], meta[(family, qry)]
    path = outdir / 'blastn' / f'{family}_{ref}__{qry}.tsv'
    counts = Counter()
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 13:
                raise ValueError(f'Malformed BLASTN line in {path}: {line[:80]}')
            qid, sid, qlen, slen, qs, qe, ss, se, identity, alen, matches, evalue, score = f
            qs, qe, ss, se, alen = map(int, (qs, qe, ss, se, alen))
            identity = float(identity)
            if qid != qmeta['sequence_id'] or sid != rmeta['sequence_id']:
                raise ValueError(f'Unexpected FASTA IDs in {path}: {qid} / {sid}')
            if int(qlen) != int(qmeta['length_bp']) or int(slen) != int(rmeta['length_bp']):
                raise ValueError(f'BLASTN length mismatch in {path}')
            a, b = sorted((ss, se)); c, d = sorted((qs, qe))
            record = dict(family=family, reference=ref, query=qry,
                          reference_chr=rmeta['chromosome'], query_chr=qmeta['chromosome'],
                          start=a, end=b, start2=c, end2=d,
                          strand='+' if (qe-qs)*(se-ss) >= 0 else '-',
                          aligned_bp=alen, matching_bp=int(matches),
                          percent_identity=identity, evalue=evalue, bitscore=score,
                          reference_genomic_start=int(rmeta['start'])+a-1,
                          reference_genomic_end=int(rmeta['start'])+b-1,
                          query_genomic_start=int(qmeta['start'])+c-1,
                          query_genomic_end=int(qmeta['start'])+d-1)
            all_rows.append(record); counts['all'] += 1
            minimum = 500 if family == 'SVMP' else 150
            if alen >= minimum and identity >= 75.0:
                filtered.append(record); counts['filtered'] += 1
    print(f'{family} {ref} vs {qry}: {counts["all"]} HSPs, {counts["filtered"]} filtered')

for name, rows in (('all_local_alignment_blocks.tsv', all_rows),
                   ('filtered_blocks_for_plot.tsv', filtered)):
    with (outdir / name).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    print(f'Output: {outdir / name} ({len(rows)} rows)')
print('NOTE: these are local HSPs, not verified orthology or nonredundant synteny blocks.')
print('===== COMPLETE =====')
PY
```

Run this as a shell script:
```sh
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Shell"

chmod +x "${ROOT}/toxin_microsynteny_local_align.sh"
bash "${ROOT}/toxin_microsynteny_local_align.sh"
```