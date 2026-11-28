## 10) Toxin gene microsynteny
In this section, we will investigate the chromosomal microsynteny and phylogenetic relationships of two major snake venom components: snake venom metalloproteinases (SVMPs) and phospholipases A2 (PLA2s). In addition to our focal species, *Gloydius ussuriensis*, we will use four comparison taxa: *G. shedaoensis*, *Crotalus adamanteus*, *Cerastes gasperettii*, and *Naja naja*.

For *G. ussuriensis*, we will use our previously curated ToxCodAn-Genome annotations. For *C. adamanteus*, *C. gasperettii*, and *N. naja*, we will leverage published toxin gene annotations and/or genomic coordinate catalogs. For *G. shedaoensis*, we will use the ToxCodAn-Genome predictions I generated, supplemented by homology-based searches and the general genome annotation.

One thing to point out for *C. adamanteus*: The assembly we used for chromosomal macrosynteny was from Hogan et al. (2024). However, Nachtigall et al. (2025) subsequently generated another chromosome-level, haplotype-resolved assembly. While these two genomes come from the same individual and the same PacBio data, the latter had higher Hi-C coverage. This is therefore preferable for our toxin gene analyses. We will use the primary assembly from Nachtigall et al. (2025).

Standardized homology searches will be conducted to assess existing gene models, identify potentially overlooked toxin homologs, and facilitate consistent toxin gene classification. For comparative microsynteny, we will primarily focus on chromosomally resolved toxin gene clusters and their conserved flanking genes.

As explained below, our preliminary investigation revealed that the SVMP gene cluster of *G. shedaoensis* is likely incompletely resolved in its current genome assembly. We will therefore exclude this species from SVMP microsynteny but retain it for PLA2 microsynteny.

### Step 0: Set up directory structure and link relevant data
The directory name we will use here is `toxin_microsynteny_v2`. This is a second version because I ended up scrapping the initial workflow I built. Run the script below to set up our directory structure and link existing data to the new workflow. This can be run directly on Mendel.
```sh
#!/bin/bash

# ================================================================
# initialize directories and link existing reference files
# ================================================================

# 1. directories
PROJECT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
ROOT="${PROJECT}/synteny/toxin_microsynteny_v2"
REF="${ROOT}/00_reference"
ASSEMBLIES="${PROJECT}/synteny/assemblies_synteny"
TOXCODAN="${PROJECT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
FUNANNOTATE="${PROJECT}/annotation/funannotate/annotate_results"
CURATED="${TOXCODAN}/all120_QC/curation/exports"

# 2. create workflow directories
mkdir -p "${ROOT}"/{00_reference,01_query,02_discovery,03_homology,04_curation,05_classification,06_microsynteny,scripts/outfiles}

for SP in G_ussuriensis G_shedaoensis C_adamanteus C_gasperettii N_naja; do
    mkdir -p "${REF}/${SP}"
done

# 3. function for safely linking existing files
link_file() {

    local SOURCE="$1"
    local DEST="$2"

    if [[ ! -f "${SOURCE}" ]]; then
        echo "MISSING: ${SOURCE}"
    elif [[ -e "${DEST}" || -L "${DEST}" ]]; then
        echo "EXISTS: ${DEST}"
    else
        ln -s "${SOURCE}" "${DEST}"
        echo "LINKED: ${DEST}"
    fi
}

# 4. G. ussuriensis
DEST="${REF}/G_ussuriensis"

link_file \
    "${PROJECT}/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa" \
    "${DEST}/genome.fa"

link_file \
    "${FUNANNOTATE}/Gloydius_ussuriensis_AMNH_21010.gff3" \
    "${DEST}/general_annotation.gff3"

for FILE in toxin_annotation.gtf toxin_annotation.tsv toxin_annotation_cds.fasta toxin_annotation_pep.fasta annotation_warning.txt annotation_removed.txt; do
    link_file "${TOXCODAN}/${FILE}" "${DEST}/${FILE}"
done

# link curated toxin annotations
for FILE in core_candidates.gtf core_candidates.cds.fasta core_candidates.pep.fasta curated_nonredundant.gtf curated_nonredundant.cds.fasta curated_nonredundant.pep.fasta curated_family_counts.tsv; do
    link_file "${CURATED}/${FILE}" "${DEST}/${FILE}"
done

# 5. G. shedaoensis
DEST="${REF}/G_shedaoensis"
SOURCE="${ASSEMBLIES}/Gloydius_shedaoensis"

link_file "${SOURCE}/Gloydius_shedaoensis.genome.fa" "${DEST}/genome.fa"
link_file "${SOURCE}/Gloydius_shedaoensis.annotation.gff3" "${DEST}/general_annotation.gff3"
link_file "${SOURCE}/Gloydius_shedaoensis.protein.faa" "${DEST}/general_proteins.faa"
link_file "${SOURCE}/Gloydius_shedaoensis.genome.fa.fai" "${DEST}/genome.fa.fai"

# 6. C. gasperettii
DEST="${REF}/C_gasperettii"
SOURCE="${ASSEMBLIES}/Cerastes_gasperettii"

link_file "${SOURCE}/Cerastes_gasperettii.genome.fa" "${DEST}/genome.fa"
link_file "${SOURCE}/Cerastes_gasperettii.annotation.gff3" "${DEST}/general_annotation.gff3"

# 7. N. naja
DEST="${REF}/N_naja"
SOURCE="${ASSEMBLIES}/Naja_naja"

link_file "${SOURCE}/Naja_naja.genome.fa" "${DEST}/genome.fa"
link_file "${SOURCE}/Naja_naja.annotation.gff3" "${DEST}/general_annotation.gff3"

# 8. inspect linked files
ls -lhR "${REF}"

echo -e "\n===== BROKEN SYMBOLIC LINKS ====="
find "${REF}" -xtype l
```

For *G. ussuriensis*, we will retain both the core and nonredundant toxin annotation sets. The core set will primarily be used to establish our initial query sequences, whereas the nonredundant set will also be examined to identify disrupted or incomplete gene models.

For the comparison species, we will retain the general genome annotations in addition to the toxin-specific annotations, as these will be needed to identify conserved flanking genes and potentially overlooked toxin homologs.

### Step 1: Data acquisition & prep
#### C. adamanteus
Run the script below to download the primary genome assembly, general genome annotation, and toxin gene annotations from Nachtigall et al. (2025).
```sh
#!/bin/bash

# 1. directories
ROOT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2/00_reference/C_adamanteus'

cd "${ROOT}"

URL='https://ndownloader.figshare.com'

# 2. file download function
download_file() {

    local DOWNLOAD_URL="$1" FILE="$2"
    local -a CHECK

    case "${FILE}" in
        *.zip) CHECK=(unzip -tq) ;;
        *.gz) CHECK=(gzip -t) ;;
    esac

    if [[ -s "${FILE}" ]] && "${CHECK[@]}" "${FILE}" > /dev/null 2>&1; then
        echo "Already downloaded and validated: ${FILE}"
        return 0
    fi

    echo "Downloading: ${FILE}"

    curl -fL --retry 2 --connect-timeout 30 "${DOWNLOAD_URL}" -o "${FILE}.part" || return 1

    [[ -s "${FILE}.part" ]] || return 1

    "${CHECK[@]}" "${FILE}.part" || return 1

    mv "${FILE}.part" "${FILE}"
}

# 3. curated toxin annotation
download_file "${URL}/articles/26526730/versions/1" "Cadam_toxin_annotation.zip"

# 4. funannotate general annotation
download_file "${URL}/articles/26526748/versions/1" "Cadam_funannotate.zip"

# 5. primary genome assembly
download_file "${URL}/files/45450430" "Cadam_primary_chromosomes.fasta.gz"

# 6. inspect downloaded files
ls -lh

unzip -l Cadam_toxin_annotation.zip
unzip -l Cadam_funannotate.zip

zcat Cadam_primary_chromosomes.fasta.gz | grep '^>'
```

Decompress the files after downloading and standardize their filenames using symbolic links.
```sh
# 1. directories
ROOT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2/00_reference/C_adamanteus'

cd "${ROOT}"

# 2. extract primary toxin annotation
unzip -n Cadam_toxin_annotation.zip Cadam_primary_toxin_annotation.gtf.gz

# 3. extract primary funannotate annotation
mapfile -t FILES < <(unzip -Z1 Cadam_funannotate.zip | grep -Ei 'primary.*\.(gff3?|gtf)(\.gz)?$')

if [[ ${#FILES[@]} -gt 0 ]]; then
    unzip -n Cadam_funannotate.zip "${FILES[@]}"
else
    echo "ERROR: Primary annotation not found"
fi

# 4. decompress primary genome and annotations
for FILE in Cadam_primary_chromosomes.fasta.gz Cadam_primary_toxin_annotation.gtf.gz "${FILES[@]}"; do

    if [[ "${FILE}" == *.gz && ! -s "${FILE%.gz}" ]]; then
        gunzip -k "${FILE}"
    fi

done

# 5. standardize reference filenames
ln -sfn Cadam_primary_chromosomes.fasta genome.fa
ln -sfn Cadam_primary_chromosomes.funannotate.gff general_annotation.gff3
ln -sfn Cadam_primary_toxin_annotation.gtf toxin_annotation.gtf

# 6. inspect files
ls -lh
grep '^>' genome.fa
head -n 10 toxin_annotation.gtf

find . -xtype l
```

One thing to keep in mind is that Nachtigall et al. (2025) generated haplotype-resolved assemblies and identified structural variation within the SVMP cluster. Our microsynteny analysis will focus on the primary assembly, and the resulting toxin gene counts should therefore be interpreted as those of the selected reference assembly.

#### C. gasperettii
Run the script below to obtain supplementary data from the original genome assembly paper by Mochales-Riaño et al. (2025).

```sh
ROOT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2'

cd "${ROOT}/00_reference/C_gasperettii"

URL='https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102647'

# 1. curated genome annotation
wget -c "${URL}/annotation/curated_annotation.gff"

# 2. coding sequences
wget -c "${URL}/annotation/Cgasperettii.cds"

# 3. protein sequences
wget -c "${URL}/annotation/coding_proteins.fa"

# 4. svmp supplementary data
wget -c "${URL}/alignment/SVMP/SVMP.csv"

# 5. pla2 supplementary data
wget -c "${URL}/alignment/PLA/PLA2.csv"

# 6. inspect files
ls -lh
```

After examining the downloaded annotation, I confirmed that the annotation we previously used for chromosomal macrosynteny corresponds to the same curated genome annotation. The previously acquired annotation also has sequence identifiers that are compatible with our local genome assembly.

We will therefore continue to use the existing `general_annotation.gff3` linked in Step 0 rather than replace it with the newly downloaded `curated_annotation.gff`. The supplementary files downloaded here will be retained as reference material for toxin gene identification and classification.

#### N. naja
Use the script below to download the toxin gene catalog from Suryamohan et al. (2020).

```sh
ROOT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2'

cd "${ROOT}/00_reference/N_naja"

wget -c \
    'https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-019-0559-8/MediaObjects/41588_2019_559_MOESM3_ESM.xlsx' \
    -O Suryamohan_2020_supplementary_tables.xlsx

ls -lh Suryamohan_2020_supplementary_tables.xlsx

file Suryamohan_2020_supplementary_tables.xlsx
```

We will primarily use Supplementary Tables 6b and 6c to identify the published SVMP and PLA2 gene models and their associated genomic coordinates.

#### G. shedaoensis
Since we do not have a separately curated toxin gene annotation for this species, I initially ran ToxCodAn-Genome using the entire genome assembly of *G. shedaoensis*.

However, the genome contains numerous unplaced scaffolds in addition to explicitly designated chromosomes. So I reran ToxCodAn-Genome using only the chromosomal sequences to determine which toxin gene models are placed on these chromosomes.

The chromosome-only genome assembly and the corresponding ToxCodAn-Genome predictions are stored in:
```sh
/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_Gshe
```

Run the script below to link these files to our current workflow.
```sh
# 1. directories
PROJECT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
REF="${PROJECT}/synteny/toxin_microsynteny_v2/00_reference"
GSHE="${PROJECT}/annotation/toxin_gene_annotation_Gshe"

cd "${REF}/G_shedaoensis"

# 2. link chromosome-only genome
ln -sfn "${GSHE}/Gloydius_shedaoensis.chromosome_only.genome.fa" chromosome_only.genome.fa

# 3. link chromosome-only toxcodan annotations
for FILE in toxin_annotation.gtf toxin_annotation_cds.fasta toxin_annotation_pep.fasta annotation_warning.txt annotation_removed.txt; do

    if [[ -f "${GSHE}/${FILE}" ]]; then
        ln -sfn "${GSHE}/${FILE}" "${FILE}"
    else
        echo "MISSING: ${GSHE}/${FILE}"
    fi

done

# 4. inspect linked files
ls -lh

find . -xtype l
```
Note that `genome.fa` refers to the original full genome assembly, whereas `chromosome_only.genome.fa` contains only explicitly designated chromosomes. We will retain both files.

#### Verify assembly-annotation compatibility across species
Before conducting homology searches or extracting toxin gene models, we need to make sure that the genome assemblies and their associated annotations use the same sequence identifiers.

The script below checks for annotation sequence IDs absent from the corresponding genome, annotation coordinates extending beyond genomic sequence boundaries, and invalid genomic coordinates.
```sh
# 1. directories
REF='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2/00_reference'

cd "${REF}"

# 2. verify assembly-annotation compatibility
python3 - <<'PY'
from pathlib import Path

species = {
    "G_ussuriensis": ["general_annotation.gff3", "core_candidates.gtf", "curated_nonredundant.gtf"],
    "G_shedaoensis": ["general_annotation.gff3"],
    "C_adamanteus": ["general_annotation.gff3", "toxin_annotation.gtf"],
    "C_gasperettii": ["general_annotation.gff3"],
    "N_naja": ["general_annotation.gff3"]
}

for sp, annotations in species.items():

    print(f"\n===== {sp} =====", flush=True)

    lengths = {}

    with open(Path(sp) / "genome.fa", "rb") as f:

        for line in f:

            if line.startswith(b">"):
                seqid = line[1:].split()[0].decode()
                lengths[seqid] = 0
            else:
                lengths[seqid] += len(line.strip())

    print("Genome sequences:", len(lengths))

    for annotation in annotations:

        ids, missing = set(), set()
        total = invalid = outside = 0

        with open(Path(sp) / annotation) as f:

            for line in f:

                if line.startswith("##FASTA"):
                    break

                if line.startswith("#") or not line.strip():
                    continue

                fields = line.rstrip().split("\t")

                if len(fields) != 9:
                    continue

                seqid = fields[0]
                start, end = int(fields[3]), int(fields[4])

                total += 1
                ids.add(seqid)

                if start < 1 or end < start:
                    invalid += 1

                if seqid not in lengths:
                    missing.add(seqid)
                elif end > lengths[seqid]:
                    outside += 1

        print("\nAnnotation:", annotation)
        print("Annotation records:", total)
        print("Annotation sequence IDs:", len(ids))
        print("Unmatched sequence IDs:", len(missing))
        print("Out-of-bound records:", outside)
        print("Invalid coordinates:", invalid)
        print("Status:", "PASS" if total and not (missing or outside or invalid) else "REVIEW")

PY
```

All eight annotation files examined here passed these compatibility checks.

| Species | Genome sequences | Unmatched sequence IDs | Out-of-bound records | Invalid coordinates |
|---|---:|---:|---:|---:|
| G. ussuriensis | 19 | 0 | 0 | 0 |
| G. shedaoensis | 710 | 0 | 0 | 0 |
| C. adamanteus | 19 | 0 | 0 | 0 |
| C. gasperettii | 222 | 0 | 0 | 0 |
| N. naja | 1,897 | 0 | 0 | 0 |

This confirms that the locally available genome assemblies and annotations are compatible in terms of sequence identifiers and genomic coordinate bounds. Note that this does not necessarily mean that toxin gene models from separate annotation catalogs correspond to one another.

#### Validate and inventory chromosome-only G. shedaoensis toxin annotations
Since the new ToxCodAn-Genome predictions for *G. shedaoensis* were generated using only explicitly designated chromosomes, we will conduct an additional compatibility check to make sure that the retained chromosomal sequences are identical to their counterparts in the full genome assembly.

We will also inspect the recovered SVMP and PLA2 gene models, including their CDS and protein lengths, genomic coordinates, start and stop codons, internal stop codons, translation consistency, and annotation warnings.
```sh
# 0. activate conda env
conda activate funannotate

# 1. directories
REF='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2/00_reference'

cd "${REF}/G_shedaoensis"

# 2. verify chromosome-only genome and toxin annotations
python - <<'PY'
import re
import hashlib
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq

# 2.1. obtain fasta sequence lengths and hashes
def fasta_stats(path, targets=None):

    records = {}
    seqid = None
    length = 0
    digest = None

    with open(path, "rb") as f:

        for line in f:

            if line.startswith(b">"):

                if digest is not None:
                    records[seqid] = (length, digest.hexdigest())

                seqid = line[1:].split()[0].decode()
                length = 0
                digest = hashlib.sha256() if targets is None or seqid in targets else None

            elif digest is not None:

                seq = b"".join(line.split()).upper()
                length += len(seq)
                digest.update(seq)

    if digest is not None:
        records[seqid] = (length, digest.hexdigest())

    return records

chromosomes = fasta_stats("chromosome_only.genome.fa")
full = fasta_stats("genome.fa", set(chromosomes))

missing = set(chromosomes) - set(full)
different = [x for x in chromosomes if x in full and chromosomes[x] != full[x]]

print("\n===== CHROMOSOME-ONLY GENOME =====")
print("Retained chromosomes:", len(chromosomes))
print("Missing sequences:", len(missing))
print("Different sequences:", len(different))

for seqid, (length, digest) in chromosomes.items():
    print(seqid, length, sep="\t")

if not chromosomes or missing or different:
    raise SystemExit("Genome compatibility check failed.")

# 2.2. verify toxin annotation coordinates
models = defaultdict(list)
total = invalid = outside = 0
unmatched = set()
annotated = set()

with open("toxin_annotation.gtf") as f:

    for line in f:

        if line.startswith("#") or not line.strip():
            continue

        r = line.rstrip().split("\t")

        if len(r) != 9:
            continue

        seqid = r[0]
        start, end = int(r[3]), int(r[4])

        total += 1
        annotated.add(seqid)

        if start < 1 or end < start:
            invalid += 1

        if seqid not in chromosomes:
            unmatched.add(seqid)
        elif end > chromosomes[seqid][0]:
            outside += 1

        if r[2] == "CDS":

            match = re.search(r'gene_id\s+"((SVMP|PLA2)-[^"]+)"', r[8])

            if match:
                models[match.group(1)].append(r)

print("\n===== ANNOTATION COMPATIBILITY =====")
print("Annotation records:", total)
print("Annotated chromosomes:", len(annotated))
print("Unmatched sequence IDs:", len(unmatched))
print("Out-of-bound records:", outside)
print("Invalid coordinates:", invalid)

if not total or unmatched or outside or invalid:
    raise SystemExit("Annotation compatibility check failed.")

# 3. read cds and protein sequences
cds = SeqIO.to_dict(SeqIO.parse("toxin_annotation_cds.fasta", "fasta"))
pep = SeqIO.to_dict(SeqIO.parse("toxin_annotation_pep.fasta", "fasta"))

# 4. inspect svmp and pla2 gene models
print("\n===== SVMP AND PLA2 INVENTORY =====")
print("Gene_ID\tChromosome\tStart\tEnd\tStrand\tCDS_segments\tGTF_bp\tFASTA_bp\tProtein_aa\tStart_codon\tStop_codon\tInternal_stops\tTranslation_match")

signatures = defaultdict(list)

for gene, records in sorted(models.items()):

    chrom = records[0][0]
    strand = records[0][6]

    start = min(int(r[3]) for r in records)
    end = max(int(r[4]) for r in records)
    length = sum(int(r[4]) - int(r[3]) + 1 for r in records)

    cds_seq = str(cds[gene].seq).upper() if gene in cds else ""
    pep_seq = str(pep[gene].seq) if gene in pep else ""

    signature = tuple(sorted((r[0], int(r[3]), int(r[4]), r[6]) for r in records))
    signatures[signature].append(gene)

    if cds_seq and len(cds_seq) % 3 == 0:

        translated = str(Seq(cds_seq).translate())
        internal = translated[:-1].count("*") if translated.endswith("*") else translated.count("*")
        match = translated.rstrip("*") == pep_seq.rstrip("*") if pep_seq else "MISSING"

    else:

        internal = "NA"
        match = "FRAME_ERROR" if cds_seq else "MISSING"

    print(gene, chrom, start, end, strand, len(records), length, len(cds_seq), len(pep_seq.rstrip("*")), cds_seq[:3], cds_seq[-3:], internal, match, sep="\t")

# 5. inspect redundant annotations
print("\n===== IDENTICAL CDS COORDINATE MODELS =====")

duplicates = [genes for genes in signatures.values() if len(genes) > 1]

for genes in duplicates:
    print(", ".join(genes))

if not duplicates:
    print("None")

# 6. summarize gene counts
print("\n===== ANNOTATION COUNTS =====")
print(dict(Counter(g.split("-")[0] for g in models)))

PY

# 7. inspect annotation warnings
grep -iE 'SVMP|PLA2' annotation_warning.txt || echo "No SVMP or PLA2 warnings"

# 8. inspect removed annotations
grep '^#' annotation_removed.txt | grep -iE 'SVMP|PLA2' || echo "No SVMP or PLA2 removal headers"
```

The chromosome-only genome retained 19 chromosomal sequences, all of which were identical to their counterparts in the full genome assembly.

The new ToxCodAn-Genome annotation contained 277 records distributed across 14 chromosomes. There were no unmatched sequence identifiers or invalid genomic coordinates.

However, only one SVMP and two PLA2 gene models were recovered from this chromosome-only run.

| Gene ID | Chromosome | Start | End | Strand | CDS length | Protein length |
|---|---|---:|---:|:---:|---:|---:|
| SVMP-1 | GWHBWDU00000013 | 16828413 | 16852744 | − | 1839 | 612 |
| PLA2-1 | GWHBWDU00000018 | 7771779 | 7773153 | − | 435 | 144 |
| PLA2-2 | GWHBWDU00000018 | 7776041 | 7777615 | + | 417 | 138 |

All three gene models contained canonical start and stop codons, lacked internal stop codons, and exhibited consistent CDS-protein translations. No coordinate-identical gene models or SVMP/PLA2 annotation warnings were detected.

The initial ToxCodAn-Genome run using the entire genome assembly recovered seven SVMP annotations and two PLA2 annotations. However, SVMP-4 and SVMP-5 from this initial run had identical genomic boundaries, strands, CDS segment counts, and CDS lengths, suggesting that these may represent redundant annotations. More importantly, the seven SVMP annotations were distributed across three assembly sequences, two of which were unplaced scaffolds (GWHBWDU00000186 and GWHBWDU00000187).

The results from the two ToxCodAn-Genome runs can be summarized like so:
| Input genome | SVMP annotations | Provisional SVMP loci | PLA2 annotations |
|---|---:|---:|---:|
| Full genome assembly | 7 | 6 | 2 |
| Chromosome-only genome | 1 | 1 | 2 |

This means that the smaller number of chromosomally resolved SVMPs could reflect assembly fragmentation or incomplete annotation rather than genuine biological contraction of the SVMP repertoire.

#### Preliminary homology searches for SVMP and PLA2 in G. shedaoensis
Given the unexpectedly small number of recovered SVMP gene models, I conducted preliminary tBLASTn searches to determine whether additional SVMP homologs occur on the chromosome containing SVMP-1. I also conducted the same analysis for the chromosome containing the two annotated PLA2 genes.

For these analyses, I used the 17 curated SVMP proteins and five curated venom PLA2 proteins from *G. ussuriensis* as queries. Note that PLA2-1 from *G. ussuriensis* was excluded because we previously reclassified this gene as a PLA2 inhibitor-like protein.

The two preliminary searches can be consolidated into the script below. This script does not generate any additional files and can be run directly on Mendel.
```sh
# 0. activate conda env
conda activate funannotate

(
set -euo pipefail

# 1. directories
PROJECT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
REF="${PROJECT}/synteny/toxin_microsynteny_v2/00_reference"

GENOME="${REF}/G_shedaoensis/genome.fa"
PROTEINS="${REF}/G_ussuriensis/core_candidates.pep.fasta"

[[ -s "${GENOME}.fai" ]] || { echo "ERROR: Genome FASTA index missing"; exit 1; }

# 2. conduct preliminary homology searches
for FAMILY in SVMP PLA2; do

    if [[ "${FAMILY}" == "SVMP" ]]; then

        CHR='GWHBWDU00000013'
        PATTERN='^>SVMP-'
        EXPECTED=17

    else

        CHR='GWHBWDU00000018'
        PATTERN='^>PLA2-(2|3|5|6|7)$'
        EXPECTED=5

    fi

    QUERY=$(awk -v pattern="${PATTERN}" '/^>/{keep=($1 ~ pattern)} keep' "${PROTEINS}")

    NQUERY=$(printf '%s\n' "${QUERY}" | grep -c '^>')

    echo -e "\n===== ${FAMILY} ====="
    echo "Target chromosome: ${CHR}"
    echo "Number of query proteins: ${NQUERY}"

    [[ "${NQUERY}" -eq "${EXPECTED}" ]] || { echo "Unexpected query count"; exit 1; }

    # 3. run tblastn
    printf '%s\n' "${QUERY}" | tblastn \
        -query - \
        -subject <(samtools faidx "${GENOME}" "${CHR}") \
        -evalue 1e-5 \
        -seg yes \
        -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen' | \
    python -c '

import sys
from collections import defaultdict

family = sys.argv[1]

known = {
    "SVMP": {
        "SVMP-1": (16828413, 16852744, "-")
    },
    "PLA2": {
        "PLA2-1": (7771779, 7773153, "-"),
        "PLA2-2": (7776041, 7777615, "+")
    }
}[family]

regions = defaultdict(list)
total = 0

# 4. parse blast hits
for line in sys.stdin:

    r = line.rstrip().split("\t")

    query = r[0]
    identity = float(r[2])
    length = int(r[3])
    qstart, qend = int(r[4]), int(r[5])
    sstart, send = int(r[6]), int(r[7])
    evalue = float(r[8])
    bitscore = float(r[9])
    qlen = int(r[10])

    start, end = sorted([sstart, send])
    strand = "+" if sstart <= send else "-"

    status = "Outside_known_loci"

    for gene, (left, right, direction) in known.items():

        if start <= right and end >= left:

            status = gene if strand == direction else gene + "_opposite_strand"
            break

    regions[(strand, status)].append(
        (start, end, query, identity, length, evalue, bitscore, qstart, qend, qlen)
    )

    total += 1

# 5. group nearby hits
print("\n===== GENOMIC HIT REGIONS =====")
print("Start\tEnd\tStrand\tQueries\tHSPs\tIdentity_pct\tAlignment_aa\tEvalue\tStatus")

results = []

for (strand, status), hits in regions.items():

    groups = []

    for hit in sorted(hits):

        if not groups or hit[0] > groups[-1]["end"] + 1000:
            groups.append({"end": hit[1], "hits": [hit]})
        else:
            groups[-1]["hits"].append(hit)
            groups[-1]["end"] = max(groups[-1]["end"], hit[1])

    for group in groups:

        records = group["hits"]

        start = min(r[0] for r in records)
        end = max(r[1] for r in records)
        queries = len({r[2] for r in records})
        best = max(records, key=lambda r: r[6])

        results.append(
            (start, end, strand, queries, len(records), best[3], best[4], best[5], status)
        )

for result in sorted(results):
    print("\t".join(map(str, result)))

# 6. inspect individual hits outside known loci
print("\n===== INDIVIDUAL HITS OUTSIDE KNOWN LOCI =====")
print("Start\tEnd\tStrand\tQuery\tIdentity_pct\tAlignment_aa\tQuery_start\tQuery_end\tQuery_length\tEvalue")

for (strand, status), hits in sorted(regions.items()):

    if status != "Outside_known_loci":
        continue

    for r in sorted(hits):

        print(r[0], r[1], strand, r[2], r[3], r[4], r[7], r[8], r[9], r[5], sep="\t")

# 7. summarize results
print("\n===== SUMMARY =====")
print("Total BLAST HSPs:", total)
print("Total genomic hit groups:", len(results))

for gene in known:
    print("Hit groups overlapping", gene, ":", sum(r[-1] == gene for r in results))

print("Hit groups outside known loci:", sum(r[-1] == "Outside_known_loci" for r in results))

' "${FAMILY}"

done

)
```
For the purpose of this preliminary investigation, nearby BLAST hits within 1 kb of one another were grouped together. Note that these groups do not necessarily correspond to individual gene loci, particularly for SVMPs, whose CDS segments may be separated by relatively long introns.

##### SVMP
The tBLASTn search recovered 170 HSPs distributed across 13 genomic hit groups. Nine of these groups overlapped the known SVMP-1 gene model, whereas four groups occurred outside this annotation.

| Start | End | Strand | Identity (%) | Alignment (aa) | E-value |
|---:|---:|:---:|---:|---:|---:|
| 7464817 | 7465518 | − | 26.606 | 218 | 3.05e-10 |
| 7493006 | 7493548 | + | 28.099 | 121 | 6.53e-07 |
| 16877633 | 16877803 | − | 71.739 | 46 | 1.72e-13 |
| 16879839 | 16879958 | − | 72.500 | 40 | 1.70e-12 |

The first two hits exhibited relatively low sequence similarity to our SVMP queries and were located far away from the known SVMP locus. However, the other two hits were located relatively close to SVMP-1 and exhibited high sequence similarity to *G. ussuriensis* SVMP proteins. These two hits occurred on the same strand and were located near the chromosome boundary.

Interestingly, both hits corresponded to the C-terminal region of *G. ussuriensis* SVMP-20.
| Query | Query start | Query end | Identity (%) | E-value |
|---|---:|---:|---:|---:|
| SVMP-20 | 514 | 553 | 72.500 | 1.70e-12 |
| SVMP-20 | 553 | 598 | 71.739 | 1.72e-13 |

The order of these two alignments was consistent with the negative strand orientation of a potential additional SVMP gene model.

##### PLA2
The tBLASTn search recovered 49 HSPs distributed across six genomic hit groups. Two of these groups corresponded to the two ToxCodAn-Genome PLA2 predictions, whereas four additional groups occurred outside these annotations.

| Start | End | Strand | Identity (%) | Alignment (aa) | E-value |
|---:|---:|:---:|---:|---:|---:|
| 5228152 | 5228289 | − | 45.652 | 46 | 1.30e-06 |
| 7766366 | 7767083 | − | 39.189 | 74 | 4.07e-08 |
| 7771779 | 7772882 | − | 72.500 | 40 | 2.02e-15 |
| 7776311 | 7777615 | + | 72.917 | 48 | 7.17e-19 |
| 7779670 | 7780517 | − | 83.333 | 60 | 1.10e-27 |
| 7784012 | 7784843 | + | 25.781 | 128 | 5.55e-08 |

The identities, alignment lengths, and E-values reported in the tables above correspond to the highest-bitscore HSP within each hit group.

Five of the six PLA2-associated hit groups were located within an approximately 18.5 kb genomic region. This suggests that the chromosome contains additional PLA2 homologs that were not recovered by ToxCodAn-Genome. The most interesting additional candidate was located at 7779670–7780517 bp. Three alignments corresponding to this region collectively covered the entire 119-aa *G. ussuriensis* PLA2-2 query.

| Query start | Query end | Identity (%) | E-value |
|---:|---:|---:|---:|
| 1 | 35 | 97.143 | 1.22e-17 |
| 29 | 74 | 80.435 | 4.07e-19 |
| 60 | 119 | 83.333 | 1.10e-27 |

This provides strong evidence for the presence of a third PLA2 homolog within the same chromosomal cluster. Note, however, that our PLA2-2 query is shorter than the other curated PLA2 proteins. Complete query coverage therefore does not necessarily indicate a complete PLA2 gene model.

#### Reconcile additional G. shedaoensis homologs with general genome annotations
Next, I examined the general genome annotation of *G. shedaoensis* to determine whether the additional SVMP and PLA2 homologs identified by tBLASTn correspond to previously annotated gene models. I also examined the genomic regions surrounding the SVMP and PLA2 clusters for assembly gaps represented by Ns.
```sh
# 0. activate conda env
conda activate funannotate

# 1. directories
PROJECT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
REF="${PROJECT}/synteny/toxin_microsynteny_v2/00_reference"

GENOME="${REF}/G_shedaoensis/genome.fa"
GFF="${REF}/G_shedaoensis/general_annotation.gff3"

# 2. inspect gene models near the svmp chromosome boundary
echo -e "\n===== SVMP REGION ====="

awk -F '\t' '$1 == "GWHBWDU00000013" && $4 <= 16880701 && $5 >= 16852745 && ($3 == "gene" || $3 == "mRNA" || $3 == "CDS" || $3 == "pseudogene")' "${GFF}"

# 3. inspect gene models within the pla2 cluster
echo -e "\n===== PLA2 REGION ====="

awk -F '\t' '$1 == "GWHBWDU00000018" && $4 <= 7790000 && $5 >= 7760000 && ($3 == "gene" || $3 == "mRNA" || $3 == "CDS" || $3 == "pseudogene")' "${GFF}"

# 4. inspect surrounding pla2 genes
echo -e "\n===== PLA2 FLANKING GENES ====="

awk -F '\t' '$1 == "GWHBWDU00000018" && $4 <= 7850000 && $5 >= 7700000 && ($3 == "gene" || $3 == "pseudogene")' "${GFF}" | sort -k4,4n

# 5. inspect assembly gaps
for REGION in GWHBWDU00000013:16852745-16880701 GWHBWDU00000018:7700000-7850000; do

    echo -e "\n===== ${REGION} ====="

    samtools faidx "${GENOME}" "${REGION}" | python -c '

import sys
import re

region = sys.argv[1]
offset = int(region.split(":")[1].split("-")[0]) - 1

seq = "".join(line.strip() for line in sys.stdin if not line.startswith(">")).upper()

print("Region length:", len(seq))
print("N bases:", seq.count("N"))
print("N percentage:", round(100 * seq.count("N") / len(seq), 4) if seq else "NA")

print("Gap_start\tGap_end\tLength")

for match in re.finditer("N+", seq):
    print(offset + match.start() + 1, offset + match.end(), len(match.group()), sep="\t")

' "${REGION}"

done
```

##### SVMP
The two high-identity SVMP hits near the chromosome boundary corresponded to an existing Liftoff gene model, Gshe_Contig7_G00267.

| Attribute | Result |
|---|---|
| Chromosome | GWHBWDU00000013 |
| Start | 16877027 |
| End | 16879975 |
| Strand | − |
| CDS segments | 3 |
| CDS length | 317 bp |
| Valid ORF | False |
| Missing start codon | True |
| Distance from chromosome boundary | 726 bp |

The two tBLASTn hits overlapped two of the three CDS segments of this gene model, providing additional support for its identity as a partial SVMP homolog. One thing to point out here is that this gene model is located only 726 bp from the chromosome boundary. Given its negative strand orientation, the missing N-terminal coding region would be expected toward the end of the assembled chromosome.

While no Ns were detected within the region between SVMP-1 and the chromosome boundary, this does not necessarily mean that the entire SVMP cluster is assembled. Rather, the assembly may have terminated within the SVMP cluster. This interpretation is also consistent with the recovery of additional SVMP gene models on two unplaced scaffolds in our initial ToxCodAn-Genome run.

We cannot establish whether these scaffolds are physical extensions of this chromosome based on the current results. However, it seems likely that the SVMP gene cluster is incompletely resolved in this genome assembly. Therefore, I decided to exclude *G. shedaoensis* from SVMP microsynteny. The partial SVMP gene model will not be classified as a pseudogene, as its incomplete coding sequence may be a consequence of assembly fragmentation rather than biological pseudogenization.

##### PLA2
All three additional homologous regions near the annotated PLA2 genes corresponded to existing Liftoff gene models. Interestingly, while the second ToxCodAn-Genome PLA2 prediction corresponded to an existing Liftoff gene model (Gshe_Contig255_G00047), the first prediction was not represented in the general genome annotation.

The local PLA2-associated gene models can be summarized like so:

| Gene ID | Start | End | Strand | CDS segments |
|---|---:|---:|:---:|---:|
| Gshe_Contig255_G00046 | 7766324 | 7768353 | − | 4 |
| PLA2-1 | 7771779 | 7773153 | − | 4 |
| PLA2-2 / Gshe_Contig255_G00047 | 7776041 | 7777615 | + | 4 |
| Gshe_Contig255_G00048 | 7779670 | 7781071 | − | 4 |
| Gshe_Contig255_G00049 | 7783517 | 7788912 | + | 9 |

All four Liftoff gene models were annotated as having valid ORFs. Furthermore, no Ns were detected within the approximately 150 kb region surrounding the PLA2 cluster (7700000–7850000 bp).

Unlike the SVMP cluster, the PLA2 cluster does not exhibit obvious evidence of truncation or unresolved sequence gaps. This does not necessarily rule out other assembly errors, including collapsed duplications.

#### Validate additional PLA2 gene models in G. shedaoensis
To independently verify the coding sequences of the PLA2-associated Liftoff gene models, I reconstructed their CDS sequences from the genome assembly and translated them into protein sequences.

Run the script below to reconstruct the four gene models. This script will print the reconstructed proteins and their associated ORF assessments without generating any additional files.
```sh
# 0. activate conda env
conda activate funannotate

# 1. directories
PROJECT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
REF="${PROJECT}/synteny/toxin_microsynteny_v2/00_reference"

GENOME="${REF}/G_shedaoensis/genome.fa"
GFF="${REF}/G_shedaoensis/general_annotation.gff3"

# 2. reconstruct cds and translate proteins
python - "${GENOME}" "${GFF}" <<'PY'

import sys
import subprocess
from collections import defaultdict
from Bio.Seq import Seq

genome, gff = sys.argv[1:3]

chrid = "GWHBWDU00000018"

targets = [
    "Gshe_Contig255_G00046.t1",
    "Gshe_Contig255_G00047.t1",
    "Gshe_Contig255_G00048.t1",
    "Gshe_Contig255_G00049.t1"
]

# 3. extract cds coordinates
models = defaultdict(list)

with open(gff) as f:

    for line in f:

        if line.startswith("##FASTA"):
            break

        if line.startswith("#"):
            continue

        r = line.rstrip().split("\t")

        if len(r) != 9 or r[0] != chrid or r[2] != "CDS":
            continue

        attrs = dict(x.split("=", 1) for x in r[8].split(";") if "=" in x)

        for parent in attrs.get("Parent", "").split(","):

            if parent in targets:
                models[parent].append((int(r[3]), int(r[4]), r[6]))

# 4. reconstruct genomic cds
for gene in targets:

    records = sorted(models[gene])

    if not records:
        print("\nMissing gene model:", gene)
        continue

    strands = {r[2] for r in records}

    if len(strands) != 1:
        print("\nInconsistent strand:", gene)
        continue

    strand = records[0][2]
    fragments = []

    for start, end, direction in records:

        fasta = subprocess.check_output(
            ["samtools", "faidx", genome, f"{chrid}:{start}-{end}"],
            text=True
        )

        seq = "".join(x.strip() for x in fasta.splitlines() if not x.startswith(">"))

        fragments.append(seq)

    cds = Seq("".join(fragments).upper())

    if strand == "-":
        cds = cds.reverse_complement()

    # 5. inspect open reading frame
    print("\n=====", gene, "=====")
    print("Strand:", strand)
    print("CDS segments:", len(records))
    print("CDS length:", len(cds))
    print("Start codon:", cds[:3])
    print("Stop codon:", cds[-3:])
    print("Divisible by 3:", len(cds) % 3 == 0)

    if len(cds) % 3 != 0:
        print("WARNING: CDS frame error")
        continue

    protein = str(cds.translate())

    terminal_stop = protein.endswith("*")
    internal_stops = protein[:-1].count("*") if terminal_stop else protein.count("*")

    intact = str(cds[:3]) == "ATG" and terminal_stop and internal_stops == 0

    print("Protein length:", len(protein.rstrip("*")))
    print("Internal stops:", internal_stops)
    print("Terminal stop:", terminal_stop)
    print("Intact ORF:", intact)

    print("\nProtein sequence:")
    print(">" + gene)
    print(protein.rstrip("*"))

PY
```

All four reconstructed gene models contained intact ORFs, consistent with the ORF assessments in the original Liftoff annotation.

| Gene ID (abbreviated) | CDS segments | CDS length | Protein length | Start codon | Stop codon | Internal stops |
|---|---:|---:|---:|---|---|---:|
| G00046 | 4 | 471 | 156 | ATG | TGA | 0 |
| G00047 | 4 | 417 | 138 | ATG | TAA | 0 |
| G00048 | 4 | 438 | 145 | ATG | TGA | 0 |
| G00049 | 9 | 1554 | 517 | ATG | TGA | 0 |

Below are some observations regarding these gene models.

1) Gshe_Contig255_G00046: This gene encodes a 156-aa protein containing several characteristic features of secretory PLA2 proteins. However, its relatively low sequence similarity to our venom PLA2 queries suggests that it may represent a divergent PLA2 family member, potentially a nonvenom homolog. Its identity will be examined during subsequent toxin gene classification.

2) Gshe_Contig255_G00047: The genomic coordinates, strand, and CDS length of this gene model correspond to those of ToxCodAn-Genome PLA2-2. These two annotations will therefore be treated as a single genomic locus.

3) Gshe_Contig255_G00048: This gene encodes a 145-aa protein and exhibits high sequence similarity to our *G. ussuriensis* PLA2 proteins. Its reconstructed coding sequence also contains an intact ORF. We will retain this gene as an additional putative venom PLA2 candidate, pending subsequent classification.

4) Gshe_Contig255_G00049: Unlike the other PLA2-associated gene models, this gene encodes a considerably longer protein (517 aa). While its N-terminal region exhibits characteristics of a secretory PLA2 protein, the remaining sequence is substantially different. This could reflect an incorrectly merged gene model, a genuine gene fusion, or another protein architecture. We will therefore flag this gene for subsequent curation rather than include the entire protein sequence in our PLA2 phylogenetic analysis.

The distant PLA2-associated BLAST hit at approximately 5.23 Mb remains unclassified and will be examined during the standardized homology searches.

So, the preliminary PLA2 inventory for *G. shedaoensis* can be summarized like so:

| Gene ID | Current classification |
|---|---|
| Gshe_Contig255_G00046 | Intact divergent PLA2 homolog |
| PLA2-1 | Intact putative venom PLA2 |
| PLA2-2 / Gshe_Contig255_G00047 | Intact putative venom PLA2 |
| Gshe_Contig255_G00048 | Additional intact putative venom PLA2 |
| Gshe_Contig255_G00049 | PLA2-associated gene model requiring curation |

This gives us three putative venom PLA2 gene candidates, one divergent PLA2 homolog, and one potentially fused gene model within the same chromosomal region.

Unlike the SVMP gene cluster, the PLA2 gene cluster of *G. shedaoensis* appears to be reasonably well resolved. We will retain *G. shedaoensis* for comparative PLA2 microsynteny, subject to identification of conserved flanking genes and final gene model validation.

#### Inspect toxin annotation formats and inventory existing toxin gene models
Before extracting SVMP and PLA2 gene models, we need to determine how each annotation identifies toxin genes. The script below will provide the necessary information to establish an extraction strategy for each species. Since we already examined the toxin gene models of *G. shedaoensis* above, this species is not included in the script below.

```sh
# 0. activate conda env
conda activate funannotate

# 1. directories
REF='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/toxin_microsynteny_v2/00_reference'

cd "${REF}"

# 2. inspect annotation formats
echo -e "\n===== G. USSURIENSIS ====="

cat G_ussuriensis/curated_family_counts.tsv

head -n 5 G_ussuriensis/core_candidates.gtf

grep -m 5 '^>' G_ussuriensis/core_candidates.cds.fasta
grep -m 5 '^>' G_ussuriensis/core_candidates.pep.fasta

echo -e "\n===== C. ADAMANTEUS ====="

head -n 5 C_adamanteus/toxin_annotation.gtf

grep -im 5 'SVMP' C_adamanteus/toxin_annotation.gtf
grep -im 5 'PLA2' C_adamanteus/toxin_annotation.gtf

echo -e "\n===== C. GASPERETTII ====="

head -n 5 C_gasperettii/general_annotation.gff3

grep -im 10 -E 'SVMP|metalloproteinase|ADAM28' C_gasperettii/general_annotation.gff3

grep -im 10 -E 'PLA2|phospholipase' C_gasperettii/general_annotation.gff3

echo -e "\n===== N. NAJA ====="

head -n 5 N_naja/general_annotation.gff3

# 3. inventory toxin gene models
python - <<'PY'

import re
from collections import defaultdict, Counter
from openpyxl import load_workbook

annotations = {
    "G_ussuriensis_core": "G_ussuriensis/core_candidates.gtf",
    "G_ussuriensis_nonredundant": "G_ussuriensis/curated_nonredundant.gtf",
    "C_adamanteus": "C_adamanteus/toxin_annotation.gtf",
    "C_gasperettii": "C_gasperettii/general_annotation.gff3"
}

for species, annotation in annotations.items():

    loci = defaultdict(list)

    with open(annotation) as f:

        for line in f:

            if line.startswith("##FASTA"):
                break

            if line.startswith("#"):
                continue

            r = line.rstrip().split("\t")

            if len(r) != 9 or r[2] != "CDS":
                continue

            if species == "C_gasperettii":

                match = re.search(r'(?:^|;)Parent=([^;]+)', r[8])

                if not match:
                    continue

                for parent in match.group(1).split(","):

                    if re.search(r'(?:^|[-_])(SVMP|PLA2)(?:[-_]|$)', parent, re.I):
                        loci[parent].append(r)

            else:

                match = re.search(r'gene_id\s+"((SVMP|PLA2)-[^"]+)"', r[8])

                if match:
                    loci[match.group(1)].append(r)

    print(f"\n===== {species} =====")
    print("Gene_ID\tChromosome\tStart\tEnd\tStrand\tCDS_segments\tCDS_bp")

    for gene, records in sorted(loci.items()):

        chrom = ",".join(sorted({r[0] for r in records}))
        strand = ",".join(sorted({r[6] for r in records}))

        start = min(int(r[3]) for r in records)
        end = max(int(r[4]) for r in records)

        length = sum(int(r[4]) - int(r[3]) + 1 for r in records)

        print(gene, chrom, start, end, strand, len(records), length, sep="\t")

    if species != "C_gasperettii":
        print("Raw annotation counts:", dict(Counter(g.split("-")[0] for g in loci)))

# 4. inspect N. naja published toxin gene catalog
wb = load_workbook(
    "N_naja/Suryamohan_2020_supplementary_tables.xlsx",
    read_only=True,
    data_only=True
)

print("\n===== N. NAJA: SUPPLEMENTARY TABLE 6B =====")
print("Scaffold\tStart\tEnd\tStrand\tGene_ID\tHomolog\tProtein_length\tORF_status\tProtein_family")

for row in wb["S6b. VG genes"].iter_rows(min_row=5, values_only=True):

    if not row[4]:
        continue

    family = str(row[24] or "")
    homolog = str(row[5] or "")
    text = (family + " " + homolog).lower()

    if any(x in text for x in ["metalloproteinase", "phospholipase a2", "phospholipase a 2", "svmp", "pla2"]):

        values = [row[i] for i in [0, 1, 2, 3, 4, 5, 6, 13, 24]]

        print("\t".join("" if x is None else str(x) for x in values))

print("\n===== N. NAJA: SUPPLEMENTARY TABLE 6C =====")

for row in wb["S6c. VG gene scaffold info"].iter_rows(values_only=True):

    text = " ".join(str(x or "") for x in row).lower()

    if any(x in text for x in ["metalloproteinase", "phospholipase", "svmp", "pla2"]):

        print("\t".join("" if x is None else str(x) for x in row))

wb.close()

PY
```

Below are species-specific observations.

##### G. ussuriensis
Our curated core set contains 17 SVMPs and five venom PLA2 candidates. The apparent recovery of an additional PLA2 annotation from the GTF is explained by the inclusion of PLA2-1, which we previously reclassified as PLA2I-like-1. Its original identifier was preserved in the GTF, but this gene should not be included in the venom PLA2 counts.

| Gene ID | CDS length (bp) | Current classification |
|---|---:|---|
| PLA2-1 | 465 | PLA2 inhibitor-like protein |
| PLA2-2 | 357 | Core PLA2 candidate |
| PLA2-3 | 423 | Core PLA2 candidate |
| PLA2-4 | 417 | Disrupted PLA2 gene model |
| PLA2-5 | 417 | Core PLA2 candidate |
| PLA2-6 | 417 | Core PLA2 candidate |
| PLA2-7 | 417 | Core PLA2 candidate |

Note that PLA2-5 and PLA2-6 have identical CDS sequences but are located at distinct genomic loci approximately 10 kb apart. We already decided to retain both annotations in the core set. ALso, all five core venom PLA2 candidates and the disrupted PLA2-4 gene model are located on chromosome 16. The PLA2 inhibitor-like gene is located on chromosome 10.

For SVMP, our core set contains 17 candidates located on chromosome 9. The nonredundant set contains 19 candidate gene models, including models that were excluded from the core set during our previous annotation curation.

##### C. adamanteus
The published toxin annotation contains 23 SVMP annotations, all of which are located on Cadam_mi-1. The PLA2 family contains five annotated gene models (PLA2-A through PLA2-E), all of which are located on Cadam_mi-8.

| Gene ID | CDS length (bp) |
|---|---:|
| PLA2-A | 417 |
| PLA2-B | 420 |
| PLA2-C | 417 |
| PLA2-D | 375 |
| PLA2-E | 471 |

These five annotations appear to include both venom-associated and nonvenom PLA2 homologs. We expect three venom PLA2 candidates, but their identities need to be verified before finalizing the toxin gene inventory. The original SVMP identifiers also contain suffixes that may correspond to differences in protein domain architecture. These will be examined during SVMP classification.

##### C. gasperettii
The curated genome annotation contains 13 SVMP gene models located on CM101933.1. An additional SVMP ancestor-like gene model (rna-SVMP-ancestor-1) is located on a separate chromosome and will be considered separately from the venom SVMP cluster.

The curated venom PLA2 cluster contains two gene models on CM101942.1.

| Gene ID | Start | End | Strand | CDS length (bp) |
|---|---:|---:|:---:|---:|
| rna-CG-PLA2-gD | 7891295 | 7892874 | − | 417 |
| rna-CG-PLA2-gC | 7901523 | 7904795 | + | 441 |

Additional PLA2-related gene models occur on other chromosomes. We will not automatically count these as venom PLA2 genes simply because their annotation identifiers contain PLA2. These additional homologs may nevertheless be useful as reference sequences for differentiating venom-associated PLA2s from other PLA2 family members.

##### N. naja
The published toxin gene catalog contains eight SVMP gene models and two group I venom PLA2 gene models. Six of the eight SVMP gene models were reported as having complete ORFs, whereas two were reported as incomplete, and all eight SVMP gene models occur on ScVE01q_1043;HRSCAF=1199.

The two venom PLA2 gene models (Nana39244 and Nana39246) occur on ScVE01q_344;HRSCAF=411. Additional PLA2 family members were also identified in the published catalog. However, these correspond to other PLA2 groups and should not be directly counted as venom group I PLA2 genes.

One issue we need to resolve here is that the genomic sequence identifiers and gene identifiers in the published toxin catalog differ from those in our locally available NCBI genome annotation.

For example, ScVE01q_344;HRSCAF=411 corresponds to the NCBI sequence SOZL01000342.1 based on the FASTA header. However, the coordinates of the two published venom PLA2 gene models do not directly correspond to the gene models identified in the local general genome annotation.

So, we will need to reconcile the published toxin gene coordinates with our local genome assembly and annotation before extracting gene models or conducting comparative microsynteny.

#### Summary of preliminary toxin gene inventories and taxon selection
Based on our preliminary investigation, the initial SVMP and PLA2 inventories can be summarized like so:

| Species | SVMP candidates | Venom PLA2 candidates | Annotation source |
|---|---:|---:|---|
| G. ussuriensis | 17 | 5 | Curated core set |
| G. shedaoensis | 1 intact + 1 partial (chromosomal) | 3 (provisional) | ToxCodAn, Liftoff, and tBLASTn |
| C. adamanteus | 23 | 3 (provisional) | Published toxin annotation |
| C. gasperettii | 13 | 2 | Published curated annotation |
| N. naja | 8 (6 complete + 2 incomplete) | 2 | Published toxin catalog |

These are preliminary annotation inventories and should not yet be interpreted as directly comparable estimates of complete genomic toxin repertoires.

For *G. ussuriensis*, the table only includes the curated core SVMP and venom PLA2 sets. Disrupted models and other homologs will be considered separately during gene model reconciliation. For *G. shedaoensis*, the additional SVMP annotations on unplaced scaffolds are not included in the chromosomal counts above. Given the apparent incomplete resolution of the SVMP cluster of *G. shedaoensis*, I will use different taxon sets for SVMP and PLA2 microsynteny.

| Species | SVMP microsynteny | PLA2 microsynteny |
|---|---|---|
| G. ussuriensis | Include | Include |
| G. shedaoensis | Exclude | Include |
| C. adamanteus | Include | Include |
| C. gasperettii | Include | Include |
| N. naja | Pending coordinate validation | Pending coordinate validation |

This gives us four species for SVMP microsynteny and five species for PLA2 microsynteny. The exclusion of *G. shedaoensis* from SVMP microsynteny is based on the apparent incomplete assembly of the specific genomic region under investigation, rather than the overall quality of the genome assembly. The PLA2 cluster appears to be considerably better resolved, and retaining this species is particularly useful given its phylogenetic proximity to our focal species. For context, four venom PLA2 isoforms have previously been reported for *G. intermedius*. However, counts obtained from venom-derived sequences should not be directly equated with genomic toxin gene copy numbers.

#### Remaining tasks before proceeding to Step 2
At this point, most reference preparation and preliminary annotation inspection have now been completed. However, there are two things we need to resolve before establishing our initial toxin sequence datasets:

1) Verify which of the five published PLA2 annotations of *C. adamanteus* correspond to venom-associated PLA2s and which correspond to other PLA2 family members.

2) Reconcile the published SVMP and PLA2 gene coordinates of *N. naja* with our locally available genome assembly and general genome annotation.

The additional PLA2 homologs identified in *G. shedaoensis* will be examined further during the standardized homology searches and gene classification steps. In particular, we will need to determine the identity of Gshe_Contig255_G00046 and resolve the structure of Gshe_Contig255_G00049. We will also identify conserved flanking genes around the SVMP and PLA2 clusters during subsequent microsynteny analyses.

### Step 2: Extract existing SVMP and PLA2 gene models and sequences

(To be conducted.)

In this step, we will extract existing SVMP and PLA2 gene models and their associated coding and protein sequences from the reference data prepared in Step 1.

These sequences will serve as starting material for standardized homology searches, which will be used to reconcile existing annotations and identify additional toxin homologs.

We will retain the original gene identifiers throughout this process to facilitate tracking of gene models across different annotation sources.

The downstream workflow will broadly consist of:

1) Extraction of existing toxin gene models and associated sequences.

2) Construction of a common toxin and related nontoxin query dataset.

3) Standardized homology searches against the selected reference genomes.

4) Reconciliation and targeted curation of toxin gene models, including identification of redundant, incomplete, and disrupted annotations.

5) Classification of SVMP and PLA2 genes based on their protein sequences, domain architectures, and phylogenetic relationships.

6) Identification of conserved flanking genes and comparative chromosomal microsynteny.

For SVMPs, we will distinguish P-I, P-II, and P-III proteins where their domain architectures can be reliably established. For PLA2s, we will distinguish the major evolutionary groups, venom-associated and nonvenom homologs, and relevant mature-protein residue classifications.

We will primarily focus on describing differences in the genomic organization of these two major toxin gene families. Given that we are comparing individual reference genome assemblies, differences in gene copy number and cluster organization will not automatically be interpreted as fixed species-level differences or evidence of evolutionary gene gain or loss.