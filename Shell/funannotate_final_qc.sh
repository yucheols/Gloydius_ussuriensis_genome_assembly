#!/bin/bash
#SBATCH --job-name=funannotate_final_qc
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# final annotation QC
#
# run after completing these steps:
#
#   1. funannotate train
#   2. funannotate predict
#   3. funannotate update
#   5. InterProScan
#   6. funannotate annotate
#
# This script performs:
#
#   - clean restart of final QC directory
#   - final Funannotate output-file detection
#   - gene/transcript/exon/CDS counts
#   - FASTA sequence statistics
#   - FASTA ID uniqueness checks
#   - longest-protein inspection
#   - genome/GFF3 contig consistency checks
#   - GFF3 coordinate validation
#   - GFF3 ID/Parent validation
#   - valid multipart-ID vs conflicting-ID distinction
#   - functional annotation coverage summary
#   - gffread reconstruction validation
#   - BUSCO protein-mode QC
#
# BUSCO lineage:
#
#   squamata_odb12.2
#
# BUSCO runs completely offline using:
#
#   /home/yshin/mendel-nas1/busco_downloads/
#       lineages/squamata_odb12.2
#
# ============================================================


# ------------------------------------------------------------
# initialize conda
# ------------------------------------------------------------

source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
set -euo pipefail


# ============================================================
# PATHS
# ============================================================

FUN_DIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
ANN_DIR="${FUN_DIR}/annotate_results"
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
OUT="${FUN_DIR}/final_annotation_QC"
COMPLEASM_LIBRARY="/home/yshin/mendel-nas1/busco_downloads/lineages"
COMPLEASM_LINEAGE="sauropsida"
COMPLEASM_ODB="odb12"
COMPLEASM_DATASET="${COMPLEASM_LIBRARY}/${COMPLEASM_LINEAGE}_${COMPLEASM_ODB}"
COMPLEASM_OUT="${OUT}/compleasm"


# ============================================================
# 0. PRE-FLIGHT CHECKS
# ============================================================

echo
echo "============================================================"
echo "0. PRE-FLIGHT CHECKS"
echo "============================================================"
echo

echo "Date:"
date
echo

echo "Host:"
hostname
echo

echo "Funannotate directory:"
echo "${FUN_DIR}"
echo

echo "Annotation directory:"
echo "${ANN_DIR}"
echo

echo "Reference genome:"
echo "${REF}"
echo

echo "compleasm lineage:"
echo "${COMPLEASM_LINEAGE}_${COMPLEASM_ODB}"
echo


# ------------------------------------------------------------
# check input directories/files
# ------------------------------------------------------------

if [[ ! -d "${ANN_DIR}" ]]; then
    echo "ERROR: annotate_results directory does not exist:"
    echo "${ANN_DIR}"
    exit 1
fi


if [[ ! -s "${REF}" ]]; then
    echo "ERROR: reference genome does not exist or is empty:"
    echo "${REF}"
    exit 1
fi


# ------------------------------------------------------------
# check local BUSCO lineage BEFORE deleting/running anything
# ------------------------------------------------------------

if [[ ! -d "${COMPLEASM_DATASET}" ]]; then
    echo "ERROR: local compleasm lineage directory does not exist:"
    echo "${COMPLEASM_DATASET}"
    exit 1
fi


if [[ ! -s "${COMPLEASM_DATASET}/dataset.cfg" ]]; then
    echo "ERROR: compleasm dataset.cfg does not exist:"
    echo "${COMPLEASM_DATASET}/dataset.cfg"
    exit 1
fi


if [[ ! -d "${COMPLEASM_DATASET}/hmms" ]]; then
    echo "ERROR: compleasm HMM directory does not exist:"
    echo "${COMPLEASM_DATASET}/hmms"
    exit 1
fi


echo "Local compleasm lineage found:"
ls -lh "${COMPLEASM_DATASET}/dataset.cfg"
echo


# ------------------------------------------------------------
# check BUSCO environment without contacting network
# ------------------------------------------------------------

conda activate compleasm


if ! command -v compleasm >/dev/null 2>&1; then
    echo "ERROR: compleasm executable not found in compleasm environment."
    exit 1
fi


COMPLEASM_EXEC=$(command -v compleasm)


echo "compleasm executable:"
echo "${COMPLEASM_EXEC}"
echo


# ------------------------------------------------------------
# obtain BUSCO count from local dataset.cfg if available
# ------------------------------------------------------------

COMPLEASM_EXPECTED=$(awk -F '=' '
    $1 ~ /^number_of_BUSCOs$/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${COMPLEASM_DATASET}/dataset.cfg" || true)


if [[ -z "${COMPLEASM_EXPECTED}" ]]; then
    COMPLEASM_EXPECTED="NA"
fi


echo "Expected Sauropsida BUSCOs:"
echo "${COMPLEASM_EXPECTED}"
echo


# ------------------------------------------------------------
# activate funannotate environment
# ------------------------------------------------------------

conda activate funannotate


if ! command -v python >/dev/null 2>&1; then
    echo "ERROR: python not found in funannotate environment."
    exit 1
fi


echo "Funannotate environment:"
echo "${CONDA_DEFAULT_ENV:-UNKNOWN}"
echo


# ============================================================
# 1. CLEAN START
# ============================================================

echo "============================================================"
echo "1. CLEANING PREVIOUS QC OUTPUT"
echo "============================================================"
echo


# safety check before rm -rf

EXPECTED_OUT="${FUN_DIR}/final_annotation_QC"


if [[ "${OUT}" != "${EXPECTED_OUT}" ]]; then
    echo "ERROR: unexpected QC output path."
    echo "Refusing to remove:"
    echo "${OUT}"
    exit 1
fi


if [[ -d "${OUT}" ]]; then
    echo "Removing previous QC directory:"
    echo "${OUT}"
    rm -rf "${OUT}"
fi


mkdir -p \
    "${OUT}" \
    "${OUT}/logs" \
    "${OUT}/compleasm" \
    "${OUT}/gffread"


echo "Fresh QC directory created:"
echo "${OUT}"
echo


# ============================================================
# 2. DETECT FINAL FUNANNOTATE OUTPUT FILES
# ============================================================

echo "============================================================"
echo "2. DETECTING FINAL FUNANNOTATE FILES"
echo "============================================================"
echo


ls -lh "${ANN_DIR}" \
    | tee "${OUT}/annotate_results_file_list.txt"


echo


# ------------------------------------------------------------
# helper function
# ------------------------------------------------------------

find_first_file() {

    local pattern
    local result

    for pattern in "$@"; do

        result=$(find "${ANN_DIR}" \
            -maxdepth 1 \
            -type f \
            -name "${pattern}" \
            | sort \
            | head -n 1 || true)

        if [[ -n "${result}" ]]; then
            echo "${result}"
            return 0
        fi

    done

    return 1
}


# ------------------------------------------------------------
# detect GFF
# ------------------------------------------------------------

GFF=$(find_first_file \
    "*.gff3" \
    "*.gff" \
    || true)


# ------------------------------------------------------------
# detect protein FASTA
# ------------------------------------------------------------

PROT=$(find_first_file \
    "*.proteins.fa" \
    "*.proteins.fasta" \
    "*protein*.fa" \
    "*protein*.fasta" \
    || true)


# ------------------------------------------------------------
# detect transcript FASTA
# ------------------------------------------------------------

TRANS=$(find_first_file \
    "*.transcripts.fa" \
    "*.transcripts.fasta" \
    "*transcript*.fa" \
    "*transcript*.fasta" \
    "*mrna*.fa" \
    "*mrna*.fasta" \
    || true)


# prevent CDS file from accidentally being selected as transcript file

if [[ -n "${TRANS}" ]] \
    && [[ "$(basename "${TRANS}")" == *cds* ]]; then
    TRANS=""
fi


# ------------------------------------------------------------
# detect CDS FASTA
# ------------------------------------------------------------

CDS=$(find_first_file \
    "*.cds-transcripts.fa" \
    "*.cds-transcripts.fasta" \
    "*cds*.fa" \
    "*cds*.fasta" \
    || true)


echo "Detected files:"
echo
echo "Reference    = ${REF}"
echo "GFF3         = ${GFF:-NOT FOUND}"
echo "Proteins     = ${PROT:-NOT FOUND}"
echo "Transcripts  = ${TRANS:-NOT FOUND}"
echo "CDS          = ${CDS:-NOT FOUND}"
echo


# ------------------------------------------------------------
# require critical files
# ------------------------------------------------------------

if [[ -z "${GFF}" ]] || [[ ! -s "${GFF}" ]]; then
    echo "ERROR: final GFF3 could not be detected."
    exit 1
fi


if [[ -z "${PROT}" ]] || [[ ! -s "${PROT}" ]]; then
    echo "ERROR: final protein FASTA could not be detected."
    exit 1
fi


# ------------------------------------------------------------
# save file paths
# ------------------------------------------------------------

{
    echo -e "file_type\tpath"
    echo -e "reference\t${REF}"
    echo -e "gff3\t${GFF}"
    echo -e "proteins\t${PROT}"
    echo -e "transcripts\t${TRANS:-NA}"
    echo -e "cds\t${CDS:-NA}"
} > "${OUT}/input_files.tsv"


cat "${OUT}/input_files.tsv"

echo


# ============================================================
# 3. BASIC FASTA COUNTS
# ============================================================

echo "============================================================"
echo "3. FASTA RECORD COUNTS"
echo "============================================================"
echo


REF_COUNT=$(grep -c '^>' "${REF}" || true)

PROT_COUNT=$(grep -c '^>' "${PROT}" || true)


if [[ -n "${TRANS}" ]] && [[ -s "${TRANS}" ]]; then
    TRANS_COUNT=$(grep -c '^>' "${TRANS}" || true)
else
    TRANS_COUNT="NA"
fi


if [[ -n "${CDS}" ]] && [[ -s "${CDS}" ]]; then
    CDS_COUNT=$(grep -c '^>' "${CDS}" || true)
else
    CDS_COUNT="NA"
fi


{
    echo -e "data_type\tcount"
    echo -e "genome_sequences\t${REF_COUNT}"
    echo -e "proteins\t${PROT_COUNT}"
    echo -e "transcripts\t${TRANS_COUNT}"
    echo -e "CDS\t${CDS_COUNT}"
} > "${OUT}/fasta_record_counts.tsv"


cat "${OUT}/fasta_record_counts.tsv"

echo


# ============================================================
# 4. GFF3 FEATURE COUNTS
# ============================================================

echo "============================================================"
echo "4. GFF3 FEATURE COUNTS"
echo "============================================================"
echo


{
    echo -e "feature\tcount"

    awk -F '\t' '
        $0 !~ /^#/ && NF >= 9 {
            count[$3]++
        }

        END {
            for (feature in count) {
                print feature "\t" count[feature]
            }
        }
    ' "${GFF}" \
        | sort -k1,1

} > "${OUT}/gff_feature_counts.tsv"


cat "${OUT}/gff_feature_counts.tsv"

echo


# ------------------------------------------------------------
# major feature counts
# ------------------------------------------------------------

GENE_COUNT=$(awk -F '\t' '
    $0 !~ /^#/ && $3=="gene" {
        n++
    }
    END {
        print n+0
    }
' "${GFF}")


MRNA_COUNT=$(awk -F '\t' '
    $0 !~ /^#/ && ($3=="mRNA" || $3=="transcript") {
        n++
    }
    END {
        print n+0
    }
' "${GFF}")


EXON_COUNT=$(awk -F '\t' '
    $0 !~ /^#/ && $3=="exon" {
        n++
    }
    END {
        print n+0
    }
' "${GFF}")


GFF_CDS_COUNT=$(awk -F '\t' '
    $0 !~ /^#/ && $3=="CDS" {
        n++
    }
    END {
        print n+0
    }
' "${GFF}")


{
    echo -e "metric\tvalue"
    echo -e "genes\t${GENE_COUNT}"
    echo -e "mRNAs\t${MRNA_COUNT}"
    echo -e "exons\t${EXON_COUNT}"
    echo -e "CDS_features\t${GFF_CDS_COUNT}"
} > "${OUT}/annotation_counts.tsv"


cat "${OUT}/annotation_counts.tsv"

echo


# ============================================================
# 5. FASTA LENGTH / INTEGRITY STATISTICS
# ============================================================

echo "============================================================"
echo "5. FASTA LENGTH AND SEQUENCE INTEGRITY"
echo "============================================================"
echo


python - \
    "${REF}" \
    "${PROT}" \
    "${TRANS:-NA}" \
    "${CDS:-NA}" \
    "${OUT}/fasta_length_statistics.tsv" \
    "${OUT}/fasta_integrity_summary.tsv" \
    "${OUT}/longest_proteins.tsv" <<'PYTHON'

import os
import sys
import statistics
from collections import Counter


ref_file = sys.argv[1]
protein_file = sys.argv[2]
transcript_file = sys.argv[3]
cds_file = sys.argv[4]

stats_file = sys.argv[5]
integrity_file = sys.argv[6]
longest_file = sys.argv[7]


def read_fasta(path):

    records = []

    header = None
    seq_parts = []

    with open(path) as handle:

        for line in handle:

            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    records.append(
                        (header, "".join(seq_parts))
                    )

                header = line[1:]
                seq_parts = []

            else:
                seq_parts.append(line)

    if header is not None:
        records.append(
            (header, "".join(seq_parts))
        )

    return records


def n50(lengths):

    if not lengths:
        return 0

    target = sum(lengths) / 2
    cumulative = 0

    for value in sorted(lengths, reverse=True):

        cumulative += value

        if cumulative >= target:
            return value

    return 0


datasets = [
    ("genome", ref_file),
    ("proteins", protein_file),
    ("transcripts", transcript_file),
    ("CDS", cds_file),
]


# ------------------------------------------------------------
# length statistics
# ------------------------------------------------------------

with open(stats_file, "w") as out:

    out.write(
        "dataset\t"
        "sequences\t"
        "total_length\t"
        "min\t"
        "median\t"
        "mean\t"
        "N50\t"
        "max\n"
    )

    for name, path in datasets:

        if path == "NA" or not os.path.isfile(path):

            out.write(
                f"{name}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n"
            )

            continue

        records = read_fasta(path)

        lengths = [
            len(seq)
            for _, seq in records
        ]

        if not lengths:

            out.write(
                f"{name}\t0\t0\tNA\tNA\tNA\tNA\tNA\n"
            )

            continue

        out.write(
            f"{name}\t"
            f"{len(lengths)}\t"
            f"{sum(lengths)}\t"
            f"{min(lengths)}\t"
            f"{statistics.median(lengths):.2f}\t"
            f"{statistics.mean(lengths):.2f}\t"
            f"{n50(lengths)}\t"
            f"{max(lengths)}\n"
        )


# ------------------------------------------------------------
# integrity statistics
# ------------------------------------------------------------

with open(integrity_file, "w") as out:

    out.write(
        "dataset\t"
        "sequences\t"
        "duplicate_primary_IDs\t"
        "empty_sequences\t"
        "invalid_character_sequences\n"
    )

    allowed = {
        "genome": set("ACGTNRYKMSWBDHVacgtnrykmswbdhv"),
        "proteins": set("ABCDEFGHIKLMNPQRSTVWXYZJUO*abcdefghiklmnpqrstvwxyzjuo"),
        "transcripts": set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv"),
        "CDS": set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv"),
    }

    for name, path in datasets:

        if path == "NA" or not os.path.isfile(path):

            out.write(
                f"{name}\tNA\tNA\tNA\tNA\n"
            )

            continue

        records = read_fasta(path)

        primary_ids = [
            header.split()[0]
            for header, _ in records
        ]

        id_counts = Counter(primary_ids)

        duplicate_ids = sum(
            1
            for count in id_counts.values()
            if count > 1
        )

        empty_sequences = sum(
            1
            for _, seq in records
            if len(seq) == 0
        )

        invalid_sequences = sum(
            1
            for _, seq in records
            if any(char not in allowed[name] for char in seq)
        )

        out.write(
            f"{name}\t"
            f"{len(records)}\t"
            f"{duplicate_ids}\t"
            f"{empty_sequences}\t"
            f"{invalid_sequences}\n"
        )


# ------------------------------------------------------------
# longest proteins
# ------------------------------------------------------------

proteins = read_fasta(protein_file)

proteins.sort(
    key=lambda x: len(x[1]),
    reverse=True
)


with open(longest_file, "w") as out:

    out.write(
        "rank\tlength_aa\tprotein_header\n"
    )

    for rank, (header, seq) in enumerate(
        proteins[:20],
        start=1
    ):

        out.write(
            f"{rank}\t"
            f"{len(seq)}\t"
            f"{header}\n"
        )

PYTHON


cat "${OUT}/fasta_length_statistics.tsv"

echo

cat "${OUT}/fasta_integrity_summary.tsv"

echo

echo "20 longest predicted proteins:"
cat "${OUT}/longest_proteins.tsv"

echo


# ============================================================
# 6. GENOME / GFF3 CONTIG CONSISTENCY
# ============================================================

echo "============================================================"
echo "6. GENOME / GFF3 CONTIG CONSISTENCY"
echo "============================================================"
echo


grep '^>' "${REF}" \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | sort -u \
    > "${OUT}/reference_contigs.txt"


awk -F '\t' '
    $0 !~ /^#/ && NF >= 9 {
        print $1
    }
' "${GFF}" \
    | sort -u \
    > "${OUT}/gff_contigs.txt"


comm -23 \
    "${OUT}/gff_contigs.txt" \
    "${OUT}/reference_contigs.txt" \
    > "${OUT}/gff_contigs_missing_from_reference.txt"


comm -13 \
    "${OUT}/gff_contigs.txt" \
    "${OUT}/reference_contigs.txt" \
    > "${OUT}/reference_contigs_without_annotations.txt"


MISSING_CONTIGS=$(wc -l < "${OUT}/gff_contigs_missing_from_reference.txt")

UNANNOTATED_CONTIGS=$(wc -l < "${OUT}/reference_contigs_without_annotations.txt")


echo "GFF contigs absent from genome:    ${MISSING_CONTIGS}"
echo "Genome contigs without annotation: ${UNANNOTATED_CONTIGS}"
echo


if [[ "${MISSING_CONTIGS}" -gt 0 ]]; then

    echo "WARNING:"
    echo "GFF3 contains sequence IDs absent from the reference:"
    echo

    cat "${OUT}/gff_contigs_missing_from_reference.txt"

    echo
fi


# ============================================================
# 7. GFF3 COORDINATE VALIDATION
# ============================================================

echo "============================================================"
echo "7. GFF3 COORDINATE VALIDATION"
echo "============================================================"
echo


python - \
    "${REF}" \
    "${GFF}" \
    "${OUT}/gff_coordinate_errors.tsv" \
    "${OUT}/gff_coordinate_summary.tsv" <<'PYTHON'

import sys


fasta = sys.argv[1]
gff = sys.argv[2]
errors_file = sys.argv[3]
summary_file = sys.argv[4]


reference_lengths = {}

sequence_id = None
sequence_length = 0


with open(fasta) as handle:

    for line in handle:

        if line.startswith(">"):

            if sequence_id is not None:
                reference_lengths[sequence_id] = sequence_length

            sequence_id = line[1:].split()[0]
            sequence_length = 0

        else:

            sequence_length += len(line.strip())


if sequence_id is not None:
    reference_lengths[sequence_id] = sequence_length


features_checked = 0
errors = 0
malformed_lines = 0
unknown_seqid = 0
invalid_start_end = 0
outside_reference = 0


with open(errors_file, "w") as out:

    out.write(
        "line\t"
        "seqid\t"
        "start\t"
        "end\t"
        "error\n"
    )

    with open(gff) as handle:

        for line_number, line in enumerate(handle, 1):

            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:

                errors += 1
                malformed_lines += 1

                out.write(
                    f"{line_number}\t"
                    "NA\tNA\tNA\t"
                    "malformed_gff_line\n"
                )

                continue

            features_checked += 1

            seqid = fields[0]

            try:

                start = int(fields[3])
                end = int(fields[4])

            except ValueError:

                errors += 1
                invalid_start_end += 1

                out.write(
                    f"{line_number}\t"
                    f"{seqid}\t"
                    f"{fields[3]}\t"
                    f"{fields[4]}\t"
                    "non_integer_coordinate\n"
                )

                continue

            if seqid not in reference_lengths:

                errors += 1
                unknown_seqid += 1

                out.write(
                    f"{line_number}\t"
                    f"{seqid}\t"
                    f"{start}\t"
                    f"{end}\t"
                    "seqid_not_in_reference\n"
                )

                continue

            if start < 1 or end < start:

                errors += 1
                invalid_start_end += 1

                out.write(
                    f"{line_number}\t"
                    f"{seqid}\t"
                    f"{start}\t"
                    f"{end}\t"
                    "invalid_start_end\n"
                )

                continue

            if end > reference_lengths[seqid]:

                errors += 1
                outside_reference += 1

                out.write(
                    f"{line_number}\t"
                    f"{seqid}\t"
                    f"{start}\t"
                    f"{end}\t"
                    f"end_exceeds_contig_length_"
                    f"{reference_lengths[seqid]}\n"
                )


with open(summary_file, "w") as out:

    out.write("metric\tvalue\n")
    out.write(f"features_checked\t{features_checked}\n")
    out.write(f"errors\t{errors}\n")
    out.write(f"malformed_lines\t{malformed_lines}\n")
    out.write(f"unknown_seqid\t{unknown_seqid}\n")
    out.write(f"invalid_start_end\t{invalid_start_end}\n")
    out.write(f"outside_reference\t{outside_reference}\n")

PYTHON


cat "${OUT}/gff_coordinate_summary.tsv"

echo


# ============================================================
# 8. GFF3 ID / PARENT VALIDATION
#
# IMPORTANT:
#
# GFF3 permits repeated IDs for parts of a discontinuous
# feature. Therefore repeated IDs are NOT automatically errors.
#
# We distinguish:
#
#   repeated_ID_names
#
# from:
#
#   conflicting_repeated_ID_names
#
# A repeated ID is considered conflicting if its occurrences
# disagree in sequence ID, feature type, strand, or Parent.
# ============================================================

echo "============================================================"
echo "8. GFF3 ID / PARENT VALIDATION"
echo "============================================================"
echo


python - \
    "${GFF}" \
    "${OUT}/gff_id_parent_details.tsv" \
    "${OUT}/gff_parent_summary.tsv" <<'PYTHON'

import sys
from collections import defaultdict


gff = sys.argv[1]
details_file = sys.argv[2]
summary_file = sys.argv[3]


id_records = defaultdict(list)

all_ids = set()

parent_records = []


def parse_attributes(text):

    attrs = {}

    for item in text.split(";"):

        item = item.strip()

        if not item:
            continue

        if "=" in item:

            key, value = item.split("=", 1)
            attrs[key] = value

    return attrs


with open(gff) as handle:

    for line_number, line in enumerate(handle, 1):

        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        seqid = fields[0]
        feature_type = fields[2]
        start = fields[3]
        end = fields[4]
        strand = fields[6]

        attrs = parse_attributes(fields[8])

        feature_id = attrs.get("ID")

        parent_text = attrs.get("Parent", "")

        parents = tuple(
            sorted(
                p
                for p in parent_text.split(",")
                if p
            )
        )

        if feature_id:

            all_ids.add(feature_id)

            id_records[feature_id].append(
                {
                    "line": line_number,
                    "seqid": seqid,
                    "type": feature_type,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "parents": parents,
                }
            )

        for parent in parents:

            parent_records.append(
                {
                    "parent": parent,
                    "line": line_number,
                    "seqid": seqid,
                    "type": feature_type,
                }
            )


repeated_ids = {
    feature_id: records
    for feature_id, records in id_records.items()
    if len(records) > 1
}


conflicting_ids = {}


for feature_id, records in repeated_ids.items():

    signatures = {
        (
            record["seqid"],
            record["type"],
            record["strand"],
            record["parents"],
        )
        for record in records
    }

    if len(signatures) > 1:

        conflicting_ids[feature_id] = records


missing_parents = [
    record
    for record in parent_records
    if record["parent"] not in all_ids
]


with open(details_file, "w") as out:

    out.write(
        "status\t"
        "ID_or_Parent\t"
        "line\t"
        "seqid\t"
        "feature\t"
        "start\t"
        "end\n"
    )

    for feature_id, records in repeated_ids.items():

        status = (
            "conflicting_repeated_ID"
            if feature_id in conflicting_ids
            else "valid_repeated_ID"
        )

        for record in records:

            out.write(
                f"{status}\t"
                f"{feature_id}\t"
                f"{record['line']}\t"
                f"{record['seqid']}\t"
                f"{record['type']}\t"
                f"{record['start']}\t"
                f"{record['end']}\n"
            )

    for record in missing_parents:

        out.write(
            "missing_Parent\t"
            f"{record['parent']}\t"
            f"{record['line']}\t"
            f"{record['seqid']}\t"
            f"{record['type']}\t"
            "NA\tNA\n"
        )


with open(summary_file, "w") as out:

    out.write("metric\tvalue\n")

    out.write(
        f"unique_IDs\t"
        f"{len(all_ids)}\n"
    )

    out.write(
        f"repeated_ID_names\t"
        f"{len(repeated_ids)}\n"
    )

    out.write(
        f"conflicting_repeated_ID_names\t"
        f"{len(conflicting_ids)}\n"
    )

    out.write(
        f"Parent_links\t"
        f"{len(parent_records)}\n"
    )

    out.write(
        f"missing_Parent_links\t"
        f"{len(missing_parents)}\n"
    )

PYTHON


cat "${OUT}/gff_parent_summary.tsv"

echo


# ============================================================
# 9. FUNCTIONAL ANNOTATION COVERAGE
# ============================================================

echo "============================================================"
echo "9. FUNCTIONAL ANNOTATION COVERAGE"
echo "============================================================"
echo


python - \
    "${GFF}" \
    "${OUT}/functional_annotation_coverage.tsv" <<'PYTHON'

import sys


gff = sys.argv[1]
output_file = sys.argv[2]


total_transcripts = 0


counts = {
    "product": 0,
    "Name": 0,
    "Dbxref": 0,
    "Ontology_term": 0,
    "GO": 0,
    "InterPro": 0,
    "Pfam": 0,
    "MEROPS": 0,
    "CAZy_or_dbCAN": 0,
    "EggNog": 0,
}


with open(gff) as handle:

    for line in handle:

        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        if fields[2] not in ("mRNA", "transcript"):
            continue

        total_transcripts += 1

        attributes = fields[8]


        if "product=" in attributes:
            counts["product"] += 1

        if "Name=" in attributes:
            counts["Name"] += 1

        if "Dbxref=" in attributes:
            counts["Dbxref"] += 1

        if "Ontology_term=" in attributes:
            counts["Ontology_term"] += 1

        if "GO:" in attributes:
            counts["GO"] += 1

        if (
            "InterPro" in attributes
            or "IPR" in attributes
        ):
            counts["InterPro"] += 1

        if (
            "Pfam" in attributes
            or "PFAM" in attributes
        ):
            counts["Pfam"] += 1

        if "MEROPS" in attributes:
            counts["MEROPS"] += 1

        if (
            "CAZy" in attributes
            or "CAZY" in attributes
            or "dbCAN" in attributes
        ):
            counts["CAZy_or_dbCAN"] += 1

        if (
            "EggNog" in attributes
            or "eggNOG" in attributes
            or "EGGNOG" in attributes
        ):
            counts["EggNog"] += 1


with open(output_file, "w") as out:

    out.write(
        "annotation_type\t"
        "transcripts_with_annotation\t"
        "total_transcripts\t"
        "percent\n"
    )

    for annotation_type, count in counts.items():

        percent = (
            100.0 * count / total_transcripts
            if total_transcripts
            else 0.0
        )

        out.write(
            f"{annotation_type}\t"
            f"{count}\t"
            f"{total_transcripts}\t"
            f"{percent:.2f}\n"
        )

PYTHON


cat "${OUT}/functional_annotation_coverage.tsv"

echo


# ============================================================
# 10. GFFREAD VALIDATION
# ============================================================

echo "============================================================"
echo "10. GFFREAD VALIDATION"
echo "============================================================"
echo


GFFREAD_STATUS="NOT_RUN"

GFFREAD_TRANS="NA"
GFFREAD_CDS="NA"
GFFREAD_PROT="NA"


if command -v gffread >/dev/null 2>&1; then

    echo "gffread executable:"
    command -v gffread
    echo


    gffread \
        "${GFF}" \
        -g "${REF}" \
        -w "${OUT}/gffread/reconstructed_transcripts.fa" \
        -x "${OUT}/gffread/reconstructed_CDS.fa" \
        -y "${OUT}/gffread/reconstructed_proteins.fa" \
        -E \
        2> "${OUT}/logs/gffread_validation.log"


    GFFREAD_TRANS=$(grep -c '^>' \
        "${OUT}/gffread/reconstructed_transcripts.fa" \
        || true)


    GFFREAD_CDS=$(grep -c '^>' \
        "${OUT}/gffread/reconstructed_CDS.fa" \
        || true)


    GFFREAD_PROT=$(grep -c '^>' \
        "${OUT}/gffread/reconstructed_proteins.fa" \
        || true)


    GFFREAD_STATUS="PASS"


    if [[ "${GFFREAD_PROT}" != "${PROT_COUNT}" ]]; then
        GFFREAD_STATUS="WARNING"
    fi


    if [[ "${CDS_COUNT}" != "NA" ]] \
        && [[ "${GFFREAD_CDS}" != "${CDS_COUNT}" ]]; then
        GFFREAD_STATUS="WARNING"
    fi


    if [[ "${TRANS_COUNT}" != "NA" ]] \
        && [[ "${GFFREAD_TRANS}" != "${TRANS_COUNT}" ]]; then
        GFFREAD_STATUS="WARNING"
    fi


    {
        echo -e "metric\tcount"
        echo -e "reconstructed_transcripts\t${GFFREAD_TRANS}"
        echo -e "reconstructed_CDS\t${GFFREAD_CDS}"
        echo -e "reconstructed_proteins\t${GFFREAD_PROT}"
        echo -e "status\t${GFFREAD_STATUS}"
    } > "${OUT}/gffread_reconstruction_counts.tsv"


else

    echo "WARNING:"
    echo "gffread was not found in the funannotate environment."
    echo

    {
        echo -e "metric\tcount"
        echo -e "reconstructed_transcripts\tNA"
        echo -e "reconstructed_CDS\tNA"
        echo -e "reconstructed_proteins\tNA"
        echo -e "status\tNOT_RUN"
    } > "${OUT}/gffread_reconstruction_counts.tsv"

fi


cat "${OUT}/gffread_reconstruction_counts.tsv"

echo


# ============================================================
# 11. EXTRACT VALUES NEEDED FOR STRUCTURAL QC
# ============================================================

COORD_ERRORS=$(awk -F '\t' '
    $1=="errors" {
        print $2
    }
' "${OUT}/gff_coordinate_summary.tsv")


CONFLICTING_IDS=$(awk -F '\t' '
    $1=="conflicting_repeated_ID_names" {
        print $2
    }
' "${OUT}/gff_parent_summary.tsv")


REPEATED_IDS=$(awk -F '\t' '
    $1=="repeated_ID_names" {
        print $2
    }
' "${OUT}/gff_parent_summary.tsv")


MISSING_PARENTS=$(awk -F '\t' '
    $1=="missing_Parent_links" {
        print $2
    }
' "${OUT}/gff_parent_summary.tsv")


# ------------------------------------------------------------
# FASTA duplicate primary IDs
# ------------------------------------------------------------

REF_DUP_IDS=$(awk -F '\t' '
    $1=="genome" {
        print $3
    }
' "${OUT}/fasta_integrity_summary.tsv")


PROT_DUP_IDS=$(awk -F '\t' '
    $1=="proteins" {
        print $3
    }
' "${OUT}/fasta_integrity_summary.tsv")


# ============================================================
# 12. STRUCTURAL QC STATUS
# ============================================================

STRUCTURAL_STATUS="PASS"


if [[ "${MISSING_CONTIGS}" -ne 0 ]] \
    || [[ "${COORD_ERRORS}" -ne 0 ]] \
    || [[ "${CONFLICTING_IDS}" -ne 0 ]] \
    || [[ "${MISSING_PARENTS}" -ne 0 ]] \
    || [[ "${REF_DUP_IDS}" -ne 0 ]] \
    || [[ "${PROT_DUP_IDS}" -ne 0 ]]; then

    STRUCTURAL_STATUS="WARNING"

fi


echo "============================================================"
echo "12. STRUCTURAL QC STATUS"
echo "============================================================"
echo

echo "Status: ${STRUCTURAL_STATUS}"
echo

echo "GFF contigs absent from genome:      ${MISSING_CONTIGS}"
echo "Coordinate errors:                   ${COORD_ERRORS}"
echo "Repeated GFF IDs:                    ${REPEATED_IDS}"
echo "Conflicting repeated GFF IDs:        ${CONFLICTING_IDS}"
echo "Missing Parent links:                ${MISSING_PARENTS}"
echo "Duplicate genome FASTA IDs:          ${REF_DUP_IDS}"
echo "Duplicate protein FASTA IDs:         ${PROT_DUP_IDS}"
echo "gffread status:                      ${GFFREAD_STATUS}"
echo


# ============================================================
# 13. BUSCO
#
# IMPORTANT:
#
# No --list-datasets call.
#
# No network lineage name.
#
# BUSCO receives the full path to the local lineage and is
# explicitly run with --offline.
# ============================================================

echo "============================================================"
echo "13. COMPLEASM ON FINAL PROTEIN SET"
echo "============================================================"
echo


conda activate compleasm


echo "Current conda environment:"
echo "${CONDA_DEFAULT_ENV:-UNKNOWN}"
echo

echo "compleasm executable:"
command -v compleasm
echo

echo "compleasm input:"
echo "${PROT}"
echo

echo "compleasm lineage:"
echo "${COMPLEASM_LINEAGE}"
echo

echo "OrthoDB:"
echo "${COMPLEASM_ODB}"
echo

echo "compleasm mode:"
echo "proteins"
echo

echo "Expected BUSCOs:"
echo "${COMPLEASM_EXPECTED}"
echo


# ------------------------------------------------------------
# BUSCO run
# ------------------------------------------------------------

rm -rf "${COMPLEASM_OUT}"

compleasm protein \
    -p "${PROT}" \
    -o "${COMPLEASM_OUT}" \
    -t "${SLURM_CPUS_PER_TASK}" \
    -l "${COMPLEASM_LINEAGE}" \
    -L "${COMPLEASM_LIBRARY}" \
    --odb "${COMPLEASM_ODB}"


echo
echo "compleasm command completed."
echo


# ============================================================
# 14. LOCATE BUSCO SUMMARY
# ============================================================

COMPLEASM_SUMMARY="${COMPLEASM_OUT}/summary.txt"


if [[ ! -s "${COMPLEASM_SUMMARY}" ]]; then

    echo "ERROR:"
    echo "compleasm completed but summary.txt was not found."
    exit 1

fi


cp \
    "${COMPLEASM_SUMMARY}" \
    "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt"


echo "compleasm summary:"
echo

cat "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt"

echo


# ============================================================
# 15. EXTRACT BUSCO METRICS
# ============================================================

COMPLEASM_SINGLE=$(awk -F ',' '
    /^S:/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt" || true)


COMPLEASM_DUPLICATED=$(awk -F ',' '
    /^D:/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt" || true)


COMPLEASM_FRAGMENTED=$(awk -F ',' '
    /^F:/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt" || true)


COMPLEASM_INTERSPACED=$(awk -F ',' '
    /^I:/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt" || true)


COMPLEASM_MISSING=$(awk -F ',' '
    /^M:/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt" || true)


COMPLEASM_TOTAL=$(awk -F ':' '
    /^N:/ {
        gsub(/[[:space:]]/, "", $2)
        print $2
        exit
    }
' "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt" || true)


if [[ -n "${COMPLEASM_SINGLE}" ]] && [[ -n "${COMPLEASM_DUPLICATED}" ]]; then
    COMPLEASM_COMPLETE=$((COMPLEASM_SINGLE + COMPLEASM_DUPLICATED))
else
    COMPLEASM_COMPLETE="NA"
fi


{
    echo -e "metric\tcount"
    echo -e "expected_BUSCOs\t${COMPLEASM_EXPECTED}"
    echo -e "complete_BUSCOs\t${COMPLEASM_COMPLETE}"
    echo -e "single_copy_BUSCOs\t${COMPLEASM_SINGLE:-NA}"
    echo -e "duplicated_BUSCOs\t${COMPLEASM_DUPLICATED:-NA}"
    echo -e "fragmented_BUSCOs\t${COMPLEASM_FRAGMENTED:-NA}"
    echo -e "interspaced_BUSCOs\t${COMPLEASM_INTERSPACED:-NA}"
    echo -e "missing_BUSCOs\t${COMPLEASM_MISSING:-NA}"
    echo -e "total_BUSCOs\t${COMPLEASM_TOTAL:-NA}"
} > "${OUT}/compleasm_counts.tsv"


cat "${OUT}/compleasm_counts.tsv"

echo


# ============================================================
# 16. FINAL QC SUMMARY
# ============================================================

echo "============================================================"
echo "16. GENERATING FINAL QC SUMMARY"
echo "============================================================"
echo


{
    echo "============================================================"
    echo "Gloydius ussuriensis"
    echo "FINAL FUNANNOTATE ANNOTATION QC"
    echo "============================================================"
    echo


    echo "INPUT FILES"
    echo "------------------------------------------------------------"
    echo

    echo "Genome:"
    echo "${REF}"
    echo

    echo "GFF3:"
    echo "${GFF}"
    echo

    echo "Proteins:"
    echo "${PROT}"
    echo

    echo "Transcripts:"
    echo "${TRANS:-NA}"
    echo

    echo "CDS:"
    echo "${CDS:-NA}"
    echo


    echo "ANNOTATION COUNTS"
    echo "------------------------------------------------------------"
    echo

    echo "Genome sequences:       ${REF_COUNT}"
    echo "Genes:                  ${GENE_COUNT}"
    echo "mRNAs/transcripts:      ${MRNA_COUNT}"
    echo "Exons:                  ${EXON_COUNT}"
    echo "CDS features:           ${GFF_CDS_COUNT}"
    echo "Protein sequences:      ${PROT_COUNT}"
    echo "Transcript sequences:   ${TRANS_COUNT}"
    echo "CDS sequences:          ${CDS_COUNT}"
    echo


    echo "STRUCTURAL VALIDATION"
    echo "------------------------------------------------------------"
    echo

    echo "Overall structural status:           ${STRUCTURAL_STATUS}"
    echo
    echo "GFF contigs absent from genome:      ${MISSING_CONTIGS}"
    echo "Genome contigs without annotation:   ${UNANNOTATED_CONTIGS}"
    echo "Coordinate errors:                   ${COORD_ERRORS}"
    echo "Repeated GFF IDs:                    ${REPEATED_IDS}"
    echo "Conflicting repeated GFF IDs:        ${CONFLICTING_IDS}"
    echo "Missing Parent links:                ${MISSING_PARENTS}"
    echo "Duplicate reference FASTA IDs:       ${REF_DUP_IDS}"
    echo "Duplicate protein FASTA IDs:         ${PROT_DUP_IDS}"
    echo


    echo "NOTE ON REPEATED GFF3 IDs"
    echo "------------------------------------------------------------"
    echo

    echo "Repeated IDs are not automatically errors."
    echo "GFF3 permits the same ID on multiple rows belonging to"
    echo "one discontinuous feature."
    echo
    echo "The relevant error metric is:"
    echo
    echo "Conflicting repeated GFF IDs = ${CONFLICTING_IDS}"
    echo


    echo "GFFREAD VALIDATION"
    echo "------------------------------------------------------------"
    echo

    echo "Status:                       ${GFFREAD_STATUS}"
    echo "Original transcripts:         ${TRANS_COUNT}"
    echo "Reconstructed transcripts:    ${GFFREAD_TRANS}"
    echo "Original CDS sequences:       ${CDS_COUNT}"
    echo "Reconstructed CDS:            ${GFFREAD_CDS}"
    echo "Original proteins:            ${PROT_COUNT}"
    echo "Reconstructed proteins:       ${GFFREAD_PROT}"
    echo


    echo "FUNCTIONAL ANNOTATION COVERAGE"
    echo "------------------------------------------------------------"
    echo

    cat "${OUT}/functional_annotation_coverage.tsv"

    echo


    echo "FASTA LENGTH STATISTICS"
    echo "------------------------------------------------------------"
    echo

    cat "${OUT}/fasta_length_statistics.tsv"

    echo


    echo "COMPLEASM"
    echo "------------------------------------------------------------"
    echo

    echo "Dataset: sauropsida_odb12"
    echo "Dataset path: ${COMPLEASM_DATASET}"
    echo "Expected BUSCOs: ${COMPLEASM_EXPECTED}"
    echo "Mode: proteins"
    echo

    cat "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt"

    echo


    echo "QC OUTPUT DIRECTORY"
    echo "------------------------------------------------------------"
    echo

    echo "${OUT}"
    echo


} > "${OUT}/FINAL_QC_SUMMARY.txt"


# ============================================================
# 17. PRINT FINAL SUMMARY
# ============================================================

cat "${OUT}/FINAL_QC_SUMMARY.txt"


echo
echo "============================================================"
echo "FINAL ANNOTATION QC COMPLETE"
echo "============================================================"
echo

echo "Main summary:"
echo "${OUT}/FINAL_QC_SUMMARY.txt"
echo

echo "Annotation counts:"
echo "${OUT}/annotation_counts.tsv"
echo

echo "GFF feature counts:"
echo "${OUT}/gff_feature_counts.tsv"
echo

echo "FASTA statistics:"
echo "${OUT}/fasta_length_statistics.tsv"
echo

echo "FASTA integrity:"
echo "${OUT}/fasta_integrity_summary.tsv"
echo

echo "Longest proteins:"
echo "${OUT}/longest_proteins.tsv"
echo

echo "Functional annotation coverage:"
echo "${OUT}/functional_annotation_coverage.tsv"
echo

echo "GFF coordinate QC:"
echo "${OUT}/gff_coordinate_summary.tsv"
echo

echo "GFF ID / Parent QC:"
echo "${OUT}/gff_parent_summary.tsv"
echo

echo "Detailed repeated-ID information:"
echo "${OUT}/gff_id_parent_details.tsv"
echo

echo "gffread reconstruction:"
echo "${OUT}/gffread_reconstruction_counts.tsv"
echo

echo "compleasm summary:"
echo "${OUT}/compleasm_final_proteins_sauropsida_odb12_summary.txt"
echo

echo "compleasm counts:"
echo "${OUT}/compleasm_counts.tsv"
echo

echo "Finished:"
date
echo