#!/bin/bash
#SBATCH --job-name=synteny_prep_proteins
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate synteny_qc
set -uo pipefail

# root directory
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"
cd "$ROOT"

# check programs
for cmd in gffread gzip grep awk find readlink; do

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $cmd" >&2
        exit 1
    fi

done

# output summary
SUMMARY="${ROOT}/annotation_protein_inventory.tsv"

printf \
"species\tgenome\tannotation\tprotein\tCDS_records\tprotein_sequences\taction\tstatus\n" \
> "$SUMMARY"


# ============================================================
# helper: choose genome
# ============================================================

find_genome() {

    local dir="$1"
    local sp="$2"
    local f=""

    # preferred standardized filename
    if [[ -s "${dir}/${sp}.genome.fa" ]]; then
        printf '%s\n' "${dir}/${sp}.genome.fa"
        return
    fi

    # otherwise search
    f=$(find "$dir" -type f \
        \( -iname '*genomic.fna' \
           -o -iname '*.genome.fa' \
           -o -iname '*.genome.fasta' \
           -o -iname '*genome*.fa' \
           -o -iname '*genome*.fasta' \) \
        -print -quit)

    printf '%s\n' "$f"
}


# ============================================================
# helper: choose annotation
# ============================================================

find_annotation() {

    local dir="$1"
    local sp="$2"
    local f=""

    # preferred standardized filename
    if [[ -s "${dir}/${sp}.annotation.gff3" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gff3"
        return
    fi

    # prefer GFF3
    f=$(find "$dir" -type f \
        \( -iname '*.gff3' -o -iname '*.gff3.gz' \) \
        -print -quit)

    if [[ -n "$f" ]]; then
        printf '%s\n' "$f"
        return
    fi

    # then GFF
    f=$(find "$dir" -type f \
        \( -iname '*.gff' -o -iname '*.gff.gz' \) \
        -print -quit)

    if [[ -n "$f" ]]; then
        printf '%s\n' "$f"
        return
    fi

    # finally GTF
    f=$(find "$dir" -type f \
        \( -iname '*.gtf' -o -iname '*.gtf.gz' \) \
        -print -quit)

    printf '%s\n' "$f"
}


# ============================================================
# helper: choose existing protein FASTA
# ============================================================

find_protein() {

    local dir="$1"
    local sp="$2"
    local f=""

    # preferred standardized filename
    if [[ -s "${dir}/${sp}.protein.faa" ]]; then
        printf '%s\n' "${dir}/${sp}.protein.faa"
        return
    fi

    # *.faa
    f=$(find "$dir" -type f \
        \( -iname '*.faa' -o -iname '*.faa.gz' \) \
        -print -quit)

    if [[ -n "$f" ]]; then
        printf '%s\n' "$f"
        return
    fi

    # protein-named FASTAs
    f=$(find "$dir" -type f \
        \( -iname '*protein*.fa' \
           -o -iname '*protein*.fasta' \
           -o -iname '*protein*.fa.gz' \
           -o -iname '*protein*.fasta.gz' \
           -o -iname '*proteins*.fa' \
           -o -iname '*proteins*.fasta' \
           -o -iname '*pep*.fa' \
           -o -iname '*pep*.fasta' \) \
        -print -quit)

    printf '%s\n' "$f"
}


# ============================================================
# process each species directory
# ============================================================

for dir in "${ROOT}"/*/; do

    [[ -d "$dir" ]] || continue

    sp=$(basename "$dir")

    echo
    echo "============================================================"
    echo "Species: ${sp}"
    echo "============================================================"


    # --------------------------------------------------------
    # find genome
    # --------------------------------------------------------

    genome=$(find_genome "$dir" "$sp")

    if [[ -z "$genome" || ! -s "$genome" ]]; then

        echo "Genome: MISSING"

        printf \
        "%s\tMISSING\tMISSING\tMISSING\tNA\tNA\tnone\tFAILED_NO_GENOME\n" \
            "$sp" >> "$SUMMARY"

        continue

    fi

    echo "Genome:"
    echo "  ${genome}"


    # --------------------------------------------------------
    # find annotation
    # --------------------------------------------------------

    annotation=$(find_annotation "$dir" "$sp")

    if [[ -z "$annotation" || ! -s "$annotation" ]]; then

        echo "Annotation: MISSING"

        protein=$(find_protein "$dir" "$sp")

        if [[ -n "$protein" && -s "$protein" ]]; then

            echo "Protein FASTA:"
            echo "  ${protein}"

            if [[ "$protein" == *.gz ]]; then

                gzip -dc "$protein" > "${dir}/${sp}.protein.faa"

            elif [[ "$protein" != "${dir}/${sp}.protein.faa" ]]; then

                ln -sfn \
                    "$(readlink -f "$protein")" \
                    "${dir}/${sp}.protein.faa"

            fi

            nprot=$(grep -c '^>' "${dir}/${sp}.protein.faa")

            printf \
            "%s\t%s\tMISSING\t%s\tNA\t%s\tprotein_exists\tWARNING_NO_ANNOTATION\n" \
                "$sp" \
                "$genome" \
                "${dir}/${sp}.protein.faa" \
                "$nprot" \
                >> "$SUMMARY"

        else

            printf \
            "%s\t%s\tMISSING\tMISSING\tNA\tNA\tnone\tFAILED_NO_ANNOTATION\n" \
                "$sp" \
                "$genome" \
                >> "$SUMMARY"

        fi

        continue

    fi

    echo "Annotation:"
    echo "  ${annotation}"


    # --------------------------------------------------------
    # prepare annotation for gffread if compressed
    # --------------------------------------------------------

    annotation_use="$annotation"
    annotation_tmp=""

    if [[ "$annotation" == *.gz ]]; then

        annotation_tmp="${dir}/.${sp}.annotation.tmp.gff"

        echo "Decompressing annotation temporarily..."

        gzip -dc "$annotation" > "$annotation_tmp"

        annotation_use="$annotation_tmp"

    fi


    # --------------------------------------------------------
    # standardized annotation symlink
    #
    # only create when annotation itself is uncompressed
    # --------------------------------------------------------

    if [[ "$annotation_use" != "${dir}/${sp}.annotation.gff3" ]]; then

        if [[ -z "$annotation_tmp" ]]; then

            ln -sfn \
                "$(readlink -f "$annotation_use")" \
                "${dir}/${sp}.annotation.gff3"

        else

            cp "$annotation_use" \
                "${dir}/${sp}.annotation.gff3"

            annotation_use="${dir}/${sp}.annotation.gff3"

            rm -f "$annotation_tmp"
            annotation_tmp=""

        fi

    else

        annotation_use="${dir}/${sp}.annotation.gff3"

    fi


    # --------------------------------------------------------
    # count CDS records
    # --------------------------------------------------------

    ncds=$(awk -F '\t' \
        '$0 !~ /^#/ && $3=="CDS"{n++} END{print n+0}' \
        "$annotation_use")

    echo "CDS records: ${ncds}"


    # --------------------------------------------------------
    # find existing protein FASTA
    # --------------------------------------------------------

    protein=$(find_protein "$dir" "$sp")

    if [[ -n "$protein" && -s "$protein" ]]; then

        echo "Protein FASTA already exists:"
        echo "  ${protein}"

        # standardize filename
        if [[ "$protein" == *.gz ]]; then

            echo "Decompressing protein FASTA..."

            gzip -dc "$protein" \
                > "${dir}/${sp}.protein.faa"

        elif [[ "$protein" != "${dir}/${sp}.protein.faa" ]]; then

            ln -sfn \
                "$(readlink -f "$protein")" \
                "${dir}/${sp}.protein.faa"

        fi

        nprot=$(grep -c '^>' "${dir}/${sp}.protein.faa")

        echo "Protein sequences: ${nprot}"

        printf \
        "%s\t%s\t%s\t%s\t%s\t%s\tprotein_exists\tOK\n" \
            "$sp" \
            "$genome" \
            "${dir}/${sp}.annotation.gff3" \
            "${dir}/${sp}.protein.faa" \
            "$ncds" \
            "$nprot" \
            >> "$SUMMARY"

        continue

    fi


    # --------------------------------------------------------
    # protein missing -> generate from GFF + genome
    # --------------------------------------------------------

    echo "Protein FASTA: MISSING"


    if [[ "$ncds" -eq 0 ]]; then

        echo "ERROR: Annotation contains no CDS features."
        echo "Cannot derive proteins."

        printf \
        "%s\t%s\t%s\tMISSING\t0\tNA\tcannot_generate\tFAILED_NO_CDS\n" \
            "$sp" \
            "$genome" \
            "${dir}/${sp}.annotation.gff3" \
            >> "$SUMMARY"

        continue

    fi


    outfile="${dir}/${sp}.protein.faa"
    tmp="${outfile}.tmp"

    rm -f "$tmp"

    echo "Generating protein FASTA with gffread..."

    if gffread \
        "${dir}/${sp}.annotation.gff3" \
        -g "$genome" \
        -y "$tmp"; then

        if [[ -s "$tmp" ]] && grep -m1 -q '^>' "$tmp"; then

            mv "$tmp" "$outfile"

            nprot=$(grep -c '^>' "$outfile")

            echo "Protein generation: OK"
            echo "Protein sequences: ${nprot}"

            printf \
            "%s\t%s\t%s\t%s\t%s\t%s\tgenerated_from_annotation\tOK\n" \
                "$sp" \
                "$genome" \
                "${dir}/${sp}.annotation.gff3" \
                "$outfile" \
                "$ncds" \
                "$nprot" \
                >> "$SUMMARY"

        else

            echo "ERROR: gffread produced no valid protein FASTA."

            rm -f "$tmp"

            printf \
            "%s\t%s\t%s\tMISSING\t%s\t0\tgffread_failed\tFAILED\n" \
                "$sp" \
                "$genome" \
                "${dir}/${sp}.annotation.gff3" \
                "$ncds" \
                >> "$SUMMARY"

        fi

    else

        echo "ERROR: gffread failed for ${sp}"

        rm -f "$tmp"

        printf \
        "%s\t%s\t%s\tMISSING\t%s\tNA\tgffread_failed\tFAILED\n" \
            "$sp" \
            "$genome" \
            "${dir}/${sp}.annotation.gff3" \
            "$ncds" \
            >> "$SUMMARY"

    fi

done


# ============================================================
# print summary
# ============================================================

echo
echo "============================================================"
echo "ANNOTATION / PROTEIN SUMMARY"
echo "============================================================"
echo

if command -v column >/dev/null 2>&1; then

    column -t -s $'\t' "$SUMMARY"

else

    cat "$SUMMARY"

fi

echo
echo "Summary written to:"
echo "  ${SUMMARY}"