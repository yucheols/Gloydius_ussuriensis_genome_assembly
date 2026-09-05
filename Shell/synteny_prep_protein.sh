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

set -euo pipefail


# paths
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"
cd "$ROOT"

SUMMARY="${ROOT}/annotation_protein_inventory.tsv"


# check required programs
for cmd in gffread gzip grep awk find readlink; do

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $cmd" >&2
        exit 1
    fi

done


# initialize summary
printf \
"species\tgenome\tannotation\tannotation_format\tprotein\tCDS_records\tprotein_sequences\taction\tstatus\n" \
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
    f=$(find "$dir" -maxdepth 2 -type f \
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

    # --------------------------------------------------------
    # preferred standardized filenames
    # --------------------------------------------------------

    if [[ -s "${dir}/${sp}.annotation.gff3" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gff3"
        return
    fi

    if [[ -s "${dir}/${sp}.annotation.gtf" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gtf"
        return
    fi

    if [[ -s "${dir}/${sp}.annotation.gff" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gff"
        return
    fi

    if [[ -s "${dir}/${sp}.annotation.gff3.gz" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gff3.gz"
        return
    fi

    if [[ -s "${dir}/${sp}.annotation.gtf.gz" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gtf.gz"
        return
    fi

    if [[ -s "${dir}/${sp}.annotation.gff.gz" ]]; then
        printf '%s\n' "${dir}/${sp}.annotation.gff.gz"
        return
    fi


    # --------------------------------------------------------
    # fallback search: GFF3
    # --------------------------------------------------------

    f=$(find "$dir" -maxdepth 2 -type f \
        \( -iname '*.gff3' -o -iname '*.gff3.gz' \) \
        ! -iname '*original*' \
        ! -iname '*old*' \
        ! -iname '*repeat_annotation*' \
        -print -quit)

    if [[ -n "$f" ]]; then
        printf '%s\n' "$f"
        return
    fi


    # --------------------------------------------------------
    # fallback search: GTF
    # --------------------------------------------------------

    f=$(find "$dir" -maxdepth 2 -type f \
        \( -iname '*.gtf' -o -iname '*.gtf.gz' \) \
        ! -iname '*original*' \
        ! -iname '*old*' \
        ! -iname '*repeat_annotation*' \
        -print -quit)

    if [[ -n "$f" ]]; then
        printf '%s\n' "$f"
        return
    fi


    # --------------------------------------------------------
    # fallback search: GFF
    # --------------------------------------------------------

    f=$(find "$dir" -maxdepth 2 -type f \
        \( -iname '*.gff' -o -iname '*.gff.gz' \) \
        ! -iname '*original*' \
        ! -iname '*old*' \
        ! -iname '*repeat_annotation*' \
        -print -quit)

    printf '%s\n' "$f"
}


# ============================================================
# helper: determine annotation format
# ============================================================

annotation_format() {

    local f="$1"

    case "$f" in

        *.gff3|*.gff3.gz)
            printf 'GFF3\n'
            ;;

        *.gtf|*.gtf.gz)
            printf 'GTF\n'
            ;;

        *.gff|*.gff.gz)
            printf 'GFF\n'
            ;;

        *)
            printf 'UNKNOWN\n'
            ;;

    esac
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
    f=$(find "$dir" -maxdepth 2 -type f \
        \( -iname '*.faa' -o -iname '*.faa.gz' \) \
        -print -quit)

    if [[ -n "$f" ]]; then
        printf '%s\n' "$f"
        return
    fi

    # protein-named FASTAs
    f=$(find "$dir" -maxdepth 2 -type f \
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
# helper: count FASTA sequences
# ============================================================

count_fasta_sequences() {

    local f="$1"

    awk '
        /^>/ {
            n++
        }
        END {
            print n+0
        }
    ' "$f"
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


    # ========================================================
    # genome
    # ========================================================

    genome=$(find_genome "$dir" "$sp")

    if [[ -z "$genome" || ! -s "$genome" ]]; then

        echo "Genome: MISSING"

        printf \
        "%s\tMISSING\tMISSING\tNA\tMISSING\tNA\tNA\tnone\tFAILED_NO_GENOME\n" \
            "$sp" \
            >> "$SUMMARY"

        continue

    fi

    echo "Genome:"
    echo "  ${genome}"


    # ========================================================
    # annotation
    # ========================================================

    annotation=$(find_annotation "$dir" "$sp")

    if [[ -z "$annotation" || ! -s "$annotation" ]]; then

        echo "Annotation: MISSING"

        protein=$(find_protein "$dir" "$sp")

        if [[ -n "$protein" && -s "$protein" ]]; then

            echo "Protein FASTA:"
            echo "  ${protein}"

            if [[ "$protein" == *.gz ]]; then

                echo "Decompressing protein FASTA..."

                gzip -dc "$protein" \
                    > "${dir}/${sp}.protein.faa"

            elif [[ "$protein" != "${dir}/${sp}.protein.faa" ]]; then

                ln -sfn \
                    "$(readlink -f "$protein")" \
                    "${dir}/${sp}.protein.faa"

            fi

            nprot=$(count_fasta_sequences "${dir}/${sp}.protein.faa")

            printf \
            "%s\t%s\tMISSING\tNA\t%s\tNA\t%s\tprotein_exists\tWARNING_NO_ANNOTATION\n" \
                "$sp" \
                "$genome" \
                "${dir}/${sp}.protein.faa" \
                "$nprot" \
                >> "$SUMMARY"

        else

            printf \
            "%s\t%s\tMISSING\tNA\tMISSING\tNA\tNA\tnone\tFAILED_NO_ANNOTATION\n" \
                "$sp" \
                "$genome" \
                >> "$SUMMARY"

        fi

        continue

    fi


    format=$(annotation_format "$annotation")

    echo "Annotation:"
    echo "  ${annotation}"

    echo "Annotation format:"
    echo "  ${format}"


    # ========================================================
    # prepare compressed annotation
    # ========================================================

    annotation_use="$annotation"
    annotation_tmp=""

    if [[ "$annotation" == *.gz ]]; then

        case "$annotation" in

            *.gff3.gz)
                annotation_tmp="${dir}/.${sp}.annotation.tmp.gff3"
                ;;

            *.gtf.gz)
                annotation_tmp="${dir}/.${sp}.annotation.tmp.gtf"
                ;;

            *.gff.gz)
                annotation_tmp="${dir}/.${sp}.annotation.tmp.gff"
                ;;

            *)
                annotation_tmp="${dir}/.${sp}.annotation.tmp"
                ;;

        esac

        echo "Decompressing annotation temporarily..."

        gzip -dc "$annotation" > "$annotation_tmp"

        annotation_use="$annotation_tmp"

    fi


    # ========================================================
    # count CDS records
    # ========================================================

    ncds=$(awk -F '\t' \
        '$0 !~ /^#/ && $3=="CDS"{n++} END{print n+0}' \
        "$annotation_use")

    echo "CDS records: ${ncds}"


    # ========================================================
    # existing protein FASTA
    # ========================================================

    protein=$(find_protein "$dir" "$sp")

    if [[ -n "$protein" && -s "$protein" ]]; then

        echo "Protein FASTA already exists:"
        echo "  ${protein}"


        # ----------------------------------------------------
        # standardize protein filename
        # ----------------------------------------------------

        if [[ "$protein" == *.gz ]]; then

            echo "Decompressing protein FASTA..."

            gzip -dc "$protein" \
                > "${dir}/${sp}.protein.faa"

        elif [[ "$protein" != "${dir}/${sp}.protein.faa" ]]; then

            ln -sfn \
                "$(readlink -f "$protein")" \
                "${dir}/${sp}.protein.faa"

        fi


        # ----------------------------------------------------
        # count proteins
        # ----------------------------------------------------

        nprot=$(count_fasta_sequences "${dir}/${sp}.protein.faa")

        echo "Protein sequences: ${nprot}"


        # ----------------------------------------------------
        # cleanup temporary annotation
        # ----------------------------------------------------

        if [[ -n "$annotation_tmp" ]]; then
            rm -f "$annotation_tmp"
        fi


        # ----------------------------------------------------
        # summary
        # ----------------------------------------------------

        printf \
        "%s\t%s\t%s\t%s\t%s\t%s\t%s\tprotein_exists\tOK\n" \
            "$sp" \
            "$genome" \
            "$annotation" \
            "$format" \
            "${dir}/${sp}.protein.faa" \
            "$ncds" \
            "$nprot" \
            >> "$SUMMARY"

        continue

    fi


    # ========================================================
    # protein missing
    # ========================================================

    echo "Protein FASTA: MISSING"


    # --------------------------------------------------------
    # annotation must contain CDS
    # --------------------------------------------------------

    if [[ "$ncds" -eq 0 ]]; then

        echo "ERROR: Annotation contains no CDS features."
        echo "Cannot derive proteins."


        if [[ -n "$annotation_tmp" ]]; then
            rm -f "$annotation_tmp"
        fi


        printf \
        "%s\t%s\t%s\t%s\tMISSING\t0\tNA\tcannot_generate\tFAILED_NO_CDS\n" \
            "$sp" \
            "$genome" \
            "$annotation" \
            "$format" \
            >> "$SUMMARY"

        continue

    fi


    # ========================================================
    # generate protein FASTA with gffread
    # ========================================================

    outfile="${dir}/${sp}.protein.faa"
    tmp="${outfile}.tmp"

    rm -f "$tmp"

    echo "Generating protein FASTA with gffread..."
    echo "Annotation used:"
    echo "  ${annotation_use}"


    if gffread \
        "$annotation_use" \
        -g "$genome" \
        -y "$tmp"; then


        # ----------------------------------------------------
        # validate generated FASTA
        # ----------------------------------------------------

        if [[ -s "$tmp" ]] && grep -m1 -q '^>' "$tmp"; then

            mv "$tmp" "$outfile"

            nprot=$(count_fasta_sequences "$outfile")

            echo "Protein generation: OK"
            echo "Protein sequences: ${nprot}"


            if [[ -n "$annotation_tmp" ]]; then
                rm -f "$annotation_tmp"
            fi


            printf \
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\tgenerated_from_annotation\tOK\n" \
                "$sp" \
                "$genome" \
                "$annotation" \
                "$format" \
                "$outfile" \
                "$ncds" \
                "$nprot" \
                >> "$SUMMARY"


        else

            echo "ERROR: gffread produced no valid protein FASTA."

            rm -f "$tmp"


            if [[ -n "$annotation_tmp" ]]; then
                rm -f "$annotation_tmp"
            fi


            printf \
            "%s\t%s\t%s\t%s\tMISSING\t%s\t0\tgffread_failed\tFAILED\n" \
                "$sp" \
                "$genome" \
                "$annotation" \
                "$format" \
                "$ncds" \
                >> "$SUMMARY"

        fi


    else

        echo "ERROR: gffread failed for ${sp}"

        rm -f "$tmp"


        if [[ -n "$annotation_tmp" ]]; then
            rm -f "$annotation_tmp"
        fi


        printf \
        "%s\t%s\t%s\t%s\tMISSING\t%s\tNA\tgffread_failed\tFAILED\n" \
            "$sp" \
            "$genome" \
            "$annotation" \
            "$format" \
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

# print when done
echo
echo "Summary written to:"
echo "  ${SUMMARY}"
echo