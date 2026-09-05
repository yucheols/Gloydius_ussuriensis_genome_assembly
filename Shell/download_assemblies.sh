#!/bin/bash
#SBATCH --job-name=download_assemblies
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=compute
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================
# download snake assemblies and available annotations/proteins
# for synteny analyses
#
# Based on synteny_samples_list(4).csv.
#
# Rules:
#   1. Download each genome according to genome_url.
#   2. If annotation_url is present, download the annotation
#      and standardize it as Species.annotation.gff3.
#   3. If protein_url is present, download the protein FASTA
#      and standardize it as Species.protein.faa.
#
# Important download behavior:
#   - Figshare ndownloader URLs are automatically converted to
#     the public Figshare API download endpoint to avoid the
#     AWS WAF challenge seen from the cluster.
#   - Direct curl downloads are restartable. Partial files are
#     retained as *.part and resumed after connection failures.
#   - Annotation/protein failures generate warnings and do not
#     stop downloads for the remaining species.
#   - Genome failures are fatal.
# ============================================================

# ============================================================
# activate conda environment
# ============================================================

source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate ncbi_datasets

set -euo pipefail


# ============================================================
# paths
# ============================================================

WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny"
OUTDIR="${WORKDIR}/assemblies_synteny"

mkdir -p "$OUTDIR"
cd "$OUTDIR"


# ============================================================
# check required programs
# ============================================================

for cmd in datasets unzip curl gzip jq find readlink tar awk grep; do

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $cmd" >&2
        exit 1
    fi

done


# ============================================================
# manifest from synteny_samples_list(4).csv
#
# fields:
#   species
#   family
#   subfamily
#   assembly_name
#   NCBI_RefSeq
#   genome_accession
#   source
#   genome_url
#   annotation_url
#   protein_url
#
# Empty values are represented by empty fields between |.
# ============================================================

readarray -t MANIFEST <<'EOF_MANIFEST'
Crotalus_adamanteus|Viperidae|Crotalinae|Cadamanteus_3dDNAHiC_1.2||GCA_039797435.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_039797435.1/||
Crotalus_viridis|Viperidae|Crotalinae|UTA_CroVir_3.0||GCA_003400415.2|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003400415.2/|http://herpmuseum.cib.ac.cn/downloads/genomes/Cvir.gene.annotation.tar.gz|
Bothrops_insularis|Viperidae|Crotalinae|ASM5582466v1||GCA_055824665.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_055824665.1/|https://figshare.com/ndownloader/files/52417235|
Gloydius_shedaoensis|Viperidae|Crotalinae|Gshe1||GWHBWDU00000000|NGDC|https://ngdc.cncb.ac.cn/gwh/Assembly/34318/show|http://herpmuseum.cib.ac.cn/downloads/genomes/Gshe.gene.annotation.tar.gz|
Vipera_berus|Viperidae|Viperinae|rVipBer3.hap1.1|GCF_964194415.1|GCA_964194415.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_964194415.1/||
Cerastes_gasperettii|Viperidae|Viperinae|ASM4652402v1||GCA_046524025.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_046524025.1/|https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102647/annotation/curated_annotation.gff|https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/102001_103000/102647/annotation/coding_proteins.fa
Naja_naja|Elapidae|Elapinae|Nana_v5||GCA_009733165.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009733165.1/|http://herpmuseum.cib.ac.cn/downloads/genomes/Nnaj.gene.annotation.tar.gz|
Elaphe_schrenckii|Colubridae|Colubrinae|ASM5023117v1||GCA_050231175.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_050231175.1/|https://figshare.com/ndownloader/files/52988567|
Candoia_aspera|Boidae|Candoiinae|rCanAsp1.hap2|GCF_035149785.1|GCA_035149785.1|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_035149785.1/||
Xenopeltis_unicolor|Xenopeltidae||Xuni1||GWHBWDW00000000|NGDC|https://ngdc.cncb.ac.cn/gwh/Assembly/34320/show|http://herpmuseum.cib.ac.cn/downloads/genomes/Xuni.gene.annotation.tar.gz|
Argyrophis_diardii|Typhlopidae||Adia1||GWHBWDT00000000|NGDC|https://ngdc.cncb.ac.cn/gwh/Assembly/34317/show|http://herpmuseum.cib.ac.cn/downloads/genomes/Adia.gene.annotation.tar.gz|
EOF_MANIFEST


# ============================================================
# save manifest used by this run
# ============================================================

MANIFEST_OUT="${OUTDIR}/assembly_manifest.tsv"

printf "species\tfamily\tsubfamily\tassembly_name\tNCBI_RefSeq\tgenome_accession\tsource\tgenome_url\tannotation_url\tprotein_url\n" \
    > "$MANIFEST_OUT"

for entry in "${MANIFEST[@]}"; do

    IFS='|' read -r \
        sp \
        family \
        subfamily \
        assembly_name \
        refseq \
        genome_acc \
        source \
        genome_url \
        annotation_url \
        protein_url <<< "$entry"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sp" \
        "$family" \
        "$subfamily" \
        "$assembly_name" \
        "$refseq" \
        "$genome_acc" \
        "$source" \
        "$genome_url" \
        "$annotation_url" \
        "$protein_url" \
        >> "$MANIFEST_OUT"

done


# ============================================================
# resolve special download URLs
#
# Figshare's regular ndownloader endpoint currently returns an
# AWS WAF challenge from the Mendel cluster. Convert it to the
# public Figshare API file-download endpoint instead.
# ============================================================

resolve_download_url() {

    local url="$1"

    if [[ "$url" =~ ^https?://(www\.)?figshare\.com/ndownloader/files/([0-9]+) ]]; then

        local file_id="${BASH_REMATCH[2]}"
        printf 'https://api.figshare.com/v2/file/download/%s\n' "$file_id"

    elif [[ "$url" =~ ^https?://ndownloader\.figshare\.com/files/([0-9]+) ]]; then

        local file_id="${BASH_REMATCH[1]}"
        printf 'https://api.figshare.com/v2/file/download/%s\n' "$file_id"

    else

        printf '%s\n' "$url"

    fi
}


# ============================================================
# robust direct URL downloader
#
# The partial file is retained after failures so the next retry
# or SLURM resubmission can continue from the existing byte.
#
# If the URL associated with a .part file changes, the old part
# is discarded to prevent resuming from the wrong resource.
# ============================================================

fetch_url() {

    local original_url="$1"
    local outfile="$2"

    local url
    local partial="${outfile}.part"
    local partial_url="${partial}.url"
    local attempt
    local rc=0

    url=$(resolve_download_url "$original_url")

    if [[ "$url" != "$original_url" ]]; then
        echo "Converted download URL:"
        echo "  original: ${original_url}"
        echo "  resolved: ${url}"
    fi

    # A completed staging file from an interrupted older script
    # cannot be distinguished reliably from an incomplete file.
    # Remove it; resumable state is kept only in *.part.
    rm -f "$outfile"

    # Do not resume a partial file created for a different URL.
    if [[ -s "$partial" ]]; then

        if [[ ! -s "$partial_url" ]] || \
           [[ "$(cat "$partial_url")" != "$url" ]]; then

            echo "Discarding partial file because its source URL changed:"
            echo "  ${partial}"

            rm -f "$partial" "$partial_url"

        fi

    fi

    printf '%s\n' "$url" > "$partial_url"

    for ((attempt=1; attempt<=10; attempt++)); do

        echo "Download attempt ${attempt}/10:"
        echo "  ${url}"

        if [[ -s "$partial" ]]; then
            echo "Resuming existing partial file:"
            echo "  ${partial}"
            echo "  current size: $(du -h "$partial" | cut -f1)"
        fi

        if curl \
            -L \
            --fail \
            --connect-timeout 60 \
            --continue-at - \
            --output "$partial" \
            "$url"; then

            if [[ ! -s "$partial" ]]; then
                echo "ERROR: Download completed but file is empty:" >&2
                echo "  ${url}" >&2
                return 1
            fi

            mv "$partial" "$outfile"
            rm -f "$partial_url"

            echo "Download complete:"
            echo "  ${outfile}"

            return 0

        else

            rc=$?

            echo "WARNING: curl failed with exit code ${rc} on attempt ${attempt}/10." >&2

            # curl exit code 33 means the remote server rejected
            # the requested byte range. Restart from byte zero.
            if [[ "$rc" -eq 33 ]]; then

                echo "Server rejected resume request; restarting this file from byte 0." >&2
                rm -f "$partial"

            fi

            if [[ "$attempt" -lt 10 ]]; then
                echo "Retrying in 10 seconds..."
                sleep 10
            fi

        fi

    done

    echo "ERROR: Download failed after 10 attempts:" >&2
    echo "  ${url}" >&2
    echo "Partial file retained, if present:" >&2
    echo "  ${partial}" >&2

    return 1
}


# ============================================================
# copy/decompress a plain or gzip-compressed file
# ============================================================

copy_maybe_gzip() {

    local infile="$1"
    local outfile="$2"
    local tmp="${outfile}.tmp"

    rm -f "$tmp"

    if gzip -t "$infile" >/dev/null 2>&1; then

        gzip -dc "$infile" > "$tmp"

    else

        cp "$infile" "$tmp"

    fi

    if [[ ! -s "$tmp" ]]; then
        echo "ERROR: Failed to create ${outfile}" >&2
        rm -f "$tmp"
        return 1
    fi

    mv "$tmp" "$outfile"
}


# ============================================================
# validation helpers
# ============================================================

validate_fasta() {

    local fasta="$1"

    if [[ ! -s "$fasta" ]] || ! grep -m1 -q '^>' "$fasta"; then
        echo "ERROR: File does not look like FASTA:" >&2
        echo "  ${fasta}" >&2
        return 1
    fi
}


validate_gff3() {

    local gff="$1"

    if [[ ! -s "$gff" ]]; then
        echo "ERROR: Empty annotation file:" >&2
        echo "  ${gff}" >&2
        return 1
    fi

    if grep -m1 -q '^##gff-version 3' "$gff"; then
        return 0
    fi

    if awk -F '\t' '
        BEGIN { ok=0 }
        /^#/ { next }
        NF == 9 && ($9 ~ /(^|;)ID=/ || $9 ~ /(^|;)Parent=/ || $9 ~ /;Parent=/) {
            ok=1
            exit
        }
        END { exit(ok ? 0 : 1) }
    ' "$gff"; then
        return 0
    fi

    echo "ERROR: Annotation does not look like GFF3:" >&2
    echo "  ${gff}" >&2

    return 1
}


# ============================================================
# download genome according to genome_url
# ============================================================

download_genome() {

    local sp="$1"
    local genome_acc="$2"
    local genome_url="$3"

    local spdir="${OUTDIR}/${sp}"
    local outfile="${spdir}/${sp}.genome.fa"
    local marker="${spdir}/.genome_source_url"

    mkdir -p "$spdir"

    # Skip only when the current standardized genome came from
    # the exact genome_url in the current manifest.
    if [[ -s "$outfile" && -s "$marker" ]] && \
       [[ "$(cat "$marker")" == "$genome_url" ]]; then

        echo "Genome already downloaded from current genome_url; skipping."
        return 0

    fi

    echo "Downloading genome according to genome_url:"
    echo "  ${genome_url}"

    # --------------------------------------------------------
    # NCBI Datasets record URL
    # --------------------------------------------------------

    if [[ "$genome_url" == *"ncbi.nlm.nih.gov/datasets/genome/"* ]]; then

        local ncbi_acc=""
        local zipfile=""
        local datadir=""
        local genome=""

        ncbi_acc=$(printf '%s\n' "$genome_url" \
            | grep -oE 'GC[AF]_[0-9]+\.[0-9]+' \
            | head -n 1 || true)

        if [[ -z "$ncbi_acc" ]]; then
            echo "ERROR: Could not extract GCA/GCF accession from genome_url:" >&2
            echo "  ${genome_url}" >&2
            return 1
        fi

        echo "Resolved NCBI accession: ${ncbi_acc}"

        zipfile="${spdir}/${ncbi_acc}.genome.zip"

        # Reuse an already downloaded NCBI genome package when possible.
        if [[ ! -s "$zipfile" ]]; then

            rm -f "${zipfile}.tmp"

            datasets download genome accession "$ncbi_acc" \
                --include genome,seq-report \
                --filename "${zipfile}.tmp" \
                --no-progressbar

            if [[ ! -s "${zipfile}.tmp" ]]; then
                echo "ERROR: NCBI genome download failed for ${sp}" >&2
                return 1
            fi

            mv "${zipfile}.tmp" "$zipfile"

        else

            echo "Existing NCBI genome ZIP found; reusing:"
            echo "  ${zipfile}"

        fi

        unzip -oq "$zipfile" -d "$spdir"

        datadir="${spdir}/ncbi_dataset/data/${ncbi_acc}"

        if [[ ! -d "$datadir" ]]; then
            echo "ERROR: Expected NCBI directory not found:" >&2
            echo "  ${datadir}" >&2
            return 1
        fi

        genome=$(find "$datadir" \
            -maxdepth 1 \
            -type f \
            \( -name '*_genomic.fna' -o -name 'genomic.fna' \) \
            -print -quit)

        if [[ -z "${genome:-}" ]]; then
            echo "ERROR: Genome FASTA not found for ${sp}" >&2
            return 1
        fi

        ln -sfn \
            "$(readlink -f "$genome")" \
            "$outfile"

    # --------------------------------------------------------
    # NGDC/GWH record URL
    #
    # genome_url is a human-facing record, so resolve the DNA
    # file through the GWH API using genome_accession.
    # --------------------------------------------------------

    elif [[ "$genome_url" == *"ngdc.cncb.ac.cn/gwh/Assembly/"* ]]; then

        local api_json="${spdir}/${genome_acc}.ngdc.json"
        local dna_url=""
        local raw="${spdir}/.${sp}.genome.download"

        if [[ -z "$genome_acc" ]]; then
            echo "ERROR: NGDC genome requires genome_accession for ${sp}" >&2
            return 1
        fi

        echo "Resolving NGDC DNA URL using accession:"
        echo "  ${genome_acc}"

        if ! fetch_url \
            "https://ngdc.cncb.ac.cn/gwh/api/public/assembly/${genome_acc}" \
            "${api_json}.tmp"; then

            return 1

        fi

        mv "${api_json}.tmp" "$api_json"

        dna_url=$(jq -r '.ftpPathDna // empty' "$api_json")

        if [[ -z "$dna_url" ]]; then
            echo "ERROR: NGDC API did not provide ftpPathDna for ${sp}" >&2
            return 1
        fi

        echo "Resolved NGDC DNA URL:"
        echo "  ${dna_url}"

        if ! fetch_url "$dna_url" "$raw"; then
            return 1
        fi

        if ! copy_maybe_gzip "$raw" "$outfile"; then
            return 1
        fi

        rm -f "$raw"

    # --------------------------------------------------------
    # direct genome download URL fallback
    # --------------------------------------------------------

    else

        local raw="${spdir}/.${sp}.genome.download"

        echo "Treating genome_url as a direct file URL."

        if ! fetch_url "$genome_url" "$raw"; then
            return 1
        fi

        if ! copy_maybe_gzip "$raw" "$outfile"; then
            return 1
        fi

        rm -f "$raw"

    fi

    if ! validate_fasta "$outfile"; then
        return 1
    fi

    printf '%s\n' "$genome_url" > "$marker"

    echo "Genome: OK"
}


# ============================================================
# find GFF/GFF3 inside an extracted annotation archive
# ============================================================

find_annotation_file() {

    local dir="$1"
    local candidate=""

    candidate=$(find "$dir" -type f \
        \( -iname '*.gff3' -o -iname '*.gff' \
           -o -iname '*.gff3.gz' -o -iname '*.gff.gz' \) \
        -print -quit)

    printf '%s\n' "$candidate"
}


# ============================================================
# download annotation according to annotation_url
#
# Handles:
#   - plain GFF/GFF3
#   - gzip-compressed GFF/GFF3
#   - tar.gz/tgz archive containing GFF/GFF3
#   - ZIP archive containing GFF/GFF3
# ============================================================

download_annotation() {

    local sp="$1"
    local annotation_url="$2"

    local spdir="${OUTDIR}/${sp}"
    local outfile="${spdir}/${sp}.annotation.gff3"
    local marker="${spdir}/.annotation_source_url"
    local raw="${spdir}/.${sp}.annotation.download"
    local tmpdir="${spdir}/.annotation_extract.$$"
    local candidate=""

    if [[ -z "$annotation_url" ]]; then
        echo "Annotation URL: not provided; skipping."
        return 0
    fi

    if [[ -s "$outfile" && -s "$marker" ]] && \
       [[ "$(cat "$marker")" == "$annotation_url" ]]; then

        echo "Annotation already downloaded from current annotation_url; skipping."
        return 0

    fi

    echo "Downloading annotation according to annotation_url:"
    echo "  ${annotation_url}"

    rm -rf "$tmpdir"
    mkdir -p "$tmpdir"

    if ! fetch_url "$annotation_url" "$raw"; then
        rm -rf "$tmpdir"
        return 1
    fi

    # --------------------------------------------------------
    # tar.gz / tgz archive
    # --------------------------------------------------------

    if tar -tzf "$raw" >/dev/null 2>&1; then

        echo "Annotation resource is a tar.gz archive."

        tar -xzf "$raw" -C "$tmpdir"
        candidate=$(find_annotation_file "$tmpdir")

    # --------------------------------------------------------
    # ZIP archive
    # --------------------------------------------------------

    elif unzip -tqq "$raw" >/dev/null 2>&1; then

        echo "Annotation resource is a ZIP archive."

        unzip -oq "$raw" -d "$tmpdir"
        candidate=$(find_annotation_file "$tmpdir")

    # --------------------------------------------------------
    # gzip-compressed single GFF/GFF3
    # --------------------------------------------------------

    elif gzip -t "$raw" >/dev/null 2>&1; then

        echo "Annotation resource is gzip-compressed."

        gzip -dc "$raw" > "${outfile}.tmp"
        candidate="${outfile}.tmp"

    # --------------------------------------------------------
    # plain annotation file
    # --------------------------------------------------------

    else

        cp "$raw" "${outfile}.tmp"
        candidate="${outfile}.tmp"

    fi

    if [[ -z "$candidate" || ! -s "$candidate" ]]; then
        echo "ERROR: Could not locate GFF/GFF3 for ${sp}" >&2
        rm -f "$raw" "${outfile}.tmp"
        rm -rf "$tmpdir"
        return 1
    fi

    # Archive member may itself be gzip-compressed.
    if [[ "$candidate" != "${outfile}.tmp" ]]; then

        if gzip -t "$candidate" >/dev/null 2>&1; then
            gzip -dc "$candidate" > "${outfile}.tmp"
        else
            cp "$candidate" "${outfile}.tmp"
        fi

    fi

    if ! validate_gff3 "${outfile}.tmp"; then
        rm -f "$raw" "${outfile}.tmp"
        rm -rf "$tmpdir"
        return 1
    fi

    mv "${outfile}.tmp" "$outfile"
    printf '%s\n' "$annotation_url" > "$marker"

    rm -f "$raw"
    rm -rf "$tmpdir"

    echo "Annotation GFF3: OK"
}


# ============================================================
# download protein FASTA according to protein_url
# ============================================================

download_protein() {

    local sp="$1"
    local protein_url="$2"

    local spdir="${OUTDIR}/${sp}"
    local outfile="${spdir}/${sp}.protein.faa"
    local marker="${spdir}/.protein_source_url"
    local raw="${spdir}/.${sp}.protein.download"

    if [[ -z "$protein_url" ]]; then
        echo "Protein URL: not provided; skipping."
        return 0
    fi

    if [[ -s "$outfile" && -s "$marker" ]] && \
       [[ "$(cat "$marker")" == "$protein_url" ]]; then

        echo "Protein FASTA already downloaded from current protein_url; skipping."
        return 0

    fi

    echo "Downloading protein FASTA according to protein_url:"
    echo "  ${protein_url}"

    if ! fetch_url "$protein_url" "$raw"; then
        return 1
    fi

    if ! copy_maybe_gzip "$raw" "$outfile"; then
        return 1
    fi

    rm -f "$raw"

    if ! validate_fasta "$outfile"; then
        rm -f "$outfile"
        return 1
    fi

    printf '%s\n' "$protein_url" > "$marker"

    echo "Protein FASTA: OK"
}


# ============================================================
# download all resources
# ============================================================

echo
echo "============================================================"
echo "starting downloads for synteny analysis"
echo "============================================================"
echo
echo "output directory:"
echo "  ${OUTDIR}"
echo


for entry in "${MANIFEST[@]}"; do

    IFS='|' read -r \
        sp \
        family \
        subfamily \
        assembly_name \
        refseq \
        genome_acc \
        source \
        genome_url \
        annotation_url \
        protein_url <<< "$entry"

    echo
    echo "============================================================"
    echo "Species:        ${sp}"
    echo "Family:         ${family}"
    echo "Subfamily:      ${subfamily:-NA}"
    echo "Assembly:       ${assembly_name}"
    echo "Genome acc.:    ${genome_acc}"
    echo "RefSeq:         ${refseq:-NA}"
    echo "Source:         ${source}"
    echo "Genome URL:     ${genome_url}"
    echo "Annotation URL: ${annotation_url:-NA}"
    echo "Protein URL:    ${protein_url:-NA}"
    echo "============================================================"

    # Genome is required. A genome failure stops the job.
    if ! download_genome \
        "$sp" \
        "$genome_acc" \
        "$genome_url"; then

        echo "ERROR: Genome download/preparation failed for ${sp}" >&2
        exit 1

    fi

    # Annotation/protein are optional according to the CSV.
    # A failed external resource does not stop later species.
    if [[ -n "$annotation_url" ]]; then

        if ! download_annotation \
            "$sp" \
            "$annotation_url"; then

            echo "WARNING: Annotation download/preparation failed for ${sp}" >&2

        fi

    else

        echo "Annotation URL: not provided; skipping."

    fi

    if [[ -n "$protein_url" ]]; then

        if ! download_protein \
            "$sp" \
            "$protein_url"; then

            echo "WARNING: Protein download/preparation failed for ${sp}" >&2

        fi

    else

        echo "Protein URL: not provided; skipping."

    fi

done


# ============================================================
# final inventory
# ============================================================

INVENTORY="${OUTDIR}/assembly_download_inventory.tsv"

printf "species\tfamily\tsubfamily\tassembly_name\tgenome_accession\tNCBI_RefSeq\tsource\tgenome\tannotation_gff3\tprotein_fasta\tannotation_url_provided\tprotein_url_provided\n" \
    > "$INVENTORY"


for entry in "${MANIFEST[@]}"; do

    IFS='|' read -r \
        sp \
        family \
        subfamily \
        assembly_name \
        refseq \
        genome_acc \
        source \
        genome_url \
        annotation_url \
        protein_url <<< "$entry"

    dir="${OUTDIR}/${sp}"

    genome_status="MISSING"
    annotation_status="MISSING"
    protein_status="MISSING"
    annotation_url_status="NO"
    protein_url_status="NO"

    [[ -s "${dir}/${sp}.genome.fa" ]] \
        && genome_status="OK"

    [[ -s "${dir}/${sp}.annotation.gff3" ]] \
        && annotation_status="OK"

    [[ -s "${dir}/${sp}.protein.faa" ]] \
        && protein_status="OK"

    [[ -n "$annotation_url" ]] \
        && annotation_url_status="YES"

    [[ -n "$protein_url" ]] \
        && protein_url_status="YES"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sp" \
        "$family" \
        "$subfamily" \
        "$assembly_name" \
        "$genome_acc" \
        "$refseq" \
        "$source" \
        "$genome_status" \
        "$annotation_status" \
        "$protein_status" \
        "$annotation_url_status" \
        "$protein_url_status" \
        >> "$INVENTORY"

done


# ============================================================
# summary
# ============================================================

echo
echo "============================================================"
echo "DOWNLOAD INVENTORY"
echo "============================================================"
echo

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$INVENTORY"
else
    cat "$INVENTORY"
fi

echo
echo "============================================================"
echo "Download run completed."
echo "============================================================"
echo
echo "Manifest:"
echo "  ${MANIFEST_OUT}"
echo
echo "Inventory:"
echo "  ${INVENTORY}"
echo
echo "Output directory:"
echo "  ${OUTDIR}"
