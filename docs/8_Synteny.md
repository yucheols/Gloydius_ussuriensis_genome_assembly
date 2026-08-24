## 8) Chromosomal synteny
### Step 1: Setup directory and prepare data downloads
First, set up a directory for synteny analyses.
```sh
# under the "G_ussuriensis_Chromo" directory
mkdir -p synteny/{raw_ncbi,assemblies,gff3,proteins,metadata,logs,scripts}
```

Activate the conda env to access NCBI Datasets CLI
```sh
conda activate ncbi_datasets
```

### Step 2: Download assemblies
There are several chromosom-level snake reference genome assemblies, especially viperid assemblies, available. We can use these to investigate the synteny across species. We already downloaded *C. adamanteus*, *C. viridis*, and *V. berus* assemblies. We will download nine additional assemblies. These additional species are: *Cerastes gasperettii*, *Sistrurus catenatus*, *Bothrops insularis*, *GLoydius shedaoensis*, *Naja naja*, *Thamnophis elegans*, *Elaphe schrenckii*, *Liasis olivaceus*, *Candoia aspera*.

First, create a .csv file containing species name, family/subfamily names, assembly name, assembly accession, data source, and url. Note that all assemblies other than *G. shedaoensis* are deposited GenBank. The *G. shedaoensis* assembly is available from National Genomics Data Center (NGDC) Genome Warehouse (GWH).

Store the .csv file in the "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/metadata" with the following name:
```sh
synteny_samples_list.csv
```

For the NCBI data download, we can use the "ncbi_datasets" conda env. But, to download data from NGDC, we need to use curl and jq. Download these two commands into the "ncbi_datasets" conda env.
```sh
conda install -n ncbi_datasets -c conda-forge \
    curl jq -y
```

Next, submit the following shell script to Mendel. This script will take the .csv file and download the assemblies that will be used in synteny analyses. 
```sh
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
# download chromosome-level snake assemblies for synteny
#
# what this script does:
#   1. download genome assemblies
#   2. download GFF3/protein/CDS files when available
#   3. create standardized filenames
#   4. tolerate assemblies with no public annotation
#   5. create an inventory of downloaded files
#   6. list species that will require structural annotation
#
# species/accessions are based on synteny_samples_list(3).csv
#
# missing annotations are NOT considered download failures.
# ============================================================


# ============================================================
# activate conda environment
# ============================================================

source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate ncbi_datasets

set -euo pipefail


# ============================================================
# set up working directories
# ============================================================

WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny"
OUTDIR="${WORKDIR}/assemblies_synteny"

mkdir -p "$OUTDIR"
cd "$OUTDIR"


# ============================================================
# check required programs
# ============================================================

for cmd in datasets unzip curl gzip jq find readlink; do

    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $cmd" >&2
        exit 1
    fi

done


# ============================================================
# assembly manifest
#
# fields:
#   species
#   family
#   subfamily
#   assembly_name
#   assembly_accession
#   downloaded_status_from_csv
#   source
#   url
#
# NA = no subfamily given in supplied CSV
# ============================================================

readarray -t MANIFEST <<'EOF'
Crotalus_adamanteus|Viperidae|Crotalinae|Cadamanteus_3dDNAHiC_1.2|GCA_039797435.1|Y|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_039797435.1/
Crotalus_viridis|Viperidae|Crotalinae|UTA_CroVir_3.0|GCA_003400415.2|Y|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_003400415.2/
Bothrops_insularis|Viperidae|Crotalinae|ASM5582466v1|GCA_055824665.1|N|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_055824665.1/
Gloydius_shedaoensis|Viperidae|Crotalinae|Gshe1|GWHBWDU00000000|N|NGDC|https://ngdc.cncb.ac.cn/gwh/Assembly/34318/show
Vipera_berus|Viperidae|Viperinae|rVipBer3.hap1.1|GCA_964194415.1|Y|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_964194415.1/
Cerastes_gasperettii|Viperidae|Viperinae|ASM4652402v1|GCA_046524025.1|N|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_046524025.1/
Naja_naja|Elapidae|Elapinae|Nana_v5|GCA_009733165.1|N|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009733165.1/
Elaphe_schrenckii|Colubridae|Colubrinae|ASM5023117v1|GCA_050231175.1|N|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_050231175.1/
Candoia_aspera|Boidae|Candoiinae|rCanAsp1.hap2|GCA_035149785.1|N|GenBank|https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_035149785.1/
Xenopeltis_unicolor|Xenopeltidae|NA|Xuni1|GWHBWDW00000000|N|NGDC|https://ngdc.cncb.ac.cn/gwh/Assembly/34320/show
Argyrophis_diardii|Typhlopidae|NA|Adia1|GWHBWDT00000000|N|NGDC|https://ngdc.cncb.ac.cn/gwh/Assembly/34317/show
EOF


# ============================================================
# save manifest used by this run
# ============================================================

MANIFEST_OUT="${OUTDIR}/assembly_manifest.tsv"

printf \
"species\tfamily\tsubfamily\tassembly_name\tassembly_accession\toriginal_downloaded_status\tsource\turl\n" \
> "$MANIFEST_OUT"


for entry in "${MANIFEST[@]}"; do

    IFS='|' read -r \
        sp \
        family \
        subfamily \
        assembly_name \
        acc \
        downloaded \
        source \
        url <<< "$entry"

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sp" \
        "$family" \
        "$subfamily" \
        "$assembly_name" \
        "$acc" \
        "$downloaded" \
        "$source" \
        "$url" \
        >> "$MANIFEST_OUT"

done


# ============================================================
# standardize NCBI files
# ============================================================

standardize_ncbi_files() {

    local sp="$1"
    local acc="$2"
    local spdir="$3"
    local datadir="$4"

    local genome=""

    # --------------------------------------------------------
    # genome
    # --------------------------------------------------------

    genome=$(find "$datadir" \
        -maxdepth 1 \
        -type f \
        \( -name '*_genomic.fna' -o -name 'genomic.fna' \) \
        -print -quit)

    if [[ -z "${genome:-}" ]]; then
        echo "ERROR: Genome FASTA not found for ${sp}" >&2
        exit 1
    fi

    ln -sfn \
        "$(readlink -f "$genome")" \
        "${spdir}/${sp}.genome.fa"

    echo "Genome: OK"


    # --------------------------------------------------------
    # GFF3
    # --------------------------------------------------------

    if [[ -s "${datadir}/genomic.gff" ]]; then

        ln -sfn \
            "$(readlink -f "${datadir}/genomic.gff")" \
            "${spdir}/${sp}.annotation.gff3"

        echo "GFF3: OK"

    else

        echo "GFF3: not available"

    fi


    # --------------------------------------------------------
    # protein FASTA
    # --------------------------------------------------------

    if [[ -s "${datadir}/protein.faa" ]]; then

        ln -sfn \
            "$(readlink -f "${datadir}/protein.faa")" \
            "${spdir}/${sp}.protein.faa"

        echo "Protein FASTA: OK"

    else

        echo "Protein FASTA: not available"

    fi


    # --------------------------------------------------------
    # CDS FASTA
    # --------------------------------------------------------

    if [[ -s "${datadir}/cds_from_genomic.fna" ]]; then

        ln -sfn \
            "$(readlink -f "${datadir}/cds_from_genomic.fna")" \
            "${spdir}/${sp}.cds.fna"

        echo "CDS FASTA: OK"

    else

        echo "CDS FASTA: not available"

    fi
}


# ============================================================
# download NCBI assembly
# ============================================================

download_ncbi() {

    local sp="$1"
    local acc="$2"

    local spdir="${OUTDIR}/${sp}"
    local zipfile="${spdir}/${acc}.zip"
    local checked="${spdir}/.ncbi_repository_checked"
    local datadir="${spdir}/ncbi_dataset/data/${acc}"

    mkdir -p "$spdir"

    echo
    echo "============================================================"
    echo "Species:   ${sp}"
    echo "Accession: ${acc}"
    echo "Source:    NCBI GenBank"
    echo "============================================================"


    # --------------------------------------------------------
    # skip if repository was previously queried successfully
    # --------------------------------------------------------

    if [[ -s "${spdir}/${sp}.genome.fa" && -e "$checked" ]]; then

        echo "Repository already checked and genome exists."
        echo "Skipping NCBI download."

        return

    fi


    # --------------------------------------------------------
    # download NCBI package
    # --------------------------------------------------------

    echo "Downloading NCBI dataset..."

    rm -f "${zipfile}.tmp"

    datasets download genome accession "$acc" \
        --include genome,gff3,protein,cds,seq-report \
        --filename "${zipfile}.tmp" \
        --no-progressbar


    if [[ ! -s "${zipfile}.tmp" ]]; then

        echo "ERROR: NCBI download failed for ${sp}" >&2
        exit 1

    fi

    mv "${zipfile}.tmp" "$zipfile"


    # --------------------------------------------------------
    # extract
    # --------------------------------------------------------

    unzip -oq "$zipfile" -d "$spdir"


    if [[ ! -d "$datadir" ]]; then

        echo "ERROR: Expected NCBI directory not found:" >&2
        echo "  ${datadir}" >&2
        exit 1

    fi


    # --------------------------------------------------------
    # standardize filenames
    # --------------------------------------------------------

    standardize_ncbi_files \
        "$sp" \
        "$acc" \
        "$spdir" \
        "$datadir"


    # mark repository as checked even when annotation is absent
    touch "$checked"

    echo "Finished ${sp}"
}


# ============================================================
# download and decompress an NGDC file
#
# handles either gzip-compressed or uncompressed downloads
# ============================================================

download_ngdc_file() {

    local url="$1"
    local outfile="$2"

    local raw="${outfile}.download"
    local tmp="${outfile}.tmp"

    rm -f "$raw" "$tmp"

    curl \
        -L \
        --fail \
        --retry 5 \
        --retry-delay 5 \
        "$url" \
        -o "$raw"


    if [[ ! -s "$raw" ]]; then

        echo "ERROR: Empty NGDC download:" >&2
        echo "  ${url}" >&2
        exit 1

    fi


    # determine whether downloaded file is gzip-compressed
    if gzip -t "$raw" >/dev/null 2>&1; then

        gzip -dc "$raw" > "$tmp"

    else

        cp "$raw" "$tmp"

    fi


    if [[ ! -s "$tmp" ]]; then

        echo "ERROR: Failed to create ${outfile}" >&2
        exit 1

    fi


    mv "$tmp" "$outfile"

    rm -f "$raw"
}


# ============================================================
# download assemblies from NGDC
# ============================================================

download_ngdc() {

    local sp="$1"
    local acc="$2"

    local spdir="${OUTDIR}/${sp}"
    local api_json="${spdir}/${acc}.json"
    local checked="${spdir}/.ngdc_repository_checked"

    mkdir -p "$spdir"

    echo
    echo "============================================================"
    echo "Species:   ${sp}"
    echo "Accession: ${acc}"
    echo "Source:    NGDC"
    echo "============================================================"


    # --------------------------------------------------------
    # skip after successful previous repository check
    # --------------------------------------------------------

    if [[ -s "${spdir}/${sp}.genome.fa" && -e "$checked" ]]; then

        echo "Repository already checked and genome exists."
        echo "Skipping NGDC download."

        return

    fi


    # --------------------------------------------------------
    # query NGDC API
    # --------------------------------------------------------

    echo "Querying NGDC API..."

    curl \
        -L \
        --fail \
        --retry 5 \
        --retry-delay 5 \
        "https://ngdc.cncb.ac.cn/gwh/api/public/assembly/${acc}" \
        -o "${api_json}.tmp"


    if [[ ! -s "${api_json}.tmp" ]]; then

        echo "ERROR: Empty NGDC API response for ${sp}" >&2
        exit 1

    fi

    mv "${api_json}.tmp" "$api_json"


    # --------------------------------------------------------
    # retrieve URLs
    # --------------------------------------------------------

    dna_url=$(jq -r '.ftpPathDna // empty' "$api_json")
    gff_url=$(jq -r '.ftpPathGff // empty' "$api_json")
    protein_url=$(jq -r '.ftpPathProtein // empty' "$api_json")
    cds_url=$(jq -r '.ftpPathCds // empty' "$api_json")


    # --------------------------------------------------------
    # genome
    # --------------------------------------------------------

    if [[ -z "$dna_url" ]]; then

        echo "ERROR: NGDC genome URL not found for ${sp}" >&2
        exit 1

    fi


    echo "Downloading genome..."

    download_ngdc_file \
        "$dna_url" \
        "${spdir}/${sp}.genome.fa"

    echo "Genome: OK"


    # --------------------------------------------------------
    # GFF3
    # --------------------------------------------------------

    if [[ -n "$gff_url" ]]; then

        echo "Downloading GFF3..."

        download_ngdc_file \
            "$gff_url" \
            "${spdir}/${sp}.annotation.gff3"

        echo "GFF3: OK"

    else

        echo "GFF3: not available"

    fi


    # --------------------------------------------------------
    # protein FASTA
    # --------------------------------------------------------

    if [[ -n "$protein_url" ]]; then

        echo "Downloading protein FASTA..."

        download_ngdc_file \
            "$protein_url" \
            "${spdir}/${sp}.protein.faa"

        echo "Protein FASTA: OK"

    else

        echo "Protein FASTA: not available"

    fi


    # --------------------------------------------------------
    # CDS FASTA
    # --------------------------------------------------------

    if [[ -n "$cds_url" ]]; then

        echo "Downloading CDS FASTA..."

        download_ngdc_file \
            "$cds_url" \
            "${spdir}/${sp}.cds.fna"

        echo "CDS FASTA: OK"

    else

        echo "CDS FASTA: not available"

    fi


    touch "$checked"

    echo "Finished ${sp}"
}


# ============================================================
# download all assemblies
# ============================================================

echo
echo "============================================================"
echo "starting assembly downloads to be used for synteny analysis"
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
        acc \
        downloaded \
        source \
        url <<< "$entry"

    echo
    echo "------------------------------------------------------------"
    echo "Species:       ${sp}"
    echo "Family:        ${family}"
    echo "Subfamily:     ${subfamily}"
    echo "Assembly:      ${assembly_name}"
    echo "Accession:     ${acc}"
    echo "CSV status:    ${downloaded}"
    echo "Source:        ${source}"
    echo "Record URL:    ${url}"
    echo "------------------------------------------------------------"


    case "$source" in

        GenBank)

            download_ncbi \
                "$sp" \
                "$acc"
            ;;


        NGDC)

            download_ngdc \
                "$sp" \
                "$acc"
            ;;


        *)

            echo "ERROR: Unknown source '${source}' for ${sp}" >&2
            exit 1
            ;;

    esac

done


# ============================================================
# final inventory
# ============================================================

INVENTORY="${OUTDIR}/assembly_download_inventory.tsv"
NEED_ANNOT="${OUTDIR}/assemblies_requiring_annotation.txt"


printf \
"species\tfamily\tsubfamily\tassembly_name\taccession\tsource\tgenome\tgff3\tprotein\tcds\n" \
> "$INVENTORY"


: > "$NEED_ANNOT"


for entry in "${MANIFEST[@]}"; do

    IFS='|' read -r \
        sp \
        family \
        subfamily \
        assembly_name \
        acc \
        downloaded \
        source \
        url <<< "$entry"

    dir="${OUTDIR}/${sp}"

    genome_status="MISSING"
    gff_status="MISSING"
    protein_status="MISSING"
    cds_status="MISSING"


    [[ -s "${dir}/${sp}.genome.fa" ]] \
        && genome_status="OK"

    [[ -s "${dir}/${sp}.annotation.gff3" ]] \
        && gff_status="OK"

    [[ -s "${dir}/${sp}.protein.faa" ]] \
        && protein_status="OK"

    [[ -s "${dir}/${sp}.cds.fna" ]] \
        && cds_status="OK"


    printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sp" \
        "$family" \
        "$subfamily" \
        "$assembly_name" \
        "$acc" \
        "$source" \
        "$genome_status" \
        "$gff_status" \
        "$protein_status" \
        "$cds_status" \
        >> "$INVENTORY"


    # GENESPACE ultimately needs gene coordinates + proteins.
    # Flag species lacking either GFF3 or protein FASTA.
    if [[ "$gff_status" != "OK" || "$protein_status" != "OK" ]]; then

        echo "$sp" >> "$NEED_ANNOT"

    fi

done


# ============================================================
# print summary
# ============================================================

echo
echo "============================================================"
echo "DOWNLOAD INVENTORY"
echo "============================================================"
echo


column -t -s $'\t' "$INVENTORY" || cat "$INVENTORY"


echo
echo "============================================================"
echo "ASSEMBLIES REQUIRING STRUCTURAL ANNOTATION"
echo "============================================================"
echo


if [[ -s "$NEED_ANNOT" ]]; then

    cat "$NEED_ANNOT"

else

    echo "None"

fi


echo
echo "============================================================"
echo "All download attempts completed."
echo "============================================================"
echo
echo "Assembly manifest:"
echo "  ${MANIFEST_OUT}"
echo
echo "Download inventory:"
echo "  ${INVENTORY}"
echo
echo "Assemblies requiring annotation:"
echo "  ${NEED_ANNOT}"
echo
echo "Output directory:"
echo "  ${OUTDIR}"
```

### Step 3: Download RNA-seq data


### Step 4: Install software
We will create several conda env for handling various steps in the synteny analyses.

##### 1) General genome data handling & qc env
```sh
conda create -n synteny_qc \
    -c conda-forge -c bioconda \
    seqkit \
    samtools \
    bedtools \
    gffread \
    agat \
    compleasm \
    minimap2 \
    mummer4
```
##### 2) Repeat masking env
```sh
conda create -n repeatmask \
    -c conda-forge -c bioconda \
    repeatmodeler \
    repeatmasker
```    

##### 3) RNA-seq mapping env
```sh
conda create -n rnaseq_map \
    -c conda-forge -c bioconda \
    fastp \
    star \
    samtools \
    stringtie \
    gffread
```

##### 4) GENESPACE env
```sh
conda create -n genespace \
    -c conda-forge -c bioconda \
    r-base \
    r-devtools \
    r-remotes \
    r-ggplot2 \
    r-igraph \
    r-dbscan \
    r-r.utils \
    r-data.table \
    bioconductor-biostrings \
    bioconductor-rtracklayer \
    orthofinder \
    diamond \
    git \
    make \
    gcc_linux-64 \
    gxx_linux-64
```

Check installation:
```sh
orthofinder --version
diamond version
R --version
```

After this is done, activate this conda env and start R:
```sh
conda activate genespace
R
```
Then install GENESPACE in it:
```R
if (!requireNamespace('devtools', quietly = T))
    install.packages('devtools')

devtools::install_github('jtlovell/GENESPACE')

library(GENESPACE)
citation('GENESPACE')
```

#####  BRAKER3 env
BRAKER3 was already installed in a conda env called "braker3" prior to this project. Check installation:
```sh
conda activate braker3

braker.pl --version
augustus --version
diamond version
samtools --version
```
