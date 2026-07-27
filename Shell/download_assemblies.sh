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

# activate download environment
source ~/.bash_profile
conda activate ncbi_datasets

# strict bash error handling
set -euo pipefail

# set working and output directories
WORKDIR=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny
OUTDIR=${WORKDIR}/assemblies_synteny

# make outdir if it's not there
mkdir -p "$OUTDIR"
cd "$OUTDIR"

# NCBI assemblies
species=(
    Crotalus_adamanteus
    Crotalus_viridis
    Bothrops_insularis
    Sistrurus_catenatus
    Vipera_berus
    Cerastes_gasperettii
    Elaphe_schrenckii
    Naja_naja
    Candoia_aspera
    Liasis_olivaceus
    Thamnophis_elegans
)

accessions=(
    GCA_039797435.1
    GCA_003400415.2
    GCA_055824665.1
    GCA_039880765.1
    GCA_964194415.1
    GCA_046524025.1
    GCA_050231175.1
    GCA_009733165.1
    GCA_035149785.1
    GCA_057929265.1
    GCA_009769535.1
)

if [[ ${#species[@]} -ne ${#accessions[@]} ]]; then
    echo "ERROR: Species and accession arrays have different lengths." >&2
    exit 1
fi

# loop through species and download assemblies
for i in "${!species[@]}"; do
    sp=${species[$i]}
    acc=${accessions[$i]}
    spdir=${OUTDIR}/${sp}

    echo "Downloading ${sp}: ${acc}"
    mkdir -p "$spdir"

    # skip assemblies already downloaded successfully
    if [[ -s "${spdir}/${sp}.genome.fa" ]]; then
        echo "Genome already exists; skipping ${sp}"
        continue
    fi

    datasets download genome accession "$acc" \
        --include genome,gff3,protein,cds,seq-report \
        --filename "${spdir}/${acc}.zip" \
        --no-progressbar

    unzip -oq "${spdir}/${acc}.zip" -d "$spdir"

    datadir="${spdir}/ncbi_dataset/data/${acc}"

    if [[ ! -d "$datadir" ]]; then
        echo "ERROR: NCBI data directory not found for ${sp}: ${datadir}" >&2
        exit 1
    fi

    genome=$(find "$datadir" -maxdepth 1 -type f \
        -name '*_genomic.fna' -print -quit)

    if [[ -z "$genome" ]]; then
        echo "ERROR: Genome FASTA not found for ${sp}" >&2
        exit 1
    fi

    ln -sfn "$(readlink -f "$genome")" \
        "${spdir}/${sp}.genome.fa"

    if [[ -f "${datadir}/genomic.gff" ]]; then
        ln -sfn "$(readlink -f "${datadir}/genomic.gff")" \
            "${spdir}/${sp}.annotation.gff3"
    fi

    if [[ -f "${datadir}/protein.faa" ]]; then
        ln -sfn "$(readlink -f "${datadir}/protein.faa")" \
            "${spdir}/${sp}.protein.faa"
    fi

    if [[ -f "${datadir}/cds_from_genomic.fna" ]]; then
        ln -sfn "$(readlink -f "${datadir}/cds_from_genomic.fna")" \
            "${spdir}/${sp}.cds.fna"
    fi

    echo "Finished ${sp}"
done

# NGDC assembly: Gloydius shedaoensis
sp=Gloydius_shedaoensis
acc=GWHBWDU00000000
spdir=${OUTDIR}/${sp}
api_json=${spdir}/${acc}.json

echo "Downloading ${sp}: ${acc}"
mkdir -p "$spdir"

if [[ ! -s "${spdir}/${sp}.genome.fa" ]]; then
    curl -L --fail --retry 5 \
        "https://ngdc.cncb.ac.cn/gwh/api/public/assembly/${acc}" \
        -o "$api_json"

    dna_url=$(jq -r '.ftpPathDna // empty' "$api_json")
    gff_url=$(jq -r '.ftpPathGff // empty' "$api_json")
    protein_url=$(jq -r '.ftpPathProtein // empty' "$api_json")
    cds_url=$(jq -r '.ftpPathCds // empty' "$api_json")

    if [[ -z "$dna_url" ]]; then
        echo "ERROR: NGDC genome download URL was not found." >&2
        exit 1
    fi

    curl -L --fail --retry 5 "$dna_url" -o "${spdir}/${sp}.genome.fa.gz"
    gzip -t "${spdir}/${sp}.genome.fa.gz"
    gzip -dc "${spdir}/${sp}.genome.fa.gz" \
        > "${spdir}/${sp}.genome.fa.tmp"
    mv "${spdir}/${sp}.genome.fa.tmp" \
        "${spdir}/${sp}.genome.fa"

    if [[ -n "$gff_url" ]]; then
        curl -L --fail --retry 5 "$gff_url" -o "${spdir}/${sp}.annotation.gff3.gz"
        gzip -t "${spdir}/${sp}.annotation.gff3.gz"
        gzip -dc "${spdir}/${sp}.annotation.gff3.gz" \
            > "${spdir}/${sp}.annotation.gff3.tmp"
        mv "${spdir}/${sp}.annotation.gff3.tmp" \
            "${spdir}/${sp}.annotation.gff3"
    fi

    if [[ -n "$protein_url" ]]; then
        curl -L --fail --retry 5 "$protein_url" -o "${spdir}/${sp}.protein.faa.gz"
        gzip -t "${spdir}/${sp}.protein.faa.gz"
        gzip -dc "${spdir}/${sp}.protein.faa.gz" \
            > "${spdir}/${sp}.protein.faa.tmp"
        mv "${spdir}/${sp}.protein.faa.tmp" \
            "${spdir}/${sp}.protein.faa"
    fi

    if [[ -n "$cds_url" ]]; then
        curl -L --fail --retry 5 "$cds_url" -o "${spdir}/${sp}.cds.fna.gz"
        gzip -t "${spdir}/${sp}.cds.fna.gz"
        gzip -dc "${spdir}/${sp}.cds.fna.gz" \
            > "${spdir}/${sp}.cds.fna.tmp"
        mv "${spdir}/${sp}.cds.fna.tmp" \
            "${spdir}/${sp}.cds.fna"
    fi
else
    echo "Genome already exists; skipping ${sp}"
fi

echo
echo "All downloads completed."
echo "Output directory: ${OUTDIR}"

# simple final check
for dir in "${OUTDIR}"/*; do
    [[ -d "$dir" ]] || continue
    sp=$(basename "$dir")

    if [[ -s "${dir}/${sp}.genome.fa" ]]; then
        echo "${sp}: genome OK"
    else
        echo "${sp}: genome MISSING"
    fi
done