## 10) Toxin gene analyses
In this section, we will analyze the microsynteny and phylogenetic relationships for two major toxin components in viperids: snake venom metalloproteases (SVMPs) and phospholipase A2 (PLA2s). We will use five comparison taxa alongside our focal species *G. ussuriensis*. These are: G. shedaoensis, C. adamanteus, C. gasperettii, and N. naja. All species except G. shedaoensis have toxin gene annotations, either as a separate toxin GTF/GFF files, or general genome annotations with toxin coordinate catalogs - we will leverage these resources for these species. For G. shedaoensis, we will use blastn to locate SVMPs and PLA2s and conduct targeted gene model reconstructions.

We will apply an exception to C. adamanteus: The assembly we used for chromosomal macrosynteny was from Hogan et al. (2024). However, Nachtigall et al. (2025) subsequently generated another chromosome-level, haplotype-resolved  assembly. While these two genomes come from the same individual and the same PacBio data, the latter had higher Hi-C coverage. This is therefore preferable to our toxin analyses.   

### Step 0: Set up directory structure and link relevant data
The directory name we will use here is `toxin_microsynteny_v2`. This is a second version because I ended up scraping the initial workflow I built. Run the script below to set up our directory structure and link existing data to the new workflow. This can be run directly on Mendel.
```sh
#!/bin/bash

# ================================================================
# initialize directories and locate local reference files
# ================================================================

PROJECT='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo'
ROOT="${PROJECT}/synteny/toxin_microsynteny_v2"
ASSEMBLIES="${PROJECT}/synteny/assemblies_synteny"
TOXCODAN="${PROJECT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"
FUNANNOTATE="${PROJECT}/annotation/funannotate/annotate_results"

# ------------------------------------------------
# 1. create workflow directories
# ------------------------------------------------
mkdir -p "${ROOT}"/{00_reference,01_query,02_discovery,03_homology,04_curation,05_classification,06_microsynteny}

for SP in G_ussuriensis G_shedaoensis C_adamanteus C_gasperettii N_naja; do
    mkdir -p "${ROOT}/00_reference/${SP}"
done

# ------------------------------------------------
# 2. function for safely linking existing files
# ------------------------------------------------
link_file() {

    local SOURCE="$1"
    local DEST="$2"

    if [[ ! -s "${SOURCE}" ]]; then

        echo "MISSING: ${SOURCE}"

    elif [[ -e "${DEST}" || -L "${DEST}" ]]; then

        echo "EXISTS: ${DEST}"

    else

        ln -s "${SOURCE}" "${DEST}"
        echo "LINKED: ${DEST}"

    fi

}

# ------------------------------------------------
# 3. Gloydius ussuriensis
# ------------------------------------------------
echo
echo "===== G. USSURIENSIS ====="

DEST="${ROOT}/00_reference/G_ussuriensis"

link_file \
    "${PROJECT}/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa" \
    "${DEST}/genome.fa"

link_file \
    "${FUNANNOTATE}/Gloydius_ussuriensis_AMNH_21010.gff3" \
    "${DEST}/general_annotation.gff3"

for FILE in \
    toxin_annotation.gtf \
    toxin_annotation.tsv \
    toxin_annotation_cds.fasta \
    toxin_annotation_pep.fasta \
    annotation_warning.txt \
    annotation_removed.txt; do

    link_file "${TOXCODAN}/${FILE}" "${DEST}/${FILE}"

done

# ------------------------------------------------
# 4. Cerastes gasperettii
# ------------------------------------------------
echo
echo "===== C. GASPERETTII ====="

DEST="${ROOT}/00_reference/C_gasperettii"

link_file \
    "${ASSEMBLIES}/Cerastes_gasperettii/Cerastes_gasperettii.genome.fa" \
    "${DEST}/genome.fa"

link_file \
    "${ASSEMBLIES}/Cerastes_gasperettii/Cerastes_gasperettii.annotation.gff3" \
    "${DEST}/general_annotation.gff3"

# ------------------------------------------------
# 5. Naja naja
# ------------------------------------------------
echo
echo "===== N. NAJA ====="

DEST="${ROOT}/00_reference/N_naja"

link_file \
    "${ASSEMBLIES}/Naja_naja/Naja_naja.genome.fa" \
    "${DEST}/genome.fa"

link_file \
    "${ASSEMBLIES}/Naja_naja/Naja_naja.annotation.gff3" \
    "${DEST}/general_annotation.gff3"

# ------------------------------------------------
# 6. locate Gloydius shedaoensis reference files
# ------------------------------------------------
echo
echo "===== G. SHEDAOENSIS REFERENCE FILES ====="

find "${PROJECT}" \
    \( -type f -o -type l \) \
    \( -iname '*shedaoensis*' -o -iname '*GWHBWDU*' \) \
    -print

# ------------------------------------------------
# 7. locate existing curated toxin annotation files
# ------------------------------------------------
echo
echo "===== G. USSURIENSIS CURATED TOXIN FILES ====="

find "${TOXCODAN}" -maxdepth 2 \
    -type f \
    \( -iname '*.gtf' -o \
       -iname '*.tsv' -o \
       -iname '*.fasta' -o \
       -iname '*warning*' -o \
       -iname '*removed*' -o \
       -iname '*core*' -o \
       -iname '*curated*' \) \
    -print

# ------------------------------------------------
# 8. print new reference directory
# ------------------------------------------------
echo
echo "===== NEW REFERENCE DIRECTORY ====="

ls -lhR "${ROOT}/00_reference"

echo
echo "Initialization complete."
```


### Step 1: Data acquisition