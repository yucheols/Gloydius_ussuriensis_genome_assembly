#!/bin/bash
#SBATCH --job-name=funannotate
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=672:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================
# 1. Activate environment and set up analyses
# ============================================================

source ~/.bash_profile
conda activate funannotate

# Important: set -u only after conda activation
set -euo pipefail

# Avoid system libstdc++ issue for SignalP / PIL / matplotlib
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# use a short temporary directory to avoid Python AF_UNIX socket path limit
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

# clean it on exit
trap 'rm -rf "$TMPDIR"' EXIT

# ============================================================
# 2. Set paths
# ============================================================

basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation"

# Earl Grey softmasked genome
genome="${basedir}/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# Main tissue RNA-seq reads
main_reads="${basedir}/paired_RNAseq_reads"

# Published venom-gland RNA-seq reads
venom_reads="${basedir}/venom_gland/trimmed_fastq"

# Working directory
workdir="${basedir}/funannotate"

# Output directories
indexdir="${workdir}/hisat2_index"
main_mapdir="${workdir}/rnaseq_main_hisat2"
venom_mapdir="${workdir}/rnaseq_venom_hisat2"
main_stringtie="${workdir}/stringtie_main"
venom_stringtie="${workdir}/stringtie_venom"
sort_tmpdir="${workdir}/tmp"

# Funannotate output
outdir="${workdir}/G_ussuriensis_funannotate"

mkdir -p \
    "$workdir" \
    "$indexdir" \
    "$main_mapdir" \
    "$venom_mapdir" \
    "$main_stringtie" \
    "$venom_stringtie" \
    "$sort_tmpdir"

cd "$workdir"

# ============================================================
# 3. Check environment and inputs
# ============================================================

echo "Checking Funannotate environment..."
funannotate check --show-versions

if [ ! -s "$genome" ]; then
    echo "ERROR: Genome file not found or empty:"
    echo "$genome"
    exit 1
fi

if [ ! -d "$main_reads" ]; then
    echo "ERROR: Main RNA-seq directory not found:"
    echo "$main_reads"
    exit 1
fi

if [ ! -d "$venom_reads" ]; then
    echo "ERROR: Venom-gland RNA-seq directory not found:"
    echo "$venom_reads"
    exit 1
fi

echo
echo "Genome:"
echo "$genome"
echo
echo "Main RNA-seq path:"
echo "$main_reads"
echo
echo "Venom-gland RNA-seq path:"
echo "$venom_reads"
echo
echo "TMPDIR:"
echo "$TMPDIR"
echo

# ============================================================
# 4. Clean and sort genome for Funannotate
# ============================================================

clean_genome="${workdir}/Gloydius_ussuriensis.clean.fa"
sorted_genome="${workdir}/Gloydius_ussuriensis.clean.sorted.fa"

if [ ! -s "$clean_genome" ]; then
    echo "Running funannotate clean..."

    funannotate clean \
        -i "$genome" \
        -o "$clean_genome" \
        --minlen 1000
else
    echo "$clean_genome already exists. Skipping clean."
fi

if [ ! -s "$sorted_genome" ]; then
    echo "Running funannotate sort..."

    funannotate sort \
        -i "$clean_genome" \
        -o "$sorted_genome" \
        --minlen 1000
else
    echo "$sorted_genome already exists. Skipping sort."
fi

annot_genome="$sorted_genome"

# ============================================================
# 5. Build HISAT2 index
# ============================================================

if [ ! -s "${indexdir}/Gus_softmasked.1.ht2" ]; then
    echo "Building HISAT2 index..."

    hisat2-build \
        "$annot_genome" \
        "${indexdir}/Gus_softmasked"
else
    echo "HISAT2 index already exists. Skipping index build."
fi

# ============================================================
# 6. Function to map paired RNA-seq reads
# ============================================================

map_rnaseq_dir () {
    local reads_dir="$1"
    local map_dir="$2"
    local sample_prefix="$3"

    echo
    echo "Mapping RNA-seq reads from:"
    echo "$reads_dir"
    echo

    mkdir -p "$map_dir"

    shopt -s nullglob

    local found_any="no"
    declare -A seen_r1

    for r1 in \
        "$reads_dir"/*_R1*.fq.gz \
        "$reads_dir"/*_R1*.fastq.gz \
        "$reads_dir"/*_1*.fq.gz \
        "$reads_dir"/*_1*.fastq.gz
    do
        [ -e "$r1" ] || continue

        # Avoid duplicate processing if multiple glob patterns overlap
        if [[ -n "${seen_r1[$r1]:-}" ]]; then
            continue
        fi
        seen_r1[$r1]=1

        local r2=""

        if [[ "$r1" == *_R1* ]]; then
            r2="${r1/_R1/_R2}"
        elif [[ "$r1" == *_1* ]]; then
            r2="${r1/_1/_2}"
        else
            echo "WARNING: Could not infer read pattern for:"
            echo "$r1"
            continue
        fi

        if [ ! -s "$r2" ]; then
            echo "WARNING: Could not find matching R2 for:"
            echo "$r1"
            echo "Expected:"
            echo "$r2"
            continue
        fi

        found_any="yes"

        local base
        local sample
        local bam
        local sam

        base=$(basename "$r1")
        base=${base%.fastq.gz}
        base=${base%.fq.gz}

        if [[ "$base" == *_R1* ]]; then
            sample="${base%%_R1*}"
        elif [[ "$base" == *_1* ]]; then
            sample="${base%%_1*}"
        else
            sample="$base"
        fi

        bam="${map_dir}/${sample_prefix}_${sample}.sorted.bam"
        sam="${map_dir}/${sample_prefix}_${sample}.sam"

        if [ -s "$bam" ]; then
            echo "BAM already exists for $sample. Skipping mapping."
            continue
        fi

        echo "Mapping sample: $sample"
        echo "R1: $r1"
        echo "R2: $r2"

        hisat2 \
            -x "${indexdir}/Gus_softmasked" \
            -1 "$r1" \
            -2 "$r2" \
            -p 32 \
            --dta \
            -S "$sam" \
            2> "${map_dir}/${sample_prefix}_${sample}.hisat2.log"

        samtools sort \
            -@ 16 \
            -T "${sort_tmpdir}/${sample_prefix}_${sample}.sorttmp" \
            -o "$bam" \
            "$sam"

        samtools index \
            -@ 8 \
            "$bam"

        samtools flagstat \
            -@ 8 \
            "$bam" \
            > "${map_dir}/${sample_prefix}_${sample}.flagstat.txt"

        rm -f "$sam"

        echo "Finished mapping sample: $sample"
        echo
    done

    shopt -u nullglob

    if [ "$found_any" = "no" ]; then
        echo "ERROR: No paired FASTQ files found in:"
        echo "$reads_dir"
        echo "Expected patterns include *_R1*.fq.gz, *_R1*.fastq.gz, *_1*.fq.gz, or *_1*.fastq.gz"
        exit 1
    fi
}

# ============================================================
# 7. Map main RNA-seq and venom-gland RNA-seq
# ============================================================

map_rnaseq_dir "$main_reads" "$main_mapdir" "main"
map_rnaseq_dir "$venom_reads" "$venom_mapdir" "venom"

# ============================================================
# 8. Function to assemble transcripts with StringTie
# ============================================================

run_stringtie_dir () {
    local map_dir="$1"
    local stringtie_dir="$2"
    local label="$3"
    local merged_gtf="$4"
    local transcript_fa="$5"

    echo
    echo "Running StringTie for $label RNA-seq BAMs..."
    echo

    mkdir -p "$stringtie_dir"

    shopt -s nullglob

    local bam_files=("$map_dir"/*.sorted.bam)

    if [ "${#bam_files[@]}" -eq 0 ]; then
        echo "ERROR: No sorted BAM files found in:"
        echo "$map_dir"
        exit 1
    fi

    for bam in "${bam_files[@]}"
    do
        local sample
        local gtf

        sample=$(basename "$bam" .sorted.bam)
        gtf="${stringtie_dir}/${sample}.gtf"

        if [ ! -s "$gtf" ]; then
            echo "Running StringTie for $sample"

            stringtie "$bam" \
                -p 32 \
                -o "$gtf"
        else
            echo "GTF already exists for $sample. Skipping StringTie."
        fi
    done

    shopt -u nullglob

    ls "$stringtie_dir"/*.gtf > "${stringtie_dir}/${label}_gtf_list.txt"

    if [ ! -s "$merged_gtf" ]; then
        echo "Merging StringTie GTFs for $label..."

        stringtie --merge \
            -p 32 \
            -o "$merged_gtf" \
            "${stringtie_dir}/${label}_gtf_list.txt"
    else
        echo "Merged GTF already exists for $label. Skipping merge."
    fi

    if [ ! -s "$transcript_fa" ]; then
        echo "Extracting transcript FASTA for $label..."

        gffread "$merged_gtf" \
            -g "$annot_genome" \
            -w "$transcript_fa"
    else
        echo "Transcript FASTA already exists for $label. Skipping gffread."
    fi

    if [ ! -s "$transcript_fa" ]; then
        echo "ERROR: Transcript FASTA was not created:"
        echo "$transcript_fa"
        exit 1
    fi

    echo "Transcript FASTA for $label:"
    echo "$transcript_fa"
}

# ============================================================
# 9. Assemble main and venom transcripts
# ============================================================

main_merged_gtf="${main_stringtie}/Gloydius_main_merged.gtf"
main_transcripts="${main_stringtie}/Gloydius_main_rnaseq_transcripts.fa"

venom_merged_gtf="${venom_stringtie}/Gloydius_venom_merged.gtf"
venom_transcripts="${venom_stringtie}/Gloydius_venom_gland_transcripts.fa"

run_stringtie_dir \
    "$main_mapdir" \
    "$main_stringtie" \
    "main" \
    "$main_merged_gtf" \
    "$main_transcripts"

run_stringtie_dir \
    "$venom_mapdir" \
    "$venom_stringtie" \
    "venom" \
    "$venom_merged_gtf" \
    "$venom_transcripts"

# ============================================================
# 10. Combine transcript evidence
# ============================================================

combined_transcripts="${workdir}/Gloydius_main_plus_venom_transcripts.fa"

if [ ! -s "$combined_transcripts" ]; then
    echo
    echo "Combining main + venom transcript evidence..."

    cat "$main_transcripts" "$venom_transcripts" > "$combined_transcripts"
else
    echo
    echo "Combined transcript evidence already exists. Skipping concat."
fi

if [ ! -s "$combined_transcripts" ]; then
    echo "ERROR: Combined transcript evidence file was not created:"
    echo "$combined_transcripts"
    exit 1
fi

n_transcripts=$(grep -c "^>" "$combined_transcripts" || true)

echo
echo "Combined transcript evidence:"
echo "$combined_transcripts"
echo "Number of transcript records: $n_transcripts"

if [ "$n_transcripts" -eq 0 ]; then
    echo "ERROR: Combined transcript FASTA has zero records."
    exit 1
fi

# ============================================================
# 11. Merge RNA BAMs for Funannotate
# ============================================================

merged_rna_bam="${workdir}/Gloydius_all_RNAseq_merged.sorted.bam"

mapfile -t RNA_BAM_ARRAY < <(
    find "$main_mapdir" "$venom_mapdir" -name "*.sorted.bam" | sort
)

if [ "${#RNA_BAM_ARRAY[@]}" -eq 0 ]; then
    echo "ERROR: No RNA-seq BAM files found for --rna_bam"
    exit 1
fi

echo
echo "RNA BAM files to merge for Funannotate:"
printf '%s\n' "${RNA_BAM_ARRAY[@]}"
echo

for bam in "${RNA_BAM_ARRAY[@]}"; do
    if [ ! -s "$bam" ]; then
        echo "ERROR: Missing or empty BAM:"
        echo "$bam"
        exit 1
    fi

    if [ ! -s "${bam}.bai" ]; then
        echo "Indexing BAM:"
        echo "$bam"
        samtools index -@ 8 "$bam"
    fi

    samtools quickcheck "$bam" || {
        echo "ERROR: samtools quickcheck failed for:"
        echo "$bam"
        exit 1
    }
done

if [ ! -s "$merged_rna_bam" ]; then
    echo
    echo "Merging all RNA-seq BAMs into one BAM for Funannotate:"
    echo "$merged_rna_bam"
    echo

    samtools merge \
        -@ 16 \
        -f \
        "$merged_rna_bam" \
        "${RNA_BAM_ARRAY[@]}"

    samtools index -@ 16 "$merged_rna_bam"
else
    echo
    echo "Merged RNA BAM already exists. Skipping merge:"
    echo "$merged_rna_bam"

    if [ ! -s "${merged_rna_bam}.bai" ]; then
        samtools index -@ 16 "$merged_rna_bam"
    fi
fi

samtools quickcheck "$merged_rna_bam" || {
    echo "ERROR: merged RNA BAM failed samtools quickcheck:"
    echo "$merged_rna_bam"
    exit 1
}

echo
echo "Merged RNA BAM evidence for Funannotate:"
echo "$merged_rna_bam"
echo

# ============================================================
# 12. Check Funannotate output directory
# ============================================================

if [ -d "$outdir" ]; then
    echo "ERROR: Funannotate output directory already exists:"
    echo "$outdir"
    echo
    echo "Because a previous run failed partway, remove it before rerunning:"
    echo "rm -rf \"$outdir\""
    exit 1
fi

# ============================================================
# 13. Run Funannotate predict
# ============================================================

echo
echo "Running funannotate predict..."
echo
echo "DIAMOND executable: $(which diamond)"
diamond version
echo

funannotate predict \
    -i "$annot_genome" \
    -o "$outdir" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --name GUS \
    --cpus 32 \
    --busco_db tetrapoda \
    --organism other \
    --max_intronlen 100000 \
    --rna_bam "$merged_rna_bam" \
    --transcript_evidence "$combined_transcripts" \
    --tmpdir "$TMPDIR" \
    --force

echo
echo "Funannotate predict finished."
echo

# ============================================================
# 14. Run Funannotate annotate
# ============================================================

echo
echo "Running funannotate annotate..."
echo

funannotate annotate \
    -i "$outdir" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --cpus 32 \
    --busco_db tetrapoda \
    --database "$FUNANNOTATE_DB" \
    --tmpdir "$TMPDIR" \
    --force

echo
echo "Funannotate annotate finished."
echo

# ============================================================
# 15. Final summary
# ============================================================

echo "Important files:"
echo
echo "Clean sorted genome:"
echo "$annot_genome"
echo
echo "Main transcript evidence:"
echo "$main_transcripts"
echo
echo "Venom-gland transcript evidence:"
echo "$venom_transcripts"
echo
echo "Combined transcript evidence used for prediction:"
echo "$combined_transcripts"
echo
echo "RNA BAM evidence used for prediction:"
echo "$merged_rna_bam"
echo
echo "Funannotate output:"
echo "$outdir"
echo
echo "Done."