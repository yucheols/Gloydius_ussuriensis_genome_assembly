#!/bin/bash
#SBATCH --job-name=funannotate
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250G
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================
# 1. Activate environment
# ============================================================

source ~/.bash_profile
conda activate funannotate
export FUNANNOTATE_DB=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db

set -euo pipefail

# ============================================================
# 2. Set paths
# ============================================================
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation"

# Earl Grey softmasked genome
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# Main tissue RNA-seq reads
main_reads="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads"

# Published venom-gland RNA-seq reads
venom_reads="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq"

# Working directory
workdir="${basedir}/funannotate"

# Output directories
indexdir="${workdir}/hisat2_index"
main_mapdir="${workdir}/rnaseq_main_hisat2"
venom_mapdir="${workdir}/rnaseq_venom_hisat2"
main_stringtie="${workdir}/stringtie_main"
venom_stringtie="${workdir}/stringtie_venom"
outdir="${workdir}/G_ussuriensis_funannotate"

mkdir -p "$workdir" "$indexdir" "$main_mapdir" "$venom_mapdir" "$main_stringtie" "$venom_stringtie"

cd "$workdir"

# ============================================================
# 3. Check inputs
# ============================================================

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

echo "Genome:"
echo "$genome"
echo
echo "Main RNA-seq path:"
echo "$main_reads"
echo
echo "Venom-gland RNA-seq path:"
echo "$venom_reads"
echo

# ============================================================
# 4. Clean and sort genome for funannotate
# ============================================================

if [ ! -s Gloydius_ussuriensis.clean.fa ]; then
    echo "Running funannotate clean..."

    funannotate clean \
        -i "$genome" \
        -o Gloydius_ussuriensis.clean.fa \
        --minlen 1000
else
    echo "Gloydius_ussuriensis.clean.fa already exists. Skipping clean."
fi

if [ ! -s Gloydius_ussuriensis.clean.sorted.fa ]; then
    echo "Running funannotate sort..."

    funannotate sort \
        -i Gloydius_ussuriensis.clean.fa \
        -o Gloydius_ussuriensis.clean.sorted.fa \
        --minlen 1000
else
    echo "Gloydius_ussuriensis.clean.sorted.fa already exists. Skipping sort."
fi

annot_genome="${workdir}/Gloydius_ussuriensis.clean.sorted.fa"

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
    reads_dir="$1"
    map_dir="$2"
    sample_prefix="$3"

    echo "Mapping RNA-seq reads from:"
    echo "$reads_dir"

    # This assumes files contain _R1 and _R2.
    # Example:
    # sample_R1.fq.gz / sample_R2.fq.gz
    # sample_R1.trim.fq.gz / sample_R2.trim.fq.gz

    found_any="no"

    for r1 in "$reads_dir"/*_R1*.fq.gz "$reads_dir"/*_R1*.fastq.gz "$reads_dir"/*_1*.fq.gz "$reads_dir"/*_1*.fastq.gz
    do
        [ -e "$r1" ] || continue

        # Try to infer R2 name
        r2="${r1/_R1/_R2}"
        r2="${r2/_1/_2}"

        if [ ! -s "$r2" ]; then
            echo "WARNING: Could not find matching R2 for:"
            echo "$r1"
            echo "Expected:"
            echo "$r2"
            continue
        fi

        found_any="yes"

        sample=$(basename "$r1")
        sample=${sample%%.fastq.gz}
        sample=${sample%%.fq.gz}
        sample=${sample/_R1*/}
        sample=${sample/_1*/}

        bam="${map_dir}/${sample_prefix}_${sample}.sorted.bam"

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
            -S "${map_dir}/${sample_prefix}_${sample}.sam" \
            2> "${map_dir}/${sample_prefix}_${sample}.hisat2.log"

        samtools sort \
            -@ 16 \
            -o "$bam" \
            "${map_dir}/${sample_prefix}_${sample}.sam"

        samtools index "$bam"

        samtools flagstat "$bam" \
            > "${map_dir}/${sample_prefix}_${sample}.flagstat.txt"

        rm "${map_dir}/${sample_prefix}_${sample}.sam"
    done

    if [ "$found_any" = "no" ]; then
        echo "ERROR: No paired FASTQ files found in:"
        echo "$reads_dir"
        echo "Check file naming. Expected patterns include *_R1*.fq.gz or *_1*.fq.gz"
        exit 1
    fi
}

# ============================================================
# 7. Map main RNA-seq and venom-gland RNA-seq
# ============================================================

map_rnaseq_dir "$main_reads" "$main_mapdir" "main"
map_rnaseq_dir "$venom_reads" "$venom_mapdir" "venom"

# ============================================================
# 8. Assemble main RNA-seq transcripts with StringTie
# ============================================================

echo "Running StringTie on main RNA-seq BAMs..."

for bam in "$main_mapdir"/*.sorted.bam
do
    sample=$(basename "$bam" .sorted.bam)
    gtf="${main_stringtie}/${sample}.gtf"

    if [ ! -s "$gtf" ]; then
        stringtie "$bam" \
            -p 32 \
            -o "$gtf"
    fi
done

ls "$main_stringtie"/*.gtf > "${main_stringtie}/main_gtf_list.txt"

if [ ! -s "${main_stringtie}/Gloydius_main_merged.gtf" ]; then
    stringtie --merge \
        -p 32 \
        -o "${main_stringtie}/Gloydius_main_merged.gtf" \
        "${main_stringtie}/main_gtf_list.txt"
fi

if [ ! -s "${main_stringtie}/Gloydius_main_rnaseq_transcripts.fa" ]; then
    gffread "${main_stringtie}/Gloydius_main_merged.gtf" \
        -g "$annot_genome" \
        -w "${main_stringtie}/Gloydius_main_rnaseq_transcripts.fa"
fi

main_transcripts="${main_stringtie}/Gloydius_main_rnaseq_transcripts.fa"

# ============================================================
# 9. Assemble venom-gland RNA-seq transcripts with StringTie
# ============================================================

echo "Running StringTie on venom-gland RNA-seq BAMs..."

for bam in "$venom_mapdir"/*.sorted.bam
do
    sample=$(basename "$bam" .sorted.bam)
    gtf="${venom_stringtie}/${sample}.gtf"

    if [ ! -s "$gtf" ]; then
        stringtie "$bam" \
            -p 32 \
            -o "$gtf"
    fi
done

ls "$venom_stringtie"/*.gtf > "${venom_stringtie}/venom_gtf_list.txt"

if [ ! -s "${venom_stringtie}/Gloydius_venom_merged.gtf" ]; then
    stringtie --merge \
        -p 32 \
        -o "${venom_stringtie}/Gloydius_venom_merged.gtf" \
        "${venom_stringtie}/venom_gtf_list.txt"
fi

if [ ! -s "${venom_stringtie}/Gloydius_venom_gland_transcripts.fa" ]; then
    gffread "${venom_stringtie}/Gloydius_venom_merged.gtf" \
        -g "$annot_genome" \
        -w "${venom_stringtie}/Gloydius_venom_gland_transcripts.fa"
fi

venom_transcripts="${venom_stringtie}/Gloydius_venom_gland_transcripts.fa"

# ============================================================
# 10. Check transcript evidence files
# ============================================================

if [ ! -s "$main_transcripts" ]; then
    echo "ERROR: Main transcript evidence file was not created:"
    echo "$main_transcripts"
    exit 1
fi

if [ ! -s "$venom_transcripts" ]; then
    echo "ERROR: Venom transcript file was not created:"
    echo "$venom_transcripts"
    exit 1
fi

echo "Main transcript evidence:"
echo "$main_transcripts"

echo "Venom-gland transcript evidence:"
echo "$venom_transcripts"

# ============================================================
# 11. Run funannotate predict
# ============================================================

if [ -d "$outdir" ]; then
    echo "ERROR: funannotate output directory already exists:"
    echo "$outdir"
    echo "Remove it or change outdir if you want to rerun."
    exit 1
fi

echo "Running funannotate predict using MAIN RNA-seq transcript evidence..."

funannotate predict \
    -i "$annot_genome" \
    -o "$outdir" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --cpus 32 \
    --busco_db tetrapoda \
    --organism other \
    --transcript_evidence "$main_transcripts" \
    --force

echo "Funannotate predict finished."

echo
echo "Important files:"
echo "Main RNA-seq transcript evidence:"
echo "$main_transcripts"
echo
echo "Venom-gland transcript evidence, for later venom-gene curation:"
echo "$venom_transcripts"
echo
echo "Funannotate output:"
echo "$outdir"