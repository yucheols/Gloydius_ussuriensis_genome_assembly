#!/bin/bash
#SBATCH --job-name=curated_hicmap
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=450G
#SBATCH --time=72:00:00
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# init conda env
source ~/.bash_profile
conda activate scaffolding

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK}"

# working directory for rebuilt Hi-C map
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/hic_rebuild"
mkdir -p "$workdir"
cd "$workdir"

# curated split assembly
FASTA="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"

# path to hic reads
hicdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined"
R1="${hicdir}/Gloydius_ussuriensis_HiC_R1.fastq.gz"
R2="${hicdir}/Gloydius_ussuriensis_HiC_R2.fastq.gz"

# output prefix
prefix="Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split"

# juicer .jar path
JUICER_TOOLS_JAR="/home/yshin/mendel-nas1/juicer_tools/juicer_tools_1.22.01.jar"

# sanity checks
echo "Checking input files: $(date)"
ls -lh "$FASTA"
ls -lh "$R1"
ls -lh "$R2"
ls -lh "$JUICER_TOOLS_JAR"

which bwa-mem2
which samtools
which yahs
which juicer
which java

java -version

echo "Indexing curated FASTA: $(date)"
samtools faidx "$FASTA"

if [ ! -e "${FASTA}.0123" ] && [ ! -e "${FASTA}.bwt.2bit.64" ]; then
    bwa-mem2 index "$FASTA"
fi

# map hic reads to curated split assembly
echo "Mapping Hi-C reads to curated split FASTA: $(date)"
bwa-mem2 mem -5SP -t "$THREADS" "$FASTA" "$R1" "$R2" | \
    samtools view -@ "$THREADS" -b -F 3340 - | \
    samtools sort -@ "$THREADS" -o "${prefix}.sorted.bam" -

samtools index "${prefix}.sorted.bam"

# run yahs
echo "Running YaHS on curated split assembly: $(date)"
yahs "$FASTA" "${prefix}.sorted.bam" -o "$prefix"

# run juicer pre to generate .hic file
echo "Running juicer pre: $(date)"
juicer pre \
    -a \
    -o "${prefix}_JBAT" \
    "${prefix}.bin" \
    "${prefix}_scaffolds_final.agp" \
    "${FASTA}.fai" \
    > "${prefix}_JBAT.log" 2>&1

echo "Making chromosome sizes file: $(date)"
grep PRE_C_SIZE "${prefix}_JBAT.log" | awk '{print $2, $3}' > "${prefix}_JBAT.chrom.sizes"

echo "Chromosome sizes:"
cat "${prefix}_JBAT.chrom.sizes"

echo "Checking JBAT contact file:"
ls -lh "${prefix}_JBAT.txt"
head -n 5 "${prefix}_JBAT.txt"

echo "Creating .hic file: $(date)"
java -Xmx400G -jar "$JUICER_TOOLS_JAR" pre \
    "${prefix}_JBAT.txt" \
    "${prefix}_JBAT.hic" \
    "${prefix}_JBAT.chrom.sizes" \
    > "${prefix}_JBAT.make_hic.stdout" \
    2> "${prefix}_JBAT.make_hic.stderr"

echo "Final files:"
ls -lh "${prefix}_JBAT.hic" "${prefix}_JBAT.assembly" "${prefix}_JBAT.chrom.sizes" "${prefix}_JBAT.txt"

# test .hic readability
echo "Testing .hic readability:"
java -jar "$JUICER_TOOLS_JAR" dump observed NONE \
    "${prefix}_JBAT.hic" \
    assembly assembly BP 1000000 | head \
    > "${prefix}_JBAT.hic_read_test.txt" || true

cat "${prefix}_JBAT.hic_read_test.txt" || true

echo "Done: $(date)"