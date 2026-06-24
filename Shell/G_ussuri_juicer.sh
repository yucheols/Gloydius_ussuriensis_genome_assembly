#!/bin/bash
#SBATCH --job-name=juicer
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# stop if anything fails
set -euo pipefail

# set dir
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/yahs_out"
prefix="Gloydius_ussuriensis_AMNH_21010_yahs"

# set vars
BIN="${workdir}/${prefix}.bin"
AGP="${workdir}/${prefix}_scaffolds_final.agp"
FASTA="${workdir}/${prefix}_scaffolds_final.fa"

# IMPORTANT: use the original PacBio draft FASTA index here,
# not the final YaHS scaffolded FASTA index
FAI="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa.fai"

JBAT="${workdir}/${prefix}_JBAT"

# standalone juicer_tools jar
JUICER_TOOLS_JAR="/home/yshin/mendel-nas1/juicer_tools/juicer_tools_1.22.01.jar"

# cd into working dir
cd "$workdir"

# sanity checks
echo "Checking files: $(date)"
ls -lh "$BIN"
ls -lh "$AGP"
ls -lh "$FASTA"
ls -lh "$FAI"
ls -lh "$JUICER_TOOLS_JAR"

which juicer
which java
java -version
which samtools

echo "Using juicer tools jar:"
echo "$JUICER_TOOLS_JAR"

# remove previous incomplete hic only
rm -f "${JBAT}.hic"

# run juicer pre from YaHS
echo "Running juicer pre: $(date)"
juicer pre \
    -a \
    -o "$JBAT" \
    "$BIN" \
    "$AGP" \
    "$FAI" \
    > "${JBAT}.log" 2>&1

echo "Finished juicer pre: $(date)"

# make chromosome sizes file from juicer pre log
echo "Making chromosome sizes file from juicer pre log: $(date)"
grep PRE_C_SIZE "${JBAT}.log" | awk '{print $2, $3}' > "${JBAT}.chrom.sizes"

echo "Chromosome sizes:"
cat "${JBAT}.chrom.sizes"

# check JBAT outputs
echo "Checking JBAT outputs:"
ls -lh "${JBAT}"*

# make sure JBAT txt exists and has contacts
echo "Checking JBAT contact file:"
ls -lh "${JBAT}.txt"
head -n 5 "${JBAT}.txt"

# run juicer tools directly with java
echo "Creating .hic file: $(date)"

java -Xmx80G -jar "$JUICER_TOOLS_JAR" pre \
    "${JBAT}.txt" \
    "${JBAT}.hic" \
    "${JBAT}.chrom.sizes" \
    > "${JBAT}.make_hic.stdout" \
    2> "${JBAT}.make_hic.stderr"

echo "Finished creating .hic file: $(date)"

# final check
echo "Final JBAT files:"
ls -lh "${JBAT}.hic" "${JBAT}.assembly" "${JBAT}.chrom.sizes" "${JBAT}.txt"

echo "Testing whether .hic is readable:"
java -jar "$JUICER_TOOLS_JAR" dump observed NONE \
    "${JBAT}.hic" \
    assembly assembly BP 1000000 | head \
    > "${JBAT}.hic_read_test.txt"

echo "Read test output:"
cat "${JBAT}.hic_read_test.txt"

echo "Done: $(date)"