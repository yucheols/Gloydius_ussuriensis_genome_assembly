#!/bin/bash
#SBATCH --job-name=iprscan
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# 1. environment
# ============================================================

source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# ============================================================
# 2. paths
# ============================================================

workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
outdir="${workdir}/G_ussuriensis_funannotate"
iprdir="/home/yshin/mendel-nas1/interproscan/interproscan-5.78-109.0"
iprscan="${iprdir}/interproscan.sh"

# ============================================================
# 3. validate InterProScan
# ============================================================

if [ ! -x "$iprscan" ]; then
    echo "ERROR: InterProScan executable not found:"
    echo "$iprscan"
    exit 1
fi

echo
echo "InterProScan:"
"$iprscan" -version

echo
echo "Java:"
java -version

echo
echo "Python:"
python3 --version


# ============================================================
# 4. validate funannotate predict output
# ============================================================

if [ ! -d "$outdir/predict_results" ]; then
    echo "ERROR: Funannotate predict_results directory not found:"
    echo "$outdir/predict_results"
    exit 1
fi

proteins=$(find "$outdir/predict_results" \
    -maxdepth 1 \
    -type f \
    -name "*.proteins.fa" \
    | head -n 1)

if [ -z "$proteins" ] || [ ! -s "$proteins" ]; then
    echo "ERROR: Funannotate protein FASTA not found."
    exit 1
fi

echo
echo "Protein FASTA:"
echo "$proteins"

echo
echo "Number of proteins:"
grep -c '^>' "$proteins"


# ============================================================
# 5. remove stale empty InterProScan result if present
# ============================================================

iprxml="${outdir}/annotate_misc/iprscan.xml"

if [ -e "$iprxml" ] && [ ! -s "$iprxml" ]; then
    echo
    echo "Removing empty previous InterProScan XML:"
    rm -f "$iprxml"
fi


# ============================================================
# 6. run InterProScan through funannotate
# ============================================================

cd "$workdir"

echo
echo "Running funannotate iprscan..."
echo

funannotate iprscan \
    -i "$outdir" \
    -m local \
    --iprscan_path "$iprscan" \
    -c 4 \
    --debug


# ============================================================
# 7. validate output
# ============================================================

if [ ! -s "$iprxml" ]; then
    echo
    echo "ERROR: InterProScan XML is missing or empty:"
    echo "$iprxml"
    exit 1
fi

echo
echo "InterProScan completed successfully."
echo

ls -lh "$iprxml"

echo
echo "XML beginning:"
head -n 3 "$iprxml"

echo
echo "Done."