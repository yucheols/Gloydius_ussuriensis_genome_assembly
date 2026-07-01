#!/bin/bash
#SBATCH --job-name=merqury_scaffold
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env 
source ~/.bash_profile
conda activate merqury

# inputs 
asm_dir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"
read_dir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz"

# parameters
K=21
THREADS=${SLURM_CPUS_PER_TASK}

# outputs 
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/merqury/final_assembly_k${K}
PREFIX=AMNH_21010_chromosome_level_k${K}

mkdir -p "$out_dir"
cd "$out_dir"

echo "[INFO] Using K=$K threads=$THREADS"
echo "[INFO] Assembly: $asm_dir"
echo "[INFO] Reads:    $read_dir"
echo "[INFO] Output:   $out_dir/$PREFIX.*"

# sanity checks
command -v meryl >/dev/null 2>&1 || { echo "[ERROR] meryl not found in PATH"; exit 1; }
command -v merqury.sh >/dev/null 2>&1 || { echo "[ERROR] merqury.sh not found in PATH"; exit 1; }

# build read k-mer database
echo "[INFO] Counting k-mers in reads with meryl..."
meryl k=$K threads=$THREADS count "$read_dir" output reads.k${K}.meryl

# run Merqury == computes QV, completeness, spectra-cn, etc.
echo "[INFO] Running Merqury..."
merqury.sh reads.k${K}.meryl "$asm_dir" "$PREFIX" 2>&1 | tee "$PREFIX.merqury.log"

echo "[INFO] Done."
echo "[INFO] Key outputs to look at:"
echo "  - ${PREFIX}.qv            (assembly QV)"
echo "  - ${PREFIX}.completeness  (k-mer completeness)"
echo "  - ${PREFIX}.spectra-cn.*  (copy-number spectra plots/data)"