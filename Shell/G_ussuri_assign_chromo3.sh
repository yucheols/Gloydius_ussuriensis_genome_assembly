#!/bin/bash
#SBATCH --job-name=assign_chromo_3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate the conda env
source ~/.bash_profile
conda activate genome_assembly  # to access minimap2

set -euo pipefail

# set dir
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Vipera_berus"
mkdir -p "$workdir"
cd "$workdir"

# set variables
QUERY="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies/Vipera_berus_GCA_964194415.1.fa"

# sanity checks
echo "Checking files:"
ls -lh "$QUERY" "$REF"
which minimap2

# run minimap2
echo "Starting minimap2: $(date)"
minimap2 -x asm20 -t "${SLURM_CPUS_PER_TASK}" "$REF" "$QUERY" \
  > Gloydius_vs_Vipera_berus.asm20.paf

echo "Finished minimap2: $(date)"
ls -lh Gloydius_vs_Vipera_berus.asm20.paf