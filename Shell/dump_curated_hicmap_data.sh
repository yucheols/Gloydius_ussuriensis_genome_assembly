#!/bin/bash
#SBATCH --job-name=dump_hicmap_data
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda env before strict unset-variable checking
source ~/.bash_profile
conda activate scaffolding

set -euo pipefail

# paths
hic="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/hic_rebuild/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split_JBAT.hic"
juicer_tools="/home/yshin/mendel-nas1/juicer_tools/juicer_tools_1.22.01.jar"
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/hic_rebuild/hic_data_dump"

mkdir -p "$outdir"

# resolution = 1 Mb
res=1000000

# JBAT .hic files often use "assembly" as the pseudo-chromosome name
chr="assembly"

java -Xmx24g -jar "$juicer_tools" dump observed NONE \
  "$hic" \
  "$chr" "$chr" \
  BP "$res" \
  "$outdir/Gloydius_wholeassembly_${res}bp_raw.tsv"

echo "Done: $outdir/Gloydius_wholeassembly_${res}bp_raw.tsv"