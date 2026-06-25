#!/bin/bash
#SBATCH --job-name=find_telomere
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=100:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate the conda env
source ~/.bash_profile
conda activate tidk

# set output dir
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/find_telomere"
mkdir -p ${outdir}
cd ${outdir}

# input scaffolded genome
GENOME="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"

# run tidk
echo "Running tidk search for vertebrate telomere repeat TTAGGG: $(date)"
tidk search \
  --string TTAGGG \
  --output Gloydius_ussuriensis_tidk \
  --dir ${outdir} \
  "$GENOME"