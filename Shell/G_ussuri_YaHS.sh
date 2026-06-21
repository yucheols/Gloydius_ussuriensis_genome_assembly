#!/bin/bash
#SBATCH --job-name=yahs_hic
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# stop if anything fails
set -euo pipefail

# set paths
indir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding"
outdir="${indir}/yahs_out"

REF="${indir}/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
BAM="${indir}/Gloydius_ussuriensis_HiC_to_draft.markdup.bam"

# make yahs outdir
mkdir -p "$outdir"
cd "$outdir"

# run yahs
yahs "$REF" "$BAM" -o Gloydius_ussuriensis_AMNH_21010_yahs