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

# activate conda env
source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# set paths
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
outdir="${workdir}/G_ussuriensis_funannotate"
iprdir="/home/yshin/mendel-nas1/interproscan/interproscan-5.78-109.0"
iprscan="${iprdir}/interproscan.sh"

# run interproscan
cd "$workdir"

funannotate iprscan \
    -i "$outdir" \
    -m local \
    --iprscan_path "$iprscan" \
    -c 4 \
    --debug