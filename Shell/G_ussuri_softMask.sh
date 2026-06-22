#!/bin/bash
#SBATCH --job-name=earlGrey_ussuri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=180:00:00
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# soft masking the draft assembly with Earl Grey // do this prior to Hi-C incorporation

# activate the conda env
source ~/.bash_profile
conda activate earlgrey

# path to draft assembly == does not contain mitogenome
path_to_asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa

# output path
outpath=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/
mkdir -p $outpath

# run Earl Grey
# -d flag == Create soft-masked genome at the end? (yes/no, Default: no)
earlGrey -g ${path_to_asm} -s Gloydius_ussuriensis -o ${outpath} -d yes -t ${SLURM_CPUS_PER_TASK}