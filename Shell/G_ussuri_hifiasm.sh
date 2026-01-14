#!/bin/sh
#SBATCH --job-name=ussuriensis_hifi
#SBATCH --nodes=1
#SBATCH --mem=512G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=12
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate mytools

# set taxon name as a variable
name="Gloydius_ussuriensis"

# run hifiasm - put results in their own directory named after the species
hifiasm -o ${name}/${name} -t 30 /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz
