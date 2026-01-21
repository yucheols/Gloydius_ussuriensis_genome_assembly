#!/bin/bash
#SBATCH --job-name=busco_ussuri
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate busco

# handle java memory issues before starting busco
export _JAVA_OPTIONS="-Xms2g -Xmx64g"

# path to assembly
path_to_asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/hifiasm

# output path
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/busco

# path to busco lineage database
path_to_busco=/home/yshin/mendel-nas1/snake_genome_ass/busco_db/lineages/sauropsida_odb12

# run busco
busco -m genome -i ${path_to_asm}/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa \
  -l ${path_to_busco} -o ${out_dir} -f --metaeuk --offline \
  --download_path /home/yshin/mendel-nas1/snake_genome_ass/busco_db