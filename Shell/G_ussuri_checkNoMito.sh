#!/bin/bash
#SBATCH --job-name=checkNoMito_ussuri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# use mitohifi to verify the mito contig is absent from the "cleaned" assembly
# if we correctly eliminated mito contig, this script should fail to run
# the mitohifi .sif image was copied from Dani's file directory

# load apptainer to run mitohifi
module load Apptainer/apptainer-1.2.5

# MitoHiFi flags
# -c: assembled contigs in .fasta (e.g. output from hifiasm) 
# -r: PacBio read .fasta
# -f: Close-related Mitogenome is fasta format
# -g: .gb file associated with the mitogenome reference provided with -c
# -o: genetic code // vertebrate = 2 
# -t: n threads

# cd into working directory (optional)
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi

# set paths for inputs (use /mendel-nas1, not /home/yshin/mendel-nas1)
contigs=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa
mito_ref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/NC_026553.1.fa
mito_gb=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/NC_026553.1.gb

# run mitohifi
apptainer exec -B $PWD --bind /mendel-nas1 mitohifi.sif \
  mitohifi.py -c ${contigs} -f ${mito_ref} -g ${mito_gb} -o 2 -t ${SLURM_CPUS_PER_TASK}