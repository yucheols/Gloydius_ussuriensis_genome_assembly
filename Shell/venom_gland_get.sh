#!/bin/bash
#SBATCH --job-name=get_venom_data
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=25:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

### commands start here ###
# always include the PATH export command so it knows where the program you are calling is 
export PATH=$PWD/sratoolkit.3.4.1-alma_linux64/bin:$PATH

# use prefetch tool to download files
prefetch SRR35908235 -O /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/ncbi_seq  # sample 1
prefetch SRR35908238 -O /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/ncbi_seq  # sample 2

# use fasterq-dump tool to convert .sra files to fastq files
fasterq-dump SRR35908235
fasterq-dump SRR35908238