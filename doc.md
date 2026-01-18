# Whole-genome assembly of the Ussuri Pitviper (*Gloydius ussuriensis*)
*Gloydius ussuriensis* PacBio HiFi genome assembly. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation and https://github.com/amandamarkee/actias-luna-genome

1. __*k*-mer analysis of raw reads using jellyfish__
2. __Genome assembly using hifiasm__
3. __BUSCO__
4. __Genome stats with QUAST__

## 2) Genome assembly using hifiasm
Use the following script to submit a hifiasm job to Mendel. The estimated coverage for this sample is very high (~87x) and the FASTQ file is very big. Use the bigmem partition and request sufficient amount of CPUs to assemble this genome.

``` 
#!/bin/bash
#SBATCH --job-name=hifi_ussuri
#SBATCH --nodes=1
#SBATCH --mem=650G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=48
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%j_%x.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%j_%x.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate mytools

# set taxon name as a variable
name="Gloydius_ussuriensis"

# run hifiasm - put results in their own directory named after the species
hifiasm -o ${name}/${name} -t ${SLURM_CPUS_PER_TASK} /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz

``` 

## 3) BUSCO