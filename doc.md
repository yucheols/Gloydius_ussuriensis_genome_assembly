# Whole-genome assembly of the Ussuri mamushi (*Gloydius ussuriensis*)
*Gloydius ussuriensis* PacBio HiFi genome assembly. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation and https://github.com/amandamarkee/actias-luna-genome

1. __*k*-mer analysis of raw reads using jellyfish__
2. __Genome assembly using hifiasm__
3. __BUSCO__
4. __Genome stats with QUAST__

## 2) Genome assembly using hifiasm
Use the following script to submit a hifiasm job to Mendel

``` 
#!/bin/sh
#SBATCH --job-name=yshin_hifiasm_G_ussuri
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/out/yshin_hifiasm_G_ussuri_%j-%x.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/err/yshin_hifiasm_G_ussuri_%j-%x.err

#conda init

source ~/.bash_profile
conda activate mytools

hifiasm -o home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/Shell/outfiles/hifiasm/out/Gloydius_ussuriensis_v1.asm -t /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis/24GUHW001.hifireads.fastq.gz
``` 

## 3) BUSCO