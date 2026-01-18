# Whole-genome assembly of the Ussuri Pitviper (*Gloydius ussuriensis*)
*Gloydius ussuriensis* PacBio HiFi genome assembly. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation and https://github.com/amandamarkee/actias-luna-genome

The individual used for this genome assembly is accessioned at the AMNH Herpetology Collections under the voucher number AMNH 21010.

__Workflow__
1. __Raw read QC with *FastQC*__
2. __*k*-mer analysis of raw reads using *jellyfish*__
3. __Draft genome assembly using *hifiasm*__
4. __Genome completeness using BUSCO__
5. __Genome assembly stats with *QUAST*__
6. __Genome annotation__
   - __RNA read QC__
   - __Adapter trimming__ 
   - __Structural annotation__
   - __Functional annotation__
7. __Scaffolding through Hi-C data incorporation__

## 1) Raw read QC with FastQC
Run QC on the raw PacBio HiFi reads using FastQC. This is only meant to be a "sanity check" and not the actual quality assessment because FastQC assumes short Illumina reads as an input.

```sh
#!/bin/bash
#SBATCH --job-name=fastQC_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --time=10:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# path to the fastq file
path_to_seq=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ

# output directory
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FastQC

# run FastQC
fastqc -o ${out_dir} ${path_to_seq}/AMNH_21010_HiFi.fastq.gz
```

## 2) *k*-mer analysis of raw reads using jellyfish
Conduct a *k*-mer count analysis on the raw reads using jellyfish. This can be useful to estimate the genome size, heterozygosity, etc. Use the following script to submit a job to the AMNH Mendel HPC cluster. The "zcat [...]" line first unpacks the HiFi fastq.gz file (without permanently extracting the content, because this file is massive), converts it into the FASTA format (dropping quality scores, etc.), and pipes that into jellyfish. The output .jf file is then fed into the "jellyfish histo" command to produce the .histo file for viewing.

```sh
#!/bin/bash
#SBATCH --job-name=kmer_ussuri
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --time=06:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# make sure the script to fail fast and loudly if something is broken
set -euo pipefail

# path to sequence read
path_to_seq=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ

# specify output directory
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/

# run jellyfish
zcat ${path_to_seq}/AMNH_21010_HiFi.fastq.gz | awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' | jellyfish count -m 21 -s 10G -t ${SLURM_CPUS_PER_TASK} -C /dev/fd/0 -o ${outdir}/Gloydius_ussuriensis_kmer.jf
jellyfish histo ${outdir}/Gloydius_ussuriensis_kmer.jf -t ${SLURM_CPUS_PER_TASK} > ${outdir}/Gloydius_ussuriensis_kmer.histo
```

The output .histo can be fed into GenomeScope 2.0 (http://genomescope.org/genomescope2.0/) to visualize the results, which look something like: 

![alt text](etc/genomescope_result_summary.png)

The results suggest:
  - Estimated haploid genome size of 1.18 Gb
  - Estimated repeat content of 21.3%
  - High homozygosity (~99%)
  - Very low read error rate (~0.15%) 


## 3) Draft genome assembly using hifiasm
Hifiasm (https://github.com/chhylp123/hifiasm) is a fast, haplotype-resolved assembler for PacBio long-read sequencing data. Use the following script to submit a hifiasm job to the Mendel cluster. The estimated coverage for this sample is very high (~87x) and the fastq.gz file of raw reads is very big (503 Gb). Use the bigmem partition and request a sufficient amount of CPUs and walltime to assemble this genome.

```sh
#!/bin/bash
#SBATCH --job-name=hifi_ussuri
#SBATCH --nodes=1
#SBATCH --mem=650G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# set taxon name as a variable
name="Gloydius_ussuriensis"

# run hifiasm - put results in their own directory named after the species
hifiasm -o ${name}/${name} -t ${SLURM_CPUS_PER_TASK} /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz
``` 

## 4) Genome completeness using BUSCO

## 5) Genome assembly stats with *QUAST*

## 6) Genome annotation
   - __RNA read QC__
   Run FastQC on raw, untrimmed reads.

```sh
   #!/bin/bash
   #SBATCH --job-name=rnaQC_ussuri
   #SBATCH --nodes=1
   #SBATCH --mem=200G
   #SBATCH --partition=compute
   #SBATCH --cpus-per-task=24
   #SBATCH --time=10:00:00
   #SBATCH --mail-type=ALL
   #SBATCH --mail-user=yshin@amnh.org
   #SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
   #SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

   # initiate conda and activate the conda environment
   source ~/.bash_profile
   conda activate genome_assembly

   # path to the fastq file
   path_to_seq=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/FASTQ

   # output directory
   out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/FastQC/pretrim

   # run FastQC
   fastqc -o ${out_dir} ${path_to_seq}/AMNH_21010_Ht_1.fastq.gz ${path_to_seq}/AMNH_21010_Ht_2.fastq.gz ${path_to_seq}/AMNH_21010_Ky_1.fastq.gz ${path_to_seq}/AMNH_21010_Ky_2.fastq.gz ${path_to_seq}/AMNH_21010_Lg_1.fastq.gz ${path_to_seq}/AMNH_21010_Lg_2.fastq.gz ${path_to_seq}/AMNH_21010_Lr_1.fastq.gz ${path_to_seq}/AMNH_21010_Lr_2.fastq.gz ${path_to_seq}/AMNH_21010_Me_1.fastq.gz ${path_to_seq}/AMNH_21010_Me_2.fastq.gz ${path_to_seq}/AMNH_21010_Skin_1.fastq.gz ${path_to_seq}/AMNH_21010_Skin_2.fastq.gz
```

   - __Adapter trimming__ 
   Use trimmomatic to trim adapters and then run FastQC on the trimmed reads.


   - __Structural annotation__
   - __Functional annotation__
## 7) Scaffolding through Hi-C data incorporation