# Whole-genome assembly of the Ussuri Pitviper (*Gloydius ussuriensis*)
*Gloydius ussuriensis* PacBio HiFi genome assembly. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation and https://github.com/amandamarkee/actias-luna-genome

The individual used for this genome assembly is accessioned at the AMNH Herpetology Collections under the voucher number AMNH 21010.

__Workflow__
1. __Raw read QC with FastQC__
2. __*k*-mer analysis of raw reads using jellyfish__
3. __Draft genome assembly using hifiasm__
4. __Genome completeness using BUSCO__
5. __Genome assembly stats with QUAST__
6. __Genome annotation__
   - __RNA read QC__
   - __Adapter trimming__ 
   - __Transcriptome assembly__
   - __Structural annotation__
   - __Functional annotation__
7. __Scaffolding through Hi-C data incorporation__
8. __Mitogenome assembly__

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
Hifiasm (https://github.com/chhylp123/hifiasm) is a fast, haplotype-resolved assembler for PacBio long-read sequencing data. Use the following script to submit a hifiasm job to the Mendel cluster. The estimated coverage for this sample is very high (~87x) and the fastq.gz file of raw reads is very big (63 GB). Use the bigmem partition and request a sufficient amount of CPUs and walltime to assemble this genome.

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

# set output directory
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/hifiasm

# set taxon name as a variable
name="Gloydius_ussuriensis"

# run hifiasm - put results in their own directory named after the species
hifiasm -o ${out_dir}/${name}_v1.asm -t ${SLURM_CPUS_PER_TASK} /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz
``` 

## 4) Genome completeness using BUSCO

## 5) Genome assembly stats with QUAST

## 6) Genome annotation
   - __RNA read QC:__
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

   - __Adapter trimming:__ 
   Use trimmomatic to trim adapters and then run FastQC on the trimmed reads.

   - __Transcriptome assembly:__

   - __Structural annotation:__
   - __Functional annotation:__

## 7) Scaffolding through Hi-C data incorporation

## 8) Mitogenome assembly
Because PacBio HiFi reads are long and highly contiguous, it is possible to assemble a full mitogenome as a bycatch. This can be done easily using some existing tools and a reference mitogenome to "fish out" the mitochondrial contigs from HiFi reads. There is already a conspecific mitogenome reference available on GenBank (NC_026553.1). We can fetch this mitogenome like so:

```txt
# create working directory for mitogenome assembly
mkdir -p /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref

# install entrez to fetch seq from ncbi
conda activate genome_assembly
conda install -y -c bioconda entrez-direct

# fetch seq and annotation
ref_acc=NC_026553.1
echo ${ref_acc}

efetch -db nucleotide -id ${ref_acc} -format fasta > ${ref_acc}.fa
efetch -db nucleotide -id ${ref_acc} -format gbwithparts > ${ref_acc}.gb
```

Then, create a conda environment for mitogenome assembly and install some packages that will be used.

```txt
# create conda environment for mitogenome assembly
conda create -n mito_assembly
conda activate mito_assembly

# install some packages
conda install -y -c bioconda minimap2 samtools flye
```

Once everything is prepared, use the pipeline below to assemble a mitogenome. 

```sh
#!/bin/bash
#SBATCH --job-name=mito_ussuri
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda environment
source ~/.bash_profile
conda activate mito_assembly

##############
#  1. setup  #
##############

# make the run fail loudly if something is broken
set -euo pipefail

# set some paths as variables
path_to_hifi=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ
path_to_mitoref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_out

# set mapping quality score
# set this to 60 to drop numts and other crap that we do not want
mq=60


################
#  2. mapping  #
################

# print this message when starting
echo "start mapping reads to reference..."

# run minimap2
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${path_to_mitoref}/NC_026553.1.fa \
  ${path_to_hifi}/AMNH_21010_HiFi.fastq.gz | \
  samtools view -@ 8 -b -F 4 | \
  samtools sort -@ 8 -T ${out_dir}/tmp.mito -o ${out_dir}/mito_map.bam
samtools index ${out_dir}/mito_map.bam

#############################
#  3. extract mapped reads  #
#############################

# print this message when starting
echo "extract mapped reads for mapping quality score >= ${mq}"

# use samtools
samtools view -h ${out_dir}/mito_map.bam | \
  awk -v mq=${mq} 'BEGIN{OFS="\t"} $1 ~ /^@/ || $5 >= mq' | \
  samtools view -b -F 4 -o ${out_dir}/mito_map.MQ.bam

samtools index ${out_dir}/mito_map.MQ.bam

# Convert to FASTQ for assembly
samtools fastq ${out_dir}/mito_map.MQ.bam | gzip > ${out_dir}/mito_reads.MQ.fastq.gz

# print out this message
echo "mito reads written as FASTQ: mito_reads.MQ.fastq.gz"


#######################
#  4. assemble reads  #
#######################

# print this message when starting
echo "assemble reads with flye..."

# run flye
flye --pacbio-hifi ${out_dir}/mito_reads.MQ.fastq.gz \
  --out-dir ${out_dir}/flye_mito -t ${SLURM_CPUS_PER_TASK} --genome-size 20k

mito_asm=${out_dir}/flye_mito/assembly.fasta

# print out this message
echo "mitochondrial genome assembly complete! assembly written to: ${mito_asm}"


###############################################################
#  5. map mito reads back to mito assembly and compute depth  #
###############################################################

# print this message when starting
echo "remap mito reads back to assembled contigs..."

# run minimap2
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${mito_asm} ${out_dir}/mito_reads.MQ.fastq.gz | \
  samtools sort -@ 8 -o ${out_dir}/mito_realign.bam
samtools index ${out_dir}/mito_realign.bam

# compute depth
echo "compute depth..."

# run samtools depth
samtools depth ${out_dir}/mito_realign.bam > ${out_dir}/mito_depth.tsv
awk '{sum+=$3; n++} END{print "mean_depth\t"sum/n"\npositions\t"n}' ${out_dir}/mito_depth.tsv > ${out_dir}/depth_summary.txt


###########################################
#  6. validate by comparing to reference  #
###########################################

# print this message when starting
echo "validate assembly by comparing to reference..."

# run minimap2
minimap2 -x asm5 ${path_to_mitoref}/NC_026553.1.fa ${mito_asm} > ${out_dir}/mito_vs_ref.paf


###################################
# print this when all is complete #
###################################
echo "mitochondrial genome assembly pipeline complete!"

```

Once this is done, run a quick sanity check:

```txt
# run this in the directory containing the assembly.fasta file (which is the "flye_mito" folder)
seqkit stats assembly.fasta
grep -c "^>" assembly.fasta
```