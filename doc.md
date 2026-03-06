# Whole-genome assembly of the Ussuri Pitviper (*Gloydius ussuriensis*)
Reference genome assembly of *Gloydius ussuriensis* using PacBio HiFi long-read sequencing, RNA-seq, and Hi-C. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation, https://github.com/amandamarkee/onigra_genome, and https://github.com/amandamarkee/actias-luna-genome. Also, thanks to Amanda Markee (AMNH IZ), Daniela Garcia (AMNH Herp), Jon Hoffman (AMNH Herp), Dylan DeBaun (AMNH Herp), Dean Bobo (AMNH ICG), and Sajesh Singh (AMNH CS) for help and discussions.

The genome sequencing was done on the PacBio Revio system (1 SMRT cell) and RNA sequencing was done on the Illumina NovaSeq X (151bp PE). The individual used for this genome assembly is accessioned at the AMNH Herpetology Collections under the field number AMNH 21010.

__Note:__ I created this documentation in the hopes that my friends and future RGGS students doing genomics may find it useful. If you have any questions about any parts of the assembly process outlined below, please don't hesitate to email me at yshin@amnh.org or post a question in the "issues" section of the repository.

### Workflow

1. __[A quick sanity check on the dataset](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#1-a-quick-sanity-check-on-the-dataset)__
2. __[*k*-mer analysis of raw reads using jellyfish](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#2-k-mer-analysis-of-raw-reads-using-jellyfish)__
3. __[Draft genome assembly using hifiasm](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#3-draft-genome-assembly-using-hifiasm)__
4. __[Contamination screening](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#4-contamination-screening)__
   - __Screening for potential non-vertebrate contaminants using blobtools__
   - __Identifying and removing mitochondrial contigs from the draft assembly__
   - __Genomewide mean sequencing coverage after contamination screening__
5. __[Draft assembly evaluation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#5-draft-assembly-evaluation)__ 
   - __Genome completeness using BUSCO__
   - __Genome assembly stats with QUAST__
   - __*k*-mer based assembly evaluation with Merqury__
6. __[Scaffolding through Hi-C data incorporation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#6-scaffolding-through-hi-c-data-incorporation)__
   - __Setup__
7. __[Genome annotation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#7-genome-annotation)__
   - __Setup__
   - __RNA read QC (pre-trimming)__
   - __Adapter trimming & post-trimming QC__ 
   - __(Pre-Hi-C) RNA alignment to draft using HiSat2__
   - __(Pre-Hi-C) Draft-guided transcriptome assembly using StringTie__
   - __(Post-Hi-C) Repeat masking (soft masking) using Earl Grey__
   - __(Post-Hi-C) Re-run HiSat2 and StringTie on the scaffolded and masked genome__
   - __(Post-Hi-C) Structural annotation__
   - __(Post-Hi-C) Functional annotation__
8. __[Mitogenome assembly](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/doc.md#8-mitogenome-assembly)__
    - __"Manual" annotation with MITOS2__
    - __Submitting mitogenome to GenBank__

## 1) A quick sanity check on the dataset
Even before doing anything, let's do a very quick sanity check on the HiFi data to check the reads we have are actually from our target species. Let's take a chunk from the HiFi FASTQ file, after cd'ing into the directory containing the fastq.gz file:
```txt
zcat AMNH_21010_HiFi.fastq.gz | head -n 2
```

When you run the line above you will get a string of sequence read printed out:
![alt text](etc/seqchunk.PNG)

Let's copy the whole chunk and paste it into NCBI BLAST. The result  will look something like this:
![alt text](etc/blast.PNG)

We can see that the top hits are from the Tiger Rattlesnake (*Crotalus tigris*). This is the expected result, although the species itself is not in the genus *Gloydius*. This is because our target species (*G. ussuriensis*) only have isolated, Sanger-sequenced mitochondrial and nuclear markers accessioned in GenBank. Since *C. tigris* is in the same subfamily as *Gloydius*, this means that our FASTQ file contains actual genome sequence reads from *G. ussuriensis*. 

Let's do the same for the RNA seq reads:
```txt
zcat AMNH_21010_Ht_1.fastq.gz | head -n 2
```
Repeat this for each tissue type and read, and you will see that all the hits come out to be snake mRNA genes. For the skin RNA reads, using the regular megablast option to optimize for highly similar sequences will print out a warning "No significant similarity found." However, switching to optimization for somewhat similar sequences (blastn) will print out hits for snake genes (from *Thamnophis* and *Candoia*).

Also, let's estimate the mean sequencing depth ("coverage") we got from the raw gzipped FASTQ file. We can get the coverage based on the number of bp in our FASTQ file divided by the estimated genome size. We don't know the genome size of *G. ussuriensis*, but we can use the genome size of related species as a proxy - here we will use the Prairie Rattlesnake (*Crotalus viridis*) reference genome for this purpose (1.3 Gb). 

```txt
# get the number of bp
zcat AMNH_21010_HiFi.fastq.gz | awk 'NR%4==2{bp+=length($0)} END{print bp/1e9 " Gb"}'
```
The output is 134.449 Gb. If we use the *C. viridis* ref genome size, the coverage would be 134.449/1.3 = 103x.

Next, let's check the average read length. We can do so by running the "seqkit stats" command on our fastq.gz file:
```txt
# cd into the directory containing the .fastq.gz file
# cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ

# activate the "genome_assembly" conda env to access seqkit / this conda env was already in place
# with seqkit installed
conda activate genome_assembly
seqkit stats AMNH_21010_HiFi.fastq.gz
```
![alt text](etc/raw_read_avg_seqlen.PNG)

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

The output .histo file can be fed into GenomeScope 2.0 (http://genomescope.org/genomescope2.0/) to visualize the results, which look something like: 

![alt text](etc/genomescope_result_summary.png)

The results suggest:
  - Estimated haploid genome size of 1.18 Gb
  - Estimated repeat content of 21.3%
  - High homozygosity (~99%)
  - Very low read error rate (~0.15%) 
  - Very high sequencing coverage/depth (~100x)

## 3) Draft genome assembly using hifiasm
Hifiasm (https://github.com/chhylp123/hifiasm) is a fast, haplotype-resolved assembler for PacBio long-read sequencing data. Use the following script to submit a hifiasm job to the Mendel cluster. The estimated coverage for this sample is very high (100x ~ 114x) and the fastq.gz file of raw reads is very big (63 GB). Use the bigmem partition and request a sufficient amount of CPUs and walltime to assemble this genome.

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
This will take approximately 22 hours to run on Mendel. The output files should look something like this:
![alt text](etc/hifiasm.PNG)

Among all these files the "bp.p_ctg.gfa" file contains the assembly graph of primary contigs and this is the file we will use in all downstream stuff. However, some packages take .fa as input, not .gfa. So we first need to convert the .gfa file of primary contigs into .fa file. Run the line below in the directory containing the hifiasm output files:

```txt
awk '$1=="S"{print ">"$2"\n"$3}' Gloydius_ussuriensis_v1.asm.bp.p_ctg.gfa > Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa
```

## 4) Contamination screening
The draft assembly likely contains mitochondrial contigs and/or potential microbial contaminants. So, it is always a good idea to check for these and "clean up" the genome before finalizing the assembly and publishing it. This section is based on https://github.com/amandamarkee/actias-luna-genome

Let's setup the workspace for this clean up step.
```
# make a directory for clean up work
mkdir genome_cleanup

# copy the draft assembly from its directory to the clean up workspace
cp Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# create a conda env for blobtools
conda create -n blobtools
conda activate blobtools

# install blobtools
conda install bioconda::blobtools

```

Let's start by getting the "pre-cleanup" stats
```
seqkit stats Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa > preclean_stats.txt
cat preclean_stats.txt
```
We can see that there are 140 contigs total.

----------------------------------------------------------------------------------------------------
### Screening for potential non-vertebrate contaminants using blobtools
We will use blobtools to identify potential non-vertebrate (e.g., microbial) contaminants. To run blobtools, we need the following:
  - nodes.dmp and names.dmp files from NCBI taxdump
  - .bam file output from minimap2 and samtools
  - .nt blast hit file from ncbi megablast run

Let's prep these files. First, get nodes.dmp and names.dmp files from NCBI taxdump.
```
# activate the blobtools conda env to use blobtools
conda activate blobtools

# create a directory to store data
mkdir data

# download NCBI taxdump and create nodes.dmp and names.dmp
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp
```
Then, if you do "ls data/" you should see something like:
![alt text](etc/lsdata.PNG)


To generate the blast hit file, we will use the Huxley cluster because it already has the database soft linked to /nas3. Blast is also available as a module on Huxley. But before running blast, there is one thing to consider - blasting a whole genome assembly againt the database will take a long time. The job will either fail due to insufficient computational resources or walltime. 

One loophole is to divide up the assembly fasta file and then running an array job (thanks Amanda for this tip!).

To divide up the assembly fasta file, navigate to the directory containing the assembled genome and run these commands in the console:

```txt
# make directory for split fasta files
mkdir asm_div

# spilit fasta
awk 'BEGIN {n_seq=0;} /^>/ \
{if(n_seq%50==0){file=sprintf("seq_%d.fa",n_seq);} \
print >> file; n_seq++; next;} { print >> file; }' < Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa

# move split fastas to asm_div directory
mv seq_* asm_div/

# put split files into a list
cd asm_div/
ls *.fa > blast.list

# check the number of lines in the blast.list file to use as an input to the array command
wc -l blast.list
```

Then run the array job as below:
```sh
#!/bin/bash
#SBATCH --job-name=blobBlst_ussuri
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --array=1-3
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/nas4/G_ussuriensis_Chromo/slurm_logs/slurm-%x_%j.out
#SBATCH --error=/home/yshin/nas4/G_ussuriensis_Chromo/slurm_logs/slurm-%x_%j.err

# this is run on the AMNH Huxley cluster
# load module
module load BLAST/2.11.0-Linux_x86_64 

# path to genome cleanup directory
wdir=/home/yshin/nas4/G_ussuriensis_Chromo/genome_cleanup

# path to DB
export BLASTDB=/nas3/database/nt_2024_10_23
blastdbcmd -db nt -info

# path to assembly
path_to_asm=${wdir}/asm_div

# output dir
out_dir=/home/yshin/nas4/G_ussuriensis_Chromo/genome_cleanup/blast_out

# assign and print array ID
i=${SLURM_ARRAY_TASK_ID}
P=$(awk "NR==$i" ${path_to_asm}/blast.list)

# run blast
blastn \
  -db nt \
  -task megablast \
  -query ${path_to_asm}/${P} \
  -out ${out_dir}/${P}.blast.out \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -evalue 1e-20 \
  -outfmt "6 qseqid staxids bitscore sseqid sskingdoms sscinames" \
  -num_threads=${SLURM_CPUS_PER_TASK}  
```

Once this job finishes running, navigate to the output directory and merge the three output files into a single .nt file to be used as an input for blobtools.
```
# the output dir is /home/yshin/nas4/G_ussuriensis_Chromo/genome_cleanup/blast_out
cat *.blast.out > Gloydius_megablast.nt
```
Now, send this file to Mendel so that it can be used as input for blobtools.

We also need mapping files as a final input. These can be made by mapping the raw fastq file back to the reference (i.e. draft assembly) using minimap2, and then converting the output SAM file into a BAM file and sorting and indexing them using samtools. Run the job script below on Mendel:

```sh
#!/bin/bash
#SBATCH --job-name=blobMapping_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate genome_assembly

# designate paths to assembly and fastq
path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# output directory
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/mapping

# mapping
# -@ 8 means "use 8 threads" 
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${path}/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa \
  ${path}/AMNH_21010_HiFi.fastq.gz | samtools view -@ 8 -b \
  | samtools sort -@ 8 -o ${outdir}/ussuri_aln_sorted.bam

# index BAM 
samtools index ${outdir}/ussuri_aln_sorted.bam 
```
Now we have all the input files required to run blobtools. Run blobtools with this script on Mendel:

```sh
#!/bin/bash
#SBATCH --job-name=blobView_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate blobtools

# set paths to different inputs
asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa
aln_sorted_bam=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/mapping/ussuri_aln_sorted.bam
megablast=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/blast_out/Gloydius_megablast.nt
nodes=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/data/nodes.dmp
names=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/data/names.dmp

# output dir
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/blob_out

# run blobtools
blobtools create -i ${asm} -b ${aln_sorted_bam} -t ${megablast} --nodes ${nodes} --names ${names} \
  -o ${outdir}/ussuri_blob_results
```

Now, run the lines below in the head node to generate output files and plots.
```
# activate blobtools conda env
conda activate blobtools

# cd into the blobtools output directory
# view blobtools results and plot
blobtools view -i ussuri_blob_results.blobDB.json
blobtools plot -i ussuri_blob_results.blobDB.json
```
This will generate a whole bunch of files:
![alt text](etc/blobtools.PNG)

Now, scp these files from the cluster to a local device for viewing.
![alt text](etc/ussuri_blob_results1.png)
![alt text](etc/ussuri_blob_results2.png)

We can see there are no obvious non-vertebrate contaminants.

----------------------------------------------------------------------------------------------------
### Identifying and removing mitochondrial contigs from the draft assembly
PacBio long-read assemblies usually contain complete mitogenomes as a sequencing bycatch, and it is a good idea to identify and remove them from primary assemblies in order to get a clean nuclear whole genome. This can be done by blasting the mitochondrial reference to the draft assembly. *G. ussuriensis* already has a couple assembled mitogenomes on GenBank. We will use one of them (NC_026553.1) as a mitocohndrial reference. Let's fetch this sequence from GenBank using entrez.

```txt
# create a directory to store mito reference
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

# move the .fa file to the "genome_cleanup" directory
cp NC_026553.1.fa /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup
```
Once you get this file, we can run blast on Mendel with the script below. Note that we are creating a blast database internally in this script. The "-outfmt 6" flag means that the output format will be tabular. 
```sh
#!/bin/bash
#SBATCH --job-name=idMito_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=32
#SBATCH --time=144:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# load blast
module load NCBI/blast-2.10.1+

# set paths to input
mito_ref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/NC_026553.1.fa
asm_fa=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa

# set output dir
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# build DB in a dedicated folder with a proper prefix name
dbdir=${outdir}/blast_db_mito
mkdir -p "${dbdir}"
dbprefix="${dbdir}/ussuri_asm"

echo "Making BLAST DB:"
echo "  asm_fa   = ${asm_fa}"
echo "  dbprefix = ${dbprefix}"

makeblastdb -in "${asm_fa}" -dbtype nucl -out "${dbprefix}"

# sanity check DB
echo "DB files:"
ls -lh "${dbprefix}".n* || { echo "ERROR: DB files not found"; exit 1; }
blastdbcmd -db "${dbprefix}" -info || { echo "ERROR: blastdbcmd failed"; exit 1; }

# run blast
blastn \
-query "${mito_ref}" \
-db "${dbprefix}" \
-out ${outdir}/mito_contigs_blast.out \
-outfmt 6 \
-evalue 1e-20 \
-num_threads ${SLURM_CPUS_PER_TASK}

# convert output into csv
tr '\t' ',' < ${outdir}/mito_contigs_blast.out > ${outdir}/mito_contigs_blast.csv
```

The output should look something like this:
![alt text](etc/mito_blast_out.PNG)

Let's break this down. 
- Column 1: "qseqid" (means "query sequence ID") - Shows the GenBank accession number of the mitogenome reference. 
- Column 2: "sseqid" (means "subject sequence ID") - Shows the contig names in the draft assembly. 
- Column 3: "pident" (percent identity) - Shows the percent match between the mitochondrial reference and contig. 
- Column 4: length - Shows the length of aligned region. 
- Column 5: mismatch - Shows the number of mismatched base pairs.
- Column 6: "gapopen" - Number of gaps in the alignment.
- Column 7: "qstart" - Start position on query (mitochondrial reference).
- Column 8: "qend" - End position on query.
- Column 9: "sstart" - Start position on subject (draft assembly). 
- Column 10: "send" - End position on subject.
- Column 11: "evalue" - evalue of 0 means that the probability the match is random is 0.
- Column 12: "bitscore" - Alignment score

Based on this, let's focus on the contig ptg000073c. This contig has a percent identity of 99% to the mitochondrial reference, the length of aligned region is 17,212bp, and mismatch is relatively small at 161bp over 17kb. The mitochondrial reference is 17,208bp. Based on this, this contig very likely contains the complete mitogenome, which is expected from PacBio long-read sequencing results. 

Let's look at this contig in more detail. First, let's check the length:
```
# make sure that you are working in the "genome_cleanup" directory
# activate conda env
conda activate genome_assembly

# check contig length
seqkit fx2tab -n -l Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa | grep -w ptg000073c
```
This will show that the contig is 68.8Kb long, which is much longer than the actual mitogenome. But this likely because one whole mitogenome was repeated several time on the same contig. For example:
![alt text](etc/mito_blast_out_head.PNG)
In the first three lines we can see that 17,212bp stretch is repeated three times, and the fourth line shows another 15,348bp region, and there are other, shorter hits. Combined, this is about 67Kb long. Also, 68.8Kb is about 17.2Kb x 4. 

We can verify this by blasting the putative mitochondrial contig to itself.
First, let's store this contig as a separate file.
```
# store mito contigs into a separate file
seqkit grep -p ptg000073c Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa \
-o Gloydius_ussuriensis_AMNH_21010_mito.fa
```

Then, use the script below to run blast:
```sh
#!/bin/bash
#SBATCH --job-name=selfalignMito_ussuri
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# load modules
module load NCBI/blast-2.10.1+

# activate conda env
source ~/.bash_profile
conda activate genome_assembly

# inputs
mito_fa=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/Gloydius_ussuriensis_AMNH_21010_mito.fa
mito_contig="ptg000073c"

# output dir
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/selfalign_mito
mkdir -p "${outdir}"

# build BLAST DB from mito contig
dbdir="${outdir}/blast_db"
mkdir -p "${dbdir}"
dbprefix="${dbdir}/${mito_contig}"

echo "Making BLAST DB:"
echo "  mito_fa   = ${mito_fa}"
echo "  dbprefix  = ${dbprefix}"
makeblastdb -in "${mito_fa}" -dbtype nucl -out "${dbprefix}"

# sanity check DB
echo "DB info:"
blastdbcmd -db "${dbprefix}" -info || { echo "ERROR: blastdbcmd failed"; exit 1; }

# self-BLAST (mito vs mito)
# use stricter settings + filter to highlight tandem repeat structure
# include qlen/slen to help interpret offsets and copy structure.
tsv="${outdir}/${mito_contig}_selfblast.tsv"
csv="${outdir}/${mito_contig}_selfblast.csv"

blastn \
  -query "${mito_fa}" \
  -db "${dbprefix}" \
  -out "${tsv}" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -evalue 1e-50 \
  -num_threads "${SLURM_CPUS_PER_TASK}" \
  -dust no \
  -soft_masking false

# remove the trivial full-length self hit (diagonal) to make repeats obvious
# (keeps alignments where query and subject coordinates differ)
filtered_tsv="${outdir}/${mito_contig}_selfblast_filtered.tsv"
awk 'BEGIN{OFS="\t"} !($4==$13 && $7==1 && $8==$13 && $9==1 && $10==$13) {print}' "${tsv}" > "${filtered_tsv}"

# 5) convert to CSV (both full and filtered)
tr '\t' ',' < "${tsv}" > "${csv}"
tr '\t' ',' < "${filtered_tsv}" > "${outdir}/${mito_contig}_selfblast_filtered.csv"

echo "Done."
echo "Outputs:"
echo "  Full TSV      : ${tsv}"
echo "  Filtered TSV  : ${filtered_tsv}"
echo "  Full CSV      : ${csv}"
echo "  Filtered CSV  : ${outdir}/${mito_contig}_selfblast_filtered.csv"
```

Let's examine the output:
```
cat ptg000073c_selfblast_filtered.tsv
```
![alt text](etc/selfalign_mito.PNG)

We can see that the alignment length here is 17,211bp. 17,211 x 4 is exactly 68,844bp. In the results, we can also see alignment lengths of 51,633bp (17,211 x 3) and 34,422bp (17,211 x 2). This basically means that the contig ptg000073c contains four tandem copies of the complete mitogenome.

Now that we identified the mitochondrial contig, let's snip it out from our draft assembly. This "no mito" assembly will be used in all downstream assembly steps.
```
# run on the head node
# set paths to input
asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa
mito_fa=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/ptg000073c.fa

# set output dir
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito

# read the contig ID from the mito file
mito_id=$(grep -m1 "^>" "$mito_fa" | sed 's/^>//' | awk '{print $1}')
echo "mito contig ID = $mito_id"

# create the no-mito assembly
seqkit grep -v -p "$mito_id" "$asm" > "$outdir/Gloydius_ussuriensis_AMNH_21010_noMito.fa"

# sanity check == contig count should drop by 1
echo "original contigs:"; grep -c "^>" "$asm"
echo "no-mito contigs:"; grep -c "^>" "$outdir/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
```
We can see that the total number of contigs went from 140 to 139 after removing the mito contig.

Also, let's store a single mitogenome copy for later use:
```
# run this in the "genome_cleanup" directory
seqkit grep -p ptg000073c Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa > ptg000073c.fa
seqkit subseq -r 1:17211 ptg000073c.fa > mito_singlecopy.fa
seqkit stats mito_singlecopy.fa
```
![alt text](etc/mito_singlecopy_seqkit.PNG)

This is exactly what we expect to see.

As one final step, let's verify that the mitochondrial contig is gone from our draft assembly. This can be done using MitoHiFi (https://github.com/marcelauliano/MitoHiFi), which is a pipeline that finds, circularises and annotates mitogenome from PacBio assemblies. Installing this can be quite tricky. So, instead of fighting my way through installation, I created a directory for mitohifi under the "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo" directory just copied the mitohifi.sif file from Dani's directory on Mendel. Also, it is possible to run MitoHiFi through the Galaxy web server (https://galaxy-main.usegalaxy.org/).
```
# make mitohifi directory
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo
mkdir -p mitohifi
cd mitohifi/

# copy mitohifi.sif file over to this directory
cp /home/dgarcia/mendel-nas1/PacBio/Helicops_angulatus_Aug2024/mitogenome_assembly/MitoHiFi/mitohifi.sif .
``` 
After this, run the script below. If we correctly eliminated the mitochondrial contig, __this run should fail because there are no mitochondrial reads to find and annotate__

```sh
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

# make output dir
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/check_no_mito
mkdir -p ${outdir}

# copy mitohifi.sif into outdir
cp /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/mitohifi.sif ${outdir}

# cd into working directory
cd ${outdir}

# set paths for inputs (use /mendel-nas1, not /home/yshin/mendel-nas1)
contigs=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa
mito_ref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref/NC_026553.1.fa
mito_gb=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/mito_ref/NC_026553.1.gb

# in-container file visibility test
apptainer exec --bind "/home/yshin:/home/yshin" --bind "/mendel-nas1:/mendel-nas1" --pwd "$PWD" mitohifi.sif \
  bash -lc 'ls -lh /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/NC_026553.1.fa'

# run mitohifi
apptainer exec --bind "/home/yshin:/home/yshin" --bind "/mendel-nas1:/mendel-nas1" --pwd "$PWD" \
  mitohifi.sif python3 /opt/MitoHiFi/src/mitohifi.py -c ${contigs} -f ${mito_ref} -g ${mito_gb} -o 2 -t ${SLURM_CPUS_PER_TASK}
``` 
The error log should look something like this:
![alt text](etc/check_nomito.PNG)

In contrast, running MitoHiFi with the assembly still containing the mitochondrial contig will produce the full result, including the annotated mitogenome.
```sh
#!/bin/bash
#SBATCH --job-name=mitoHiFi_ussuri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# use mitohifi to identify and assemble a mitogenome from pacbio assembly
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
contigs=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa
mito_ref=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/NC_026553.1.fa
mito_gb=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/NC_026553.1.gb

# in-container file visibility test
apptainer exec --bind "/home/yshin:/home/yshin" --bind "/mendel-nas1:/mendel-nas1" --pwd "$PWD" mitohifi.sif \
  bash -lc 'ls -lh /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi/NC_026553.1.fa'

# run mitohifi
apptainer exec --bind "/home/yshin:/home/yshin" --bind "/mendel-nas1:/mendel-nas1" --pwd "$PWD" \
  mitohifi.sif python3 /opt/MitoHiFi/src/mitohifi.py -c ${contigs} -f ${mito_ref} -g ${mito_gb} -o 2 -t ${SLURM_CPUS_PER_TASK}
```
![alt text](etc/mitohifi.PNG)

Let's scp this file into a local device and have a closer look.
```
# run this locally in "/home/yshin/Gloydius_ussuriensis_genome_assembly/outfiles" directory
mkdir mitohifi
scp -r yshin@mendel.sdmz.amnh.org:/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mitohifi .
```

----------------------------------------------------------------------------------------------------
### Genomewide mean sequencing coverage after contamination screening
The coverage we calculated above from the raw fastq file represents the coverage based on total sequencing yield (i.e., amount of bases sequenced), prior to contamination screening. This value is going to be somewhat different from the mean sequencing coverage based on the reads mapped back to the assembly cleaned of contaminants and mitochondrial contigs. To calculate the mean sequencing coverage, let's first minimap2 to map the reads to the "no mito" assembly to generate and index a bam file. We will then calculate the coverage from this bam file:
```sh
#!/bin/bash
#SBATCH --job-name=calcCov_ussuri
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env for genome assembly
source ~/.bash_profile
conda activate genome_assembly

# set paths
no_mito_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa
reads_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/AMNH_21010_HiFi.fastq.gz   
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/cov
mkdir -p ${out_dir}

# index fasta
samtools faidx ${no_mito_dir}

# map HiFi reads
echo "start mapping....."
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${no_mito_dir} ${reads_dir} \
| samtools sort -@ 8 -o ${out_dir}/no_mito_hifi.bam

samtools index ${out_dir}/no_mito_hifi.bam

# activate a separate conda env for the newer samtools version
source ~/.bash_profile
conda activate samtools

echo "activated a separate conda env for the newer samtools version....."

# calculate mean genome wide coverage 
echo "calculate mean genome wide coverage....."
samtools coverage ${out_dir}/no_mito_hifi.bam \
| awk 'NR>1 {sum += $3*$7; len += $3} END {print "mean depth of coverage =", sum/len "x"}'

# calculate how much of the genome (%) is covered >= 20x, 30x, and 50x
echo "calculate how much of the genome (%) is covered >= 20x, 30x, and 50x....."
samtools depth -a ${out_dir}/no_mito_hifi.bam \
| awk '{
  total++;
  if($3>=20) c20++;
  if($3>=30) c30++;
  if($3>=50) c50++;
} END {
  print "≥20× =", 100*c20/total "%";
  print "≥30× =", 100*c30/total "%";
  print "≥50× =", 100*c50/total "%";
}'

echo "done"
```
The result shows:
- 83.3x genome-wide mean coverage 
- 99% of the genome has ≥20× coverage
- 97.4% of the genome has ≥30× coverage
- 85.2% of the genome has ≥50× coverage


## 5) Draft assembly evaluation 
### Genome completeness using BUSCO
Now let's assess the completeness of our draft assembly output from hifiasm. BUSCO (Benchmarking Universal Single-Copy Orthologs) is a common metric to assess genome completeness. It uses a lineage-specific dataset to search for the presence/absence of highly conserved genes for that lineage in your genome assembly. We will use *compleasm* (https://github.com/huangnengCSU/compleasm) to assess genome completeness. This provides a faster alternative to the regular BUSCO package for large genome assemblies.

First, create a conda environment for *compleasm* and install the package:
```txt
# create a conda environment for compleasm and install it
conda create -n compleasm -c conda-forge -c bioconda compleasm
conda activate compleasm
compleasm -h
```

Run compleasm on Mendel with this script:
```sh
#!/bin/bash
#SBATCH --job-name=compleasm_ussuri
#SBATCH --nodes=1
#SBATCH --mem=200G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate compleasm

# set paths as variables
path_to_asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito
out_path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/compleasm

# run compleasm
compleasm run -a ${path_to_asm}/Gloydius_ussuriensis_AMNH_21010_noMito.fa -o ${out_path} \
  -t ${SLURM_CPUS_PER_TASK} -l sauropsida --odb odb12
```

The script above should take only about 35 minutes to run. The results are:
![alt text](etc/busco_compleasm.PNG)

Compleasm output can be visualized in R:
```r
######  visualize compleasem output

# clean working environment
rm(list = ls(all.names = T))
gc()

# load packages
library(ggplot2)


######   compleasm result categories

# S (Single Copy Complete Genes): The BUSCO genes that can be entirely aligned in the assembly, with only one copy present.
# D (Duplicated Complete Genes): The BUSCO genes that can be completely aligned in the assembly, with more than one copy present.
# F (Fragmented Genes, subclass 1): The BUSCO genes which only a portion of the gene is present in the assembly, and the rest of the gene cannot be aligned.
# I (Fragmented Genes, subclass 2): The BUSCO genes in which a section of the gene aligns to one position in the assembly,
#                                   while the remaining part aligns to another position.
# M (Missing Genes): The BUSCO genes with no alignment present in the assembly.


######   compleasm result dataframe
compl_res <- data.frame(category = c('S', 'M', 'D', 'F'),
                        count = c(5764, 285, 50, 19))

print(compl_res)

######   calculate proportion
compl_res$percentage <- round((compl_res$count / sum(compl_res$count)) * 100, digits = 2) 
print(compl_res)

######   plotting order
compl_res$category <- factor(compl_res$category, levels = c('S', 'D', 'F', 'M'))


######   plot
ggplot(compl_res, aes(x = 2, y = percentage, fill = category)) +
  geom_col(width = 1, color = 'white') +
  xlim(0.1, 3.5) +
  coord_polar(theta = 'y') +
  scale_fill_manual(values = c(S = '#B3E2CD',
                               D = '#FDCDAC',
                               F = "#FFF2AE",
                               M = '#CBD5E8')) +
  theme_void() +
  theme(legend.position = 'none')

######  export plot
ggsave('Rplots/compleasm_result.png', width = 21, height = 20, dpi = 800, units = 'cm')
```
----------------------------------------------------------------------------------------------------
### Genome assembly stats with QUAST
QUAST is a quality assessment tool for genome assemblies. Installing QUAST in the "genome_assembly" conda environment is not possible because of python version clashes - we need to create a separate conda environment for this package.

Install QUAST:
```txt
conda create -n quast
conda activate quast
conda install bioconda::quast
``` 

When this is done, run QUAST with the script below:

```sh
#!/bin/bash
#SBATCH --job-name=quast_ussuri
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the quast conda environment
source ~/.bash_profile
conda activate quast

# path to assembly
path_to_asm=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito

# output directory
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/quast
mkdir -p ${out_dir}

# run quast
quast.py -t ${SLURM_CPUS_PER_TASK} ${path_to_asm}/Gloydius_ussuriensis_AMNH_21010_noMito.fa -o ${out_dir} 
```
----------------------------------------------------------------------------------------------------
### *k*-mer based assembly evaluation with Merqury
Let's use Merqury as another assembly evaluation tool. Similar to jellyfish, this one is also based on *k*-mers. But unlike jellyfish we ran above (which was on raw read fastq), we are running merqury on the assembly. Basically, jellyfish just counts *k*-mers from on your reads. On the other hand, merqury is a framework to evaluate assemblies based on *k*-mer comparisons. It basically compares the assembly to reads. In other words, it checks whether the assembly faithfully represent the sequences that exist in the reads, using *k*-mers as evidence. It does so by computing assembly QV (sequence accuracy), *k*-mer completeness, etc. 

The genome size estimate based on jellyfish was 1.18Gb, but the size estimate from QUAST was 1.59Gb.
By comparing read k-mers and assembly k-mers, Merqury can tell you whether the genome size based on jellyfish output is an underestimation of the true genome size.

Let's begin by creating a separate conda environment for merqury. We also need to install one of the dependencies (meryl) separately in that conda environment.
```
# create a new conda env for merqury
conda create -n merqury
conda activate merqury

# install meryl
conda install bioconda::meryl==1.4.1

# install merqury
conda install bioconda::merqury
```
Then run Merqury like so:
```sh
#!/bin/bash
#SBATCH --job-name=merqury_ussuri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env 
source ~/.bash_profile
conda activate merqury

# inputs 
asm_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa
read_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/FASTQ/AMNH_21010_HiFi.fastq.gz

# parameters
K=21
THREADS=${SLURM_CPUS_PER_TASK}

# outputs 
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/merqury/noMito_k${K}
PREFIX=AMNH_21010_noMito_k${K}

mkdir -p "$out_dir"
cd "$out_dir"

echo "[INFO] Using K=$K threads=$THREADS"
echo "[INFO] Assembly: $asm_dir"
echo "[INFO] Reads:    $read_dir"
echo "[INFO] Output:   $out_dir/$PREFIX.*"

# sanity checks
command -v meryl >/dev/null 2>&1 || { echo "[ERROR] meryl not found in PATH"; exit 1; }
command -v merqury.sh >/dev/null 2>&1 || { echo "[ERROR] merqury.sh not found in PATH"; exit 1; }

# build read k-mer database
echo "[INFO] Counting k-mers in reads with meryl..."
meryl k=$K threads=$THREADS count "$read_dir" output reads.k${K}.meryl

# run Merqury == computes QV, completeness, spectra-cn, etc.
echo "[INFO] Running Merqury..."
merqury.sh reads.k${K}.meryl "$asm_dir" "$PREFIX" 2>&1 | tee "$PREFIX.merqury.log"

echo "[INFO] Done."
echo "[INFO] Key outputs to look at:"
echo "  - ${PREFIX}.qv            (assembly QV)"
echo "  - ${PREFIX}.completeness  (k-mer completeness)"
echo "  - ${PREFIX}.spectra-cn.*  (copy-number spectra plots/data)"
```
Check the results when the job finishes running. Let's check the .qv file first.
```
# cd into "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/merqury/noMito_k21" dirctory"
cat AMNH_21010_noMito_k21.qv 
```
![alt text](etc/qv.PNG)

Let's break it down:
- Column 1: Assembly name
- Column 2: Error *k*-mers
- Column 3: Total bases
- Column 4: QV score
- Column 5: Error rate

We can see the QV score is very high (> 70) and error rate is very low. This means that the estimated genome size based on QUAST is unlikely to be an overestimation due to sequencing artifacts/errors. 

Let's also check the *k*-mer completeness.
```
cat AMNH_21010_noMito_k21.completeness.stats
```
![alt text](etc/merqury_completeness_stats.PNG)
- Column 1: Assembly name
- Column 2: Subset ("all" means that completeness was calculated from all read *k*-mers)
- Column 3: Total read *k*-mers found in the assembly
- Column 4: Total *k*-mers found in the reads
- Column 5: *k*-mer completeness (%)

So, the *k*-mer completeness is very high. The ~9.4% of the *k*-mers in the reads are not present in the assembly likely due to: 1) Repeat complexity, 2) *k*-mers specific to alternative alleles at heterozygous loci not present in collapsed primary assembly, 3) Filtering based on read coverage. This does not mean that ~9.4% of the genome is missing. For these reasons, the haploid genome size estimated by jellyfish (1.18Gb) is likely and underestimation of the true genome size (1.6Gb) due to repeats and heterozygosity.

## 6) Scaffolding through Hi-C data incorporation
### Setup
```
# create a directory for scaffolding
mkdir -p /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding

# clone YaHS from github to the assembly project directory (/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo)
git clone https://github.com/c-zhou/yahs.git
cd yahs
make
```
## 7) Genome annotation
### Setup
Let's create new conda environments for packages to be used in genome annotation. Trimmomatic will be used for trimming Illumina adaptera. The funannotation package provides an automated pipeline for gene prediction, annotation, and comparison. The Earl Grey package automates transposable element annotation.
```
### create a conda environment and install funannotate
# add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create a conda env for trimmomatic and install it
conda create -n trimmomatic -c conda-forge -c bioconda trimmomatic

# create a conda env for funannotate and install it
conda create -n funannotate "python>=3.6,<3.9" funannotate
```
----------------------------------------------------------------------------------------------------
### RNA read QC (pre-trimming)
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
----------------------------------------------------------------------------------------------------
### Adapter trimming & post-trimming QC 
Use trimmomatic to trim adapters and then run FastQC on the trimmed reads. The RNA sequencing was done on Illumina NovaSeq X in a paired-end mode. We will use trimmomatic to trim the Illumina adapters. Since we did paired end sequencing on six different tissues, there are a total of 12 FASTQ files. Repeating trimmomatic independently for each tissue type is not very effective. Instead, I wrote a simple for loop to do the trimming in one go:

```sh
#!/bin/bash
#SBATCH --job-name=adapterTrim_ussuri
#SBATCH --nodes=1
#SBATCH --mem=100G
#SBATCH --partition=compute
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate trimmomatic

# paths to input forward & reverse reads, adapters, and output trimmed reads 
path_to_seq=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/FASTQ
adapters=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/custom_adapters.fa
out_path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/trimmed

# loop through the read files in the directory and run trimmomatic
# print this before looping
echo "start adapter trimming..."

# loop through the file in the directory and run trimmomatic
for f_read in ${path_to_seq}/*_1.fastq.gz; do
  
  # echo forward read
  echo "found a forward read: ${f_read##*/}"

  # designate reverse read
  r_read=${f_read/_1.fastq.gz/_2.fastq.gz}
  echo "found a corresponding reverse read: ${r_read##*/}"

  # print out a message on the type of tissue being processed
  tissue=${f_read%_1.fastq.gz}
  echo "Start adapter trimming ${tissue##*/} reads..."

  # run trimmomatic
  trimmomatic PE -threads ${SLURM_CPUS_PER_TASK} -phred33 \
    -Xmx80g \
    ${f_read} ${r_read} \
    ${tissue}_R1_paired.fastq.gz ${tissue}_R1_unpaired.fastq.gz \
    ${tissue}_R2_paired.fastq.gz ${tissue}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
done

# move trimmed files to the output directory
echo "move trimmed files to output dir..."

mv ${path_to_seq}/*_R1_paired.fastq.gz ${out_path}
mv ${path_to_seq}/*_R2_paired.fastq.gz ${out_path}
mv ${path_to_seq}/*_R1_unpaired.fastq.gz ${out_path}
mv ${path_to_seq}/*_R2_unpaired.fastq.gz ${out_path}

echo "all files moved to output dir"

# print this at the end
echo "Trimming on all tissue types finished successfully"
```
This run will result in a total of 24 files, two files (paired & unpaired) for each read. Now run FastQC again on trimmed read files:

```sh
#!/bin/bash
#SBATCH --job-name=rnaQC_posttrim
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

# paths to the trimmed fastq files
path_to_trimmed=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/trimmed

# output directory
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/FastQC/posttrim

# run FastQC
for file in ${path_to_trimmed}/*_paired.fastq.gz; do
  echo "running FastQC on $(basename "${file}")..."
  fastqc -o ${out_dir} -t ${SLURM_CPUS_PER_TASK} ${file}
  echo "FastQC on ${file##*/} completed"
done

# print when completed
echo "FastQC on all files completed"
```
Let's compare pre-trimming and post-trimming FastQC results, focusing on adapter content.
Let's open up the results for heart RNA reads:
![alt text](etc/adapters.png)
We can see that the adapters are basically gone after trimming.
   
It can be annoying to look through all the different .html files containing QC results for each tissue type. MultiQC (https://github.com/MultiQC/MultiQC) is a really neat tool that enables the user to merge outputs from different bioinformatics software to generate one, clean QC output.

Install MultiQC:
```txt
conda create -n multiqc -c conda-forge multiqc
conda activate multiqc
multiqc --help
```

You only need to supply MultiQC some options and path to the files you want to merge (e.g., path to multiple FastQC outputs).
```sh
# set directories as variables
outdir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/FastQC/multiqc
indir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/FastQC

# run multiqc == can run this on the head node because it runs very fast
multiqc -o ${outdir} --filename "posttrim_QC" ${indir}/posttrim
```
![alt text](etc/multiqc.PNG)

The MultiQC output gives a single, neatly organized .html file:
![alt text](etc/multiqc_out.PNG)

----------------------------------------------------------------------------------------------------
### (Pre-Hi-C) RNA alignment to draft using HiSat2
Now that we have our RNA-seq reads trimmed of adapter sequences and QCed, we can go ahead and align these reads to our draft genome to obtain mapping quality, etc. 

Let's create a conda environment for HiSat2 and install the package in it:
```
conda create -n hisat2 -c bioconda hisat2
conda activate hisat2
```

Also, let's copy paired reads over to our annotation directory. This is not entirely necessary, but I did it to keep the files neatly organized for different tasks.
```
# create a directory for annotation under the main assembly directory
# which is "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo
mkdir -p annotation
mkdir -p annotation/paired_RNAseq_reads

# copy the files
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/RNAseq/trimmed
cp *_R1_paired.fastq.gz /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads
cp *_R2_paired.fastq.gz /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads

# check files
ls /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads
```
![alt text](etc/paired_reads.PNG)

Now, let's run HiSat2 on Mendel with the script below
```sh
#!/bin/bash
#SBATCH --job-name=hisat2_preHiC
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=180:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
export PATH=/home/yshin/mendel-nas1/miniconda3/bin:$PATH
conda activate hisat2

# set paths for draft genome and RNA reads
draft_fa=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa
rna_path=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads

# output directory; create these dirs if they dont exist
outpath=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/hisat2_preHiC
mkdir -p ${outpath}
mkdir -p ${outpath}/{index,bam,logs}

# make index prefix 
index_prefix=${outpath}/index/Gloydius_ussuriensis_AMNH_21010_noMito

# build hisat2 index, skip if it already exist
if [[ ! -e ${index_prefix}.1.ht2 && ! -e ${index_prefix}.1.ht2l ]]; then
  echo "[INFO] Building HISAT2 index..."
  hisat2-build -p ${SLURM_CPUS_PER_TASK} ${draft_fa} ${index_prefix}
else
  echo "[INFO] HISAT2 index exists. Skipping."
fi

# map RNA reads to draft genome
# print out message
echo "[INFO] Mapping RNA-seq reads in: ${rna_path}"

# prevent loop from running if there are no files with matching names
shopt -s nullglob

# loop over all R1 files, get sample names, and get R2 file names
for R1 in ${rna_path}/AMNH_21010_*_R1_paired.fastq.gz; do
  sample=$(basename ${R1} _R1_paired.fastq.gz)
  R2=${rna_path}/${sample}_R2_paired.fastq.gz

  # skip if R2 is missing
  if [[ ! -f ${R2} ]]; then
    echo "[WARN] Missing R2 for ${sample}. Skipping."
    continue
  fi
  
  # print out some helpful message
  echo "[INFO] Mapping: ${sample}"

  # define output files
  bam=${outpath}/bam/${sample}.sorted.bam
  log=${outpath}/logs/${sample}.hisat2.log

  # run hisat with parameters specified above and then sort using samtools
  # --dta: means "downstream transcriptome assembly"; alignment report optimized for StringTie
  # -x: draft genome
  # -1: read 1 file
  # -2: read 2 file
  # 2> ${log}: direct & store output stats into log file
  # | samtools sort: pipe hisat2 output directly into samtools sort command
  # -: means "take piped input"
  # samtools index: index the sorted bam file
  hisat2 -p ${SLURM_CPUS_PER_TASK} --dta \
    -x ${index_prefix} \
    -1 ${R1} -2 ${R2} \
    2> ${log} \
  | samtools sort -@ ${SLURM_CPUS_PER_TASK} -o ${bam} -

  samtools index ${bam}
done

# print when done
echo "[DONE] HISAT2 mapping complete."
echo "[INFO] BAMs in: ${outpath}/bam"
echo "[INFO] Logs in: ${outpath}/logs"
```

We can open up the log files to check the overall alignment rate per tissue type. This is as follows:
  - Heart: 64.79%
  - Lung: 78.11%
  - Kidney: 82.91%
  - Liver: 87.25%
  - Muscle: 88.53%
  - Skin: 85.89%

----------------------------------------------------------------------------------------------------
### (Pre-Hi-C) Draft-guided transcriptome assembly using StringTie
Now. let's use StringTie to conduct transcriptome assembly guided by the draft genome. Let's create a conda environment for StringTie and install the package:
```
conda create -n stringtie -c bioconda stringtie
```

Next, run StringTie on Mendel with the script below. 
```sh
#!/bin/bash
#SBATCH --job-name=stringtie_preHiC
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G
#SBATCH --time=180:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
export PATH=/home/yshin/mendel-nas1/miniconda3/bin:$PATH
conda activate stringtie

# dir to bam files
bam_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/hisat2_preHiC/bam

# output dir
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/StringTie_preHiC
mkdir -p ${out_dir}
mkdir -p ${out_dir}/gtf ${out_dir}/merge 

# run stringtie // loop through file per tissue and assemble per tissue
# each loop takes one bam file (RNA alignment), assembles a transcriptome, and outputs a gtf file for a given tisue type
for bam in ${bam_dir}/*.sorted.bam; do
  sample=$(basename "$bam". sorted.bam)
  echo "assembling transcripts from ${sample}..."
  stringtie "$bam" -p ${SLURM_CPUS_PER_TASK} -o ${out_dir}/gtf/${sample}.gtf
  echo "done assembling transcripts from ${sample}..."
done

# merge gtf per tissue into a single gtf
echo "merging transcripts"
ls ${out_dir}/gtf/*.gtf > ${out_dir}/merge/gtf_list.txt
stringtie --merge -p ${SLURM_CPUS_PER_TASK} -o ${out_dir}/merge/merged.gtf ${out_dir}/merge/gtf_list.txt

# print when done
echo "DONE"
```
Now, let's see how many transcripts, exons, and genes we've got from this run:
```
# check the number of transcripts (counting the number of unique transcript_id) and number of exons
cut -f3 merged.gtf | sort | uniq -c | sort -nr | head

# check the number of genes (counting the number of unique gene_id)
awk '$3=="transcript" {match($0,/gene_id "[^"]+"/,a); if(a[0]!=""){gsub(/gene_id "|"/,"",a[0]); print a[0]}}' merged.gtf | sort -u | wc -l
```
The outputs will show that we have:
  - 381,058 exons
  - 37,147 transcripts
  - 20,527 genes

These results seem to be in great shape in terms of transcript evidence-building. We will put a pause here, finish Hi-C incorporation, and move on to final annotation once our assembly is at the scaffold level. 

----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Repeat masking (soft masking) using Earl Grey
Before moving on to final annotation, it is necessary to annotate the repeats and soft mask the genome. We can do this using the Earl Grey pipeline (https://github.com/TobyBaril/EarlGrey), which is a fully-automated pipeline for transposable element (TE)/repeat annotation.

Let's start by creating a separate conda environment and install Earl Grey: 

```
# create a conda environment for Earl Grey and install it
conda create -n earlgrey -c conda-forge -c bioconda earlgrey=7.0.1
```

Once this is done, run Earl Grey with the script below on Mendel:
```sh
#!/bin/bash
#SBATCH --job-name=earlGrey_ussuri
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=180:00:00
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/outfiles/slurm-%x_%j.err

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
```
----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Re-run HiSat2 and StringTie on the scaffolded and masked genome
----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Structural annotation
----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Functional annotation

## 8) Mitogenome assembly
### "Manual" annotation with MITOS2
Running mitofifi will automatically give you the annotated mitogenome. This is an alternative method for annotating the mitogenome that is more time consuming (but I guess it is still valuabe for the purpose of learning). 

In the "Contamination screening" section above, we identified and stored the mitochondrial contig into a separate fasta file containing a single copy mitogenome ("mito_singlecopy.fa"). Now let's annotate this file.

First, let's check the fasta header (we are in the "genome_cleanup" directory):
```txt
head mito_singlecopy.fa 
```
![alt text](etc/head_mitogenome_fa.PNG)

We can see that the header is the contig name. This is not ideal; let's change the header with species name:
```txt
sed 's/^>.*/>Gloydius_ussuriensis_mitogenome/' mito_singlecopy.fa > ussuri_mitogenome_ann.fa
head ussuri_mitogenome_ann.fa 
```
![alt text](etc/head_mitogenome_ann_fa.PNG)

Now, scp this file to your local device.

__*The following steps are (mostly) run on a local device*__

Now that this step is done, we can annotate the assembled mitogenome using MITOS2 (https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu%2Frepos%2Fiuc%2Fmitos2%2Fmitos2%2F2.1.3%20galaxy0).

Use the following settings:
1) Input file: ussuri_mitogenome_ann.fasta
2) In file uploads, use "auto-detect" and do not specify the reference
3) Use "Vertebrate (2)" genetic code and "RefSeq63 Metazoa" for reference data
4) Select BED, GFF, and nucleotide FASTA, and zipped raw results as outputs
5) Turn on email notification
6) Turn on "Attempt to re-use jobs with identical parameters?"

MITOS2 can be run through the run on the Galaxy Server. The window looks something like this:
![alt text](etc/mitos2.PNG)

After MITOS2 run is complete, run a quick sanity check on the output GFF (annotation) file. First, download all output files in the "mitos2" directory in your local device and cd into it (in my case the directory in my local device is: /home/yshin/Gloydius_ussuriensis_genome_assembly/outfiles/mito_assembly/mitos2)
```txt
grep -v "^#" "Galaxy7-[MITOS2 on dataset 5_ GFF].gff" | cut -f3 | sort | uniq -c
```
![alt text](etc/check_mitos2_gff1.PNG)

Looks good! Also check if the names of all 13 protein-coding genes:
```
grep -v "^#" "Galaxy7-[MITOS2 on dataset 5_ GFF].gff" | awk -F'\t' '$3=="gene"{print $9}' | \
  grep -Eo 'cox1|cox2|cox3|cob|cytb|atp6|atp8|nad1|nad2|nad3|nad4l|nad4|nad5|nad6' | sort | uniq -c
```
![alt text](etc/check_mitos2_gff2.PNG)

Sweet. Now let's rotate the mitogenome FASTA based on the location of tRNA-Phe (trnF)
```txt
grep -i -E "trnF|tRNA-Phe|phenylalanine" "Galaxy7-[MITOS2 on dataset 5_ GFF].gff" 
```
![alt text](etc/check_mitos2_gff3.PNG)

We can see that trnF starts at position 15284 and ends at position 15347. __BUT,__ because the orientation of this sequence is (-), this means that the actual 5' end is at position 15347 and the 3' is at position 15284. 

Now we can rotate the sequence based on this information. To do so, first install the "rotate" package under the mendel-nas1 directory:

```txt
###  this is run on the cluster
# install
git clone https://github.com/richarddurbin/rotate.git ; cd rotate ; make

# to run, cd into the folder containing rotate and then run:
./rotate
```

Once the package is installed, run it on the head node. Since the start position orientation is on the (-) strand, make position 15347 as the start position. In addition, because the strand orientation is (-), make sure to reverse complement when you rotate (use the -rc flag).

```txt
### run rotate with the input sequence
# path to seq 
path_to_seq=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

# output dir
out_dir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/mt_genome

# run rotate. run this in the rotate directory under mendel-nas1/
./rotate -x 15347 -rc ${path_to_seq}/ussuri_mitogenome_ann.fa > ${out_dir}/ussuri_mt_rotated.fasta

# check stats
conda activate genome_assembly
seqkit stats ${out_dir}/ussuri_mt_rotated.fasta
```
The output should look something like this:
![alt text](etc/rotate_checkstat.PNG)

Awesome. Now, scp the rotated fasta the directory in your local device and re-run MITOS2 to get GFF of the rotated FASTA. Once this is done, run some final checks

```txt
###  this is run on a local device
# cd into mito assembly directory
cd outfiles/mito_assembly/rotated/

# make directory for the final output files
mkdir final

# copy the "ussuri_mt_rotated.fasta" and the .gff file to the "final" folder
cp ussuri_mt_rotated.fasta final/
cp *.gff final/

# change the gff file name
mv 'Galaxy28-[MITOS2 on dataset 26_ GFF].gff' annotation_rotated_final.gff

# cd into final folder and change file names and confirm gene content again
cd final/
grep -v "^#" *.gff | cut -f3 | sort | uniq -c
```
![alt text](etc/final_check_mitos2_gff1.PNG)

```
# check trnF position
grep -i -E "gene_trnF|Name=trnF|trnF" *.gff
```
![alt text](etc/final_check_mitos2_gff2.PNG)

The trnF is now shown to span 1 bp - 64 bp.This is exactly what we expect from a correctly rotated sequence. Check the sequence stats in the local device using seqkit (install it if you need to).

```txt
seqkit stats ussuri_mt_rotated.fasta
```
![alt text](etc/final_mt_seqstat.PNG)

Awesome! this can now be imported into Geneious for visualization.
![alt text](etc/mt_geneious_vis.PNG)

__Note:__ The steps outlined above already have a lot of manual work. But then it needs further manual curation to adjust annotation. But we don't have to go through all of this for this sample because MitoHiFi run we did in the "Contamination screening" step will automatically produce fully annotated mitogenome that requires minimal manual data processing.
![alt text](etc/mitohifi_geneious.PNG)

From here, we only need to manually annotate the replication origin and two D-loop regions. I did this by extracting these sequences from the reference mitogenome and aligning them to the MitoHiFi assembly. The final output looks like this (after some name edits):
![alt text](etc/fullMt_MitoHiFi.PNG)

Note that this genome is also 17,211bp long, like the assembly generated through the steps in this subsection. But overall, MitoHiFi output produced a cleaner output. This is because the assembly we annotated with MITOS2 had some issues determining gene boundaries in some places. Let's look at this by comparing the MitoHiFi assembly, MITOS2 assembly, and the reference mitogenome:
![alt text](etc/gene_boundaries.PNG)

Finally, recall that our mitogenome was contained in the contig ptg000073c and that this contig was 68,844 bp in length. Since our mitogenome is 17.211 bp in length, this mean that the contig has four complete mitogenomes stitched back to back four times. Let's make sure this is actually the case by splitting up this contig into four equal chunks and see if each of them is indeed a complete mitogenome.

First, let's cd into the "genome_cleanup" directory on Mendel, which is where the ptg000073c contig fasta file is located in. Activate the "genome_assembly" conda env and use seqkit to check the length of this contig:
```txt
# in the "genome_cleanup" directory
conda activate genome_assembly
seqkit stats ptg000073c.fa
```
![alt text](etc/mito_dup.PNG)

Now, let's create a "mito_chunks" directory under the current directory to store the split files, and use the samtools faidx command to slice this contig into four chunks, each exactlty 17,211 bp in length. 
```txt
# make the directory
mkdir mito_chunks

# copy the contig file into this directory
cp ptg000073c.fa mito_chunks

# cd into the mito_chunks directory
cd mito_chunks/

# index the fasta first
samtools faidx ptg000073c.fa

# set contig name as variable
c_name="ptg000073c"

# split
samtools faidx ptg000073c.fa ${c_name}:1-17211     > mito_part1.fa
samtools faidx ptg000073c.fa ${c_name}:17212-34422 > mito_part2.fa
samtools faidx ptg000073c.fa ${c_name}:34423-51633 > mito_part3.fa
samtools faidx ptg000073c.fa ${c_name}:51634-68844 > mito_part4.fa
```
Let's use seqkit to verify the size of each chunk
```txt
seqkit stats mito_part*
```
![alt text](etc/mito_chunk_size.PNG)

scp this folder to our local device and run blast on each of them. The results will show that they all match our *G. ussuriensis* reference mitogenome with > 99% percent identity with 100% query cover. This confirms that mitogenomes were stitched four times on to the contig ptg000073c. This likely happened because the assembler failed to recognize the circular nature of the mitogenome.

----------------------------------------------------------------------------------------------------
### Submitting mitogenome to GenBank
Let's submit the completed mitogenome assembly to GenBank. Prior to this step. I loaded the final mitogenome output from mitohifi on to Geneious Prime and manually annotated the two D-loops and replication origin. I did this by extracting these sequences from the conspecific reference mitogenome and mapping them to the assembled mitogenome using the "Map to reference" tool. I then downloaded this as a fasta file. This file is stored in a new directory at "/home/yshin/Gloydius_ussuriensis_genome_assembly/outfiles/mitogenome_GenBank_submission"

Let's also copy the "final_mitogenome.gb" file from the mitohifi output directory to this directory:
```txt
# cd into mitohifi output directory
cd /home/yshin/Gloydius_ussuriensis_genome_assembly/outfiles/mitohifi

# copy .gb file over to mitogenome submission directory
cp final_mitogenome.gb /home/yshin/Gloydius_ussuriensis_genome_assembly/outfiles/mitogenome_GenBank_submission

# cd into mitogenome submission directory
cd /home/yshin/Gloydius_ussuriensis_genome_assembly/outfiles/mitogenome_GenBank_submission
```
