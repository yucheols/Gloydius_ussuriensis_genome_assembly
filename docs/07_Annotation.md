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
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

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
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

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
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

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
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

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
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

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
### (Pre-Hi-C) Venom gland transcriptome data
Park et al. (2026) produced venom gland transcriptome data for the three *Gloydius* species found in South Korea. This can be useful for placing expressed venom gene candidates to specific scaffolds.

First, create a new directory under the annotation directory
```sh
mkdir -p venom_gland venom_gland/ncbi_seq venom_gland/fastq
```
Then, download NCBI SRA Toolkit on Mendel. Go to the following link (https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit) and copy the link for AlmaLinux 64 bit architecture. Then:
```sh
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/sratoolkit.3.4.1-alma_linux64.tar.gz
```
You should see a tar.gz file created under the annotation directory. Untar this file like so:
```sh
tar -xzvf sratoolkit.3.4.1-alma_linux64.tar.gz
rm *.tar.gz
```
Then add the executable to the PATH variable:
```sh
export PATH=$PWD/sratoolkit.3.4.1-alma_linux64/bin:$PATH
```

Then submit the script below to SLURM.
```sh
#!/bin/bash
#SBATCH --job-name=get_venom_data
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=50G
#SBATCH --time=25:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

set -euo pipefail

### commands start here ###

# path to SRA toolkit
export PATH=$PWD/sratoolkit.3.4.1-alma_linux64/bin:$PATH

# set directories
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland"
sradir="${basedir}/ncbi_seq"
fastqdir="${basedir}/fastq"
tmpdir="${basedir}/tmp"

mkdir -p "$sradir" "$fastqdir" "$tmpdir"

# download SRA files
prefetch SRR35908235 -O "$sradir"
prefetch SRR35908238 -O "$sradir"

# convert .sra files to FASTQ files
fasterq-dump "$sradir/SRR35908235/SRR35908235.sra" \
  --split-files \
  --threads 12 \
  --temp "$tmpdir" \
  -O "$fastqdir"

fasterq-dump "$sradir/SRR35908238/SRR35908238.sra" \
  --split-files \
  --threads 12 \
  --temp "$tmpdir" \
  -O "$fastqdir"

# gzip FASTQ files
gzip "$fastqdir"/*.fastq
```

After this, run a quick sanity check from the fastq directory.
```sh
# verify the files have sequences
zcat SRR35908235_1.fastq.gz | head -n 4
zcat SRR35908235_2.fastq.gz | head -n 4

zcat SRR35908238_1.fastq.gz | head -n 4
zcat SRR35908238_2.fastq.gz | head -n 4

# check file size == R1/R2 files should have the same size
ls -lh
```
Then QC the files using FastQC and MultiQC.
```sh
#!/bin/bash
#SBATCH --job-name=venom_gland_qc
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

set -euo pipefail

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# set directory
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland"
cd "$basedir"

# make output dir
mkdir -p qc_fastqc

# run fastqc
fastqc fastq/*.fastq.gz \
  -o qc_fastqc \
  -t 12

# activate multiqc conda env 
conda activate multiqc

# run multiqc
multiqc qc_fastqc \
  -o qc_fastqc
```

The MultiQC report shows considerable adapter content. Use fastp for automatic adapter detection and trimming for Illumina paired-end reads.
```sh
#!/bin/bash
#SBATCH --job-name=venom_trim
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

source ~/.bash_profile
conda activate scaffolding

# set directories
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland"
indir="${basedir}/fastq"
outdir="${basedir}/trimmed_fastq"
reportdir="${basedir}/fastp_reports"

mkdir -p "$outdir" "$reportdir"

ACCESSIONS=(
  SRR35908235
  SRR35908238
)

for acc in "${ACCESSIONS[@]}"; do
  echo "Trimming ${acc}"

  fastp \
    -i "${indir}/${acc}_1.fastq.gz" \
    -I "${indir}/${acc}_2.fastq.gz" \
    -o "${outdir}/${acc}_1.trimmed.fastq.gz" \
    -O "${outdir}/${acc}_2.trimmed.fastq.gz" \
    --detect_adapter_for_pe \
    --thread 8 \
    --html "${reportdir}/${acc}.fastp.html" \
    --json "${reportdir}/${acc}.fastp.json"

  echo "Finished ${acc}"
done
```

After this is done, run FastQC and MultiQC again on trimmed reads.
```sh
#!/bin/bash
#SBATCH --job-name=venom_gland_posttrim_qc
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate genome_assembly

# set directory
basedir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland"
cd "$basedir"

# make output dir
mkdir -p postrrim_qc_fastqc

# run fastqc
fastqc trimmed_fastq/*.fastq.gz \
  -o postrrim_qc_fastqc \
  -t 12

# activate multiqc conda env 
conda activate multiqc

# run multiqc
multiqc postrrim_qc_fastqc \
  -o postrrim_qc_fastqc
```

----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Repeat masking (soft masking) using Earl Grey
Before moving on to final annotation, it is necessary to annotate the repeats and soft mask the genome. We can do this using the Earl Grey pipeline (https://github.com/TobyBaril/EarlGrey), which is a fully-automated pipeline for transposable element (TE)/repeat annotation.

Let's start by creating a separate conda environment and install Earl Grey: 

```
# create a conda environment for Earl Grey and install it
conda create -n earlgrey -c conda-forge -c bioconda earlgrey=7.0.1
```

Once this is done, run Earl Grey with the script below on Mendel (make sure the input genome is the final, curated chromosome-level assembly):
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
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# soft masking the scaffolded assembly with Earl Grey

# activate the conda env
source ~/.bash_profile
conda activate earlgrey

# set variables
GENOME="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
SPECIES="Gloydius_ussuriensis"

# output path
outpath=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked
mkdir -p $outpath

# run Earl Grey
echo "Starting EarlGrey: $(date)"

# -d flag == Create soft-masked genome at the end? (yes/no, Default: no)
earlGrey -g ${GENOME} -s ${SPECIES} -o ${outpath} -d yes -t ${SLURM_CPUS_PER_TASK}

echo "Finished EarlGrey: $(date)"

# list outputs
echo "EarlGrey output files:"
find "$outdir" -maxdepth 3 -type f | head -n 50

echo "Done: $(date)"
```
----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Annotation using funannotate
The general annotation workflow is as follows:
```
1. funannotate train: organ tissue RNA-seq data + publicly available venom gland RNA-seq data
2. funannotate predict
3. funannotate update
4. funannotate fix
5. InterProScan on proteins predicted by funannotate predict
6. funannotate annotate
7. Final QC / manual curation
```

Let's start by symlinking the softmasked fasta for annotation.
```sh
# symlink softmasked fasta output from earl grey
ln -s /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta \
  Gloydius_ussuriensis.softmasked.fa
```
Then activate the funannotate conda env and install gffread.
```sh
conda activate funannotate
conda install bioconda::gffread
```

####  Setup funannotate 1: GeneMark
funannotate also requires an initial database setup step. Run setup as below:
```sh
# create db path
db_path="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/an
notation/funannotate_db"

mkdir -p ${db_path}

# check dependencies
funannotate check --show-versions

# setup db
funannotate setup -d ${db_path} -i all --wget

# run this once the setup is done
echo 'export FUNANNOTATE_DB=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db' >> ~/.bash_profile

# also set up BUSCO lineage db 
# this version of funannotate do not support the sauropsida lineage DB; use tetrapoda db instead 
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
funannotate setup -d "$FUNANNOTATE_DB" -b tetrapoda --wget
```
After this, if you run "funannotate check --show-versions", you'll see the current configuration is missing some key dependencies, such as eggnog mapper ("ERROR: emapper.py not installed"), GeneMark ("ERROR: gmes_petap.pl not installed"), and signalp ("ERROR: signalp not installed"). 

GeneMark cannot be installed through conda. We need to manually download the .tar.gz file from the GeneMark website and scp it to Mendel.
```
# on Mendel
mkdir -p "/home/yshin/mendel-nas1/gmes_linux_64"
cd "/home/yshin/mendel-nas1/gmes_linux_64"

# in local device
scp gmes_linux_64_4.tar.gz yshin@mendel.sdmz.amnh.org:/home/yshin/mendel-nas1/gmes_linux_64

# unzip on Mendel
tar -xvzf gmes_linux_64_4.tar.gz
```

Now, set the GeneMark path:
```
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"
```

Also, download the key file from the website (Linux 64 bit) and scp it to the Mendel home directory. After that, run:
```
zcat gm_key_64.gz > ~/.gm_key
chmod 600 ~/.gm_key
```

Then, run a full check
```sh
conda activate funannotate

export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"

funannotate check --show-versions
```
Now you shouldn't see the GeneMark error.

####  Setup funannotate 2: eggNOG-mapper
```
# in the funannotate conda env 
conda install -c conda-forge -c bioconda eggnog-mapper

# download eggnog db
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"
mkdir -p "$EGGNOG_DATA_DIR"
cd "$EGGNOG_DATA_DIR"

BASE="http://eggnog5.embl.de/download/emapperdb-5.0.2"

wget -c -O eggnog.db.gz "$BASE/eggnog.db.gz"
wget -c -O eggnog.taxa.tar.gz "$BASE/eggnog.taxa.tar.gz"
wget -c -O eggnog_proteins.dmnd.gz "$BASE/eggnog_proteins.dmnd.gz"

# decompress
gunzip -f eggnog.db.gz
gunzip -f eggnog_proteins.dmnd.gz
tar -xvzf eggnog.taxa.tar.gz
rm -f eggnog.taxa.tar.gz

# check installation
which emapper.py
emapper.py --version
```

####  Setup funannotate 3: SignalP
Install signalp from the software website and scp it to Mendel. After that:
```
# make install directory
mkdir -p /home/yshin/mendel-nas1/signalp

# unpack
tar -xvzf /home/yshin/mendel-nas1/signalp-6.0i.fast.tar.gz -C /home/yshin/mendel-nas1/signalp

# cd into the package dir
/home/yshin/mendel-nas1/signalp/signalp6_fast/signalp-6-package

# put pip cache and temp files on nas, not /home, to prevent disk quota error
mkdir -p /home/yshin/mendel-nas1/tmp/pip_cache
mkdir -p /home/yshin/mendel-nas1/tmp/pip_tmp

export PIP_CACHE_DIR="/home/yshin/mendel-nas1/tmp/pip_cache"
export TMPDIR="/home/yshin/mendel-nas1/tmp/pip_tmp"

# install
python -m pip install --no-cache-dir .

# symlink for funannotate
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

cd "$CONDA_PREFIX/bin"
ln -sf signalp6 signalp
```

Then run a final check. This should print "All 37 external dependencies are installed" at the end.
```
funannotate check --show-versions
```
Also, copy the "distilled_model_signalp6.pt" file into the signalp installation "model_weights" dir.
```sh
cp -av \
/home/yshin/mendel-nas1/signalp/signalp6_fast/signalp-6-package/models/distilled_model_signalp6.pt \
/home/yshin/mendel-nas1/miniconda3/envs/funannotate/lib/python3.8/site-packages/signalp/model_weights/
```

Check whether the file is correctly copied.
```sh
ls -lh \
/home/yshin/mendel-nas1/miniconda3/envs/funannotate/lib/python3.8/site-packages/signalp/model_weights/
```

####  Install InterProScan
First, create a new directory to store interproscan
```sh
cd /home/yshin/mendel-nas1

mkdir -p interproscan
cd interproscan
```
Next, install interproscan and checksum files.
```sh
# download interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz

# download checksum
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz.md5
```
Check download:
```
md5sum -c interproscan-5.78-109.0-64-bit.tar.gz.md5
```

Extract the file while preserving permissions:
```sh
# extract
tar -pxvzf interproscan-5.78-109.0-64-bit.tar.gz

# check script
ls -lh /home/yshin/mendel-nas1/interproscan/interproscan-5.78-109.0/interproscan.sh
```

Also check the runtime requirements. What we want is Java 11 or above.
```sh
java -version
python3 --version
perl -version | head
echo "$JAVA_HOME"
```

####  funannotate step 1: Train
Run with the script below.
```sh
#!/bin/bash
#SBATCH --job-name=funannotate_train
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=120:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate environment and set up analyses
source ~/.bash_profile
conda activate funannotate
set -euo pipefail

# make TransDecoder utilities available to PASA
export PATH="$CONDA_PREFIX/opt/transdecoder:$CONDA_PREFIX/opt/transdecoder/util:$PATH"

# avoid system library conflicts
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# use a short temporary directory and clean it up on exit
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT

# set paths for input and output
# softmakesked genome assembly
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# path to paired-end RNA-seq reads from organ tissues and venom gland
tissue_reads="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/paired_RNAseq_reads"
venom_reads="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq"

# make output dir for all funannotate results
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate/"
mkdir -p "$outdir"

# Train to align RNA-seq data, run Trinity, and then run PASA
funannotate train -i ${genome} -o ${outdir} \
    --left ${tissue_reads}/AMNH_21010_Ht_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lr_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Ky_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Me_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lg_R1_paired.fastq.gz ${tissue_reads}/AMNH_21010_Skin_R1_paired.fastq.gz ${venom_reads}/SRR35908235_1.trimmed.fastq.gz \
    --right ${tissue_reads}/AMNH_21010_Ht_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lr_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Ky_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Me_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Lg_R2_paired.fastq.gz ${tissue_reads}/AMNH_21010_Skin_R2_paired.fastq.gz ${venom_reads}/SRR35908235_2.trimmed.fastq.gz \
    --species "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --stranded RF \
    --max_intronlen 100000 \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --memory 250G \
    --no_trimmomatic
```

####  funannotate step 2: Predict
```sh
#!/bin/bash
#!/bin/bash
#SBATCH --job-name=funannotate_predict
#SBATCH --nodes=1
#SBATCH --partition=bigmem
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# funannotate predict
#
# GlimmerHMM is disabled because the previous run produced
# 431,192 GlimmerHMM predictions, causing an excessively large
# EVM input and contributing to an OOM failure.
# ============================================================

# activate environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail

# avoid system library conflicts
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="${GENEMARK_PATH}:${PATH}"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# temporary directory
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT


### set paths for input and output
# input genome
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# output directory
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"

# run funannotate predict
echo
echo "============================================================"
echo "Starting funannotate predict"
echo "GlimmerHMM weight = 0 (disabled)"
echo "============================================================"
echo

funannotate predict \
    -i "${genome}" \
    -o "${outdir}" \
    --species "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --busco_db tetrapoda \
    --organism other \
    --busco_seed_species Taeniopygia_guttata \
    --max_intronlen 100000 \
    --repeats2evm \
    --weights glimmerhmm:0 \
    --cpus "${SLURM_CPUS_PER_TASK}"

# print when finished
echo
echo "============================================================"
echo "funannotate predict completed successfully"
echo "End time: $(date)"
echo "============================================================"
```

####  funannotate step 3: Update
```sh
#!/bin/bash
#SBATCH --job-name funannotate_update
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate environment and set up analyses
source ~/.bash_profile
conda activate funannotate
set -euo pipefail

# avoid system library conflicts
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# funannotate databases and external tools
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"
export GENEMARK_PATH="/home/yshin/mendel-nas1/gmes_linux_64/gmes_linux_64_4"
export PATH="$GENEMARK_PATH:$PATH"
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"

# use a short temporary directory and clean it up on exit
export TMPDIR="/tmp/yshin_fun_${SLURM_JOB_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"
mkdir -p "$TMPDIR"

trap 'rm -rf "$TMPDIR"' EXIT

# output directory for funannotate results
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate/"

# run funannotate update
funannotate update -i ${outdir} --cpus $SLURM_CPUS_PER_TASK
```

####  funannotate step 4: InterProScan
Run interproscan through "funannotate iprscan" using the script below:
```sh
#!/bin/bash
#SBATCH --job-name=funannotate_iprscan
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=720:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# ------------------------------------------------------------
# set paths
# ------------------------------------------------------------

FUN_DIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
IPR_DIR="/home/yshin/mendel-nas1/interproscan/interproscan-5.78-109.0"
IPRSCAN="${IPR_DIR}/interproscan.sh"


# ------------------------------------------------------------
# run InterProScan through funannotate
# ------------------------------------------------------------

cd "${FUN_DIR}"

echo "Starting InterProScan: $(date)"

funannotate iprscan \
    -i "${FUN_DIR}" \
    -m local \
    --iprscan_path "${IPRSCAN}" \
    -c 4


# ------------------------------------------------------------
# verify output
# ------------------------------------------------------------

IPR_XML="${FUN_DIR}/annotate_misc/iprscan.xml"

echo
echo "============================================================"

if [[ -s "${IPR_XML}" ]]; then
    echo "InterProScan completed successfully."
    echo "Output:"
    ls -lh "${IPR_XML}"
else
    echo "ERROR: Expected InterProScan XML was not produced:"
    echo "${IPR_XML}"
    exit 1
fi


# ------------------------------------------------------------
# print when finished
# ------------------------------------------------------------

echo "Finished: $(date)"
echo "============================================================"
```

####  funannotate step 5: Annotate 
Using the outputs from "funannotate predict" and interproscan runs, now run funannotate annotate step:
```sh
#!/bin/bash
#SBATCH --job-name=funannotate_annotate
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=350G
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ------------------------------------------------------------
# activate conda environment
# ------------------------------------------------------------

source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# ------------------------------------------------------------
# environment
# ------------------------------------------------------------

# avoid system libstdc++ conflict
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

# Funannotate databases
export FUNANNOTATE_DB="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate_db"

# EggNOG database
export EGGNOG_DATA_DIR="/home/yshin/mendel-nas1/eggnog_db"


# ------------------------------------------------------------
# temporary directory
# ------------------------------------------------------------

export TMPDIR="/tmp/yshin_funannotate_${SLURM_JOB_ID}"
export TEMP="${TMPDIR}"
export TMP="${TMPDIR}"

mkdir -p "${TMPDIR}"
trap 'rm -rf "${TMPDIR}"' EXIT


# ------------------------------------------------------------
# set paths
# ------------------------------------------------------------

FUN_DIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"
IPR_XML="${FUN_DIR}/annotate_misc/iprscan.xml"
PROTEINS="${FUN_DIR}/update_results/Gloydius_ussuriensis_AMNH_21010.proteins.fa"


# ------------------------------------------------------------
# run Funannotate annotate
# ------------------------------------------------------------

cd "${FUN_DIR}"

echo "============================================================"
echo "Starting funannotate annotate"
echo "Date: $(date)"
echo "============================================================"
echo

funannotate annotate \
    -i "${FUN_DIR}" \
    -s "Gloydius ussuriensis" \
    --isolate AMNH_21010 \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --busco_db tetrapoda \
    --database "${FUNANNOTATE_DB}" \
    --tmpdir "${TMPDIR}"


# ------------------------------------------------------------
# verify annotation output
# ------------------------------------------------------------

ANNOTATE_DIR="${FUN_DIR}/annotate_results"

echo
echo "============================================================"
echo "Checking Funannotate annotation output"
echo "============================================================"
echo

if [[ ! -d "${ANNOTATE_DIR}" ]]; then
    echo "ERROR: annotate_results directory was not created:"
    echo "${ANNOTATE_DIR}"
    exit 1
fi

echo "annotate_results created successfully:"
ls -lh "${ANNOTATE_DIR}"

echo
echo "============================================================"
echo "Funannotate annotate completed"
echo "Finished: $(date)"
echo "============================================================"
```

#### Post-annotation checks
Run the following script to QC the annotation results:
```sh
funannotate_final_qc.sh
```

__NOTE:__ This script uses a newer compleasm version, which is installed in a separate conda env created as:
```sh
conda create -n compleasm_v0.2.9 \
    -c conda-forge \
    -c bioconda \
    compleasm=0.2.9
```

----------------------------------------------------------------------------------------------------
### (Post-Hi-C) Toxin gene annotation using ToxCodAn-Genome
Snake toxin gene annotation can be tricky. [...]

Install ToxCodAn-Genome and its dependencies.
```sh
# create conda env and install dependencies
conda create -n ToxcodanGenome -c bioconda python biopython pandas blast exonerate miniprot gffread hisat2 samtools stringtie trinity spades

# install ToxCodAn-Genome from github and add to PATH
cd mendel-nas1/
git clone https://github.com/pedronachtigall/ToxCodAn-Genome.git
echo "export PATH=$PATH:$PWD/ToxCodAn-Genome/bin/" >> ~/.bashrc

# test
conda activate ToxcodanGenome
cd ToxCodAn-Genome/bin/
python toxcodan-genome.py -h

# apply permission to all executables
chmod +x /home/yshin/mendel-nas1/ToxCodAn-Genome/bin/*
```
Also, create a separate conda env and install ToxCodAn to be used for transcriptome annotation.
```sh
# setup conda env
conda create -n Toxcodan \
  --strict-channel-priority \
  -c conda-forge \
  -c bioconda \
  python=3.10 \
  codan=1.2 \
  blast \
  hmmer

# git clone toxcodan from repo
git clone https://github.com/pedronachtigall/ToxCodAn.git
export PATH=$PATH:path/to/ToxCodAn/bin/

# unzip models (in "/home/yshin/mendel-nas1/ToxCodAn" dir)
unzip models.zip

# download the SignalP-4.1 from its website, decompress, and add it to PATH:
tar -xzf signalp-4.1g.Linux.tar.gz
export PATH="$PATH:/home/yshin/mendel-nas1/signalp-4.1"

### signalp 4.1 setup
# cd into dir
cd /home/yshin/mendel-nas1/signalp-4.1

# confirm that the required Perl module exists
ls -l lib/FASTA.pm
$ENV{SIGNALP} = '/home/yshin/mendel-nas1/signalp-4.1/';

# back up and fix the script
cp signalp signalp.original

sed -i \
  's|/usr/opt/www/pub/CBS/services/SignalP-4.1/signalp-4.1|/home/yshin/mendel-nas1/signalp-4.1|g' \
  signalp

# verify that the correct path is now embedded
grep -n 'ENV{SIGNALP}' signalp

### apply "execution permission" to all bin executables
chmod 777 /home/yshin/mendel-nas1/ToxCodAn/bin/*
```

First, run TRassembly.py included in the ToxCodAn-Genome package for venom gland transcriptome assembly.

Before running TRassembly, make sure the jellyfish k-mer counter is installed and accessible by Trinity (which is part of the TRassembly script). In my case, jellyfish was not installed in the ToxcodanGenome conda env and this caused TRassembly script to crash as soon as it entered the de novo transcriptome assembly step. Also, install bowtie/bowtie2 and salmon alongside jellyfish.
```sh
# activate conda env
conda activate ToxcodanGenome

# install
conda install --freeze-installed \
    -c conda-forge \
    -c bioconda \
    bioconda::kmer-jellyfish \
    bioconda::bowtie2 \
    bioconda::bowtie \
    bioconda::salmon
```
Now, run the TRassembly.py script:
```sh
#!/bin/bash
#SBATCH --job-name=TRassembly
#SBATCH --nodes=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env 
# TRassembly.py is a part of ToxCodAn-Genome
source ~/.bash_profile
conda activate ToxcodanGenome

# force UTF-8 output in non-interactive SLURM jobs
# locale available on all standard Linux nodes
export LANG=C
export LC_ALL=C

# force Python itself to use UTF-8
export PYTHONIOENCODING=UTF-8
export PYTHONUTF8=1

set -euo pipefail

# set TRassembly path and cd into it
dir_TRassembly="/home/yshin/mendel-nas1/ToxCodAn-Genome/bin"
cd ${dir_TRassembly}

# output dir
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added"
mkdir -p ${outdir}
mkdir -p ${outdir}/SRR35908235_TRassembly

# input dir
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
venom_read_1="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq/SRR35908235_1.trimmed.fastq.gz"
venom_read_2="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq/SRR35908235_2.trimmed.fastq.gz"

# run TRassembly.py to assemble venom gland RNA-seq reads
python TRassembly.py \
    -g ${genome} \
    -r ${venom_read_1},${venom_read_2} \
    -o ${outdir}/SRR35908235_TRassembly \
    -c ${SLURM_CPUS_PER_TASK} \
    -M 250G

# echo the following message when done
echo "TRassembly completed successfully."
echo "The output files are located in ${outdir}/SRR35908235_TRassembly"
ls ${outdir}/SRR35908235_TRassembly
```

Next, run toxcodan.py, which is the main script of ToxCodAn (*NOT* the same as ToxCodAn-Genome!!), for transcriptome annotation:
```sh
#!/bin/bash
#SBATCH --job-name=ToxCodAn
#SBATCH --nodes=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env for mapping
source ~/.bash_profile
conda activate Toxcodan

# force UTF-8 output in non-interactive SLURM jobs
# locale available on all standard Linux nodes
export LANG=C
export LC_ALL=C

# force Python itself to use UTF-8
export PYTHONIOENCODING=UTF-8
export PYTHONUTF8=1

set -euo pipefail

# set toxcodan and signalp paths
dir_toxcodan="/home/yshin/mendel-nas1/ToxCodAn/bin"
signalp_dir="/home/yshin/mendel-nas1/signalp-4.1"

# make toxcodan and signalp available in the PATH
export PATH="${dir_toxcodan}:${signalp_dir}:${PATH}"
hash -r

# set input paths
models_dir="/home/yshin/mendel-nas1/ToxCodAn/models"

# output dir for all toxin gene annotation steps
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added"
mkdir -p ${outdir}

# output dir for ToxCodAn
mkdir -p ${outdir}/SRR35908235_ToxCodAn

# run ToxCodAn to annotate the transcriptome
python "${dir_toxcodan}/toxcodan.py" \
    -s SRR35908235 \
    -t ${outdir}/SRR35908235_TRassembly/transcripts.fasta \
    -m ${models_dir} \
    -o ${outdir}/SRR35908235_ToxCodAn \
    -c ${SLURM_CPUS_PER_TASK}

cat ${outdir}/SRR35908235_ToxCodAn/SRR35908235_Toxins_cds_RedundancyFiltered.fasta ${outdir}/SRR35908235_ToxCodAn/SRR35908235_PutativeToxins_cds_SPfiltered.fasta > ${outdir}/G_ussuriensis_VG_toxins.toxcodan.fasta

echo "ToxCodAn completed successfully."
echo "The output files are located in ${outdir}/SRR35908235_ToxCodAn"
ls ${outdir}/SRR35908235_ToxCodAn       
```

Finally, run toxcodan-genome.py (the main script of ToxCodAn-Genome) for toxin gene annotation of the genome:
```sh
#!/bin/bash
#SBATCH --job-name=ToxCodAn-Genome
#SBATCH --nodes=1
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=300:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env for mapping
source ~/.bash_profile
conda activate ToxcodanGenome

# force UTF-8 output in non-interactive SLURM jobs
# locale available on all standard Linux nodes
export LANG=C
export LC_ALL=C

# force Python itself to use UTF-8
export PYTHONIOENCODING=UTF-8
export PYTHONUTF8=1

set -euo pipefail

# set toxcodan-genome path
dir_toxcodan_genome="/home/yshin/mendel-nas1/ToxCodAn-Genome/bin"

# allow ToxCodAn-Genome to find helper scripts
export PATH="${dir_toxcodan_genome}:${PATH}"

# set input paths
genome="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
db_dir="/home/yshin/mendel-nas1/ToxCodAn-Genome/Databases/Viperidae_db_appended.fasta"

# output dir for all toxin gene annotation steps
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added"
mkdir -p ${outdir}

# output dir for ToxCodAn-Genome
mkdir -p ${outdir}/SRR35908235_ToxCodAn-Genome

# run ToxCodAn-Genome
python ${dir_toxcodan_genome}/toxcodan-genome.py \
    -g ${genome} \
    -d ${db_dir} \
    -C ${outdir}/G_ussuriensis_VG_toxins.toxcodan.fasta \
    -o ${outdir}/SRR35908235_ToxCodAn-Genome \
    -c ${SLURM_CPUS_PER_TASK}

echo "ToxCodAn-Genome completed successfully."
echo "The output files are located in ${outdir}/SRR35908235_ToxCodAn-Genome"
ls ${outdir}/SRR35908235_ToxCodAn-Genome
```

After this finishes running, convert the output annotation .gtf file to a .tsv file using "fromCDStoGENE.py" script.
```sh
# from "/home/yshin/mendel-nas1/ToxCodAn-Genome/bin" dir
python3 fromCDStoGENE.py /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/toxin_annotation.gtf  /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/toxin_annotation.tsv
```  
After the .tsv is created, plot annotation using "PlotToxinLoci.R" script. First, install R and required packages on the conda env.
```sh
conda install conda-forge::r-base
conda install -c conda-forge r-tidyverse r-ggplot2 r-gggenes r-ggrepel r-argparse
```
Then, run "PlotToxinLoci.R"
```sh
# make a plot dir
mkdir -p \
  /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/plots

# cd to toxcodan genome dir
cd /home/yshin/mendel-nas1/ToxCodAn-Genome/bin/

# run plotting script
Rscript \
  PlotToxinLoci.R \
  -i /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/toxin_annotation.tsv \
  -o /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome/plots/Gloydius_ussuriensis_toxin_annotation
```

### (Post-Hi-C) Finalizing non-toxin gene annotation

### (Post-Hi-C) Finalizing toxin gene annotation
```sh
#cd into output dir
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome

# warning / removed counts == 7 warnings and 48 removed
wc -l annotation_warning.txt
wc -l annotation_removed.txt

# GTF feature types
grep -v '^#' toxin_annotation.gtf | cut -f3 | sort | uniq -c

# CDS count
grep -c '^>' toxin_annotation_cds.fasta

# protein count
grep -c '^>' toxin_annotation_pep.fasta

# check TSV
head -20 toxin_annotation.tsv
```

Next, check the annotations that are removed and flagged with warnings:
```sh
cat annotation_warning.txt
cat annotation_removed.txt
```

The removed candidate gene models were:
```sh
- Ficolin-Partial-4
- SVSP-Partial-8
- CTL-Partial-7
- CTL-Partial-10
- SVMP-Partial-5
```

The annotations flagged for manual inspection were:
```sh
- SVSP-1
- SVSP-16
- Ficolin-4
- PLA2-4
- CTL-10
- HYAL-4
- CTL-16
```

So, there are five gene models what were removed, and seven that were retained with warnings. Now, let's inspect these flagged models.
```sh
awk -F '\t' '
$1=="SVSP-1" ||
$1=="SVSP-16" ||
$1=="Ficolin-4" ||
$1=="PLA2-4" ||
$1=="CTL-10" ||
$1=="HYAL-4" ||
$1=="CTL-16"
' toxin_annotation.tsv
```

The output will look like this:
```sh
|toxin_ID   |  toxin_family    |      chromosome   |     start  |       end   |  strand|
|-----------|------------------|-------------------|------------|-------------|--------|
|SVSP-1     |          SVSP    |   G_ussuri_chr1   |  16190589  |  16212840   |       -|
|SVSP-16    |          SVSP    |  G_ussuri_chr10   |  18308689  |  18312041   |       -|
|Ficolin-4  |       Ficolin    |  G_ussuri_chr15   |  10789984  |  10801121   |       -|
|PLA2-4     |          PLA2    |  G_ussuri_chr16   |    285580  |    286640   |       +|
|CTL-10     |           CTL    |   G_ussuri_chr6   |   1351446  |   1355701   |       +|
|HYAL-4     |          HYAL    |   G_ussuri_chr6   |  72034421  |  72039204   |       +|
|CTL-16     |           CTL    |   G_ussuri_chrZ   | 119509275  | 119512397   |       +|
```

Now, let's verify that all seven are present in both the CDS and protein FASTA outputs, then inspect their protein lengths.
```sh
echo "===== WARNED CDS HEADERS ====="

grep '^>' toxin_annotation_cds.fasta | \
grep -E 'SVSP-1($|[[:space:]])|SVSP-16($|[[:space:]])|Ficolin-4($|[[:space:]])|PLA2-4($|[[:space:]])|CTL-10($|[[:space:]])|HYAL-4($|[[:space:]])|CTL-16($|[[:space:]])'

echo
echo "===== WARNED PROTEIN HEADERS ====="

grep '^>' toxin_annotation_pep.fasta | \
grep -E 'SVSP-1($|[[:space:]])|SVSP-16($|[[:space:]])|Ficolin-4($|[[:space:]])|PLA2-4($|[[:space:]])|CTL-10($|[[:space:]])|HYAL-4($|[[:space:]])|CTL-16($|[[:space:]])'
```

Then get the protein lengths for the full toxin set, with the seven warning models labeled:
```sh
awk '
BEGIN {
    warn["SVSP-1"]=1
    warn["SVSP-16"]=1
    warn["Ficolin-4"]=1
    warn["PLA2-4"]=1
    warn["CTL-10"]=1
    warn["HYAL-4"]=1
    warn["CTL-16"]=1
}
/^>/ {
    if (id != "") {
        status = (id in warn ? "WARNING" : "")
        print id, length(seq), status
    }

    id = substr($1, 2)
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    if (id != "") {
        status = (id in warn ? "WARNING" : "")
        print id, length(seq), status
    }
}
' OFS='\t' toxin_annotation_pep.fasta | column -t
```

The outputs of above scripts will show that all seven warned models are retained consistently in both the CDS and protein FASTA outputs, and that the protein lengths by themselves do not make any of the seven look obviously truncated. Therefore, retain all seven for now and move to sequence-integrity QC. We need to check whether their CDSs have proper starts/stops and whether the translated proteins contain internal stop codons.

Run this first:
```sh
awk '
BEGIN {
    warn["SVSP-1"]=1
    warn["SVSP-16"]=1
    warn["Ficolin-4"]=1
    warn["PLA2-4"]=1
    warn["CTL-10"]=1
    warn["HYAL-4"]=1
    warn["CTL-16"]=1
}
/^>/ {
    if (id != "" && id in warn) {
        start = substr(seq,1,3)
        stop = substr(seq,length(seq)-2,3)
        print id, length(seq), length(seq)%3, start, stop
    }
    id = substr($1,2)
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    if (id != "" && id in warn) {
        start = substr(seq,1,3)
        stop = substr(seq,length(seq)-2,3)
        print id, length(seq), length(seq)%3, start, stop
    }
}
' OFS='\t' toxin_annotation_cds.fasta | column -t
```

The output will look like this:
```sh
gene        CDS_length    length_mod_3    first_codon    last_codon
SVSP-1             798               0            ATG           tag
SVSP-16            810               0            atg           tga
Ficolin-4         1113               0            ATG           TAG
PLA2-4             417               0            ATG           TAA
CTL-10             531               0            ATG           TGA
HYAL-4            1626               0            ATG           TGA
CTL-16             585               0            ATG           TAA
```

Also run:
```sh
awk '
BEGIN {
    warn["SVSP-1"]=1
    warn["SVSP-16"]=1
    warn["Ficolin-4"]=1
    warn["PLA2-4"]=1
    warn["CTL-10"]=1
    warn["HYAL-4"]=1
    warn["CTL-16"]=1
}
/^>/ {
    if (id != "" && id in warn) {
        tmp = seq
        nstop = gsub(/\*/, "", tmp)

        tmp = seq
        nx = gsub(/[Xx]/, "", tmp)

        first = substr(seq,1,1)
        last = substr(seq,length(seq),1)

        print id, length(seq), first, last, nstop, nx
    }
    id = substr($1,2)
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    if (id != "" && id in warn) {
        tmp = seq
        nstop = gsub(/\*/, "", tmp)

        tmp = seq
        nx = gsub(/[Xx]/, "", tmp)

        first = substr(seq,1,1)
        last = substr(seq,length(seq),1)

        print id, length(seq), first, last, nstop, nx
    }
}
' OFS='\t' toxin_annotation_pep.fasta | column -t
```

The output will look like this:
```sh
gene        protein_length    first_AA    last_AA    number_of_*    number_of_X
SVSP-1                 266           M          *              2              0
SVSP-16                270           M          *              7              0
Ficolin-4              371           M          *              3              0
PLA2-4                 139           M          *              2              0
CTL-10                 177           M          *              2              0
HYAL-4                 542           M          *              2              0
CTL-16                 195           M          *              3              0
```

This reveals some important clues:
 1) CDS length is divisible by 3 for every model.
 2) All begin with ATG; case differences are irrelevant.
 3) All end with a canonical stop codon.
 4) No proteins contain X.
 5) BUT, every warned protein contains at least one premature/internal stop codon (*) in addition to the normal terminal stop.

SVSP-16 is immediately the most suspicious because it has six internal stops. The others may still represent pseudogenes, truncated paralogs, or gene-model boundary issues.

Next, let's locate the exact amino-acid positions of those internal stops.

Run this:
```sh
awk '
BEGIN {
    warn["SVSP-1"]=1
    warn["SVSP-16"]=1
    warn["Ficolin-4"]=1
    warn["PLA2-4"]=1
    warn["CTL-10"]=1
    warn["HYAL-4"]=1
    warn["CTL-16"]=1
}
function report() {
    if (id != "" && id in warn) {
        for (i=1; i<=length(seq); i++) {
            if (substr(seq,i,1)=="*") {
                status = (i==length(seq) ? "terminal" : "INTERNAL")
                print id, i, status
            }
        }
    }
}
/^>/ {
    report()
    id=substr($1,2)
    seq=""
    next
}
{
    seq=seq $0
}
END {
    report()
}
' OFS='\t' toxin_annotation_pep.fasta | column -t
```

The output will look like this:
```sh
SVSP-1     265  INTERNAL
SVSP-1     266  terminal
SVSP-16    56   INTERNAL
SVSP-16    59   INTERNAL
SVSP-16    68   INTERNAL
SVSP-16    70   INTERNAL
SVSP-16    72   INTERNAL
SVSP-16    93   INTERNAL
SVSP-16    270  terminal
Ficolin-4  334  INTERNAL
Ficolin-4  349  INTERNAL
Ficolin-4  371  terminal
PLA2-4     113  INTERNAL
PLA2-4     139  terminal
CTL-10     143  INTERNAL
CTL-10     177  terminal
HYAL-4     443  INTERNAL
HYAL-4     542  terminal
CTL-16     75   INTERNAL
CTL-16     121  INTERNAL
CTL-16     195  terminal
```

The output suggests that SVSP-1 is the least alarming: its internal stop is at aa 265 of 266, immediately before the terminal stop. Ficolin-4, PLA2-4, CTL-10, and HYAL-4 also have relatively late disruptions. CTL-16 is more concerning because its two stops are at aa 75 and 121, and SVSP-16 is the strongest pseudogene/error candidate because it has six clustered stops between aa 56–93.

The next step is to map these premature stops back onto the genomic CDS exons. That will tell us whether they sit inside ordinary exons or near problematic splice boundaries.

Run:
```sh
python3 - <<'PY'
import re
from collections import defaultdict

warned = {
    "SVSP-1",
    "SVSP-16",
    "Ficolin-4",
    "PLA2-4",
    "CTL-10",
    "HYAL-4",
    "CTL-16"
}

# read CDS FASTA
cds = {}
current = None

with open("toxin_annotation_cds.fasta") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current = line[1:].split()[0]
            cds[current] = []
        elif current is not None:
            cds[current].append(line)

cds = {k: "".join(v).upper() for k, v in cds.items()}

# read protein FASTA
pep = {}
current = None

with open("toxin_annotation_pep.fasta") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            current = line[1:].split()[0]
            pep[current] = []
        elif current is not None:
            pep[current].append(line)

pep = {k: "".join(v) for k, v in pep.items()}

# read GTF CDS features
features = defaultdict(list)

with open("toxin_annotation.gtf") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip().split("\t")

        if len(fields) < 9 or fields[2] != "CDS":
            continue

        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        phase = fields[7]
        attributes = fields[8]

        m = re.search(r'gene_id "([^"]+)"', attributes)
        if not m:
            continue

        gene = m.group(1)

        if gene not in warned:
            continue

        features[gene].append({
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "phase": phase
        })

print(
    "gene",
    "aa_position",
    "protein_length",
    "codon",
    "chromosome",
    "genomic_positions",
    "exon_rank",
    "strand",
    sep="\t"
)

for gene in sorted(warned):

    sequence = cds[gene]
    protein = pep[gene]
    segs = features[gene]

    strand = segs[0]["strand"]

    if strand == "+":
        segs = sorted(segs, key=lambda x: x["start"])
    else:
        segs = sorted(segs, key=lambda x: x["start"], reverse=True)

    nt_to_genome = {}
    nt_to_exon = {}

    cds_pos = 1

    for exon_rank, seg in enumerate(segs, start=1):

        if strand == "+":
            coords = range(seg["start"], seg["end"] + 1)
        else:
            coords = range(seg["end"], seg["start"] - 1, -1)

        for genomic_pos in coords:
            nt_to_genome[cds_pos] = genomic_pos
            nt_to_exon[cds_pos] = exon_rank
            cds_pos += 1

    internal_stops = [
        i + 1
        for i, aa in enumerate(protein)
        if aa == "*" and i + 1 < len(protein)
    ]

    for aa_pos in internal_stops:

        nt_start = (aa_pos - 1) * 3 + 1
        nt_positions = [nt_start, nt_start + 1, nt_start + 2]

        codon = sequence[nt_start - 1:nt_start + 2]

        genomic_positions = [
            str(nt_to_genome.get(x, "NA"))
            for x in nt_positions
        ]

        exon_ranks = sorted({
            nt_to_exon.get(x, "NA")
            for x in nt_positions
        }, key=str)

        print(
            gene,
            aa_pos,
            len(protein),
            codon,
            segs[0]["chrom"],
            ",".join(genomic_positions),
            ",".join(map(str, exon_ranks)),
            strand,
            sep="\t"
        )
PY
```

The output of the above script confirms that every premature stop codon is contained entirely within a single CDS exon; none is produced by a codon spanning an exon junction.

The patterns are:
- SVSP-1: one stop in exon 5, essentially at the very end of the protein.
- SVSP-16: all six stops cluster in exon 1.
- Ficolin-4: both stops are in exon 9.
- PLA2-4: one stop in exon 4.
- CTL-10: one stop in exon 4.
- HYAL-4: one stop in exon 4.
- CTL-16: both stops are in exon 3.

So SVSP-16 and CTL-16 remain the most suspicious structurally, whereas SVSP-1 may simply have a terminal-boundary issue.

Next, inspect the complete CDS/exon structures for these seven models. 

Run:
```sh
for gene in SVSP-1 SVSP-16 Ficolin-4 PLA2-4 CTL-10 HYAL-4 CTL-16
do
    echo
    echo "===== ${gene} ====="
    grep "gene_id \"${gene}\"" toxin_annotation.gtf
done
```

The output of this script will provide some provisional results: SVSP-1 has its stop essentially at the end of the final CDS exon, so a terminal-boundary issue is most likely. Ficolin-4, PLA2-4, and CTL-10 also have stops in their final CDS exon, but farther upstream. HYAL-4 is more concerning because its stop is in exon 4 of 5. CTL-16 is clearly more suspicious because both stops occur in exon 3 of 5. SVSP-16 is the strongest disruption candidate because all six premature stops occur early in exon 1.

Now we should bring in the evidence ToxCodAn used to construct these models. Rather than searching by final names, search matched_regions.gtf by genomic overlap:
```sh
awk -F '\t' '
($1=="G_ussuri_chr1"  && $4<=16212840  && $5>=16190589) ||
($1=="G_ussuri_chr10" && $4<=18312041  && $5>=18308689) ||
($1=="G_ussuri_chr15" && $4<=10801121  && $5>=10789984) ||
($1=="G_ussuri_chr16" && $4<=286640    && $5>=285580) ||
($1=="G_ussuri_chr6"  && $4<=1355701   && $5>=1351446) ||
($1=="G_ussuri_chr6"  && $4<=72039204  && $5>=72034421) ||
($1=="G_ussuri_chrZ"  && $4<=119512397 && $5>=119509275)
' matched_regions.gtf
```

Then, separately, identify the matched-GTF files corresponding to those loci:
```sh
for region in \
"G_ussuri_chr1--161905" \
"G_ussuri_chr10--183086" \
"G_ussuri_chr15--107899" \
"G_ussuri_chr16--285580" \
"G_ussuri_chr6--1351446" \
"G_ussuri_chr6--72034421" \
"G_ussuri_chrZ--119509275"
do
    echo
    echo "===== ${region} ====="
    ls matched_GTFs | grep "$region"
done
```

The output will show that three loci have clear transcript-derived support in matched_regions.gtf:
 - SVSP-1 has multiple independent transcript assemblies (TRINITY, STRG, and NODE) covering essentially the same   locus, from 16190592–16212840.   
 - SVSP-16 has extensive transcript evidence: Trinity, several StringTie (STRG) transcripts, and numerous NODE assemblies all map to this locus.
 - PLA2-4 has several transcript-derived matches, including a very high-coverage NODE_11840... hit from 285580–286562, plus another NODE, STRG, and Trinity transcript.

By contrast, the evidence shown for Ficolin-4, CTL-10, HYAL-4, and CTL-16 is toxin-database homology (Viperidae_Toxin_*) rather than transcript-derived evidence.

There are also two particularly interesting boundary clues.

For SVSP-1, the transcript-supported region begins at 16190592, exactly where its premature stop codon sits on the minus strand (16190594–16190592). The final ToxCodAn model extends three nucleotides farther to 16190589. That makes a terminal-boundary artifact quite plausible rather than a genuine disrupted gene.

For PLA2-4, the strongest transcript-derived hit ends at 286562, exactly where the premature TGA is located (286560–286562). Meanwhile, homology models extend the locus through 286640. Again, that strongly suggests we may be looking at competing terminal boundaries rather than automatically at a pseudogene.

SVSP-16 is different. It has abundant transcript evidence, yet six premature stops occur early in exon 1. So expression alone does not rescue it; this remains our strongest candidate for a transcribed pseudogene or erroneous CDS model.

Another thing to note, CTL-16 (which is based on database homology matches) being located on the Z chromosome is rather suspicious given toxin gene annotations in other snake species. So this one also deserves careful inspection. 

The next thing to do is to inspect the actual seven matched_GTFs files, because those tell us which model ToxCodAn selected as the best gene structure.

Run:
```sh
for f in \
matched_GTFs/G_ussuri_chr1--16190592-16212840.gtf \
matched_GTFs/G_ussuri_chr10--18308689-18313198.gtf \
matched_GTFs/G_ussuri_chr15--10789984-10801121.gtf \
matched_GTFs/G_ussuri_chr16--285580-286562.gtf \
matched_GTFs/G_ussuri_chr6--1351446-1355701.gtf \
matched_GTFs/G_ussuri_chr6--72034421-72039204.gtf \
matched_GTFs/G_ussuri_chrZ--119509275-119512397.gtf
do
    echo
    echo "===== $f ====="
    cat "$f"
done
```

The output can be summarized like so:
```txt
1) SVSP-1: Selected best hit = STRG.279.1; internal codon located at aa 265/266. This gene CDS is supported by RNA transcript-supported CDS, but penultimate codon is a stop. This annotation is mildly suspicious; could be a tiny terminal modeling issue.

2) SVSP-16: Selected best hit =	NODE_12577;	it is supported by transcript evidence but has 6 stops, spanning aa 56–93. This annotation is	highly suspicious; likely disrupted/transcribed pseudogene or bad model.

3) Ficolin-4: Selected best hit = Viperidae_Toxin_745; internal stop codon located at	aa 334, 349/371. Therefore, this annotation is a homology-only model with two late stops.

4) PLA2-4: Selected best hit = Viperidae_Toxin_1149; internal stop codon located at aa 113/139. This annotation is based on a homology model, but independent RNA evidence exists and one transcript ends at the internal stop.

5) CTL-10: Selected best hit = Viperidae_Toxin_515;	internal stop codon located at aa 143/177. This is a homology-only disrupted model.

6) HYAL-4: Selected best hit = Viperidae_Toxin_686; internal stop codon located at aa 443/542, exon 4/5. This is more concerning because substantial CDS follows the internal stop codon.

7) CTL-16: Selected best hit = Viperidae_Toxin_739;	internal stop codon located at aa 75, 121/195, exon 3/5. This is highly suspicious, especially given its location on the Z chromosome.
```

For SVSP-1, the selected RNA-derived model itself contains exactly the same five CDS exons as the final ToxCodAn annotation. So the stop at aa 265 is not something ToxCodAn added after choosing the RNA model. The transcript-supported model itself produces that stop.

SVSP-16 is even clearer. ToxCodAn selected a transcript-derived NODE model, and the final exon coordinates are identical. Yet six stop codons occur in its first coding exon. For now, we will put SVSP-16 in the pseudogene/problematic-model bucket, pending sequence validation.

For PLA2-4, the situation is more nuanced. The selected model is actually the database-derived Viperidae_Toxin_1149_PLA2, extending through 286640. But independently, we have an RNA-derived NODE match ending at 286562, exactly where the apparent internal TGA occurs. So PLA2-4 deserves specific RNA inspection before we classify it. There may genuinely be a shorter expressed transcript ending at that stop.

CTL-16 is looking increasingly questionable: database-only best hit, two internal stops in the middle of the gene, substantial coding sequence downstream, and its location on the Z chromosome. The Z location alone is not disqualifying, but combined with everything else it raises the possibility that this is a degraded/general C-type lectin-like locus that ToxCodAn matched to its toxin database rather than a functional venom CTL.

Before going to raw-read validation, The cleanest next check would be to inspect the genomic context. Let's see whether these warning genes sit with the other members of their toxin families or are genomic outliers. 

Run:
```sh
awk -F '\t' '
BEGIN {
    OFS="\t"
}
$2=="SVSP" ||
$2=="Ficolin" ||
$2=="PLA2" ||
$2=="CTL" ||
$2=="SNACLEC" ||
$2=="HYAL" {
    print $1,$2,$3,$4,$5,$6
}
' toxin_annotation.tsv | sort -k2,2 -k3,3V -k4,4n | column -t
```

The output will look like this:
```sh
> ' toxin_annotation.tsv | sort -k2,2 -k3,3V -k4,4n | column -t
CTL-1      CTL      G_ussuri_chr1   111265878  111272415  -
CTL-2      CTL      G_ussuri_chr5   95628854   95631475   -
CTL-3      CTL      G_ussuri_chr6   448197     450286     +
CTL-4      CTL      G_ussuri_chr6   867856     872443     +
CTL-5      CTL      G_ussuri_chr6   1016550    1022773    -
CTL-6      CTL      G_ussuri_chr6   1046375    1053967    +
CTL-7      CTL      G_ussuri_chr6   1292124    1325647    +
CTL-8      CTL      G_ussuri_chr6   1298499    1305085    +
CTL-9      CTL      G_ussuri_chr6   1332993    1343087    -
CTL-10     CTL      G_ussuri_chr6   1351446    1355701    +
CTL-11     CTL      G_ussuri_chr6   1365922    1371014    -
CTL-12     CTL      G_ussuri_chr6   1374967    1377991    +
CTL-13     CTL      G_ussuri_chr6   1416996    1422999    -
CTL-14     CTL      G_ussuri_chr6   1427637    1432276    +
CTL-15     CTL      G_ussuri_chrZ   53703611   53706064   -
CTL-16     CTL      G_ussuri_chrZ   119509275  119512397  +
CTL-17     CTL      G_ussuri_chrZ   119519092  119521350  -
CTL-18     CTL      G_ussuri_chrZ   133429644  133435453  -
Ficolin-1  Ficolin  G_ussuri_chr1   107463308  107476757  +
Ficolin-2  Ficolin  G_ussuri_chr1   107493750  107504245  +
Ficolin-3  Ficolin  G_ussuri_chr1   107517199  107526160  +
Ficolin-4  Ficolin  G_ussuri_chr15  10789984   10801121   -
HYAL-1     HYAL     G_ussuri_chr2   202028395  202036528  -
HYAL-2     HYAL     G_ussuri_chr6   71969588   71980263   +
HYAL-3     HYAL     G_ussuri_chr6   72004412   72016075   +
HYAL-4     HYAL     G_ussuri_chr6   72034421   72039204   +
PLA2-1     PLA2     G_ussuri_chr10  20681231   20685339   +
PLA2-2     PLA2     G_ussuri_chr16  278191     278997     +
PLA2-3     PLA2     G_ussuri_chr16  281098     282719     -
PLA2-4     PLA2     G_ussuri_chr16  285580     286640     +
PLA2-5     PLA2     G_ussuri_chr16  293721     295345     +
PLA2-6     PLA2     G_ussuri_chr16  302420     304043     +
PLA2-7     PLA2     G_ussuri_chr16  311102     312734     +
SNACLEC-1  SNACLEC  G_ussuri_chr6   1385725    1389251    -
SNACLEC-2  SNACLEC  G_ussuri_chr6   1404163    1409516    +
SVSP-1     SVSP     G_ussuri_chr1   16190589   16212840   -
SVSP-21    SVSP     G_ussuri_chr2   24809373   24816201   +
SVSP-22    SVSP     G_ussuri_chr2   137498143  137506778  +
SVSP-2     SVSP     G_ussuri_chr10  17846269   17853935   -
SVSP-3     SVSP     G_ussuri_chr10  17897988   17916218   -
SVSP-4     SVSP     G_ussuri_chr10  17933543   17943451   -
SVSP-5     SVSP     G_ussuri_chr10  17934046   17943451   -
SVSP-6     SVSP     G_ussuri_chr10  17988255   17993030   -
SVSP-7     SVSP     G_ussuri_chr10  18042288   18046724   -
SVSP-8     SVSP     G_ussuri_chr10  18091907   18096890   -
SVSP-9     SVSP     G_ussuri_chr10  18108802   18115307   -
SVSP-10    SVSP     G_ussuri_chr10  18138584   18149575   -
SVSP-11    SVSP     G_ussuri_chr10  18155360   18162289   -
SVSP-12    SVSP     G_ussuri_chr10  18170088   18175096   -
SVSP-13    SVSP     G_ussuri_chr10  18182232   18199625   -
SVSP-14    SVSP     G_ussuri_chr10  18222991   18230561   -
SVSP-15    SVSP     G_ussuri_chr10  18278654   18282773   -
SVSP-16    SVSP     G_ussuri_chr10  18308689   18312041   -
SVSP-17    SVSP     G_ussuri_chr10  18319874   18325066   -
SVSP-18    SVSP     G_ussuri_chr10  18341314   18354244   -
SVSP-19    SVSP     G_ussuri_chr10  18372679   18386109   -
SVSP-20    SVSP     G_ussuri_chr10  18401912   18415458   -
```
This shows that the location of CTL-16 on the Z chromosome is not inherently suspicious after all. There are four CTL annotations on Z. More importantly, CTL-16 and CTL-17 are only ~6.7 kb apart, so CTL-16 sits next to another CTL-family member. That strongly supports this being a genuine CTL-family locus rather than a random spurious homology hit. What remains suspicious is its functionality: CTL-16 alone has two internal stops at aa 75 and 121.

The same context is informative for the other warnings:
 - SVSP-16 sits right inside the large chr10 SVSP tandem cluster (SVSP-2 through SVSP-20). So its identity as an SVSP-family locus is very convincing; the six internal stops make it look more like a pseudogenized SVSP paralog than a false annotation.
 - PLA2-4 sits squarely inside the chr16 PLA2 cluster (PLA2-2 through PLA2-7). Again, family identity is strong; the issue is whether it is intact or pseudogenized.
 - HYAL-4 is immediately downstream of HYAL-2 and HYAL-3 on chr6, so its HYAL identity is also well supported.
 - CTL-10 is embedded inside the huge chr6 CTL/SNACLEC cluster, so it is almost certainly a real CTL-family locus.
 - Ficolin-4 is the odd one: Ficolin-1–3 form a cluster on chr1, whereas Ficolin-4 is isolated on chr15. That doesn't make it false, but it deserves a higher bar.
 - SVSP-1 is likewise isolated on chr1 while most SVSPs form the chr10 cluster.

So at this point, there are two separate questions:

__1) Is this really a member of the toxin-gene family?__
For SVSP-16, PLA2-4, CTL-10, HYAL-4, and probably CTL-16, genomic context says yes, very plausibly.

__2) Is the locus an intact functional toxin gene?__
This is still unresolved, and the internal stops are the problem.

For CTL-16 specifically, the next clean check should be to compare it directly with its nearby Z-linked paralog CTL-17.

Run:
```sh
for gene in CTL-15 CTL-16 CTL-17 CTL-18
do
    echo
    echo "===== ${gene} GTF ====="
    grep "gene_id \"${gene}\"" toxin_annotation.gtf

    echo
    echo "===== ${gene} PROTEIN ====="
    awk -v gene="$gene" '
    /^>/ {
        if (p) exit
        if ($1==">"gene) {
            p=1
            print
            next
        }
    }
    p {print}
    ' toxin_annotation_pep.fasta
done
```

This comparison makes CTL-16 substantially more suspicious as a toxin annotation, and for a reason beyond the two premature stops.

For example, CTL-17 looks exactly like a secreted venom C-type lectin/SNACLEC-like protein: it is ~155 aa in length and begins with a strongly hydrophobic N-terminal signal peptide:
```txt
MRPLAYSSLWLLGFFMAVSQAAF...
```
Then it transitions directly into a cysteine-rich lectin-like mature protein.

CTL-16 looks quite different; the conspicuous hydrophobic stretch is not right at the N-terminus; it appears around aa ~50–70. This looks much more like a transmembrane helix than a conventional secretory signal peptide. 

Combined with:
- 195 aa length, longer than the nearby intact CTL-17;
- two internal stop codons at aa 75 and 121;
- five-exon structure;
- database-only rather than transcript-derived ToxCodAn support;
- location next to another CTL-like gene on the Z;

The most likely conclusion is that the annotated CTL-16 is a genuine C-type lectin-domain-family locus, but probably not a functional secreted venom CTL. It could be a degraded membrane-associated CTL-like paralog/pseudogene that happens to match the toxin database.

Another interesting pattern is that CTL-15 and CTL-18 also don't immediately resemble canonical secreted venom CTLs at their N-termini, whereas CTL-17 does. So the Z may contain a broader CTL-domain gene family, only some members of which should necessarily be considered venom toxins.

The next useful check is very focused: inspect whether the autosomal chr6 CTLs have the expected short N-terminal signal peptides. Let's print all CTL proteins in one compact form, first 40 aa plus total length:
```sh
awk '
function print_record() {
    if (id ~ /^CTL-/) {
        print id, length(seq), substr(seq,1,40)
    }
}

/^>/ {
    print_record()

    id = substr($1,2)
    seq = ""
    next
}

{
    seq = seq $0
}

END {
    print_record()
}
' OFS='\t' toxin_annotation_pep.fasta | column -t
```

The output looks like this:
```sh
CTL-1   178  MCSWLFICAFLGITFMQDVSAQTCLCAQGFCASGWVQYKS
CTL-2   167  MRLILLSSLLLLGTFLANGHETDPDGHVLNSLMLTVLRLQ
CTL-3   168  MASFSLGLFLGLVLAPLAANAQADPCPSGWLAYNGHCYGY
CTL-4   165  MFLITCFIFGVLGSLTWAGPSVRTVCPSGTFVYKDGSQWY
CTL-5   167  MFLITCFIFGVLGSLTWAGPSVRTVCPSGTFAYKDRNQWY
CTL-6   159  MGRLIFMSFGLLVVFLSLSGVKGCCCPRDWLPMNGLCYKV
CTL-7   147  MGRFIFVSFGLLVVFLSLGGTGADCPSGWSSFEGHCYRVF
CTL-8   139  MGRFIFMSFGLLVVFLSLRIGGNHHCPFGWFFHNASCYKL
CTL-9   101  MGRFIFVSFGLLVLFLSLSGTAADCPSGWSSYEGHCYQPF
CTL-10  177  MGRFIFVSFGLLVVFLSLSGTGADCPSDWSSYEGHCYKTF
CTL-11  159  MGRFIFVSFSLLVVFLSLSGAEDDFDCPSDWPSNGRFCYK
CTL-12  101  MGRLIFVSFGLLVVFLSLSGTGADCPSDWSSYEGHCYRIF
CTL-13  159  MGRFIFVSFGLLVLFLSLSGTGADFDCPPGWSAYDRYCYQ
CTL-14  151  MGRFIFVSFGLLVVAISLSGTGAGFCCPMGWSSYDQHCYR
CTL-15  161  MIEELKEQLDNLSQEVALLKEKQALQTVCLKGTKINLKCV
CTL-16  195  MLLHICLRLYKATAQKVREYLRGIQLDSPPLLSIQKEKTI
CTL-17  155  MRPLAYSSLWLLGFFMAVSQAAFCPREWLRYQRNCYGLFH
CTL-18  165  MKLTQAQLETMEQNSRALNCEIIEMKSNGSKSGCCPRGWL
```
This comparison is quite revealing. The chr6 CTLs (CTL-3 through CTL-14) are extremely consistent: nearly all begin with a strongly hydrophobic N-terminal signal peptide. CTL-17 on Z also fits this pattern. But CTL-15, CTL-16, and CTL-18 do not match this pattern.

For CTL-16, this is especially problematic because it additionally has two internal stops and a hydrophobic stretch appearing much later in the sequence. So we would now treat this annotation as very unlikely to represent an intact secreted venom CTL. More likely it is a CTL-like paralog/pseudogene or a non-toxin CTL-family locus that ToxCodAn picked up based database homology.

But this also uncovers something broader: CTL-15 and CTL-18 deserve review too, even though ToxCodAn did not flag them. The automated warning set is therefore not sufficient by itself for final curation.

The next useful step is to inspect whether the stop codons in the warned models are genuinely encoded in the genome or could reflect assembly/model errors. Let's check the reference bases first at the exact stop positions.

Run:
```sh
samtools faidx \
/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta \
G_ussuri_chr6:1355590-1355605 \
G_ussuri_chrZ:119511508-119511662 \
G_ussuri_chr15:10790040-10790105 \
G_ussuri_chr6:72038492-72038508 \
G_ussuri_chr16:286553-286569 \
G_ussuri_chr1:16190585-16190600 \
G_ussuri_chr10:18311755-18311885 \
G_ussuri_chrZ:53705975-53706070 \
G_ussuri_chrZ:133435360-133435460
```

This confirms that the premature stops are genuinely present in the reference assembly sequence; they are not artifacts introduced during translation or FASTA extraction. Also, for CTL-15 and CTL-18, their unusual N-termini are genuinely encoded by the reference assembly and are not caused by the protein FASTA extraction.

Now, we still need to distinguish two possibilities: real biological disruption/pseudogenization versus an assembly/base-call error at the reference locus. Given our HiFi assembly quality, widespread assembly errors are unlikely, especially for models like SVSP-16 with six separate stops, but checking the underlying genomic reads is worthwhile before we label genes as pseudogenes.

At this stage, this is the provisional status for different toxin gene annotations:
  1) SVSP-1:   genuine SVSP locus; one very late internal stop; needs read validation
  2) SVSP-16:  genuine SVSP locus; strong pseudogene candidate
  3) Ficolin-4: homology-only isolated locus; questionable functional gene
  4) PLA2-4:   genuine PLA2-cluster locus; one internal stop; needs read/RNA evaluation
  5) CTL-10:   genuine CTL-cluster locus; likely disrupted if stop confirmed
  6) HYAL-4:   genuine HYAL-cluster locus; likely disrupted if stop confirmed
  7) CTL-15:   intact, transcript-supported CTL-like gene; no obvious signal peptide; toxin status questionable
  8) CTL-16:   genuine CTL-like locus; strong nonfunctional/pseudogene candidate; two internal stops; no obvious  
               signal peptide
  9) CTL-17:   intact; no internal stops; conventional secreted CTL-like architecture
 10) CTL-18:   intact, transcript-supported CTL-like gene; no obvious signal peptide; toxin status questionable

Now, let's move on to genomic read validation. The purpose here is very specific: determine whether each premature stop codon is genuinely supported by the underlying genomic reads or whether it could reflect an assembly/base-calling error. CTL-15, CTL-17, and CTL-18 do not need stop-codon validation because they have no internal stops.

First, create a BED file containing all 14 internal stop codons:
```sh
printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
'G_ussuri_chr1'  '16190591'  '16190594'  'SVSP-1_stop265'      '0' '-' \
'G_ussuri_chr10' '18311873'  '18311876'  'SVSP-16_stop56'     '0' '-' \
'G_ussuri_chr10' '18311864'  '18311867'  'SVSP-16_stop59'     '0' '-' \
'G_ussuri_chr10' '18311837'  '18311840'  'SVSP-16_stop68'     '0' '-' \
'G_ussuri_chr10' '18311831'  '18311834'  'SVSP-16_stop70'     '0' '-' \
'G_ussuri_chr10' '18311825'  '18311828'  'SVSP-16_stop72'     '0' '-' \
'G_ussuri_chr10' '18311762'  '18311765'  'SVSP-16_stop93'     '0' '-' \
'G_ussuri_chr15' '10790094'  '10790097'  'Ficolin-4_stop334'  '0' '-' \
'G_ussuri_chr15' '10790049'  '10790052'  'Ficolin-4_stop349'  '0' '-' \
'G_ussuri_chr16' '286559'    '286562'    'PLA2-4_stop113'     '0' '+' \
'G_ussuri_chr6'  '1355596'   '1355599'   'CTL-10_stop143'     '0' '+' \
'G_ussuri_chr6'  '72038498'  '72038501'  'HYAL-4_stop443'     '0' '+' \
'G_ussuri_chrZ'  '119511514' '119511517' 'CTL-16_stop75'      '0' '+' \
'G_ussuri_chrZ'  '119511652' '119511655' 'CTL-16_stop121'     '0' '+' \
> toxin_internal_stops.bed
```

Check:
```sh
column -t toxin_internal_stops.bed
```

Next, locate the genomic-read alignment we can use to validate internal-stop codons. There is a AMNH_21010 illumina WGS BAM file under the population genomics project separate from the genome assembly project directory:
```sh
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

BAM="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam/AMNH_21010.markdup.bam"
```

Confirm that the BAM is intact, indexed, and uses the same chromosome naming scheme:
```sh
samtools quickcheck -v "$BAM"

ls -lh "$BAM" "$BAM".bai 2>/dev/null

samtools idxstats "$BAM" | head -20
```

Now, let's generate a per-base pileup of the Illumina reads specifically at the 14 internal stop codons we identified:
```sh
samtools mpileup \
    -f "$REF" \
    -l toxin_internal_stops.bed \
    -q 20 \
    -Q 20 \
    -aa \
    "$BAM" \
    > toxin_internal_stops_AMNH_21010.mpileup
```

Inspect the pileup:
```sh
ls -lh toxin_internal_stops_AMNH_21010.mpileup

wc -l toxin_internal_stops_AMNH_21010.mpileup
column -t toxin_internal_stops_AMNH_21010.mpileup
```

The result of the "column -t" line should look something like this:
```sh
G_ussuri_chr1   16190592   t  19  ..,.,,..,..,.......       JJJJJJJJJJJJJJJJJJJ
G_ussuri_chr1   16190593   c  19  ..,.,,..,..,.......       JJJJJJJJJJJJJJJJJJJ
G_ussuri_chr1   16190594   A  20  ..,.,,..,..,.......^].    JJJJJJJJJJJJJJJJJJJB
G_ussuri_chr6   1355597    T  9   ,.,..,.,.                 JsJJJJJJJ
G_ussuri_chr6   1355598    A  9   ,.,..,.,.                 JsJJJJJJJ
G_ussuri_chr6   1355599    G  9   ,$.,..,.,.                EsJJJJJJJ
G_ussuri_chr6   72038499   T  11  .,..,,,,...               sJssJJJJJJJ
G_ussuri_chr6   72038500   G  12  .,..,,,,..,.              sJssJJJJJJ9J
G_ussuri_chr6   72038501   A  12  .$,..,,,,..,.             nJssJJJJJJ9J
G_ussuri_chrZ   119511515  T  15  ,,..,,...,,,..,           JJJ9JJJJJJsJJJJ
G_ussuri_chrZ   119511516  G  15  ,,..,,...,,,..,           JJJJJJJJJJsJJJJ
G_ussuri_chrZ   119511517  A  15  ,,..,,...,,,..,           JJJJJJJJJJsJJJJ
G_ussuri_chrZ   119511653  T  15  ,,.,.,..,,..,,,           JJJJJ9JJJJJJJJJ
G_ussuri_chrZ   119511654  A  15  ,,.,.,..,,..,,,           JJJJJJJ9JJJJJJJ
G_ussuri_chrZ   119511655  A  15  ,,.,.,..,,..,,,           JJJJJJJJJJJJJJJ
```
Here, . and , mean the Illumina base agrees with the reference. The ^] and $ symbols are just read-start/read-end markers and are not mismatches. The genomic read validation basically tells us that the Illumina reads agree with the reference sequence for essentially every stop position. Therefore, these premature stops are not reference assembly base-call errors.

But Ficolin-4 has one interesting exception. At the second Ficolin-4 stop codon (chr15:10790050), the reference base is a T, read depth for that position is 16x. Five of those reads match the reference base, while the other 11 have an alternative base, C (this is printed as .Ccc,,c.Cc.cCCCc). Because Ficolin-4 is on the minus strand, the reference allele produces a canonical stop (TGA), whereas the genomic T→C alternative converts the coding-strand codon to TGG, which is a Trp. So this appears to be a heterozygous stop-loss polymorphism at aa349. 

However, this does not rescue Ficolin-4, because its earlier stop at aa334 is essentially completely supported by ~20× reads. Thus its current CDS still contains a genuine premature stop regardless of the aa349 allele.

__At this stage, the assessment for flagged annotations can be summarized like this:__
1) SVSP-1:
   genuine SVSP locus; very late premature stop at aa265;
   stop strongly supported by genomic Illumina reads;
   assembly error ruled out; functional status still requires RNA/model evaluation

2) SVSP-16:
   genuine SVSP locus;
   six premature stops at aa56–93;
   all six strongly supported by genomic reads;
   very strong pseudogene/disrupted paralog candidate

3) Ficolin-4:
   isolated homology-derived Ficolin locus;
   premature stop at aa334 strongly supported;
   second stop at aa349 appears heterozygous stop/TGG;
   strong disrupted/nonfunctional candidate under current model

4) PLA2-4:
   genuine PLA2-cluster locus;
   premature TGA at aa113 strongly supported by genomic reads;
   requires RNA/model evaluation because transcript evidence exists

5) CTL-10:
   genuine CTL-cluster locus;
   premature TAG at aa143 strongly supported by genomic reads;
   likely disrupted under current model

6) HYAL-4:
   genuine HYAL-cluster locus;
   premature TGA at aa443 strongly supported by genomic reads;
   likely disrupted under current model

7) CTL-15:
   intact transcript-supported CTL-like gene;
   no internal stops;
   no obvious N-terminal signal peptide;
   toxin status uncertain

8) CTL-16:
   genuine CTL-like locus;
   two premature stops strongly supported by genomic reads;
   no obvious N-terminal signal peptide;
   very strong nonfunctional/pseudogene candidate

9) CTL-17:
   intact;
   no internal stops;
   conventional secreted CTL-like architecture;
   intact, potentially secreted CTL-like gene; venom-toxin identity unconfirmed.

10) CTL-18:
    intact transcript-supported CTL-like gene;
    no internal stops;
    no obvious N-terminal signal peptide;
    toxin status uncertain 

Now, let's move on to RNA read validation. For each of the 10 loci, we'll assess:
 1) RNA expression: Are reads mapped to the gene, particularly in venom-gland RNA-seq?
 2) Exon structure: Do splice-junction reads support the predicted exon boundaries?
 3) Alternative transcripts: Is there an alternative isoform that avoids the premature stop codon?
 4) Stop-codon support: Do RNA reads support the premature stop or a different allele?
 5) Toxin identity: For CTL-15, CTL-17, and CTL-18, is there evidence for venom-gland expression rather than merely general transcription?

We'll prioritize venom-gland datasets SRR35908235 and SRR35908238, but the six tissue RNA-seq datasets can provide useful comparative expression evidence. There is one important limitation in our sampling design: because the venom gland RNA-seq data were from different specimens than the reference individual, they may have different alleles from the reference individual. Therefore, a non-stop allele in RNA would not necessarily indicate an assembly error.

First, let's create a BED file containing all 10 loci. This can be generated directly from toxin_annotation.tsv to avoid introducing coordinate errors. Run the script below from the "SRR35908235_ToxCodAn-Genome" directory.
```sh
python3 - <<'PY'
targets = [
    "SVSP-1",
    "SVSP-16",
    "Ficolin-4",
    "PLA2-4",
    "CTL-10",
    "HYAL-4",
    "CTL-15",
    "CTL-16",
    "CTL-17",
    "CTL-18"
]

records = {}

with open("toxin_annotation.tsv") as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if fields[0] in targets:
            records[fields[0]] = fields

missing = set(targets) - set(records)

if missing:
    raise ValueError(f"Missing genes: {sorted(missing)}")

with open("toxin_RNA_validation_targets.bed", "w") as out:
    for gene in targets:
        fields = records[gene]

        chrom = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[5]

        # BED uses 0-based, half-open coordinates
        out.write(
            f"{chrom}\t{start-1}\t{end}\t{gene}\t0\t{strand}\n"
        )

print(f"Created BED containing {len(targets)} loci.")
PY
```

Verify:
```sh
column -t toxin_RNA_validation_targets.bed
```

The venom gland RNA-seq fastq files are stored under this directory:
```sh
/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq
```

Next, let's check the location of RNA alignments generated in previous steps:
```sh
find /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation \
    -type f \
    \( -name "*.bam" -o -name "*.cram" -o -name "*.sam" \) \
    | sort
```

Let's verify the existing RNA alignments:
```sh
RNA_BAM="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.sorted.bam"

conda activate funannotate

echo "===== FILES ====="
ls -lh "$RNA_BAM" "${RNA_BAM}.bai" "${RNA_BAM%.bam}.bai" 2>/dev/null

echo
echo "===== BAM INTEGRITY ====="
samtools quickcheck -v "$RNA_BAM" && echo "PASS"

echo
echo "===== ALIGNMENT HEADER ====="
samtools view -H "$RNA_BAM" \
    | grep -E '^@HD|^@RG|^@PG' \
    | cut -c1-250 \
    | head -20

echo
echo "===== REFERENCE CHROMOSOMES ====="
samtools view -H "$RNA_BAM" \
    | grep '^@SQ' \
    | head -20
```
This will show us that we have a usable, indexed, splice-aware RNA-seq BAM.

Next, verify the original input reads. The script below will display the complete original HISAT2 command, but it will also show that HISAT2 recorded its inputs as temporary pipes (-1 /tmp/2326903.inpipe1 -2 /tmp/2326903.inpipe2).
```sh
samtools view -H "$RNA_BAM" \
    | awk -F '\t' '
        $1=="@PG" && $2=="ID:hisat2" {
            for (i=1; i<=NF; i++) {
                if ($i ~ /^CL:/) {
                    print substr($i,4)
                }
            }
        }
    ' \
    | fold -s -w 120
```

Next, let's check the alignment statistics and read identifiers:
```sh
echo "===== HISAT2 SUMMARY ====="
cat /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.hisat2_summary.txt

echo
echo "===== BAM READ IDENTIFIERS ====="
samtools view "$RNA_BAM" | head -n 5 | cut -f1

echo
echo "===== SRR35908235 FASTQ IDENTIFIER ====="
zcat /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq/SRR35908235_1.trimmed.fastq.gz | head -n 1

echo
echo "===== SRR35908238 FASTQ IDENTIFIER ====="
zcat /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/venom_gland/trimmed_fastq/SRR35908238_1.trimmed.fastq.gz | head -n 1
```

The `samtools view "$RNA_BAM" | head -n 5 | cut -f1` line will print out:
```sh
SRR35908235.1905510
SRR35908235.4122839
SRR35908235.1905510
SRR35908235.4122839
SRR35908235.9683651
```
This means that the SRA accession SRR35908235 was successfully aligned to the reference genome. The overall alignment level is 87.69%, which means that the existing BAM is suitable for RNA read-based toxin gene validation. One important consideration is the relatively high multimapping rate (22.38%). This is likely because SVSP, PLA2, and CTL genes form tandem repeats.

Now, let's check RNA read support across all 10 loci. First, let's quantify the number of RNA alignments overlapping each candidate gene.
```sh
# ============================================================
# Preliminary venom-gland RNA support
# Dataset: SRR35908235
# Targets: 10 manually reviewed toxin loci
# ============================================================

RNA_BAM="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.sorted.bam"

OUT="toxin_RNA_validation_read_counts.tsv"

# Output header
printf 'Gene\tChromosome\tStart\tEnd\tPrimary_alignments\tMAPQ20_alignments\n' > "$OUT"

# Process each locus
while IFS='|' read -r gene chrom start end; do

    region="${chrom}:${start}-${end}"

    # Count primary mapped alignments
    total=$(samtools view \
        -c \
        -F 0x904 \
        "$RNA_BAM" \
        "$region") || break

    # Count primary mapped alignments with MAPQ >= 20
    highconf=$(samtools view \
        -c \
        -F 0x904 \
        -q 20 \
        "$RNA_BAM" \
        "$region") || break

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$gene" "$chrom" "$start" "$end" "$total" "$highconf" \
        >> "$OUT"

done <<'EOF'
SVSP-1|G_ussuri_chr1|16190589|16212840
SVSP-16|G_ussuri_chr10|18308689|18312041
Ficolin-4|G_ussuri_chr15|10789984|10801121
PLA2-4|G_ussuri_chr16|285580|286640
CTL-10|G_ussuri_chr6|1351446|1355701
HYAL-4|G_ussuri_chr6|72034421|72039204
CTL-15|G_ussuri_chrZ|53703611|53706064
CTL-16|G_ussuri_chrZ|119509275|119512397
CTL-17|G_ussuri_chrZ|119519092|119521350
CTL-18|G_ussuri_chrZ|133429644|133435453
EOF

# Display results
echo "===== RNA SUPPORT SUMMARY ====="
column -t "$OUT"

echo
echo "===== ROW COUNT ====="
wc -l "$OUT"
```
THe results will show that seven of the ten loci have detectable venom-gland RNA alignments, but the strength and specificity of that evidence vary substantially. The most striking patterns are the strong RNA signal at the disrupted SVSP-16 and CTL-10 loci, the extensive ambiguous mapping at PLA2-4, and the very weak RNA support for the apparently intact Z-linked CTL-17.


```sh
Gene	  All primary	MAPQ ≥20	High-MAPQ fraction
SVSP-1	        6,098	   6,096	            100.0%
SVSP-16	      147,861	  82,277	             55.6%
Ficolin-4	       85	      85	              100%
PLA2-4	       58,132	     476	              0.8%
CTL-10	       98,510	  46,502	             47.2%
HYAL-4	            0	       0	                 —
CTL-15	           43	      43	              100%
CTL-16	            0	       0	                 —
CTL-17	           11	      11	              100%
CTL-18	           10	       9	               90%
```

Based on this result, the assessment can now be updated to:
 1) SVSP-1: Strong and nearly entirely high-MAPQ RNA coverage across the interval. This is encouraging evidence of locus-specific transcription, although we still need to check whether the predicted CDS containing the late premature stop is the transcript actually expressed.

 2) SVSP-16: Very substantial RNA signal, including more than 82,000 high-MAPQ alignments. This is particularly interesting because its predicted protein contains six premature stops. It could be a transcribed pseudogene, or the current gene model could be incorrect. RNA abundance alone cannot distinguish these explanations.

 3) Ficolin-4: Modest RNA support with no apparent low-MAPQ contribution. Because this locus was identified primarily through toxin-sequence homology, we should determine whether these reads cover its predicted exons and support splice junctions.

 4) PLA2-4: This is the most concerning result for mapping specificity. Only 0.8% of overlapping primary alignments have MAPQ ≥20. The chromosome 16 PLA2 cluster contains highly similar paralogs, so apparent transcription at this locus may largely reflect ambiguous assignment. We should be particularly cautious about interpreting its RNA support.

 5) CTL-10: Substantial RNA coverage, including approximately 46,500 high-MAPQ alignments. Like SVSP-16, it may represent a transcribed disrupted gene or a gene whose current model needs revision.

 6) HYAL-4 and CTL-16: No RNA alignments were detected within their intervals. This provides no transcriptional support in SRR35908235, but does not independently establish pseudogene status.

 7) CTL-15, CTL-17, and CTL-18: All have detectable but sparse RNA support. For CTL-17 specifically, 11 high-MAPQ alignments do not establish strong venom-gland expression, despite its intact ORF and apparent signal peptide. Its toxin status should remain unresolved.

The next step now is to examine the actual CDS-level RNA coverage. The key remaining question is whether RNA reads support the specific CDS structures predicted by ToxCodAn, especially those containing premature stops. One important point for eventual annotation curation is that a premature-stop-containing gene can be strongly transcribed yet encode no functional toxin, while weak venom-gland expression does not rule out functionality. We'll keep those distinctions explicit throughout the analysis.

We can extract the predicted CDS segments for all 10 loci directly from `toxin_annotation.gtf`. We'll use the existing SRR35908235 HISAT2 BAM and calculate the mean RNA read depth, percentage of CDS bases covered by at least one read, number of covered bases, and corresponding gene-wide CDS coverage.

To reduce ambiguous paralog assignments, we'll initially use MAPQ ≥20 and base quality ≥20, excluding secondary, supplementary, duplicate, and QC-failed alignments. Run the python script below:
```sh
python3 /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/toxin_RNA_exon_coverage.py
```

...and inspect the outputs:
```sh
echo "===== GENE-LEVEL CDS COVERAGE ====="
column -t toxin_RNA_validation_CDS_summary.tsv

echo
echo "===== CDS SEGMENT COVERAGE ====="
column -t toxin_RNA_validation_CDS_exons.tsv
```

The output of `toxin_RNA_validation_CDS_summary.tsv` looks like this:
```sh
Gene       Chromosome      Strand  CDS_segments  Total_CDS_bp  Covered_bases  Coverage_percent  Mean_depth
SVSP-1     G_ussuri_chr1   -       5             798           798            100.0             979.45
SVSP-16    G_ussuri_chr10  -       4             810           810            100.0             486.74
Ficolin-4  G_ussuri_chr15  -       9             1113          671            60.29             1.33
PLA2-4     G_ussuri_chr16  +       4             417           417            100.0             37.09
CTL-10     G_ussuri_chr6   +       4             531           172            32.39             2116.58
HYAL-4     G_ussuri_chr6   +       5             1626          0              0.0               0.0
CTL-15     G_ussuri_chrZ   -       2             483           483            100.0             10.85
CTL-16     G_ussuri_chrZ   +       5             585           0              0.0               0.0
CTL-17     G_ussuri_chrZ   -       5             465           184            39.57             0.54
CTL-18     G_ussuri_chrZ   -       4             495           449            90.71             1.88
```

This result show that some disrupted toxin genes have excellent RNA coverage, whereas others have incomplete or no support for their predicted coding sequences.In particular, SVSP-1 and SVSP-16 have exceptionally strong coverage across every predicted CDS segment. Conversely, CTL-10 has extremely uneven coverage, which raises a serious question about its predicted gene structure.

At this point, these are the assessments we can make:
 1) SVSP-1 and SVSP-16: strong transcription despite premature stops. For SVSP-1, 100% of the CDS is covered and all 5 CDS segments have 507–1,176× mean depth. For SVSP-16, 100% of the CDS is covered and all 4 CDS segments have 155–1,129× mean depth. This is interesting because both predicted proteins contain internal stop codons, and genomic reads support those stops. For SVSP-16, all six premature stops occur within CDS segment 1, which has a mean RNA depth of approximately 202×. However, coverage alone cannot establish whether the RNA reads support the same mature transcript or an alternative splicing arrangement. We should not yet classify these as conclusively expressed pseudogenes.

 2) CTL-10: This is the most suspicious gene model; the predicted CDS consists of four segments, but nearly all RNA coverage is concentrated in the first two. CDS 3 has no detectable coverage, and CDS 4 has almost none. The premature stop at amino acid 143 lies in CDS segment 4, which is poorly supported by venom-gland RNA. This suggests that the highly abundant transcript reads at the beginning of CTL-10 might not support the complete ToxCodAn gene model. A truncated transcript, alternative gene model, or incorrect joining of CDS segments are plausible explanations.

 3) PLA2-4: The coverage is good, but mapping specificity remains a concern. All four CDS segments have 100% coverage, including CDS 4, which contains the premature stop at amino acid 113. However, there are only 476 MAPQ ≥20 alignments among 58,132 primary alignments overlapping this locus. The substantial multimapping makes it important to verify that the stop-containing exon and predicted splice junctions are supported by reads that genuinely distinguish PLA2-4 from adjacent paralogs.

 4) Ficolin-4: Only the later coding segments are supported. The first four CDS segments have zero coverage. The remaining five segments show low coverage. Both premature stops occur in CDS 9, which has approximately 91% positional coverage but less than 1× mean depth. This is weak support for a complete nine-segment transcript and reinforces the need for manual gene-model curation.

 5) CTL-15, CTL-16, CTL-17, and CTL-18: The Z-linked candidates show different patterns. For CTL-15, all coding bases are covered, with mean depth 10.85×. This provides reasonably good transcriptional evidence, although its venom-toxin identity remains uncertain. For CTL-16, there is no coverage across any CDS segment. Its two premature stops remain genomically supported, but it has no detectable venom-gland transcription in this dataset. For CTL-17, only the first three CDS segments have coverage; the last two have none. This is weak evidence for a complete transcript and does not justify calling it a confirmed venom toxin. For CTL-18, four segments have broad but low-depth coverage. Transcription is plausible, but venom-specific function remains unestablished.

The next now would be splice-junction validation for all 10 loci. This will test whether RNA reads actually connect adjacent CDS segments predicted by ToxCodAn. Based on our `toxin_annotation.gtf.`, there are 47 CDS segments across 10 genes, which means that there are 37 predicted connections between CDS.

Run this:
```sh
python3 /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/toxin_RNA_splice_junctions.py
```

The script will spit out:
```sh
SVSP-1: 4/4 predicted junctions with >=3 MAPQ20 alignments
SVSP-16: 1/3 predicted junctions with >=3 MAPQ20 alignments
Ficolin-4: 0/8 predicted junctions with >=3 MAPQ20 alignments
PLA2-4: 3/3 predicted junctions with >=3 MAPQ20 alignments
CTL-10: 1/3 predicted junctions with >=3 MAPQ20 alignments
HYAL-4: 0/4 predicted junctions with >=3 MAPQ20 alignments
CTL-15: 1/1 predicted junctions with >=3 MAPQ20 alignments
CTL-16: 0/4 predicted junctions with >=3 MAPQ20 alignments
CTL-17: 0/4 predicted junctions with >=3 MAPQ20 alignments
CTL-18: 1/3 predicted junctions with >=3 MAPQ20 alignments

Completed splice-junction analysis.
Expected junctions: 37
Alternative junctions: 1
```

Let's inspect junction support as well:
```sh
column -t toxin_RNA_expected_junctions.tsv
```

The outputs are highly informative. We can now distinguish toxin loci with RNA-supported coding structures from those whose predicted exon connections are poorly supported. The strongest result is SVSP-1: all four predicted splice junctions have hundreds of high-confidence RNA reads supporting them. On the other hand, SVSP-16 and CTL-10 show abundant transcription but incomplete support for their predicted gene structures. 

We can now summarize the results like so (≥3 reads cutoff means ≥3 high quality reads):
```txt
|Gene	    |  Junctions supported by ≥3 reads |   Interpretation                                            |
|-----------|----------------------------------|-------------------------------------------------------------| 
|SVSP-1	    |                              4/4 |   Excellent support for the predicted CDS structure         |
|SVSP-16	|                              1/3 |   Strong expression, but incomplete structural support      |
|Ficolin-4	|                              0/8 |   No junction support for the predicted nine-CDS model      |
|PLA2-4	    |                              3/3 |   All predicted junctions supported, last one weakly        |
|CTL-10	    |                              1/3 |   Strong evidence for first connection only                 |
|HYAL-4	    |                              0/4 |   No RNA evidence in this dataset                           |
|CTL-15	    |                              1/1 |   Both CDS segments connected by RNA reads                  |
|CTL-16	    |                              0/4 |   No RNA evidence in this dataset                           |
|CTL-17	    |                              0/4 |   Insufficient RNA evidence to validate the structure       |
|CTL-18	    |                              1/3 |   Sparse support for all three junctions                    |
```

Based on this output, SVSP-1 is likely to be a genuine transcript because we now have three independent observations: complete CDS coverage, all four predicted splice junctions supported, and genomic confirmation of the internal stop. The premature stop is at amino acid 265 of a 266-amino-acid predicted translation. That means translation would terminate only one codon earlier than the annotated terminal stop.

For SVSP-16, there is abundant transcription, but the complete gene model is questionable. Because despite hundreds-fold RNA coverage across all four CDS segments, only the last connection is supported. This is a significant discrepancy. It suggests the segments may belong to different isoforms, or that ToxCodAn has incorrectly connected parts of the locus. The six premature stops in CDS1 remain concerning, but we first need to establish which transcripts are actually produced.

The PLA2-4 annotation has considerably stronger evidence than initially expected. All three predicted junctions have RNA support, although the final connection is relatively weak. Because the PLA2 paralogs are highly similar, we should still check that these reads are uniquely assignable.

CTL-10 annotation has a probable gene-model problem. The first connection has 7,045 high-MAPQ supporting reads, whereas the other two have zero. Combined with our earlier RNA coverage results, this basically tells us that the first two CDS segments are highly expressed, but the third and fourth have virtually no support. This is much more consistent with an incomplete or incorrectly joined predicted model than with a confidently validated four-CDS transcript. Crucially, the premature stop at amino acid 143 occurs in CDS4. We therefore have no convincing RNA support that the abundant transcript actually includes this stop-containing segment.

For Z-linked CTL genes:
CTL-15 has complete CDS coverage and five reads supporting its sole predicted splice junction. This is reasonable evidence for an expressed two-CDS gene, although toxin identity remains unconfirmed.

CTL-16 has neither RNA coverage nor junction support, consistent with its provisional classification as a disrupted, possibly nonexpressed CTL-like gene.

CTL-17 has only 11 overlapping RNA alignments and zero detected splice-junction reads. Its intact ORF and apparent signal peptide remain plausible structural features, but we cannot currently validate a complete expressed transcript or identify it as a venom toxin.

CTL-18 has three predicted junctions that have at least one supporting high-MAPQ read (1, 1, and 4). The evidence is weak, but consistent with low-level transcription of its predicted model.

The next step is to inspect the alternative splice junctions:
```sh
echo "===== ALTERNATIVE SPLICE JUNCTION ====="

column -t toxin_RNA_alternative_junctions.tsv
```

This will show that SVSP-1 has one splice junction, but the annotated junction has a substantially stronger support than the alternative junction (887 vs. 16 alignments spanning the splice junction). Possible explanations include genuine low-frequency alternative splicing, transcription from another SVSP paralog, or a minority of incorrectly aligned reads.

Now, let's examine the RNA nucleotide evidence at all 14 premature-stop codons. Because RNA pileup cannot distinguish whether a read comes from the predicted full-length transcript or an alternative isoform, we will interpret RNA nucleotide support together with the splice-junction results, especially for SVSP-16 and CTL-10. 

Let's run this as a SLURM job, because running this on the head node turned out to be time-consuming.
```sh
#!/bin/bash
#SBATCH --job-name=toxin_RNA_stop_validation
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate funannotate

set -euo pipefail


# directories
ROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
WORKDIR="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/SRR35908235_ToxCodAn-Genome"

# input files
REF="${ROOT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
RNA_BAM="${ROOT}/annotation/toxin_gene_annotation_RNAdata_added/venom_RNA_mapping/venom_genome.sorted.bam"
BED="${WORKDIR}/toxin_internal_stops.bed"

# output
OUT="${WORKDIR}/toxin_internal_stops_SRR35908235.mpileup"
TMP="${OUT}.tmp.${SLURM_JOB_ID:-$$}"

cd "$WORKDIR"

trap 'rm -f "$TMP"' EXIT

# ------------------------------------------------------------
# RNA pileup using indexed genomic regions
#
# BED: 0-based, half-open
# mpileup -r: 1-based, inclusive
#
# -a: include zero-depth positions in requested regions
# -q 20: minimum mapping quality
# -Q 20: minimum base quality
# -d 100000: maximum depth
# ------------------------------------------------------------

: > "$TMP"

count=0

while IFS=$'\t' read -r chrom start end gene score strand; do

    region="${chrom}:$((start + 1))-${end}"

    echo "Processing ${gene}: ${region}"

    samtools mpileup \
        -f "$REF" \
        -r "$region" \
        -q 20 \
        -Q 20 \
        -d 100000 \
        -a \
        "$RNA_BAM" \
        >> "$TMP"

    count=$((count + 1))

done < "$BED"


# validate output
lines=$(wc -l < "$TMP")

echo
echo "===== RNA STOP VALIDATION ====="
echo "Stop codons processed: $count"
echo "Pileup positions:       $lines"

if [[ "$count" -ne 14 || "$lines" -ne 42 ]]; then
    echo "ERROR: Expected 14 codons and 42 positions." >&2
    exit 1
fi

mv -f "$TMP" "$OUT"

echo "Output: $OUT"
echo "Status: PASS"
```