## 11) Demographic history
### Step 1: Environment setup and software installation
Set up the analysis directory:
```sh
BASE="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/"

mkdir -p "$BASE"/{00_env,01_reference,02_bams,03_qc,04_vcf,05_masks,06_smc,07_models,08_plots}

cd "$BASE"
```

Now install SMC++ for demographic inference. Mendel already has apptainer as a module, so that we can just pull SMC. The "apptainer pull" command will create a SIF file and installation check will show that the software is correctly installed and launches properly.
```sh
# load appatainer as a module
module load Apptainer/apptainer-1.2.5
apptainer --version

# download smc
cd 00_env

apptainer pull \
    smcpp_latest.sif \
    docker://terhorst/smcpp:latest

# check installation
apptainer exec smcpp_latest.sif smc++ --help    
apptainer exec smcpp_latest.sif smc++ version

# save smc version
apptainer exec smcpp_latest.sif smc++ version > smc_version.txt
```

Now, lets create a conda environment called "smc_tools" to contain software that will be used for BAM/VCF/QC work.
```sh
# from the "00_env" dir
source ~/.bash_profile

conda create -n smc_tools \
    -c conda-forge \
    -c bioconda \
    --strict-channel-priority \
    bwa \
    samtools \
    bcftools \
    htslib \
    bedtools \
    mosdepth \
    vcftools \
    seqkit \
    -y

# activate the conda env
conda activate smc_tools

# verify installation
command -v bwa
command -v samtools
command -v bcftools
command -v bgzip
command -v tabix
command -v bedtools
command -v mosdepth

samtools --version | head -n 2
bcftools --version | head -n 2
bedtools --version
mosdepth --version
```

Record the environment:
```sh
conda env export --no-builds > smc_tools.yml
conda list > smc_tools_conda_list.txt
```

### Step 2: Reference genome preparation
Create a symlink for the softmasked, chromosome-assigned reference genome.
```sh
# cd into relevant dir
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/

# reference genome path
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# make symlink
ln -s "$REF" 01_reference/G_ussuriensis.softmasked.fa

# record a checksum
md5sum "$REF" > 01_reference/G_ussuriensis.softmasked.fa.md5
```

Now, index the genome fasta file:
```sh
# make fasta index file
conda activate smc_tools
samtools faidx "$REF"

# check fasta index file
head -n 30 "${REF}.fai"
```

Also make a contig file:
```sh
# make contig file
cut -f1,2 "${REF}.fai" \
    > 01_reference/reference_contigs.tsv

# view
column -t 01_reference/reference_contigs.tsv
```

In our softmasked reference genome, all chromosomes other than ZW should be regarded as autosomes, because the mitogenome and unplaced scaffolds have been removed.

Create a list of autosomes:
```sh
# from the "demography" dir
cut -f1 01_reference/G_ussuriensis.softmasked.fa.fai \
    | grep -v -E '^G_ussuri_chr[ZW]$' \
    > 01_reference/autosomes.txt

# check
cat 01_reference/autosomes.txt
wc -l 01_reference/autosomes.txt
``` 
Next, create an autosome size file; this will be useful for BED operations.
```sh
awk '
NR==FNR {
    keep[$1]=1
    next
}
($1 in keep) {
    print $1 "\t" $2
}
' \
01_reference/autosomes.txt \
01_reference/reference_contigs.tsv \
> 01_reference/autosomes.genome

# check
column -t 01_reference/autosome.genome
``` 

Create a BED spanning every autosome:
```sh
awk '{print $1 "\t0\t" $2}' \
    01_reference/autosomes.genome \
    > 01_reference/autosomes.bed

# check; BED coordinates start at 0
head 01_reference/autosomes.bed
```

Verify total sequence length for the autosomes:
```sh
awk '{sum += $3-$2} END {
    printf "Autosomal bases: %d\n", sum
    printf "Autosomal Gb: %.6f\n", sum/1000000000
}' 01_reference/autosomes.bed
```
The output is 1.38 Gb.

### Step 3: Prepare mainland sample BAMs
We will use 8 low-coverage whole genome samples out of 12 we used for sex chromosome assignment validation. This is because 4 of those samples were collected from island populations that may have different demographic histories. Also, the BAM files for these  samples already exist because I generated them as part of a broader population genomics study on this species. Therefore, we will simply symlink these files rather than recreating them.  

In the "demography" dir, create a list of mainland samples to use:
```
nano mainland_samples.txt
```