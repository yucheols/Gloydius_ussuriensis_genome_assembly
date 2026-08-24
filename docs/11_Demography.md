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
cat mainland_samples.txt
```

Now set the BAM dir and symlink both BAM and index files for each mainland sample:
```sh
BAMDIR="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam"
while read -r sample; do

    ln -s \
        "${BAMDIR}/${sample}.markdup.bam" \
        "02_bams/${sample}.markdup.bam"

    ln -s \
        "${BAMDIR}/${sample}.markdup.bam.bai" \
        "02_bams/${sample}.markdup.bam.bai"

done < mainland_samples.txt

# check
ls -lh 02_bams/
```

### Step 4: BAM validation and autosomal coverage QC
From the "demography" dir, run:
```sh
# make dir
mkdir -p 03_qc/{bam_validation,flagstat,idxstats,mosdepth}

# create a text file containing the full paths to your 8 mainland BAM files
> 02_bams/mainland_bams.txt

while read -r sample; do
    readlink -f "02_bams/${sample}.markdup.bam"
done < mainland_samples.txt \
    > 02_bams/mainland_bams.txt

# check
cat 02_bams/mainland_bams.txt
wc -l 02_bams/mainland_bams.txt
```

Then run below to automate BAM validation across all files:
```sh
REFDICT="01_reference/reference_contigs.tsv"
OUT="03_qc/bam_validation"

mkdir -p "$OUT"

printf "sample\tquickcheck\tsort_order\treference_match\tread_group_SM\tindex\n" \
    > "$OUT/bam_validation_summary.tsv"

while read -r bam; do

    sample=$(basename "$bam" .markdup.bam)

    echo "Checking $sample ..."

    # -------------------------
    # BAM integrity
    # -------------------------
    if samtools quickcheck -v "$bam"; then
        quickcheck="PASS"
    else
        quickcheck="FAIL"
    fi

    # -------------------------
    # sorting
    # -------------------------
    sort_order=$(
        samtools view -H "$bam" |
        awk -F'\t' '
        /^@HD/ {
            for(i=1;i<=NF;i++) {
                if($i ~ /^SO:/) {
                    print substr($i,4)
                }
            }
        }'
    )

    [[ -z "$sort_order" ]] && sort_order="UNKNOWN"

    # -------------------------
    # BAM reference dictionary
    # -------------------------
    samtools view -H "$bam" |
    awk -F'\t' '
    /^@SQ/ {
        sn="";
        ln="";
        for(i=1;i<=NF;i++) {
            if($i ~ /^SN:/) sn=substr($i,4);
            if($i ~ /^LN:/) ln=substr($i,4);
        }
        print sn "\t" ln
    }' > "$OUT/${sample}.reference_contigs.tsv"

    if diff -q \
        "$REFDICT" \
        "$OUT/${sample}.reference_contigs.tsv" \
        >/dev/null 2>&1
    then
        reference_match="PASS"
    else
        reference_match="FAIL"
    fi

    # -------------------------
    # read-group sample name
    # -------------------------
    sm=$(
        samtools view -H "$bam" |
        awk -F'\t' '
        /^@RG/ {
            for(i=1;i<=NF;i++) {
                if($i ~ /^SM:/) {
                    print substr($i,4)
                }
            }
        }' |
        sort -u |
        paste -sd ',' -
    )

    [[ -z "$sm" ]] && sm="NONE"

    # -------------------------
    # BAM index
    # -------------------------
    if [[ -f "${bam}.bai" || -f "${bam%.bam}.bai" ]]; then
        index_status="PASS"
    else
        index_status="FAIL"
    fi

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" \
        "$quickcheck" \
        "$sort_order" \
        "$reference_match" \
        "$sm" \
        "$index_status" \
        >> "$OUT/bam_validation_summary.tsv"

done < 02_bams/mainland_bams.txt
```

Then check - should see PASS for all samples.
```sh
column -t -s $'\t' \
    03_qc/bam_validation/bam_validation_summary.tsv
```

Next, generate flagstat and idxstats for all eight
```sh
while read -r bam; do

    sample=$(basename "$bam" .markdup.bam)

    echo "Processing $sample ..."

    samtools flagstat -@ 8 "$bam" \
        > "03_qc/flagstat/${sample}.flagstat.txt"

    samtools idxstats "$bam" \
        > "03_qc/idxstats/${sample}.idxstats.tsv"

done < 02_bams/mainland_bams.txt
```

Get a quick look at mapping, pairing, and duplicates:
```sh
for f in 03_qc/flagstat/*.flagstat.txt; do

    echo "===== $(basename "$f" .flagstat.txt) ====="

    grep -E \
        'mapped \(|properly paired|duplicates' \
        "$f"

    echo

done
```
We can see the mapping rates are extremely high across all 8 samples (> 99%). Properly paired reads are also high (96.3–98.1%)

Now run autosomal coverage QC with mosdepth. Run this as an array job on Mendel:
```sh
#!/bin/bash
#SBATCH --job-name=smc_mosdepth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-8%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

set -euo pipefail

# --------------------------------------------------
# activate conda
# --------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

# --------------------------------------------------
# directories
# --------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"

SAMPLEFILE="${WORKDIR}/mainland_samples.txt"
BAMDIR="${WORKDIR}/02_bams"
AUTOSOMES="${WORKDIR}/01_reference/autosomes.bed"
OUTDIR="${WORKDIR}/03_qc/mosdepth"

mkdir -p "$OUTDIR"

# --------------------------------------------------
# get sample from array index
# --------------------------------------------------
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLEFILE")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

BAM="${BAMDIR}/${sample}.markdup.bam"

if [[ ! -f "$BAM" ]]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

if [[ ! -f "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found: ${BAM}.bai"
    exit 1
fi

echo "Sample: $sample"
echo "BAM: $BAM"
echo "Autosome BED: $AUTOSOMES"
echo "Starting: $(date)"

# --------------------------------------------------
# autosomal coverage
# --------------------------------------------------
mosdepth \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --mapq 30 \
    --by "$AUTOSOMES" \
    --no-per-base \
    "${OUTDIR}/${sample}" \
    "$BAM"

echo "Finished: $(date)"
```