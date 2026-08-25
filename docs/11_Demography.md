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

When this job finishes running, summarize the output into a single table using the script below:
```sh
# run from "demography" dir
printf "sample\tautosomal_bp\tmean_autosomal_depth\n" \
    > 03_qc/mosdepth/autosomal_mean_depth.tsv

while read -r sample; do

    zcat "03_qc/mosdepth/${sample}.regions.bed.gz" |
    awk -v sample="$sample" '
    BEGIN {OFS="\t"}
    {
        len=$3-$2
        bp += len
        covsum += len*$4
    }
    END {
        print sample,bp,covsum/bp
    }'

done < mainland_samples.txt \
    >> 03_qc/mosdepth/autosomal_mean_depth.tsv
```

Check:
```sh
column -t -s $'\t' \
    03_qc/mosdepth/autosomal_mean_depth.tsv
```
The output should look like this:
```sh
sample      autosomal_bp  mean_autosomal_depth
AMNH_21010  1380185643    11.8822
AMNH_21128  1380185643    11.6865
AMNH_21147  1380185643    11.2356
AMNH_21161  1380185643    11.9259
AMNH_21162  1380185643    12.5469
AMNH_21164  1380185643    11.958
AMNH_21172  1380185643    13.0561
AMNH_21185  1380185643    11.4408
```

Also make a chromosome-by-chromosome summary table:
```sh
printf "sample\tchromosome\tstart\tend\tmean_depth\n" \
    > 03_qc/mosdepth/autosomal_chromosome_depth.tsv

while read -r sample; do

    zcat "03_qc/mosdepth/${sample}.regions.bed.gz" |
    awk -v sample="$sample" '
    BEGIN {OFS="\t"}
    {
        print sample,$1,$2,$3,$4
    }'

done < mainland_samples.txt \
    >> 03_qc/mosdepth/autosomal_chromosome_depth.tsv
```
Check:
```sh
column -t -s $'\t' \
    03_qc/mosdepth/autosomal_chromosome_depth.tsv \
    | head -n 40
```

The output should look like this:
```sh
sample      chromosome      start  end        mean_depth
AMNH_21010  G_ussuri_chr1   0      340623792  11.62
AMNH_21010  G_ussuri_chr2   0      272076671  11.87
AMNH_21010  G_ussuri_chr3   0      203076995  11.43
AMNH_21010  G_ussuri_chr4   0      122805809  11.28
AMNH_21010  G_ussuri_chr5   0      99259852   11.51
AMNH_21010  G_ussuri_chr6   0      93675145   11.51
AMNH_21010  G_ussuri_chr7   0      79146131   11.41
AMNH_21010  G_ussuri_chr9   0      24199037   13.34
AMNH_21010  G_ussuri_chr10  0      27181165   13.79
AMNH_21010  G_ussuri_chr11  0      17175000   14.51
AMNH_21010  G_ussuri_chr12  0      21009094   13.60
AMNH_21010  G_ussuri_chr13  0      17611041   14.37
AMNH_21010  G_ussuri_chr14  0      19255245   13.71
AMNH_21010  G_ussuri_chr15  0      11289862   14.34
AMNH_21010  G_ussuri_chr16  0      8920326    14.31
AMNH_21010  G_ussuri_chr17  0      13149378   14.65
AMNH_21010  G_ussuri_chr18  0      9731100    15.42
AMNH_21128  G_ussuri_chr1   0      340623792  11.65
AMNH_21128  G_ussuri_chr2   0      272076671  11.69
AMNH_21128  G_ussuri_chr3   0      203076995  11.60
AMNH_21128  G_ussuri_chr4   0      122805809  11.23
AMNH_21128  G_ussuri_chr5   0      99259852   11.44
AMNH_21128  G_ussuri_chr6   0      93675145   11.29
AMNH_21128  G_ussuri_chr7   0      79146131   11.22
AMNH_21128  G_ussuri_chr9   0      24199037   12.42
AMNH_21128  G_ussuri_chr10  0      27181165   12.53
AMNH_21128  G_ussuri_chr11  0      17175000   13.32
AMNH_21128  G_ussuri_chr12  0      21009094   12.65
AMNH_21128  G_ussuri_chr13  0      17611041   12.99
AMNH_21128  G_ussuri_chr14  0      19255245   12.94
AMNH_21128  G_ussuri_chr15  0      11289862   12.45
AMNH_21128  G_ussuri_chr16  0      8920326    12.50
AMNH_21128  G_ussuri_chr17  0      13149378   12.92
AMNH_21128  G_ussuri_chr18  0      9731100    13.28
AMNH_21147  G_ussuri_chr1   0      340623792  11.22
AMNH_21147  G_ussuri_chr2   0      272076671  11.25
AMNH_21147  G_ussuri_chr3   0      203076995  11.14
AMNH_21147  G_ussuri_chr4   0      122805809  10.69
AMNH_21147  G_ussuri_chr5   0      99259852   11.01
```

These coverage results are very good for SMC++. Across the eight mainland samples, mean autosomal depth ranges only from 11.24× to 13.06×, with an overall mean of about 11.97×.

### Step 5: Set depth cutoff
mosdepth generated files like:
```sh
03_qc/mosdepth/AMNH_21010.mosdepth.global.dist.txt
```

First, extract converage quantile across all samples:
```sh
printf "sample\tq2.5\tq5\tmedian\tq95\tq97.5\n" \
    > 03_qc/mosdepth/autosomal_depth_quantiles.tsv

while read -r sample; do

    awk -v sample="$sample" '
    $1=="total" {

        depth=$2
        frac=$3

        # mosdepth frac = proportion of bases with coverage >= depth.
        # Rows run from high depth toward depth 0.

        if (q025=="" && frac >= 0.975) q025=depth
        if (q05==""  && frac >= 0.950) q05=depth
        if (q50==""  && frac >= 0.500) q50=depth
        if (q95==""  && frac >= 0.050) q95=depth
        if (q975=="" && frac >= 0.025) q975=depth
    }

    END {
        print sample "\t" q025 "\t" q05 "\t" q50 "\t" q95 "\t" q975
    }' \
    "03_qc/mosdepth/${sample}.mosdepth.global.dist.txt"

done < mainland_samples.txt \
    >> 03_qc/mosdepth/autosomal_depth_quantiles.tsv
```

Check:
```sh
column -t -s $'\t' \
    03_qc/mosdepth/autosomal_depth_quantiles.tsv
```

The output should look like:
```sh
sample      q2.5  q5  median  q95  q97.5
AMNH_21010  0     2   12      21   23
AMNH_21128  0     0   12      21   22
AMNH_21147  0     0   12      20   22
AMNH_21161  0     0   12      21   23
AMNH_21162  0     0   13      22   23
AMNH_21164  0     0   13      22   23
AMNH_21172  0     0   14      23   25
AMNH_21185  0     0   12      21   23
```
Thus, the median depths are 12–14×, while the upper 97.5th-percentile depths are 22–25×. Next, let's also calculate what percentage of the autosomal genome has at least 1×, 3×, 5×, 7×, and 10× coverage in each sample.

```sh
printf "sample\tge1\tge3\tge5\tge7\tge10\n" \
    > 03_qc/mosdepth/autosomal_callable_fractions.tsv

while read -r sample; do

    file="03_qc/mosdepth/${sample}.mosdepth.global.dist.txt"

    ge1=$(awk '$1=="total" && $2==1  {print $3}' "$file")
    ge3=$(awk '$1=="total" && $2==3  {print $3}' "$file")
    ge5=$(awk '$1=="total" && $2==5  {print $3}' "$file")
    ge7=$(awk '$1=="total" && $2==7  {print $3}' "$file")
    ge10=$(awk '$1=="total" && $2==10 {print $3}' "$file")

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$sample" "$ge1" "$ge3" "$ge5" "$ge7" "$ge10"

done < mainland_samples.txt \
    >> 03_qc/mosdepth/autosomal_callable_fractions.tsv
```
Check:
```sh
column -t -s $'\t' \
    03_qc/mosdepth/autosomal_callable_fractions.tsv
```

The output should look like this:
```sh
sample      ge1   ge3   ge5   ge7   ge10
AMNH_21010  0.97  0.94  0.90  0.84  0.67
AMNH_21128  0.93  0.90  0.87  0.83  0.69
AMNH_21147  0.93  0.90  0.86  0.81  0.66
AMNH_21161  0.93  0.90  0.87  0.82  0.70
AMNH_21162  0.93  0.90  0.88  0.84  0.74
AMNH_21164  0.93  0.89  0.86  0.81  0.69
AMNH_21172  0.93  0.90  0.87  0.83  0.73
AMNH_21185  0.93  0.89  0.84  0.78  0.64
```

At 5× (ge5), each individual retains 84–90% of autosomal bases. Thus, it is sensible to set 5× as the minimum callable depth. Also, the highest 97.5th-percentile depth value was 25x. So, set this as the maximum callable depth.

### Step 6: Create actual callable-region BEDs
The summary files we generated above tell us how much of the genome passes certain depth thresholds, but they don't tell us where those bases are. SMC++ needs that spatial information.

For that, we will have to run a second mosdepth jon with a --quantize mode: this will collapse adjacent bases into coverage categories such as low, callable, and high coverage.
We will define 0–4× as low coverage, 5–25× as callable, and 26× as high coverage.

Use this script:
```sh
#!/bin/bash
#SBATCH --job-name=smc_callability
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-8%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ============================================================
# generate depth-based callability masks for SMC++
#
# coverage classes:
#   0x       = NO_COVERAGE
#   1-4x     = LOW_COVERAGE
#   5-25x    = CALLABLE
#   >=26x    = HIGH_COVERAGE
#
# one sample is processed per SLURM array task.
# ============================================================

# ------------------------------------------------------------
# activate conda environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
SAMPLEFILE="${WORKDIR}/mainland_samples.txt"
BAMDIR="${WORKDIR}/02_bams"
OUTDIR="${WORKDIR}/03_qc/mosdepth_callability"

mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# get sample corresponding to SLURM array task
# ------------------------------------------------------------
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLEFILE")

if [[ -z "$sample" ]]; then
    echo "ERROR: No sample found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


# ------------------------------------------------------------
# BAM
# ------------------------------------------------------------
BAM="${BAMDIR}/${sample}.markdup.bam"

if [[ ! -e "$BAM" ]]; then
    echo "ERROR: BAM not found:"
    echo "$BAM"
    exit 1
fi

if [[ ! -e "${BAM}.bai" ]]; then
    echo "ERROR: BAM index not found:"
    echo "${BAM}.bai"
    exit 1
fi


# ------------------------------------------------------------
# print job information
# ------------------------------------------------------------

echo "============================================================"
echo "SMC++ callability mask"
echo "============================================================"
echo "Sample:           $sample"
echo "BAM:              $BAM"
echo "Output directory: $OUTDIR"
echo "Array task:       ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "Host:             $(hostname)"
echo "Start time:       $(date)"
echo "============================================================"


# ------------------------------------------------------------
# define labels for mosdepth quantization
#
# --quantize 0:1:5:26:
#
# bins:
#   Q0 = 0x
#   Q1 = 1-4x
#   Q2 = 5-25x
#   Q3 = >=26x
# ------------------------------------------------------------
export MOSDEPTH_Q0="NO_COVERAGE"
export MOSDEPTH_Q1="LOW_COVERAGE"
export MOSDEPTH_Q2="CALLABLE"
export MOSDEPTH_Q3="HIGH_COVERAGE"


# ------------------------------------------------------------
# run mosdepth
#
# --mapq 30:
#   ignore reads with mapping quality <30
#
# --flag 1796:
#   exclude unmapped, secondary, QC-failed, and duplicate reads
#
# --no-per-base:
#   do not produce huge per-base BED files
#
# --quantize:
#   merge adjacent bases belonging to the same coverage class
# ------------------------------------------------------------
mosdepth \
    --threads "${SLURM_CPUS_PER_TASK}" \
    --mapq 30 \
    --flag 1796 \
    --no-per-base \
    --quantize 0:1:5:26: \
    "${OUTDIR}/${sample}.callability" \
    "$BAM"


# ------------------------------------------------------------
# check expected output
# ------------------------------------------------------------
QUANT="${OUTDIR}/${sample}.callability.quantized.bed.gz"

if [[ ! -s "$QUANT" ]]; then
    echo "ERROR: Expected quantized BED was not produced:"
    echo "$QUANT"
    exit 1
fi


# ------------------------------------------------------------
# report class counts
# ------------------------------------------------------------
echo
echo "Coverage-class intervals:"
zcat "$QUANT" |
awk '
{
    n[$4]++
    bp[$4] += $3-$2
}
END {
    for (x in n) {
        printf "%-15s intervals=%d bp=%d\n", x, n[x], bp[x]
    }
}'


echo
echo "============================================================"
echo "Completed sample: $sample"
echo "Finish time:      $(date)"
echo "============================================================"
```

Confirm the output:
```sh
ls -lh 03_qc/mosdepth_callability/*.quantized.bed.gz
```
And make a mask directory:
```sh
mkdir -p 05_masks/per_sample
```

Now, convert those eight quantized BED files into autosomal callable BEDs. We will extract callable autosomal regions for each individual
```sh
conda activate smc_tools

while read -r sample; do

    echo "Processing $sample ..."

    zcat "03_qc/mosdepth_callability/${sample}.callability.quantized.bed.gz" |
    awk '
        NR==FNR {
            keep[$1]=1
            next
        }
        ($1 in keep) && $4=="CALLABLE" {
            print $1 "\t" $2 "\t" $3
        }
    ' 01_reference/autosomes.txt - |
    bedtools sort \
        -g 01_reference/autosomes.genome \
        -i - \
        > "05_masks/per_sample/${sample}.callable.autosomes.bed"

done < mainland_samples.txt
```

Now make a summary table containing information on how much callable sequence each sample has:
```sh
printf "sample\tcallable_bp\tcallable_fraction\n" \
    > 05_masks/per_sample_callable_summary.tsv

while read -r sample; do

    bp=$(
        awk '{sum += $3-$2} END {print sum}' \
        "05_masks/per_sample/${sample}.callable.autosomes.bed"
    )

    awk -v sample="$sample" \
        -v bp="$bp" \
        -v total=1380185643 \
        'BEGIN {
            printf "%s\t%d\t%.4f\n", sample,bp,bp/total
        }'

done < mainland_samples.txt \
    >> 05_masks/per_sample_callable_summary.tsv
```

Check the output:
```sh
column -t -s $'\t' \
    05_masks/per_sample_callable_summary.tsv
```

The output should look like:
```sh
sample      callable_bp  callable_fraction
AMNH_21010  1231746411   0.8924
AMNH_21128  1197101199   0.8673
AMNH_21147  1185567014   0.8590
AMNH_21161  1187123569   0.8601
AMNH_21162  1202950360   0.8716
AMNH_21164  1172894304   0.8498
AMNH_21172  1170033533   0.8477
AMNH_21185  1153023324   0.8354
```

Next, find regions callable in all 8 samples:
```sh
> 05_masks/callable_bed_files.txt

while read -r sample; do
    echo "05_masks/per_sample/${sample}.callable.autosomes.bed"
done < mainland_samples.txt \
    > 05_masks/callable_bed_files.txt
```

And load the file names into a bash array:
```sh
mapfile -t CALLABLE_FILES < 05_masks/callable_bed_files.txt
```

Now that this is done, find intervals callable in all 8 samples:
```sh
bedtools multiinter \
    -i "${CALLABLE_FILES[@]}" |
awk '$4==8 {
    print $1 "\t" $2 "\t" $3
}' \
> 05_masks/shared_callable_depth.autosomes.bed
```

Also check how much of the genome has 5–25× high-MAPQ coverage in every one of the eight individuals:
```sh
CALLABLE_BP=$(
    awk '{sum += $3-$2} END {print sum}' \
    05_masks/shared_callable_depth.autosomes.bed
)

awk -v bp="$CALLABLE_BP" \
    -v total=1380185643 \
    'BEGIN {
        printf "Shared callable bp: %d\n",bp
        printf "Shared callable Gb: %.3f\n",bp/1e9
        printf "Shared callable fraction: %.4f\n",bp/total
        printf "Shared callable percent: %.2f%%\n",(bp/total)*100
    }'
```

The output should be:
```sh
Shared callable bp: 953290613
Shared callable Gb: 0.953
Shared callable fraction: 0.6907
Shared callable percent: 69.07%
```
This means that 69.07% of our autosomal genome is callable in all eight individuals simultaneously under the 5–25×, MAPQ ≥30 criterion.

### Step 7: Make the depth-based SMC++ exclusion mask
First, merge and sort the callable intervals:
```sh
bedtools sort \
    -g 01_reference/autosomes.genome \
    -i 05_masks/shared_callable_depth.autosomes.bed |
bedtools merge -i - \
    > 05_masks/shared_callable_depth.autosomes.merged.bed
```

Then take the complement (mask), which contains what we do not want SMC++ to use:
```sh
bedtools complement \
    -i 05_masks/shared_callable_depth.autosomes.merged.bed \
    -g 01_reference/autosomes.genome \
    > 05_masks/noncallable_depth.autosomes.bed
```

Check the size of this mask file. This should be 0.427 Gb (30.93% masked percentage).
```sh
awk '
{
    bp += $3-$2
}
END {
    printf "Masked bp: %d\n", bp
    printf "Masked Gb: %.3f\n", bp/1e9
    printf "Masked percent: %.2f%%\n", (bp/1380185643)*100
}' 05_masks/noncallable_depth.autosomes.bed
```

Also run this sanity check. The expected output is "Total: 1380185643"
```sh
CALLABLE=$(
    awk '{x += $3-$2} END {print x}' \
    05_masks/shared_callable_depth.autosomes.merged.bed
)

MASKED=$(
    awk '{x += $3-$2} END {print x}' \
    05_masks/noncallable_depth.autosomes.bed
)

echo "Callable: $CALLABLE"
echo "Masked:   $MASKED"
echo "Total:    $((CALLABLE + MASKED))"
```

Now, we need to make sure that every N region in the reference is also masked.
```sh
REF="01_reference/G_ussuriensis.softmasked.fa"
awk '
BEGIN {OFS="\t"}

/^>/ {
    if (inN) {
        print seqname,start,pos
        inN=0
    }

    seqname=substr($1,2)
    pos=0
    next
}

{
    line=toupper($0)

    for (i=1; i<=length(line); i++) {

        base=substr(line,i,1)

        if (base=="N" && !inN) {
            start=pos+i-1
            inN=1
        }

        if (base!="N" && inN) {
            print seqname,start,pos+i-1
            inN=0
        }
    }

    pos += length(line)
}

END {
    if (inN) {
        print seqname,start,pos
    }
}
' "$REF" \
> 05_masks/reference_N.all.bed
```

Restrict this to 17 autosomes in the assembly:
```sh
awk '
NR==FNR {
    keep[$1]=1
    next
}
($1 in keep)
' \
01_reference/autosomes.txt \
05_masks/reference_N.all.bed \
> 05_masks/reference_N.autosomes.bed
```

Check how much autosomal sequence is N:
```sh
awk '
{
    bp += $3-$2
}
END {
    printf "Autosomal N bp: %d\n",bp
    printf "Autosomal N Mb: %.3f\n",bp/1e6
    printf "Autosomal N percent: %.4f%%\n",(bp/1380185643)*100
}' 05_masks/reference_N.autosomes.bed
```

And check whether any of those N bases are currently considered callable:
```sh
bedtools intersect \
    -a 05_masks/reference_N.autosomes.bed \
    -b 05_masks/shared_callable_depth.autosomes.merged.bed \
    > 05_masks/reference_N_callable_overlap.bed
```           

Check if the autosomal Ns and callable autosomal Ns overlap. The output should be 0 because reference N regions should have no valid read coverage and therefore should already have failed the 5–25× criterion.
```sh
awk '
{
    bp += $3-$2
}
END {
    printf "N bases incorrectly callable: %d\n",bp+0
}' 05_masks/reference_N_callable_overlap.bed
```

Our primary SMC++ mask is now:
```sh
05_masks/noncallable_depth.autosomes.bed
```

Run a final integrity check:
```sh
bedtools sort \
    -g 01_reference/autosomes.genome \
    -i 05_masks/noncallable_depth.autosomes.bed |
bedtools merge -i - \
    > 05_masks/smcpp_exclude.primary.autosomes.sorted.bed
```

Confirm the size is still exactly the same as expected (426895030, 30.93%):
```sh
awk '
{
    bp += $3-$2
}
END {
    printf "Final masked bp: %d\n",bp
    printf "Final masked percent: %.2f%%\n",(bp/1380185643)*100
}' 05_masks/smcpp_exclude.primary.autosomes.sorted.bed
```

### Step 8: Joint variant calling using BCFtools
The script will call one autosome per SLURM array task. However, all 8 individuals are called together in every task. So this is still joint calling; we're only splitting the 1.38-Gb genome by chromosome to make the computation manageable. This script will generate a relatively complete raw joint VCF; filtering comes immediately afterward.

```sh
#!/bin/bash
#SBATCH --job-name=smc_jvarcall
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ============================================================
# Joint variant calling for SMC++
#
# Eight mainland G. ussuriensis individuals are jointly called
# One autosome is processed per SLURM array task
#
# This is RAW variant calling
# Final SNP/DP/QUAL/missingness filtering is done afterward
# ============================================================

# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
REF="${WORKDIR}/01_reference/G_ussuriensis.softmasked.fa"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
BAMLIST="${WORKDIR}/02_bams/mainland_bams.txt"
OUTDIR="${WORKDIR}/04_vcf/raw"
STATDIR="${WORKDIR}/04_vcf/raw_stats"

mkdir -p "$OUTDIR"
mkdir -p "$STATDIR"


# ------------------------------------------------------------
# check input files
# ------------------------------------------------------------
if [[ ! -s "$REF" ]]; then
    echo "ERROR: Reference not found:"
    echo "$REF"
    exit 1
fi

if [[ ! -s "${REF}.fai" ]]; then
    echo "ERROR: Reference FASTA index not found:"
    echo "${REF}.fai"
    exit 1
fi

if [[ ! -s "$CHRFILE" ]]; then
    echo "ERROR: Autosome list not found:"
    echo "$CHRFILE"
    exit 1
fi

if [[ ! -s "$BAMLIST" ]]; then
    echo "ERROR: BAM list not found:"
    echo "$BAMLIST"
    exit 1
fi

# ------------------------------------------------------------
# verify number of BAMs
# ------------------------------------------------------------
NBAMS=$(grep -cve '^[[:space:]]*$' "$BAMLIST")

if [[ "$NBAMS" -ne 8 ]]; then
    echo "ERROR: Expected 8 BAMs but found $NBAMS"
    echo "BAM list: $BAMLIST"
    exit 1
fi

echo "Found $NBAMS BAMs in BAM list."


# ------------------------------------------------------------
# get chromosome for this array task
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")

if [[ -z "$CHR" ]]; then
    echo "ERROR: No chromosome found for array task:"
    echo "${SLURM_ARRAY_TASK_ID}"
    exit 1
fi


# ------------------------------------------------------------
# output files
# ------------------------------------------------------------
VCF="${OUTDIR}/${CHR}.raw.vcf.gz"
STATS="${STATDIR}/${CHR}.raw.stats.txt"


# ------------------------------------------------------------
# print job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Joint variant calling"
echo "============================================================"
echo "Chromosome:       $CHR"
echo "Number of BAMs:   $NBAMS"
echo "Reference:        $REF"
echo "BAM list:         $BAMLIST"
echo "Output:           $VCF"
echo "Array task:       ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "Start time:       $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# check BAMs before calling
# ------------------------------------------------------------
while read -r bam; do

    if [[ ! -s "$bam" ]]; then
        echo "ERROR: BAM not found:"
        echo "$bam"
        exit 1
    fi

    samtools quickcheck "$bam"

done < "$BAMLIST"


# ------------------------------------------------------------
# joint variant calling
#
# mpileup:
#
#   -f       reference genome
#   -r       current chromosome
#   -b       BAM list containing all 8 individuals
#   -q 30    require mapping quality >=30
#   -Q 20    require base quality >=20
#   -d 1000  allow up to 1000 reads/sample/site
#
#   --ff:
#     exclude unmapped, secondary, QC-failed,
#     and duplicate-marked reads
#
#   -a:
#     retain per-individual depth (DP)
#     and allele depth (AD)
#
# call:
#
#   -m       multiallelic calling model
#   -v       output variant sites only
#   -f GQ    include genotype quality
#
# no final SNP or depth filtering is performed here.
# ------------------------------------------------------------

bcftools mpileup \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f "$REF" \
    -r "$CHR" \
    -b "$BAMLIST" \
    -q 30 \
    -Q 20 \
    -d 1000 \
    --ff UNMAP,SECONDARY,QCFAIL,DUP \
    -a FORMAT/DP,FORMAT/AD \
    -Ou \
| \
bcftools call \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -m \
    -v \
    -f GQ \
    -Oz \
    -o "$VCF"


# ------------------------------------------------------------
# index VCF
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$VCF"


# ------------------------------------------------------------
# generate raw VCF statistics
# ------------------------------------------------------------
bcftools stats "$VCF" > "$STATS"


# ------------------------------------------------------------
# quick output validation
# ------------------------------------------------------------
echo
echo "Samples in VCF:"

bcftools query -l "$VCF"

NSAMPLES=$(bcftools query -l "$VCF" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: VCF contains $NSAMPLES samples instead of 8."
    exit 1
fi


# ------------------------------------------------------------
# count raw variants
# ------------------------------------------------------------
NVAR=$(bcftools view -H "$VCF" | wc -l)

echo
echo "Raw variant records: $NVAR"


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Completed chromosome: $CHR"
echo "Raw variants:         $NVAR"
echo "Output:               $VCF"
echo "Finish time:          $(date)"
echo "============================================================"
```

Once joint variant calling finishes running, QC the raw VCF. First confirm all 17 chromosome VCFs and indexes exist:
```sh
ls 04_vcf/raw/*.raw.vcf.gz | wc -l
ls 04_vcf/raw/*.raw.vcf.gz.tbi | wc -l
```
Both should return 17.

Then confirm the eight sample names in one VCF:
```sh
bcftools query -l \
    04_vcf/raw/G_ussuri_chr1.raw.vcf.gz
```

Next, verify that the fields we need survived variant calling. These fields include DP, AD, GT, GQ.
```sh
bcftools view -h \
    04_vcf/raw/G_ussuri_chr1.raw.vcf.gz |
grep -E '^##FORMAT=<ID=(GT|DP|AD|GQ),'
```

Now, let's summarize raw variants across all 17 autosomes. First, create:
```sh
printf "chromosome\trecords\tsnps\tindels\tmultiallelic\n" \
    > 04_vcf/raw_stats/raw_variant_summary.tsv
```
And then:
```sh
for f in 04_vcf/raw_stats/*.raw.stats.txt; do

    chr=$(basename "$f" .raw.stats.txt)

    records=$(awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' "$f")

    snps=$(awk -F'\t' \
        '$1=="SN" && $3=="number of SNPs:" {print $4}' "$f")

    indels=$(awk -F'\t' \
        '$1=="SN" && $3=="number of indels:" {print $4}' "$f")

    multi=$(awk -F'\t' \
        '$1=="SN" && $3=="number of multiallelic sites:" {print $4}' "$f")

    printf "%s\t%s\t%s\t%s\t%s\n" \
        "$chr" "$records" "$snps" "$indels" "$multi"

done >> 04_vcf/raw_stats/raw_variant_summary.tsv
```

View it:
```sh
column -t -s $'\t' \
    04_vcf/raw_stats/raw_variant_summary.tsv
```

Next, get totals:
```sh
awk '
NR>1 {
    records += $2
    snps += $3
    indels += $4
    multi += $5
}
END {
    printf "Total records:             %d\n",records
    printf "Total SNPs:                %d\n",snps
    printf "Total indels:              %d\n",indels
    printf "Total multiallelic sites:  %d\n",multi
}' 04_vcf/raw_stats/raw_variant_summary.tsv
```

### Step 9: VCF filtering - Stage 1 (quality filtering)
This step takes the raw joint VCF ans keeps only shared-callable regions, SNPs, and biallelic sites. It also requires site QUAL >= 20, and consider DP < 5, DP > 25, and GQ < 20 as missing.

Below is the overall workflow:
```sh
1) raw joint VCF
2) normalize alleles
3) retain biallelic SNPs
4) genotype DP/GQ filtering
5) site QUAL filtering
6) callable-region restriction
7) final VCF QC
8) SMC++ vcf2smc
```

Run the first filtering stage with the script below:
```sh
#!/bin/bash
#SBATCH --job-name=smc_vcf_filt_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ================================================================
# SMC++ VCF filtering - Stage 1
#
# 1) starting data:
#   raw joint-called autosomal VCFs
#
# 2) site filters:
#   shared callable genome only
#   SNPs only
#   biallelic sites only
#   QUAL >= 20
#
# 3) genotype filters:
#   DP 5-25
#   GQ >= 20
#
# genotypes failing DP/GQ are changed to ./.
#
# this script DOES NOT remove sites based on missingness just yet.
# ================================================================

# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
INDIR="${WORKDIR}/04_vcf/raw"
OUTDIR="${WORKDIR}/04_vcf/filtered_stage1"
STATDIR="${WORKDIR}/04_vcf/filtered_stage1_stats"
CALLABLE="${WORKDIR}/05_masks/shared_callable_depth.autosomes.merged.bed"
RAWSTATDIR="${WORKDIR}/04_vcf/raw_stats"

mkdir -p "$OUTDIR"
mkdir -p "$STATDIR"


# ------------------------------------------------------------
# get chromosome
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")


# ------------------------------------------------------------
# files
# ------------------------------------------------------------
RAW="${INDIR}/${CHR}.raw.vcf.gz"
OUT="${OUTDIR}/${CHR}.filtered_stage1.vcf.gz"
STATS="${STATDIR}/${CHR}.filtered_stage1.stats.txt"
RAWSTATS="${RAWSTATDIR}/${CHR}.raw.stats.txt"


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------

echo
echo "============================================================"
echo "SMC++ VCF filtering - Stage 1"
echo "============================================================"
echo "Chromosome:     $CHR"
echo "Raw VCF:        $RAW"
echo "Callable BED:   $CALLABLE"
echo "Output:         $OUT"
echo "Array task:     ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:           ${SLURM_CPUS_PER_TASK}"
echo "Host:           $(hostname)"
echo "Start:          $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# FILTERING
#
# 1) Step 1:
#   Remove INFO/MQ because its raw header definition causes
#   the bcftools MQ sanity warning and we do not use this tag.
#
# 2) Step 2:
#   Retain positions within shared callable regions.
#
# 3)  Step 3:
#   Retain SNPs only.
#
# 4) Step 4:
#   Retain exactly two alleles (REF + one ALT).
#
# 5) Step 5:
#   Require site QUAL >=20.
#
# 6) Step 6:
#   Change an individual genotype to missing (./.) if:
#
#       DP <5
#       DP >25
#       DP missing
#       GQ <20
#       GQ missing
#
# IMPORTANT:
# Single "|" operators are intentional. For FORMAT fields,
# this evaluates the conditions at the individual-sample level.
# ------------------------------------------------------------
bcftools annotate \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -x INFO/MQ \
    -Ou \
    "$RAW" \
| \
bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -T "$CALLABLE" \
    --targets-overlap 0 \
    -v snps \
    -m2 \
    -M2 \
    -i 'QUAL>=20' \
    -Ou \
| \
bcftools +setGT \
    -Ou \
    -- \
    -t q \
    -n . \
    -i 'FMT/DP="." | FMT/DP<5 | FMT/DP>25 | FMT/GQ="." | FMT/GQ<20' \
| \
bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -Oz \
    -o "$OUT"


# ------------------------------------------------------------
# index filtered VCF
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$OUT"


# ------------------------------------------------------------
# statistics
#
# -s - tells bcftools stats to generate per-sample statistics
# for all samples.
# ------------------------------------------------------------
bcftools stats \
    -s - \
    "$OUT" \
    > "$STATS"


# ------------------------------------------------------------
# validate number of samples
# ------------------------------------------------------------
NSAMPLES=$(bcftools query -l "$OUT" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: Filtered VCF contains $NSAMPLES samples"
    exit 1
fi


# ------------------------------------------------------------
# obtain raw and filtered record counts from stats
# ------------------------------------------------------------
if [[ -s "$RAWSTATS" ]]; then

    RAW_N=$(
        awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "$RAWSTATS"
    )

else

    RAW_N="NA"

fi


FILTERED_N=$(
    awk -F'\t' \
    '$1=="SN" && $3=="number of records:" {print $4}' \
    "$STATS"
)


# ------------------------------------------------------------
# confirm output contains only biallelic SNPs
# ------------------------------------------------------------
NON_SNP=$(
    bcftools view \
        -H \
        -V snps \
        "$OUT" |
    wc -l
)

MULTI=$(
    bcftools view \
        -H \
        -m3 \
        "$OUT" |
    wc -l
)

if [[ "$NON_SNP" -ne 0 ]]; then
    echo "ERROR: $NON_SNP non-SNP records remain"
    exit 1
fi

if [[ "$MULTI" -ne 0 ]]; then
    echo "ERROR: $MULTI multiallelic records remain"
    exit 1
fi


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Completed chromosome:   $CHR"
echo "Raw records:            $RAW_N"
echo "Stage-1 records:        $FILTERED_N"
echo "Samples:                $NSAMPLES"
echo "Non-SNP records:        $NON_SNP"
echo "Multiallelic records:   $MULTI"
echo "Output:                 $OUT"
echo "Finish:                 $(date)"
echo "============================================================"
```

Make a summary file of the output:
```sh
printf "chromosome\trecords\tsnps\tmultiallelic\n" \
    > 04_vcf/filtered_stage1_stats/stage1_variant_summary.tsv

for f in 04_vcf/filtered_stage1_stats/*.filtered_stage1.stats.txt; do

    chr=$(basename "$f" .filtered_stage1.stats.txt)

    records=$(awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' "$f")

    snps=$(awk -F'\t' \
        '$1=="SN" && $3=="number of SNPs:" {print $4}' "$f")

    multi=$(awk -F'\t' \
        '$1=="SN" && $3=="number of multiallelic sites:" {print $4}' "$f")

    printf "%s\t%s\t%s\t%s\n" \
        "$chr" "$records" "$snps" "$multi"

done >> 04_vcf/filtered_stage1_stats/stage1_variant_summary.tsv
```

View the summary file:
```sh
column -t -s $'\t' \
    04_vcf/filtered_stage1_stats/stage1_variant_summary.tsv
```

Get totals:
```sh
awk '
NR>1 {
    records += $2
    snps += $3
    multi += $4
}
END {
    printf "Stage-1 records:            %d\n",records
    printf "Stage-1 SNPs:               %d\n",snps
    printf "Multiallelic sites:         %d\n",multi
}' 04_vcf/filtered_stage1_stats/stage1_variant_summary.tsv
```

Next, quantify missing genotypes:
```sh
printf "sample\tcalled_genotypes\tmissing_genotypes\ttotal_genotypes\tmissing_fraction\n" \
    > 04_vcf/filtered_stage1_stats/per_sample_missingness.tsv

for sample in $(cat mainland_samples.txt); do

    called=0
    missing=0

    for vcf in 04_vcf/filtered_stage1/*.filtered_stage1.vcf.gz; do

        vals=$(
            bcftools query \
                -s "$sample" \
                -f '[%GT\n]' \
                "$vcf" |
            awk '
            {
                total++

                if ($1=="./." || $1==".|." || $1==".")
                    missing++
                else
                    called++
            }
            END {
                print called+0,missing+0,total+0
            }'
        )

        c=$(echo "$vals" | awk '{print $1}')
        m=$(echo "$vals" | awk '{print $2}')

        called=$((called + c))
        missing=$((missing + m))

    done

    total=$((called + missing))

    awk \
        -v sample="$sample" \
        -v called="$called" \
        -v missing="$missing" \
        -v total="$total" \
        'BEGIN {
            printf "%s\t%d\t%d\t%d\t%.6f\n",
                   sample,called,missing,total,missing/total
        }'

done >> 04_vcf/filtered_stage1_stats/per_sample_missingness.tsv
```

View the file:
```sh
column -t -s $'\t' \
    04_vcf/filtered_stage1_stats/per_sample_missingness.tsv
```

This will show that per-sample missingness is only about ~4.2% across 21,177,169 SNPs.

Now look at site missingness distribution. This shows how many SNPs have 0, 1, 2, etc. missing samples.
```sh
printf "n_missing_samples\tn_sites\n" \
    > 04_vcf/filtered_stage1_stats/site_missingness_distribution.tsv

for vcf in 04_vcf/filtered_stage1/*.filtered_stage1.vcf.gz; do

    bcftools query \
        -f '[%GT\t]\n' \
        "$vcf" |
    awk '
    {
        missing=0

        for (i=1; i<=NF; i++) {
            if ($i=="./." || $i==".|." || $i==".")
                missing++
        }

        count[missing]++
    }

    END {
        for (m in count)
            print m,count[m]
    }'

done |
awk '
{
    total[$1]+=$2
}
END {
    for (m=0; m<=8; m++)
        print m "\t" total[m]+0
}' \
>> 04_vcf/filtered_stage1_stats/site_missingness_distribution.tsv
```

View the file:
```sh
column -t -s $'\t' \
    04_vcf/filtered_stage1_stats/site_missingness_distribution.tsv
```

The output should look like this:
```sh
n_missing_samples  n_sites
0                  17630051
1                  2362631
2                  659589
3                  275319
4                  129912
5                  64057
6                  32844
7                  16485
8                  6281
```

This is a very clean distribution. Keeping sites with 0 or 1 missing sample and removing sites with 2 or more missing samples could be a reasonable rule for stage 2 SNP filtering. This also means that the second filtering step will require at least 7 of 8 individuals called at each SNP. So stage 2 filtering would retain 94.41% of the stage 1 SNPs.

### Step 10: VCF filtering - Stage 2 (site missingness filtering)
This step asks whether the SNPs retained from the first filtering stage still have enough individuals with reliable data.

Use the script below for filtering step 2:
```sh
#!/bin/bash
#SBATCH --job-name=smc_vcf_filt_2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ================================================================
# SMC++ VCF filtering - Stage 2
#
# starting data:
#   stage-1 filtered autosomal VCFs
#
# stage-1 already applied:
#   shared callable genome
#   SNPs only
#   biallelic sites only
#   QUAL >= 20
#   genotype DP 5-25
#   genotype GQ >= 20
#   failed genotypes set to ./.
#
# stage-2 filter:
#   retain sites with no more than 1 missing sample
#
# therefore:
#   0 missing samples -> KEEP
#   1 missing sample  -> KEEP
#   >=2 missing       -> REMOVE
#
# this requires at least 7 of 8 samples called at each SNP.
# ================================================================


# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
INDIR="${WORKDIR}/04_vcf/filtered_stage1"
STAGE1STATDIR="${WORKDIR}/04_vcf/filtered_stage1_stats"
OUTDIR="${WORKDIR}/04_vcf/filtered_stage2"
STATDIR="${WORKDIR}/04_vcf/filtered_stage2_stats"

mkdir -p "$OUTDIR"
mkdir -p "$STATDIR"


# ------------------------------------------------------------
# get chromosome
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")


# ------------------------------------------------------------
# files
# ------------------------------------------------------------
IN="${INDIR}/${CHR}.filtered_stage1.vcf.gz"
OUT="${OUTDIR}/${CHR}.filtered_stage2.vcf.gz"
STAGE1STATS="${STAGE1STATDIR}/${CHR}.filtered_stage1.stats.txt"
STATS="${STATDIR}/${CHR}.filtered_stage2.stats.txt"

# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ VCF filtering - Stage 2"
echo "============================================================"
echo "Chromosome:       $CHR"
echo "Input:            $IN"
echo "Output:           $OUT"
echo "Missingness rule: N_MISSING <= 1"
echo "Minimum called:   7 of 8 samples"
echo "Array task:       ${SLURM_ARRAY_TASK_ID}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK}"
echo "Host:             $(hostname)"
echo "Start:            $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# stage-2 filtering
#
# N_MISSING is calculated by bcftools from the GT fields.
#
# keep sites where:
#
#     N_MISSING <= 1
#
# this retains sites with:
#     8/8 called
#     7/8 called
#
# and removes sites with:
#     <=6/8 called
# ------------------------------------------------------------
bcftools view \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -i 'N_MISSING<=1' \
    -Oz \
    -o "$OUT" \
    "$IN"


# ------------------------------------------------------------
# index output
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$OUT"


# ------------------------------------------------------------
# generate statistics
# ------------------------------------------------------------
bcftools stats \
    -s - \
    "$OUT" \
    > "$STATS"


# ------------------------------------------------------------
# validate number of samples
# ------------------------------------------------------------
NSAMPLES=$(bcftools query -l "$OUT" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: Stage-2 VCF contains $NSAMPLES samples"
    exit 1
fi


# ------------------------------------------------------------
# get Stage-1 record count
# ------------------------------------------------------------
if [[ -s "$STAGE1STATS" ]]; then

    STAGE1_N=$(
        awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "$STAGE1STATS"
    )

else

    STAGE1_N="NA"

fi


# ------------------------------------------------------------
# get Stage-2 record count
# ------------------------------------------------------------
STAGE2_N=$(
    awk -F'\t' \
    '$1=="SN" && $3=="number of records:" {print $4}' \
    "$STATS"
)


# ------------------------------------------------------------
# confirm no sites with >=2 missing samples remain
# ------------------------------------------------------------
BAD_MISSING=$(
    bcftools view \
        -H \
        -i 'N_MISSING>1' \
        "$OUT" |
    wc -l
)

if [[ "$BAD_MISSING" -ne 0 ]]; then
    echo "ERROR: $BAD_MISSING sites with >1 missing sample remain"
    exit 1
fi


# ------------------------------------------------------------
# confirm output still contains only SNPs
# ------------------------------------------------------------
NON_SNP=$(
    bcftools view \
        -H \
        -V snps \
        "$OUT" |
    wc -l
)

if [[ "$NON_SNP" -ne 0 ]]; then
    echo "ERROR: $NON_SNP non-SNP records remain"
    exit 1
fi


# ------------------------------------------------------------
# confirm output remains biallelic
# ------------------------------------------------------------
MULTI=$(
    bcftools view \
        -H \
        -m3 \
        "$OUT" |
    wc -l
)

if [[ "$MULTI" -ne 0 ]]; then
    echo "ERROR: $MULTI multiallelic records remain"
    exit 1
fi


# ------------------------------------------------------------
# calculate number and percent removed
# ------------------------------------------------------------
if [[ "$STAGE1_N" != "NA" && "$STAGE1_N" -gt 0 ]]; then

    REMOVED=$((STAGE1_N - STAGE2_N))

    PCT_RETAINED=$(
        awk \
            -v kept="$STAGE2_N" \
            -v total="$STAGE1_N" \
            'BEGIN {printf "%.2f",100*kept/total}'
    )

    PCT_REMOVED=$(
        awk \
            -v removed="$REMOVED" \
            -v total="$STAGE1_N" \
            'BEGIN {printf "%.2f",100*removed/total}'
    )

else

    REMOVED="NA"
    PCT_RETAINED="NA"
    PCT_REMOVED="NA"

fi


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Completed chromosome:     $CHR"
echo "Stage-1 records:          $STAGE1_N"
echo "Stage-2 records:          $STAGE2_N"
echo "Removed:                  $REMOVED"
echo "Percent retained:         ${PCT_RETAINED}%"
echo "Percent removed:          ${PCT_REMOVED}%"
echo "Samples:                  $NSAMPLES"
echo "Sites with >1 missing:    $BAD_MISSING"
echo "Non-SNP records:          $NON_SNP"
echo "Multiallelic records:     $MULTI"
echo "Output:                   $OUT"
echo "Finish:                   $(date)"
echo "============================================================"
```

Summarize the retained SNPs:
```sh
printf "chromosome\tstage1_snps\tstage2_snps\tremoved\tpercent_retained\n" \
    > 04_vcf/filtered_stage2_stats/stage2_summary.tsv

for chr in $(cat 01_reference/autosomes.txt); do

    s1=$(
        awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "04_vcf/filtered_stage1_stats/${chr}.filtered_stage1.stats.txt"
    )

    s2=$(
        awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "04_vcf/filtered_stage2_stats/${chr}.filtered_stage2.stats.txt"
    )

    awk \
        -v chr="$chr" \
        -v s1="$s1" \
        -v s2="$s2" \
        'BEGIN {
            printf "%s\t%d\t%d\t%d\t%.2f\n",
                   chr,s1,s2,s1-s2,100*s2/s1
        }'

done >> 04_vcf/filtered_stage2_stats/stage2_summary.tsv
```

View the file:
```sh
column -t -s $'\t' \
    04_vcf/filtered_stage2_stats/stage2_summary.tsv
```

Get genome-wide total:
```sh
awk '
NR>1 {
    s1 += $2
    s2 += $3
}
END {
    printf "Stage-1 SNPs:       %d\n",s1
    printf "Stage-2 SNPs:       %d\n",s2
    printf "Removed:            %d\n",s1-s2
    printf "Percent retained:   %.2f%%\n",100*s2/s1
}' 04_vcf/filtered_stage2_stats/stage2_summary.tsv
```

The results will show that 19,992,682 biallelic SNPs were retained, or 94.41% of the SNPs identified in the first filtering stage.


### Step 11: Final VCF QC
We will now concatenate the 17 stage-2 chromosome VCFs into one final autosomal VCF and conduct QC.
```sh
#!/bin/bash
#SBATCH --job-name=smc_vcf_concat
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%j.err

# ================================================================
# Concatenate Stage-2 autosomal VCFs and generate final QC stats
# ================================================================

source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate smc_tools

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"

CHRFILE="${WORKDIR}/01_reference/autosomes.txt"

INDIR="${WORKDIR}/04_vcf/filtered_stage2"

OUTDIR="${WORKDIR}/04_vcf/final"
QCDIR="${WORKDIR}/04_vcf/final_qc"

VCFLIST="${OUTDIR}/stage2_vcf_list.txt"

FINALVCF="${OUTDIR}/G_ussuriensis_mainland.autosomes.filtered.vcf.gz"

STATS="${QCDIR}/G_ussuriensis_mainland.autosomes.filtered.stats.txt"

mkdir -p "$OUTDIR"
mkdir -p "$QCDIR"


# ------------------------------------------------------------
# create VCF list in chromosome order
# ------------------------------------------------------------
> "$VCFLIST"

while read -r chr; do

    VCF="${INDIR}/${chr}.filtered_stage2.vcf.gz"

    if [[ ! -s "$VCF" ]]; then
        echo "ERROR: Missing VCF:"
        echo "$VCF"
        exit 1
    fi

    if [[ ! -s "${VCF}.tbi" && ! -s "${VCF}.csi" ]]; then
        echo "ERROR: Missing VCF index:"
        echo "$VCF"
        exit 1
    fi

    echo "$VCF" >> "$VCFLIST"

done < "$CHRFILE"


# ------------------------------------------------------------
# verify 17 chromosome VCFs
# ------------------------------------------------------------
NVCF=$(wc -l < "$VCFLIST")

if [[ "$NVCF" -ne 17 ]]; then
    echo "ERROR: Expected 17 VCFs but found $NVCF"
    exit 1
fi


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo "============================================================"
echo "SMC++ final autosomal VCF"
echo "============================================================"
echo "Number of chromosomes: $NVCF"
echo "Output:                $FINALVCF"
echo "Start:                 $(date)"
echo "============================================================"


# ------------------------------------------------------------
# concatenate chromosomes
# ------------------------------------------------------------
bcftools concat \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f "$VCFLIST" \
    -Oz \
    -o "$FINALVCF"


# ------------------------------------------------------------
# index
# ------------------------------------------------------------
bcftools index \
    --threads "${SLURM_CPUS_PER_TASK}" \
    -f \
    -t \
    "$FINALVCF"


# ------------------------------------------------------------
# final statistics including per-sample statistics
# ------------------------------------------------------------
bcftools stats \
    -s - \
    "$FINALVCF" \
    > "$STATS"


# ------------------------------------------------------------
# validation
# ------------------------------------------------------------
NSAMPLES=$(bcftools query -l "$FINALVCF" | wc -l)

if [[ "$NSAMPLES" -ne 8 ]]; then
    echo "ERROR: Expected 8 samples but found $NSAMPLES"
    exit 1
fi


NRECORDS=$(
    awk -F'\t' \
        '$1=="SN" && $3=="number of records:" {print $4}' \
        "$STATS"
)


NON_SNP=$(
    bcftools view \
        -H \
        -V snps \
        "$FINALVCF" |
    wc -l
)


MULTI=$(
    bcftools view \
        -H \
        -m3 \
        "$FINALVCF" |
    wc -l
)


BAD_MISSING=$(
    bcftools view \
        -H \
        -i 'N_MISSING>1' \
        "$FINALVCF" |
    wc -l
)


# ------------------------------------------------------------
# final report
# ------------------------------------------------------------
echo
echo "============================================================"
echo "FINAL VCF VALIDATION"
echo "============================================================"
echo "Records:                  $NRECORDS"
echo "Samples:                  $NSAMPLES"
echo "Non-SNP records:          $NON_SNP"
echo "Multiallelic records:     $MULTI"
echo "Sites with >1 missing:    $BAD_MISSING"
echo "Final VCF:                $FINALVCF"
echo "Stats:                    $STATS"
echo "Finish:                   $(date)"
echo "============================================================"
```

Now, verify the concatenated VCF and run final QC before SMC++.
```sh
grep -A20 "FINAL VCF VALIDATION" \
    outfiles/slurm-smc_vcf_concat_*.out
```

The output should print:
```sh
FINAL VCF VALIDATION
============================================================
Records:                  19992682
Samples:                  8
Non-SNP records:          0
Multiallelic records:     0
Sites with >1 missing:    0
Final VCF:                /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/04_vcf/final/G_ussuriensis_mainland.autosomes.filtered.vcf.gz
Stats:                    /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/04_vcf/final_qc/G_ussuriensis_mainland.autosomes.filtered.stats.txt
Finish:                   Mon Aug 24 19:52:01 EDT 2026
============================================================
```

Confirm the final VCF and index exists:
```sh
ls -lh \
    04_vcf/final/G_ussuriensis_mainland.autosomes.filtered.vcf.gz*
```

Now use the bcftools stats file that the job already generated and get a series of summary values:
```sh
# define path to file
STATS="04_vcf/final_qc/G_ussuriensis_mainland.autosomes.filtered.stats.txt"

# overall summary
grep '^SN' "$STATS"

# transition/transversion summary
grep '^TSTV' "$STATS"

# per-sample QC
grep '^# PSC' "$STATS"
grep '^PSC' "$STATS"
```

Finally, get SNP counts per chromosome directly from the concatenated VCF:
```sh
printf "chromosome\tsnps\n" \
    > 04_vcf/final_qc/snps_per_chromosome.tsv

while read -r chr; do

    n=$(
        bcftools view \
            -H \
            -r "$chr" \
            04_vcf/final/G_ussuriensis_mainland.autosomes.filtered.vcf.gz |
        wc -l
    )

    printf "%s\t%d\n" "$chr" "$n"

done < 01_reference/autosomes.txt \
    >> 04_vcf/final_qc/snps_per_chromosome.tsv
```

Look at the file:
```sh
column -t -s $'\t' \
    04_vcf/final_qc/snps_per_chromosome.tsv
```

Our samples look very good as a SMC++ input:
```sh
| Sample     | Avg. DP | Heterozygous SNPs | Missing SNP genotypes |
| ---------- | ------: | ----------------: | --------------------: |
| AMNH_21010 |   13.5× |         4,239,268 |                 1.30% |
| AMNH_21128 |   13.3× |         4,932,001 |                 1.46% |
| AMNH_21147 |   12.9× |         4,800,337 |                 1.94% |
| AMNH_21161 |   13.6× |         4,831,175 |                 1.44% |
| AMNH_21162 |   14.1× |         4,818,188 |                 1.24% |
| AMNH_21164 |   13.7× |         4,676,233 |                 1.53% |
| AMNH_21172 |   14.9× |         4,719,096 |                 1.08% |
| AMNH_21185 |   13.3× |         4,473,157 |                 1.82% |
```

Here's a summary of input samples:
```sh
1) 8 mainland individuals
2) 17 autosomes
3) 1.380 Gb total autosomal reference
4) 953.3 Mb shared callable sequence
5) 69.07% callable
6) 19,992,682 high-quality biallelic SNPs
7) QUAL >= 20
8) MAPQ >= 30 during calling
9) base quality >= 20
10) genotype: DP 5–25; GQ >= 20
11) site: at least 7/8 individuals called
12) Ti/Tv = 2.42
```

### Step 12: VCF to SMC file conversion
First, we need to convert the final VCF file into SMC++ input file, which can be done using smc++ vcf2smc command.

```sh
#!/bin/bash
#SBATCH --job-name=smc_vcf2smc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-17%4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography/outfiles/slurm-%x_%A_%a.err

# ================================================================
# convert filtered autosomal VCF to SMC++ input
#
# one autosome per SLURM array task.
#
# population: Mainland G. ussuriensis
#
# samples: 8 mainland individuals
# 
# distinguished individual (-d): AMNH_21172
#    this sample has excellent coverage and the lowest final missingness
#
# because the VCF is unphased, both distinguished lineages
# are taken from the same individual:
#
#   -d AMNH_21172 AMNH_21172
#
# the depth/callability exclusion mask is supplied with --mask.
# ================================================================

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
WORKDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography"
SIF="${WORKDIR}/00_env/smcpp_latest.sif"
VCF="${WORKDIR}/04_vcf/final/G_ussuriensis_mainland.autosomes.filtered.vcf.gz"
CHRFILE="${WORKDIR}/01_reference/autosomes.txt"
MASK="${WORKDIR}/05_masks/smcpp_exclude.primary.autosomes.sorted.bed"
MASKDIR="${WORKDIR}/05_masks/by_chrom"
OUTDIR="${WORKDIR}/06_smc/primary"

mkdir -p "$MASKDIR"
mkdir -p "$OUTDIR"


# ------------------------------------------------------------
# samples
# ------------------------------------------------------------
POP="Mainland"
SAMPLES="AMNH_21010,AMNH_21128,AMNH_21147,AMNH_21161,AMNH_21162,AMNH_21164,AMNH_21172,AMNH_21185"
DIST="AMNH_21172"


# ------------------------------------------------------------
# load Apptainer
# ------------------------------------------------------------
module load Apptainer/apptainer-1.2.5


# ------------------------------------------------------------
# get chromosome
# ------------------------------------------------------------
CHR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CHRFILE")


# ------------------------------------------------------------
# chromosome-specific mask
#
# SMC++ processes one chromosome at a time, so extract only
# the mask intervals belonging to the current chromosome.
# ------------------------------------------------------------
CHRMASK="${MASKDIR}/${CHR}.smcpp_exclude.bed"

awk -v chr="$CHR" '
BEGIN {OFS="\t"}
$1==chr {
    print $1,$2,$3
}
' "$MASK" > "$CHRMASK"


# ------------------------------------------------------------
# output
# ------------------------------------------------------------
OUT="${OUTDIR}/${CHR}.${DIST}.smc.gz"


# ------------------------------------------------------------
# job information
# ------------------------------------------------------------
echo
echo "============================================================"
echo "SMC++ vcf2smc"
echo "============================================================"
echo "Chromosome:             $CHR"
echo "Population:             $POP"
echo "Distinguished sample:   $DIST"
echo "VCF:                    $VCF"
echo "Mask:                   $CHRMASK"
echo "Container:              $SIF"
echo "Output:                 $OUT"
echo "Host:                   $(hostname)"
echo "Array task:             ${SLURM_ARRAY_TASK_ID}"
echo "Start:                  $(date)"
echo "============================================================"
echo


# ------------------------------------------------------------
# report chromosome mask size
# ------------------------------------------------------------
MASKED_BP=$(
    awk '
    {
        bp += $3-$2
    }
    END {
        print bp+0
    }' "$CHRMASK"
)

echo "Masked bp on ${CHR}: $MASKED_BP"
echo


# ------------------------------------------------------------
# verify container
# ------------------------------------------------------------
echo "Checking SMC++ container..."

apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$SIF" \
    smc++ --help \
    >/dev/null


# ------------------------------------------------------------
# convert VCF to SMC++ format
#
# -d DIST DIST:
#   for unphased data, both distinguished lineages are drawn
#   from the same diploid individual.
#
# --mask:
#   positions in this BED are marked missing rather than being
#   interpreted as long homozygous sequence.
# ------------------------------------------------------------
apptainer exec \
    --cleanenv \
    --bind "${WORKDIR}:${WORKDIR}" \
    "$SIF" \
    smc++ vcf2smc \
    -d "$DIST" "$DIST" \
    --mask "$CHRMASK" \
    "$VCF" \
    "$OUT" \
    "$CHR" \
    "${POP}:${SAMPLES}"


# ------------------------------------------------------------
# finish
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Completed chromosome:   $CHR"
echo "Distinguished sample:   $DIST"
echo "Output:                 $OUT"
echo "Output size:            $(du -h "$OUT" | cut -f1)"
echo "Finish:                 $(date)"
echo "============================================================"
```