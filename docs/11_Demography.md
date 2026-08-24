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

#### Step 5: Set depth cutoff
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

#### Step 6: Create actual callable-region BEDs
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

#### Step 7: Make the depth-based SMC++ exclusion mask
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

#### Step 8: