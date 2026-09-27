# 11) Demographic history using PSMC
We will use PSMC to infer the demographic history of *G. ussuriensis.* Our reference individual (AMNH 21010) has both chromosome-level assembly and Illumina reads. PSMC requires a whole-genome diploid consensus sequence from one individual as an input. The assembly fasta largely represents one sequence state at heterozygous sites. PSMC needs the spatial pattern of heterozygous versus homozygous regions along the genome. This is where the Illumina reads from the same individual is useful. By mapping the Illumina WGS reads back to the reference assembly, we can identify homozygous vs heterozygous sites, and from there we can construct the diploid consensus needed for PSMC.

### Step 1: Setup
Set up the basic directory structure.
```sh
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"

mkdir -p "${PSMCROOT}"/{00_setup,01_depth_qc,02_consensus,03_psmc,04_bootstrap,05_plots,scripts,slurm_logs}
``` 
... and point to existing input files:
```sh
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"

REF="${PROJECT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

BAM="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam/AMNH_21010.markdup.bam"
BAI="${BAM}.bai"

REGIONS="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/metadata/autosome_regions.txt"
```

Let's create a conda env for PSMC and install software we will need:
```sh
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh

conda create \
    -n psmc \
    -c conda-forge \
    -c bioconda \
    --strict-channel-priority \
    psmc \
    samtools \
    bcftools \
    htslib \
    seqtk \
    bedtools \
    perl \
    gnuplot \
    mosdepth \
    -y

conda activate psmc
```
Then, check the installation:
```sh
which psmc
which fq2psmcfa
which splitfa
which psmc_plot.pl
which samtools
which bcftools
which seqtk
which bedtools
```
Also, save the env information:
```sh
conda env export --no-builds \
    > "${PSMCROOT}/00_setup/psmc_environment.yml"

conda list \
    > "${PSMCROOT}/00_setup/psmc_conda_packages.txt"
```

Verify the input files:
```sh
ls -lh \
    "${REF}" \
    "${BAM}" \
    "${BAI}" \
    "${REGIONS}"
```

Copy the autosome definition file to the "demography_PSMC/00_setup" directory and create the plain chromosome-name list: 
```sh
cp "${REGIONS}" ./autosome_regions.txt
sed 's/:$//' autosome_regions.txt > autosome_contigs.txt
```

There should be 17 defined chromosomes:
```sh
cat autosome_contigs.txt
wc -l autosome_contigs.txt
```

Check whether the reference already has a fasta index:
```sh
ls -lh "${REF}.fai"
```

Make sure that all 17 autosome names exist in the reference:
```sh
while read chr; do
    if grep -q -w "^${chr}" "${REF}.fai"; then
        echo "FOUND: ${chr}"
    else
        echo "MISSING: ${chr}"
    fi
done < autosome_contigs.txt
```

Now, generate the autosome .bed file:
```sh
awk 'NR==FNR {keep[$1]=1; next} ($1 in keep) {print $1 "\t0\t" $2}' \
    autosome_contigs.txt \
    "${REF}.fai" \
    > autosomes.bed
```
... and calculate the total autosomal sequence length:
```SH
awk '{
    sum += $3-$2
}
END {
    printf "Autosomes: %d\n", NR
    printf "Total autosomal bases: %.0f\n", sum
    printf "Total autosomal Mb: %.3f\n", sum/1e6
    printf "Total autosomal Gb: %.3f\n", sum/1e9
}' autosomes.bed
```

### Step 2: Calculate autosomal depth
Use mosdepth to calculate autosomal depth. With a BED supplied via --by flag, this will produce a region-specific depth distribution needed for PSMC.
```sh
#!/bin/bash
#SBATCH --job-name=psmc_mosdepth
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda environment
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate psmc

set -euo pipefail

# setup paths
PSMCROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC"

BAM="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam/AMNH_21010.markdup.bam"
AUTOSOMES="${PSMCROOT}/00_setup/autosomes.bed"
OUTDIR="${PSMCROOT}/01_depth_qc"

# create output directory and cd into it
mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# run mosdepth to calculate depth of coverage
# -n          do not write gigantic per-base output
# -t 8        use 8 BAM decompression threads
# -Q 30       require MAPQ ≥ 30
# -F 3844     exclude unmapped, secondary, QC-fail, duplicate, and supplementary alignments
# --by BED    calculate statistics specifically for our 17 autosomes

mosdepth \
    -n \
    -t 8 \
    -Q 30 \
    -F 3844 \
    --by "${AUTOSOMES}" \
    AMNH_21010.autosomes \
    "${BAM}"
```

Run these after the job finishes running:
```sh
cat AMNH_21010.autosomes.mosdepth.summary.txt
head -20 AMNH_21010.autosomes.mosdepth.region.dist.txt
tail -20 AMNH_21010.autosomes.mosdepth.region.dist.txt
```
Also run this to extract the overall autosomal mean depth directly from the summary:
```sh
awk '$1=="total_region" {print "Mean autosomal depth:", $4}' \
    AMNH_21010.autosomes.mosdepth.summary.txt
```

This will show that the mean autosomal coverage is about ~12x. Here are other key numbers to keep in mind for setting analysis parameters (e.g. depth filters, etc.):
```sh
- Mean autosomal depth: 11.88×
- Median depth: approximately 12× (≥12× = 53%, ≥13× = 46%)
- Bases with any coverage: about 97%
- Bases at ≥4×: about 92%
- Bases at ≥5×: about 90%
```

PSMC documentation recommends to set the min depth (-d) to 1/3 of the average depth and max depth (-D) to 2x of average depth. Therefore, our parameters would become:
```sh
-d = 4
-D = 24
```

### Step 3: Build diploid consensus sequence
Run the script below:
```sh
#!/bin/bash
#SBATCH --job-name=psmc_consensus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# ============================================================================
# generate diploid consensus FASTQ for PSMC
# individual: AMNH_21010
#
# input: Illumina WGS mapped to the AMNH_21010 chromosome-level assembly
# regions: 17 autosomes only
# filtering:
#   minimum mapping quality = 30
#   minimum base quality    = 20
#   minimum depth           = 4x
#   maximum depth           = 24x
#
# output: diploid consensus FASTQ for downstream fq2psmcfa
# ============================================================================


# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate psmc

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"

REF="${PROJECT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"
BAM="/home/yshin/mendel-nas1/ussuri_popgen/WGS_mapping/bam/AMNH_21010.markdup.bam"
AUTOSOMES="${PSMCROOT}/00_setup/autosomes.bed"
OUTDIR="${PSMCROOT}/02_consensus"

mkdir -p "${OUTDIR}"


# ------------------------------------------------------------
# parameters
# ------------------------------------------------------------
SAMPLE="AMNH_21010"

MIN_MQ=30
MIN_BQ=20

MIN_DP=4
MAX_DP=24

CONSENSUS_MQ=30

# small region used to verify that current bcftools output
# is compatible with vcfutils.pl vcf2fq
TEST_REGION="G_ussuri_chr18:1-1000000"


# ------------------------------------------------------------
# output names
# ------------------------------------------------------------
PREFIX="${OUTDIR}/${SAMPLE}.autosomes.dp${MIN_DP}-${MAX_DP}"
TEST_VCF="${OUTDIR}/${SAMPLE}.preflight.vcf"
TEST_FQ="${OUTDIR}/${SAMPLE}.preflight.fq"
DIPLOID_FQ="${PREFIX}.diploid.fq.gz"


# ------------------------------------------------------------
# check inputs
# ------------------------------------------------------------

echo "============================================================"
echo "PSMC diploid consensus generation"
echo "============================================================"
echo "Sample:        ${SAMPLE}"
echo "Reference:     ${REF}"
echo "BAM:           ${BAM}"
echo "Autosomes:     ${AUTOSOMES}"
echo "MAPQ cutoff:   ${MIN_MQ}"
echo "BaseQ cutoff:  ${MIN_BQ}"
echo "Depth range:   ${MIN_DP}-${MAX_DP}x"
echo "============================================================"
echo

for FILE in "${REF}" "${REF}.fai" "${BAM}" "${BAM}.bai" "${AUTOSOMES}"; do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: missing or empty input: ${FILE}" >&2
        exit 1
    fi
done

for PROGRAM in bcftools samtools vcfutils.pl fq2psmcfa; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} not found in PATH" >&2
        exit 1
    fi
done

samtools quickcheck -v "${BAM}"


# ------------------------------------------------------------
# record software versions
# ------------------------------------------------------------
{
    echo "bcftools"
    bcftools --version | head -2
    echo
    echo "samtools"
    samtools --version | head -2
    echo
    echo "vcfutils.pl"
    command -v vcfutils.pl
    echo
    echo "fq2psmcfa"
    command -v fq2psmcfa
} > "${OUTDIR}/software_versions.txt"


# ------------------------------------------------------------
# preflight test
#
# vcfutils.pl vcf2fq was written around the original bcftools
# consensus-caller annotations. test a small autosomal region
# before processing the complete 1.38-Gb autosomal genome.
# ------------------------------------------------------------
echo
echo "Running preflight test on ${TEST_REGION} ..."
echo

bcftools mpileup \
    -f "${REF}" \
    -r "${TEST_REGION}" \
    -q "${MIN_MQ}" \
    -Q "${MIN_BQ}" \
    --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
    -Ou \
    "${BAM}" \
| bcftools call \
    -c \
    -Ov \
    -o "${TEST_VCF}"


# ------------------------------------------------------------
# verify single sample
# ------------------------------------------------------------
N_SAMPLES=$(bcftools query -l "${TEST_VCF}" | wc -l)

if [[ "${N_SAMPLES}" -ne 1 ]]; then
    echo "ERROR: expected one sample in VCF but found ${N_SAMPLES}" >&2
    bcftools query -l "${TEST_VCF}" >&2
    exit 1
fi

echo "VCF sample:"
bcftools query -l "${TEST_VCF}"


# ------------------------------------------------------------
# verify fields required by vcf2fq
# ------------------------------------------------------------
for TAG in DP MQ FQ; do
    if ! grep -q "##INFO=<ID=${TAG}," "${TEST_VCF}"; then
        echo "ERROR: INFO/${TAG} is absent from bcftools output." >&2
        echo "Current bcftools output is not directly compatible with the" >&2
        echo "legacy vcfutils.pl vcf2fq workflow." >&2
        exit 1
    fi
done


# ------------------------------------------------------------
# verify vcf2fq conversion
# ------------------------------------------------------------
vcfutils.pl vcf2fq \
    -d "${MIN_DP}" \
    -D "${MAX_DP}" \
    -Q "${CONSENSUS_MQ}" \
    "${TEST_VCF}" \
    > "${TEST_FQ}"

if [[ ! -s "${TEST_FQ}" ]]; then
    echo "ERROR: preflight diploid FASTQ is empty" >&2
    exit 1
fi

echo
echo "Preflight test passed."
echo


# ------------------------------------------------------------
# generate full autosomal diploid consensus
#
# IMPORTANT:
# Do not use bcftools call -v here.
# vcfutils.pl vcf2fq requires an all-sites VCF so that both
# homozygous and heterozygous regions are represented.
# ------------------------------------------------------------
echo "Generating full autosomal diploid consensus ..."
echo

bcftools mpileup \
    -f "${REF}" \
    -R "${AUTOSOMES}" \
    -q "${MIN_MQ}" \
    -Q "${MIN_BQ}" \
    --skip-any-set UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY \
    -Ou \
    "${BAM}" \
| bcftools call \
    -c \
    -Ov \
| vcfutils.pl vcf2fq \
    -d "${MIN_DP}" \
    -D "${MAX_DP}" \
    -Q "${CONSENSUS_MQ}" \
| gzip -c \
    > "${DIPLOID_FQ}"


# ------------------------------------------------------------
# validate output
# ------------------------------------------------------------
gzip -t "${DIPLOID_FQ}"

if [[ ! -s "${DIPLOID_FQ}" ]]; then
    echo "ERROR: diploid consensus FASTQ is empty" >&2
    exit 1
fi


# ------------------------------------------------------------
# report
# ------------------------------------------------------------
echo
echo "============================================================"
echo "Consensus generation completed"
echo "============================================================"
echo "Output:"
echo "${DIPLOID_FQ}"
echo
ls -lh "${DIPLOID_FQ}"
echo "============================================================"
```

Check the output:
```sh
# check file integrity
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC/02_consensus

FQ="AMNH_21010.autosomes.dp4-24.diploid.fq.gz"

ls -lh "${FQ}"
gzip -t "${FQ}" && echo "gzip integrity: OK"

# make sure all 17 autosomes are present
seqtk seq -A "${FQ}" \
    | grep '^>' \
    | sed 's/^>//'

seqtk seq -A "${FQ}" \
    | grep -c '^>'
```

Also quantify total sequence, callable uppercase bases, heterozygous IUPAC bases, and masked lowercase sequence:
```sh
seqtk seq -A "${FQ}" \
| awk '
BEGIN {
    total=0
    acgt=0
    het=0
    lower=0
    upperN=0
}
!/^>/ {
    s=$0
    total += length(s)

    x=s
    acgt += gsub(/[ACGT]/, "", x)

    x=s
    het += gsub(/[MRWSYK]/, "", x)

    x=s
    lower += gsub(/[a-z]/, "", x)

    x=s
    upperN += gsub(/N/, "", x)
}
END {
    callable = acgt + het

    printf "Total bases:            %d\n", total
    printf "Uppercase A/C/G/T:      %d\n", acgt
    printf "Uppercase heterozygous: %d\n", het
    printf "Lowercase masked:       %d\n", lower
    printf "Uppercase N:            %d\n", upperN
    printf "Callable bases:         %d\n", callable
    printf "Callable fraction:      %.4f (%.2f%%)\n", callable/total, 100*callable/total
    printf "Het/callable:           %.8f\n", het/callable
}'
```

This will show:
```sh
Total bases:            1380174220
Uppercase A/C/G/T:      1217013984
Uppercase heterozygous: 5652966
Lowercase masked:       157417510
Uppercase N:            89760
Callable bases:         1222666950
Callable fraction:      0.8859 (88.59%)
Het/callable:           0.00462347
```

Now, let's make a BED of lowercase regions from the softmasked reference. This is needed because the soft-masked reference already marks repetitive sequence in lowercase, but the consensus-generation step does not reliably preserve that original lowercase mask.

The problem is that vcfutils.pl vcf2fq reconstructs the diploid consensus from the VCF rather than copying the original FASTA character-for-character. The lowercase block (= repeats) can emerge from consensus generation as callable uppercase sequence if it passes the depth/MQ criteria. In short, the vcf2fq code builds sequence from the VCF records and applies its own lowercase/uppercase logic based on quality and depth, not based on the original EarlGrey mask.

Let's run the script below:
```sh
# set paths
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"
REF="${PROJECT}/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# run python script
python /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/extract_lowercase_bed.py \
    "${REF}" \
    "${PSMCROOT}/00_setup/autosome_contigs.txt" \
    > "${PSMCROOT}/00_setup/autosome_repeat_mask.bed"
```

Now, check how much of the autosomes is softmasked:
```sh
awk '{
    n++
    bp += $3-$2
}
END {
    printf "Repeat intervals: %d\n", n
    printf "Repeat-masked bp: %d\n", bp
    printf "Repeat-masked Mb: %.3f\n", bp/1e6
    printf "Fraction of autosomes: %.2f%%\n", 100*bp/1380185643
}' "${PSMCROOT}/00_setup/autosome_repeat_mask.bed"
```

The script above should print:
```sh
Repeat intervals: 1324681
Repeat-masked bp: 635512154
Repeat-masked Mb: 635.512
Fraction of autosomes: 46.05%
``` 

Now, let's apply this mask to the diploid consensus:
```sh
cd "${PSMCROOT}/02_consensus"
FQ="AMNH_21010.autosomes.dp4-24.diploid.fq.gz"
MASK="${PSMCROOT}/00_setup/autosome_repeat_mask.bed"
OUT="AMNH_21010.autosomes.dp4-24.repeatmasked.diploid.fq.gz"

seqtk seq \
    -M "${MASK}" \
    "${FQ}" \
| gzip -c \
> "${OUT}"
``` 

Verify the output:
```sh
gzip -t "${OUT}" && echo "repeat-masked FASTQ: OK"
ls -lh "${OUT}"
```

Quantify the final usable sequence:
```sh
seqtk seq -A "${OUT}" \
| awk '
BEGIN {
    total=0
    acgt=0
    het=0
    lower=0
    upperN=0
}
!/^>/ {
    s=$0
    total += length(s)

    x=s
    acgt += gsub(/[ACGT]/, "", x)

    x=s
    het += gsub(/[MRWSYK]/, "", x)

    x=s
    lower += gsub(/[a-z]/, "", x)

    x=s
    upperN += gsub(/N/, "", x)
}
END {
    callable = acgt + het

    printf "Total bases:            %d\n", total
    printf "Uppercase A/C/G/T:      %d\n", acgt
    printf "Uppercase heterozygous: %d\n", het
    printf "Lowercase masked:       %d\n", lower
    printf "Uppercase N:            %d\n", upperN
    printf "Callable bases:         %d\n", callable
    printf "Callable fraction:      %.4f (%.2f%%)\n", callable/total, 100*callable/total
    printf "Het/callable:           %.8f\n", het/callable
}'
```
This will print the following:
```sh
Total bases:            1380174220
Uppercase A/C/G/T:      695456361
Uppercase heterozygous: 2867028
Lowercase masked:       681817017
Uppercase N:            33814
Callable bases:         698323389
Callable fraction:      0.5060 (50.60%)
Het/callable:           0.00410559
```

So the tradeoff is:
```sh
1) Without repeat mask:
callable = 88.59%
het/site = 0.00462
more sequence, but greater risk of repeat-driven false heterozygosity

2) With EarlGrey mask:
callable = 50.60%
het/site = 0.00411
less sequence, but cleaner uniquely interpretable sequence
```

### Step 4: Convert the repeat-masked diploid FASTQ into PSMCFA
```sh
# cd into dir
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC/03_psmc

# set paths
PSMCROOT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC"
FQ="${PSMCROOT}/02_consensus/AMNH_21010.autosomes.dp4-24.repeatmasked.diploid.fq.gz"
OUT="AMNH_21010.autosomes.dp4-24.repeatmasked.psmcfa"

# convert
fq2psmcfa -q20 "${FQ}" > "${OUT}"

# check the output
ls -lh "${OUT}"
grep -c '^>' "${OUT}"
```

Quantify the PSMC bins:
```sh
awk '
BEGIN {
    total=0
    t=0
    k=0
    n=0
    seqs=0
}
(/^>/) {
    seqs++
    next
}
{
    total += length($0)

    x=$0
    t += gsub(/T/, "", x)

    x=$0
    k += gsub(/K/, "", x)

    x=$0
    n += gsub(/N/, "", x)
}
END {
    callable=t+k

    printf "Sequences:              %d\n", seqs
    printf "Total 100-bp bins:      %d\n", total
    printf "T bins:                 %d\n", t
    printf "K bins:                 %d\n", k
    printf "N bins:                 %d\n", n
    printf "Callable bins:          %d\n", callable
    printf "Callable-bin fraction:  %.4f (%.2f%%)\n", callable/total, 100*callable/total
    printf "Heterozygous-bin frac:  %.6f (%.3f%%)\n", k/callable, 100*k/callable
}' "${OUT}"
```

This will print this:
```sh
Sequences:              17
Total 100-bp bins:      13801750
T bins:                 6230264
K bins:                 1970294
N bins:                 5601192
Callable bins:          8200558
Callable-bin fraction:  0.5942 (59.42%)
Heterozygous-bin frac:  0.240263 (24.026%)
```
Therefore, after repeat masking you still have 8,200,558 callable 100-bp bins, corresponding to roughly 820 Mb of usable PSMC sequence.

### Step 5: Run PSMC
```sh
#!/bin/bash
#SBATCH --job-name=psmc_exploratory
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=300G
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC/slurm_logs/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/demography_PSMC/slurm_logs/slurm-%x_%j.err

# ================================================================
# PSMC demographic inference
#
# individual: AMNH_21010
#
# input:
#   Repeat-masked diploid consensus
#   Autosomes only
#   Illumina depth filter = 4-24x
#
# initial parameterization:
#   -N25
#   -t15
#   -r5
#   -p "4+25*2+4+6"
#
# NOTE:
# This is the initial exploratory fit. The interval pattern will
# be evaluated after the run before bootstrapping.
# ================================================================


# ------------------------------------------------------------
# activate environment
# ------------------------------------------------------------
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate psmc

set -euo pipefail


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
PROJECT="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo"
PSMCROOT="${PROJECT}/demography_PSMC"
INDIR="${PSMCROOT}/03_psmc"
INPUT="${INDIR}/AMNH_21010.autosomes.dp4-24.repeatmasked.psmcfa"
OUTPUT="${INDIR}/AMNH_21010.autosomes.dp4-24.repeatmasked.psmc"


# ------------------------------------------------------------
# check input
# ------------------------------------------------------------
if [[ ! -s "${INPUT}" ]]; then
    echo "ERROR: missing or empty PSMCFA:"
    echo "${INPUT}"
    exit 1
fi


# ------------------------------------------------------------
# run PSMC
# ------------------------------------------------------------
echo "============================================================"
echo "Running PSMC"
echo "============================================================"
echo "Input:  ${INPUT}"
echo "Output: ${OUTPUT}"
echo

psmc \
    -N25 \
    -t15 \
    -r5 \
    -p "4+25*2+4+6" \
    -o "${OUTPUT}" \
    "${INPUT}"


# ------------------------------------------------------------
# validate output
# ------------------------------------------------------------
if [[ ! -s "${OUTPUT}" ]]; then
    echo "ERROR: PSMC output is empty"
    exit 1
fi

echo
echo "============================================================"
echo "PSMC completed"
echo "============================================================"

ls -lh "${OUTPUT}"

echo
echo "Final iteration:"
grep '^RD' "${OUTPUT}" | tail -1

echo
echo "Final theta/rho:"
grep '^TR' "${OUTPUT}" | tail -1

echo
echo "Final inferred recombination summary:"
grep 'n_recomb' "${OUTPUT}" | tail -1

echo "============================================================"
```