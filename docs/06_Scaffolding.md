## 6) Scaffolding through Hi-C data incorporation
### Hi-C sequencing overview
Hi-C sequencing was done at Texas A&M Institute for Genome Sciences and Society (TIGGS) on three lanes of Illumina NovaSeq X Plus, using the same blood sample used for PacBio sequencing. The sequencing libraries were prepared using the Dovetail Omni-C kit.

### Setup
Let's create a directory for scaffolding and install YaHS, which is a scaffolding tool for Hi-C data.
```sh
# create a directory for scaffolding
mkdir -p /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding
```

Let's also create a conda env for scaffolding on Mendel and install YaHS.
```sh
# create conda env for scaffolding
conda create -n scaffolding
conda activate scaffolding

# install YaHS
conda install bioconda::yahs
```

### Combine sequencing reads across lanes
The sample was sequenced on three Illuminal lanes. So, there are six files total (two reads x three lanes). 
![alt text](/etc/hic_reads.png)

Let's combine reads across lanes.
```sh
# run under the scaffolding directory in Mendel
mkdir -p combined

cat \
  lane_1/26231TGS_Gloydius-ussuriensis-21010_S5_L001_R1_001.fastq.gz \
  lane_2/26231TGS_Gloydius-ussuriensis-21010_S5_L002_R1_001.fastq.gz \
  lane_3/26231TGS_Gloydius-ussuriensis-21010_S5_L005_R1_001.fastq.gz \
  > combined/Gloydius_ussuriensis_HiC_R1.fastq.gz

cat \
  lane_1/26231TGS_Gloydius-ussuriensis-21010_S5_L001_R2_001.fastq.gz \
  lane_2/26231TGS_Gloydius-ussuriensis-21010_S5_L002_R2_001.fastq.gz \
  lane_3/26231TGS_Gloydius-ussuriensis-21010_S5_L005_R2_001.fastq.gz \
  > combined/Gloydius_ussuriensis_HiC_R2.fastq.gz
```
Also note that it is not necessary to adapter/quality trim reads for Dovetail Omni-C libraries. Proceed to Hi-C mapping and scaffolding.

### Map Hi-C reads to the draft genome
The next step is to prepare the PacBio draft genome and map Hi-C reads to it. Let's symlink the draft directory into the scaffolding directory.

```sh
# under the scaffolding directory on Mendel
mkdir -p draft
ln -s /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa
```

To index the draft genome and to map Hi-C reads to it, we need to install BWA. Also install samtools.

```sh
# in the scaffolding conda env
conda install bioconda::bwa
conda install bioconda::samtools
```

Then, run the script below to index the draft genome.
```sh
#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# set directory
indir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/draft

# run bwa
bwa index ${indir}/Gloydius_ussuriensis_AMNH_21010_noMito.fa
samtools faidx ${indir}/Gloydius_ussuriensis_AMNH_21010_noMito.fa
```

After the above script finishes running, map Hi-C reads to the draft using bwa-mem.
```sh
#!/bin/bash
#SBATCH --job-name=hic_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=300G
#SBATCH --time=72:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# make sure the job will stop if any step fails
set -euo pipefail

# set directories
indir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined"
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding"
tmpdir="${outdir}/tmp_hic_map_${SLURM_JOB_ID}"

# make directories for temp files
sort1_tmp="${tmpdir}/sort1"
namesort_tmp="${tmpdir}/namesort"
coordsort_tmp="${tmpdir}/coordsort"

mkdir -p "${sort1_tmp}" "${namesort_tmp}" "${coordsort_tmp}"

# make a directory for bwa logs
mkdir -p "${outdir}/logs"

# set variables
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
R1="${indir}/Gloydius_ussuriensis_HiC_R1.fastq.gz"
R2="${indir}/Gloydius_ussuriensis_HiC_R2.fastq.gz"
BAM_PREFIX="${outdir}/Gloydius_ussuriensis_HiC_to_draft"

# basic checks
echo "Starting BWA mapping: $(date)"
echo "REF = ${REF}"
echo "R1  = ${R1}"
echo "R2  = ${R2}"
echo "TMP = ${tmpdir}"

# quota/filesystem reporting
quota -s || true
df -h "${outdir}" || true

# flag if inputs are missing
[[ -s "${REF}" ]] || { echo "ERROR: missing reference ${REF}"; exit 1; }
[[ -s "${R1}" ]] || { echo "ERROR: missing R1 ${R1}"; exit 1; }
[[ -s "${R2}" ]] || { echo "ERROR: missing R2 ${R2}"; exit 1; }

# run bwa mem
bwa mem -5SP -t ${SLURM_CPUS_PER_TASK} "${REF}" "${R1}" "${R2}" 2> "${outdir}/logs/bwa_mem.log" | \
  samtools view -@ 16 -bS - | \
  samtools sort -@ 16 -m 4G -T "${sort1_tmp}/sort1" -o "${BAM_PREFIX}.sorted.bam" -

samtools index "${BAM_PREFIX}.sorted.bam"

# name-sort for fixmate
echo "Starting name-sort: $(date)"
ls -ld "${namesort_tmp}"
df -h "${namesort_tmp}" || true
quota -s || true

samtools sort -@ 16 -m 4G -T "${namesort_tmp}/namesort" -n \
  -o "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.sorted.bam"

# add mate information
samtools fixmate -@ 16 -m \
  "${BAM_PREFIX}.namesort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

# coordinate-sort again
samtools sort -@ 16 -m 4G -T "${coordsort_tmp}/coordsort" \
  -o "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.fixmate.bam"

# mark duplicates
samtools markdup -@ 16 \
  "${BAM_PREFIX}.fixmate.coordsort.bam" \
  "${BAM_PREFIX}.markdup.bam"

samtools index "${BAM_PREFIX}.markdup.bam"

# produce basic mapping summary
samtools flagstat "${BAM_PREFIX}.markdup.bam" > "${BAM_PREFIX}.markdup.flagstat.txt"

# check final BAM integrity
samtools quickcheck -v "${BAM_PREFIX}.markdup.bam"

# clean temp files only if everything succeeded
rm -rf "${tmpdir}"

echo "Finished Hi-C mapping: $(date)"
ls -lh "${BAM_PREFIX}.markdup.bam" "${BAM_PREFIX}.markdup.bam.bai" "${BAM_PREFIX}.markdup.flagstat.txt"
``` 
The above script has the following steps:
  1) bwa-mem maps the trimmed Hi-C reads to the draft HiFi genome and outputs alignments as SAM text
  2) samtools view reads that SAM from the pipe and converts it to BAM
  3) samtools sort coordinate-sorts the BAM and writes a sorted bam file
  4) samtools index creates an index for the coordinate-sorted BAM.
  5) samtools sort -n creates a name-sorted BAM (read pairs are placed next to each other)
  6) samtools fixmate adds/corrects mate-pair information needed for duplicate marking
  7) samtools sort coordinate-sorts the BAM again
  8) samtools markdup marks duplicate alignments in the BAM
  9) samtools index indexes the duplicate-marked BAM
  10) samtools flagstat summarizes the final duplicate-marked BAM

### Scaffolding with YaHS
After mapping is done, run scaffolding with YaHS. YaHS takes optional arguments, a contig fasta file (i.e., draft genome), and one of several possible Hi-C input format (e.g., bed, bam, pa5, bin).

```sh
#!/bin/bash
#SBATCH --job-name=yahs_hic
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# stop if anything fails
set -euo pipefail

# set paths
indir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding"
outdir="${indir}/yahs_out"

REF="${indir}/draft/Gloydius_ussuriensis_AMNH_21010_noMito.fa"
BAM="${indir}/Gloydius_ussuriensis_HiC_to_draft.markdup.bam"

# make yahs outdir
mkdir -p "$outdir"
cd "$outdir"

# run yahs
yahs "$REF" "$BAM" -o Gloydius_ussuriensis_AMNH_21010_yahs
```

Submit the above script with the following SLURM dependency (the script will start running once the mapping job is done).
```sh
sbatch --dependency=afterok:10822420 G_ussuri_YaHS.sh
```
Check the assembly stats in the log file after YaHS finishes running.


![alt text](etc/yahs_out.png)

We can see that the scaffold N50 = 203.1 Mb and L50 = 3, which is much more contiguous than our PacBio N50 = 127.5 Mb and L50 = 3.

The key output here is "Gloydius_ussuriensis_AMNH_21010_yahs_scaffolds_final.fa." Let's index this file. This will generate a ".fai" file.
```
# activate the scaffolding conda env first to access samtools
samtools faidx Gloydius_ussuriensis_AMNH_21010_yahs_scaffolds_final.fa
```
Let's inspect scaffold sizes.
```
cut -f1,2 Gloydius_ussuriensis_AMNH_21010_yahs_scaffolds_final.fa.fai | sort -k2,2nr | head -n 30
```
This will show that there are 9 very large scaffolds (57–341 Mb), several medium scaffolds (8–27 Mb), and many small leftover scaffolds.

To see if these scaffolds are correctly oriented and joined, we need to look at the Hi-C contact map.

### Hi-C contact map visualization with Juicer/Juicer Tools
Juicer is already installed together with YaHS. Also install Juicer Tools.
```
# in the scaffolding conda env
# install java
conda install -c conda-forge openjdk=11

# make dir for juicer tools
mkdir -p /home/yshin/mendel-nas1/juicer_tools

# install jucer tools
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar
```

Run juicer/juicer tools as below. One thing to note is that I'm feeding it the original PacBio draft .fai file. I did this because the job crashed when I gave it the scaffold .fai file because juicer looked for the original contig name. 
```sh
#!/bin/bash
#SBATCH --job-name=juicer
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=450G
#SBATCH --time=24:00:00
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# initiate conda and activate the conda environment
source ~/.bash_profile
conda activate scaffolding

# stop if anything fails
set -euo pipefail

# set dir
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/yahs_out"
prefix="Gloydius_ussuriensis_AMNH_21010_yahs"

# set vars
BIN="${workdir}/${prefix}.bin"
AGP="${workdir}/${prefix}_scaffolds_final.agp"
FASTA="${workdir}/${prefix}_scaffolds_final.fa"

# IMPORTANT: use the original PacBio draft FASTA index here,
# not the final YaHS scaffolded FASTA index
FAI="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/PacBio_Revio/no_mito/Gloydius_ussuriensis_AMNH_21010_noMito.fa.fai"

JBAT="${workdir}/${prefix}_JBAT"

# standalone juicer_tools jar
JUICER_TOOLS_JAR="/home/yshin/mendel-nas1/juicer_tools/juicer_tools_1.22.01.jar"

# cd into working dir
cd "$workdir"

# sanity checks
echo "Checking files: $(date)"
ls -lh "$BIN"
ls -lh "$AGP"
ls -lh "$FASTA"
ls -lh "$FAI"
ls -lh "$JUICER_TOOLS_JAR"

which juicer
which java
java -version
which samtools

echo "Using juicer tools jar:"
echo "$JUICER_TOOLS_JAR"

# remove previous incomplete hic only
rm -f "${JBAT}.hic"

# run juicer pre from YaHS
echo "Running juicer pre: $(date)"
juicer pre \
    -a \
    -o "$JBAT" \
    "$BIN" \
    "$AGP" \
    "$FAI" \
    > "${JBAT}.log" 2>&1

echo "Finished juicer pre: $(date)"

# make chromosome sizes file from juicer pre log
echo "Making chromosome sizes file from juicer pre log: $(date)"
grep PRE_C_SIZE "${JBAT}.log" | awk '{print $2, $3}' > "${JBAT}.chrom.sizes"

echo "Chromosome sizes:"
cat "${JBAT}.chrom.sizes"

# check JBAT outputs
echo "Checking JBAT outputs:"
ls -lh "${JBAT}"*

# make sure JBAT txt exists and has contacts
echo "Checking JBAT contact file:"
ls -lh "${JBAT}.txt"
head -n 5 "${JBAT}.txt"

# run juicer tools directly with java
echo "Creating .hic file: $(date)"

java -Xmx400G -jar "$JUICER_TOOLS_JAR" pre \
    "${JBAT}.txt" \
    "${JBAT}.hic" \
    "${JBAT}.chrom.sizes" \
    > "${JBAT}.make_hic.stdout" \
    2> "${JBAT}.make_hic.stderr"

echo "Finished creating .hic file: $(date)"

# final check
echo "Final JBAT files:"
ls -lh "${JBAT}.hic" "${JBAT}.assembly" "${JBAT}.chrom.sizes" "${JBAT}.txt"

echo "Testing whether .hic is readable:"
java -jar "$JUICER_TOOLS_JAR" dump observed NONE \
    "${JBAT}.hic" \
    assembly assembly BP 1000000 | head \
    > "${JBAT}.hic_read_test.txt"

echo "Read test output:"
cat "${JBAT}.hic_read_test.txt"

echo "Done: $(date)"
```
I added a final .hic file validation step because initial runs without this chunk ran without errors but produced unreadable .hic file (this was due to a Java heap space limit; I increased mem and switched to the bigmem node). 

This will generate .hic and .assembly files. scp these to a local directory and inspect them in Juicebox GUI.

The Hi-C contact map showed a strong diagonal signal and about nine distinct square blocks (likely macrochromosomes). Since we don't see any evidence of major misjoins (e.g., off-diagonal signals), we can proceed to the annnotation step.

Let's create the final assembly folder and copy the final FASTA, AGP, Hi-C files, and QC summaries:
```
# make the final assembly dir
mkdir -p /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly

# copy files
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly

cp ../scaffolding/yahs_out/Gloydius_ussuriensis_AMNH_21010_yahs_scaffolds_final.fa \
   Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa

cp ../scaffolding/yahs_out/Gloydius_ussuriensis_AMNH_21010_yahs_scaffolds_final.agp \
   Gloydius_ussuriensis_AMNH_21010_chromosome_level.agp

cp ../scaffolding/yahs_out/Gloydius_ussuriensis_AMNH_21010_yahs_JBAT.hic \
   Gloydius_ussuriensis_AMNH_21010_chromosome_level_JBAT.hic

cp ../scaffolding/yahs_out/Gloydius_ussuriensis_AMNH_21010_yahs_JBAT.assembly \
   Gloydius_ussuriensis_AMNH_21010_chromosome_level_JBAT.assembly
```

Then, index the final assembly FASTA
```
conda activate scaffolding
samtools faidx Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa
```

Create a scaffold size table
```
cut -f1,2 Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa.fai \
  | sort -k2,2nr \
  > Gloydius_ussuriensis_AMNH_21010_chromosome_level.scaffold_sizes.tsv
```

Generate checksums:
```
md5sum Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa \
       Gloydius_ussuriensis_AMNH_21010_chromosome_level.agp \
       Gloydius_ussuriensis_AMNH_21010_chromosome_level_JBAT.hic \
       Gloydius_ussuriensis_AMNH_21010_chromosome_level_JBAT.assembly \
       > Gloydius_ussuriensis_AMNH_21010_final_assembly.md5
```

### Assignment of scaffolds to chromosomes and manual assembly curation
Our Hi-C genome is currently at the scaffold level, i.e., they are not assigned to chromosomes yet. To do so, we will use the Eastern Diamondback (*Crotalus adamanteus*) reference genome assembly (assembled from a female [ZW]), the Prairie Rattlesnake (*Crotalus viridis*) reference genome (assembled from a male [ZZ]), and the Adder (*Vipera berus*) reference genome (assembled from a female [ZW]).  

Let's download these genomes. Below is for the *C.adamanteus* genome. Repeat for the *C.viridis* and *v. berus* genomes.
```sh
# cd into dir to store the assembly
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies

# create a conda env and install NCBI Datasets CLI
conda create -n ncbi_datasets -c conda-forge ncbi-datasets-cli unzip seqkit -y
conda activate ncbi_datasets

# download the genome
datasets download genome accession GCA_039797435.1 \
  --include genome,seq-report \
  --filename Crotalus_adamanteus_GCA_039797435.1.zip

# unzip the genome
unzip Crotalus_adamanteus_GCA_039797435.1.zip -d Crotalus_adamanteus_GCA_039797435.1

# symlink the genome
ln -s $(find Crotalus_adamanteus_GCA_039797435.1 -type f -name "*.fna" | head -n 1) \
      Crotalus_adamanteus_GCA_039797435.1.fa
```

Then index the genome and check chromosome names
```
# index
conda activate scaffolding
samtools faidx Crotalus_adamanteus_GCA_039797435.1.fa

# check names
cut -f1,2 Crotalus_adamanteus_GCA_039797435.1.fa.fai \
  | sort -k2,2nr \
  | head -n 40
```

Now proceed with whole-genome alignment using minimap2
```sh
#!/bin/bash
#SBATCH --job-name=assign_chromo
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate the conda env
source ~/.bash_profile
conda activate genome_assembly  # to access minimap2

set -euo pipefail

# set dir
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Crotalus_adamanteus"
mkdir -p "$workdir"
cd "$workdir"

# set variables
QUERY="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"
REF="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies/Crotalus_adamanteus_GCA_039797435.1.fa"

# sanity checks
echo "Checking files:"
ls -lh "$QUERY" "$REF"
which minimap2

# run minimap2
echo "Starting minimap2: $(date)"
minimap2 -x asm20 -t "${SLURM_CPUS_PER_TASK}" "$REF" "$QUERY" \
  > Gloydius_vs_Crotalus_adamanteus.asm20.paf

echo "Finished minimap2: $(date)"
ls -lh Gloydius_vs_Crotalus_adamanteus.asm20.paf
```

This will generate a .paf file. Summarize this into a scaffold-to-reference chromosome assignment table by running the Python script below.
```py
#!/usr/bin/env python3

import sys
from collections import defaultdict

paf = sys.argv[1]

# qname -> qlen
query_lengths = {}

# qname -> ref -> list of query intervals
intervals = defaultdict(lambda: defaultdict(list))

def merge_intervals(ints):
    if not ints:
        return []
    ints = sorted(ints)
    merged = [ints[0]]
    for start, end in ints[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged

def interval_total(ints):
    return sum(end - start for start, end in ints)

with open(paf) as f:
    for line in f:
        if not line.strip():
            continue

        fields = line.rstrip("\n").split("\t")

        qname = fields[0]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])

        rname = fields[5]
        aln_len = int(fields[10])
        mapq = int(fields[11])

        # filter small/noisy alignments
        if aln_len < 10000:
            continue
        if mapq < 5:
            continue

        query_lengths[qname] = qlen
        intervals[qname][rname].append((qstart, qend))

print(
    "query_scaffold\tquery_length\tbest_reference_chromosome\t"
    "best_unique_query_bp\ttotal_unique_query_bp\t"
    "best_fraction_of_aligned\tbest_fraction_of_scaffold"
)

for qname in sorted(query_lengths, key=lambda x: query_lengths[x], reverse=True):
    ref_intervals = intervals[qname]

    if not ref_intervals:
        continue

    ref_unique = {}
    all_intervals = []

    for rname, ints in ref_intervals.items():
        merged = merge_intervals(ints)
        bp = interval_total(merged)
        ref_unique[rname] = bp
        all_intervals.extend(ints)

    total_unique = interval_total(merge_intervals(all_intervals))
    best_ref, best_bp = max(ref_unique.items(), key=lambda x: x[1])

    best_fraction_of_aligned = best_bp / total_unique if total_unique else 0
    best_fraction_of_scaffold = best_bp / query_lengths[qname] if query_lengths[qname] else 0

    print(
        f"{qname}\t{query_lengths[qname]}\t{best_ref}\t"
        f"{best_bp}\t{total_unique}\t"
        f"{best_fraction_of_aligned:.4f}\t{best_fraction_of_scaffold:.4f}"
    )
```

Run it like this:
```
# go to the Python scripts dir
python chromo_assign_summary.py /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Crotalus_adamanteus/Gloydius_vs_Crotalus_adamanteus.asm20.paf > /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Crotalus_adamanteus/Gloydius_to_Crotalus_adamanteus_scaffold_assignment.tsv
```

Repeat this for *C. adamanteus*, *C. viridis*, and *V. berus* genomes. This will show that *C. adamanteus* and *C. viridis* chromosome 18 and *V. berus* chromosome 16 are "missing" from the *G. ussuriensis* assembly. To veify whether these chromosomes are truly missing from the assembly, or whether two chromosome were misjoined, let's inspect the raw. paf file and the Hi-C contact map. Although the first look at the contact map did not show any obvious misjoins, it's worth checking again. 

First, for the *C. adamanteus* .paf file, run:
```
awk '$6=="CM077934.1" {bp[$1]+=$11; n[$1]++}END {  for (s in bp) print s, bp[s], n[s]}' Gloydius_vs_Crotalus_adamanteus.asm20.paf | sort -k2,2nr | head -30
```
This command basically shows which *G. ussuriensis* scaffolds align to *C. adamanteus* chromosome 18, and how much sequence from each scaffold aligns there. In the .paf file, $1 is the *G. ussuriensis* scaffold, $6 is the *C. adamanteus* chromosome, and $11 is the aligned length. The output will show that *G. ussuriensis* scaffold 11 has ~11.06 Mb aligned to *C. adamanteus* chromosome 18.

Now, run below for the *C. viridis* .paf file.
```
awk '$6=="CM012322.1" {bp[$1]+=$11; n[$1]++}
END {
  for (s in bp) print s, bp[s], n[s]
}' Gloydius_vs_Crotalus_viridis.asm20.paf | sort -k2,2nr | head -30
```
The output will show that *G. ussuriensis* scaffold 11 has ~6.52 Mb aligned to *C. viridis* chromosome 18.

Lastly, run this for the *V. berus* .paf file:
```
awk '$6=="OZ077574.1" {
    bp[$1]+=$11;
    n[$1]++;
}
END {
    for (s in bp) print s, bp[s], n[s]
}' Gloydius_vs_Vipera_berus.asm20.paf | sort -k2,2nr | head -30
```
The output will show that *G. ussuriensis* scaffold 11 has ~11.37 Mb aligned to *V. berus* chromosome 16.

Since the *G. ussuriensis* scaffold 11 had highest matches to chromosome 11 in *C. adamanteus* and *C. viridis*, and chromosome 12 in *V. berus*, this suggests that two chromosomes were likely misjoined in the assembly. 

When you zoom into scaffold 11 in the Hi-C contact map, you will see the block corresponding to chromosome 11 actually contains two distinct chunks.

![alt_text](/etc/scaff_11.png)

Now, let's identify where to break this misjoin. I'm using *C. adamanteus* and *V. berus* here because the output from *C. viridis* paf file was noisier (but it didn't contract the results from the other two species).

Run below for *C. adamanteus*
```
awk '$1=="scaffold_11" && ($6=="CM077927.1" || $6=="CM077934.1") {
    ref=($6=="CM077927.1" ? "C_adamanteus_Chr11_like" : "C_adamanteus_Chr18_like");
    print $3, $4, ref, $6, $8, $9, $11, $5
}' Gloydius_vs_Crotalus_adamanteus.asm20.paf | sort -k1,1n | column -t
```            
This will identify where on scaffold 11 the transition happens from *C. adamanteus* "chr11 like" sequence to C. adamanteus "chr18 like" sequence. The output will show that chr11-like ends at position 17,166,186 and chr18-like starts at position 17,183,323.


And for *V. berus*
```
awk '$1=="scaffold_11" && ($6=="OZ077570.1" || $6=="OZ077574.1") {
    ref=($6=="OZ077570.1" ? "V_berus_Chr12_like" : "V_berus_Chr16_like");
    print $3, $4, ref, $6, $8, $9, $11, $5
}' Gloydius_vs_Vipera_berus.asm20.paf | sort -k1,1n | column -t
```
This will show that in scaffold 11, *V. berus* chr12-like portion ends at position 17,164,430 and chr16-like portion starts at position 17,183,534.

From the *C. adamanteus* result, we can see that there is a small gap between the transition point. We can take the midpoint of this gap ((17,166,186 + 17,183,323) / 2 = 17,174,754.5) to break this misjoin, rounded to 17,175,000.

Let's go to the final assmebly dir and create a directory for the curated assembly
```
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly

mkdir -p curated
cd curated
```
Then, run this:
```
# set variables
ASM="../Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"
CUT=17175000

OUT="Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
UNPLACED="Gloydius_ussuriensis_AMNH_21010_unplaced_lowconf_scaffolds.fa"
MAP="Gloydius_ussuriensis_AMNH_21010_chromosome_renaming.tsv"

# index & check assembly
conda activate scaffolding
samtools faidx "$ASM"
cut -f1,2 "$ASM.fai" | column -t | head -40
```
Run the below script from the "final_assembly" dir.
```sh
#!/bin/bash
set -euo pipefail

# Input final assembly
ASM="Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"

# Split point for scaffold_11
CUT=17175000

# Output directory
OUTDIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated"
mkdir -p "$OUTDIR"

# Output files
OUT="${OUTDIR}/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
UNPLACED="${OUTDIR}/Gloydius_ussuriensis_AMNH_21010_unplaced_lowconf_scaffolds.fa"
MAP="${OUTDIR}/Gloydius_ussuriensis_AMNH_21010_chromosome_rename.tsv"

# ------------------------------------------------------------
# Index input assembly
# ------------------------------------------------------------

samtools faidx "$ASM"

# Get scaffold_11 length
len11=$(awk '$1=="scaffold_11" {print $2}' "${ASM}.fai")

echo "scaffold_11 length: $len11"
echo "split coordinate:   $CUT"
echo "Chr11-like length:  $CUT"
echo "Chr18-like length:  $((len11 - CUT))"

# Remove previous outputs if present
rm -f "$OUT" "$OUT.fai" "$UNPLACED" "$UNPLACED.fai" "$MAP"

# ------------------------------------------------------------
# Write renaming / coordinate tracking table
# ------------------------------------------------------------

cat > "$MAP" <<EOF
old_scaffold	old_start	old_end	new_name	final_call
scaffold_1	1	end	G_ussuri_chr1	Chr1
scaffold_2	1	end	G_ussuri_chr2	Chr2
scaffold_3	1	end	G_ussuri_chr3	Chr3
scaffold_5	1	end	G_ussuri_chr4	Chr4
scaffold_6	1	end	G_ussuri_chr5	Chr5
scaffold_7	1	end	G_ussuri_chr6	Chr6
scaffold_8	1	end	G_ussuri_chr7	Chr7
scaffold_4	1	end	G_ussuri_chrZ	Z
scaffold_9	1	end	G_ussuri_chrW	W
scaffold_12	1	end	G_ussuri_chr9	Chr9
scaffold_10	1	end	G_ussuri_chr10	Chr10
scaffold_11	1	17175000	G_ussuri_chr11	Chr11
scaffold_13	1	end	G_ussuri_chr12	Chr12
scaffold_15	1	end	G_ussuri_chr13	Chr13
scaffold_14	1	end	G_ussuri_chr14	Chr14
scaffold_17	1	end	G_ussuri_chr15	Chr15
scaffold_18	1	end	G_ussuri_chr16	Chr16
scaffold_16	1	end	G_ussuri_chr17	Chr17
scaffold_11	17175001	end	G_ussuri_chr18	Chr18
EOF

# ------------------------------------------------------------
# Helper function: extract full scaffold and rename FASTA header
# ------------------------------------------------------------

extract_full () {
    old="$1"
    new="$2"

    if ! grep -q -w "$old" "${ASM}.fai"; then
        echo "ERROR: $old not found in ${ASM}.fai"
        exit 1
    fi

    samtools faidx "$ASM" "$old" | sed "1s/.*/>${new}/" >> "$OUT"
}

# ------------------------------------------------------------
# Build curated chromosome-only FASTA
# ------------------------------------------------------------

# Large autosomes
extract_full scaffold_1  G_ussuri_chr1
extract_full scaffold_2  G_ussuri_chr2
extract_full scaffold_3  G_ussuri_chr3
extract_full scaffold_5  G_ussuri_chr4
extract_full scaffold_6  G_ussuri_chr5
extract_full scaffold_7  G_ussuri_chr6
extract_full scaffold_8  G_ussuri_chr7

# Sex chromosomes
extract_full scaffold_4  G_ussuri_chrZ
extract_full scaffold_9  G_ussuri_chrW

# Smaller autosomes
extract_full scaffold_12 G_ussuri_chr9
extract_full scaffold_10 G_ussuri_chr10

# Split scaffold_11: first piece = Chr11
samtools faidx "$ASM" "scaffold_11:1-${CUT}" | \
    sed "1s/.*/>G_ussuri_chr11/" >> "$OUT"

extract_full scaffold_13 G_ussuri_chr12
extract_full scaffold_15 G_ussuri_chr13
extract_full scaffold_14 G_ussuri_chr14
extract_full scaffold_17 G_ussuri_chr15
extract_full scaffold_18 G_ussuri_chr16
extract_full scaffold_16 G_ussuri_chr17

# Split scaffold_11: second piece = Chr18
samtools faidx "$ASM" "scaffold_11:$((CUT + 1))-${len11}" | \
    sed "1s/.*/>G_ussuri_chr18/" >> "$OUT"

# Index curated FASTA
samtools faidx "$OUT"

# ------------------------------------------------------------
# Save low-confidence / unplaced scaffolds separately
# ------------------------------------------------------------

HIGHCONF="${OUTDIR}/highconf_scaffolds.txt"
ALL="${OUTDIR}/all_scaffolds.txt"
LOWCONF="${OUTDIR}/lowconf_scaffolds.txt"

cat > "$HIGHCONF" <<EOF
scaffold_1
scaffold_2
scaffold_3
scaffold_4
scaffold_5
scaffold_6
scaffold_7
scaffold_8
scaffold_9
scaffold_10
scaffold_11
scaffold_12
scaffold_13
scaffold_14
scaffold_15
scaffold_16
scaffold_17
scaffold_18
EOF

cut -f1 "${ASM}.fai" > "$ALL"

grep -vxFf "$HIGHCONF" "$ALL" > "$LOWCONF"

while read -r s; do
    samtools faidx "$ASM" "$s" >> "$UNPLACED"
done < "$LOWCONF"

samtools faidx "$UNPLACED"

# ------------------------------------------------------------
# Final checks
# ------------------------------------------------------------

echo
echo "Curated chromosome FASTA:"
echo "$OUT"
echo

echo "Number of curated chromosome sequences:"
grep -c "^>" "$OUT"

echo
echo "Curated chromosome sizes:"
cut -f1,2 "$OUT.fai" | column -t

echo
echo "Split scaffold_11 pieces:"
grep -E "G_ussuri_chr11|G_ussuri_chr18" "$OUT.fai" | column -t

echo
echo "Unplaced / low-confidence scaffolds:"
echo "$UNPLACED"
echo "Number of unplaced sequences:"
grep -c "^>" "$UNPLACED"

echo
echo "Renaming table:"
echo "$MAP"

echo
echo "Done."
```

After this, run a quick check on the curated assembly:
```
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated

OUT="Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"
UNPLACED="Gloydius_ussuriensis_AMNH_21010_unplaced_lowconf_scaffolds.fa"

echo "Number of curated chromosome sequences:"
grep -c "^>" "$OUT"

echo "Curated chromosome sizes:"
cut -f1,2 "$OUT.fai" | column -t

echo "Split scaffold_11 pieces:"
grep -E "G_ussuri_chr11|G_ussuri_chr18" "$OUT.fai" | column -t

echo "Number of unplaced/low-confidence scaffolds:"
grep -c "^>" "$UNPLACED" 
```

Then make a checksum
```
md5sum Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa \
  > Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa.md5

md5sum Gloydius_ussuriensis_AMNH_21010_unplaced_lowconf_scaffolds.fa \
  > Gloydius_ussuriensis_AMNH_21010_unplaced_lowconf_scaffolds.fa.md5
```

Also check the size of the curated assembly
```
awk '{sum+=$2} END{print "curated_bp", sum}' \
Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa.fai
```
After this, run a final synteny check with the *C. adamanteus* genome and re-build the Hi-C contact map.

When you re-run minimap2 with the curated assembly and *C. adamanteus* genome, you will see the output shows excellent chromosome assignment summary, with 1:1 matches between each *G. ussuriensis* scaffold and named *C. adamanteus* chromosome. 

To rebuild the Hi-C contact map, run the script below:
```sh
#!/bin/bash
#SBATCH --job-name=curated_hicmap
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=450G
#SBATCH --time=72:00:00
#SBATCH --partition=bigmem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# init conda env
source ~/.bash_profile
conda activate scaffolding

set -euo pipefail

THREADS="${SLURM_CPUS_PER_TASK}"

# working directory for rebuilt Hi-C map
workdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/hic_rebuild"
mkdir -p "$workdir"
cd "$workdir"

# curated split assembly
FASTA="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"

# path to hic reads
hicdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scaffolding/combined"
R1="${hicdir}/Gloydius_ussuriensis_HiC_R1.fastq.gz"
R2="${hicdir}/Gloydius_ussuriensis_HiC_R2.fastq.gz"

# output prefix
prefix="Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split"

# juicer .jar path
JUICER_TOOLS_JAR="/home/yshin/mendel-nas1/juicer_tools/juicer_tools_1.22.01.jar"

# sanity checks
echo "Checking input files: $(date)"
ls -lh "$FASTA"
ls -lh "$R1"
ls -lh "$R2"
ls -lh "$JUICER_TOOLS_JAR"

which bwa-mem2
which samtools
which yahs
which juicer
which java

java -version

echo "Indexing curated FASTA: $(date)"
samtools faidx "$FASTA"

if [ ! -e "${FASTA}.0123" ] && [ ! -e "${FASTA}.bwt.2bit.64" ]; then
    bwa-mem2 index "$FASTA"
fi

# map hic reads to curated split assembly
echo "Mapping Hi-C reads to curated split FASTA: $(date)"
bwa-mem2 mem -5SP -t "$THREADS" "$FASTA" "$R1" "$R2" | \
    samtools view -@ "$THREADS" -b -F 3340 - | \
    samtools sort -@ "$THREADS" -o "${prefix}.sorted.bam" -

samtools index "${prefix}.sorted.bam"

# run yahs
echo "Running YaHS on curated split assembly: $(date)"
yahs "$FASTA" "${prefix}.sorted.bam" -o "$prefix"

# run juicer pre to generate .hic file
echo "Running juicer pre: $(date)"
juicer pre \
    -a \
    -o "${prefix}_JBAT" \
    "${prefix}.bin" \
    "${prefix}_scaffolds_final.agp" \
    "${FASTA}.fai" \
    > "${prefix}_JBAT.log" 2>&1

echo "Making chromosome sizes file: $(date)"
grep PRE_C_SIZE "${prefix}_JBAT.log" | awk '{print $2, $3}' > "${prefix}_JBAT.chrom.sizes"

echo "Chromosome sizes:"
cat "${prefix}_JBAT.chrom.sizes"

echo "Checking JBAT contact file:"
ls -lh "${prefix}_JBAT.txt"
head -n 5 "${prefix}_JBAT.txt"

echo "Creating .hic file: $(date)"
java -Xmx400G -jar "$JUICER_TOOLS_JAR" pre \
    "${prefix}_JBAT.txt" \
    "${prefix}_JBAT.hic" \
    "${prefix}_JBAT.chrom.sizes" \
    > "${prefix}_JBAT.make_hic.stdout" \
    2> "${prefix}_JBAT.make_hic.stderr"

echo "Final files:"
ls -lh "${prefix}_JBAT.hic" "${prefix}_JBAT.assembly" "${prefix}_JBAT.chrom.sizes" "${prefix}_JBAT.txt"

# test .hic readability
echo "Testing .hic readability:"
java -jar "$JUICER_TOOLS_JAR" dump observed NONE \
    "${prefix}_JBAT.hic" \
    assembly assembly BP 1000000 | head \
    > "${prefix}_JBAT.hic_read_test.txt" || true

cat "${prefix}_JBAT.hic_read_test.txt" || true

echo "Done: $(date)"
```

Before running the script, also install BWA-MEM2 in the scaffolding conda env:
```
conda install bioconda::bwa-mem2
```

While the above script is running, run __one final QC__ on the curated assembly by re-running QUAST, compleasm, and Merqury. 

Create the final QC dir and re-run:
```
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated
mkdir final_QC final_QC/QUAST final_QC/compleasm final_QC/Merqury
```

### Sex chromosome validation based on sex-specific read coverage patterns
Our synteny-based chromosome assignment, using *C. adamanteus* (ZW), *C. viridis* (ZZ), and *V. berus* (ZW) assemblies, identified Z and W chromosomes in the *G. ussuriensis* assembly. As an independent validation of this assignment, we will use low-coverage whole-genome resequencing data from 4 males and 8 females and map sequencing reads from these samples to the reference *G. ussuriensis* assembly. 


Create a directory for this analysis.
```
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly
mkdir -p sex_chr_coverage
cd sex_chr_coverage
```

And create a low coverage WGS sample metadata file. 
```
# create the metadata file
cat > sample_metadata.tsv <<'EOF'
sample	sex	locality
AMNH_21010	F	Hwacheon
AMNH_21050	M	Jeju
AMNH_21060	M	Jeju
AMNH_21070	F	Udo
AMNH_21073	F	Udo
AMNH_21128	F	Gumi
AMNH_21147	M	Daejeon
AMNH_21161	F	Sancheong
AMNH_21162	F	Sancheong
AMNH_21164	F	Yeoncheon
AMNH_21172	F	Busan
AMNH_21185	M	Yangyang
EOF

# use printf so bash inserts real tab characters
printf "sample\tsex\tlocality\n" > sample_metadata.tsv
printf "AMNH_21010\tF\tHwacheon\n" >> sample_metadata.tsv
printf "AMNH_21050\tM\tJeju\n" >> sample_metadata.tsv
printf "AMNH_21060\tM\tJeju\n" >> sample_metadata.tsv
printf "AMNH_21070\tF\tUdo\n" >> sample_metadata.tsv
printf "AMNH_21073\tF\tUdo\n" >> sample_metadata.tsv
printf "AMNH_21128\tF\tGumi\n" >> sample_metadata.tsv
printf "AMNH_21147\tM\tDaejeon\n" >> sample_metadata.tsv
printf "AMNH_21161\tF\tSancheong\n" >> sample_metadata.tsv
printf "AMNH_21162\tF\tSancheong\n" >> sample_metadata.tsv
printf "AMNH_21164\tF\tYeoncheon\n" >> sample_metadata.tsv
printf "AMNH_21172\tF\tBusan\n" >> sample_metadata.tsv
printf "AMNH_21185\tM\tYangyang\n" >> sample_metadata.tsv
```

Then, prep the *G. ussuriensis* reference assembly.
```sh
# in the "sex_chr_coverage" dir
# set ref genome as a var
ref="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"

# take the first two columns from the fasta index file and save it to a new file called "ref.genome"
cut -f1,2 "${ref}.fai" > ref.genome
```

check the chromosome names in the assembly fasta.
```sh
awk '{print NR, $1, $2}' "${ref}.fai" | sort -k3,3nr | head -30
```

create a chromosome name classification file.
```
# make chromosome class file
cat > chr_sets.tsv <<'EOF'
chrom	class
G_ussuri_chr1	AUTO
G_ussuri_chr2	AUTO
G_ussuri_chr3	AUTO
G_ussuri_chr4	AUTO
G_ussuri_chr5	AUTO
G_ussuri_chr6	AUTO
G_ussuri_chr7	AUTO
G_ussuri_chrZ	Z
G_ussuri_chrW	W
G_ussuri_chr9	AUTO
G_ussuri_chr10	AUTO
G_ussuri_chr11	AUTO
G_ussuri_chr12	AUTO
G_ussuri_chr13	AUTO
G_ussuri_chr14	AUTO
G_ussuri_chr15	AUTO
G_ussuri_chr16	AUTO
G_ussuri_chr17	AUTO
G_ussuri_chr18	AUTO
EOF

# use printf so bash inserts real tab characters
printf "chrom\tclass\n" > chr_sets.tsv
printf "G_ussuri_chr1\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr2\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr3\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr4\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr5\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr6\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr7\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chrZ\tZ\n" >> chr_sets.tsv
printf "G_ussuri_chrW\tW\n" >> chr_sets.tsv
printf "G_ussuri_chr9\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr10\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr11\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr12\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr13\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr14\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr15\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr16\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr17\tAUTO\n" >> chr_sets.tsv
printf "G_ussuri_chr18\tAUTO\n" >> chr_sets.tsv
```

Next, make 100 kb windows using bedtools. First, install bedtools and mosdepth into the scaffolding conda env:
```
conda install bioconda::bedtools
conda install bioconda::mosdepth
```
Then run:
```
bedtools makewindows -g ref.genome -w 100000 > windows_100kb.bed
```

After this, map reads and calculate window coverage:
```sh
#!/bin/bash
#SBATCH --job-name=sex_chr_depth
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --partition=compute
#SBATCH --array=1-12
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env
source ~/.bash_profile
conda activate scaffolding

set -euo pipefail

# project dir
project="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/sex_chr_coverage"

# chromosome level assembly
ref="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/curated/Gloydius_ussuriensis_AMNH_21010_curated_scaffold11_split.fa"

# low coverage whole genome sequencing reads
reads="/home/yshin/mendel-nas1/ussuri_popgen/Genome_seq/04.Trimmed"

# metadata and window files
metadata="${project}/sample_metadata.tsv"
windows="${project}/windows_100kb.bed"

# slurm threads
threads="${SLURM_CPUS_PER_TASK}"

# chooses which sample an array task should process
sample=$(awk -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n+1 {print $1}' "$metadata")

r1="${reads}/${sample}_R1.trimmed.fq.gz"
r2="${reads}/${sample}_R2.trimmed.fq.gz"

echo "Sample: $sample"
echo "R1: $r1"
echo "R2: $r2"

# spit out errors if read files are missing
if [[ ! -s "$r1" ]]; then
    echo "ERROR: missing R1: $r1"
    exit 1
fi

if [[ ! -s "$r2" ]]; then
    echo "ERROR: missing R2: $r2"
    exit 1
fi

# create output & other directories
mkdir -p "${project}/bam" "${project}/depth" "${project}/tmp"

tmpdir="${project}/tmp/${sample}"
mkdir -p "$tmpdir"

if command -v bwa-mem2 >/dev/null 2>&1; then
    aligner="bwa-mem2 mem"
else
    aligner="bwa mem"
fi

echo "Using aligner: $aligner"

# Name-sort for duplicate marking
$aligner -t "$threads" \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA" \
    "$ref" "$r1" "$r2" \
    | samtools sort -n -@ "$threads" -T "${tmpdir}/${sample}.name" \
      -o "${tmpdir}/${sample}.name.bam" -

# Fix mate information
samtools fixmate -m \
    "${tmpdir}/${sample}.name.bam" \
    "${tmpdir}/${sample}.fixmate.bam"

# Coordinate sort
samtools sort -@ "$threads" -T "${tmpdir}/${sample}.coord" \
    -o "${tmpdir}/${sample}.coord.bam" \
    "${tmpdir}/${sample}.fixmate.bam"

# Remove duplicates
samtools markdup -r -@ "$threads" \
    "${tmpdir}/${sample}.coord.bam" \
    "${project}/bam/${sample}.markdup.bam"

samtools index -@ "$threads" "${project}/bam/${sample}.markdup.bam"

# Windowed depth.
# MAPQ 30 keeps mainly confidently mapped reads.
mosdepth -t "$threads" \
    --mapq 30 \
    --by "$windows" \
    "${project}/depth/${sample}" \
    "${project}/bam/${sample}.markdup.bam"

# Basic mapping statistics
samtools flagstat -@ "$threads" \
    "${project}/bam/${sample}.markdup.bam" \
    > "${project}/bam/${sample}.flagstat.txt"

rm -rf "$tmpdir"

echo "Finished sample: $sample"
```

After this finishes running, summarize and plot normalized Z/W/autosome coverage in R.
```r
########## summarize normalized Z/W/autosome coverage in R

# clean working env
rm(list = ls(all.names = T))
gc()

# load packages
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# local R project directory
project <- '/home/yshin/Gloydius_ussuriensis_genome_assembly/R'

# local input and output folders
data_dir <- file.path(project, 'Rdata', 'sex_chr_coverage')
plot_dir <- file.path(project, 'Rplots')

# load metadata and chromosome set files
meta <- fread(file.path(data_dir, 'sample_metadata.tsv'))
chr_sets <- fread(file.path(data_dir, 'chr_sets.tsv'))

# get list of mosdepth depth files
depth_files <- list.files(file.path(data_dir, 'depth'), pattern = '\\.regions\\.bed\\.gz$', full.names = T)

if (length(depth_files) == 0) {
  stop('No mosdepth .regions.bed.gz files found')
}

cat('Number of depth files found:', length(depth_files), '\n')

# read depth files
depth_list <- lapply(depth_files, function(f) {
  x <- fread(f, col.names = c('chrom', 'start', 'end', 'depth'))
  x$sample <- sub('\\.regions\\.bed\\.gz$', '', basename(f))
  x
})

depth <- rbindlist(depth_list)
print(depth)

# keep only chromosomes listed in chr_sets.tsv
depth <- depth %>%
  inner_join(chr_sets, by = 'chrom') %>%
  left_join(meta, by = 'sample')

head(depth)

# check for missing metadata
missing_meta <- depth %>%
  filter(is.na(sex) | is.na(locality)) %>%
  distinct(sample)

if (nrow(missing_meta) > 0) {
  print(missing_meta)
  stop('Some depth files do not have matching sample metadata')
}

# autosome median depth per sample
auto_median <- depth %>%
  filter(class == 'AUTO') %>%
  group_by(sample) %>%
  summarise(auto_median_depth = median(depth, na.rm = T), .groups = 'drop')

head(auto_median)

# normalize depth by autosomal median
depth_norm <- depth %>%
  left_join(auto_median, by = 'sample') %>%
  mutate(norm_depth = depth / auto_median_depth)

head(depth_norm)

# chromosome-level summary
chr_summary <- depth_norm %>%
  group_by(sample, sex, locality, chrom, class) %>%
  summarise(median_norm_depth = median(norm_depth, na.rm = T),
            mean_norm_depth = mean(norm_depth, na.rm = T),
            n_windows = n(),
            .groups = 'drop')

head(chr_summary)

# class-level summary: AUTO, Z, W
class_summary <- depth_norm %>%
  group_by(sample, sex, locality, class) %>%
  summarise(median_norm_depth = median(norm_depth, na.rm = T), 
            mean_norm_depth = mean(norm_depth, na.rm = T),
            n_windows = n(),
            .groups = 'drop')

class_wide <- class_summary %>%
  select(sample, sex, locality, class, median_norm_depth) %>%
  pivot_wider(names_from = class, values_from = median_norm_depth)

# write output summary tables to Rdata
write.table(chr_summary, file = file.path(data_dir, 'sexchr_chromosome_depth_summary.tsv'), sep = '\t', quote = F, row.names = F)
write.table(class_summary, file = file.path(data_dir, 'sexchr_class_depth_summary.tsv'), sep = '\t', quote = F, row.names = F)
write.table(class_wide,file = file.path(data_dir, 'sexchr_class_depth_wide.tsv'), sep = '\t', quote = F, row.names = F)

# plot 1: sex-specific Z/W coverage validation
class_wide$sex <- recode(class_wide$sex, 'F' = 'Female', 'M' = 'Male')

p1 <- ggplot(class_wide, aes(x = Z, y = W, shape = sex, fill = sex)) +
  geom_point(size = 4.5, color = 'black', stroke = 1.5) +
  geom_vline(xintercept = 0.75, linetype = 'dashed') +
  geom_hline(yintercept = 0.15, linetype = 'dashed') +
  scale_fill_manual(values = c(Female = 'salmon', Male = 'cornflowerblue')) +
  scale_shape_manual(values = c(Female = 21, Male = 24)) +
  labs(x = 'Normalized Z depth', y = 'Normalized W depth') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 14))

print(p1)

# plot 2: autosome-normalized coverage by chromosome class per sample
class_summary$sample <- factor(class_summary$sample,
                               levels = c('AMNH_21010', 'AMNH_21070', 'AMNH_21073', 'AMNH_21128',
                                          'AMNH_21161', 'AMNH_21162', 'AMNH_21164', 'AMNH_21172',
                                          'AMNH_21050', 'AMNH_21060', 'AMNH_21147', 'AMNH_21185'))

class_summary$class <- recode(class_summary$class, 'AUTO' = 'Autosome', 'Z' = 'Z', 'W' = 'W')
class_summary$class <- factor(class_summary$class, levels = c('Autosome', 'Z', 'W'))


p2 <- ggplot(class_summary, aes(x = sample, y = median_norm_depth, fill = class)) +
  geom_col(position = 'dodge') +
  coord_flip() +
  scale_fill_manual(values = c(Autosome = 'grey80',
                               Z = 'steelblue3',
                               W = 'goldenrod2')) +
  labs(x = 'Sample', y = 'Median normalized depth') +
  theme_bw() +
  theme(legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 20)),
        axis.text = element_text(size = 14))

print(p2)

# plot 3: known sex comparison by chromosome class
#p3 <- class_summary %>%
#  filter(class %in% c('AUTO', 'Z', 'W')) %>%
#  ggplot(aes(x = sex, y = median_norm_depth)) +
#  geom_boxplot(outlier.shape = NA) +
#  geom_jitter(width = 0.15, height = 0, size = 2) +
#  facet_wrap(~class, scales = 'free_y') +
#  labs(x = 'Known sex', y = 'Median normalized depth') +
#  theme_bw()

#print(p3)

# save plots
ggsave(filename = file.path(plot_dir, 'Z_vs_W_normalized_depth.png'), plot = p1, width = 7, height = 6, dpi = 800)
ggsave(filename = file.path(plot_dir, 'AUTO_Z_W_depth_by_sample.png'), plot = p2, width = 8, height = 7, dpi = 800)
#ggsave(filename = file.path(plot_dir, 'known_sex_AUTO_Z_W_validation.png'), plot = p3, width = 7, height = 5, dpi = 300)
```
