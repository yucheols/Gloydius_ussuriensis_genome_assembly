#!/bin/bash
#SBATCH --job-name=genomeCleanup
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

# set some variables
indir=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup

#####  step 1 ::: map hifi reads to assembly
echo "Step 1 == Mapping hifi reads to assembly..."

# map
minimap2 -t ${SLURM_CPUS_PER_TASK} -ax map-hifi ${indir}/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa \
  ${indir}/AMNH_21010_HiFi.fastq.gz | samtools sort -@ 8 -o hifi_to_asm.bam

samtools index hifi_to_asm.bam

# mean depth per contig
samtools coverage hifi_to_asm.bam > contig_depth.tsv

echo "Step 1 == completed"


#####  step 2 ::: identify mitochondrial contigs
echo "Step 2 == Check for mitochondrial hits..."

# run minimap using the G. ussuriensis mitochondrial reference
minimap2 -x asm5 ${indir}/NC_026553.1.fa ${indir}/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa > mt_hits.paf
cut -f6 mt_hits.paf | sort | uniq > mt_contigs.txt

echo "Step 2 == completeed"


#####  step 3 ::: screen for bacterial & other contaminants
echo "Step 3 == Check for contaminants..."

# run diamond
diamond blastx -q ${indir}/Gloydius_ussuriensis_v1.asm.bp.p_ctg.fa \
  -d nr -o diamond_result.tsv -f 6 qseqid staxids bitscore evalue -k 1 -p ${SLURM_CPUS_PER_TASK}

# output contaminants
awk '
($5 <= 1e-20 && $4 >= 200) &&
($3 ~ /Bacteria|Proteobacteria|Firmicutes|Actinobacteria|Fungi|Ascomycota|Basidiomycota/)
{ print $1 }
' diamond_tax.tsv | sort -u > diamond_hits.txt


