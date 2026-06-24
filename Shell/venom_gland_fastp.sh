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