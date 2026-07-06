#!/bin/bash
#SBATCH --job-name=calcCov_scaffold
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --cpus-per-task=32
#SBATCH --partition=compute
#SBATCH --time=20:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env for mapping
source ~/.bash_profile
conda activate genome_assembly

set -euo pipefail

# set paths
ref_fa="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/Gloydius_ussuriensis_AMNH_21010_chromosome_level.fa"
reads_fq="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/genome_cleanup/AMNH_21010_HiFi.fastq.gz"
out_dir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/final_assembly/scaffold_level_cov"

mkdir -p "${out_dir}"
mkdir -p "${out_dir}/tmp"

prefix=scaffold_level

# index fasta
samtools faidx "${ref_fa}"

# map HiFi reads
echo "start mapping....."

minimap2 -t "${SLURM_CPUS_PER_TASK}" -ax map-hifi "${ref_fa}" "${reads_fq}" \
| samtools sort -@ 8 -T "${out_dir}/tmp/${prefix}" -o "${out_dir}/${prefix}.bam"

samtools index "${out_dir}/${prefix}.bam"

# activate a separate conda env for newer samtools version
source ~/.bash_profile
conda activate samtools

echo "activated separate conda env for newer samtools version....."

# calculate mean genome-wide coverage
echo "calculate mean genome-wide coverage....."

samtools coverage "${out_dir}/${prefix}.bam" \
| awk 'NR>1 {
    len_i = $3 - $2 + 1;
    sum += len_i * $7;
    len += len_i;
  }
  END {
    print "mean depth of coverage =", sum/len "x";
  }' \
| tee "${out_dir}/${prefix}.mean_depth.txt"

# calculate how much of the genome (%) is covered >= 20x, 30x, and 50x
echo "calculate how much of the genome (%) is covered >= 20x, 30x, and 50x....."

samtools depth -a "${out_dir}/${prefix}.bam" \
| awk '{
    total++;
    if ($3 >= 20) c20++;
    if ($3 >= 30) c30++;
    if ($3 >= 50) c50++;
  }
  END {
    print ">=20x =", 100*c20/total "%";
    print ">=30x =", 100*c30/total "%";
    print ">=50x =", 100*c50/total "%";
  }' \
| tee