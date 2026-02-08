#!/bin/bash
#SBATCH --job-name=blobBlst_ussuri
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=32
#SBATCH --array=1-3
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/nas4/G_ussuriensis_Chromo/slurm_logs/slurm-%x_%j.out
#SBATCH --error=/home/yshin/nas4/G_ussuriensis_Chromo/slurm_logs/slurm-%x_%j.err

# this is run on the AMNH Huxley cluster
# load module
module load BLAST/2.11.0-Linux_x86_64 

# path to genome cleanup directory
wdir=/home/yshin/nas4/G_ussuriensis_Chromo/genome_cleanup

# path to DB
export BLASTDB=/nas3/database/nt_2024_10_23
blastdbcmd -db nt -info

# path to assembly
path_to_asm=${wdir}/asm_div

# output dir
out_dir=/home/yshin/nas4/G_ussuriensis_Chromo/genome_cleanup/blast_out

# assign and print array ID
i=${SLURM_ARRAY_TASK_ID}
P=$(awk "NR==$i" ${path_to_asm}/blast.list)

# run blast
blastn \
  -db nt \
  -task megablast \
  -query ${path_to_asm}/${P} \
  -out ${out_dir}/${P}.blast.out \
  -max_target_seqs 10 \
  -max_hsps 1 \
  -evalue 1e-20 \
  -outfmt "6 qseqid staxids bitscore sseqid sskingdoms sscinames" \
  -num_threads=${SLURM_CPUS_PER_TASK}  