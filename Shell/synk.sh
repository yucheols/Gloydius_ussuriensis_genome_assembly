#!/bin/sh
#SBATCH --job-name=synk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=20:00:00
#SBATCH --partition=compute
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err

# activate conda env
source /home/yshin/mendel-nas1/miniconda3/etc/profile.d/conda.sh
conda activate synk

# synk.py path
SYNK="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Synk/synk.py"

# G. ussuriensis assembly
ussuri_asm="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.softmasked.fasta"

# comaprison genomes base dir
comp_dir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Synk/chromosomes_fa"

# outdir
outdir="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Synk/outputs"
mkdir -p ${outdir}

# run synk
echo "Running synk: $(date)"

python ${SYNK} \
  --main_name G_ussuri \
  --main_assembly ${ussuri_asm} \
  --compare A_dia=${comp_dir}/Argyrophis_diardii.chromosomes.fa  \
  --compare B_ins=${comp_dir}/Bothrops_insularis.chromosomes.fa \
  --compare C_asp=${comp_dir}/Candoia_aspera.chromosomes.fa  \
  --compare C_gas=${comp_dir}/Cerastes_gasperettii.chromosomes.fa \
  --compare C_ada=${comp_dir}/Crotalus_adamanteus.chromosomes.fa \
  --compare C_vir=${comp_dir}/Crotalus_viridis.chromosomes.fa  \
  --compare E_sch=${comp_dir}/Elaphe_schrenckii.chromosomes.fa \
  --compare G_she=${comp_dir}/Gloydius_shedaoensis.chromosomes.fa \
  --compare N_naj=${comp_dir}/Naja_naja.chromosomes.fa \
  --compare V_ber=${comp_dir}/Vipera_berus.chromosomes.fa \
  --compare X_uni=${comp_dir}/Xenopeltis_unicolor.chromosomes.fa \
  --lineage sauropsida \
  --threads ${SLURM_CPUS_PER_TASK:-8} \
  --reuse_compleasm \
  --outdir ${outdir} \
  --plot

echo "Synk run completed: $(date)"