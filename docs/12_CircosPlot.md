## 12) Circos plot
### Step 1: Setup
Set up the directories:
```sh
CIRCOS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/circos"

mkdir -p \
    "${CIRCOS}/01_windows" \
    "${CIRCOS}/02_tracks" \
    "${CIRCOS}/03_plots"
```

Run the following Python script to prepare tracks for generating a circos plot. This script will generate chromosome lengths, GC %, repeat %, and gene density at 1 Mb window.
```sh
# activate conda env to access Python
conda activate funannotate

# run the script
python \
    /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/prep_circos_tracks.py
```

After this, scp the three output files to a local directory and run the R script below to generate a circos plot:
```r

```