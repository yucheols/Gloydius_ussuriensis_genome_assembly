## 8) Chromosomal synteny
### Step 1: Setup directory and prepare data downloads
First, set up a directory for synteny analyses.
```sh
# under the "G_ussuriensis_Chromo" directory
mkdir -p synteny/{raw_ncbi,assemblies,gff3,proteins,metadata,logs,scripts}
```

Activate the conda env to access NCBI Datasets CLI
```sh
conda activate ncbi_datasets
```

### Step 2: Download assemblies
There are several chromosom-level snake reference genome assemblies, especially viperid assemblies, available. We can use these to investigate the synteny across species. Along with *C. adamanteus*, *C. viridis*, and *V. berus* assemblies, we will download eight additional assemblies. These additional species are: *Cerastes gasperettii*, *Bothrops insularis*, *Gloydius shedaoensis*, *Naja naja*, *Elaphe schrenckii*, *Candoia aspera*, *Xenopeltis unicolor*, *Argyrophis diardii*.

First, create a .csv file containing species name, family/subfamily names, assembly name, assembly accession, data source, and url. Note that all assemblies other than *G. shedaoensis* are deposited GenBank. The *G. shedaoensis* assembly is available from National Genomics Data Center (NGDC) Genome Warehouse (GWH).

Store the .csv file in the "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/metadata" with the following name:
```sh
synteny_samples_list.csv
```

For the NCBI data download, we can use the "ncbi_datasets" conda env. But, to download data from NGDC, we need to use curl and jq. Download these two commands into the "ncbi_datasets" conda env.
```sh
conda install -n ncbi_datasets -c conda-forge \
    curl jq -y
```

Next, submit the following shell script to Mendel. This script will take the .csv file of assmbly manifest and download the assemblies that will be used in synteny analyses. The script applies several conditions based on where genome assembly, gff, and protein fasta files are available. First, it downloads the genome assembly fasta based on the "genome_url" column. It then downloads the gff file based on RefSeq accession if the annotation is available there. If not, the script will attempt to download the gff file from an external url provided in the "annotation_url" column. Same applies to the protein fasta.
```sh
download_assemblies.sh
```

For *C. adamanteus* let's separately download .gff and protein fasta. We are doing this because access information for these two files is not very visible on NCBI.

First, download the protein fasta:
```sh
cd /home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny/Crotalus_adamanteus

wget \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/797/435/GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2/GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_protein.faa.gz
```

Check download:
```sh
ls -lh GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_protein.faa.gz
gzip -t GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_protein.faa.gz
zgrep -c '^>' GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_protein.faa.gz
```

Then give it a standardiezed file name:
```sh
gunzip -c \
    GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_protein.faa.gz \
    > Crotalus_adamanteus.protein.faa
```

Do the same for .gff:
```sh
# download
wget \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/039/797/435/GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2/GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_genomic.gff.gz

# standardize file name 
gunzip -c \
    GCA_039797435.1_Cadamanteus_3dDNAHiC_1.2_genomic.gff.gz \
    > Crotalus_adamanteus.annotation.gff3
```

Then check the number of genes and CDS contained in this file:
```sh
head Crotalus_adamanteus.annotation.gff3

awk -F '\t' '$3=="gene"{n++} END{print "genes:",n}' \
    Crotalus_adamanteus.annotation.gff3

awk -F '\t' '$3=="CDS"{n++} END{print "CDS:",n}' \
    Crotalus_adamanteus.annotation.gff3
```

Also, download gff and protein fasta separately for V. berus and C. aspera, because the assembly download script somehow omitted these files.

For V. berus:
```sh
# from the assemblies_synteny/Vipera_berus dir
# download gff and extract
wget -qO- \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/194/415/GCF_964194415.1_rVipBer3.hap1.1/GCF_964194415.1_rVipBer3.hap1.1_genomic.gff.gz \
    | gunzip -c \
    > Vipera_berus.annotation.gff3  

# download protein fasta
wget -qO- \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/964/194/415/GCF_964194415.1_rVipBer3.hap1.1/GCF_964194415.1_rVipBer3.hap1.1_protein.faa.gz \
    | gunzip -c \
    > Vipera_berus.protein.faa
```
Then verify downloads:
```sh
ls -lh \
    Vipera_berus.annotation.gff3 \
    Vipera_berus.protein.faa

awk -F '\t' '$3=="gene"{n++} END{print "genes:",n+0}' \
    Vipera_berus.annotation.gff3

grep -c '^>' \
    Vipera_berus.protein.faa
```

Do the same for C. aspera:
```sh
# from the assemblies_synteny/Candoia_aspera dir
# download gff and extract
wget -qO- \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/149/785/GCF_035149785.1_rCanAsp1.hap2/GCF_035149785.1_rCanAsp1.hap2_genomic.gff.gz \
    | gunzip -c \
    > Candoia_aspera.annotation.gff3

# download protein fasta
wget -qO- \
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/035/149/785/GCF_035149785.1_rCanAsp1.hap2/GCF_035149785.1_rCanAsp1.hap2_protein.faa.gz \
    | gunzip -c \
    > Candoia_aspera.protein.faa    
```

Again, verify downloads:
```sh
awk -F '\t' '$3=="gene"{n++} END{print "genes:",n+0}' \
    Candoia_aspera.annotation.gff3

awk -F '\t' '$3=="CDS"{n++} END{print "CDS:",n+0}' \
    Candoia_aspera.annotation.gff3

grep -c '^>' Candoia_aspera.protein.faa
```

### Step 3: Install software
We will create several conda env for handling various steps in the synteny analyses.

##### General genome data handling & qc env
```sh
conda create -n synteny_qc \
    -c conda-forge -c bioconda \
    seqkit \
    samtools \
    bedtools \
    gffread \
    agat \
    compleasm \
    minimap2 \
    mummer4
```

##### GENESPACE env
```sh
conda create -n genespace \
    -c conda-forge -c bioconda \
    r-base \
    r-devtools \
    r-remotes \
    r-ggplot2 \
    r-igraph \
    r-dbscan \
    r-r.utils \
    r-data.table \
    bioconductor-biostrings \
    bioconductor-rtracklayer \
    orthofinder \
    diamond \
    git \
    make \
    gcc_linux-64 \
    gxx_linux-64
```

Check installation:
```sh
orthofinder --version
diamond version
R --version
```

After this is done, activate this conda env and start R:
```sh
conda activate genespace
R
```
Then install GENESPACE in it:
```R
if (!requireNamespace('devtools', quietly = T))
    install.packages('devtools')

devtools::install_github('jtlovell/GENESPACE')

library(GENESPACE)
citation('GENESPACE')
```

### Step 5: Convert .gff into protein fasta
First, check whether sequence names in the GFF exactly match the FASTA headers:
```sh
# from the synteny/assemblies_synteny dir
for dir in */; do

    sp=${dir%/}

    gff="${dir}/${sp}.annotation.gff3"
    fa="${dir}/${sp}.genome.fa"

    echo
    echo "============================================================"
    echo "$sp"
    echo "============================================================"

    if [[ ! -s "$gff" ]]; then
        echo "GFF: MISSING"
        continue
    fi

    if [[ ! -s "$fa" ]]; then
        echo "Genome: MISSING"
        continue
    fi

    echo "--- GFF sequence IDs ---"

    awk -F '\t' '
        $0 !~ /^#/ && NF >= 9 {print $1}
    ' "$gff" \
        | sort -u \
        | head -10

    echo
    echo "--- FASTA sequence IDs ---"

    grep '^>' "$fa" \
        | sed 's/^>//' \
        | cut -d' ' -f1 \
        | head -10

done
```

Before changing