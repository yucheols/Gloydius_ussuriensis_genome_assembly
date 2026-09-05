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

Before changing anything, let's verify matches between gff IDs and fasta IDs. For cases where both genome and gff files were downloaded directly from GenBank (C. adamanteus. C. aspera, and V. berus), the two IDs should match already.

Verify like this:
```sh
# from the assemblies_synteny dir
for sp in \
    Candoia_aspera \
    Crotalus_adamanteus \
    Vipera_berus
do

    echo
    echo "===== ${sp} ====="

    awk -F '\t' \
        '$0 !~ /^#/ && NF >= 9 {print $1}' \
        "${sp}/${sp}.annotation.gff3" \
        | sort -u \
        > "${sp}/gff_seqids.txt"

    grep '^>' "${sp}/${sp}.genome.fa" \
        | sed 's/^>//' \
        | awk '{print $1}' \
        | sort -u \
        > "${sp}/fasta_seqids.txt"

    echo -n "GFF IDs:   "
    wc -l < "${sp}/gff_seqids.txt"

    echo -n "FASTA IDs: "
    wc -l < "${sp}/fasta_seqids.txt"

    echo "GFF IDs absent from FASTA:"
    comm -23 \
        "${sp}/gff_seqids.txt" \
        "${sp}/fasta_seqids.txt"

done
```

We expect the output to look like:
```sh
===== Candoia_aspera =====
GFF IDs:   149
FASTA IDs: 149
GFF IDs absent from FASTA:

===== Crotalus_adamanteus =====
GFF IDs:   27
FASTA IDs: 27
GFF IDs absent from FASTA:

===== Vipera_berus =====
GFF IDs:   478
FASTA IDs: 478
GFF IDs absent from FASTA:
```

#### Fixing slight ID mismatches
Crotalus viridis and Naja naja have slight mismatches between fasta and gff IDs (e.g., Cvir_CM012306.1 and CM012306.1 in C. viridis). For C. viridis, create a corrected GFF:
```sh
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ {
        print
        next
    }
    {
        sub(/^Cvir_/, "", $1)
        print
    }
' Crotalus_viridis/Crotalus_viridis.annotation.gff3 \
> Crotalus_viridis/Crotalus_viridis.annotation.seqids_fixed.gff3
```

...and verify the fix:
```sh
awk -F '\t' \
    '$0 !~ /^#/ {print $1}' \
    Crotalus_viridis/Crotalus_viridis.annotation.seqids_fixed.gff3 \
    | sort -u \
    > Crotalus_viridis/gff_fixed_seqids.txt

grep '^>' Crotalus_viridis/Crotalus_viridis.genome.fa \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | sort -u \
    > Crotalus_viridis/fasta_seqids.txt

comm -23 \
    Crotalus_viridis/gff_fixed_seqids.txt \
    Crotalus_viridis/fasta_seqids.txt
```
We expect these chunks to not print anything. Let's also preserve the original gff and promote the fixed gff for downstream use by changing file names:
```sh
mv Crotalus_viridis.annotation.gff3 \
   Crotalus_viridis.annotation.original.gff3

mv Crotalus_viridis.annotation.seqids_fixed.gff3 \
   Crotalus_viridis.annotation.gff3

mkdir fasta_gff_ids_check
mv *.txt *.tsv *.idclean.gff3 *.original.gff3 fasta_gff_ids_check/
```


Repeat the same for N. naja:
```sh
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ {
        print
        next
    }
    {
        sub(/^Nnaj_/, "", $1)
        sub(/_np121212$/, "", $1)
        print
    }
' Naja_naja/Naja_naja.annotation.gff3 \
> Naja_naja/Naja_naja.annotation.seqids_fixed.gff3
```
...and again verify and check nothing prints to the console.
```sh
awk -F '\t' \
    '$0 !~ /^#/ {print $1}' \
    Naja_naja/Naja_naja.annotation.seqids_fixed.gff3 \
    | sort -u \
    > Naja_naja/gff_fixed_seqids.txt

grep '^>' Naja_naja/Naja_naja.genome.fa \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | sort -u \
    > Naja_naja/fasta_seqids.txt

comm -23 \
    Naja_naja/gff_fixed_seqids.txt \
    Naja_naja/fasta_seqids.txt
```

Change the file name for downstream use:
```sh
mv Naja_naja.annotation.gff3 \
   Naja_naja.annotation.original.gff3

mv Naja_naja.annotation.seqids_fixed.gff3 \
   Naja_naja.annotation.gff3

mkdir fasta_gff_ids_check
mv *.txt *.tsv *.idclean.gff3 *.original.gff3 fasta_gff_ids_check/
```

#### Assemblies downloaded from NGDC
For the assemblies downloaded from NGDC, run the following. cat -A will expose tabs and the complete chromosome labels in the FASTA headers.
```sh
for sp in \
    Argyrophis_diardii \
    Gloydius_shedaoensis \
    Xenopeltis_unicolor
do

    echo
    echo "===== ${sp} ====="

    grep '^>' "${sp}/${sp}.genome.fa" \
        | head -30 \
        | cat -A

done
```

The output will tell us that, for example, Gshe_Chr01 is GWHBWDW00000001.

##### Gloydius shedaoensis
Run this script first to build a seq id mapping directly from fasta headers:
```sh
cd Gloydius_shedaoensis

awk -F '\t' '
    /^>/ {
        acc=$1
        sub(/^>/, "", acc)

        for(i=1;i<=NF;i++) {
            if($i ~ /^OriSeqID=/) {
                orig=$i
                sub(/^OriSeqID=/, "", orig)
                print orig "\t" acc
            }
        }
    }
' Gloydius_shedaoensis.genome.fa \
> seqid_map.tsv
```

Check:
```sh
head seqid_map.tsv
```
....and we should see:
```sh
Gshe_Chr01      GWHBWDU00000001
Gshe_Chr02      GWHBWDU00000002
Gshe_Chr03      GWHBWDU00000003
Gshe_Chr04      GWHBWDU00000004
Gshe_Chr05      GWHBWDU00000005
Gshe_Chr06      GWHBWDU00000006
Gshe_Chr07      GWHBWDU00000007
Gshe_Chr08      GWHBWDU00000008
Gshe_Chr09      GWHBWDU00000009
Gshe_Chr10      GWHBWDU00000010
```

Now verify every GFF seqid has a mapping:
```sh
awk -F '\t' '$0 !~ /^#/ {print $1}' \
    Gloydius_shedaoensis.annotation.gff3 \
    | sort -u \
    > gff_ids.txt

cut -f1 seqid_map.tsv | sort -u > mapped_ids.txt

comm -23 gff_ids.txt mapped_ids.txt
```

Here, no output is a good sign. But the last chunk will print something like: 
```sh
Gshe_Contig5_ERROPOS18435237+
```

This name is suspicious. Let's run a series of validation steps:
```sh
# check whether the assembly contains Gshe_Contig5
grep -F 'OriSeqID=Gshe_Contig5' Gloydius_shedaoensis.genome.fa
grep -F $'Gshe_Contig5\t' seqid_map.tsv

# inspect every GFF record using the weird ID
grep -n -F 'Gshe_Contig5_ERROPOS18435237+' \
    Gloydius_shedaoensis.annotation.gff3
```

This will show that the assembly actually contains a contig named Gshe_Contig5_ERROPOS18435237 which maps to GWHBWDU00000198. gff file is instead calling the same sequence Gshe_Contig5_ERROPOS18435237+

First, strip the trailing plus sign and create a gff file with cleaned id:
```sh
sed 's/^Gshe_Contig5_ERROPOS18435237+/Gshe_Contig5_ERROPOS18435237/' \
    Gloydius_shedaoensis.annotation.gff3 \
    > Gloydius_shedaoensis.annotation.idclean.gff3
```

Now confirm that all original GFF IDs can be mapped. We want no output:
```sh
awk -F '\t' '$0 !~ /^#/ {print $1}' \
    Gloydius_shedaoensis.annotation.idclean.gff3 \
    | sort -u \
    > gff_ids_clean.txt

cut -f1 seqid_map.tsv | sort -u > mapped_ids.txt

comm -23 gff_ids_clean.txt mapped_ids.txt
```

The empty output from the last chunk is especially informative, It means that means every sequence ID used by the cleaned G. shedaoensis gff now has a corresponding entry in the FASTA-derived mapping.

Now, convert the cleaned gff from the original Gshe_* names to the GWH accession names:
```sh
awk -F '\t' 'BEGIN{OFS="\t"}

    NR==FNR {
        map[$1]=$2
        next
    }

    /^#/ {
        print
        next
    }

    {
        if(!($1 in map)) {
            print "ERROR: unmapped sequence ID: " $1 > "/dev/stderr"
            exit 1
        }

        $1=map[$1]
        print
    }

' seqid_map.tsv \
  Gloydius_shedaoensis.annotation.idclean.gff3 \
> Gloydius_shedaoensis.annotation.seqids_fixed.gff3
```

Then run the final compatibility check against the actual genome fasta:
```sh
awk -F '\t' '$0 !~ /^#/ {print $1}' \
    Gloydius_shedaoensis.annotation.seqids_fixed.gff3 \
    | sort -u \
    > fixed_gff_ids.txt

grep '^>' Gloydius_shedaoensis.genome.fa \
    | sed 's/^>//' \
    | awk '{print $1}' \
    | sort -u \
    > fasta_ids.txt

comm -23 fixed_gff_ids.txt fasta_ids.txt
```
Again, we want no output.

Let's preserve the original and promote the fixed file for downstream use by changing the file names:
```sh
mv Gloydius_shedaoensis.annotation.gff3 \
   Gloydius_shedaoensis.annotation.original.gff3

mv Gloydius_shedaoensis.annotation.seqids_fixed.gff3 \
   Gloydius_shedaoensis.annotation.gff3

mkdir fasta_gff_ids_check
mv *.txt *.tsv *.idclean.gff3 *.original.gff3 fasta_gff_ids_check/
```

##### Xenopeltis unicolor
The id fix for this species is straightforward 