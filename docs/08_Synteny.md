## 8) Chromosomal synteny
### Step 1: Setup directory and prepare data downloads
First, set up a directory for synteny analyses.
```sh
# under the "G_ussuriensis_Chromo" directory
mkdir -p synteny/
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
__NOTE:__ OrthoFinder v3 and above are NOT compatible with GENESPACE. The version used here is v2.5.5.
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
    'orthofinder=2.5.5' \
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

### Step 4: GFF-FASTA sequence ID matching
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

#### Crotalus viridis
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

#### Gloydius shedaoensis
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

#### Xenopeltis unicolor
The id fix for this species is straightforward because the original sequence names are preserved in OriSeqID=.

First, build the original seq name - GWH accession mapping:
```sh
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
' Xenopeltis_unicolor.genome.fa \
> seqid_map.tsv
```

Check the output. We can see that the mapping is clean, with chromosome scale contigs clearly labeled.
```sh
head -30 seqid_map.tsv
```

Now, run this:
```sh
awk -F '\t' '$0 !~ /^#/ {print $1}' \
    Xenopeltis_unicolor.annotation.gff3 \
    | sort -u \
    > gff_ids.txt

cut -f1 seqid_map.tsv \
    | sort -u \
    > mapped_ids.txt

comm -23 gff_ids.txt mapped_ids.txt
```

The last chunk will print this:
```sh
Xuni_Chr01
Xuni_Chr02
Xuni_Chr03
Xuni_Chr04
Xuni_Chr05
Xuni_Chr06
Xuni_Chr07
Xuni_Chr08
Xuni_Chr09
Xuni_Chr10
Xuni_Chr11
Xuni_Chr12
Xuni_Chr13
Xuni_Chr14
Xuni_Chr15
Xuni_Chr16
Xuni_Chr17
Xuni_Chr18
```

This is a simple, clean naming mismatch. The fasta ID is using Xhai, whereas gff ID is using Xuni. Let's check whether changing Xuni to Xhai resolves every remaining gff sequence ID.
```sh
awk -F '\t' '$0 !~ /^#/ {
    x=$1
    sub(/^Xuni_/, "Xhai_", x)
    print x
}' Xenopeltis_unicolor.annotation.gff3 \
    | sort -u \
    > gff_ids_transformed.txt

comm -23 gff_ids_transformed.txt mapped_ids.txt
```

This will not print any outputs, which means that the prefix difference is the entire problem.

Now, create an intermediate cleaned gff:
```sh
awk -F '\t' 'BEGIN{OFS="\t"}
    /^#/ {
        print
        next
    }
    {
        sub(/^Xuni_/, "Xhai_", $1)
        print
    }
' Xenopeltis_unicolor.annotation.gff3 \
> Xenopeltis_unicolor.annotation.idclean.gff3
```

Then convert all those original names to the actual GWH accessions:
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
  Xenopeltis_unicolor.annotation.idclean.gff3 \
> Xenopeltis_unicolor.annotation.seqids_fixed.gff3
```

Now, convert filenames for downstream use:
```sh
mv Xenopeltis_unicolor.annotation.gff3 \
   Xenopeltis_unicolor.annotation.original.gff3

mv Xenopeltis_unicolor.annotation.seqids_fixed.gff3 \
   Xenopeltis_unicolor.annotation.gff3

mkdir fasta_gff_ids_check
mv *.txt *.tsv *.idclean.gff3 *.original.gff3 fasta_gff_ids_check/
```

#### Argyrophis diardii
Resolve in the same way as above.

#### Bothrops insularis
For this species, before trying to map gff Binsu_* chromosome names to fasta CM1480*.1, let's check the NCBI sequence report:
```sh
# inspect one record
head -1 $(find . -name 'sequence_report.jsonl' | head -1) | jq .

# show full fasta header
grep '^>' Bothrops_insularis.genome.fa \
    | head -25 \
    | cat -A
```
The FASTA headers already expose the chromosome identities, and the numbering strongly suggests the published Binsu_* aliases correspond directly to those chromosomes such that:
```sh
Binsu_ma-1  -> CM148030.1   chromosome 1
Binsu_ma-2  -> CM148031.1   chromosome 2
Binsu_ma-3  -> CM148032.1   chromosome 3
Binsu_ma-4  -> CM148033.1   chromosome 4
Binsu_ma-5  -> CM148034.1   chromosome 5
Binsu_ma-6  -> CM148035.1   chromosome 6
Binsu_ma-7  -> CM148036.1   chromosome 7

Binsu_mi-1  -> CM148037.1   chromosome 9
Binsu_mi-2  -> CM148038.1   chromosome 10
...
Binsu_mi-10 -> CM148046.1   chromosome 18

Binsu_Z     -> CM148047.1   chromosome Z
```

Now, let's get every sequence ID in the gff and its maximum feature coordinate:
```sh
awk -F '\t' '
    $0 !~ /^#/ && NF >= 5 {
        if($5 > max[$1])
            max[$1]=$5
    }
    END {
        for(x in max)
            print x, max[x]
    }
' Bothrops_insularis.annotation.gff3 \
| sort -V
```

Also get the chromosome accessions and lengths from the fasta:
```sh
samtools faidx Bothrops_insularis.genome.fa

awk '
    NR <= 18 {
        print $1, $2
    }
' Bothrops_insularis.genome.fa.fai
```

The outputs are enough to confidently resolve naming discrepencies.
```sh
GFF             max feature    FASTA accession    FASTA length
Binsu_ma-1      352151586      CM148030.1         352152128
Binsu_ma-2      279938784      CM148031.1         279938784
Binsu_ma-3      206955912      CM148032.1         206955912
Binsu_ma-4      122534957      CM148033.1         122534993
Binsu_ma-5      100144331      CM148034.1         100144331
Binsu_ma-6       88848289      CM148035.1          88848289
Binsu_ma-7       82184342      CM148036.1          82184716

Binsu_mi-1       24375632      CM148037.1          24375633
Binsu_mi-2       22207381      CM148038.1          22207440
Binsu_mi-3       19746646      CM148039.1          19746646
Binsu_mi-4       18317824      CM148040.1          18318580
Binsu_mi-5       17755719      CM148041.1          17755719
Binsu_mi-6       17602135      CM148042.1          17602211
Binsu_mi-7       14071518      CM148043.1          14071536
Binsu_mi-8       13870028      CM148044.1          13870031
Binsu_mi-9       13660605      CM148045.1          13662181
Binsu_mi-10      11461905      CM148046.1          11461905

Binsu_Z         138394576      CM148047.1         138394576
```

Directly create the mapping file:
```sh
printf "Binsu_ma-1\tCM148030.1\n\
Binsu_ma-2\tCM148031.1\n\
Binsu_ma-3\tCM148032.1\n\
Binsu_ma-4\tCM148033.1\n\
Binsu_ma-5\tCM148034.1\n\
Binsu_ma-6\tCM148035.1\n\
Binsu_ma-7\tCM148036.1\n\
Binsu_mi-1\tCM148037.1\n\
Binsu_mi-2\tCM148038.1\n\
Binsu_mi-3\tCM148039.1\n\
Binsu_mi-4\tCM148040.1\n\
Binsu_mi-5\tCM148041.1\n\
Binsu_mi-6\tCM148042.1\n\
Binsu_mi-7\tCM148043.1\n\
Binsu_mi-8\tCM148044.1\n\
Binsu_mi-9\tCM148045.1\n\
Binsu_mi-10\tCM148046.1\n\
Binsu_Z\tCM148047.1\n" > seqid_map.tsv
```

... and verify the script below prints nothing.
```sh
awk -F '\t' '$0 !~ /^#/ {print $1}' \
    Bothrops_insularis.annotation.gff3 \
    | sort -u \
    > gff_ids.txt

cut -f1 seqid_map.tsv | sort -u > mapped_ids.txt

comm -23 gff_ids.txt mapped_ids.txt
```

Now, create fixed gff:
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
  Bothrops_insularis.annotation.gff3 \
> Bothrops_insularis.annotation.seqids_fixed.gff3
```

Change file names:
```sh
mv Bothrops_insularis.annotation.gff3 \
   Bothrops_insularis.annotation.original.gff3

mv Bothrops_insularis.annotation.seqids_fixed.gff3 \
   Bothrops_insularis.annotation.gff3

mkdir fasta_gff_ids_check
mv *.txt *.tsv *.original.gff3 fasta_gff_ids_check/
```

#### Cerastes gasperettii
Use the same strategy as Bothrops.

#### Elaphe schrenckii
Same deal as above.

### Step 5: Convert annotations into protein fasta
Some comparison taxa had protein fasta already downloaded alongside the genome fasta and annotation. However, some species only had genome fasta and annotation files. Run the slurm script below. This script will look for annotation files in the order GFF3 - GFF - GTF, and either generate the protein fasta or reuse it if it already existed.
```sh
synteny_prep_protein.sh
```

### Step 6: Setup for GENESPACE analyses
```sh
DIR="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

mkdir -p \
    "${DIR}/bed" \
    "${DIR}/peptide" \
    "${DIR}/mapping" \
```

Next, we will prepare a representative protein set that will be used as an input for GENESPACE analyses. Our annotation has multiple transcript/protein isoforms for some genes, while for macrosynteny we generally want one locus per protein. GENESPACE compares genes across genomes, and feeding every isoform would make one biological gene look like several nearby genes. Therefore, we will run the Python script below to select one representative (i.e. the longest) protein isoform per gene.
```sh
make_genespace_representative.py
```

The script will read the final GFF3:
```sh
Gloydius_ussuriensis_AMNH_21010.gff3
```

Which has lines like:
```sh
gene  ID=FUN_000001
mRNA  ID=FUN_000001-T1;Parent=FUN_000001
mRNA  ID=FUN_000001-T2;Parent=FUN_000001
```

From this, it will learn that:
```sh
FUN_000001-T1 = FUN_000001
FUN_000001-T2 = FUN_000001
```

Then it will read the protein fasta and select the longest protein isoform.

Run the script like this:
```sh
# activate conda env to access python
conda activate synteny_qc

# set paths
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"
FUN="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/funannotate"

# run script
python3 "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/make_genespace_representative.py" \
    "${FUN}/annotate_results/Gloydius_ussuriensis_AMNH_21010.gff3" \
    "${FUN}/annotate_results/Gloydius_ussuriensis_AMNH_21010.proteins.fa" \
    "${GS}/bed/Gloydius_ussuriensis.bed" \
    "${GS}/peptide/Gloydius_ussuriensis.fa" \
    "${GS}/mapping/Gloydius_ussuriensis.longest_isoform.tsv"
```

The output will print like this:
```sh
============================================================
GENESPACE representative-protein preparation
============================================================
GFF genes:                         23,157
GFF transcripts:                   28,000
Input proteins:                    28,000
Genes represented by proteins:     22,534
GENESPACE loci written:            22,534
Genes without protein:              623
Proteins without transcript map:    0
Representative genes without coord: 0
============================================================
```

Let's also run a BED-fasta identity check:
```sh
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

echo "===== counts ====="

printf "BED loci: "
wc -l < "${GS}/bed/Gloydius_ussuriensis.bed"

printf "Protein sequences: "
grep -c '^>' "${GS}/peptide/Gloydius_ussuriensis.fa"

printf "Mapping rows: "
awk 'NR>1' "${GS}/mapping/Gloydius_ussuriensis.longest_isoform.tsv" | wc -l


echo
echo "===== duplicate BED IDs ====="

cut -f4 "${GS}/bed/Gloydius_ussuriensis.bed" \
    | sort \
    | uniq -d \
    | head


echo
echo "===== duplicate FASTA IDs ====="

grep '^>' "${GS}/peptide/Gloydius_ussuriensis.fa" \
    | sed 's/^>//; s/[[:space:]].*$//' \
    | sort \
    | uniq -d \
    | head


echo
echo "===== BED IDs absent from protein FASTA ====="

comm -23 \
    <(cut -f4 "${GS}/bed/Gloydius_ussuriensis.bed" | sort) \
    <(grep '^>' "${GS}/peptide/Gloydius_ussuriensis.fa" | sed 's/^>//; s/[[:space:]].*$//' | sort) \
    | head


echo
echo "===== protein IDs absent from BED ====="

comm -13 \
    <(cut -f4 "${GS}/bed/Gloydius_ussuriensis.bed" | sort) \
    <(grep '^>' "${GS}/peptide/Gloydius_ussuriensis.fa" | sed 's/^>//; s/[[:space:]].*$//' | sort) \
    | head
```

This will show:
```sh
BED loci:          22,534
Protein sequences: 22,534
Mapping rows:      22,534
```
This means one gene locus - one BED entry - one representative protein match.

Now, we need to do the same for all comparison taxa. Let's inspect the exact annotation/protein filenames and formats in the comparison folders.
```sh
BASE="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"

for d in "${BASE}"/*; do

    [ -d "${d}" ] || continue

    echo
    echo "============================================================"
    echo "$(basename "${d}")"
    echo "============================================================"

    find "${d}" \
        -maxdepth 1 \
        -type f \
        \( -iname "*.gff" \
        -o -iname "*.gff3" \
        -o -iname "*.gtf" \
        -o -iname "*.faa" \
        -o -iname "*.fa" \
        -o -iname "*.fasta" \
        -o -iname "*.pep" \) \
        -printf '%f\n' \
        | sort

done
```

Let's also check how the protein IDs correspond to the annotation IDs. Check the protein headers first:
```sh
BASE="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"

for d in "${BASE}"/*; do

    [ -d "${d}" ] || continue

    echo
    echo "============================================================"
    echo "$(basename "${d}")"
    echo "============================================================"

    for f in "${d}"/*.faa; do

        [ -f "${f}" ] || continue

        echo "FILE: $(basename "${f}")"
        grep '^>' "${f}" | head -3

    done

done
```

Do the same for the annotations:
```sh
BASE="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/assemblies_synteny"

for d in "${BASE}"/*; do

    [ -d "${d}" ] || continue

    echo
    echo "============================================================"
    echo "$(basename "${d}")"
    echo "============================================================"

    for f in "${d}"/*.gff3 "${d}"/*.gtf; do

        [ -f "${f}" ] || continue

        echo "FILE: $(basename "${f}")"

        awk -F'\t' '
            $0 !~ /^#/ {
                print
                n++
                if (n == 8) exit
            }
        ' "${f}"

    done

done
```

Running these two scripts will reveal that the annotation/protein relationships differ substantially among species. Therefore, let's run audit the ID mapping first:
```sh
audit_comparison_annotations.py
``` 

Run the script:
```sh
# set path
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

# run
python3 "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/audit_comparison_annotations.py" \
    | tee "${GS}/comparison_annotation_ID_audit.txt"
```

This will show that all 11 comparison protein FASTAs can be mapped back to their annotations cleanly. The comparison dataset contains three mapping patterns, and we can use the batch converter script below to create a clean GENESPACE input for the comparison taxa.
```sh
# set path
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

# run
python3 "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/Python/make_comparison_genespace_inputs.py" \
    | tee "${GS}/comparison_genespace_preprocessing.log"
``` 

The results look clean. All 11 comparison genomes passed the critical mapping test, with 0 unmapped proteins and 0 ambiguous mappings.

Before moving on to the actual GENESPACE run, let's do one global QC across all 12 genomes (11 comparison taxa genome + 1 *G. ussuriensis*):
```sh
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

printf "species\tBED\tproteins\tdupBED\tdupFASTA\tBED_not_FASTA\tFASTA_not_BED\tbad_coords\tseqIDs\n"

for bed in "${GS}"/bed/*.bed; do

    sp=$(basename "${bed}" .bed)
    pep="${GS}/peptide/${sp}.fa"

    nbed=$(wc -l < "${bed}")
    npep=$(grep -c '^>' "${pep}")

    dupbed=$(
        cut -f4 "${bed}" |
        sort |
        uniq -d |
        wc -l
    )

    dupfa=$(
        grep '^>' "${pep}" |
        sed 's/^>//; s/[[:space:]].*$//' |
        sort |
        uniq -d |
        wc -l
    )

    bed_not_fasta=$(
        comm -23 \
            <(cut -f4 "${bed}" | sort) \
            <(grep '^>' "${pep}" |
              sed 's/^>//; s/[[:space:]].*$//' |
              sort) |
        wc -l
    )

    fasta_not_bed=$(
        comm -13 \
            <(cut -f4 "${bed}" | sort) \
            <(grep '^>' "${pep}" |
              sed 's/^>//; s/[[:space:]].*$//' |
              sort) |
        wc -l
    )

    badcoords=$(
        awk '
            $2 < 0 || $3 <= $2 {n++}
            END {print n+0}
        ' "${bed}"
    )

    nseq=$(
        cut -f1 "${bed}" |
        sort -u |
        wc -l
    )

    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${sp}" \
        "${nbed}" \
        "${npep}" \
        "${dupbed}" \
        "${dupfa}" \
        "${bed_not_fasta}" \
        "${fasta_not_bed}" \
        "${badcoords}" \
        "${nseq}"

done | column -t
```  
The results will show:
```sh
BED count = protein count
duplicate BED IDs = 0
duplicate FASTA IDs = 0
BED IDs missing from FASTA = 0
FASTA IDs missing from BED = 0
bad coordinates = 0
```
This means that the GENESPACE input prep has been successful.

### Step 7: GENESPACE run
__NOTE:__ Before running GENESPACE, I noticed that *C. viridis* had lots of duplicated transcript coordinates that needed further correction. Therefore, I excluded this species from the initial GENESPACE run.

First, let's check whether the required packages are installed:
```sh
conda activate genespace

command -v Rscript
command -v orthofinder
command -v diamond
command -v MCScanX_h

Rscript -e 'cat("GENESPACE ", as.character(packageVersion("GENESPACE")), "\n")'
```

This will show that MCScanX_h is not installed. Let's install it in the genespace conda env:
```sh
conda install -c conda-forge -c bioconda mcscanx
command -v MCScanX_h
MCScanX_h 2>&1 | head
```

After this, create an R script to run GENESPACE:
```r
# ============================================================
# GENESPACE macrosynteny analysis
#
# This run uses 11 snake genomes
# NOTE: Crotalus viridis excluded pending annotation cleanup
#
# This script is run on AMNH Mendel HPC
# ============================================================

library(GENESPACE)


# ------------------------------------------------------------
# working directory
# ------------------------------------------------------------

wd <- '/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE'


# ------------------------------------------------------------
# genomes
#
# Order is approximately phylogenetic for convenient plotting
# Crotalus_viridis is deliberately excluded
# ------------------------------------------------------------

genomeIDs <- c('Argyrophis_diardii',
               'Xenopeltis_unicolor',
               'Candoia_aspera',
               'Elaphe_schrenckii',
               'Naja_naja',
               'Cerastes_gasperettii',
               'Vipera_berus',
               'Bothrops_insularis',
               'Crotalus_adamanteus',
               'Gloydius_shedaoensis',
               'Gloydius_ussuriensis')


# ------------------------------------------------------------
# computational resources
# ------------------------------------------------------------

nCores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK', unset = '48'))


# ------------------------------------------------------------
# path to external software
# ------------------------------------------------------------

path2orthofinder <- '/home/yshin/mendel-nas1/miniconda3/envs/genespace/bin/orthofinder'
path2diamond <- '/home/yshin/mendel-nas1/miniconda3/envs/genespace/bin/diamond'

# GENESPACE expects the directory containing MCScanX_h
path2mcscanx <- '/home/yshin/mendel-nas1/miniconda3/envs/genespace/bin'


# ------------------------------------------------------------
# print run information
# ------------------------------------------------------------

cat('\n')
cat('============================================================\n')
cat('GENESPACE genome macrosynteny analysis\n')
cat('============================================================\n')

cat('GENESPACE version: ', as.character(packageVersion('GENESPACE')), '\n', sep = '')

cat('Working directory: ', wd, '\n', sep = '')
cat('Cores: ', nCores, '\n', sep = '')

cat('\nGenome IDs:\n')
print(genomeIDs)

cat('\nExternal software:\n')
cat('OrthoFinder: ', path2orthofinder, '\n', sep = '')
cat('DIAMOND:     ', path2diamond, '\n', sep = '')
cat('MCScanX dir: ', path2mcscanx, '\n', sep = '')


# ------------------------------------------------------------
# check required input files
# ------------------------------------------------------------

bedFiles <- file.path(wd, 'bed', paste0(genomeIDs, '.bed'))
pepFiles <- file.path(wd, 'peptide', paste0(genomeIDs, '.fa'))


if (!all(file.exists(bedFiles))) {
  
  stop(paste('Missing BED files:',
             paste(bedFiles[!file.exists(bedFiles)],
                   collapse = '\n')))
}


if (!all(file.exists(pepFiles))) {
  
  stop(paste('Missing peptide files:',
             paste(pepFiles[!file.exists(pepFiles)],
                   collapse = '\n')))
}


# explicit safeguard

if ('Crotalus_viridis' %in% genomeIDs) {
  stop('Crotalus_viridis should not be included in this run.')
}


cat('\nPASS: all 11 BED and peptide files found.\n')


# ------------------------------------------------------------
# initialize GENESPACE
#
# ploidy = 1:
# each assembly represents one haploid chromosome complement
#
# useHOGs = T:
# use hierarchical orthogroups from OrthoFinder
#
# other synteny parameters remain at GENESPACE defaults
# ------------------------------------------------------------

gpar <- init_genespace(wd = wd, genomeIDs = genomeIDs, ploidy = 1,
                       path2orthofinder = path2orthofinder,
                       path2diamond = path2diamond,
                       path2mcscanx = path2mcscanx,
                       useHOGs = T,
                       nCores = nCores,
                       dotplots = 'check')


# ------------------------------------------------------------
# save initialization parameters
# ------------------------------------------------------------

saveRDS(gpar, file.path(wd, 'GENESPACE_parameters_11snake.rds'))


cat('\n')
cat('============================================================\n')
cat('GENESPACE initialization passed\n')
cat('============================================================\n\n')


# ------------------------------------------------------------
# run complete GENESPACE pipeline
# ------------------------------------------------------------

out <- run_genespace(gsParam = gpar)


# ------------------------------------------------------------
# save final object
# ------------------------------------------------------------

saveRDS(out, file.path(wd, 'GENESPACE_results_11snake.rds'))


# ------------------------------------------------------------
# session information
# ------------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(wd, 'GENESPACE_sessionInfo_11snake.txt'))


cat('\n')
cat('============================================================\n')
cat('GENESPACE RUN COMPLETED\n')
cat('============================================================\n')
```

Run this script as a SLURM job:
```sh
#!/bin/bash
#SBATCH --job-name=genespace
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=168:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yshin@amnh.org
#SBATCH --output=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.out
#SBATCH --error=/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/outfiles/slurm-%x_%j.err


# ============================================================
# GENESPACE macrosynteny
# 11 snake genomes
# ============================================================

# activate conda environment
source ~/.bash_profile
conda activate genespace

set -euo pipefail

# set path
GS="/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE"

# print helpful information
echo "============================================================"
echo "GENESPACE genome macrosynteny"
echo "============================================================"
echo "Date:   $(date)"
echo "Node:   $(hostname)"
echo "Job ID: ${SLURM_JOB_ID}"
echo "CPUs:   ${SLURM_CPUS_PER_TASK}"
echo


# ------------------------------------------------------------
# software
# ------------------------------------------------------------

echo "===== SOFTWARE ====="

echo "Rscript:"
command -v Rscript

echo "OrthoFinder:"
command -v orthofinder

echo "DIAMOND:"
command -v diamond

echo "MCScanX_h:"
command -v MCScanX_h

echo

Rscript -e '
cat(
    "GENESPACE ",
    as.character(packageVersion("GENESPACE")),
    "\n"
)
'


# ------------------------------------------------------------
# input summary
# ------------------------------------------------------------

echo
echo "===== INPUT FILES ====="

for sp in \
    Argyrophis_diardii \
    Xenopeltis_unicolor \
    Candoia_aspera \
    Elaphe_schrenckii \
    Naja_naja \
    Cerastes_gasperettii \
    Vipera_berus \
    Bothrops_insularis \
    Crotalus_adamanteus \
    Gloydius_shedaoensis \
    Gloydius_ussuriensis
do

    bed="${GS}/bed/${sp}.bed"
    pep="${GS}/peptide/${sp}.fa"

    if [[ ! -s "${bed}" ]]; then
        echo "ERROR: missing BED: ${bed}"
        exit 1
    fi

    if [[ ! -s "${pep}" ]]; then
        echo "ERROR: missing peptide FASTA: ${pep}"
        exit 1
    fi

    nbed=$(wc -l < "${bed}")
    npep=$(grep -c '^>' "${pep}")

    echo -e "${sp}\tBED=${nbed}\tproteins=${npep}"

done


# ------------------------------------------------------------
# run
# ------------------------------------------------------------

echo "============================================================"
echo "STARTING GENESPACE"
echo "============================================================"

Rscript \
    "/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/scripts/R/genespace.R"


echo
echo "============================================================"
echo "GENESPACE FINISHED"
echo "Date: $(date)"
echo "============================================================"
```
__NOTE:__ The first GENESPACE run failed during DIAMOND database creation. I ran the script below to see if the protein files contained any invalid amino acid caharcters:
```sh
GS='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE'

for sp in \
    Argyrophis_diardii \
    Xenopeltis_unicolor \
    Candoia_aspera \
    Elaphe_schrenckii \
    Naja_naja \
    Cerastes_gasperettii \
    Vipera_berus \
    Bothrops_insularis \
    Crotalus_adamanteus \
    Gloydius_shedaoensis \
    Gloydius_ussuriensis
do

    echo "===== ${sp} ====="

    grep -v '^>' "${GS}/peptide/${sp}.fa" \
        | tr -d '[:space:]ABCDEFGHIKLMNPQRSTVWXYZ*abcdefghiklmnpqrstvwxyz' \
        | fold -w1 \
        | sort \
        | uniq -c

done
```
This revealed that A. diardii had 465 periods (.), X. unicolor had 371, B. insiularis had 3, and G. shedaoensis had 1151. The periods are not recognized by DIAMOND and this is what caused the job failure. The "U" characters in C. aspera and V. berus are harmless. Let's replace periods with "X", which deisgnates ambiguous amino acids:
```sh
GS='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE'
for sp in \
    Argyrophis_diardii \
    Xenopeltis_unicolor \
    Bothrops_insularis \
    Gloydius_shedaoensis
do

    sed -i '/^>/! s/\./X/g' \
        "${GS}/peptide/${sp}.fa"

done
```