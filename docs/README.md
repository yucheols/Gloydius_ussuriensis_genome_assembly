# Whole-genome assembly of the Ussuri Pitviper (*Gloydius ussuriensis*)
Reference genome assembly of *Gloydius ussuriensis* using PacBio HiFi long-read sequencing, RNA-seq, and Hi-C. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation, https://github.com/amandamarkee/onigra_genome, and https://github.com/amandamarkee/actias-luna-genome. Also, thanks to Amanda Markee (AMNH IZ), Daniela Garcia (AMNH Herp), Jon Hoffman (AMNH Herp), Dylan DeBaun (AMNH Herp), Dean Bobo (AMNH ICG), and Sajesh Singh (AMNH CS) for help and discussions.

The genome sequencing was done on the PacBio Revio system (1 SMRT cell) and RNA sequencing was done on the Illumina NovaSeq X (151bp PE). The individual used for this genome assembly is accessioned at the AMNH Herpetology Collections under the field number AMNH 21010.

__Note:__ I created this documentation in the hopes that my friends and future RGGS students doing genomics may find it useful. If you have any questions about any parts of the assembly process outlined below, please don't hesitate to email me at yshin@amnh.org or post a question in the "issues" section of the repository.

### Workflow

1. __[A quick sanity check on the dataset](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/01_DataCheck.md#1-a-quick-sanity-check-on-the-dataset)__
2. __[*k*-mer analysis of raw reads using jellyfish](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/02_KmerJellyfish.md#2-k-mer-analysis-of-raw-reads-using-jellyfish)__
3. __[Draft genome assembly using hifiasm](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/03_DraftAssembly.md#3-draft-genome-assembly-using-hifiasm)__
4. __[Contamination screening](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/04_ContamCheck.md#4-contamination-screening)__
   - __Screening for potential non-vertebrate contaminants using blobtools__
   - __Identifying and removing mitochondrial contigs from the draft assembly__
   - __Genomewide mean sequencing coverage after contamination screening__
5. __[Draft assembly evaluation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/05_DraftEval.md#5-draft-assembly-evaluation)__ 
   - __Genome completeness using BUSCO__
   - __Genome assembly stats with QUAST__
   - __*k*-mer based assembly evaluation with Merqury__
6. __[Scaffolding through Hi-C data incorporation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/06_Scaffolding.md#6-scaffolding-through-hi-c-data-incorporation)__
   - __Hi-C sequencing overview__
   - __Setup__
   - __Combine sequencing reads across lanes__
   - __Map Hi-C reads to the draft genome__
   - __Scaffolding with YaHS__
   - __Hi-C contact map visualization with Juicer/Juicer Tools__
   - __Assignment of scaffolds to chromosomes and manual assembly curation__
   - __Sex chromosome validation based on sex-specific read coverage patterns__
7. __[Genome annotation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/07_Annotation.md#7-genome-annotation)__
   - __Setup__
   - __RNA read QC (pre-trimming)__
   - __Adapter trimming & post-trimming QC__ 
   - __(Pre-Hi-C) RNA alignment to draft using HiSat2__
   - __(Pre-Hi-C) Draft-guided transcriptome assembly using StringTie__
   - __(Pre-Hi-C) Venom gland transcriptome data__
   - __(Post-Hi-C) Repeat masking (soft masking) using Earl Grey__
   - __(Post-Hi-C) Annotation using funannotate__
   - __(Post-Hi-C) Toxin gene annotation using ToxCodAn-Genome__
8. __[Chromosomal synteny](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/08_Synteny.md#8-chromosomal-synteny)__
   - __Step 1: Setup directory and prepare data downloads__
   - __Step 2: Download assemblies__
9. __[Mitogenome assembly](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/09_MitoAssembly.md#9-mitogenome-assembly)__
    - __"Manual" annotation with MITOS2__
    - __Manual curation of the mitogenome assembled from MitoHiFi__
    - __Submitting mitogenome to GenBank__
10. __[Telomere identification](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/10_Telomere.md#10-telomere-identification)__

11. __[Demographic history](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/11_Demography.md#11-demographic-history)__
    - __Step 1: Environment setup and software installation__
    - __Step 2: Reference genome preparation__
    - __Step 3: Prepare mainland sample BAMs__
    - __Step 4: BAM validation and autosomal coverage QC__
    - __Step 5: Set depth cutoff__
    - __Step 6: Create actual callable-region BEDs__
    - __Step 7: Make the depth-based SMC++ exclusion mask__
    - __Step 8: Step 8: Joint variant calling using BCFtools__