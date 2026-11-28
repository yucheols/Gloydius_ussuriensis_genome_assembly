# Whole-genome assembly of the Ussuri Pitviper (*Gloydius ussuriensis*)
Reference genome assembly of *Gloydius ussuriensis* using PacBio HiFi long-read sequencing, RNA-seq, and Hi-C. Workflow adapted from: https://github.com/danielagarciacobos4/PacBio_GenomeAssembly_annotation, https://github.com/amandamarkee/onigra_genome, and https://github.com/amandamarkee/actias-luna-genome. Also, thanks to Amanda Markee (AMNH IZ), Rebecca Lopez-Anido (AMNH-IZ), Daniela Garcia (AMNH Herp), Jon Hoffman (AMNH Herp), Dylan DeBaun (AMNH Herp), Dean Bobo (AMNH ICG), and Sajesh Singh (AMNH CS) for help and discussions.

The genome sequencing was done on the PacBio Revio system (1 SMRT cell) and RNA sequencing was done on the Illumina NovaSeq X (151bp PE). The Hi-C sequencing was done on three lanes of Illumina NovaSeq X Plus (150bp PE). Additional short-read whole genome sequencing data used for sex chromosome validation were generated on four lanes of Illumina NovaSeq X Plus (150bp PE). The individual used for this genome assembly is accessioned at the AMNH Herpetology Collections under the catalog number AMNH-R-179997 (field number AMNH 21010).

__Note:__ I created this documentation in the hopes that my friends and future RGGS students doing genomics may find it useful. If you have any questions about any parts of the assembly process outlined below, please don't hesitate to email me at yshin@amnh.org or post a question in the "issues" section of the repository.

### Workflow
1. __[A quick sanity check on the dataset](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/01_DataCheck.md#1-a-quick-sanity-check-on-the-dataset)__
2. __[*k*-mer analysis of raw reads using jellyfish](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/02_KmerJellyfish.md#2-k-mer-analysis-of-raw-reads-using-jellyfish)__
3. __[Draft genome assembly using hifiasm](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/03_DraftAssembly.md#3-draft-genome-assembly-using-hifiasm)__
4. __[Contamination screening](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/04_ContamCheck.md#4-contamination-screening)__
5. __[Draft assembly evaluation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/05_DraftEval.md#5-draft-assembly-evaluation)__ 
6. __[Scaffolding through Hi-C data incorporation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/06_Scaffolding.md#6-scaffolding-through-hi-c-data-incorporation)__
7. __[Telomere identification](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/07_Telomere.md#7-telomere-identification)__
8. __[Genome annotation](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/08_Annotation.md#8-genome-annotation)__
9. __[Curation of annotated toxin genes](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/09_ToxinAnnotationEval.md#9-post-hi-c-finalizing-toxin-gene-annotation)__
10. __[Toxin gene microsynteny](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/10_ToxinGeneAnalyses.md#10-toxin-gene-evolution)__
11. __[Toxin gene phylogeny](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/11_ToxinPhylogeny.md#11-toxin-gene-phylogeny)__
12. __[Chromosomal synteny](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/12_ChromosomeSynteny.md#12-chromosomal-macrosynteny)__  
13. __[Historical demography](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/13_Demography.md#13-demographic-history-using-psmc)__
14. __[Circos plot](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/14_CircosPlot.md#14-circos-plot)__
15. __[Repeat landscape plot](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/15_PlotRepeats.md#15-plotting-the-repeat-landscape)__
16. __[Mitogenome assembly](https://github.com/yucheols/Gloydius_ussuriensis_genome_assembly/blob/main/docs/16_MitoAssembly.md#16-mitogenome-assembly)__


