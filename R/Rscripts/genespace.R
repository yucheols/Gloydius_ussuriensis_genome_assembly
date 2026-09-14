# ============================================================
# GENESPACE macrosynteny analysis
#
# 11 snake genomes
#
# Input preparation:
#   - one longest representative protein per biological locus
#   - chromosome-scale pseudomolecules only
#   - unplaced/unlocalized sequences excluded
#   - explicitly identified W chromosomes excluded
#   - Z chromosomes retained
#   - numerically labelled potential sex chromosomes retained
#
# NOTE:
#   Crotalus viridis excluded pending annotation cleanup
#
# This script is run on AMNH Mendel HPC
# ============================================================


# ------------------------------------------------------------
# clean R environment
# ------------------------------------------------------------

rm(list = ls(all.names = T))
gc()


# ------------------------------------------------------------
# load GENESPACE
# ------------------------------------------------------------

library(GENESPACE)


# ------------------------------------------------------------
# working directory
# ------------------------------------------------------------

wd <- '/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/GENESPACE'


# ------------------------------------------------------------
# genomes
#
# Order is approximately phylogenetic for convenient plotting.
#
# Crotalus_viridis is deliberately excluded because its
# annotation contains extensive duplicated GeneWise models
# requiring separate cleanup.
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
# paths to external software
# ------------------------------------------------------------

path2orthofinder <- '/home/yshin/mendel-nas1/miniconda3/envs/genespace/bin/orthofinder'
path2diamond <- '/home/yshin/mendel-nas1/miniconda3/envs/genespace/bin/diamond'

# GENESPACE expects the directory containing MCScanX_h
path2mcscanx <- '/home/yshin/mendel-nas1/miniconda3/envs/genespace/bin'


# ------------------------------------------------------------
# output file names
#
# Distinguish this chromosome-only run from earlier
# exploratory GENESPACE analyses.
# ------------------------------------------------------------

parameterFile <- file.path(wd, 'GENESPACE_parameters_11snake_chrOnly.rds')
resultFile <- file.path(wd, 'GENESPACE_results_11snake_chrOnly.rds')
sessionFile <- file.path(wd, 'GENESPACE_sessionInfo_11snake_chrOnly.txt')


# ============================================================
# print run information
# ============================================================

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


# ============================================================
# basic safeguards
# ============================================================

if (length(genomeIDs) != 11) {
  stop(paste('Expected 11 genomes, found', length(genomeIDs)))
}

if (anyDuplicated(genomeIDs)) {
  stop('Duplicate genome IDs found in genomeIDs.')
}

if ('Crotalus_viridis' %in% genomeIDs) {
  stop('Crotalus_viridis should not be included in this run.')
}

if (!dir.exists(wd)) {
  stop(paste('Working directory does not exist:', wd))
}


# ============================================================
# check external programs
# ============================================================

if (!file.exists(path2orthofinder)) {
  stop(paste('OrthoFinder executable not found:', path2orthofinder))
}

if (!file.exists(path2diamond)) {
  stop(paste('DIAMOND executable not found:', path2diamond))
}


mcscanxExecutable <- file.path(path2mcscanx, 'MCScanX_h')

if (!file.exists(mcscanxExecutable)) {
  stop(paste('MCScanX_h executable not found:', mcscanxExecutable))
}

cat('\nPASS: external software paths found.\n')


# ============================================================
# check required GENESPACE input files
# ============================================================

bedFiles <- file.path(wd, 'bed', paste0(genomeIDs, '.bed'))
pepFiles <- file.path(wd, 'peptide', paste0(genomeIDs, '.fa'))

missingBeds <- bedFiles[!file.exists(bedFiles)]
missingPeptides <- pepFiles[!file.exists(pepFiles)]


if (length(missingBeds) > 0) {
  stop(paste('Missing BED files:\n', paste(missingBeds, collapse = '\n')))
}


if (length(missingPeptides) > 0) {
  stop(paste('Missing peptide files:\n', paste(missingPeptides, collapse = '\n')))
}

cat('PASS: all 11 BED and peptide files found.\n')


# ============================================================
# check BED files
# ============================================================

bedData <- lapply(bedFiles, function(x) {
                    read.delim(x, header = F, stringsAsFactors = F, sep = '\t')
                    })

names(bedData) <- genomeIDs


# ------------------------------------------------------------
# BED files must contain exactly four columns
# ------------------------------------------------------------

badBedColumns <- names(bedData)[vapply(bedData, ncol, integer(1)) != 4]

if (length(badBedColumns) > 0) {
  stop(paste('BED files do not contain exactly four columns:',
             paste(badBedColumns, collapse = ', ')))
}

cat('PASS: all BED files contain four columns.\n')


# ============================================================
# verify explicit W chromosomes are absent
# ============================================================

wCheck <- vapply(bedData,
                 function(x) {
                   any(x[[1]] == 'chrW')
                   },
                 logical(1))


if (any(wCheck)) {
  
  stop(paste('Explicit chrW found in BED file(s):', 
             paste(names(wCheck)[wCheck], collapse = ', ')))
}

cat('PASS: no explicitly labelled chrW remains in BED files.\n')


# ============================================================
# report chromosome counts and locus counts
# ============================================================

cat('\n')
cat('============================================================\n')
cat('GENESPACE INPUT SUMMARY\n')
cat('============================================================\n')

for (species in genomeIDs) {
  
  dat <- bedData[[species]]
  chromosomeCount <- length(unique(dat[[1]]))
  locusCount <- nrow(dat)
  
  cat(sprintf('%-25s chromosomes=%2d  loci=%6d\n', species, chromosomeCount, locusCount))
}


# ============================================================
# initialize GENESPACE
#
# ploidy = 1:
#   each genome assembly represents one haploid chromosome
#   complement for the purposes of GENESPACE.
#
# useHOGs = T:
#   use hierarchical orthogroups generated by OrthoFinder.
#
# orthofinderInBlk = F:
#   run OrthoFinder normally rather than using GENESPACE
#   syntenic blocks to constrain OrthoFinder.
#
# onlySameChrs = F:
#   do NOT restrict synteny searches to chromosomes carrying
#   the same chromosome label. This is important because the
#   source chromosome numbers are not being treated as
#   pre-assigned homologous chromosomes across species.
#
# dotplots = 'check':
#   generate diagnostic dotplots where practical.
#
# Other GENESPACE synteny parameters remain at defaults.
# ============================================================

gpar <- init_genespace(wd = wd,
                       genomeIDs = genomeIDs,
                       ploidy = 1,
                       path2orthofinder = path2orthofinder,
                       path2diamond = path2diamond,
                       path2mcscanx = path2mcscanx,
                       orthofinderInBlk = F,
                       useHOGs = T,
                       onlySameChrs = F,
                       nCores = nCores,
                       dotplots = 'check')


# ============================================================
# save initialization parameters
# ============================================================

saveRDS(gpar, parameterFile)

cat('\n')
cat('============================================================\n')
cat('GENESPACE initialization passed\n')
cat('============================================================\n')

cat('Parameters saved to:\n', parameterFile, '\n', sep = '')


# ============================================================
# save session information before long analysis
# ============================================================

writeLines(capture.output(sessionInfo()), sessionFile)


# ============================================================
# run complete GENESPACE pipeline
# ============================================================

cat('\n')
cat('============================================================\n')
cat('STARTING GENESPACE PIPELINE\n')
cat('============================================================\n\n')


out <- run_genespace(gsParam = gpar)


# ============================================================
# save final GENESPACE object
# ============================================================

saveRDS(out, resultFile)


# ============================================================
# update session information
# ============================================================

writeLines(capture.output(sessionInfo()), sessionFile)


# ============================================================
# final report
# ============================================================

cat('\n')
cat('============================================================\n')
cat('GENESPACE RUN COMPLETED\n')
cat('============================================================\n')

cat('Final GENESPACE object:\n', resultFile, '\n', sep = '')
cat('Session information:\n', sessionFile, '\n', sep = '')

cat('\n')
