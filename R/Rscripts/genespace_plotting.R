# ============================================================
# Curate GENESPACE macrosynteny figure
#
# Reference:
#   Gloydius ussuriensis
#
# Display:
#   Autosomes + Z chromosome
#   W chromosome EXCLUDED
#
# Input files:
#   GENESPACE_results_11snake.rds
#   Gloydius_ussuriensis_phasedBlks.csv
#   combBed.txt
#
# This script DOES NOT rerun GENESPACE.
# It only replots the completed analysis.
# ============================================================


# ------------------------------------------------------------
# clean workspace
# ------------------------------------------------------------

rm(list = ls(all.names = T))
gc()


# ------------------------------------------------------------
# load packages
# ------------------------------------------------------------

library(GENESPACE)
library(data.table)
library(ggplot2)


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------

gs_dir <- '/home/yshin/Gloydius_ussuriensis_genome_assembly/R/Rdata/synteny'
out_dir <- 'Rplots/'


# ------------------------------------------------------------
# input files
# ------------------------------------------------------------

gs_rds <- file.path(
  gs_dir,
  'GENESPACE_results_11snake.rds'
)

comb_bed <- file.path(
  gs_dir,
  'combBed.txt'
)

phased_blocks <- file.path(
  gs_dir,
  'Gloydius_ussuriensis_phasedBlks.csv'
)


# ------------------------------------------------------------
# check inputs
# ------------------------------------------------------------

input_files <- c(
  gs_rds,
  comb_bed,
  phased_blocks
)

if (!all(file.exists(input_files))) {
  
  stop(
    paste0(
      'Missing input file(s):\n',
      paste(
        input_files[!file.exists(input_files)],
        collapse = '\n'
      )
    )
  )
  
}

cat('\nAll input files found.\n')


# ------------------------------------------------------------
# load completed GENESPACE object
# ------------------------------------------------------------

gs <- readRDS(gs_rds)


# ------------------------------------------------------------
# redirect GENESPACE results path to local files
#
# plot_riparian() reads:
#
#   combBed.txt
#
# from gs$paths$results
# ------------------------------------------------------------

gs$paths$results <- gs_dir


# ------------------------------------------------------------
# load phased syntenic blocks
# ------------------------------------------------------------

blks <- fread(phased_blocks)

cat('\n========================================\n')
cat('Original phased blocks\n')
cat('========================================\n')

cat('\nNumber of phased blocks:\n')
print(nrow(blks))

cat('\nReference genome(s):\n')
print(unique(blks$refGenome))


# ------------------------------------------------------------
# check expected columns
# ------------------------------------------------------------

required_cols <- c(
  'genome1',
  'genome2',
  'chr1',
  'chr2',
  'refChr',
  'refGenome'
)

missing_cols <- setdiff(
  required_cols,
  colnames(blks)
)

if (length(missing_cols) > 0) {
  
  stop(
    paste(
      'Missing columns in phased block file:',
      paste(
        missing_cols,
        collapse = ', '
      )
    )
  )
  
}


# ------------------------------------------------------------
# check how many W-associated blocks exist
# ------------------------------------------------------------

n_w_blocks <- sum(
  blks$refChr == 'G_ussuri_chrW',
  na.rm = T
)

cat('\nW-associated phased blocks before filtering:\n')
print(n_w_blocks)


# ------------------------------------------------------------
# exclude G. ussuriensis W chromosome
#
# refChr identifies the G. ussuriensis reference chromosome
# assigned to each phased syntenic block.
#
# Therefore, filtering refChr removes:
#
#   1. the G. ussuriensis W chromosome
#   2. all syntenic ribbons corresponding to W
#      throughout the other genomes
#
# This is different from simply hiding the W label.
# ------------------------------------------------------------

blks_noW <- blks[
  is.na(refChr) |
    refChr != 'G_ussuri_chrW'
]


# ------------------------------------------------------------
# verify W removal
# ------------------------------------------------------------

if ('G_ussuri_chrW' %in% blks_noW$refChr) {
  
  stop(
    'W chromosome is still present after filtering.'
  )
  
}

cat('\n========================================\n')
cat('W chromosome exclusion\n')
cat('========================================\n')

cat('\nBlocks before filtering:\n')
print(nrow(blks))

cat('\nBlocks after filtering:\n')
print(nrow(blks_noW))

cat('\nBlocks removed:\n')
print(
  nrow(blks) -
    nrow(blks_noW)
)

cat('\nW chromosome successfully excluded.\n')


# ------------------------------------------------------------
# create separate riparian directory for no-W plotting
#
# The original phased block file remains untouched.
# ------------------------------------------------------------

noW_rip_dir <- file.path(
  gs_dir,
  'riparian_noW'
)

dir.create(
  noW_rip_dir,
  recursive = T,
  showWarnings = F
)


# ------------------------------------------------------------
# GENESPACE expects this exact filename
# ------------------------------------------------------------

noW_phased_blocks <- file.path(
  noW_rip_dir,
  'Gloydius_ussuriensis_phasedBlks.csv'
)


# ------------------------------------------------------------
# write filtered phased-block file
# ------------------------------------------------------------

fwrite(
  blks_noW,
  noW_phased_blocks
)

cat('\nFiltered phased-block file written to:\n')
cat(noW_phased_blocks, '\n')


# ------------------------------------------------------------
# redirect GENESPACE riparian path
#
# forceRecalcBlocks = F below tells GENESPACE to use
# this filtered file rather than recalculate blocks.
# ------------------------------------------------------------

gs$paths$riparian <- noW_rip_dir


# ------------------------------------------------------------
# confirm redirected paths
# ------------------------------------------------------------

cat('\nGENESPACE results path:\n')
print(gs$paths$results)

cat('\nGENESPACE riparian path:\n')
print(gs$paths$riparian)

cat('\ncombBed exists:\n')
print(
  file.exists(
    file.path(
      gs$paths$results,
      'combBed.txt'
    )
  )
)

cat('\nno-W phased blocks exist:\n')
print(
  file.exists(
    file.path(
      gs$paths$riparian,
      'Gloydius_ussuriensis_phasedBlks.csv'
    )
  )
)


# ------------------------------------------------------------
# species order
#
# first genome = bottom of plot
# last genome  = top of plot
#
# Therefore, from top -> bottom:
#
# G. ussuriensis
# G. shedaoensis
# C. adamanteus
# Bothrops
# Vipera
# Cerastes
# Naja
# Elaphe
# Candoia
# Xenopeltis
# Argyrophis
# ------------------------------------------------------------

genomeIDs <- c(
  'Argyrophis_diardii',
  'Xenopeltis_unicolor',
  'Candoia_aspera',
  'Elaphe_schrenckii',
  'Naja_naja',
  'Cerastes_gasperettii',
  'Vipera_berus',
  'Bothrops_insularis',
  'Crotalus_adamanteus',
  'Gloydius_shedaoensis',
  'Gloydius_ussuriensis'
)


# ------------------------------------------------------------
# make sure all requested genomes remain represented
# ------------------------------------------------------------

block_genomes <- unique(
  c(
    blks_noW$genome1,
    blks_noW$genome2
  )
)

missing_genomes <- setdiff(
  genomeIDs,
  block_genomes
)

if (length(missing_genomes) > 0) {
  
  stop(
    paste(
      'These genomes are absent from phased blocks:',
      paste(
        missing_genomes,
        collapse = ', '
      )
    )
  )
  
}


# ------------------------------------------------------------
# determine reference chromosome order
#
# W has already been removed.
#
# Expected final reference chromosomes:
#
# G_ussuri_chr1
# ...
# G_ussuri_chr18
# G_ussuri_chrZ
# ------------------------------------------------------------

ref_chrs <- unique(
  blks_noW[
    refGenome == 'Gloydius_ussuriensis',
    refChr
  ]
)

ref_chrs <- ref_chrs[
  !is.na(ref_chrs)
]


cat('\n========================================\n')
cat('Reference chromosomes retained\n')
cat('========================================\n\n')

print(ref_chrs)


# ------------------------------------------------------------
# autosomes / numerically named chromosomes
# ------------------------------------------------------------

numeric_chrs <- ref_chrs[
  grepl(
    '^G_ussuri_chr[0-9]+$',
    ref_chrs
  )
]

numeric_chrs <- numeric_chrs[
  order(
    as.integer(
      sub(
        '^G_ussuri_chr',
        '',
        numeric_chrs
      )
    )
  )
]


# ------------------------------------------------------------
# Z chromosome
# ------------------------------------------------------------

z_chr <- intersect(
  'G_ussuri_chrZ',
  ref_chrs
)


# ------------------------------------------------------------
# anything unexpected
# ------------------------------------------------------------

other_chrs <- setdiff(
  ref_chrs,
  c(
    numeric_chrs,
    z_chr
  )
)


# ------------------------------------------------------------
# final reference chromosome order
#
# W deliberately excluded
# ------------------------------------------------------------

ref_chr_order <- c(
  numeric_chrs,
  sort(other_chrs),
  z_chr
)


# ------------------------------------------------------------
# final safety check
# ------------------------------------------------------------

if ('G_ussuri_chrW' %in% ref_chr_order) {
  
  stop(
    'W chromosome unexpectedly present in ref_chr_order.'
  )
  
}


cat('\nFinal reference chromosome plotting order:\n')
print(ref_chr_order)


# ------------------------------------------------------------
# chromosome label function
#
# Only G. ussuriensis receives chromosome labels.
#
# Examples:
#
# G_ussuri_chr1 -> 1
# G_ussuri_chr12 -> 12
# G_ussuri_chrZ -> Z
#
# W is excluded.
# ------------------------------------------------------------

chr_lab_fun <- function(x) {
  
  x <- sub(
    '^G_ussuri_chr',
    '',
    x
  )
  
  toupper(x)
  
}


# ------------------------------------------------------------
# plotting parameters
# ------------------------------------------------------------

# Minimum displayed chromosome/scaffold span.
#
# This DOES NOT alter the synteny analysis.
#
# It only filters short sequence fragments from the
# visualization.
# ------------------------------------------------------------

min_chr_len <- 5e6


# ------------------------------------------------------------
# final figure dimensions
# ------------------------------------------------------------

plot_width <- 14
plot_height <- 8.5


# ------------------------------------------------------------
# plot theme
# ------------------------------------------------------------

rip_theme <- theme(
  
  panel.background = element_rect(
    fill = 'white',
    color = NA
  ),
  
  plot.background = element_rect(
    fill = 'white',
    color = NA
  ),
  
  panel.grid = element_blank(),
  
  panel.border = element_blank(),
  
  axis.ticks = element_blank(),
  
  axis.text.x = element_blank(),
  
  axis.text.y = element_text(
    size = 10,
    face = 'italic',
    color = 'black'
  ),
  
  axis.title.x = element_text(
    size = 10,
    color = 'black',
    margin = margin(
      t = 8
    )
  ),
  
  plot.margin = margin(
    t = 10,
    r = 15,
    b = 10,
    l = 10
  )
  
)


# ------------------------------------------------------------
# temporary graphics device
#
# GENESPACE uses device dimensions when calculating
# chromosome polygon height.
# ------------------------------------------------------------

tmp_pdf <- tempfile(
  fileext = '.pdf'
)

pdf(
  file = tmp_pdf,
  width = plot_width,
  height = plot_height
)


# ------------------------------------------------------------
# generate curated riparian plot
#
# IMPORTANT:
#
# forceRecalcBlocks = F
#
# This ensures GENESPACE reads our filtered no-W phased
# block file instead of reconstructing the original blocks.
# ------------------------------------------------------------

rip <- plot_riparian(
  
  gsParam = gs,
  
  genomeIDs = genomeIDs,
  
  refGenome = 'Gloydius_ussuriensis',
  
  useOrder = F,
  
  forceRecalcBlocks = F,
  
  labelTheseGenomes = 'Gloydius_ussuriensis',
  
  customRefChrOrder = ref_chr_order,
  
  minChrLen2plot = min_chr_len,
  
  reorderBySynteny = T,
  
  syntenyWeight = 0.5,
  
  gapProp = 0.004,
  
  scaleGapSize = 0.20,
  
  scaleBraidGap = 0.15,
  
  braidAlpha = 0.55,
  
  chrFill = 'white',
  
  chrBorderCol = 'black',
  
  chrBorderLwd = 0.35,
  
  chrLabFontSize = 7,
  
  chrExpand = 0.55,
  
  chrLabFun = chr_lab_fun,
  
  xlabel = 'Chromosomes scaled by physical position',
  
  addThemes = rip_theme,
  
  verbose = T
  
)


dev.off()

unlink(tmp_pdf)


# ------------------------------------------------------------
# extract ggplot object
# ------------------------------------------------------------

p <- rip$plotData$ggplotObj


# ------------------------------------------------------------
# extract chromosome plotting data
# ------------------------------------------------------------

chr_plot_data <- as.data.table(
  rip$plotData$sourceData$chromosomes
)


# ------------------------------------------------------------
# verify W is absent from actual plotted chromosomes
# ------------------------------------------------------------

w_in_plot <- chr_plot_data[
  genome == 'Gloydius_ussuriensis' &
    chr == 'G_ussuri_chrW'
]

if (nrow(w_in_plot) > 0) {
  
  stop(
    'W chromosome is unexpectedly present in final plot data.'
  )
  
}

cat('\nW chromosome absent from final plotting data: PASS\n')


# ------------------------------------------------------------
# verify no W-assigned ribbons survived
# ------------------------------------------------------------

if (
  'refChr' %in% colnames(rip$blks) &&
  any(
    rip$blks$refChr == 'G_ussuri_chrW',
    na.rm = T
  )
) {
  
  stop(
    'W-associated syntenic blocks unexpectedly remain.'
  )
  
}

cat('W-associated syntenic ribbons absent: PASS\n')


# ------------------------------------------------------------
# clean species labels
#
# Gloydius_ussuriensis
#
# becomes
#
# Gloydius ussuriensis
#
# and is displayed in italics.
# ------------------------------------------------------------

species_plot_data <- chr_plot_data[
  !duplicated(genome)
]

setorder(
  species_plot_data,
  y1
)

species_breaks <- (
  species_plot_data$y1 +
    species_plot_data$y2
) / 2

species_labels <- gsub(
  '_',
  ' ',
  species_plot_data$genome
)


# ------------------------------------------------------------
# replace y-axis labels
# ------------------------------------------------------------

p <- p +
  
  scale_y_continuous(
    
    breaks = species_breaks,
    
    labels = species_labels,
    
    expand = c(
      0.01,
      0.01
    ),
    
    name = NULL
    
  ) +
  
  theme(
    
    axis.text.y = element_text(
      size = 10,
      face = 'italic',
      color = 'black'
    )
    
  )


# ------------------------------------------------------------
# output filenames
# ------------------------------------------------------------

pdf_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_noW_curated_v1.pdf'
)

png_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_noW_curated_v1.png'
)

rds_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_noW_curated_v1.rds'
)

chr_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_noW_curated_v1_chromosomes.tsv'
)


# ------------------------------------------------------------
# save PDF
# ------------------------------------------------------------

ggsave(
  
  filename = pdf_out,
  
  plot = p,
  
  width = plot_width,
  
  height = plot_height,
  
  units = 'in',
  
  device = 'pdf'
  
)


# ------------------------------------------------------------
# save high-resolution PNG
# ------------------------------------------------------------

ggsave(
  
  filename = png_out,
  
  plot = p,
  
  width = plot_width,
  
  height = plot_height,
  
  units = 'in',
  
  dpi = 300
  
)


# ------------------------------------------------------------
# save complete plotting object
# ------------------------------------------------------------

saveRDS(
  rip,
  rds_out
)


# ------------------------------------------------------------
# save chromosome/scaffold plotting table
# ------------------------------------------------------------

fwrite(
  chr_plot_data,
  chr_out,
  sep = '\t'
)


# ------------------------------------------------------------
# summarize sequences shown per genome
# ------------------------------------------------------------

plot_summary <- chr_plot_data[
  ,
  .(
    n_sequences_plotted = uniqueN(chr)
  ),
  by = genome
]

plot_summary[
  ,
  plot_order := match(
    genome,
    genomeIDs
  )
]

setorder(
  plot_summary,
  plot_order
)

plot_summary[
  ,
  plot_order := NULL
]


cat('\n========================================\n')
cat('Sequences retained in curated plot\n')
cat('========================================\n\n')

print(plot_summary)


# ------------------------------------------------------------
# report G. ussuriensis chromosomes shown
# ------------------------------------------------------------

gussuri_chrs_plotted <- chr_plot_data[
  genome == 'Gloydius_ussuriensis',
  unique(chr)
]

cat('\n========================================\n')
cat('G. ussuriensis chromosomes plotted\n')
cat('========================================\n\n')

print(gussuri_chrs_plotted)


# ------------------------------------------------------------
# final output report
# ------------------------------------------------------------

cat('\n========================================\n')
cat('CURATED NO-W FIGURE COMPLETE\n')
cat('========================================\n')

cat('\nPDF:\n')
cat(pdf_out, '\n')

cat('\nPNG:\n')
cat(png_out, '\n')

cat('\nRDS:\n')
cat(rds_out, '\n')

cat('\nChromosome table:\n')
cat(chr_out, '\n')

cat('\nFiltered phased blocks:\n')
cat(noW_phased_blocks, '\n')
