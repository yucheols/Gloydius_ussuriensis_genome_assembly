# ============================================================
# Curate GENESPACE macrosynteny figure
#
# Reference:
#   Gloydius ussuriensis
#
# Display:
#   - chromosome labels for all species
#   - chromosome order = Macro -> Z -> Micro
#   - chromosome bars = light grey
#   - thicker chromosome bars
#   - EXTRA-LARGE species labels
#   - EXTRA-LARGE chromosome labels
#   - chromosome widths scaled by representative gene order
#   - Naja naja microchromosomes shown as numbers only
#   - selected chromosomes inverted for cleaner visualization
#   - asterisks added by GENESPACE to inverted chromosomes removed
#
# This script replots a completed GENESPACE run.
# It does NOT rerun GENESPACE.
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
out_dir <- 'Rplots'

dir.create(
  out_dir,
  recursive = T,
  showWarnings = F
)


# ------------------------------------------------------------
# input files
# ------------------------------------------------------------

gs_rds <- file.path(
  gs_dir,
  'GENESPACE_results_11snake_chrOnly.rds'
)

phased_blocks <- file.path(
  gs_dir,
  'Gloydius_ussuriensis_phasedBlks.csv'
)

comb_bed <- file.path(
  gs_dir,
  'combBed.txt'
)


req_files <- c(
  gs_rds,
  phased_blocks,
  comb_bed
)


if (!all(file.exists(req_files))) {
  
  stop(
    paste0(
      'Missing input file(s):\n',
      paste(
        req_files[!file.exists(req_files)],
        collapse = '\n'
      )
    )
  )
  
}


# ------------------------------------------------------------
# load GENESPACE objects
# ------------------------------------------------------------

gs <- readRDS(
  gs_rds
)

blks <- fread(
  phased_blocks
)


# redirect internal paths to local files

gs$paths$results <- gs_dir
gs$paths$riparian <- gs_dir


# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------

is_w_chr <- function(x) {
  
  grepl(
    '(^|_)chrW$|(^|_)W$',
    x,
    ignore.case = T
  )
  
}


get_chr_num <- function(x) {
  
  suppressWarnings(
    as.integer(
      sub(
        '^.*chr',
        '',
        x,
        ignore.case = T
      )
    )
  )
  
}


# ------------------------------------------------------------
# safety check: explicit W absent
# ------------------------------------------------------------

if (
  any(
    is_w_chr(
      blks$refChr
    ),
    na.rm = T
  )
) {
  
  stop(
    'Explicit W chromosome detected in phased block data.'
  )
  
}


# ------------------------------------------------------------
# genome order
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
# verify genomes present
# ------------------------------------------------------------

block_genomes <- unique(
  c(
    blks$genome1,
    blks$genome2
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


# ============================================================
# reference chromosome order
#
# Macro -> Z -> Micro
# ============================================================

ref_chrs <- unique(
  blks[
    refGenome == 'Gloydius_ussuriensis',
    refChr
  ]
)


ref_chrs <- ref_chrs[
  !is.na(ref_chrs)
]


# ------------------------------------------------------------
# macrochromosomes
# ------------------------------------------------------------

macro_chrs <- ref_chrs[
  grepl(
    '(^|_)chr[1-7]$',
    ref_chrs,
    ignore.case = T
  )
]


macro_chrs <- macro_chrs[
  order(
    get_chr_num(
      macro_chrs
    )
  )
]


# ------------------------------------------------------------
# Z chromosome
# ------------------------------------------------------------

z_chr <- ref_chrs[
  grepl(
    '(^|_)chrZ$',
    ref_chrs,
    ignore.case = T
  )
]


# ------------------------------------------------------------
# microchromosomes
# ------------------------------------------------------------

micro_chrs <- ref_chrs[
  grepl(
    '(^|_)chr(9|1[0-8])$',
    ref_chrs,
    ignore.case = T
  )
]


micro_chrs <- micro_chrs[
  order(
    get_chr_num(
      micro_chrs
    )
  )
]


# ------------------------------------------------------------
# final reference chromosome order
# ------------------------------------------------------------

ref_chr_order <- c(
  macro_chrs,
  z_chr,
  micro_chrs
)


if (!setequal(
  ref_chr_order,
  ref_chrs
)) {
  
  stop(
    'Reference chromosome classification is incomplete.'
  )
  
}


# ============================================================
# chromosome label function
# ============================================================

chr_lab_fun <- function(x) {
  
  y <- x
  
  
  y <- sub(
    '^.*_chr',
    'chr',
    y,
    ignore.case = T
  )
  
  
  y <- gsub(
    '^chromosome',
    'chr',
    y,
    ignore.case = T
  )
  
  
  y <- gsub(
    '^scaffold',
    '',
    y,
    ignore.case = T
  )
  
  
  y <- gsub(
    '^scaf',
    '',
    y,
    ignore.case = T
  )
  
  
  y <- gsub(
    '^lg',
    '',
    y,
    ignore.case = T
  )
  
  
  y <- gsub(
    '^chr',
    '',
    y,
    ignore.case = T
  )
  
  
  # Naja microchromosomes:
  # MIC_1 -> 1
  # MIC_10 -> 10
  
  y <- gsub(
    '^MIC[_\\.:-]*',
    '',
    y,
    ignore.case = T
  )
  
  
  y <- gsub(
    '^0+',
    '',
    y
  )
  
  
  ifelse(
    y == '',
    x,
    toupper(y)
  )
  
}


# ============================================================
# selected chromosomes to invert
# ============================================================

invertTheseChrs <- data.frame(
  
  genome = c(
    
    # Argyrophis diardii
    'Argyrophis_diardii',
    'Argyrophis_diardii',
    'Argyrophis_diardii',
    'Argyrophis_diardii',
    'Argyrophis_diardii',
    'Argyrophis_diardii',
    'Argyrophis_diardii',
    
    # Bothrops insularis
    'Bothrops_insularis',
    'Bothrops_insularis',
    'Bothrops_insularis',
    
    # Candoia aspera
    'Candoia_aspera',
    'Candoia_aspera',
    'Candoia_aspera',
    
    # Cerastes gasperettii
    'Cerastes_gasperettii',
    'Cerastes_gasperettii',
    
    # Crotalus adamanteus
    'Crotalus_adamanteus',
    'Crotalus_adamanteus',
    'Crotalus_adamanteus',
    
    # Elaphe schrenckii
    'Elaphe_schrenckii',
    'Elaphe_schrenckii',
    'Elaphe_schrenckii',
    'Elaphe_schrenckii',
    
    # Gloydius ussuriensis
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    'Gloydius_ussuriensis',
    
    # Naja naja
    'Naja_naja',
    'Naja_naja',
    'Naja_naja',
    
    # Vipera berus
    'Vipera_berus',
    'Vipera_berus',
    'Vipera_berus',
    'Vipera_berus',
    
    # Xenopeltis unicolor
    'Xenopeltis_unicolor',
    'Xenopeltis_unicolor',
    'Xenopeltis_unicolor'
    
  ),
  
  chr = c(
    
    # Argyrophis diardii
    'chr1',
    'chr3',
    'chr4',
    'chr5',
    'chr7',
    'chr8',
    'chr9',
    
    # Bothrops insularis
    'chr1',
    'chr2',
    'chr7',
    
    # Candoia aspera
    'chr1',
    'chr4',
    'chr8',
    
    # Cerastes gasperettii
    'chr2',
    'chr4',
    
    # Crotalus adamanteus
    'chr3',
    'chr4',
    'chr6',
    
    # Elaphe schrenckii
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    
    # Gloydius ussuriensis
    'chr2',
    'chr3',
    'chr5',
    'chr9',
    'chr13',
    'chr15',
    'chr16',
    'chr17',
    
    # Naja naja
    'chr1',
    'chr4',
    'chr7',
    
    # Vipera berus
    'chr2',
    'chr5',
    'chr6',
    'chrZ',
    
    # Xenopeltis unicolor
    'chr1',
    'chr3',
    'chr4'
    
  ),
  
  stringsAsFactors = F
  
)

# ============================================================
# plotting parameters
# ============================================================

min_chr_genes <- 1


# ------------------------------------------------------------
# larger canvas for larger labels
# ------------------------------------------------------------

plot_width <- 24
plot_height <- 13.5


# ------------------------------------------------------------
# label sizes
#
# Increase these two values further if desired.
# ------------------------------------------------------------

species_font_size <- 22
chromosome_font_size <- 14
axis_title_font_size <- 17


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
  
  # extra-large species labels
  
  axis.text.y = element_text(
    size = species_font_size,
    face = 'italic',
    color = 'black'
  ),
  
  axis.title.x = element_text(
    size = axis_title_font_size,
    color = 'black',
    margin = margin(
      t = 14
    )
  ),
  
  plot.margin = margin(
    t = 18,
    r = 28,
    b = 18,
    l = 22
  )
  
)


# ============================================================
# temporary graphics device
# ============================================================

tmp_pdf <- tempfile(
  fileext = '.pdf'
)


pdf(
  file = tmp_pdf,
  width = plot_width,
  height = plot_height
)


# ============================================================
# generate riparian plot
# ============================================================

rip <- plot_riparian(
  
  gsParam = gs,
  
  genomeIDs = genomeIDs,
  
  refGenome = 'Gloydius_ussuriensis',
  
  # scale chromosome widths by representative gene order
  
  useOrder = T,
  
  forceRecalcBlocks = F,
  
  # chromosome labels for all species
  
  labelTheseGenomes = genomeIDs,
  
  # Macro -> Z -> Micro
  
  customRefChrOrder = ref_chr_order,
  
  # useOrder = T:
  # minimum length is in gene-order units
  
  minChrLen2plot = min_chr_genes,
  
  reorderBySynteny = T,
  
  syntenyWeight = 0.5,
  
  gapProp = 0.004,
  
  scaleGapSize = 0.20,
  
  scaleBraidGap = 0.15,
  
  braidAlpha = 0.55,
  
  # chromosome appearance
  
  chrFill = 'grey90',
  
  chrBorderCol = 'black',
  
  chrBorderLwd = 0.40,
  
  # EXTRA-LARGE chromosome labels
  
  chrLabFontSize = chromosome_font_size,
  
  # thick chromosome bars
  
  chrExpand = 1.25,
  
  chrLabFun = chr_lab_fun,
  
  # flip selected chromosomes
  
  invertTheseChrs = invertTheseChrs,
  
  xlabel = 'Chromosomes scaled by representative gene order',
  
  addThemes = rip_theme,
  
  verbose = T
  
)


dev.off()


unlink(
  tmp_pdf
)


# ============================================================
# extract ggplot object
# ============================================================

p <- rip$plotData$ggplotObj


# ============================================================
# remove asterisks added by GENESPACE to inverted chromosomes
# ============================================================

for (i in seq_along(p$layers)) {
  
  layer_data <- p$layers[[i]]$data
  
  
  if (
    !is.null(layer_data) &&
    is.data.frame(layer_data) &&
    nrow(layer_data) > 0
  ) {
    
    char_cols <- names(layer_data)[
      vapply(
        layer_data,
        is.character,
        logical(1)
      )
    ]
    
    
    for (j in char_cols) {
      
      layer_data[[j]] <- sub(
        '\\*$',
        '',
        layer_data[[j]]
      )
      
    }
    
    
    p$layers[[i]]$data <- layer_data
    
  }
  
}


# ------------------------------------------------------------
# remove asterisks from source chromosome labels
# ------------------------------------------------------------

if (
  'chrLab' %in%
  colnames(
    rip$plotData$sourceData$chromosomes
  )
) {
  
  rip$plotData$sourceData$chromosomes$chrLab <- sub(
    '\\*$',
    '',
    rip$plotData$sourceData$chromosomes$chrLab
  )
  
}


# ============================================================
# chromosome plotting data
# ============================================================

chr_plot_data <- as.data.table(
  rip$plotData$sourceData$chromosomes
)


# ------------------------------------------------------------
# verify W absent
# ------------------------------------------------------------

if (
  any(
    is_w_chr(
      chr_plot_data$chr
    ),
    na.rm = T
  )
) {
  
  stop(
    'Explicit W chromosome still present in final plot.'
  )
  
}


# ============================================================
# species labels
# ============================================================

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
# replace y-axis species labels
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
      size = species_font_size,
      face = 'italic',
      color = 'black'
    ),
    
    axis.title.x = element_text(
      size = axis_title_font_size,
      color = 'black'
    )
    
  )


# ============================================================
# output files
# ============================================================

pdf_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_macroZmicro_XLlabels_v7.pdf'
)


png_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_macroZmicro_XLlabels_v7.png'
)


rds_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_macroZmicro_XLlabels_v7.rds'
)


chr_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_macroZmicro_XLlabels_v7_chromosomes.tsv'
)


# ============================================================
# save outputs
# ============================================================

ggsave(
  
  filename = pdf_out,
  
  plot = p,
  
  width = plot_width,
  
  height = plot_height,
  
  units = 'in',
  
  device = 'pdf'
  
)


ggsave(
  
  filename = png_out,
  
  plot = p,
  
  width = plot_width,
  
  height = plot_height,
  
  units = 'in',
  
  dpi = 300
  
)


saveRDS(
  rip,
  rds_out
)


fwrite(
  chr_plot_data,
  chr_out,
  sep = '\t'
)


# ============================================================
# summary
# ============================================================

plot_summary <- chr_plot_data[
  ,
  .(
    n_chromosomes_plotted = uniqueN(chr)
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
cat('Chromosomes retained in curated plot\n')
cat('========================================\n\n')

print(
  plot_summary
)


cat('\n========================================\n')
cat('CURATED MACROSYNTENY FIGURE COMPLETE\n')
cat('========================================\n')


cat('\nPDF:\n')
cat(
  pdf_out,
  '\n'
)


cat('\nPNG:\n')
cat(
  png_out,
  '\n'
)


cat('\nRDS:\n')
cat(
  rds_out,
  '\n'
)


cat('\nChromosome table:\n')
cat(
  chr_out,
  '\n'
)


cat('\nInverted chromosomes:\n')

print(
  invertTheseChrs
)