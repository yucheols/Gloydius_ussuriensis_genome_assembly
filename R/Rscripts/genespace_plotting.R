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
#   - extra-large species labels
#   - extra-large chromosome labels
#   - chromosome widths scaled by representative gene order
#   - Naja naja microchromosomes shown as numbers only
#   - selected chromosomes inverted for cleaner visualization
#   - asterisks added by GENESPACE to inverted chromosomes removed
#   - "5000 genes" above scale bar
#   - scale description to LEFT of scale bar
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
        req_files[
          !file.exists(req_files)
        ],
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


# ============================================================
# helper functions
# ============================================================

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


# ============================================================
# genome order
#
# first = bottom
# last  = top
# ============================================================

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


if (
  !setequal(
    ref_chr_order,
    ref_chrs
  )
) {
  
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
    'chr7',
    
    # Crotalus adamanteus
    'chr3',
    'chr4',
    'chr6',
    
    # Elaphe schrenckii
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    
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
# canvas
#
# Slightly wider to accommodate long scale annotation
# on the LEFT side of scale bar.
# ------------------------------------------------------------

plot_width <- 27

plot_height <- 13.5


# ------------------------------------------------------------
# font sizes
# ------------------------------------------------------------

species_font_size <- 22

chromosome_font_size <- 14

scale_font_size <- species_font_size


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
    size = species_font_size,
    face = 'italic',
    color = 'black'
  ),
  
  # remove old x-axis title
  
  axis.title.x = element_blank(),
  
  plot.margin = margin(
    t = 38,
    r = 35,
    b = 18,
    l = 25
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
  
  useOrder = T,
  
  forceRecalcBlocks = F,
  
  labelTheseGenomes = genomeIDs,
  
  customRefChrOrder = ref_chr_order,
  
  minChrLen2plot = min_chr_genes,
  
  reorderBySynteny = T,
  
  syntenyWeight = 0.5,
  
  gapProp = 0.004,
  
  scaleGapSize = 0.20,
  
  scaleBraidGap = 0.15,
  
  braidAlpha = 0.55,
  
  chrFill = 'grey90',
  
  chrBorderCol = 'black',
  
  chrBorderLwd = 0.40,
  
  chrLabFontSize = chromosome_font_size,
  
  chrExpand = 1.25,
  
  chrLabFun = chr_lab_fun,
  
  invertTheseChrs = invertTheseChrs,
  
  # remove GENESPACE's bottom description
  
  xlabel = NULL,
  
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
# remove asterisks from manually inverted chromosomes
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
# remove asterisks from source chromosome data
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
# locate GENESPACE scale bar
# ============================================================

scale_layer_index <- NA_integer_

scale_bar_data <- NULL


for (i in seq_along(p$layers)) {
  
  layer_data <- p$layers[[i]]$data
  
  
  if (
    is.null(layer_data) ||
    !is.data.frame(layer_data) ||
    nrow(layer_data) == 0
  ) {
    
    next
    
  }
  
  
  required_scale_cols <- c(
    'line',
    'x',
    'xend',
    'y',
    'yend'
  )
  
  
  if (
    all(
      required_scale_cols %in%
      names(layer_data)
    ) &&
    all(
      c(
        'left',
        'right',
        'mid'
      ) %in%
      layer_data$line
    )
  ) {
    
    scale_layer_index <- i
    
    scale_bar_data <- as.data.table(
      layer_data
    )
    
    break
    
  }
  
}


if (is.na(scale_layer_index)) {
  
  stop(
    'Could not locate GENESPACE scale-bar layer.'
  )
  
}


# ============================================================
# scale-bar geometry
# ============================================================

scale_mid <- scale_bar_data[
  line == 'mid'
]


scale_left <- min(
  scale_mid$x,
  scale_mid$xend
)


scale_right <- max(
  scale_mid$x,
  scale_mid$xend
)


scale_mid_x <- (
  scale_left +
    scale_right
) / 2


scale_mid_y <- scale_mid$y[1]


scale_top <- max(
  c(
    scale_bar_data$y,
    scale_bar_data$yend
  ),
  na.rm = T
)


scale_bottom <- min(
  c(
    scale_bar_data$y,
    scale_bar_data$yend
  ),
  na.rm = T
)


scale_bar_height <- (
  scale_top -
    scale_bottom
)


# ------------------------------------------------------------
# overall x span
# ------------------------------------------------------------

plot_x_min <- min(
  chr_plot_data$x1,
  na.rm = T
)


plot_x_max <- max(
  chr_plot_data$x2,
  na.rm = T
)


plot_x_span <- (
  plot_x_max -
    plot_x_min
)


# ============================================================
# hide original GENESPACE scale-bar text
# ============================================================

for (i in seq_along(p$layers)) {
  
  layer_data <- p$layers[[i]]$data
  
  
  if (
    is.null(layer_data) ||
    !is.data.frame(layer_data) ||
    nrow(layer_data) == 0
  ) {
    
    next
    
  }
  
  
  if (
    'line' %in% names(layer_data) &&
    nrow(layer_data) == 1 &&
    all(
      layer_data$line == 'mid'
    ) &&
    inherits(
      p$layers[[i]]$geom,
      'GeomText'
    )
  ) {
    
    p$layers[[i]]$aes_params$alpha <- 0
    
  }
  
}


# ============================================================
# custom scale annotations
#
# Layout:
#
# Chromosomes scaled by the number of genes    5000 genes
#                                               |---------|
#
# More precisely:
#
# description = LEFT of bar
# 5000 genes  = ABOVE bar
# ============================================================


# ------------------------------------------------------------
# "5000 genes" above scale bar
# ------------------------------------------------------------

scale_value_y <- (
  scale_top +
    0.65 * scale_bar_height
)


# ------------------------------------------------------------
# scale description immediately LEFT of scale bar
#
# hjust = 1 means the RIGHT EDGE of the text terminates
# at scale_description_x.
# ------------------------------------------------------------

scale_description_x <- (
  scale_left -
    0.025 * plot_x_span
)


scale_description_y <- scale_mid_y


# ============================================================
# add custom scale annotations
# ============================================================

p <- p +
  
  # ----------------------------------------------------------
# value above scale bar
# ----------------------------------------------------------

annotate(
  
  geom = 'text',
  
  x = scale_mid_x,
  
  y = scale_value_y,
  
  label = '5000 genes',
  
  size = scale_font_size / ggplot2::.pt,
  
  hjust = 0.5,
  
  vjust = 0.5,
  
  color = 'black'
  
) +
  
  
  # ----------------------------------------------------------
# description LEFT of scale bar
# ----------------------------------------------------------

annotate(
  
  geom = 'text',
  
  x = scale_description_x,
  
  y = scale_mid_y,
  
  label = 'Chromosomes scaled by the number of genes',
  
  size = scale_font_size / ggplot2::.pt,
  
  hjust = 1,
  
  vjust = 0.5,
  
  color = 'black'
  
)


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
  
  labs(
    x = NULL
  ) +
  
  theme(
    
    axis.text.y = element_text(
      size = species_font_size,
      face = 'italic',
      color = 'black'
    ),
    
    axis.title.x = element_blank()
    
  )


# ============================================================
# output files
# ============================================================

pdf_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_scaleLeft_v9.pdf'
)


png_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_scaleLeft_v9.png'
)


rds_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_scaleLeft_v9.rds'
)


chr_out <- file.path(
  out_dir,
  'Gloydius_ussuriensis_macrosynteny_geneScaled_scaleLeft_v9_chromosomes.tsv'
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
cat('CUSTOM SCALE ANNOTATION\n')
cat('========================================\n\n')


cat(
  'Scale description:\n'
)

cat(
  'Chromosomes scaled by the number of genes\n\n'
)


cat(
  'Scale description position:\n'
)

cat(
  'LEFT of scale bar\n\n'
)


cat(
  'Scale value:\n'
)

cat(
  '5000 genes above scale bar\n'
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