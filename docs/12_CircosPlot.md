## 12) Circos plot
Circos plots provide a useful way to visualize chromosome lengths, gene density, GC content, repeat content, etc. 

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
# ============================================================
# Gloydius ussuriensis genome circos plot
#
# Tracks:
#   chromosome ideogram
#   GC %
#   Repeat %
#   gene density
#
# Input windows:
#   100 kb non-overlapping windows
#
# Customizations:
#   - white overall figure background
#   - macrochromosomes alternate dark/light gray
#   - microchromosomes alternate dark/light green
#   - chrZ = steelblue3
#   - chrW = goldenrod2
#   - chromosome ideograms have true rounded ends
#   - chromosome bars are inset to prevent overlap
#   - microchromosomes have additional spacing
#   - GC/Repeat backgrounds alternate light/darker gray
#   - GC/Repeat values above genome-wide median are red
#   - gene density shown as red heatmap
#   - chromosome labels omitted
#   - larger axis text
#   - larger chromosome ticks
#   - microchromosomes show upper-bound axis label only
#
# ============================================================


# ------------------------------------------------------------
# clean working environment
# ------------------------------------------------------------

rm(list = ls(all.names = T))
gc()


# ------------------------------------------------------------
# load package
# ------------------------------------------------------------

library(circlize)


# ------------------------------------------------------------
# set paths
# ------------------------------------------------------------

root <- paste0('/home/yshin/Gloydius_ussuriensis_genome_assembly/R/', 'Rdata/circos')
track_file <- file.path(root, 'G_ussuriensis.100kb.circos_tracks.tsv')
length_file <- file.path(root, 'G_ussuriensis.chromosome_lengths.tsv')
class_file <- file.path(root, 'G_ussuriensis.chromosome_classes.tsv')

outdir <- '/home/yshin/Gloydius_ussuriensis_genome_assembly/R/Rplots'
out_png <- file.path(outdir, 'G_ussuriensis_circos_v1_100kb_custom.png')


# ------------------------------------------------------------
# check input files
# ------------------------------------------------------------

if (!file.exists(track_file)) {
  stop(paste('Cannot find track file:', track_file))
}

if (!file.exists(length_file)) {
  stop(paste('Cannot find chromosome length file:', length_file))
}

if (!file.exists(class_file)) {
  stop(paste('Cannot find chromosome class file:', class_file))
  
}


# ------------------------------------------------------------
# read data
# ------------------------------------------------------------

dat <- read.delim(track_file, header = T, stringsAsFactors = F, check.names = F)
chrom <- read.delim(length_file, header = T, stringsAsFactors = F, check.names = F)
chr_class <- read.delim(class_file, header = T, stringsAsFactors = F, check.names = F)


# ------------------------------------------------------------
# define chromosome plotting order
# ------------------------------------------------------------

chr_order <- c(
  'G_ussuri_chr1',
  'G_ussuri_chr2',
  'G_ussuri_chr3',
  'G_ussuri_chr4',
  'G_ussuri_chr5',
  'G_ussuri_chr6',
  'G_ussuri_chr7',
  'G_ussuri_chrZ',
  'G_ussuri_chrW',
  'G_ussuri_chr9',
  'G_ussuri_chr10',
  'G_ussuri_chr11',
  'G_ussuri_chr12',
  'G_ussuri_chr13',
  'G_ussuri_chr14',
  'G_ussuri_chr15',
  'G_ussuri_chr16',
  'G_ussuri_chr17',
  'G_ussuri_chr18'
)


# ------------------------------------------------------------
# confirm chromosomes exist
# ------------------------------------------------------------

missing_length_chr <- setdiff(chr_order, chrom$chr)
missing_class_chr <- setdiff(chr_order, chr_class$chr)

if (length(missing_length_chr) > 0) {
  stop(paste('Missing chromosomes from length table:',
             paste(missing_length_chr, collapse = ', ')))
}

if (length(missing_class_chr) > 0) {
  stop(paste('Missing chromosomes from class table:',
             paste(missing_class_chr, collapse = ', ')))
}


# ------------------------------------------------------------
# arrange chromosome metadata
# ------------------------------------------------------------

chrom <- chrom[match(chr_order, chrom$chr), ]
chr_class <- chr_class[match(chr_order, chr_class$chr), ]

if (!all(chrom$chr == chr_class$chr)) {
  stop('Chromosome length and chromosome class tables do not match.')
}

chrom$display <- chr_class$display
chrom$class <- chr_class$class


# ------------------------------------------------------------
# order track data
# ------------------------------------------------------------

dat <- dat[dat$chr %in% chr_order, ]
dat$chr <- factor(dat$chr, levels = chr_order)
dat <- dat[order(dat$chr, dat$start), ]
dat$mid <- (dat$start + dat$end) / 2


# ------------------------------------------------------------
# basic sanity checks
# ------------------------------------------------------------

cat('\n')
cat('============================================================\n')
cat('Circos input summary\n')
cat('============================================================\n\n')

cat('Chromosomes:', nrow(chrom), '\n')
cat('Windows:', nrow(dat), '\n')
cat('Genes represented:', sum(dat$gene_count, na.rm = T), '\n')
cat('Mean GC %:', mean(dat$GC_pct, na.rm = T), '\n')
cat('Mean repeat %:', mean(dat$Repeat_pct, na.rm = T), '\n\n')


# ------------------------------------------------------------
# chromosome colors
# ------------------------------------------------------------

macro_dark <- '#8f8f8f'
macro_light <- '#bdbdbd'

micro_dark <- '#4c9a67'
micro_light <- '#7fbc96'

z_col <- 'steelblue3'
w_col <- 'goldenrod2'

chrom$color <- NA_character_
chrom$shade <- NA_character_


# ------------------------------------------------------------
# macrochromosomes
# ------------------------------------------------------------

macro_idx <- which(chrom$class == 'Macrochromosome')
macro_shades <- rep(c('dark', 'light'), length.out = length(macro_idx))

chrom$shade[macro_idx] <- macro_shades
chrom$color[macro_idx] <- ifelse(macro_shades == 'dark', macro_dark, macro_light)


# ------------------------------------------------------------
# microchromosomes
# ------------------------------------------------------------

micro_idx <- which(chrom$class == 'Microchromosome')
micro_shades <- rep(c('dark', 'light'), length.out = length(micro_idx))

chrom$shade[micro_idx] <- micro_shades
chrom$color[micro_idx] <- ifelse(micro_shades == 'dark', micro_dark, micro_light)


# ------------------------------------------------------------
# Z chromosome
# ------------------------------------------------------------

z_idx <- which(chrom$chr == 'G_ussuri_chrZ')

chrom$color[z_idx] <- z_col
chrom$shade[z_idx] <- 'dark'


# ------------------------------------------------------------
# W chromosome
# ------------------------------------------------------------

w_idx <- which(chrom$chr == 'G_ussuri_chrW')

chrom$color[w_idx] <- w_col
chrom$shade[w_idx] <- 'light'


# ------------------------------------------------------------
# alternating GC/Repeat background colors
# ------------------------------------------------------------

track_bg_dark <- '#e1e1e1'
track_bg_light <- '#f0f0f0'

chrom$track_bg <- ifelse(chrom$shade == 'dark', track_bg_dark, track_bg_light)


# ------------------------------------------------------------
# other plot colors
# ------------------------------------------------------------

# pure white overall canvas
figure_bg <- 'white'

line_col <- '#222222'

highlight_col <- '#ef4b35'

gene_low <- '#fff5f0'
gene_mid <- '#fcae91'
gene_high <- '#e34a33'


# ------------------------------------------------------------
# genome-wide median thresholds
# ------------------------------------------------------------

gc_median <- median(dat$GC_pct, na.rm = T)
repeat_median <- median(dat$Repeat_pct, na.rm = T)


# ------------------------------------------------------------
# robust plotting limits
# ------------------------------------------------------------

gc_ylim <- quantile(dat$GC_pct, probs = c(0.01, 0.99), na.rm = T)
repeat_ylim <- quantile(dat$Repeat_pct, probs = c(0.01, 0.99), na.rm = T)

dat$GC_plot <- pmin(pmax(dat$GC_pct, gc_ylim[1]), gc_ylim[2])
dat$Repeat_plot <- pmin(pmax(dat$Repeat_pct, repeat_ylim[1]), repeat_ylim[2])


# ------------------------------------------------------------
# gene density color scale
# ------------------------------------------------------------

gene_cap <- quantile(dat$gene_count, probs = 0.99, na.rm = T)

if (
  !is.finite(gene_cap) |
  gene_cap <= 0
) {
  
  gene_cap <- max(dat$gene_count, na.rm = T)
  
}

if (gene_cap <= 0) {
  
  gene_cap <- 1
  
}

gene_col_fun <- colorRamp2(c(0, gene_cap / 2, gene_cap), c(gene_low, gene_mid, gene_high))


# ------------------------------------------------------------
# chromosome tick spacing
# ------------------------------------------------------------

get_tick_step <- function(length_bp) {
  
  if (length_bp >= 250e6) {
    return(50e6)
  }
  
  if (length_bp >= 150e6) {
    return(25e6)
  }
  
  if (length_bp >= 70e6) {
    return(20e6)
  }
  
  if (length_bp >= 30e6) {
    return(10e6)
  }
  
  return(5e6)
  
}


# ------------------------------------------------------------
# tick positions
#
# microchromosomes:
#   start and end ticks only
#
# larger chromosomes:
#   regular Mb intervals
# ------------------------------------------------------------

get_tick_positions <- function(length_bp, chr_class) {
  
  if (chr_class == 'Microchromosome') {
    
    return(c(0, length_bp))
    
  }
  
  return(seq(0, length_bp, by = get_tick_step(length_bp)))
  
}


# ------------------------------------------------------------
# tick labels
#
# microchromosomes:
#   lower-bound tick remains
#   lower-bound label is suppressed
#   upper chromosome length is shown
#
# larger chromosomes:
#   all major labels are shown
# ------------------------------------------------------------

get_tick_labels <- function(
    tick_at,
    length_bp,
    chr_class
) {
  
  if (chr_class == 'Microchromosome') {
    
    return(c('', format(round(length_bp / 1e6, 1), trim = T, scientific = F)))
  }
  
  return(format(round(tick_at / 1e6, 1), trim = T, scientific = F))
}


# ------------------------------------------------------------
# axis text size
# ------------------------------------------------------------

get_tick_cex <- function(
    chr_class,
    length_bp
) {
  
  if (chr_class == 'Microchromosome') {
    return(0.50)
  }
  
  if (length_bp < 70e6) {
    return(0.54)
  }
  
  if (length_bp < 150e6) {
    return(0.58)
  }
  
  return(0.62)
  
}


# ------------------------------------------------------------
# chromosome inset
# ------------------------------------------------------------

get_cap_fraction <- function(length_bp) {
  
  if (length_bp < 30e6) {
    return(0.085)
  }
  
  if (length_bp < 70e6) {
    return(0.055)
  }
  
  if (length_bp < 150e6) {
    return(0.032)
  }
  
  return(0.025)
  
}


# ------------------------------------------------------------
# define sector gaps
# ------------------------------------------------------------

gap_after <- rep(2.2, length(chr_order))
names(gap_after) <- chr_order


# ------------------------------------------------------------
# give microchromosomes more space
# ------------------------------------------------------------

micro_chr <- c(
  'G_ussuri_chr9',
  'G_ussuri_chr10',
  'G_ussuri_chr11',
  'G_ussuri_chr12',
  'G_ussuri_chr13',
  'G_ussuri_chr14',
  'G_ussuri_chr15',
  'G_ussuri_chr16',
  'G_ussuri_chr17',
  'G_ussuri_chr18'
)

gap_after[micro_chr] <- 2.8


# ------------------------------------------------------------
# larger gaps between chromosome classes
# ------------------------------------------------------------

gap_after['G_ussuri_chr7'] <- 4.5
gap_after['G_ussuri_chrW'] <- 5.0


# ------------------------------------------------------------
# final sector gap
# ------------------------------------------------------------

gap_after['G_ussuri_chr18'] <- 6.0


# ------------------------------------------------------------
# main Circos drawing function
# ------------------------------------------------------------

draw_circos <- function() {
  
  circos.clear()
  
  
  # ----------------------------------------------------------
  # global Circos parameters
  # ----------------------------------------------------------
  
  circos.par(start.degree = 90,
             gap.after = gap_after, track.margin = c(0.004, 0.004),
             cell.padding = c(0, 0, 0, 0),
             points.overflow.warning = F,
             canvas.xlim = c(-1.22, 1.22),
             canvas.ylim = c(-1.22, 1.22))
  
  
  # ----------------------------------------------------------
  # initialize chromosomes
  # ----------------------------------------------------------
  
  circos.initialize(factors = chrom$chr, xlim = cbind(rep(0, nrow(chrom)), chrom$length))
  
  
  # ==========================================================
  # OUTER CHROMOSOME TRACK
  #
  # chromosome ideogram + Mb axis
  # ==========================================================
  
  circos.trackPlotRegion(ylim = c(0, 1),
                         track.height = 0.068,
                         bg.border = NA, 
                         panel.fun = function(x, y) {
                           chr <- CELL_META$sector.index
                           idx <- match(chr, chrom$chr)
                           chr_length <- chrom$length[idx]
                           chr_color <- chrom$color[idx]
                           chr_type <- chrom$class[idx]
                           tick_cex <- get_tick_cex(chr_type, chr_length)
      
      
      # ------------------------------------------------------
      # shorten chromosome bar for rounded caps
      # ------------------------------------------------------
      
      cap_fraction <- get_cap_fraction(chr_length)
      cap_pad <- chr_length * cap_fraction
      
      x1 <- CELL_META$xlim[1] + cap_pad
      x2 <- CELL_META$xlim[2] - cap_pad
      
      
      if (x2 <= x1) {
        
        x1 <- CELL_META$xlim[1] + chr_length * 0.10
        x2 <- CELL_META$xlim[2] - chr_length * 0.10
      }
      
      
      # ------------------------------------------------------
      # rounded chromosome bar
      # ------------------------------------------------------
      
      x_arc <- seq(x1, x2, length.out = 500)
      y_arc <- rep(0.50, length(x_arc))
      chrom_lwd <- 15
      
      circos.lines(x_arc, y_arc, col = chr_color, lwd = chrom_lwd)
      
      
      # ------------------------------------------------------
      # Mb axis
      # ------------------------------------------------------
      
      tick_at <- get_tick_positions(chr_length, chr_type)
      tick_labels <- get_tick_labels(tick_at, chr_length, chr_type)
      
      circos.axis(h = 'top', major.at = tick_at, labels = tick_labels,
                  labels.cex = tick_cex, labels.facing = 'clockwise',
                  labels.niceFacing = T, major.tick.length = 0.10,
                  minor.ticks = 4, lwd = 0.80)
    }
  )
  
  
  # ==========================================================
  # GC TRACK
  # ==========================================================
  
  circos.trackPlotRegion(ylim = gc_ylim, track.height = 0.105,
                         bg.col = NA, bg.border = NA,
                         panel.fun = function(x, y) {
                           chr <- CELL_META$sector.index
                           idx <- match(chr, chrom$chr)
                           sector_bg <- chrom$track_bg[idx]
                           circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                                       CELL_META$xlim[2], CELL_META$ylim[2], 
                                       col = sector_bg, border = NA)
      
      tmp <- dat[dat$chr == chr, ]
      tmp <- tmp[order(tmp$start), ]
      
      circos.lines(tmp$mid, tmp$GC_plot, col = line_col, lwd = 0.40)
      gc_highlight <- ifelse(tmp$GC_pct > gc_median, tmp$GC_plot, NA)
      
      circos.lines(tmp$mid, gc_highlight, col = highlight_col, lwd = 0.65)
    }
  )
  
  
  # ==========================================================
  # REPEAT TRACK
  # ==========================================================
  
  circos.trackPlotRegion(ylim = repeat_ylim, track.height = 0.105,
                         bg.col = NA, bg.border = NA,
                         panel.fun = function(x, y) {
                           chr <- CELL_META$sector.index
                           idx <- match(chr, chrom$chr)
                           sector_bg <- chrom$track_bg[idx]
                           
                           circos.rect(CELL_META$xlim[1], CELL_META$ylim[1],
                                       CELL_META$xlim[2], CELL_META$ylim[2],
                                       col = sector_bg, border = NA)
                           
                           tmp <- dat[dat$chr == chr, ]
                           tmp <- tmp[order(tmp$start), ]
                           
                           circos.lines(tmp$mid, tmp$Repeat_plot, col = line_col, lwd = 0.40)
                           repeat_highlight <- ifelse(tmp$Repeat_pct > repeat_median, tmp$Repeat_plot, NA)
                           circos.lines(tmp$mid, repeat_highlight, col = highlight_col, lwd = 0.65)
    }
  )
  
  
  # ==========================================================
  # GENE DENSITY TRACK
  # ==========================================================
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.058, 
                         bg.col = gene_low, bg.border = NA,
                         panel.fun = function(x, y) {
                           
                           chr <- CELL_META$sector.index
                           tmp <- dat[dat$chr == chr, ]
                           tmp <- tmp[order(tmp$start), ]
                           
                           gene_values <- pmin(tmp$gene_count, gene_cap)
                           gene_colors <- gene_col_fun(gene_values)
                           
                           circos.rect(tmp$start, 0, tmp$end, 1, col = gene_colors, border = NA)
    }
  )
  
  
  # ----------------------------------------------------------
  # legend
  # ----------------------------------------------------------
  
  legend(x = 'bottom', legend = c('Macrochromosome', 'Microchromosome', 'Chr Z', 'Chr W'),
         col = c(macro_dark, micro_dark, z_col, w_col),
         pch = 15, pt.cex = 1.6, cex = 0.70, horiz = T,
         bty = 'n', inset = c(0, -0.075), xpd = NA)
  
  circos.clear()
  
}


# ============================================================
# save PNG
# ============================================================

png(filename = out_png, width = 3600, height = 3600, res = 360, bg = figure_bg)
par(mar = c(5.5, 2, 2, 2), bg = figure_bg, lend = 'round', ljoin = 'round')

draw_circos()
dev.off()


# ============================================================
# final report
# ============================================================

cat('\n')
cat('============================================================\n')
cat('Circos plot completed\n')
cat('============================================================\n\n')

cat('Genome-wide median GC:', round(gc_median, 3), '%\n')
cat('Genome-wide median Repeat:', round(repeat_median, 3), '%\n')
cat('Gene density 99th percentile:', round(gene_cap, 3), 'genes / 100 kb\n')
cat('\n')

cat('PNG:\n', out_png, '\n\n', sep = '')
```