# ============================================================
#
# circos plot
#
# chromosome ideograms, GC %, repeat %, gene density,
# and genomic positions of eight major toxin gene families
#
# ============================================================

# clean workspace
rm(list = ls(all.names = T))
gc()

# load packages
library(circlize)
library(dplyr)

# ------------------------------------------------------------
# 1. directories and input files
# ------------------------------------------------------------

root <- '/home/yshin/Gloydius_ussuriensis_genome_assembly/R/Rdata/circos'
outdir <- '/home/yshin/Gloydius_ussuriensis_genome_assembly/R/Rplots'

track_file <- file.path(root, 'G_ussuriensis.100kb.circos_tracks.tsv')
length_file <- file.path(root, 'G_ussuriensis.chromosome_lengths.tsv')
class_file <- file.path(root, 'G_ussuriensis.chromosome_classes.tsv')
toxin_file <- file.path(root, 'G_ussuriensis.toxin_gene_locations.tsv')
out_png <- file.path(outdir, 'G_ussuriensis_circos_toxin_families.png')

dir.create(outdir, recursive = T, showWarnings = F)

for (f in c(track_file, length_file, class_file, toxin_file)) {
  if (!file.exists(f)) stop(paste('Input file not found:', f))
}

# ------------------------------------------------------------
# 2. plotting settings
# ------------------------------------------------------------

families_to_plot <- c('SVSP', 'SVMP', 'CTL', 'LAO', 'PLA2', '3FTx', 'CRISP', 'KUN')
expected_counts <- c(SVSP = 20, SVMP = 17, CTL = 16, LAO = 10, PLA2 = 5, '3FTx' = 5, CRISP = 4, KUN = 4)

# cluster grouping is only for visualizing nearby toxin gene loci
# plot one thin radial line per local family cluster to avoid thick patches
cluster_gap_bp <- 1000000

toxin_line_lwd <- 2.5
# vertical extent of toxin ticks within the chromosome ideogram (0 to 1)
toxin_tick_bottom <- 0.15
toxin_tick_top <- 0.85

# ------------------------------------------------------------
# 3. read input data
# ------------------------------------------------------------

dat <- read.delim(track_file, header = T, stringsAsFactors = F, check.names = F)
chrom <- read.delim(length_file, header = T, stringsAsFactors = F, check.names = F)
chr_class <- read.delim(class_file, header = T, stringsAsFactors = F, check.names = F)
toxin <- read.delim(toxin_file, header = T, stringsAsFactors = F, check.names = F)

# ------------------------------------------------------------
# 4. chromosome order and metadata
# ------------------------------------------------------------

chr_order <- paste0('G_ussuri_chr', c(1:7, 'Z', 'W', 9:18))

if (!all(chr_order %in% chrom$chr)) stop('Missing chromosomes from chromosome length table.')
if (!all(chr_order %in% chr_class$chr)) stop('Missing chromosomes from chromosome class table.')

chrom <- chrom[match(chr_order, chrom$chr), ]
chr_class <- chr_class[match(chr_order, chr_class$chr), ]

if (!all(chrom$chr == chr_class$chr)) stop('Chromosome metadata tables do not match.')

chrom$display <- chr_class$display
chrom$class <- chr_class$class

# ------------------------------------------------------------
# 5. validate and prepare toxin annotation data
# ------------------------------------------------------------

required_cols <- c('chr', 'start', 'end', 'family', 'gene_id')
if (!all(required_cols %in% names(toxin))) stop('Missing columns from toxin input table.')

toxin <- toxin %>%
  filter(family %in% families_to_plot) %>%
  select(chr, start, end, family, gene_id) %>%
  distinct() %>%
  arrange(match(chr, chr_order), start, end, gene_id) %>%
  mutate(midpoint = (start + end) / 2)

if (nrow(toxin) == 0) stop('No toxin loci match selected families.')
if (anyDuplicated(toxin$gene_id)) stop('Duplicate gene IDs found in toxin coordinate file.')
if (any(!toxin$chr %in% chr_order)) stop('Toxin input contains unrecognized chromosome names.')

chr_lengths <- setNames(as.numeric(chrom$length), chrom$chr)
if (any(!is.finite(toxin$start) | !is.finite(toxin$end) | toxin$start < 1 | toxin$end < toxin$start | toxin$end > chr_lengths[toxin$chr])) stop('Toxin coordinates are invalid or exceed chromosome lengths.')

# validate all nonredundant annotations before excluding Z-linked CTLs
# input now contains core candidates; validate the 81 core loci before removing Z-linked CTLs
observed_counts <- table(factor(toxin$family, levels = families_to_plot))
if (!all(as.integer(observed_counts) == unname(expected_counts[families_to_plot]))) stop('Toxin family counts differ from the validated core set.')

# exclude Z-linked CTLs from the figure only; source TSV is unchanged
n_z_ctl <- sum(toxin$chr == 'G_ussuri_chrZ' & toxin$family == 'CTL')
if (n_z_ctl != 3) stop('Expected 3 Z-linked core CTLs before plotting.')

toxin <- toxin %>%
  filter(!(chr == 'G_ussuri_chrZ' & family == 'CTL'))

if (nrow(toxin) != 78) stop('Expected 78 core toxin loci after excluding Z-linked CTLs.')

plotted_counts <- table(factor(toxin$family, levels = families_to_plot))
print(data.frame(family = families_to_plot, plotted = as.integer(plotted_counts)), row.names = F)
cat('Z-linked CTLs excluded from plot:', n_z_ctl, '\n')

# ------------------------------------------------------------
# 6. group nearby loci for thin family-position ticks
# ------------------------------------------------------------

make_clusters <- function(x, gap_bp) {
  
  groups <- split(x, interaction(x$chr, x$family, drop = T))
  
  result <- lapply(groups, function(z) {
    
    z <- z[order(z$start, z$end), ]
    
    new_cluster <- c(T, z$start[-1] > (head(cummax(z$end), -1) + gap_bp))
    
    z$cluster_id <- cumsum(new_cluster)
    
    z %>%
      group_by(cluster_id) %>%
      summarise(chr = first(chr), family = first(family), start = min(start), end = max(end), n_genes = n_distinct(gene_id), .groups = 'drop')
    
  })
  
  bind_rows(result) %>%
    mutate(midpoint = (start + end) / 2)
  
}

toxin_clusters <- make_clusters(toxin, cluster_gap_bp) %>%
  arrange(match(chr, chr_order), midpoint)

# ------------------------------------------------------------
# 7. chromosome and track colors
# ------------------------------------------------------------

macro_dark <- '#8f8f8f'
macro_light <- '#bdbdbd'
micro_dark <- '#4c9a67'
micro_light <- '#7fbc96'
z_col <- 'steelblue3'
w_col <- 'goldenrod2'

chrom$color <- NA_character_
chrom$shade <- NA_character_

macro_idx <- which(chrom$class == 'Macrochromosome')
micro_idx <- which(chrom$class == 'Microchromosome')

chrom$shade[macro_idx] <- rep(c('dark', 'light'), length.out = length(macro_idx))
chrom$shade[micro_idx] <- rep(c('dark', 'light'), length.out = length(micro_idx))
chrom$color[macro_idx] <- ifelse(chrom$shade[macro_idx] == 'dark', macro_dark, macro_light)
chrom$color[micro_idx] <- ifelse(chrom$shade[micro_idx] == 'dark', micro_dark, micro_light)
chrom$color[chrom$chr == 'G_ussuri_chrZ'] <- z_col
chrom$shade[chrom$chr == 'G_ussuri_chrZ'] <- 'dark'
chrom$color[chrom$chr == 'G_ussuri_chrW'] <- w_col
chrom$shade[chrom$chr == 'G_ussuri_chrW'] <- 'light'

chrom$track_bg <- ifelse(chrom$shade == 'dark', '#e1e1e1', '#f0f0f0')

figure_bg <- 'white'
line_col <- '#222222'
highlight_col <- '#ef4b35'
gene_low <- '#fff5f0'
gene_mid <- '#fcae91'
gene_high <- '#e34a33'

# colors for thin toxin-family position ticks
family_colors <- c(SVSP = '#B17A63', SVMP = '#168A72', CTL = '#253A70', LAO = '#C6B323', PLA2 = '#87384E', '3FTx' = '#438BC4', CRISP = '#D09E43', KUN = '#7960AE')

# ------------------------------------------------------------
# 8. prepare GC, repeat and gene density tracks
# ------------------------------------------------------------

dat <- dat %>%
  filter(chr %in% chr_order) %>%
  arrange(match(chr, chr_order), start) %>%
  mutate(mid = (start + end) / 2)

gc_median <- median(dat$GC_pct, na.rm = T)
repeat_median <- median(dat$Repeat_pct, na.rm = T)

gc_ylim <- quantile(dat$GC_pct, probs = c(0.01, 0.99), na.rm = T)
repeat_ylim <- quantile(dat$Repeat_pct, probs = c(0.01, 0.99), na.rm = T)

dat$GC_plot <- pmin(pmax(dat$GC_pct, gc_ylim[1]), gc_ylim[2])
dat$Repeat_plot <- pmin(pmax(dat$Repeat_pct, repeat_ylim[1]), repeat_ylim[2])

gene_cap <- as.numeric(quantile(dat$gene_count, probs = 0.99, na.rm = T))
if (!is.finite(gene_cap) | gene_cap <= 0) gene_cap <- max(dat$gene_count, na.rm = T)
if (!is.finite(gene_cap) | gene_cap <= 0) gene_cap <- 1

gene_col_fun <- colorRamp2(c(0, gene_cap / 2, gene_cap), c(gene_low, gene_mid, gene_high))

# ------------------------------------------------------------
# 9. chromosome axes and gaps
# ------------------------------------------------------------

get_tick_step <- function(length_bp) {
  if (length_bp >= 250e6) return(50e6)
  if (length_bp >= 150e6) return(25e6)
  if (length_bp >= 70e6) return(20e6)
  if (length_bp >= 30e6) return(10e6)
  return(5e6)
}

get_tick_positions <- function(length_bp, chr_type) {
  if (chr_type == 'Microchromosome') return(c(0, length_bp))
  seq(0, length_bp, by = get_tick_step(length_bp))
}

get_tick_labels <- function(tick_at, length_bp, chr_type) {
  if (chr_type == 'Microchromosome') return(c('', format(round(length_bp / 1e6, 1), trim = T, scientific = F)))
  format(round(tick_at / 1e6, 1), trim = T, scientific = F)
}

get_tick_cex <- function(chr_type, length_bp) {
  if (chr_type == 'Microchromosome') return(0.50)
  if (length_bp < 70e6) return(0.54)
  if (length_bp < 150e6) return(0.58)
  return(0.62)
}

get_cap_fraction <- function(length_bp) {
  if (length_bp < 30e6) return(0.085)
  if (length_bp < 70e6) return(0.055)
  if (length_bp < 150e6) return(0.032)
  return(0.025)
}


# chromosome sector gaps
gap_after <- rep(2.2, length(chr_order))
names(gap_after) <- chr_order

# microchromosome spacing
gap_after[paste0('G_ussuri_chr', 9:18)] <- 2.8

# larger gaps between chromosome classes
gap_after['G_ussuri_chr7'] <- 4.5
gap_after['G_ussuri_chrW'] <- 5.0

# final sector gap
gap_after['G_ussuri_chr18'] <- 6.0



# ------------------------------------------------------------
# 11. draw Circos plot
# ------------------------------------------------------------

draw_circos <- function() {
  
  circos.clear()
  
  
  # global Circos parameters
  circos.par(start.degree = 90, gap.after = gap_after, track.margin = c(0.004, 0.004), cell.padding = c(0, 0, 0, 0), points.overflow.warning = F, canvas.xlim = c(-1.22, 1.22), canvas.ylim = c(-1.22, 1.22))
  circos.initialize(factors = factor(chrom$chr, levels = chr_order), xlim = cbind(rep(0, nrow(chrom)), chrom$length))
  
  # ----------------------------------------------------------
  # chromosome ideogram and toxin family highlights
  # ----------------------------------------------------------
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.068, bg.border = NA, panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    idx <- match(chr, chrom$chr)
    chr_length <- chrom$length[idx]
    chr_type <- chrom$class[idx]
    cap_pad <- chr_length * get_cap_fraction(chr_length)
    x1 <- CELL_META$xlim[1] + cap_pad
    x2 <- CELL_META$xlim[2] - cap_pad
    
    # original rounded chromosome ideogram
    circos.lines(seq(x1, x2, length.out = 500), rep(0.5, 500), col = chrom$color[idx], lwd = 15)
    
    # one thin radial tick at each local toxin-family cluster midpoint
    # avoids filled patches when many tandem toxin genes occur close together
    ctmp <- toxin_clusters[toxin_clusters$chr == chr, ]
    
    if (nrow(ctmp) > 0) {
      
      for (i in seq_len(nrow(ctmp))) {
        
        circos.segments(ctmp$midpoint[i], toxin_tick_bottom, ctmp$midpoint[i], toxin_tick_top, col = unname(family_colors[ctmp$family[i]]), lwd = toxin_line_lwd)
        
      }
      
    }
    
    tick_at <- get_tick_positions(chr_length, chr_type)
    tick_labels <- get_tick_labels(tick_at, chr_length, chr_type)
    
    circos.axis(h = 'top', major.at = tick_at, labels = tick_labels, labels.cex = get_tick_cex(chr_type, chr_length), labels.facing = 'clockwise', labels.niceFacing = T, major.tick.length = 0.10, minor.ticks = 4, lwd = 0.80)
    
  })
  
  # ----------------------------------------------------------
  # GC content
  # ----------------------------------------------------------
  
  circos.trackPlotRegion(ylim = gc_ylim, track.height = 0.105, bg.col = NA, bg.border = NA, panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    idx <- match(chr, chrom$chr)
    
    circos.rect(CELL_META$xlim[1], CELL_META$ylim[1], CELL_META$xlim[2], CELL_META$ylim[2], col = chrom$track_bg[idx], border = NA)
    
    tmp <- dat[dat$chr == chr, ]
    
    circos.lines(tmp$mid, tmp$GC_plot, col = line_col, lwd = 0.40)
    circos.lines(tmp$mid, ifelse(tmp$GC_pct > gc_median, tmp$GC_plot, NA), col = highlight_col, lwd = 0.65)
    
  })
  
  # ----------------------------------------------------------
  # repeat content
  # ----------------------------------------------------------
  
  circos.trackPlotRegion(ylim = repeat_ylim, track.height = 0.105, bg.col = NA, bg.border = NA, panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    idx <- match(chr, chrom$chr)
    
    circos.rect(CELL_META$xlim[1], CELL_META$ylim[1], CELL_META$xlim[2], CELL_META$ylim[2], col = chrom$track_bg[idx], border = NA)
    
    tmp <- dat[dat$chr == chr, ]
    
    circos.lines(tmp$mid, tmp$Repeat_plot, col = line_col, lwd = 0.40)
    circos.lines(tmp$mid, ifelse(tmp$Repeat_pct > repeat_median, tmp$Repeat_plot, NA), col = highlight_col, lwd = 0.65)
    
  })
  
  # ----------------------------------------------------------
  # gene density heatmap
  # ----------------------------------------------------------
  
  circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.058, bg.col = gene_low, bg.border = NA, panel.fun = function(x, y) {
    
    chr <- CELL_META$sector.index
    tmp <- dat[dat$chr == chr, ]
    
    gene_values <- pmin(tmp$gene_count, gene_cap)
    circos.rect(tmp$start, 0, tmp$end, 1, col = gene_col_fun(gene_values), border = NA)
    
  })
  
  # original chromosome class legend
  legend(x = 'bottom', legend = c('Macrochromosome', 'Microchromosome', 'Chr Z', 'Chr W'), col = c(macro_dark, micro_dark, z_col, w_col), pch = 15, pt.cex = 1.6, cex = 0.70, horiz = T, bty = 'n', inset = c(0, -0.075), xpd = NA)
  
  circos.clear()
  
}

# ============================================================
# 12. save PNG
# ============================================================

png(filename = out_png, width = 3600, height = 3600, res = 360, bg = figure_bg)
par(mar = c(5.5, 2, 2, 2), bg = figure_bg, lend = 'round', ljoin = 'round')

draw_circos()

dev.off()

# ============================================================
# 13. final report
# ============================================================

cat('\n')
cat('============================================================\n')
cat('G. ussuriensis Circos toxin family plot completed\n')
cat('============================================================\n\n')

cat('Toxin loci plotted:', nrow(toxin), '\n')
cat('Toxin family clusters marked:', nrow(toxin_clusters), '\n')
cat('Z-linked CTLs excluded:', n_z_ctl, '\n')
cat('GC median:', round(gc_median, 3), '%\n')
cat('Repeat median:', round(repeat_median, 3), '%\n')
cat('Gene density cap:', round(gene_cap, 3), 'genes / 100 kb\n')
cat('PNG:', out_png, '\n\n')

