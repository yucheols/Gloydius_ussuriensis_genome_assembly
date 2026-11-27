# ============================================================
#
# PLA2 comparative microsynteny
# G. ussuriensis, C. adamanteus, C. gasperettii, N. naja
#
# ============================================================

# clean workspace
rm(list = ls(all.names = T))
gc()

# load packages
library(ggplot2)
library(gggenes)
library(dplyr)
library(readr)
library(stringr)

# ------------------------------------------------------------
# 1. directories and plotting options
# ------------------------------------------------------------

root <- '//wsl.localhost/Ubuntu/home/yshin/Gloydius_ussuriensis_genome_assembly'
if (!dir.exists(root)) root <- '/home/yshin/Gloydius_ussuriensis_genome_assembly'

datadir <- file.path(root, 'R', 'Rdata', 'toxin_microsynteny')
outdir <- file.path(root, 'R', 'Rplots')
dir.create(outdir, recursive = T, showWarnings = F)

interval_file <- file.path(datadir, 'anchor_bounded_regions.tsv')
gene_file <- file.path(datadir, 'comparative_all_gene_models_in_regions.tsv')

for (f in c(interval_file, gene_file)) {
  if (!file.exists(f)) stop('Missing input file: ', f)
}

# No species are flipped for PLA2
FLIP_SPECIES <- character(0)

species_order <- c('Gloydius_ussuriensis', 'Crotalus_adamanteus', 'Cerastes_gasperettii', 'Naja_naja')
species_labels <- c('G. ussuriensis', 'C. adamanteus', 'C. gasperettii', 'N. naja')
names(species_labels) <- species_order

colors <- c('Flanking genes' = '#9FA6AD', 'PLA2s' = '#D87348', 'PLA2G2F' = '#8E629C', 'Other PLA2-related' = '#B49AC0')

# ------------------------------------------------------------
# 2. import genomic intervals and gene catalog
# ------------------------------------------------------------

read_table <- function(path) {
  read_tsv(path, col_types = cols(.default = col_character()), show_col_types = F, progress = F)
}

intervals <- read_table(interval_file)
genes <- read_table(gene_file)

required <- list(intervals = c('family', 'species', 'chromosome', 'start', 'end'), genes = c('family', 'species', 'chr', 'gene_id', 'start', 'end', 'strand', 'plot_role', 'plot_include', 'anchor_symbol'))

for (nm in names(required)) {
  missing <- setdiff(required[[nm]], names(get(nm)))
  if (length(missing)) stop(nm, ' missing columns: ', paste(missing, collapse = ', '))
}

intervals <- intervals %>%
  filter(family == 'PLA2', species %in% species_order) %>%
  mutate(start = as.numeric(start), end = as.numeric(end), length_bp = end - start + 1, flip = species %in% FLIP_SPECIES)

if (nrow(intervals) != 4) stop('Expected four PLA2 interval entries.')
if (anyDuplicated(intervals[c('family', 'species')])) stop('Duplicated PLA2 interval entry.')
if (any(is.na(intervals$length_bp) | intervals$length_bp <= 0)) stop('Invalid PLA2 interval coordinates.')

genes <- genes %>%
  filter(family == 'PLA2', species %in% species_order, plot_include == '1') %>%
  mutate(start = as.numeric(start), end = as.numeric(end)) %>%
  left_join(intervals %>% select(family, species, interval_chr = chromosome, interval_start = start, interval_end = end, length_bp, flip), by = c('family', 'species'))

if (any(is.na(genes$interval_start))) stop('Some gene records lack a PLA2 interval.')
if (any(is.na(genes$chr) | genes$chr != genes$interval_chr)) stop('Gene and interval chromosomes disagree.')

genes <- genes %>%
  filter(end >= interval_start, start <= interval_end) %>%
  mutate(local_start = pmax(start, interval_start) - interval_start + 1, local_end = pmin(end, interval_end) - interval_start + 1, plot_start = if_else(flip, length_bp - local_end + 1, local_start), plot_end = if_else(flip, length_bp - local_start + 1, local_end), forward = if_else(flip, strand == '-', strand == '+'))

# ------------------------------------------------------------
# 3. classify plotted genes
# ------------------------------------------------------------

anchor_roles <- c('conserved_anchor_candidate', 'provisional_VWA5B1_anchor_group')

# PLA2G2F is identified using the anchor_symbol field.
# All four currently identified PLA2G2F candidates are
# displayed in purple.
#
# PLA2-1 / PLA2I-like-1 remains excluded from enzymatic
# PLA2 counts if present in the catalog.
#
# Other ambiguous PLA2-related models remain separate.
# No missing flanking anchors are inferred.

genes <- genes %>%
  mutate(
    is_inhibitor_like = species == 'Gloydius_ussuriensis' & gene_id %in% c('PLA2-1', 'PLA2I-like-1'),
    is_pla2g2f = str_to_upper(coalesce(anchor_symbol, '')) == 'PLA2G2F',
    plot_class = case_when(
      is_pla2g2f ~ 'PLA2G2F',
      is_inhibitor_like ~ 'Other PLA2-related',
      plot_role == 'ambiguous_anchor_or_venom_PLA2' ~ 'Other PLA2-related',
      plot_role %in% anchor_roles & !is.na(anchor_symbol) & anchor_symbol != '' ~ 'Flanking genes',
      plot_role == 'putative_venom_toxin' ~ 'PLA2s',
      plot_role == 'PLA2_like_candidate' ~ 'PLA2s',
      TRUE ~ 'Other annotated gene'
    ),
    anchor_label = case_when(
      is_pla2g2f ~ 'PLA2G2F',
      is_inhibitor_like ~ 'PLA2I-like',
      plot_role == 'ambiguous_anchor_or_venom_PLA2' ~ 'PLA2-like?',
      plot_role %in% anchor_roles ~ anchor_symbol,
      TRUE ~ NA_character_
    )
  )

# ------------------------------------------------------------
# 4. remove non-plotted classes
# ------------------------------------------------------------

genes <- genes %>%
  filter(plot_class != 'Other annotated gene')

if (!all(genes$plot_class %in% names(colors))) {
  stop('Missing plot colors for: ', paste(setdiff(unique(genes$plot_class), names(colors)), collapse = ', '))
}

if (!nrow(genes)) stop('No PLA2 genes selected for plotting.')

# Verify PLA2G2F classification

cat('===== PLA2G2F GENE MODELS =====\n')

pla2g2f_check <- genes %>%
  filter(plot_class == 'PLA2G2F')

print(pla2g2f_check %>% select(species, gene_id, plot_role, anchor_symbol, start, end), n = Inf)

# Verify inhibitor-like candidate classification

inhibitor_check <- genes %>%
  filter(is_inhibitor_like)

cat('===== GLOYDIUS PLA2 INHIBITOR-LIKE CHECK =====\n')

if (nrow(inhibitor_check)) {
  print(inhibitor_check %>% select(species, gene_id, plot_class), n = Inf)
} else {
  warning('PLA2-1 / PLA2I-like-1 was not identified among the included gene models. Check its gene_id in the input table.')
}

# ------------------------------------------------------------
# 5. audit PLA2 candidate models
# ------------------------------------------------------------

# Keep all distinct annotated toxin candidates.
# Unlike SVMP, PLA2 candidates are not automatically collapsed
# solely because their genomic intervals overlap.

pla2_candidates <- genes %>%
  filter(plot_class == 'PLA2s') %>%
  arrange(species, start, end)

cat('===== PLA2 CANDIDATE MODELS =====\n')
print(pla2_candidates %>% select(species, gene_id, chr, start, end, plot_role), n = Inf)

# Inspect overlapping toxin models without deleting them

overlap_check <- pla2_candidates %>%
  group_by(species) %>%
  arrange(start, end, .by_group = T) %>%
  mutate(previous_max_end = lag(cummax(end)), overlaps_previous = !is.na(previous_max_end) & start <= previous_max_end) %>%
  ungroup()

if (any(overlap_check$overlaps_previous)) {
  cat('===== POTENTIALLY OVERLAPPING PLA2 MODELS =====\n')
  print(overlap_check %>% filter(overlaps_previous) %>% select(species, gene_id, start, end), n = Inf)
}

# ------------------------------------------------------------
# 6. align species using relative genomic intervals
# ------------------------------------------------------------

# Unlike the SVMP region, conserved anchor pairs have not
# been identified for every species in the PLA2 region.
#
# Each extracted PLA2 interval is therefore scaled to 0-1.
# Internal spacing within each interval is retained
# proportionally.
#
# This is a relative-position display, not a shared bp scale.

genes <- genes %>%
  mutate(plot_start_bp = plot_start, plot_end_bp = plot_end, plot_start = plot_start_bp / length_bp, plot_end = plot_end_bp / length_bp, center = (plot_start + plot_end) / 2)

intervals <- intervals %>%
  mutate(track_start = 0, track_end = 1)

y_map <- c('Gloydius_ussuriensis' = 3.1, 'Crotalus_adamanteus' = 2.4, 'Cerastes_gasperettii' = 1.7, 'Naja_naja' = 1.0)

genes <- genes %>%
  mutate(y = unname(y_map[species]))

intervals <- intervals %>%
  mutate(y = unname(y_map[species]))

gene_models <- genes

cat('===== PLA2 REGION ALIGNMENT =====\n')
print(intervals %>% select(species, chromosome, start, end, length_bp, flip, track_start, track_end), n = Inf)

# ------------------------------------------------------------
# 7. audit identified flanking anchors
# ------------------------------------------------------------

expected_anchors <- c('OTUD3', 'MUL1', 'VWA5B1', 'PLA2G2F')

anchor_presence <- expand.grid(species = species_order, anchor_symbol = expected_anchors, stringsAsFactors = F) %>%
  as_tibble() %>%
  left_join(genes %>% filter(!is.na(anchor_symbol), anchor_symbol %in% expected_anchors, plot_role %in% anchor_roles) %>% count(species, anchor_symbol, name = 'identified_models'), by = c('species', 'anchor_symbol')) %>%
  mutate(identified_models = coalesce(identified_models, 0L), status = if_else(identified_models > 0, 'identified', 'not_identified_in_catalog'))

cat('===== PLA2 FLANKING ANCHOR STATUS =====\n')
print(anchor_presence, n = Inf)

write_tsv(anchor_presence, file.path(outdir, 'PLA2_anchor_presence_revised.tsv'))

# ------------------------------------------------------------
# 8. draw plot
# ------------------------------------------------------------

plot_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = 'bold', size = 18),
    plot.title.position = 'plot',
    axis.text.y = element_text(face = 'italic', size = 15),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.text = element_text(size = 13),
    legend.title = element_blank(),
    legend.key.width = grid::unit(1.0, 'cm'),
    legend.key.height = grid::unit(0.75, 'cm'),
    legend.spacing.x = grid::unit(0.6, 'cm'),
    plot.margin = margin(10, 18, 8, 12)
  )

# chromosome midlines

p <- ggplot() +
  geom_segment(data = intervals, aes(x = track_start, xend = track_end, y = y, yend = y), color = '#929BA3', linewidth = 1.5)

# ------------------------------------------------------------
# arrows with flush arrowheads and black outlines
# ------------------------------------------------------------

if (nrow(gene_models) > 0) {
  p <- p + gggenes::geom_gene_arrow(data = gene_models, aes(xmin = plot_start, xmax = plot_end, y = y, fill = plot_class, forward = forward), arrowhead_width = grid::unit(1.3, 'mm'), arrowhead_height = grid::unit(6.0, 'mm'), arrow_body_height = grid::unit(6.0, 'mm'), colour = 'black', linewidth = 0.7)
}

plot_range <- range(c(intervals$track_start, intervals$track_end, genes$plot_start, genes$plot_end), na.rm = T)

p <- p +
  scale_fill_manual(values = colors, breaks = c('Flanking genes', 'PLA2s', 'PLA2G2F', 'Other PLA2-related'), name = NULL) +
  scale_x_continuous(breaks = NULL, limits = c(min(-0.04, plot_range[1] - 0.02), max(1.04, plot_range[2] + 0.02)), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(3.1, 2.4, 1.7, 1.0), labels = species_labels[species_order], limits = c(0.5, 3.6), expand = c(0, 0)) +
  labs(title = '(B)', x = NULL, y = NULL) +
  plot_theme +
  guides(fill = guide_legend(nrow = 1, byrow = T, title = NULL))

print(p)

# ------------------------------------------------------------
# 9. export outputs
# ------------------------------------------------------------

plot_prefix <- file.path(outdir, 'PLA2_gene_overview_revised')

ggsave(paste0(plot_prefix, '.png'), p, width = 14, height = 5.2, dpi = 300)
ggsave(paste0(plot_prefix, '.pdf'), p, width = 14, height = 5.2)

write_tsv(genes %>% select(family, species, chr, gene_id, start, end, strand, plot_role, plot_class, anchor_symbol, anchor_label, is_inhibitor_like, plot_start_bp, plot_end_bp, plot_start, plot_end, forward), file.path(outdir, 'PLA2_displayed_gene_models_revised.tsv'))

cat('===== PLA2 PLOT SUMMARY =====\n')
print(genes %>% count(species, plot_class), n = Inf)

cat('Output plot:', paste0(plot_prefix, '.png'), '\n')
cat('Output PDF:', paste0(plot_prefix, '.pdf'), '\n')
cat('Output table:', file.path(outdir, 'PLA2_displayed_gene_models_revised.tsv'), '\n')
cat('Anchor audit:', file.path(outdir, 'PLA2_anchor_presence_revised.tsv'), '\n')
cat('===== COMPLETE =====\n')
