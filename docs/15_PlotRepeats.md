## 15) Plotting the repeat landscape
We already ran Earl Grey prior to annotation and this process annotated the repeat and TE elements in our genome assembly. This process also generated some useful plots as pdf files. Let's scp these files and have a look:
```sh
# first make a directory to put these files in
mkdir /home/yshin/Gloydius_ussuriensis_genome_assembly/outputs/repeat_landscape
cd /home/yshin/Gloydius_ussuriensis_genome_assembly/outputs/repeat_landscape

# scp the files
scp yshin@mendel.sdmz.amnh.org:/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis_classification_landscape.pdf .

scp yshin@mendel.sdmz.amnh.org:/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.summaryPie.pdf .
```

Let's recreate these plots in R. scp the required files locally:
```sh
# cd into destination dir
cd /home/yshin/Gloydius_ussuriensis_genome_assembly/outputs/repeat_landscape

# mendel Earl Grey directory
REMOTE='/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/annotation/soft_masked/Gloydius_ussuriensis_EarlGrey'

# download repeat divergence GFF
scp -C "yshin@mendel.sdmz.amnh.org:${REMOTE}/Gloydius_ussuriensis_RepeatLandscape/Gloydius_ussuriensis.filteredRepeats.withDivergence.gff" .

# download repeat composition table
scp -C "yshin@mendel.sdmz.amnh.org:${REMOTE}/Gloydius_ussuriensis_summaryFiles/Gloydius_ussuriensis.highLevelCount.txt" .
```

...and use the R script below to generate the plot:
```sh
# ============================================================
# 
# plot G. ussuriensis repeat landscape
# 
# ============================================================

# clean workspace
rm(list = ls(all.names = T))
gc()

# load packages
library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(cowplot)

# ------------------------------------------------------------
# 1. directories and input files (Windows RStudio -> WSL)
# ------------------------------------------------------------

root <- '//wsl.localhost/Ubuntu/home/yshin/Gloydius_ussuriensis_genome_assembly/outputs/repeat_landscape'
outdir <- file.path(root, 'plots')
dir.create(outdir, recursive = T, showWarnings = F)

gff_file <- file.path(root, 'Gloydius_ussuriensis.filteredRepeats.withDivergence.gff')
summary_file <- file.path(root, 'Gloydius_ussuriensis.highLevelCount.txt')

if (!file.exists(gff_file)) stop(paste('GFF not found:', gff_file))
if (!file.exists(summary_file)) stop(paste('Summary table not found:', summary_file))


# ------------------------------------------------------------
# 2. classifications and colors matching the screenshot
# ------------------------------------------------------------

class_order <- c('DNA', 'Rolling Circle', 'Penelope', 'LINE', 'SINE', 'LTR', 'Other', 'Unclassified', 'Non-Repeat')

repeat_colors <- c('DNA' = '#FF6B4A',
                   'Rolling Circle' = '#FF9A13',
                   'Penelope' = '#D2A6F7',
                   'LINE' = '#37C9EF',
                   'SINE' = '#F578CB',
                   'LTR' = '#5FD76E',
                   'Other' = '#FAA4C2',
                   'Unclassified' = '#646A77',
                   'Non-Repeat' = '#D4DAE5')

normalize_class <- function(x) {
  x <- as.character(x)
  case_when(grepl('Penelope|^PLE', x, ignore.case = T) ~ 'Penelope',
            grepl('^DNA', x, ignore.case = T) ~ 'DNA',
            grepl('^RC|^Rolling Circle', x, ignore.case = T) ~ 'Rolling Circle',
            grepl('^LINE', x, ignore.case = T) ~ 'LINE',
            grepl('^SINE', x, ignore.case = T) ~ 'SINE',
            grepl('^LTR', x, ignore.case = T) ~ 'LTR',
            grepl('^Unknown|^Unclassified|^Retroposon|^Unspecified', x, ignore.case = T) ~ 'Unclassified',
            grepl('^Non-Repeat', x, ignore.case = T) ~ 'Non-Repeat',
            T ~ 'Other')
}


# ------------------------------------------------------------
# 3. repeat composition: non-overlapping genome fractions
# ------------------------------------------------------------

counts <- read_tsv(summary_file, show_col_types = F)
required <- c('TE Classification', 'Coverage (bp)', 'Genome Size')
if (!all(required %in% names(counts))) stop(paste('Unexpected table columns:', paste(names(counts), collapse = ', ')))

genome_size <- unique(as.numeric(counts[['Genome Size']][!is.na(counts[['Genome Size']])]))
if (length(genome_size) != 1 || !is.finite(genome_size) || genome_size <= 0) stop('Invalid genome size.')

composition <- counts %>%
  transmute(raw_class = as.character(.data[['TE Classification']]), bp = as.numeric(.data[['Coverage (bp)']])) %>%
  filter(!grepl('-nested$', raw_class, ignore.case = T),
         !raw_class %in% c('Total Interspersed Repeat', 'Non-Repeat')) %>%
  mutate(Class = normalize_class(raw_class)) %>%
  group_by(Class) %>%
  summarise(bp = sum(bp, na.rm = T), .groups = 'drop')

reported_interspersed <- counts[['Coverage (bp)']][counts[['TE Classification']] == 'Total Interspersed Repeat']
calculated_interspersed <- sum(composition$bp[composition$Class != 'Other'])
if (length(reported_interspersed) != 1 || abs(calculated_interspersed - reported_interspersed) > 1) stop('Interspersed repeat coverage does not match Earl Grey table.')

nonrepeat_bp <- genome_size - sum(composition$bp)
if (nonrepeat_bp < 0) stop('Repeat coverage exceeds genome size.')

composition <- bind_rows(composition, tibble(Class = 'Non-Repeat', bp = nonrepeat_bp)) %>%
  mutate(pct = 100 * bp / genome_size,
         Class = factor(Class, levels = class_order)) %>%
  arrange(Class)

if (abs(sum(composition$pct) - 100) > 1e-6) stop('Repeat composition does not sum to 100%.')


# ------------------------------------------------------------
# 4. read Earl Grey divergence GFF
# ------------------------------------------------------------

gff <- read_tsv(gff_file,
                col_names = c('seqid', 'method', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'),
                col_types = cols(.default = col_character(), start = col_double(), end = col_double()),
                comment = '#', show_col_types = F)

kimura_text <- str_match(gff$attributes, regex('(?:^|;)KIMURA80=([^;]+)', ignore_case = T))[, 2]
gff$kimura <- suppressWarnings(as.numeric(kimura_text))
if (all(is.na(gff$kimura))) stop('No KIMURA80 values found in divergence GFF.')

# Earl Grey: sum base pairs in each 0.01-divergence bin
landscape <- gff %>%
  mutate(Class = normalize_class(type)) %>%
  filter(is.finite(kimura), kimura >= 0, kimura <= 0.5, end >= start) %>%
  transmute(Class, kimura_bin = round(kimura, 2), bp = end - start + 1) %>%
  group_by(Class, kimura_bin) %>%
  summarise(bp = sum(bp), .groups = 'drop') %>%
  mutate(Class = factor(Class, levels = class_order))

if (nrow(landscape) == 0) stop('No valid divergence data were found.')

# check repeat classifications
cat('\n===== REPEAT LANDSCAPE CLASSIFICATIONS =====\n')

landscape_summary <- landscape %>%
  group_by(Class) %>%
  summarise(total_bp = sum(bp), .groups = 'drop') %>%
  arrange(desc(total_bp))

print(landscape_summary, n = Inf)

if (!'Unclassified' %in% as.character(landscape$Class)) {
  warning('Unclassified repeats were not detected in the divergence GFF.')
}


# ------------------------------------------------------------
# 5. landscape (behind the pie inset)
# ------------------------------------------------------------

plot_landscape <- ggplot(landscape, aes(x = kimura_bin, y = bp, fill = Class)) +
  geom_col(width = 0.0098, colour = NA) +
  scale_fill_manual(values = repeat_colors, limits = class_order, drop = F) +
  scale_x_reverse(limits = c(0.505, -0.005), breaks = seq(0, 0.5, 0.1), labels = function(x) sprintf('%.1f', x), expand = c(0, 0)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 5), labels = function(x) ifelse(x == 0, '0', formatC(x, format = 'e', digits = 0)), expand = expansion(mult = c(0, 0.04))) +
  labs(x = 'Kimura 2-parameter distance', y = 'Base pairs') +
  theme_classic(base_size = 11) +
  theme(axis.line = element_line(linewidth = 0.45),
        axis.text = element_text(size = 14),
        axis.ticks = element_line(linewidth = 0.4),
        axis.ticks.length = grid::unit(0.13, 'cm'),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15)),
        legend.position = 'none',
        plot.margin = margin(5, 0, 0, 5))


# ------------------------------------------------------------
# 6. pie chart (no donut, matching the screenshot)
# ------------------------------------------------------------

plot_pie <- ggplot(composition, aes(x = 1, y = pct, fill = Class)) +
  geom_col(width = 1, colour = NA) +
  coord_polar(theta = 'y', start = 0) +
  scale_fill_manual(values = repeat_colors, limits = class_order, drop = F) +
  theme_void() +
  theme(legend.position = 'none',
        plot.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        panel.background = element_rect(fill = 'transparent', colour = NA))




# ------------------------------------------------------------
# 7. assemble figure
# ------------------------------------------------------------

fig <- ggdraw() +
  draw_plot(plot_landscape, x = 0.035, y = 0.055, width = 0.595, height = 0.855) +
  draw_plot(plot_pie, x = 0.085, y = 0.500, width = 0.230, height = 0.410) +
  draw_label('TE Classification', x = 0.660, y = 0.800, hjust = 0, fontface = 'bold', size = 14)

# right-hand classification legend with percentages
for (i in seq_along(class_order)) {
  this_class <- class_order[i]
  this_pct <- composition$pct[as.character(composition$Class) == this_class]
  if (length(this_pct) == 0) next
  
  y <- 0.735 - (i - 1) * 0.055
  
  fig <- fig +
    draw_grob(grid::rectGrob(gp = grid::gpar(fill = repeat_colors[[this_class]], col = NA)), x = 0.663, y = y - 0.018, width = 0.020, height = 0.036) +
    draw_label(sprintf('%s (%.1f%%)', this_class, this_pct), x = 0.695, y = y, hjust = 0, size = 11)
}

fig <- fig +
  theme(plot.background = element_rect(fill = 'white', colour = NA),
        panel.background = element_rect(fill = 'white', colour = NA))

print(fig)


# ------------------------------------------------------------
# 8. save images and supporting tables
# ------------------------------------------------------------

ggsave(file.path(outdir, 'G_ussuriensis_repeat_landscape.png'), plot = fig, width = 15, height = 6, dpi = 400, bg = 'white')

write_csv(landscape, file.path(outdir, 'G_ussuriensis_landscape_data.csv'))
write_csv(composition, file.path(outdir, 'G_ussuriensis_composition_data.csv'))

cat('\nGenome size:', format(genome_size, big.mark = ','), 'bp\n')
cat('Total repeat content:', sprintf('%.2f%%', 100 - composition$pct[composition$Class == 'Non-Repeat']), '\n')
cat('Composition total:', sprintf('%.5f%%', sum(composition$pct)), '\n')
cat('Output directory:', outdir, '\n')
```