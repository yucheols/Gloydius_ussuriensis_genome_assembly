# ============================================================
# plot PSMC demographic history with bootstrap replicates
#
# PSMC run info:
#   primary N500 estimate
#   100 N500 bootstrap replicates
#
# scaling:
#   mutation rate   = 1.25e-8 / site / generation
#   generation time = 3 years
#
# NOTE:
# psmc_plot.pl output:
#   column 1 = years before present
#   column 2 = Ne / 10^4
# ============================================================


# ------------------------------------------------------------
# clean workspace
# ------------------------------------------------------------

rm(list = ls(all.names = T))
gc()


# ------------------------------------------------------------
# packages
# ------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(scales)


# ------------------------------------------------------------
# paths
# ------------------------------------------------------------

root <- '/home/yshin/Gloydius_ussuriensis_genome_assembly/R/Rdata/psmc_bootstrap'
prefix <- file.path(root, 'AMNH_21010_PSMC_bootstrap')


# ------------------------------------------------------------
# plotting range
#
# x-axis is plotted on a log10 scale.
# y-axis is plotted in units of 10^4 Ne.
#
# plot_max_ne controls the y-axis crop.
# ------------------------------------------------------------

plot_min_year <- 2e4
plot_max_year <- 2e6
plot_max_ne <- 50


# ------------------------------------------------------------
# function to read psmc_plot.pl output
# ------------------------------------------------------------

read_psmc <- function(file, replicate, type) {
  
  x <- read.table(file, header = F, stringsAsFactors = F)
  
  if (ncol(x) < 2) {
    stop(paste('Unexpected PSMC plot file:', file))
  }
  
  data.frame(years = x[[1]],
             Ne = x[[2]],
             replicate = replicate,
             type = type)
}


# ------------------------------------------------------------
# primary PSMC estimate
# ------------------------------------------------------------

primary <- read_psmc(paste0(prefix, '.0.txt'),
                     replicate = 0,
                     type = 'Primary')


# ------------------------------------------------------------
# bootstrap replicates
# ------------------------------------------------------------

bootstrap_list <- lapply(1:100, function(i) {
  
  read_psmc(paste0(prefix, '.', i, '.txt'),
            replicate = i,
            type = 'Bootstrap')
})

bootstrap <- bind_rows(bootstrap_list)


# ------------------------------------------------------------
# restrict data to plotting interval
# ------------------------------------------------------------

primary_plot <- primary %>%
  filter(years >= plot_min_year,
         years <= plot_max_year)

bootstrap_plot <- bootstrap %>%
  filter(years >= plot_min_year,
         years <= plot_max_year)


# ------------------------------------------------------------
# colors
# ------------------------------------------------------------

bootstrap_color <- '#91CBC8'
primary_color <- '#4FA7A0'


# ------------------------------------------------------------
# plot
#
# Bootstrap replicates are shown individually as thin,
# transparent turquoise step curves.
#
# The primary N500 estimate is shown as a darker turquoise
# step curve on top.
# ------------------------------------------------------------

p <- ggplot() +
  geom_step(data = bootstrap_plot,
            aes(x = years, y = Ne, group = replicate), color = bootstrap_color, linewidth = 0.6, alpha = 0.15, direction = 'hv') +
  geom_step(data = primary_plot, aes(x = years, y = Ne), color = primary_color, linewidth = 1.6, direction = 'hv') +
  scale_x_log10(limits = c(plot_min_year, plot_max_year),
                breaks = 10^(4:6),
                labels = c(expression(10^4), expression(10^5), expression(10^6))) +
  scale_y_continuous(limits = c(0, plot_max_ne), expand = expansion(mult = c(0, 0.02))) +
  labs(x = 'Time since present (years)', y = expression('Effective population size (x' * 10^4 * ')')) +
  theme_classic(base_size = 14) +
  theme(axis.title = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 13),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = 0.7),
        panel.border = element_rect(fill = NA, color = 'black', linewidth = 1.0),
        plot.margin = margin(10, 15, 10, 10))


# ------------------------------------------------------------
# display
# ------------------------------------------------------------

print(p)


# ------------------------------------------------------------
# save
# ------------------------------------------------------------

ggsave(filename = file.path('/home/yshin/Gloydius_ussuriensis_genome_assembly/R/Rplots/AMNH_21010_PSMC_bootstrap.png'),
       plot = p,
       width = 6,
       height = 5,
       dpi = 800)
