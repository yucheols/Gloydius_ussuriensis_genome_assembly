
if (!requireNamespace("RIdeogram", quietly=TRUE)) {
  if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
  devtools::install_github("TickingClock1992/RIdeogram")
}
library(RIdeogram)
if (!requireNamespace("rsvg", quietly=TRUE)) install.packages("rsvg")
library(rsvg)

setwd("/home/yshin/mendel-nas1/snake_genome_ass/G_ussuriensis_Chromo/synteny/Synk/outputs/pairwise/G_ussuri_vs_C_gas")
print(getwd())
print(list.files())

kary <- read.table("dual_karyotype.txt", header=TRUE)
synt <- read.table("final_synteny.txt", header=TRUE, colClasses = c("numeric", "integer", "integer", "numeric", "integer", "integer", "character"))

ideogram(karyotype=kary, synteny=synt)
rsvg_png("chromosome.svg", "chromosome.png", width=1000)
