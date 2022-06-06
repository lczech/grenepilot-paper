# Basic plotting
library(data.table)
library(tidyr)
library(ggplot2)
library("GGally")
library(cowplot)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# square size for all plots
ssize=12

print("reading data")
seeds <- fread( "frequency.csv", sep="\t", select=c(
#seeds <- fread( "head.csv", sep="\t", select=c(
    "CHROM_POS", "S1.FREQ", "S2.FREQ", "S3.FREQ", "S4.FREQ", "S5.FREQ", "S6.FREQ", "S7.FREQ", "S8.FREQ"
))

# turn into minor freq
#seeds <- 1 - seeds
seeds$S1.FREQ <- 1 - seeds$S1.FREQ
seeds$S2.FREQ <- 1 - seeds$S2.FREQ
seeds$S3.FREQ <- 1 - seeds$S3.FREQ
seeds$S4.FREQ <- 1 - seeds$S4.FREQ
seeds$S5.FREQ <- 1 - seeds$S5.FREQ
seeds$S6.FREQ <- 1 - seeds$S6.FREQ
seeds$S7.FREQ <- 1 - seeds$S7.FREQ
seeds$S8.FREQ <- 1 - seeds$S8.FREQ
#seeds <- as.data.frame(seeds)
print(head(seeds))

#data <- gather(seeds, Sample, Freq, c(
#    `S1.FREQ`, `S2.FREQ`, `S3.FREQ`, `S4.FREQ`, `S5.FREQ`, `S6.FREQ`, `S7.FREQ`, `S8.FREQ`
#), factor_key=TRUE)
#print(head(data))

bim <- fread( "/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/515g.bim" , sep=" ", header=FALSE )
colnames(bim)[which(names(bim) == "V2")] <- "CHROM_POS"

print("plotting all")

#ggpairs(seeds, columns=`S1.FREQ`:`S8.FREQ`)
ggpairs(seeds, columns=c(
    "S1.FREQ", "S2.FREQ", "S3.FREQ", "S4.FREQ", "S5.FREQ", "S6.FREQ", "S7.FREQ", "S8.FREQ"
))
ggsave("seeds_all_freqs-all.png", bg="white", width=ssize, height=ssize)

print("plotting bona fide")

#bf <- inner_join( seeds, bim, by="CHROM_POS")
#print(head(bf))
bf <- seeds[seeds$CHROM_POS %in% as.list(bim$CHROM_POS) ,]
print(head(bf))

#ggpairs(bf, columns=`S1.FREQ`:`S8.FREQ`)
ggpairs(bf, columns=c(
    "S1.FREQ", "S2.FREQ", "S3.FREQ", "S4.FREQ", "S5.FREQ", "S6.FREQ", "S7.FREQ", "S8.FREQ"
))
ggsave("seeds_all_freqs-bona-fide.png", bg="white", width=ssize, height=ssize)
