setwd("~/safedata/ath_evo/grenepilot_lucas/freq_correlation_individuals/")
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(tidyr)
library(RColorBrewer)
library(data.table)
library(dplyr)

# -------------------------------------------------------------------------
#    Load data
# -------------------------------------------------------------------------

print("Load data")

# List of bona fide snps
# bim <- fread( "/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/1001g/515g.bim" , sep=" ", header=FALSE )
bim <- fread( "~/safedata/1001g/515g.bim" , sep=" ", header=FALSE )
bim <- fread( "~/safedata/1001g/1001gbi.bim" , sep=" ", header=FALSE )
colnames(bim)[which(names(bim) == "V2")] <- "CHROM_POS"

seeds <- fread("../freq_seeds_vs_1001g/frequency.csv", sep="\t")
#seeds <- read.csv( "../freq_seeds_vs_1001g/frequency.csv", sep="\t" )
#seeds <- read.csv( "../freq_seeds_vs_1001g/head_seeds.csv", sep="\t" )
# seeds_mer <- dplyr::inner_join(flpools, bim, by=c("CHROM_POS"="V2"))
# seeds_mer<-inner_join(seeds,bim,by=c("CHROM_POS"="CHROM_POS"))


# Replace the `TOTAL` in the column names in that table by `Seeds` to make it intuitive.
colnames(seeds) <- gsub("TOTAL", "Seeds", colnames(seeds))

flpools <- fread( "frequency-flowerpools-mpileup.csv", sep="\t" )
# fpools_mer <- inner_join(flpools,bim,by=c("CHROM_POS"="CHROM_POS"))
#flpools <- read.csv( "frequency-flowerpools-mpileup.csv", sep="\t" )
#flpools <- read.csv( "head-frequency-flowerpools-mpileup.csv", sep="\t" )


flpools<-fpools_mer
seeds<-seeds_mer

#        { "S1",  "Flowerpool1001und2" },
#        { "S2",  "Flowerpool50A" },
#        { "S3",  "Flowerpool50B" },
#        { "S4",  "Flowerpool25A" },
#        { "S5",  "Flowerpool25B" },
#        { "S6",  "Flowerpool10B2" },
#        { "S7",  "Flowerpool5B" },
#        { "S8",  "Plantpool100x" },
#        { "S9",  "Plantpool50a" },
#        { "S10", "Plantpool50b" },
#        { "S11", "Plantpool25a" },
#        { "S12", "Plantpool25b" },
#        { "S13", "Plantpool10b1" },
#        { "S14", "Plantpool10b2" },
#        { "S15", "Plantpool5x" },
#        { "S16", "PlantpoolB12345" }

smpsnums = c(
    "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", 
    "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16"
)

smpnames = c(
    "Flowerpool1001und2", "Flowerpool50A", "Flowerpool50B", "Flowerpool25A", 
    "Flowerpool25B", "Flowerpool10B2", "Flowerpool5B", "Plantpool100x", 
    "Plantpool50a", "Plantpool50b", "Plantpool25a", "Plantpool25b", 
    "Plantpool10b1", "Plantpool10b2", "Plantpool5x", "PlantpoolB12345"
)

# Replace column names in the flpools to have nice readable names.
# We need to do this backwards, to not replace S16 by S1 etc
for( i in length(smpsnums):1 ) {
    colnames(flpools) <- gsub(smpsnums[i], smpnames[i], colnames(flpools))
}

# Merge using chromosome and position
joined <- merge( seeds, flpools, by.x="CHROM_POS", by.y="CHROM_POS", all=TRUE)
joined <- as.data.frame(joined)
#print(head(joined))
# 
# # Remove low cov and invariant sites using seeds.
# # We later also remove low cov and invarants from the flower and leaf pools,
# # but need to do this on a per sample basis, to just remove their particular low cov and invar rows.
# joined <- joined[(joined$Seeds.COV >= 10 & joined$Seeds.REF_CNT > 2 & joined$Seeds.ALT_CNT > 2),]
# #print(head(joined))

# Joined with the 11 M snps

joined_mer <- inner_join( joined, bim, by=c("CHROM_POS"='V2'))

# -------------------------------------------------------------------------
#    make histograms
# -------------------------------------------------------------------------
#if(FALSE){
for( i in 1:16 ) {
    print(paste0( "  At ", smpnames[i] ))

    # Make a copy of the data, containing only the columns we are interested in,
    # so that we can remove the rows that we are not interested in.
    base <- joined[, "CHROM_POS", drop=FALSE]
    base$COV <- joined[, paste0(smpnames[i], ".COV")]
    #base$SEED_FREQ <- joined[, "TOTAL.FREQ.X"]
    #base$SMP_FREQ  <- joined[, paste0("S", i, ".FREQ")]

    # Now drop rows with low cov or invariants in the sample.
    #base <- base[( base$COV >= 25 & base$SMP_FREQ > 0 & base$SMP_FREQ < 1 ),]
    #print(head(base))

    # Full data
    make_cov_hist( base, paste0(smpnames[i], "-all") )

    # Subset of only 1001g bona fide snps
    bf <- inner_join( base, bim, by="CHROM_POS")
    make_cov_hist( bf, paste0(smpnames[i], "-bona-fide") )
}

# -------------------------------------------------------------------------
#    Compute correlations
# -------------------------------------------------------------------------
source("functions.R")

# Plot all
make_corr_plot_quality( joined_mer, "all" )


# print("Plot bona fide")
# 
# # Subset for bona fide snps, and only plot those
# #bf <- inner_join( joined, bim, by="CHROM_POS")
# bf <- joined[joined$CHROM_POS %in% as.list(bim$CHROM_POS) ,]
# make_corr_plot( bf, "bona-fide" )

