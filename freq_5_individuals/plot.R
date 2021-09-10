#Basic plotting
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# Input data
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_5_individuals/frequency.csv"
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_5_individuals/head.csv"
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_5_individuals/test.csv"
data = read.table( infile, sep="\t", header=TRUE )
#head(data)

# Translate indices to original sample names.
# Our input is based on pileup, which does not have sample names,
# so we need this step here.
colnames = c(
  "Flowerpool5B", "Plantpool5x"
  #"Flowerpool1001und2", "Flowerpool50A", "Flowerpool50B", "Flowerpool25A",
  #"Flowerpool25B", "Flowerpool10B2", "Flowerpool5B", "Plantpool100x",
  #"Plantpool50a", "Plantpool50b", "Plantpool25a", "Plantpool25b",
  #"Plantpool10b1", "Plantpool10b2", "Plantpool5x", "PlantpoolB12345"
)

for( i in 1:2 ) {

    print(paste0("At ",i, ": ", colnames[i]))

    # ggplot is a bit stubborn in selecting columns names with variables,
    # so instead we copy the data to a column with a fixed name first...
    # But even this is a hassle, and we also need to take into account
    # that R renames our columns with a leading X, as they start with digits.
    # Also, use this opportunity to convert fom ref freq to MAF.
    # We make a copy of the data, so that we can remove rows without altering
    # the original data. Of course, we also need a trick to make a copy of a df...
    cpy <- data[, "CHROM", drop=FALSE]
    cpy$POS <- data$POS
    cpy$FREQ <- 1.0 - data[, paste0("X", i, ".FREQ")]
    cpy$COV <- data[, paste0("X", i, ".COV")]
    #cpy <- cpy[cpy$CHROM != "chloroplast" & cpy$CHROM != "mitochondria", ]

    #print(head(cpy))
    #print(tail(cpy))

    # We want to ignore frequencies below a threshold, and outside of a good coverage
    cutoff = 0.1
    mincov = 100
    maxcov = 200
    cpy <- cpy[cpy$FREQ > cutoff & cpy$CHROM != "chloroplast" & cpy$CHROM != "mitochondria", ]
    #cpy <- cpy[cpy$CHROM != "chloroplast" & cpy$CHROM != "mitochondria", ]
    #cpy <- cpy[cpy$FREQ > cutoff, ]
    cpy <- cpy[cpy$COV >= mincov & cpy$COV <= maxcov, ]
    cpy <- na.omit(cpy)

    # For some reason, facet_wrap does not work if it is at the end... so put it after the geom_point
    ggplot(cpy, aes(x=POS, y=FREQ)) + 
        geom_hline(yintercept = 0.5, color="gray") + 
        geom_hline(yintercept = cutoff, color="red") +
        geom_point(aes(color=COV)) + 
	facet_wrap( ~ CHROM, ncol=1) +
        scale_color_viridis(direction=-1, trans="log10") +
        ylim(0, 1) +
        xlab("Position") +
        ylab("MAF")
        labs(title=colnames[i])
	#facet_wrap( ~ CHROM, ncol=1)
        #facet_wrap(vars(CHROM), ncol=1)

    # Somehow, these are saved with a transparent background, 
    # so we have to set it explicitly to white here.
    # No idea why we did not need this in other plots...
    ggsave(paste0(colnames[i],"-dot.png"), bg="white", width=16, height=12)

    ggplot(cpy, aes(x=POS, y=FREQ)) +
        geom_hline(yintercept = 0.5, color="gray") +
        geom_hline(yintercept = cutoff, color="red") +
        #geom_hex(aes(colour = ..count..), bins=50) +
	geom_bin_2d(bins=c(150,30)) +
	facet_wrap( ~ CHROM, ncol=1) +
	scale_fill_viridis(direction=-1, trans="log10") +
	#scale_color_viridis() +
	ylim(0, 1) +
        xlab("Position") +
        ylab("MAF")
        labs(title=colnames[i])

    ggsave(paste0(colnames[i],"-hex.png"), bg="white", width=16, height=12)

    ggplot(cpy, aes(x=FREQ)) + 
        geom_vline( xintercept=cutoff, color="red") +
        geom_histogram(bins=50) +
        geom_vline( xintercept=0.5, color="gray") +
	xlab("Frequency") +
	#xlim(0, 1) +
	#geom_vline( xintercept=cutoff, color="red") +
	labs(title=colnames[i])

    ggsave(paste0(colnames[i],"-hist.png"), bg="white")

    #cpy <- cpy[cpy$FREQ > cutoff, ]

    #ggplot(cpy, aes(x=FREQ)) +
    #    geom_histogram(bins=50) +
    #    xlab("Frequency") +
    #    xlim(0, 1) +
    #    #geom_vline( xintercept=cutoff, color="red") +
    #    labs(title=colnames[i])

    #ggsave(paste0(colnames[i],"-hist-cut.png"), bg="white")

}
