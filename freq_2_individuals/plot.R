#Basic plotting
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

print("Reading table")

# Input data.
# infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/frequency.csv"
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/head.csv"
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/test.csv"
data = read.table( infile, sep="\t", header=TRUE )
# head(data)

# Translate indices to original sample names.
# Our input is based on pileup, which does not have sample names,
# so we need this step here.
colnames = c(
    "COL0L1B",
    "COLRUM_POOL1",
    "COLRUM_POOL2",
    "COLRUM_POOL3",
    "GRENE231MM1",
    "GRENE231MM4",
    "GRENE231MM8",
    "UM20L1A",
    "TWOFH1",
    "TWOFH3",
    "TWOFH5"

    #"COL0L1B_CKDL210018053-1a-AK10689-AK9012_HGKVHCCX2_L1",
    #"COLRUM_POOL1_CKDL210018053-1a-AK34070-AK9050_HGKVHCCX2_L1",
    #"COLRUM_POOL2_CKDL210018053-1a-7UDI1539-AK8191_HGKVHCCX2_L1",
    #"COLRUM_POOL3_CKDL210018053-1a-7UDI1971-AK20024_HGKVHCCX2_L1",
    #"GRENE231MM1_CKDL210018053-1a-AK33987-AK9050_HGKVHCCX2_L1",
    #"GRENE231MM4_CKDL210018053-1a-AK28238-AK8191_HGKVHCCX2_L1",
    #"GRENE231MM8_CKDL210018053-1a-AK9787-AK20024_HGKVHCCX2_L1",
    #"UM20L1A_CKDL210018053-1a-AK34068-AK8944_HGKVHCCX2_L1",
    #"TWOFH1_CKDL210018053-1a-AK31398-AK8965_HGKVHCCX2_L1",
    #"TWOFH3_CKDL210018053-1a-AK7838-AK24534_HGKVHCCX2_L1",
    #"TWOFH5_CKDL210018053-1a-AK34069-AK26896_HGKVHCCX2_L1"
)

for( i in 1:11 ) {

    print(paste0("At ",i, ": ", colnames[i]))
    print(paste0(colnames[i], ".FREQ"))

    # ggplot is a bit stubborn in selecting columns names with variables,
    # so instead we copy the data to a column with a fixed name first...
    # But even this is a hassle, and we also need to take into account
    # that R renames our columns with a leading X, as they start with digits.
    # Also, use this opportunity to convert fom ref freq to MAF.
    # We make a copy of the data, so that we can remove rows without altering
    # the original data. Of course, we also need a trick to make a copy of a df...
    cpy <- data[, "CHROM", drop=FALSE]
    cpy$POS <- data$POS
    cpy$FREQ <- 1.0 - data[, paste0(colnames[i], ".FREQ")]
    cpy$COV <- data[, paste0(colnames[i], ".COV")]
    #cpy <- cpy[cpy$CHROM != "chloroplast" & cpy$CHROM != "mitochondria", ]

    #print(head(cpy))
    #print(tail(cpy))

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- cpy[cpy$FREQ > cutoff & cpy$CHROM != "chloroplast" & cpy$CHROM != "mitochondria", ]
    #cpy <- cpy[cpy$CHROM != "chloroplast" & cpy$CHROM != "mitochondria", ]
    #cpy <- cpy[cpy$FREQ > cutoff, ]
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
