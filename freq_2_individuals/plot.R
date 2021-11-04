#Basic plotting
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

print("Reading table")

# Input data. We use a prefix to distinguish between all and the mapq60 filtered data.
# Need to comment/uncomment as needed. Bit ugly, but easier for now.
#infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/frequency-all.csv"
#prefix="all"
infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/frequency-mapq60.csv"
prefix="mapq60"

# Test data
# infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/head.csv"
# infile="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/ath_evo/grenepilot_lucas/freq_2_individuals/test.csv"

# Read the data
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

# =================================================================================================
#     Plot
# =================================================================================================

make_plots <- function(target_dir, data, mytitle, cutoff){
    dir.create(file.path(target_dir), showWarnings = FALSE)

    with_cutoff <- cutoff > 0

    if(TRUE) {

    # -------------------------------------------------------------------------
    #     Dot plot along the genome
    # -------------------------------------------------------------------------

    # For some reason, facet_wrap does not work if it is at the end... so put it after the geom_point
    myplot <- ggplot(data, aes(x=POS, y=FREQ)) +
        geom_hline(yintercept = 0.5, color="gray") +
        geom_point(aes(color=COV)) +
        # geom_hline(yintercept = cutoff, color="red") +
        facet_wrap( ~ CHROM, ncol=1) +
        scale_color_viridis(direction=-1, trans="log10") +
        ylim(0, 1) +
        xlab("Position") +
        ylab("MAF")
        labs(title=mytitle)
        #facet_wrap( ~ CHROM, ncol=1)
        #facet_wrap(vars(CHROM), ncol=1)

    if(with_cutoff){
        myplot <- myplot + geom_hline(yintercept = cutoff, color="red")
    }

    # Somehow, these are saved with a transparent background,
    # so we have to set it explicitly to white here.
    # No idea why we did not need this in other plots...
    # Also, ggsave just insists on print out the size... that clutters our output.
    suppressMessages(ggsave(paste0(target_dir, mytitle,"-dot.png"), bg="white", width=16, height=12))

    # -------------------------------------------------------------------------
    #     Frequency heatmap along the genome
    # -------------------------------------------------------------------------

    myplot <- ggplot(data, aes(x=POS, y=FREQ)) +
        geom_hline(yintercept = 0.5, color="gray") +
        # geom_hline(yintercept = cutoff, color="red") +
        #geom_hex(aes(colour = ..count..), bins=50) +
        geom_bin_2d(bins=c(150,30)) +
        facet_wrap( ~ CHROM, ncol=1) +
        scale_fill_viridis(direction=-1, trans="log10") +
        #scale_color_viridis() +
        ylim(0, 1) +
        xlab("Position") +
        ylab("MAF")
        labs(title=mytitle)

    if(with_cutoff){
        myplot <- myplot + geom_hline(yintercept = cutoff, color="red")
    }

    suppressMessages(ggsave(paste0(target_dir, mytitle,"-heat.png"), bg="white", width=16, height=12))

    # -------------------------------------------------------------------------
    #     Histogram of frequencies
    # -------------------------------------------------------------------------

    myplot <- ggplot(data, aes(x=FREQ)) +
        # geom_vline( xintercept=cutoff, color="red") +
        geom_histogram(bins=50) +
        geom_vline( xintercept=0.5, color="gray") +
        xlab("Frequency") +
        #xlim(0, 1) +
        #geom_vline( xintercept=cutoff, color="red") +
        labs(title=mytitle)

    if(with_cutoff){
        myplot <- myplot + geom_vline( xintercept=cutoff, color="red")
    }

    suppressMessages(ggsave(paste0(target_dir, mytitle,"-hist.png"), bg="white"))

    # -------------------------------------------------------------------------
    #     Histogram of frequencies, log scaled
    # -------------------------------------------------------------------------

    myplot <- ggplot(data, aes(x=FREQ)) +
        # geom_vline( xintercept=cutoff, color="red") +
        geom_histogram(bins=50) +
        geom_vline( xintercept=0.5, color="gray") +
        # coord_trans(y = "log10") +
        scale_y_continuous(trans = "log10") +
        xlab("Frequency") +
        #xlim(0, 1) +
        #geom_vline( xintercept=cutoff, color="red") +
        labs(title=mytitle)

    if(with_cutoff){
        myplot <- myplot + geom_vline( xintercept=cutoff, color="red")
    }

    suppressMessages(ggsave(paste0(target_dir, mytitle,"-hist-log.png"), bg="white"))

    #data <- data[data$FREQ > cutoff, ]

    #ggplot(data, aes(x=FREQ)) +
    #    geom_histogram(bins=50) +
    #    xlab("Frequency") +
    #    xlim(0, 1) +
    #    #geom_vline( xintercept=cutoff, color="red") +
    #    labs(title=mytitle)

    #ggsave(paste0(target_dir, mytitle,"-hist-cut.png"), bg="white")

    }

    # -------------------------------------------------------------------------
    #     Histogram of Coverages
    # -------------------------------------------------------------------------

    # get the data that is not larger than 1000 cov, as we do not want these.
    cov_1000 <- data[data$COV <= 250, ]

    # Histogram of coverages
    ggplot(cov_1000, aes(x=COV)) +
        geom_histogram(bins=50) +
        xlab("Coverage") +
        labs(title=mytitle)

    suppressMessages(ggsave(paste0(target_dir, mytitle,"-cov-hist.png"), bg="white"))

}

# =================================================================================================
#     Main loop over sequencing libraries
# =================================================================================================

for( i in 1:11 ) {

    print(paste0("At ",i, ": ", colnames[i]))
    # print(paste0(colnames[i], ".FREQ"))

    # ggplot is a bit stubborn in selecting columns names with variables,
    # so instead we copy the data to a column with a fixed name first...
    # Also, use this opportunity to convert fom ref freq to MAF.
    # We make a copy of the data, so that we can remove rows without altering
    # the original data. Of course, we also need a trick to make a copy of a df...
    base <- data[, "CHROM", drop=FALSE]
    base$POS <- data$POS
    base$FREQ <- 1.0 - data[, paste0(colnames[i], ".FREQ")]
    base$COV <- data[, paste0(colnames[i], ".COV")]
    base <- base[base$CHROM != "chloroplast" & base$CHROM != "mitochondria", ]

    #print(head(cpy))
    #print(tail(cpy))

    # -------------------------------------------------------------------------
    #     All data
    # -------------------------------------------------------------------------

    print(paste0("    no-cutoff: ", nrow(base)))
    cutoff = 0
    make_plots(paste0(prefix, "-no-cutoff/"), base, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With coverage > 25
    # -------------------------------------------------------------------------

    We want to ignore frequencies below a threshold
    cutoff = 0
    cpy <- base[base$COV >= 25, ]
    cpy <- na.omit(cpy)

    print(paste0("    coverage-25: ", nrow(cpy)))
    make_plots(paste0(prefix, "-coverage-25/"), cpy, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With cutoff at 0.1
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- base[base$FREQ > cutoff, ]
    cpy <- na.omit(cpy)
    
    print(paste0("    cutoff-0.1: ", nrow(cpy)))
    make_plots(paste0(prefix, "-cutoff-0.1/"), cpy, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With coverage 100-200
    # -------------------------------------------------------------------------

    We want to ignore frequencies below a threshold
    cutoff = 0
    cpy <- base[base$COV >= 100 & base$COV <= 200, ]
    cpy <- na.omit(cpy)

    print(paste0("    coverage-100-200: ", nrow(cpy)))
    make_plots(paste0(prefix, "-coverage-100-200/"), cpy, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With cutoff at 0.1 and coverage 100-200
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- base[base$FREQ > cutoff, ]
    cpy <- cpy[cpy$COV >= 100 & cpy$COV <= 200, ]
    cpy <- na.omit(cpy)
    
    print(paste0("    cutoff-0.1 coverage-100-200: ", nrow(cpy)))
    make_plots(paste0(prefix, "-cutoff-0.1_coverage-100-200/"), cpy, colnames[i], cutoff)

}

warnings()
