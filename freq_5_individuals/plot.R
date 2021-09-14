#Basic plotting
library(ggplot2)
library(cowplot)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

print("Reading table")

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

snames = c(
    "S7", "S15"
)

# =================================================================================================
#     Plot
# =================================================================================================

make_plots <- function(target_dir, data, mytitle, cutoff){
    dir.create(file.path(target_dir), showWarnings = FALSE)

    with_cutoff <- cutoff > 0

    # -------------------------------------------------------------------------
    #     Dot plot along the genome
    # -------------------------------------------------------------------------

    # For some reason, facet_wrap does not work if it is at the end... so put it after the geom_point
    myplot <- ggplot(data, aes(x=POS, y=FREQ)) +
        geom_hline(yintercept = 0.2, color="gray") +
        geom_hline(yintercept = 0.4, color="gray") +
        geom_hline(yintercept = 0.6, color="gray") +
        geom_hline(yintercept = 0.8, color="gray") +
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
        geom_hline(yintercept = 0.2, color="gray") +
        geom_hline(yintercept = 0.4, color="gray") +
        geom_hline(yintercept = 0.6, color="gray") +
        geom_hline(yintercept = 0.8, color="gray") +
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
        geom_histogram(bins=100) +
        geom_vline( xintercept=0.2, color="gray") +
        geom_vline( xintercept=0.4, color="gray") +
        geom_vline( xintercept=0.6, color="gray") +
        geom_vline( xintercept=0.8, color="gray") +
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
        geom_histogram(bins=100) +
        geom_vline( xintercept=0.2, color="gray") +
        geom_vline( xintercept=0.4, color="gray") +
        geom_vline( xintercept=0.6, color="gray") +
        geom_vline( xintercept=0.8, color="gray") +
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
    #    labs(title=colnames[i])

    #ggsave(paste0(target_dir, colnames[i],"-hist-cut.png"), bg="white")
}

# =================================================================================================
#     Main loop over sequencing libraries
# =================================================================================================

for( i in 1:2 ) {

    print(paste0("At ",i, ": ", colnames[i]))
    # print(paste0(colnames[i], ".FREQ"))

    # ggplot is a bit stubborn in selecting columns names with variables,
    # so instead we copy the data to a column with a fixed name first...
    # Also, use this opportunity to convert fom ref freq to MAF.
    # We make a copy of the data, so that we can remove rows without altering
    # the original data. Of course, we also need a trick to make a copy of a df...
    base <- data[, "CHROM", drop=FALSE]
    base$POS <- data$POS
    base$FREQ <- 1.0 - data[, paste0(snames[i], ".FREQ")]
    base$COV <- data[, paste0(snames[i], ".COV")]
    base <- base[base$CHROM != "chloroplast" & base$CHROM != "mitochondria", ]

    #print(head(cpy))
    #print(tail(cpy))

    # -------------------------------------------------------------------------
    #     All data
    # -------------------------------------------------------------------------

    print(paste0("    no-cutoff: ", nrow(base)))
    cutoff = 0
    make_plots("no-cutoff/", base, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With cutoff at 0.1
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- base[base$FREQ > cutoff, ]
    cpy <- na.omit(cpy)

    print(paste0("    cutoff-0.1: ", nrow(cpy)))
    make_plots("cutoff-0.1/", cpy, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With coverage 100-200
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0
    cpy <- base[base$COV >= 100 & base$COV <= 200, ]
    cpy <- na.omit(cpy)

    print(paste0("    coverage-100-200: ", nrow(cpy)))
    make_plots("coverage-100-200/", cpy, colnames[i], cutoff)

    # -------------------------------------------------------------------------
    #     With cutoff at 0.1 and coverage 100-200
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- base[base$FREQ > cutoff, ]
    cpy <- cpy[cpy$COV >= 100 & cpy$COV <= 200, ]
    cpy <- na.omit(cpy)

    print(paste0("    cutoff-0.1 coverage-100-200: ", nrow(cpy)))
    make_plots("cutoff-0.1_coverage-100-200/", cpy, colnames[i], cutoff)

}

warnings()
