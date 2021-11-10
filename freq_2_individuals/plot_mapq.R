#Basic plotting
library(ggplot2)
library(cowplot)
library(data.table)
#library(tidyverse)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# =================================================================================================
#     Plot
# =================================================================================================

make_plots <- function(target_dir, data, mytitle, cutoff = 0){
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
#     Main loop over input files
# =================================================================================================

filenames = c(
    "q0",
    "q20",
    "q40",
    "q60",
    "q0-f",
    "q20-f",
    "q40-f",
    "q60-f",
    "Q30-q0",
    "Q30-q20",
    "Q30-q40",
    "Q30-q60",
    "Q30-q0-f",
    "Q30-q20-f",
    "Q30-q40-f",
    "Q30-q60-f"
)

for( i in 1:16 ) {
    
    print(paste0("At ",i, ": ", filenames[i])) 
    prefix <- filenames[i]
    base <- fread( paste0( "data/frequency-", filenames[i], ".csv" ), sep="\t", header=TRUE )

    base$FREQ <- 1.0 - base[, "S1.FREQ"]
    colnames(base)[which(names(base) == "S1.COV")] <- "COV"
    base <- base[base$CHROM != "chloroplast" & base$CHROM != "mitochondria", ]

    # -------------------------------------------------------------------------
    #     All data
    # -------------------------------------------------------------------------

    print(paste0("    no-cutoff: ", nrow(base)))
    cutoff = 0
    make_plots(paste0(prefix, "-no-cutoff/"), base, prefix, cutoff)

    # -------------------------------------------------------------------------
    #     With coverage > 25
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0
    cpy <- base[base$COV >= 25, ]
    cpy <- na.omit(cpy)

    print(paste0("    coverage-25: ", nrow(cpy)))
    make_plots(paste0(prefix, "-coverage-25/"), cpy, prefix, cutoff)

    # -------------------------------------------------------------------------
    #     With cutoff at 0.1
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- base[base$FREQ > cutoff, ]
    cpy <- na.omit(cpy)
    
    print(paste0("    cutoff-0.1: ", nrow(cpy)))
    make_plots(paste0(prefix, "-cutoff-0.1/"), cpy, prefix, cutoff)

    # -------------------------------------------------------------------------
    #     With coverage 100-200
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0
    cpy <- base[base$COV >= 100 & base$COV <= 200, ]
    cpy <- na.omit(cpy)

    print(paste0("    coverage-100-200: ", nrow(cpy)))
    make_plots(paste0(prefix, "-coverage-100-200/"), cpy, prefix, cutoff)

    # -------------------------------------------------------------------------
    #     With cutoff at 0.1 and coverage 100-200
    # -------------------------------------------------------------------------

    # We want to ignore frequencies below a threshold
    cutoff = 0.1
    cpy <- base[base$FREQ > cutoff, ]
    cpy <- cpy[cpy$COV >= 100 & cpy$COV <= 200, ]
    cpy <- na.omit(cpy)
    
    print(paste0("    cutoff-0.1 coverage-100-200: ", nrow(cpy)))
    make_plots(paste0(prefix, "-cutoff-0.1_coverage-100-200/"), cpy, prefix, cutoff)

}

