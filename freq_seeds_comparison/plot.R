# Basic plotting
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
suppressMessages(library(viridis))

# R is really weird concerning arrays/lists and index-based access...
# See https://stackoverflow.com/a/43534546/4184258

print("Reading tables")

# Input files
vcf_freq_files <- c(
    "bcftools-c-filtered",
    "bcftools-m-filtered",
    "freebayes-kv-filtered",
    "freebayes-nkv-filtered",
    "haplotypecaller-kv-filtered",
    "haplotypecaller-nkv-filtered"
)
pileup_freq_file = "frequency-ath-evo-seeds-pileups-all-merged-samples.csv"

# Read all vcf freq tables.
# We do it all at once, in preparation for later, when we probably want to make
# one big facet plot with all of them in one figure.
vcf_freq_data <- list()
i <- 1
for (f in vcf_freq_files) {
    fn <- paste("frequency-ath-evo-seeds-", f, ".csv", sep="")
    vcf_freq_data[[i]] <- read.table(fn, sep="\t", header=TRUE)
    i <- i+1
}
# head(vcf_freq_data)
# head(vcf_freq_data[[1]])

# Read pileup freq table.
pileup_freq_data <- read.table(pileup_freq_file, sep="\t", header=TRUE)
# head(pileup_freq_data)

# =================================================================================================
#     Plot
# =================================================================================================

make_plots <- function(target_dir, joined_data, mytitle){
    dir.create(file.path(target_dir), showWarnings = FALSE)

    # Get standard deviation of the diff
    stddev <- sd(joined_data$DIFF)

    # allow to skip this part
    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Scatter plot seeds vs caller
        # -------------------------------------------------------------------------

        # Plot the scatter plot of the frequencies.
        # We want to use different dot sizes to see how this affects the plot.
        #for( ds in c(0.01, 0.03, 0.1, 0.3, 1.0)) {
        # for( ds in c(0.1)) {
        # }
        ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
            geom_vline( xintercept=0, color="gray") +
            geom_hline( yintercept=0, color="gray") +
            geom_point(alpha = 1/80, size=0.1) +
            # geom_point(alpha = 1/100, size=ds) +
            #geom_point(alpha = 1/100, size=ds, aes(color=COV)) +
            #scale_color_viridis(trans = "log10", direction=-1) +
            geom_abline(slope = 1, colour="blue") +
            coord_fixed() +
            xlab("Pileup MAF") +
            ylab("VCF MAF") +
            labs(title=mytitle) +
            annotate("text", x = 0.02, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

        # We save different sizes, as the pixel density changes.
        suppressMessages(ggsave(paste0(target_dir, mytitle, "-scatter.png"), bg="white"))
        #ggsave(paste(f,"-large.png", sep=""), width=16, height=15)

        # -------------------------------------------------------------------------
        #     Hex scatter plot seeds vs caller
        # -------------------------------------------------------------------------

        # geom hex produces ugly white dots in the output,
        # so we have to force it to draw borders of the same color
        # around each hexagon... see https://github.com/tidyverse/ggplot2/issues/2157
        # and https://stackoverflow.com/questions/52006888/ggplot2-geom-hex-white-border
        ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
            geom_hex(aes(colour = ..count..), bins=50) +
            #scale_fill_viridis(direction=-1) +
            scale_fill_viridis() +
            scale_color_viridis() +
            coord_fixed() +
            xlab("Pileup MAF") +
            ylab("VCF MAF") +
            labs(title=mytitle)

        suppressMessages(ggsave(paste0(target_dir, mytitle,"-hex.png"), width=9, height=7, bg="white"))

        ggplot(joined_data, aes(x=TOTAL.FREQ.x, y=TOTAL.FREQ.y)) +
            geom_hex(aes(colour = ..count..), bins=50) +
            #scale_fill_viridis(trans = "log10", direction=-1) +
            scale_fill_viridis(trans = "log10") +
            scale_color_viridis(trans="log10") +
            coord_fixed() +
            xlab("Pileup MAF") +
            ylab("VCF MAF") +
            labs(title=mytitle)

        suppressMessages(ggsave(paste0(target_dir, mytitle,"-hex-log.png"), width=9, height=7, bg="white"))

    }

    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Histogram freq diff seeds and caller
        # -------------------------------------------------------------------------

        # Histogram of diff between the two
        ggplot(joined_data, aes(x=DIFF)) +
            geom_histogram(bins=50) +
            xlab("Frequency Diff") +
            labs(title=mytitle) +
            annotate("text", x = -1, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

        suppressMessages(ggsave(paste0(target_dir, mytitle,"-hist.png"), bg="white"))

        # Histogram of diff between the two, log scaled
        ggplot(joined_data, aes(x=DIFF)) +
            geom_histogram(bins=50) +
            scale_y_continuous(trans = "log10") +
            xlab("Frequency Diff") +
            labs(title=mytitle) +
            annotate("text", x = -1, y = Inf, hjust = 0, vjust=1, label = paste0("SD: ", format(stddev, digits=3)))

        suppressMessages(ggsave(paste0(target_dir, mytitle,"-hist-log.png"), bg="white"))

    }

    if(TRUE) {

        # -------------------------------------------------------------------------
        #     Histogram frequency of caller
        # -------------------------------------------------------------------------

        # Histogram of frequencies
        ggplot(joined_data, aes(x=TOTAL.FREQ.y)) +
            geom_histogram(bins=50) +
            xlab("Frequency") +
            labs(title=mytitle)

        suppressMessages(ggsave(paste0(target_dir, mytitle,"-freq-hist.png"), bg="white"))
    }
}

# =================================================================================================
#     Main loop over tables
# =================================================================================================

# Go through all freq tables and plot their frequency vs the pileup.
i <- 1
for (f in vcf_freq_files) {

    print(paste0("At: ", f))

    # Make a joined table that only contains the rows of SNPs
    # that occur in both the pileup and the vcf.
    joined_data <- merge(
        pileup_freq_data, vcf_freq_data[[i]],
        by.x = "CHROM_POS", by.y = "CHROM_POS", all=FALSE
    )
    joined_data <- na.omit(joined_data)
    joined_data$COV <- joined_data$TOTAL.COV.x + joined_data$TOTAL.COV.y
    joined_data$SEED.COV <- joined_data$TOTAL.COV.x

    # Change ref freq to minor freq
    joined_data$TOTAL.FREQ.x <- 1.0 - joined_data$TOTAL.FREQ.x
    joined_data$TOTAL.FREQ.y <- 1.0 - joined_data$TOTAL.FREQ.y

    joined_data$DIFF <- joined_data$TOTAL.FREQ.y - joined_data$TOTAL.FREQ.x

    # -------------------------------------------------------------------------
    #     All data
    # -------------------------------------------------------------------------

    print(paste0("    all: ", nrow(joined_data)))
    make_plots("all/", joined_data, f)

    # -------------------------------------------------------------------------
    #     COV >= 50
    # -------------------------------------------------------------------------

    joined_data <- joined_data[joined_data$SEED.COV >= 50, ]
    print(paste0("    coverage-50: ", nrow(joined_data)))
    make_plots("coverage-50/", joined_data, f)

    # -------------------------------------------------------------------------
    #     COV >= 100 & COV <= 200
    # -------------------------------------------------------------------------

    joined_data <- joined_data[joined_data$SEED.COV >= 100 & joined_data$SEED.COV <= 200, ]
    print(paste0("    coverage-100-200: ", nrow(joined_data)))
    make_plots("coverage-100-200/", joined_data, f)

    i <- i+1
}

warnings()
